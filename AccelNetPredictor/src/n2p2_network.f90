module n2p2_network
    use iso_fortran_env, only: real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use aenet_network, only: atomic_network, ACTIVATION_LINEAR, ACTIVATION_TANH, &
        ACTIVATION_LOGISTIC, ACTIVATION_SOFTPLUS, ACTIVATION_RELU, ACTIVATION_GAUSSIAN, &
        ACTIVATION_COSINE, ACTIVATION_REVERSE_LOGISTIC, ACTIVATION_EXPONENTIAL, ACTIVATION_HARMONIC
    use accelnet_setup, only: descriptor_setup
    use accelnet_descriptor_models, only: add_behler
    use accelnet_behler, only: behler_config, initialize_behler_config, add_g2, add_g4, add_g5, add_extended, compact_subtype
    implicit none
    private

    integer, parameter :: LINE_LENGTH = 4096, MAX_LAYERS = 64

    type :: symmetry_function
        integer :: kind = 0, neighbor1 = 0, neighbor2 = 0
        character(len=3) :: subtype=""
        real(real64) :: left=0, right=0
        real(real64) :: eta = 0.0_real64, shift = 0.0_real64
        real(real64) :: lambda = 0.0_real64, zeta = 0.0_real64, cutoff = 0.0_real64
    end type symmetry_function

    type :: symmetry_function_list
        type(symmetry_function), allocatable :: values(:)
    end type symmetry_function_list

    type :: network_topology
        integer :: hidden_layers = -1
        integer :: hidden_nodes(MAX_LAYERS) = 0
        character(len=1) :: activations(MAX_LAYERS) = ""
        integer :: node_count = 0, activation_count = 0
        logical :: hidden_layers_set = .false.
        logical :: nodes_set = .false.
        logical :: activations_set = .false.
    end type network_topology

    type :: n2p2_settings
        character(len=16), allocatable :: species(:)
        type(symmetry_function_list), allocatable :: functions(:)
        type(network_topology) :: global_topology
        type(network_topology), allocatable :: topologies(:)
        integer :: cutoff_type = -1
        real(real64) :: cutoff_alpha = 0.0_real64
        logical :: scale = .false., center = .false., sigma_scale = .false.
        logical :: normalize_nodes = .false., has_elements = .false.
        real(real64) :: scale_min = 0.0_real64, scale_max = 1.0_real64
        real(real64) :: mean_energy = 0.0_real64, conv_energy = 1.0_real64
        real(real64) :: conv_length = 1.0_real64
        logical :: has_mean_energy = .false., has_conv_energy = .false., has_conv_length = .false.
        real(real64), allocatable :: atomic_references(:)
    end type n2p2_settings

    public :: load_n2p2_model, load_n2p2_setups, atomic_number

contains

    subroutine load_n2p2_model(directory, networks, setups, species_names, status, message)
        integer, optional, intent(out) :: status
        character(len=*), optional, intent(out) :: message
        character(len=512) :: error
        character(len=*), intent(in) :: directory
        type(atomic_network), allocatable, intent(out) :: networks(:)
        type(descriptor_setup), allocatable, intent(out) :: setups(:)
        character(len=16), allocatable, intent(out) :: species_names(:)
        type(n2p2_settings) :: settings
        character(len=LINE_LENGTH) :: input_file, scaling_file, weight_file
        integer :: species

        error = ''
        if (present(status)) status = 0
        if (present(message)) message = ''
        if (len_trim(directory) > LINE_LENGTH-32) then
            error = 'n2p2 directory path is too long'
            go to 900
        end if
        input_file = join_path(directory, "input.nn")
        scaling_file = join_path(directory, "scaling.data")
        call read_validated_settings(trim(input_file), settings, error)
        if (len_trim(error) /= 0) go to 900
        species_names = settings%species
        allocate(networks(size(species_names)), setups(size(species_names)))
        do species = 1, size(species_names)
            call sort_functions(settings%functions(species)%values)
            weight_file = join_path(directory, weight_filename(species_names(species)))
            call build_network(settings, species, trim(weight_file), trim(scaling_file), networks(species), error)
            if (len_trim(error) /= 0) go to 900
            call build_setup(settings, species, setups(species))
        end do
        return
900     if (allocated(setups)) deallocate(setups)
        if (allocated(species_names)) deallocate(species_names)
        if (allocated(networks)) deallocate(networks)
        if (present(message)) message = error
        if (present(status)) then
            status = 1
        else
            error stop trim(error)
        end if
    end subroutine load_n2p2_model

    subroutine load_n2p2_setups(directory, setups, species_names, status, message)
        integer, optional, intent(out) :: status
        character(len=*), optional, intent(out) :: message
        character(len=512) :: error
        character(len=*), intent(in) :: directory
        type(descriptor_setup), allocatable, intent(out) :: setups(:)
        character(len=16), allocatable, intent(out) :: species_names(:)
        type(n2p2_settings) :: settings
        character(len=LINE_LENGTH) :: input_file
        integer :: species

        error = ''
        if (present(status)) status = 0
        if (present(message)) message = ''
        if (len_trim(directory) > LINE_LENGTH-32) then
            error = 'n2p2 directory path is too long'
            go to 900
        end if
        input_file = join_path(directory, "input.nn")
        call read_validated_settings(trim(input_file), settings, error)
        if (len_trim(error) /= 0) go to 900
        species_names = settings%species
        allocate(setups(size(species_names)))
        do species = 1, size(species_names)
            call sort_functions(settings%functions(species)%values)
            call build_setup(settings, species, setups(species))
        end do
        return
900     if (allocated(setups)) deallocate(setups)
        if (allocated(species_names)) deallocate(species_names)
        if (present(message)) message = error
        if (present(status)) then
            status = 1
        else
            error stop trim(error)
        end if
    end subroutine load_n2p2_setups

    subroutine read_validated_settings(filename, settings, error)
        character(len=*), intent(out) :: error
        character(len=*), intent(in) :: filename
        type(n2p2_settings), intent(out) :: settings
        error = ''
        call read_settings(filename, settings, error)
        if (len_trim(error) /= 0) return
        call resolve_topologies(settings)
        call validate_settings(settings, error)
        if (len_trim(error) /= 0) return
    end subroutine read_validated_settings

    subroutine read_settings(filename, settings, error)
        character(len=*), intent(out) :: error
        logical :: opened
        character(len=*), intent(in) :: filename
        type(n2p2_settings), intent(out) :: settings
        character(len=LINE_LENGTH) :: line, key, rest, topology_values
        character(len=16) :: central, neighbor1, neighbor2
        type(symmetry_function) :: sf
        integer :: unit, ios, count, central_index
        real(real64) :: reference

        error = ''
        opened = .false.
        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            error = "cannot open n2p2 input.nn"
            go to 900
        end if
        opened = .true.
        do
            read(unit, "(A)", iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) then
                error = 'I/O error reading n2p2 model file'
                go to 900
            end if
            call clean_line(line)
            if (len_trim(line) == 0) cycle
            call split_key(line, key, rest)
            select case(trim(key))
            case("number_of_elements")
                if (allocated(settings%species)) then
                    error = 'duplicate n2p2 number_of_elements'
                    go to 900
                end if
                count = 0
                read(rest, *, iostat=ios) count
                if (ios /= 0 .or. count < 1 .or. count > 118) then
                    error = "invalid n2p2 number_of_elements"
                    go to 900
                end if
                allocate(settings%species(count), settings%functions(count), settings%topologies(count), &
            settings%atomic_references(count))
                settings%atomic_references = 0.0_real64
                settings%species = ""
            case("elements")
                if (.not. allocated(settings%species)) then
                    error = "n2p2 elements must follow number_of_elements"
                    go to 900
                end if
                if (settings%has_elements) then
                    error = 'duplicate n2p2 elements'
                    go to 900
                end if
                read(rest, *, iostat=ios) settings%species
                if (ios /= 0) then
                    error = "invalid n2p2 elements"
                    go to 900
                end if
                do count = 1,size(settings%species)
                    if (atomic_number(settings%species(count)) == 0) then
                        error = 'unknown chemical element in n2p2 model'
                        go to 900
                    end if
                    if (any(settings%species(:count-1) == settings%species(count))) then
                        error = 'duplicate chemical element in n2p2 model'
                        go to 900
                    end if
                end do
                call sort_species_by_atomic_number(settings%species)
                settings%has_elements = .true.
            case("nnp_type")
                if (trim(lowercase(adjustl(rest))) /= "2g" .and. trim(lowercase(adjustl(rest))) /= "2g-hdnnp" &
            .and. trim(adjustl(rest)) /= "2") then
                    error = "AccelNetPredictor supports only n2p2 2G models"
                    go to 900
                end if
            case("cutoff_type")
                read(rest, *, iostat=ios) settings%cutoff_type, settings%cutoff_alpha
                if (ios /= 0) then
                    settings%cutoff_alpha = 0.0_real64
                    read(rest, *, iostat=ios) settings%cutoff_type
                end if
                if (ios /= 0) then
                    error = "invalid n2p2 cutoff_type"
                    go to 900
                end if
            case("cutoff_alpha")
                read(rest, *, iostat=ios) settings%cutoff_alpha
                if (ios /= 0) then
                    error = "invalid n2p2 cutoff_alpha"
                    go to 900
                end if
            case("global_hidden_layers_short")
                read(rest, *, iostat=ios) settings%global_topology%hidden_layers
                if (ios /= 0) then
                    error = "invalid n2p2 global_hidden_layers_short"
                    go to 900
                end if
                settings%global_topology%hidden_layers_set = .true.
            case("global_nodes_short")
                call parse_integer_list(rest, settings%global_topology%hidden_nodes, &
            settings%global_topology%node_count, error)
                if (len_trim(error) /= 0) go to 900
                settings%global_topology%nodes_set = .true.
            case("global_activation_short")
                call parse_character_list(rest, settings%global_topology%activations, &
            settings%global_topology%activation_count, error)
                if (len_trim(error) /= 0) go to 900
                settings%global_topology%activations_set = .true.
            case("element_hidden_layers_short")
                if (.not. settings%has_elements) then
                    error = "n2p2 element topology precedes elements"
                    go to 900
                end if
                call split_key(rest, central, topology_values)
                central_index = species_index(central, settings%species)
                if (central_index == 0) then
                    error = "unknown species in n2p2 element topology"
                    go to 900
                end if
                read(topology_values, *, iostat=ios) settings%topologies(central_index)%hidden_layers
                if (ios /= 0) then
                    error = "invalid n2p2 element_hidden_layers_short"
                    go to 900
                end if
                settings%topologies(central_index)%hidden_layers_set = .true.
            case("element_nodes_short")
                if (.not. settings%has_elements) then
                    error = "n2p2 element topology precedes elements"
                    go to 900
                end if
                call split_key(rest, central, topology_values)
                central_index = species_index(central, settings%species)
                if (central_index == 0) then
                    error = "unknown species in n2p2 element topology"
                    go to 900
                end if
                call parse_integer_list(topology_values, settings%topologies(central_index)%hidden_nodes, &
            settings%topologies(central_index)%node_count, error)
                if (len_trim(error) /= 0) go to 900
                settings%topologies(central_index)%nodes_set = .true.
            case("element_activation_short")
                if (.not. settings%has_elements) then
                    error = "n2p2 element topology precedes elements"
                    go to 900
                end if
                call split_key(rest, central, topology_values)
                central_index = species_index(central, settings%species)
                if (central_index == 0) then
                    error = "unknown species in n2p2 element topology"
                    go to 900
                end if
                call parse_character_list(topology_values, settings%topologies(central_index)%activations, &
            settings%topologies(central_index)%activation_count, error)
                if (len_trim(error) /= 0) go to 900
                settings%topologies(central_index)%activations_set = .true.
            case("normalize_nodes")
                settings%normalize_nodes = .true.
            case("scale_symmetry_functions")
                settings%scale = .true.
            case("center_symmetry_functions")
                settings%center = .true.
            case("scale_symmetry_functions_sigma")
                settings%sigma_scale = .true.
            case("scale_min_short")
                read(rest, *, iostat=ios) settings%scale_min
                if (ios /= 0) then
                    error = "invalid n2p2 scale_min_short"
                    go to 900
                end if
            case("scale_max_short")
                read(rest, *, iostat=ios) settings%scale_max
                if (ios /= 0) then
                    error = "invalid n2p2 scale_max_short"
                    go to 900
                end if
            case("mean_energy")
                read(rest, *, iostat=ios) settings%mean_energy
                if (ios /= 0) then
                    error = "invalid n2p2 mean_energy"
                    go to 900
                end if
                settings%has_mean_energy = .true.
            case("conv_energy")
                read(rest, *, iostat=ios) settings%conv_energy
                if (ios /= 0) then
                    error = "invalid n2p2 conv_energy"
                    go to 900
                end if
                settings%has_conv_energy = .true.
            case("conv_length")
                read(rest, *, iostat=ios) settings%conv_length
                if (ios /= 0) then
                    error = "invalid n2p2 conv_length"
                    go to 900
                end if
                settings%has_conv_length = .true.
            case("atom_energy")
                if (.not. settings%has_elements) then
                    error = "n2p2 atom_energy precedes elements"
                    go to 900
                end if
                read(rest, *, iostat=ios) central, reference
                if (ios /= 0) then
                    error = "invalid n2p2 atom_energy"
                    go to 900
                end if
                central_index = species_index(central, settings%species)
                if (central_index == 0) then
                    error = "unknown species in n2p2 atom_energy"
                    go to 900
                end if
                settings%atomic_references(central_index) = reference
            case("symfunction_short")
                if (.not. settings%has_elements) then
                    error = "n2p2 symmetry functions precede elements"
                    go to 900
                end if
                sf = symmetry_function()
                central = ""; neighbor1 = ""; neighbor2 = ""
                read(rest, *, iostat=ios) central, sf%kind
                if (ios /= 0) then
                    error = "invalid n2p2 symfunction_short"
                    go to 900
                end if
                central_index = species_index(central, settings%species)
                if (central_index == 0) then
                    error = "unknown central species in n2p2 symmetry function"
                    go to 900
                end if
                select case(sf%kind)
                case(2)
                    read(rest, *, iostat=ios) central, sf%kind, neighbor1, sf%eta, sf%shift, sf%cutoff
                    sf%neighbor1 = species_index(neighbor1, settings%species)
                case(3, 9)
                    read(rest, *, iostat=ios) central, sf%kind, neighbor1, neighbor2, sf%eta, sf%lambda, sf%zeta, &
            sf%cutoff, sf%shift
                    if (ios /= 0) then
                        sf%shift = 0.0_real64
                        read(rest, *, iostat=ios) central, sf%kind, neighbor1, neighbor2, sf%eta, sf%lambda, &
            sf%zeta, sf%cutoff
                    end if
                    sf%neighbor1 = species_index(neighbor1, settings%species)
                    sf%neighbor2 = species_index(neighbor2, settings%species)
                    if (sf%neighbor1 > sf%neighbor2) then
                        count = sf%neighbor1; sf%neighbor1 = sf%neighbor2; sf%neighbor2 = count
                    end if
                case(12)
                    read(rest,*,iostat=ios) central,sf%kind,sf%eta,sf%shift,sf%cutoff
                case(13)
                    read(rest,*,iostat=ios) central,sf%kind,sf%eta,sf%shift,sf%lambda,sf%zeta,sf%cutoff
                case(20)
                    read(rest,*,iostat=ios) central,sf%kind,neighbor1,sf%shift,sf%cutoff,sf%subtype
                    sf%neighbor1=species_index(neighbor1,settings%species)
                case(21,22)
                    read(rest,*,iostat=ios) central,sf%kind,neighbor1,neighbor2,sf%shift,sf%cutoff, &
            sf%left,sf%right,sf%subtype
                    sf%neighbor1=species_index(neighbor1,settings%species)
                    sf%neighbor2=species_index(neighbor2,settings%species)
                    if (sf%neighbor1 > sf%neighbor2) then
                        count=sf%neighbor1; sf%neighbor1=sf%neighbor2; sf%neighbor2=count
                    end if
                case(23)
                    read(rest,*,iostat=ios) central,sf%kind,sf%shift,sf%cutoff,sf%subtype
                case(24,25)
                    read(rest,*,iostat=ios) central,sf%kind,sf%shift,sf%cutoff,sf%left,sf%right,sf%subtype
                case default
                    error = "unsupported n2p2 symmetry-function type"
                    go to 900
                end select
                if (ios /= 0) then
                    error = "invalid n2p2 symmetry-function parameters"
                    go to 900
                end if
                if (sf%kind < 12 .or. (sf%kind >= 20 .and. sf%kind <= 22)) then
                    if (sf%neighbor1 == 0) then
                        error = 'unknown n2p2 neighbor species'
                        go to 900
                    end if
                    if (sf%kind /= 2 .and. sf%kind /= 20 .and. sf%neighbor2 == 0) then
                        error = 'unknown n2p2 angular species'
                        go to 900
                    end if
                end if
                call append_function(settings%functions(central_index), sf)
            case default
                if (index(trim(key), "global_") == 1 .or. index(trim(key), "element_") == 1 .or. index(trim(key), &
            "symfunction_") == 1 .or. index(trim(key), "nnp_type") == 1) then
                    error = "unsupported n2p2 model setting"
                    go to 900
                end if
            end select
        end do
        close(unit)
        opened = .false.
900     if (opened) close(unit)
    end subroutine read_settings

    subroutine resolve_topologies(settings)
        type(n2p2_settings), intent(inout) :: settings
        integer :: species

        if (.not. allocated(settings%topologies)) return
        do species = 1, size(settings%topologies)
            if (.not. settings%topologies(species)%hidden_layers_set .and. &
                settings%global_topology%hidden_layers_set) then
                settings%topologies(species)%hidden_layers = settings%global_topology%hidden_layers
                settings%topologies(species)%hidden_layers_set = .true.
            end if
            if (.not. settings%topologies(species)%nodes_set .and. settings%global_topology%nodes_set) then
                settings%topologies(species)%hidden_nodes = settings%global_topology%hidden_nodes
                settings%topologies(species)%node_count = settings%global_topology%node_count
                settings%topologies(species)%nodes_set = .true.
            end if
            if (.not. settings%topologies(species)%activations_set .and. &
                settings%global_topology%activations_set) then
                settings%topologies(species)%activations = settings%global_topology%activations
                settings%topologies(species)%activation_count = settings%global_topology%activation_count
                settings%topologies(species)%activations_set = .true.
            end if
        end do
    end subroutine resolve_topologies

    subroutine validate_settings(settings, error)
        type(n2p2_settings), intent(in) :: settings
        character(len=*), intent(out) :: error
        integer :: species,j,n
        type(symmetry_function) :: sf
        error = ''
        if (.not. settings%has_elements) then
            error = 'n2p2 input.nn has no elements'
            return
        end if
        if (.not. all(ieee_is_finite([settings%cutoff_alpha,settings%scale_min,settings%scale_max, &
            settings%mean_energy,settings%conv_energy,settings%conv_length])) .or. &
            .not. all(ieee_is_finite(settings%atomic_references))) then
            error = 'nonfinite n2p2 normalization/cutoff parameters'
            return
        end if
        if (settings%cutoff_type < 0 .or. settings%cutoff_type > 9 .or. &
            settings%cutoff_alpha < 0 .or. settings%cutoff_alpha >= 1 .or. &
            (settings%cutoff_type == 9 .and. settings%cutoff_alpha <= 0)) then
            error = 'invalid n2p2 cutoff parameters'
            return
        end if
        if (settings%sigma_scale .and. (settings%scale .or. settings%center)) then
            error = 'invalid n2p2 symmetry-function scaling combination'
            return
        end if
        if ((settings%scale .or. settings%sigma_scale) .and. settings%scale_max <= settings%scale_min) then
            error = 'invalid n2p2 scaling range'
            return
        end if
        if ((settings%has_mean_energy .or. settings%has_conv_energy .or. settings%has_conv_length) .and. &
            .not. (settings%has_mean_energy .and. settings%has_conv_energy .and. settings%has_conv_length)) then
            error = 'n2p2 normalization requires mean_energy, conv_energy, and conv_length'
            return
        end if
        if (settings%conv_energy == 0 .or. settings%conv_length <= 0) then
            error = 'n2p2 conv_energy must be nonzero and conv_length positive'
            return
        end if
        do species = 1,size(settings%species)
            if (.not. allocated(settings%functions(species)%values)) then
                error = 'n2p2 model has no symmetry functions for one or more species'
                return
            end if
            associate(topology => settings%topologies(species))
            n = topology%hidden_layers
            if (.not. topology%hidden_layers_set .or. .not. topology%nodes_set .or. &
                .not. topology%activations_set .or. n < 0 .or. n >= MAX_LAYERS .or. &
                topology%node_count /= n .or. topology%activation_count /= n+1) then
                error = 'incomplete or inconsistent n2p2 network topology'
                return
            end if
            if (any(topology%hidden_nodes(:n) <= 0)) then
                error = 'n2p2 hidden nodes must be positive'
                return
            end if
            do j=1,n+1
                if (activation_code(topology%activations(j)) < 0) then
                    error = 'unsupported n2p2 activation function'
                    return
                end if
            end do
            end associate
            do j=1,size(settings%functions(species)%values)
                sf = settings%functions(species)%values(j)
                call validate_function(sf,error)
                if (len_trim(error) /= 0) return
            end do
        end do
    end subroutine validate_settings

    subroutine validate_function(sf,error)
        type(symmetry_function), intent(in) :: sf
        character(len=*), intent(out) :: error
        real(real64) :: middle
        error = ''
        if (.not. all(ieee_is_finite([sf%cutoff,sf%eta,sf%shift,sf%lambda,sf%zeta,sf%left,sf%right]))) then
            error = 'nonfinite n2p2 symmetry-function parameters'
        elseif (sf%cutoff <= 0) then
            error = 'n2p2 symmetry-function cutoff must be positive'
        elseif (sf%kind == 3 .or. sf%kind == 9 .or. sf%kind == 13) then
            if (abs(sf%lambda) > 1 .or. sf%zeta < 1 .or. sf%zeta > real(huge(1)-1,real64)) &
                error = 'invalid n2p2 angular power'
        elseif (sf%kind >= 20) then
            if (.not. any(sf%subtype == [character(len=3) :: 'e','p1','p2','p3','p4', &
                'p1a','p2a','p3a','p4a'])) then
                error = 'invalid n2p2 compact subtype'
            elseif (sf%shift >= sf%cutoff) then
                error = 'empty n2p2 compact radial interval'
            elseif (sf%kind /= 20 .and. sf%kind /= 23) then
                middle = 0.5_real64*(sf%left+sf%right)
                if (sf%left >= sf%right .or. sf%right-sf%left > 360 .or. &
                    (sf%left < 0 .and. middle /= 0) .or. (sf%right > 180 .and. middle /= 180)) &
                    error = 'invalid n2p2 compact angle interval'
            end if
        end if
    end subroutine validate_function

    subroutine build_network(settings, species, weight_file, scaling_file, network, error)
        character(len=*), intent(out) :: error
        type(n2p2_settings), intent(in) :: settings
        integer, intent(in) :: species
        character(len=*), intent(in) :: weight_file, scaling_file
        type(atomic_network), intent(out) :: network
        integer :: layer, nweights, offset, nsf
        integer(int64) :: next_offset
        type(network_topology) :: topology

        error = ''
        nsf = size(settings%functions(species)%values)
        topology = settings%topologies(species)
        network%atomtype = settings%species(species)
        network%description = "Imported n2p2 2G-HDNNP model"
        network%descriptor_name = "Behler2011"
        if (any(settings%functions(species)%values%kind >= 12)) network%descriptor_name="n2p2_extended"
        network%minimum_radius = 0.1_real64
        network%maximum_radius = maxval(settings%functions(species)%values%cutoff)
        network%nlayers = topology%hidden_layers + 2
        allocate(network%nodes(network%nlayers), network%activation(network%nlayers - 1), &
            network%weight_offsets(network%nlayers))
        network%nodes(1) = nsf
        if (topology%hidden_layers > 0) network%nodes(2:network%nlayers - 1) = topology%hidden_nodes(1:topology%hidden_layers)
        network%nodes(network%nlayers) = 1
        do layer = 1, network%nlayers - 1
            network%activation(layer) = activation_code(topology%activations(layer))
        end do
        network%maxnodes = maxval(network%nodes)
        offset = 0
        do layer = 1, network%nlayers - 1
            network%weight_offsets(layer) = offset
            next_offset = int(offset,int64)+(int(network%nodes(layer),int64)+1)*network%nodes(layer+1)
            if (next_offset > huge(offset)) then
                error = 'n2p2 network connection count exceeds integer range'
                return
            end if
            offset = int(next_offset)
        end do
        network%weight_offsets(network%nlayers) = offset
        nweights = offset
        allocate(network%weights(nweights))
        call read_weights(weight_file, network%weights, error)
        if (len_trim(error) /= 0) return
        if (settings%normalize_nodes) call normalize_network_connections(network)
        call set_descriptor_metadata(settings, species, network)
        call read_scaling(scaling_file, settings, species, network%descriptor_shift, network%descriptor_scale, error)
        if (len_trim(error) /= 0) return
        network%species_names = settings%species
        network%environment_names = settings%species
        network%atomic_references = settings%atomic_references
        network%energy_scale = settings%conv_energy
        network%energy_shift = settings%mean_energy
    end subroutine build_network

    subroutine normalize_network_connections(network)
        type(atomic_network), intent(inout) :: network
        integer :: layer, first, last

        do layer = 1, network%nlayers - 1
            first = network%weight_offsets(layer) + 1
            last = network%weight_offsets(layer + 1)
            network%weights(first:last) = network%weights(first:last)/real(network%nodes(layer), real64)
        end do
    end subroutine normalize_network_connections

    subroutine set_descriptor_metadata(settings, species, network)
        type(n2p2_settings), intent(in) :: settings
        integer, intent(in) :: species
        type(atomic_network), intent(inout) :: network
        integer :: i, n
        type(symmetry_function) :: sf
        n = size(settings%functions(species)%values)
        allocate(network%descriptor_kinds(n), &
            network%descriptor_parameters(merge(9,7,trim(network%descriptor_name) == "n2p2_extended"), n), &
                 network%descriptor_environments(2, n))
        network%descriptor_parameters = 0.0_real64
        if (size(network%descriptor_parameters,1) == 9) then
            network%descriptor_parameters(8,:)=real(settings%cutoff_type,real64)
            network%descriptor_parameters(9,:)=settings%cutoff_alpha
        end if
        network%descriptor_cutoff_type = settings%cutoff_type
        network%descriptor_cutoff_alpha = settings%cutoff_alpha
        network%descriptor_environments = 0
        do i = 1, n
            sf = settings%functions(species)%values(i)
            network%descriptor_environments(1, i) = sf%neighbor1
            network%descriptor_environments(2, i) = sf%neighbor2
            select case(sf%kind)
            case(2)
                network%descriptor_kinds(i) = 2
                network%descriptor_parameters(1:3, i) = [sf%cutoff, sf%shift, sf%eta]
            case(3)
                network%descriptor_kinds(i) = 4
                network%descriptor_parameters(1:4, i) = [sf%cutoff, sf%lambda, sf%zeta, sf%eta]
            case(9)
                network%descriptor_kinds(i) = 5
                network%descriptor_parameters(1:4, i) = [sf%cutoff, sf%lambda, sf%zeta, sf%eta]
            end select
            if (sf%kind >= 12) then
                network%descriptor_kinds(i)=sf%kind
                network%descriptor_parameters(1:7,i)=extended_parameters(sf,settings%cutoff_alpha)
                cycle
            end if
            network%descriptor_parameters(5, i) = real(settings%cutoff_type, real64)
            network%descriptor_parameters(6, i) = settings%cutoff_alpha
            if (sf%kind == 3 .or. sf%kind == 9) network%descriptor_parameters(7, i) = sf%shift
        end do
    end subroutine set_descriptor_metadata

    function extended_parameters(sf,alpha) result(p)
        type(symmetry_function), intent(in) :: sf
        real(real64), intent(in) :: alpha
        real(real64) :: p(7)
        p=0
        if (sf%kind < 20) then
            p=[sf%cutoff,sf%eta,sf%shift,sf%lambda,sf%zeta,0.0_real64,alpha]
        else
            p(1:5)=[sf%cutoff,sf%shift,sf%left,sf%right,real(compact_subtype(sf%subtype),real64)]
        end if
    end function

    subroutine build_setup(settings, species, setup)
        type(n2p2_settings), intent(in) :: settings
        integer, intent(in) :: species
        type(descriptor_setup), intent(out) :: setup
        type(behler_config) :: config
        type(symmetry_function) :: sf
        integer :: i
        setup%central_species = settings%species(species)
        setup%description = "Imported n2p2 2G-HDNNP model"
        setup%minimum_distance = 0.1_real64
        setup%central_global_species = species
        setup%environment_species = settings%species
        allocate(setup%global_to_local(size(settings%species)))
        setup%global_to_local = [(i, i=1, size(settings%species))]
        call initialize_behler_config(config, size(settings%species), &
            settings%cutoff_type, settings%cutoff_alpha)
        config%species_weights=[(real(atomic_number(settings%species(i)),real64),i=1,size(settings%species))]
        do i = 1, size(settings%functions(species)%values)
            sf = settings%functions(species)%values(i)
            select case(sf%kind)
            case(2)
                call add_g2(config, sf%neighbor1, sf%cutoff, sf%shift, sf%eta)
            case(3)
                call add_g4(config, sf%neighbor1, sf%neighbor2, sf%cutoff, sf%lambda, sf%zeta, sf%eta, sf%shift)
            case(9)
                call add_g5(config, sf%neighbor1, sf%neighbor2, sf%cutoff, sf%lambda, sf%zeta, sf%eta, sf%shift)
            case default
                call add_extended(config,sf%kind,sf%neighbor1,sf%neighbor2, &
                    extended_parameters(sf,settings%cutoff_alpha))
            end select
        end do
        call add_behler(setup%model, config)
    end subroutine build_setup

    subroutine read_weights(filename, weights, error)
        character(len=*), intent(out) :: error
        logical :: opened
        character(len=*), intent(in) :: filename
        real(real64), intent(out) :: weights(:)
        character(len=LINE_LENGTH) :: line
        integer :: unit, ios, count
        real(real64) :: value
        error = ''
        opened = .false.
        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            error = "cannot open n2p2 weights file"
            go to 900
        end if
        opened = .true.
        count = 0
        do
            read(unit, "(A)", iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) then
                error = 'I/O error reading n2p2 model file'
                go to 900
            end if
            call clean_line(line)
            if (len_trim(line) == 0) cycle
            read(line, *, iostat=ios) value
            if (ios /= 0) then
                error = "invalid n2p2 weights file"
                go to 900
            end if
            count = count + 1
            if (count > size(weights)) then
                error = "n2p2 weights file has too many connections"
                go to 900
            end if
            if (.not. ieee_is_finite(value)) then
                error = 'nonfinite n2p2 weight'
                go to 900
            end if
            weights(count) = value
        end do
        close(unit)
        opened = .false.
        if (count /= size(weights)) then
            error = "n2p2 weights file connection count does not match topology"
            go to 900
        end if
900     if (opened) close(unit)
    end subroutine read_weights

    subroutine read_scaling(filename, settings, requested_species, shift, scale_values, error)
        character(len=*), intent(out) :: error
        logical :: opened
        character(len=*), intent(in) :: filename
        type(n2p2_settings), intent(in) :: settings
        integer, intent(in) :: requested_species
        real(real64), allocatable, intent(out) :: shift(:), scale_values(:)
        character(len=LINE_LENGTH) :: line
        integer :: unit, ios, species, index, count, n
        real(real64) :: minimum, maximum, mean, sigma, factor
        logical, allocatable :: seen(:)
        error = ''
        opened = .false.
        n = size(settings%functions(requested_species)%values)
        allocate(shift(n), scale_values(n), seen(n))
        shift = 0.0_real64; scale_values = 1.0_real64; seen = .false.
        if (.not. settings%scale .and. .not. settings%center .and. .not. settings%sigma_scale) return
        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            error = "cannot open n2p2 scaling.data"
            go to 900
        end if
        opened = .true.
        count = 0
        do
            read(unit, "(A)", iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) then
                error = 'I/O error reading n2p2 model file'
                go to 900
            end if
            call clean_line(line)
            if (len_trim(line) == 0) cycle
            read(line, *, iostat=ios) species, index, minimum, maximum, mean, sigma
            if (ios /= 0) then
                error = "invalid n2p2 scaling.data"
                go to 900
            end if
            if (species < 1 .or. species > size(settings%species) .or. &
                .not. all(ieee_is_finite([minimum,maximum,mean,sigma]))) then
                error = 'invalid n2p2 scaling data values'
                go to 900
            end if
            if (species /= requested_species) cycle
            if (index < 1 .or. index > n) then
                error = 'invalid n2p2 scaling index'
                go to 900
            end if
            if (seen(index)) then
                error = "invalid n2p2 scaling index"
                go to 900
            end if
            seen(index) = .true.; count = count + 1
            if (settings%sigma_scale) then
                if (sigma <= 0.0_real64) then
                    error = "zero sigma in n2p2 scaling.data"
                    go to 900
                end if
                factor = (settings%scale_max - settings%scale_min)/sigma
                scale_values(index) = factor
                shift(index) = mean - settings%scale_min/factor
            elseif (settings%scale) then
                if (maximum <= minimum) then
                    error = "zero range in n2p2 scaling.data"
                    go to 900
                end if
                factor = (settings%scale_max - settings%scale_min)/(maximum - minimum)
                scale_values(index) = factor
                if (settings%center) then
                    shift(index) = mean - settings%scale_min/factor
                else
                    shift(index) = minimum - settings%scale_min/factor
                end if
            else
                shift(index) = mean
            end if
        end do
        close(unit)
        opened = .false.
        if (.not. all(ieee_is_finite(shift)) .or. .not. all(ieee_is_finite(scale_values))) then
            error = 'nonfinite n2p2 descriptor scaling'
            go to 900
        end if
        if (count /= n .or. .not. all(seen)) then
            error = "n2p2 scaling.data is incomplete"
            go to 900
        end if
900     if (opened) close(unit)
    end subroutine read_scaling

    subroutine append_function(list, value)
        type(symmetry_function_list), intent(inout) :: list
        type(symmetry_function), intent(in) :: value
        if (allocated(list%values)) then
            list%values = [list%values, value]
        else
            list%values = [value]
        end if
    end subroutine append_function

    subroutine sort_functions(values)
        type(symmetry_function), intent(inout) :: values(:)
        type(symmetry_function) :: candidate
        integer :: i, j
        do i = 2, size(values)
            candidate = values(i); j = i - 1
            do while (j >= 1)
                if (.not. function_less(candidate, values(j))) exit
                values(j + 1) = values(j); j = j - 1
            end do
            values(j + 1) = candidate
        end do
    end subroutine sort_functions

    pure logical function function_less(a, b) result(less)
        type(symmetry_function), intent(in) :: a, b
        less = .false.
        if (a%kind /= b%kind) then; less = a%kind < b%kind; return; end if
        if (a%kind >= 20) then
            if (a%subtype /= b%subtype) then; less=a%subtype < b%subtype; return; end if
            if (a%neighbor1 /= b%neighbor1) then; less=a%neighbor1 < b%neighbor1; return; end if
            if (a%neighbor2 /= b%neighbor2) then; less=a%neighbor2 < b%neighbor2; return; end if
            if (a%cutoff /= b%cutoff) then; less=a%cutoff < b%cutoff; return; end if
            if (a%shift /= b%shift) then; less=a%shift < b%shift; return; end if
            if (a%left /= b%left) then; less=a%left < b%left; return; end if
            less=a%right < b%right; return
        end if
        if (a%cutoff /= b%cutoff) then; less = a%cutoff < b%cutoff; return; end if
        if (a%eta /= b%eta) then; less = a%eta < b%eta; return; end if
        if (a%shift /= b%shift) then; less = a%shift < b%shift; return; end if
        if (a%zeta /= b%zeta) then; less = a%zeta < b%zeta; return; end if
        if (a%lambda /= b%lambda) then; less = a%lambda < b%lambda; return; end if
        if (a%neighbor1 /= b%neighbor1) then; less = a%neighbor1 < b%neighbor1; return; end if
        less = a%neighbor2 < b%neighbor2
    end function function_less

    integer function activation_code(name) result(code)
        character(len=1), intent(in) :: name
        select case(name)
        case("l"); code = ACTIVATION_LINEAR
        case("t"); code = ACTIVATION_TANH
        case("s"); code = ACTIVATION_LOGISTIC
        case("p"); code = ACTIVATION_SOFTPLUS
        case("r"); code = ACTIVATION_RELU
        case("g"); code = ACTIVATION_GAUSSIAN
        case("c"); code = ACTIVATION_COSINE
        case("S"); code = ACTIVATION_REVERSE_LOGISTIC
        case("e"); code = ACTIVATION_EXPONENTIAL
        case("h"); code = ACTIVATION_HARMONIC
        case default; code = -1
        end select
    end function activation_code

    integer function species_index(name, species_names) result(found)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: species_names(:)
        integer :: i
        found = 0
        do i = 1, size(species_names)
            if (trim(name) == trim(species_names(i))) then
                found = i; return
            end if
        end do
    end function species_index

    function weight_filename(symbol) result(filename)
        character(len=*), intent(in) :: symbol
        character(len=32) :: filename
        integer :: number
        number = atomic_number(symbol)
        write(filename, '("weights.",I3.3,".data")') number
    end function weight_filename

    integer function atomic_number(symbol) result(number)
        character(len=*), intent(in) :: symbol
        character(len=2), parameter :: elements(118) = [character(len=2) :: &
            "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar", &
            "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", &
            "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe", &
            "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu", &
            "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra", &
            "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db", &
            "Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"]
        integer :: i
        number = 0
        do i = 1, size(elements)
            if (trim(symbol) == trim(elements(i))) then
                number = i; return
            end if
        end do
    end function atomic_number

    subroutine sort_species_by_atomic_number(species)
        character(len=16), intent(inout) :: species(:)
        character(len=16) :: candidate
        integer :: i, j
        do i = 2, size(species)
            candidate = species(i); j = i - 1
            do while (j >= 1)
                if (atomic_number(species(j)) <= atomic_number(candidate)) exit
                species(j + 1) = species(j); j = j - 1
            end do
            species(j + 1) = candidate
        end do
    end subroutine sort_species_by_atomic_number

    function join_path(directory, leaf) result(path)
        character(len=*), intent(in) :: directory, leaf
        character(len=LINE_LENGTH) :: path
        integer :: n
        n = len_trim(directory)
        if (n == 0) then
            path = trim(leaf)
        elseif (directory(n:n) == "/") then
            path = trim(directory)//trim(leaf)
        else
            path = trim(directory)//"/"//trim(leaf)
        end if
    end function join_path

    subroutine clean_line(line)
        character(len=*), intent(inout) :: line
        integer :: comment
        comment = index(line, "#")
        if (comment > 0) line(comment:) = " "
        line = adjustl(line)
    end subroutine clean_line

    subroutine split_key(line, key, rest)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: key, rest
        integer :: separator
        separator = scan(trim(line), " "//achar(9))
        if (separator == 0) then
            key = trim(line); rest = ""
        else
            key = line(:separator - 1); rest = adjustl(line(separator + 1:))
        end if
    end subroutine split_key

    subroutine parse_integer_list(line, values, count, error)
        character(len=*), intent(out) :: error
        character(len=*), intent(in) :: line
        integer, intent(out) :: values(:), count
        character(len=LINE_LENGTH) :: work
        integer :: position, ios, value
        error = ''
        work = adjustl(line); count = 0
        do while (len_trim(work) > 0)
            position = scan(trim(work), " "//achar(9))
            if (position == 0) position = len_trim(work) + 1
            read(work(:position - 1), *, iostat=ios) value
            if (ios /= 0 .or. count == size(values)) then
                error = "invalid n2p2 integer list"
                return
            end if
            count = count + 1; values(count) = value
            if (position > len_trim(work)) exit
            work = adjustl(work(position + 1:))
        end do
    end subroutine parse_integer_list

    subroutine parse_character_list(line, values, count, error)
        character(len=*), intent(out) :: error
        character(len=*), intent(in) :: line
        character(len=1), intent(out) :: values(:)
        integer, intent(out) :: count
        character(len=LINE_LENGTH) :: work
        integer :: position
        error = ''
        work = adjustl(line); count = 0
        do while (len_trim(work) > 0)
            position = scan(trim(work), " "//achar(9))
            if (position == 0) position = len_trim(work) + 1
            if (count == size(values) .or. position /= 2) then
                error = "invalid n2p2 activation list"
                return
            end if
            count = count + 1; values(count) = work(1:1)
            if (position > len_trim(work)) exit
            work = adjustl(work(position + 1:))
        end do
    end subroutine parse_character_list

    pure function lowercase(input) result(output)
        character(len=*), intent(in) :: input
        character(len=len(input)) :: output
        integer :: i, code
        output = input
        do i = 1, len(input)
            code = iachar(input(i:i))
            if (code >= iachar("A") .and. code <= iachar("Z")) output(i:i) = achar(code + 32)
        end do
    end function lowercase

end module n2p2_network
