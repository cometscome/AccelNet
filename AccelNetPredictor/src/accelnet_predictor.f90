module accelnet_predictor
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: atomic_structure, neighbor_data, descriptor_config, &
        initialize_config, read_xsf, build_neighbor_list
    use accelnet_descriptor_models, only: evaluate_model_values, evaluate_model_values_derivatives, &
        contract_model_derivatives, model_supports_direct_contraction, &
        add_chebyshev, add_lj, add_behler
    use accelnet_lj, only: lj_config, initialize_lj_config
    use accelnet_behler, only: behler_config, initialize_behler_config, &
        add_g1, add_g2, add_g3, add_g4, add_g5
    use accelnet_setup, only: descriptor_setup, read_accelnet_setup_set
    use aenet_network, only: atomic_network, read_aenet_network_ascii, read_aenet_network
    use n2p2_network, only: load_n2p2_model
    implicit none
    private

    type, public :: predictor_model
        type(descriptor_setup), allocatable :: setups(:)
        type(atomic_network), allocatable :: networks(:)
        character(len=16), allocatable :: species_names(:)
        real(real64) :: maximum_cutoff = 0.0_real64
        real(real64) :: minimum_distance = 1.0_real64
    contains
        procedure :: reload_from_setups
        procedure :: reload_from_networks
        procedure :: reload_from_n2p2
        generic :: reload => reload_from_setups, reload_from_networks
        procedure :: predict_energy_file
        procedure :: predict_energy_structure
        generic :: predict_energy => predict_energy_file, predict_energy_structure
        procedure :: predict_energy_forces_file
        procedure :: predict_energy_forces_structure
        generic :: predict_energy_forces => predict_energy_forces_file, predict_energy_forces_structure
    end type predictor_model
    public :: load_predictor, load_predictor_from_networks, load_predictor_from_network_data, &
              load_predictor_from_n2p2
    public :: reload_predictor, reload_predictor_from_networks, reload_predictor_from_n2p2

contains
    subroutine load_predictor_from_n2p2(model_directory, model)
        character(len=*), intent(in) :: model_directory
        class(predictor_model), intent(out) :: model
        integer :: species
        call load_n2p2_model(model_directory, model%networks, model%setups, model%species_names)
        do species = 1, size(model%setups)
            model%maximum_cutoff = max(model%maximum_cutoff, model%setups(species)%model%maximum_cutoff)
            model%minimum_distance = min(model%minimum_distance, model%setups(species)%minimum_distance)
        end do
    end subroutine load_predictor_from_n2p2

    subroutine load_predictor(setup_files, network_files, model)
        character(len=*), intent(in) :: setup_files(:), network_files(:)
        class(predictor_model), intent(out) :: model
        integer :: species, descriptor_count
        if (size(setup_files) /= size(network_files)) error stop "one setup and network are required per species"
        call read_accelnet_setup_set(setup_files, model%setups, model%species_names)
        allocate(model%networks(size(network_files)))
        do species = 1, size(network_files)
            call read_aenet_network_ascii(trim(network_files(species)), model%networks(species))
            if (trim(model%networks(species)%atomtype) /= trim(model%species_names(species))) &
                error stop "setup and network species ordering differs"
            descriptor_count = model%setups(species)%model%num_descriptors()
            if (descriptor_count /= model%networks(species)%nodes(1)) error stop "setup/network descriptor mismatch"
            model%maximum_cutoff = max(model%maximum_cutoff, model%setups(species)%model%maximum_cutoff)
            model%minimum_distance = min(model%minimum_distance, model%setups(species)%minimum_distance)
        end do
    end subroutine load_predictor

    subroutine load_predictor_from_networks(network_files, model, chebyshev_version)
        character(len=*), intent(in) :: network_files(:)
        class(predictor_model), intent(out) :: model
        integer, intent(in), optional :: chebyshev_version
        integer :: species, global, local, version
        version = 0
        if (present(chebyshev_version)) version = chebyshev_version
        if (version /= 0 .and. version /= 1 .and. version /= 10) &
            error stop "Chebyshev version must be 0, 1, or 10"
        allocate(model%networks(size(network_files)), model%setups(size(network_files)))
        do species = 1, size(network_files)
            call read_aenet_network(trim(network_files(species)), model%networks(species))
        end do
        allocate(model%species_names(size(model%networks(1)%species_names)))
        model%species_names = model%networks(1)%species_names
        if (size(model%species_names) /= size(network_files)) error stop "one network is required per species"
        do species = 1, size(network_files)
            if (trim(model%networks(species)%atomtype) /= trim(model%species_names(species))) &
                error stop "network ordering must follow embedded global species ordering"
            call setup_from_network(model%networks(species), species, version, model%setups(species))
            allocate(model%setups(species)%global_to_local(size(model%species_names)))
            do global = 1, size(model%species_names)
                model%setups(species)%global_to_local(global) = 0
                do local = 1, size(model%networks(species)%environment_names)
                    if (trim(model%species_names(global)) == &
                        trim(model%networks(species)%environment_names(local))) &
                        model%setups(species)%global_to_local(global) = local
                end do
                if (model%setups(species)%global_to_local(global) == 0) &
                    error stop "embedded environment species are incomplete"
            end do
            model%maximum_cutoff = max(model%maximum_cutoff, model%setups(species)%model%maximum_cutoff)
            model%minimum_distance = min(model%minimum_distance, model%setups(species)%minimum_distance)
        end do
    end subroutine load_predictor_from_networks

    subroutine load_predictor_from_network_data(networks, model, chebyshev_version)
        type(atomic_network), intent(in) :: networks(:)
        class(predictor_model), intent(out) :: model
        integer, intent(in), optional :: chebyshev_version
        integer :: species, global, local, version
        version = 0
        if (present(chebyshev_version)) version = chebyshev_version
        if (version /= 0 .and. version /= 1 .and. version /= 10) &
            error stop "Chebyshev version must be 0, 1, or 10"
        if (size(networks) < 1) error stop "at least one network is required"
        allocate(model%networks(size(networks)), model%setups(size(networks)))
        model%networks = networks
        allocate(model%species_names(size(model%networks(1)%species_names)))
        model%species_names = model%networks(1)%species_names
        if (size(model%species_names) /= size(networks)) error stop "one network is required per species"
        do species = 1, size(networks)
            if (trim(model%networks(species)%atomtype) /= trim(model%species_names(species))) &
                error stop "network ordering must follow embedded global species ordering"
            call setup_from_network(model%networks(species), species, version, model%setups(species))
            allocate(model%setups(species)%global_to_local(size(model%species_names)))
            do global = 1, size(model%species_names)
                model%setups(species)%global_to_local(global) = 0
                do local = 1, size(model%networks(species)%environment_names)
                    if (trim(model%species_names(global)) == &
                        trim(model%networks(species)%environment_names(local))) &
                        model%setups(species)%global_to_local(global) = local
                end do
                if (model%setups(species)%global_to_local(global) == 0) &
                    error stop "embedded environment species are incomplete"
            end do
            model%maximum_cutoff = max(model%maximum_cutoff, model%setups(species)%model%maximum_cutoff)
            model%minimum_distance = min(model%minimum_distance, model%setups(species)%minimum_distance)
        end do
    end subroutine load_predictor_from_network_data

    subroutine reload_predictor(model, setup_files, network_files)
        class(predictor_model), intent(inout) :: model
        character(len=*), intent(in) :: setup_files(:), network_files(:)
        call load_predictor(setup_files, network_files, model)
    end subroutine reload_predictor

    subroutine reload_predictor_from_networks(model, network_files, chebyshev_version)
        class(predictor_model), intent(inout) :: model
        character(len=*), intent(in) :: network_files(:)
        integer, intent(in), optional :: chebyshev_version
        if (present(chebyshev_version)) then
            call load_predictor_from_networks(network_files, model, chebyshev_version)
        else
            call load_predictor_from_networks(network_files, model)
        end if
    end subroutine reload_predictor_from_networks

    subroutine reload_predictor_from_n2p2(model, model_directory)
        class(predictor_model), intent(inout) :: model
        character(len=*), intent(in) :: model_directory
        call load_predictor_from_n2p2(model_directory, model)
    end subroutine reload_predictor_from_n2p2

    subroutine reload_from_setups(self, setup_files, network_files)
        class(predictor_model), intent(inout) :: self
        character(len=*), intent(in) :: setup_files(:), network_files(:)
        call reload_predictor(self, setup_files, network_files)
    end subroutine reload_from_setups

    subroutine reload_from_networks(self, network_files, chebyshev_version)
        class(predictor_model), intent(inout) :: self
        character(len=*), intent(in) :: network_files(:)
        integer, intent(in), optional :: chebyshev_version
        if (present(chebyshev_version)) then
            call reload_predictor_from_networks(self, network_files, chebyshev_version)
        else
            call reload_predictor_from_networks(self, network_files)
        end if
    end subroutine reload_from_networks

    subroutine reload_from_n2p2(self, model_directory)
        class(predictor_model), intent(inout) :: self
        character(len=*), intent(in) :: model_directory
        call reload_predictor_from_n2p2(self, model_directory)
    end subroutine reload_from_n2p2

    subroutine setup_from_network(network, central_species, chebyshev_version, setup)
        type(atomic_network), intent(in) :: network
        integer, intent(in) :: central_species, chebyshev_version
        type(descriptor_setup), intent(out) :: setup
        type(descriptor_config) :: chebyshev
        type(lj_config) :: lj
        type(behler_config) :: behler
        character(len=100) :: name
        setup%central_species = network%atomtype
        setup%description = network%description
        setup%minimum_distance = network%minimum_radius
        setup%central_global_species = central_species
        allocate(setup%environment_species(size(network%environment_names)))
        setup%environment_species = network%environment_names
        name = lowercase(trim(network%descriptor_name))
        select case(trim(name))
        case("chebyshev")
            call initialize_config(chebyshev, size(network%environment_names), &
                network%descriptor_parameters(1, 1), nint(network%descriptor_parameters(2, 1)), &
                network%descriptor_parameters(3, 1), nint(network%descriptor_parameters(4, 1)), &
                version=chebyshev_version, central_type_index=central_species)
            call add_chebyshev(setup%model, chebyshev)
        case("lj")
            call initialize_lj_config(lj, size(network%environment_names), &
                network%descriptor_parameters(1, 1))
            call add_lj(setup%model, lj)
        case("behler2011")
            call initialize_behler_config(behler, size(network%environment_names))
            call add_aenet_ordered_behler_functions(behler, network)
            call add_behler(setup%model, behler)
        case default
            error stop "unsupported embedded descriptor type"
        end select
        if (setup%model%num_descriptors() /= network%nodes(1)) &
            error stop "embedded descriptor metadata has wrong dimension"
    end subroutine setup_from_network

    subroutine add_aenet_ordered_behler_functions(config, network)
        type(behler_config), intent(inout) :: config
        type(atomic_network), intent(in) :: network
        integer :: descriptor, kind, first_species, second_species, added
        integer :: low_species, high_species

        ! aenet's symmfunc module stores Behler functions by environment
        ! species and function kind, independently of their order in the
        ! setup/network metadata.  The NN inputs and scaling arrays follow
        ! this internal order, so reproduce it when rebuilding a setup from
        ! an embedded native-network descriptor block.
        do descriptor = 1, size(network%descriptor_kinds)
            kind = network%descriptor_kinds(descriptor)
            first_species = network%descriptor_environments(1, descriptor)
            second_species = network%descriptor_environments(2, descriptor)
            if (kind < 1 .or. kind > 5) error stop "unsupported embedded Behler descriptor kind"
            if (first_species < 1 .or. first_species > config%num_species) &
                error stop "embedded Behler descriptor species is out of range"
            if (kind >= 4 .and. (second_species < 1 .or. second_species > config%num_species)) &
                error stop "embedded Behler angular descriptor species is out of range"
        end do
        added = 0
        do first_species = 1, config%num_species
            do kind = 1, 3
                do descriptor = 1, size(network%descriptor_kinds)
                    if (network%descriptor_kinds(descriptor) /= kind .or. &
                        network%descriptor_environments(1, descriptor) /= first_species) cycle
                    select case(kind)
                    case(1)
                        call add_g1(config, first_species, network%descriptor_parameters(1, descriptor))
                    case(2)
                        call add_g2(config, first_species, network%descriptor_parameters(1, descriptor), &
                            network%descriptor_parameters(2, descriptor), &
                            network%descriptor_parameters(3, descriptor))
                    case(3)
                        call add_g3(config, first_species, network%descriptor_parameters(1, descriptor), &
                            network%descriptor_parameters(2, descriptor))
                    end select
                    added = added + 1
                end do
            end do
            do second_species = first_species, config%num_species
                do kind = 4, 5
                    do descriptor = 1, size(network%descriptor_kinds)
                        low_species = min(network%descriptor_environments(1, descriptor), &
                                          network%descriptor_environments(2, descriptor))
                        high_species = max(network%descriptor_environments(1, descriptor), &
                                           network%descriptor_environments(2, descriptor))
                        if (network%descriptor_kinds(descriptor) /= kind .or. &
                            low_species /= first_species .or. high_species /= second_species) cycle
                        if (kind == 4) then
                            call add_g4(config, first_species, second_species, &
                                network%descriptor_parameters(1, descriptor), &
                                network%descriptor_parameters(2, descriptor), &
                                network%descriptor_parameters(3, descriptor), &
                                network%descriptor_parameters(4, descriptor))
                        else
                            call add_g5(config, first_species, second_species, &
                                network%descriptor_parameters(1, descriptor), &
                                network%descriptor_parameters(2, descriptor), &
                                network%descriptor_parameters(3, descriptor), &
                                network%descriptor_parameters(4, descriptor))
                        end if
                        added = added + 1
                    end do
                end do
            end do
        end do
        if (added /= size(network%descriptor_kinds)) &
            error stop "embedded Behler descriptor ordering is incomplete"
    end subroutine add_aenet_ordered_behler_functions

    pure function lowercase(input) result(output)
        character(len=*), intent(in) :: input
        character(len=len(input)) :: output
        integer :: i, code
        output = input
        do i = 1, len(input)
            code = iachar(input(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) output(i:i) = achar(code + 32)
        end do
    end function lowercase

    subroutine predict_energy_file(self, xsf_file, total_energy)
        class(predictor_model), intent(in) :: self
        character(len=*), intent(in) :: xsf_file
        real(real64), intent(out) :: total_energy
        type(atomic_structure) :: structure
        call read_xsf(trim(xsf_file), self%species_names, structure)
        call self%predict_energy(structure, total_energy)
    end subroutine predict_energy_file

    subroutine predict_energy_structure(self, structure, total_energy)
        class(predictor_model), intent(in) :: self
        type(atomic_structure), intent(in) :: structure
        real(real64), intent(out) :: total_energy
        type(neighbor_data) :: neighbors
        real(real64), allocatable :: descriptor(:), normalized(:)
        integer, allocatable :: global_neighbors(:), local_neighbors(:)
        integer :: atom, species, first, last, n, maximum_dimension, maximum_neighbors
        real(real64) :: atomic_energy, cohesive
        call build_neighbor_list(structure, self%maximum_cutoff, neighbors, self%minimum_distance)
        maximum_dimension = 0
        do species = 1, size(self%networks)
            maximum_dimension = max(maximum_dimension, self%networks(species)%nodes(1))
        end do
        maximum_neighbors = maxval(neighbors%offsets(2:) - neighbors%offsets(:structure%natoms))
        allocate(descriptor(maximum_dimension), normalized(maximum_dimension), &
                 global_neighbors(maximum_neighbors), local_neighbors(maximum_neighbors))
        cohesive = 0.0_real64
        do atom = 1, structure%natoms
            species = structure%species(atom)
            first = neighbors%offsets(atom); last = neighbors%offsets(atom + 1) - 1; n = max(0, last - first + 1)
            if (n > 0) global_neighbors(1:n) = structure%species(neighbors%atom_indices(first:last))
            call self%setups(species)%map_species(global_neighbors(1:n), local_neighbors(1:n))
            call evaluate_model_values(self%setups(species)%model, neighbors%displacements(:, first:last), &
                local_neighbors(1:n), descriptor)
            normalized(1:self%networks(species)%nodes(1)) = &
                (descriptor(1:self%networks(species)%nodes(1)) - self%networks(species)%descriptor_shift) * &
                self%networks(species)%descriptor_scale
            call self%networks(species)%evaluate(normalized(1:self%networks(species)%nodes(1)), atomic_energy)
            cohesive = cohesive + atomic_energy
        end do
        total_energy = cohesive/self%networks(1)%energy_scale + structure%natoms*self%networks(1)%energy_shift
        do atom = 1, structure%natoms
            total_energy = total_energy + self%networks(1)%atomic_references(structure%species(atom))
        end do
    end subroutine predict_energy_structure

    subroutine predict_energy_forces_file(self, xsf_file, total_energy, forces)
        class(predictor_model), intent(in) :: self
        character(len=*), intent(in) :: xsf_file
        real(real64), intent(out) :: total_energy
        real(real64), allocatable, intent(out) :: forces(:, :)
        type(atomic_structure) :: structure
        call read_xsf(trim(xsf_file), self%species_names, structure)
        allocate(forces(3, structure%natoms))
        call self%predict_energy_forces(structure, total_energy, forces)
    end subroutine predict_energy_forces_file

    subroutine predict_energy_forces_structure(self, structure, total_energy, forces)
        class(predictor_model), intent(in) :: self
        type(atomic_structure), intent(in) :: structure
        real(real64), intent(out) :: total_energy
        real(real64), intent(out) :: forces(:, :)
        type(neighbor_data) :: neighbors
        real(real64), allocatable :: descriptor(:), normalized(:), gradient(:), contributions(:)
        real(real64), allocatable :: derivative_center(:, :), derivative_neighbors(:, :, :)
        real(real64), allocatable :: contracted_neighbors(:, :)
        real(real64) :: contracted_center(3)
        integer, allocatable :: global_neighbors(:), local_neighbors(:)
        integer :: atom, species, first, last, n, maximum_dimension, maximum_neighbors
        integer :: neighbor, target, coefficient
        real(real64) :: atomic_energy, cohesive
        logical :: use_direct_contraction
        if (size(forces, 1) /= 3 .or. size(forces, 2) /= structure%natoms) &
            error stop "forces must have shape (3, structure%natoms)"
        call build_neighbor_list(structure, self%maximum_cutoff, neighbors, self%minimum_distance)
        maximum_dimension = 0
        do species = 1, size(self%networks)
            maximum_dimension = max(maximum_dimension, self%networks(species)%nodes(1))
        end do
        use_direct_contraction = .true.
        do species = 1, size(self%setups)
            use_direct_contraction = use_direct_contraction .and. &
                model_supports_direct_contraction(self%setups(species)%model)
        end do
        maximum_neighbors = maxval(neighbors%offsets(2:) - neighbors%offsets(:structure%natoms))
        allocate(descriptor(maximum_dimension), normalized(maximum_dimension), gradient(maximum_dimension), &
                 contributions(maximum_dimension))
        if (use_direct_contraction) then
            allocate(contracted_neighbors(3, maximum_neighbors))
        else
            allocate(derivative_center(3, maximum_dimension), &
                     derivative_neighbors(3, maximum_dimension, maximum_neighbors))
        end if
        allocate(global_neighbors(maximum_neighbors), local_neighbors(maximum_neighbors))
        forces = 0.0_real64; cohesive = 0.0_real64
        do atom = 1, structure%natoms
            species = structure%species(atom)
            first = neighbors%offsets(atom); last = neighbors%offsets(atom + 1) - 1; n = max(0, last - first + 1)
            if (n > 0) global_neighbors(1:n) = structure%species(neighbors%atom_indices(first:last))
            call self%setups(species)%map_species(global_neighbors(1:n), local_neighbors(1:n))
            if (use_direct_contraction) then
                call evaluate_model_values(self%setups(species)%model, &
                    neighbors%displacements(:, first:last), local_neighbors(1:n), descriptor)
            else
                call evaluate_model_values_derivatives(self%setups(species)%model, &
                    neighbors%displacements(:, first:last), local_neighbors(1:n), descriptor, &
                    derivative_center, derivative_neighbors(:, :, 1:n))
            end if
            normalized(1:self%networks(species)%nodes(1)) = &
                (descriptor(1:self%networks(species)%nodes(1)) - self%networks(species)%descriptor_shift) * &
                self%networks(species)%descriptor_scale
            call self%networks(species)%input_gradient(normalized(1:self%networks(species)%nodes(1)), &
                atomic_energy, gradient(1:self%networks(species)%nodes(1)))
            cohesive = cohesive + atomic_energy
            contributions(1:self%networks(species)%nodes(1)) = &
                -gradient(1:self%networks(species)%nodes(1))*self%networks(species)%descriptor_scale / &
                self%networks(species)%energy_scale
            if (use_direct_contraction) then
                call contract_model_derivatives(self%setups(species)%model, &
                    neighbors%displacements(:, first:last), local_neighbors(1:n), &
                    contributions(1:self%networks(species)%nodes(1)), contracted_center, &
                    contracted_neighbors(:, 1:n))
                forces(:, atom) = forces(:, atom) + contracted_center
                do neighbor = 1, n
                    target = neighbors%atom_indices(first + neighbor - 1)
                    forces(:, target) = forces(:, target) + contracted_neighbors(:, neighbor)
                end do
            else
                do coefficient = 1, self%networks(species)%nodes(1)
                    forces(:, atom) = forces(:, atom) + &
                        contributions(coefficient)*derivative_center(:, coefficient)
                end do
                do neighbor = 1, n
                    target = neighbors%atom_indices(first + neighbor - 1)
                    do coefficient = 1, self%networks(species)%nodes(1)
                        forces(:, target) = forces(:, target) + &
                            contributions(coefficient)*derivative_neighbors(:, coefficient, neighbor)
                    end do
                end do
            end if
        end do
        total_energy = cohesive/self%networks(1)%energy_scale + structure%natoms*self%networks(1)%energy_shift
        do atom = 1, structure%natoms
            total_energy = total_energy + self%networks(1)%atomic_references(structure%species(atom))
        end do
    end subroutine predict_energy_forces_structure
end module accelnet_predictor
