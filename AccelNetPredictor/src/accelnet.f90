module accelnet
    use iso_c_binding, only: c_bool, c_char, c_double, c_f_pointer, c_int, c_null_char, c_ptr
    use iso_fortran_env, only: real64
    use accelnet_batch_target_serial, only: target_model, target_workspace, evaluate_batch_target, &
        evaluate_atomic_energy_target
    use accelnet_descriptors, only: atomic_structure, neighbor_data, descriptor_config, &
        build_neighbor_list, initialize_config, chebyshev_values, &
        CHEBYSHEV_EVALUATION_AUTO, CHEBYSHEV_EVALUATION_DIRECT, CHEBYSHEV_EVALUATION_MOMENT
    use accelnet_descriptor_models, only: add_chebyshev
    use accelnet_legacy_lcl, only: lcl_nmax_nbdist
    use accelnet_predictor, only: predictor_model, load_predictor_from_network_data, &
        load_predictor_from_n2p2
    use aenet_network, only: atomic_network, read_aenet_network, read_aenet_network_ascii
    use accelnet_behler, only: G5_EVALUATION_AUTO, G5_EVALUATION_DIRECT, &
        G5_EVALUATION_MOMENT, G5_EVALUATION_MOMENT_FORCE
    implicit none
    private
    save

    integer, parameter :: TYPE_LENGTH = 16
    integer, parameter :: PATH_LENGTH = 1024
    real(real64), parameter :: PI_COMPAT = 3.1415926535897932384626433832795_real64

    integer(c_int), bind(C, name="ACCELNET_OK"), public :: ACCELNET_OK = 0_c_int
    integer(c_int), bind(C, name="ACCELNET_ERR_INIT"), public :: ACCELNET_ERR_INIT = 1_c_int
    integer(c_int), bind(C, name="ACCELNET_ERR_MALLOC"), public :: ACCELNET_ERR_MALLOC = 2_c_int
    integer(c_int), bind(C, name="ACCELNET_ERR_IO"), public :: ACCELNET_ERR_IO = 3_c_int
    integer(c_int), bind(C, name="ACCELNET_ERR_TYPE"), public :: ACCELNET_ERR_TYPE = 4_c_int
    integer(c_int), bind(C, name="ACCELNET_ERR_ARGUMENT"), public :: ACCELNET_ERR_ARGUMENT = 5_c_int
    integer(c_int), bind(C, name="ACCELNET_TYPELEN"), public :: ACCELNET_TYPELEN = TYPE_LENGTH
    integer(c_int), bind(C, name="ACCELNET_PATHLEN"), public :: ACCELNET_PATHLEN = PATH_LENGTH
    logical(c_bool), bind(C, name="ACCELNET_TRUE"), public :: ACCELNET_TRUE = .true._c_bool
    logical(c_bool), bind(C, name="ACCELNET_FALSE"), public :: ACCELNET_FALSE = .false._c_bool
    integer(c_int), bind(C, name="ACCELNET_G5_AUTO"), public :: ACCELNET_G5_AUTO = G5_EVALUATION_AUTO
    integer(c_int), bind(C, name="ACCELNET_G5_DIRECT"), public :: ACCELNET_G5_DIRECT = G5_EVALUATION_DIRECT
    integer(c_int), bind(C, name="ACCELNET_G5_MOMENT"), public :: ACCELNET_G5_MOMENT = G5_EVALUATION_MOMENT
    integer(c_int), bind(C, name="ACCELNET_G5_MOMENT_FORCE"), public :: ACCELNET_G5_MOMENT_FORCE = &
        G5_EVALUATION_MOMENT_FORCE
    integer(c_int), bind(C, name="ACCELNET_CHEBYSHEV_AUTO"), public :: ACCELNET_CHEBYSHEV_AUTO = &
        CHEBYSHEV_EVALUATION_AUTO
    integer(c_int), bind(C, name="ACCELNET_CHEBYSHEV_DIRECT"), public :: ACCELNET_CHEBYSHEV_DIRECT = &
        CHEBYSHEV_EVALUATION_DIRECT
    integer(c_int), bind(C, name="ACCELNET_CHEBYSHEV_MOMENT"), public :: ACCELNET_CHEBYSHEV_MOMENT = &
        CHEBYSHEV_EVALUATION_MOMENT

    integer(c_int), bind(C, name="accelnet_nsf_max"), public :: accelnet_nsf_max = 0_c_int
    integer(c_int), bind(C, name="accelnet_nnb_max"), public :: accelnet_nnb_max = 0_c_int
    real(c_double), bind(C, name="accelnet_Rc_min"), public :: accelnet_Rc_min = 0.0_c_double
    real(c_double), bind(C, name="accelnet_Rc_max"), public :: accelnet_Rc_max = 0.0_c_double

    logical :: is_initialized = .false.
    logical :: is_loaded = .false.
    integer :: number_of_types = 0
    integer :: chebyshev_version = 0
    integer :: chebyshev_evaluation_mode = CHEBYSHEV_EVALUATION_AUTO
    character(len=TYPE_LENGTH), allocatable :: atom_types(:)
    type(atomic_network), allocatable :: pending_networks(:)
    logical, allocatable :: potential_loaded(:)
    type(predictor_model), allocatable :: global_model
    type(target_model) :: global_batch_model
    type(target_workspace) :: global_shared_work
    logical :: batch_cache_valid = .false.
    integer :: batch_packed_status = 1
    integer, allocatable :: atomic_species(:),atomic_indices(:)
    real(real64), allocatable :: atomic_displacements(:,:),atomic_forces(:,:)

    logical :: neighbor_list_initialized = .false.
    type(atomic_structure) :: neighbor_structure
    type(neighbor_data) :: global_neighbors

    logical :: sfb_initialized = .false.
    type(descriptor_config) :: sfb_config
    type(target_model) :: sfb_model
    type(target_workspace) :: sfb_work
    real(real64) :: sfb_radial_cutoff = 0.0_real64

    public :: accelnet_init, accelnet_init_n2p2, accelnet_final, accelnet_all_loaded
    public :: accelnet_load_potential, accelnet_load_n2p2, accelnet_print_info
    public :: accelnet_set_chebyshev_version
    public :: accelnet_set_chebyshev_evaluation, accelnet_get_chebyshev_evaluation
    public :: accelnet_set_g5_evaluation
    public :: accelnet_atomic_energy, accelnet_atomic_energy_and_forces
    public :: accelnet_atomic_energy_and_forces_virial
    public :: accelnet_batch_energy_and_forces
    public :: accelnet_convert_atom_types, accelnet_free_atom_energy
    public :: accelnet_nbl_init, accelnet_nbl_final, accelnet_nbl_neighbors
    public :: accelnet_sfb_init, accelnet_sfb_final, accelnet_sfb_nvalues
    public :: accelnet_sfb_eval, accelnet_sfb_reconstruct_radial

    ! The C wrappers have distinct Fortran names. Explicit visibility also
    ! avoids hidden C binding labels with GNU Fortran 16.2 (GCC PR126872).
    public :: accelnet_init_c, accelnet_init_n2p2_c
    public :: accelnet_load_potential_c, accelnet_load_potential_ascii_c, accelnet_load_n2p2_c
    public :: accelnet_convert_atom_types_c, accelnet_sfb_init_c

contains

    subroutine accelnet_init(species, stat)
        character(len=*), intent(in) :: species(:)
        integer, intent(out) :: stat
        integer :: allocation_status
        stat = ACCELNET_OK
        if (is_initialized .or. size(species) < 1) then
            stat = ACCELNET_ERR_INIT
            return
        end if
        number_of_types = size(species)
        allocate(atom_types(number_of_types), pending_networks(number_of_types), &
                 potential_loaded(number_of_types), stat=allocation_status)
        if (allocation_status /= 0) then
            stat = ACCELNET_ERR_MALLOC
            number_of_types = 0
            return
        end if
        atom_types = species
        potential_loaded = .false.
        chebyshev_version = 0
        chebyshev_evaluation_mode = CHEBYSHEV_EVALUATION_AUTO
        is_initialized = .true.
        is_loaded = .false.
        call reset_public_ranges()
    end subroutine accelnet_init

    subroutine accelnet_init_c(ntypes, species_pointers, stat) bind(C, name="accelnet_init")
        integer(c_int), value, intent(in) :: ntypes
        type(c_ptr), intent(in) :: species_pointers(ntypes)
        integer(c_int), intent(out) :: stat
        character(len=TYPE_LENGTH) :: species(ntypes)
        integer :: i, local_status
        do i = 1, ntypes
            call copy_c_pointer_string(species_pointers(i), species(i))
        end do
        call accelnet_init(species, local_status)
        stat = local_status
    end subroutine accelnet_init_c

    subroutine accelnet_init_n2p2(directory, stat)
        character(len=*), intent(in) :: directory
        integer, intent(out) :: stat
        type(predictor_model) :: loaded_model
        character(len=PATH_LENGTH) :: input_file
        logical :: exists
        integer :: species

        stat = ACCELNET_OK
        if (is_initialized) then
            stat = ACCELNET_ERR_INIT
            return
        end if
        input_file = trim(directory)//"/input.nn"
        inquire(file=trim(input_file), exist=exists)
        if (.not. exists) then
            stat = ACCELNET_ERR_IO
            return
        end if
        call load_predictor_from_n2p2(trim(directory), loaded_model)
        call accelnet_init(loaded_model%species_names, stat)
        if (stat /= ACCELNET_OK) return
        allocate(global_model)
        global_model = loaded_model
        pending_networks = global_model%networks
        potential_loaded = .true.
        is_loaded = .true.
        batch_cache_valid = .false.
        accelnet_Rc_min = global_model%minimum_distance
        accelnet_Rc_max = global_model%maximum_cutoff
        accelnet_nsf_max = maxval([(global_model%networks(species)%nodes(1), &
                                  species=1,number_of_types)])
        accelnet_nnb_max = lcl_nmax_nbdist(real(accelnet_Rc_min, real64), &
                                           real(accelnet_Rc_max, real64))
    end subroutine accelnet_init_n2p2

    subroutine accelnet_init_n2p2_c(directory, stat) bind(C, name="accelnet_init_n2p2")
        character(c_char), intent(in) :: directory(*)
        integer(c_int), intent(out) :: stat
        character(len=PATH_LENGTH) :: local_directory
        integer :: local_status
        call copy_c_string(directory, local_directory)
        call accelnet_init_n2p2(trim(local_directory), local_status)
        stat = local_status
    end subroutine accelnet_init_n2p2_c

    subroutine accelnet_final(stat) bind(C)
        integer(c_int), intent(out) :: stat
        stat = ACCELNET_OK
        if (neighbor_list_initialized) call accelnet_nbl_final()
        call global_shared_work%release()
        call global_batch_model%release()
        if (allocated(atomic_species)) deallocate(atomic_species,atomic_indices,atomic_displacements,atomic_forces)
        batch_cache_valid = .false.
        if (allocated(global_model)) deallocate(global_model)
        if (allocated(pending_networks)) deallocate(pending_networks)
        if (allocated(potential_loaded)) deallocate(potential_loaded)
        if (allocated(atom_types)) deallocate(atom_types)
        number_of_types = 0
        is_initialized = .false.
        is_loaded = .false.
        chebyshev_evaluation_mode = CHEBYSHEV_EVALUATION_AUTO
        call reset_public_ranges()
    end subroutine accelnet_final

    subroutine reset_public_ranges()
        accelnet_nsf_max = 0
        accelnet_nnb_max = 0
        accelnet_Rc_min = 0.0_c_double
        accelnet_Rc_max = 0.0_c_double
    end subroutine reset_public_ranges

    subroutine accelnet_set_chebyshev_version(version, stat) bind(C)
        integer(c_int), value, intent(in) :: version
        integer(c_int), intent(out) :: stat
        stat = ACCELNET_OK
        if (.not. is_initialized .or. is_loaded) then
            stat = ACCELNET_ERR_INIT
        else if (version /= 0 .and. version /= 1 .and. version /= 10) then
            stat = ACCELNET_ERR_ARGUMENT
        else
            chebyshev_version = version
        end if
    end subroutine accelnet_set_chebyshev_version

    subroutine accelnet_set_chebyshev_evaluation(mode, stat) bind(C)
        integer(c_int), value, intent(in) :: mode
        integer(c_int), intent(out) :: stat
        stat = ACCELNET_OK
        if (.not. is_initialized) then
            stat = ACCELNET_ERR_INIT
        else if (mode < CHEBYSHEV_EVALUATION_AUTO .or. mode > CHEBYSHEV_EVALUATION_MOMENT) then
            stat = ACCELNET_ERR_ARGUMENT
        else
            chebyshev_evaluation_mode = int(mode)
            batch_cache_valid = .false.
            if (is_loaded .and. allocated(global_model)) &
                call global_model%set_chebyshev_evaluation(chebyshev_evaluation_mode)
        end if
    end subroutine accelnet_set_chebyshev_evaluation

    integer(c_int) function accelnet_get_chebyshev_evaluation() bind(C)
        accelnet_get_chebyshev_evaluation = int(chebyshev_evaluation_mode, c_int)
    end function accelnet_get_chebyshev_evaluation

    subroutine accelnet_set_g5_evaluation(mode, stat) bind(C)
        integer(c_int), value, intent(in) :: mode
        integer(c_int), intent(out) :: stat
        stat = ACCELNET_OK
        if (.not. is_loaded .or. .not. allocated(global_model)) then
            stat = ACCELNET_ERR_INIT
        else if (mode < G5_EVALUATION_AUTO .or. mode > G5_EVALUATION_MOMENT_FORCE) then
            stat = ACCELNET_ERR_ARGUMENT
        else
            call global_model%set_g5_evaluation(int(mode))
            batch_cache_valid = .false.
        end if
    end subroutine accelnet_set_g5_evaluation

    subroutine accelnet_load_potential(type_id, filename, stat, is_ascii)
        integer, intent(in) :: type_id
        character(len=*), intent(in) :: filename
        integer, intent(out) :: stat
        logical, intent(in), optional :: is_ascii
        logical :: exists, read_ascii
        integer :: species, global_type
        stat = ACCELNET_OK
        if (.not. is_initialized) then
            stat = ACCELNET_ERR_INIT
            return
        end if
        if (type_id < 1 .or. type_id > number_of_types) then
            stat = ACCELNET_ERR_TYPE
            return
        end if
        inquire(file=trim(filename), exist=exists)
        if (.not. exists) then
            stat = ACCELNET_ERR_IO
            return
        end if
        read_ascii = .false.
        if (present(is_ascii)) read_ascii = is_ascii
        if (read_ascii) then
            call read_aenet_network_ascii(trim(filename), pending_networks(type_id))
        else
            call read_aenet_network(trim(filename), pending_networks(type_id))
        end if
        if (trim(pending_networks(type_id)%atomtype) /= trim(atom_types(type_id))) then
            stat = ACCELNET_ERR_TYPE
            return
        end if
        potential_loaded(type_id) = .true.
        if (all(potential_loaded)) then
            do species = 1, number_of_types
                if (size(pending_networks(species)%species_names) /= number_of_types .or. &
                    any([(trim(pending_networks(species)%species_names(global_type)) /= trim(atom_types(global_type)), &
                          global_type=1,number_of_types)])) then
                    stat = ACCELNET_ERR_TYPE
                    return
                end if
            end do
            if (allocated(global_model)) deallocate(global_model)
            allocate(global_model)
            call load_predictor_from_network_data(pending_networks, global_model, chebyshev_version)
            call global_model%set_chebyshev_evaluation(chebyshev_evaluation_mode)
            is_loaded = .true.
            batch_cache_valid = .false.
            accelnet_Rc_min = global_model%minimum_distance
            accelnet_Rc_max = global_model%maximum_cutoff
            accelnet_nsf_max = maxval([(global_model%networks(species)%nodes(1), &
                                      species=1,number_of_types)])
            accelnet_nnb_max = lcl_nmax_nbdist(real(accelnet_Rc_min, real64), &
                                               real(accelnet_Rc_max, real64))
        end if
    end subroutine accelnet_load_potential

    subroutine accelnet_load_potential_c(type_id, filename, stat) bind(C, name="accelnet_load_potential")
        integer(c_int), value, intent(in) :: type_id
        character(c_char), intent(in) :: filename(*)
        integer(c_int), intent(out) :: stat
        character(len=PATH_LENGTH) :: local_filename
        integer :: local_status
        call copy_c_string(filename, local_filename)
        call accelnet_load_potential(type_id, trim(local_filename), local_status)
        stat = local_status
    end subroutine accelnet_load_potential_c

    subroutine accelnet_load_n2p2(directory, stat)
        character(len=*), intent(in) :: directory
        integer, intent(out) :: stat
        character(len=PATH_LENGTH) :: input_file
        logical :: exists
        integer :: species, global_type

        stat = ACCELNET_OK
        if (.not. is_initialized .or. is_loaded .or. any(potential_loaded)) then
            stat = ACCELNET_ERR_INIT
            return
        end if
        input_file = trim(directory)//"/input.nn"
        inquire(file=trim(input_file), exist=exists)
        if (.not. exists) then
            stat = ACCELNET_ERR_IO
            return
        end if
        if (allocated(global_model)) deallocate(global_model)
        allocate(global_model)
        call load_predictor_from_n2p2(trim(directory), global_model)
        if (size(global_model%species_names) /= number_of_types .or. &
            any([(trim(global_model%species_names(global_type)) /= trim(atom_types(global_type)), &
                  global_type=1,number_of_types)])) then
            deallocate(global_model)
            stat = ACCELNET_ERR_TYPE
            return
        end if
        pending_networks = global_model%networks
        potential_loaded = .true.
        call global_model%set_chebyshev_evaluation(chebyshev_evaluation_mode)
        is_loaded = .true.
        batch_cache_valid = .false.
        accelnet_Rc_min = global_model%minimum_distance
        accelnet_Rc_max = global_model%maximum_cutoff
        accelnet_nsf_max = maxval([(global_model%networks(species)%nodes(1), &
                                  species=1,number_of_types)])
        accelnet_nnb_max = lcl_nmax_nbdist(real(accelnet_Rc_min, real64), &
                                           real(accelnet_Rc_max, real64))
    end subroutine accelnet_load_n2p2

    subroutine accelnet_load_n2p2_c(directory, stat) bind(C, name="accelnet_load_n2p2")
        character(c_char), intent(in) :: directory(*)
        integer(c_int), intent(out) :: stat
        character(len=PATH_LENGTH) :: local_directory
        integer :: local_status
        call copy_c_string(directory, local_directory)
        call accelnet_load_n2p2(trim(local_directory), local_status)
        stat = local_status
    end subroutine accelnet_load_n2p2_c

    subroutine accelnet_load_potential_ascii_c(type_id, filename, stat) &
        bind(C, name="accelnet_load_potential_ascii")
        integer(c_int), value, intent(in) :: type_id
        character(c_char), intent(in) :: filename(*)
        integer(c_int), intent(out) :: stat
        character(len=PATH_LENGTH) :: local_filename
        integer :: local_status
        call copy_c_string(filename, local_filename)
        call accelnet_load_potential(type_id, trim(local_filename), local_status, is_ascii=.true.)
        stat = local_status
    end subroutine accelnet_load_potential_ascii_c

    logical(c_bool) function accelnet_all_loaded() bind(C)
        accelnet_all_loaded = is_initialized .and. is_loaded
    end function accelnet_all_loaded

    subroutine accelnet_print_info() bind(C)
        integer :: species
        if (.not. is_initialized) then
            write(*, "(A)") "AccelNet API is not initialized."
            return
        end if
        write(*, "(A,I0)") "AccelNet species: ", number_of_types
        write(*, "(2A)") "Chebyshev evaluation: ", &
            trim(chebyshev_evaluation_name(chebyshev_evaluation_mode))
        do species = 1, number_of_types
            write(*, "(I0,2A,L1)") species, " ", trim(atom_types(species)), potential_loaded(species)
        end do
        if (is_loaded) then
            write(*, "(A,ES14.6)") "Minimum radius: ", accelnet_Rc_min
            write(*, "(A,ES14.6)") "Maximum cutoff: ", accelnet_Rc_max
            write(*, "(A,I0)") "Maximum descriptors: ", accelnet_nsf_max
            write(*, "(A,I0)") "Maximum neighbors: ", accelnet_nnb_max
        end if
    end subroutine accelnet_print_info

    pure function chebyshev_evaluation_name(mode) result(name)
        integer, intent(in) :: mode
        character(len=6) :: name
        select case (mode)
        case (CHEBYSHEV_EVALUATION_DIRECT)
            name = "direct"
        case (CHEBYSHEV_EVALUATION_MOMENT)
            name = "moment"
        case default
            name = "auto"
        end select
    end function chebyshev_evaluation_name

    real(c_double) function accelnet_free_atom_energy(type_id) bind(C)
        integer(c_int), value, intent(in) :: type_id
        accelnet_free_atom_energy = 0.0_c_double
        if (.not. is_loaded) return
        if (type_id < 1 .or. type_id > number_of_types) return
        accelnet_free_atom_energy = global_model%networks(type_id)%atomic_references(type_id)
    end function accelnet_free_atom_energy

    subroutine accelnet_convert_atom_types(species_in, type_id_in, type_id_out, stat)
        character(len=*), intent(in) :: species_in(:)
        integer, intent(in) :: type_id_in(:)
        integer, intent(out) :: type_id_out(:)
        integer, intent(out) :: stat
        integer :: atom, source_type, target_type
        stat = ACCELNET_OK
        if (.not. is_initialized) then
            stat = ACCELNET_ERR_INIT
            return
        end if
        if (size(type_id_out) /= size(type_id_in)) then
            stat = ACCELNET_ERR_ARGUMENT
            return
        end if
        do atom = 1, size(type_id_in)
            source_type = type_id_in(atom)
            if (source_type < 1 .or. source_type > size(species_in)) then
                stat = ACCELNET_ERR_TYPE
                return
            end if
            target_type = find_species(species_in(source_type), atom_types)
            if (target_type == 0) then
                stat = ACCELNET_ERR_TYPE
                return
            end if
            type_id_out(atom) = target_type
        end do
    end subroutine accelnet_convert_atom_types

    subroutine accelnet_convert_atom_types_c(ntypes_in, species_pointers, natoms, type_id_in, &
                                             type_id_out, stat) bind(C, name="accelnet_convert_atom_types")
        integer(c_int), value, intent(in) :: ntypes_in, natoms
        type(c_ptr), intent(in) :: species_pointers(ntypes_in)
        integer(c_int), intent(in) :: type_id_in(natoms)
        integer(c_int), intent(out) :: type_id_out(natoms)
        integer(c_int), intent(out) :: stat
        character(len=TYPE_LENGTH) :: species_in(ntypes_in)
        integer :: i, local_status
        do i = 1, ntypes_in
            call copy_c_pointer_string(species_pointers(i), species_in(i))
        end do
        call accelnet_convert_atom_types(species_in, type_id_in, type_id_out, local_status)
        stat = local_status
    end subroutine accelnet_convert_atom_types_c

    subroutine accelnet_atomic_energy(coo_i, type_i, n_j, coo_j, type_j, energy_i, stat) bind(C)
        real(c_double), intent(in) :: coo_i(3)
        integer(c_int), value, intent(in) :: type_i, n_j
        real(c_double), intent(in) :: coo_j(3, n_j)
        integer(c_int), intent(in) :: type_j(n_j)
        real(c_double), intent(out) :: energy_i
        integer(c_int), intent(out) :: stat
        energy_i=0
        call validate_atomic_arguments(type_i,n_j,type_j,stat)
        if (stat/=ACCELNET_OK) return
        call prepare_atomic_model(stat)
        if (stat/=ACCELNET_OK) return
        call evaluate_atomic_energy_target(global_batch_model,type_i,type_j,coo_i,coo_j,energy_i, &
            global_shared_work,status=stat,reuse_model=.true.)
        if (stat/=0) then
            stat=ACCELNET_ERR_ARGUMENT
            return
        end if
    end subroutine accelnet_atomic_energy

    subroutine accelnet_atomic_energy_and_forces(coo_i, type_i, index_i, n_j, coo_j, type_j, &
                                                 index_j, natoms, energy_i, forces, stat) bind(C)
        integer(c_int), value, intent(in) :: type_i, index_i, n_j, natoms
        real(c_double), intent(in) :: coo_i(3), coo_j(3,n_j)
        integer(c_int), intent(in) :: type_j(n_j), index_j(n_j)
        real(c_double), intent(out) :: energy_i
        real(c_double), intent(inout) :: forces(3,natoms)
        integer(c_int), intent(out) :: stat
        call atomic_energy_forces(coo_i, type_i, index_i, n_j, coo_j, type_j, &
                                  index_j, natoms, energy_i, forces, stat)
    end subroutine accelnet_atomic_energy_and_forces

    subroutine accelnet_atomic_energy_and_forces_virial(coo_i, type_i, index_i, n_j, coo_j, type_j, &
                                                        index_j, natoms, energy_i, forces, virial, stat) bind(C)
        integer(c_int), value, intent(in) :: type_i, index_i, n_j, natoms
        real(c_double), intent(in) :: coo_i(3), coo_j(3,n_j)
        integer(c_int), intent(in) :: type_j(n_j), index_j(n_j)
        real(c_double), intent(out) :: energy_i
        real(c_double), intent(inout) :: forces(3,natoms), virial(3,3)
        integer(c_int), intent(out) :: stat
        call atomic_energy_forces(coo_i, type_i, index_i, n_j, coo_j, type_j, &
                                  index_j, natoms, energy_i, forces, stat, virial)
    end subroutine accelnet_atomic_energy_and_forces_virial

    ! CSR batch on the same loaded model as the atomic C API. Indices are
    ! one-based; targets include ghosts/images. All forces are additive.
    subroutine accelnet_batch_energy_and_forces(natoms,nrows,nedges,species,centers,offsets,indices, &
            displacements,energies,forces,stat) bind(C)
        integer(c_int), value, intent(in) :: natoms,nrows,nedges
        integer(c_int), intent(in) :: species(natoms),centers(nrows),offsets(nrows+1),indices(nedges)
        real(c_double), intent(in) :: displacements(3,nedges)
        real(c_double), intent(out) :: energies(nrows)
        real(c_double), intent(inout) :: forces(3,natoms)
        integer(c_int), intent(out) :: stat
        stat=ACCELNET_ERR_INIT
        if (.not. is_loaded .or. .not. allocated(global_model)) return
        stat=ACCELNET_ERR_ARGUMENT
        if (natoms<0 .or. nrows<0 .or. nedges<0) return
        if (any(species<1) .or. any(species>number_of_types)) return
        if (any(centers<1) .or. any(centers>natoms)) return
        if (any(indices<1) .or. any(indices>natoms)) return
        if (offsets(1)/=1 .or. offsets(nrows+1)/=nedges+1) return
        if (any(offsets<1) .or. any(offsets>nedges+1)) return
        if (any(offsets(2:)<offsets(:nrows))) return
        ! This model is private to the C API: only loading and mode setters
        ! can change it. Reuse packed metadata across small LAMMPS chunks;
        ! each such mutation invalidates the cache above.
        if (.not. batch_cache_valid) then
            call global_shared_work%release()
            call global_batch_model%initialize(global_model,use_host=.true.,status=batch_packed_status)
            batch_cache_valid = .true.
        end if
        if (batch_packed_status /= 0) return
        call evaluate_batch_target(global_batch_model,species,centers,offsets,indices,displacements, &
            energies,forces,global_shared_work,status=stat,reuse_model=.true.)
        if (stat/=0) stat=ACCELNET_ERR_ARGUMENT
    end subroutine

    subroutine atomic_energy_forces(coo_i, type_i, index_i, n_j, coo_j, type_j, &
                                    index_j, natoms, energy_i, forces, stat, virial)
        real(c_double), intent(in) :: coo_i(3)
        integer(c_int), value, intent(in) :: type_i, index_i, n_j, natoms
        real(c_double), intent(in) :: coo_j(3,n_j)
        integer(c_int), intent(in) :: type_j(n_j), index_j(n_j)
        real(c_double), intent(out) :: energy_i
        real(c_double), intent(inout) :: forces(3,natoms)
        real(c_double), intent(inout), optional :: virial(3,3)
        integer(c_int), intent(out) :: stat
        real(real64) :: energies(1)
        integer :: j
        energy_i=0
        call validate_atomic_arguments(type_i,n_j,type_j,stat)
        if (stat/=ACCELNET_OK) return
        if (natoms<1 .or. index_i<1 .or. index_i>natoms .or. &
            any(index_j<1) .or. any(index_j>natoms)) then
            stat=ACCELNET_ERR_ARGUMENT
            return
        end if
        call prepare_atomic_environment(coo_i,type_i,n_j,coo_j,type_j,stat)
        if (stat/=ACCELNET_OK) return
        call evaluate_batch_target(global_batch_model,atomic_species(:n_j+1),[1],[1,n_j+1], &
            atomic_indices(:n_j),atomic_displacements(:,:n_j),energies,atomic_forces(:,:n_j+1), &
            global_shared_work,virial,status=stat,reuse_model=.true.)
        if (stat/=0) then
            stat=ACCELNET_ERR_ARGUMENT
            return
        end if
        energy_i=energies(1)
        forces(:,index_i)=forces(:,index_i)+atomic_forces(:,1)
        do j=1,n_j
            forces(:,index_j(j))=forces(:,index_j(j))+atomic_forces(:,j+1)
        end do
    end subroutine atomic_energy_forces

    ! Map the local environment to one CSR row. Keep image slots distinct until
    ! final scatter so repeated physical indices and additive force semantics hold.
    subroutine prepare_atomic_environment(coo_i,type_i,n_j,coo_j,type_j,stat)
        real(real64), intent(in) :: coo_i(3),coo_j(3,n_j)
        integer, intent(in) :: type_i,n_j,type_j(n_j)
        integer(c_int), intent(out) :: stat
        integer :: j,allocation_status
        logical :: grow
        call validate_atomic_arguments(type_i,n_j,type_j,stat)
        if (stat/=ACCELNET_OK) return
        grow=.not.allocated(atomic_species)
        if (.not.grow) grow=size(atomic_species)<n_j+1
        if (grow) then
            if (allocated(atomic_species)) &
                deallocate(atomic_species,atomic_indices,atomic_displacements,atomic_forces)
            allocate(atomic_species(n_j+1),atomic_indices(n_j),atomic_displacements(3,n_j), &
                atomic_forces(3,n_j+1),stat=allocation_status)
            if (allocation_status/=0) then
                stat=ACCELNET_ERR_MALLOC
                return
            end if
        end if
        call prepare_atomic_model(stat)
        if (stat/=ACCELNET_OK) return
        atomic_species(1)=type_i; atomic_species(2:n_j+1)=type_j
        do j=1,n_j
            atomic_indices(j)=j+1
            atomic_displacements(:,j)=coo_j(:,j)-coo_i
        end do
        atomic_forces(:,:n_j+1)=0
    end subroutine

    subroutine prepare_atomic_model(stat)
        integer(c_int), intent(out) :: stat
        stat=ACCELNET_OK
        if (.not.batch_cache_valid) then
            call global_shared_work%release()
            call global_batch_model%initialize(global_model,use_host=.true.,status=batch_packed_status)
            batch_cache_valid=.true.
        end if
        if (batch_packed_status/=0) stat=ACCELNET_ERR_ARGUMENT
    end subroutine

    subroutine validate_atomic_arguments(type_i, n_j, type_j, stat)
        integer, intent(in) :: type_i, n_j, type_j(n_j)
        integer(c_int), intent(out) :: stat
        stat = ACCELNET_OK
        if (.not. is_loaded) then
            stat = ACCELNET_ERR_INIT
        else if (type_i < 1 .or. type_i > number_of_types .or. &
                 any(type_j < 1) .or. any(type_j > number_of_types)) then
            stat = ACCELNET_ERR_TYPE
        else if (n_j < 0) then
            stat = ACCELNET_ERR_ARGUMENT
        end if
    end subroutine validate_atomic_arguments

    subroutine accelnet_nbl_init(lattice, natoms, species, coordinates, cartesian, pbc) bind(C)
        real(c_double), intent(in) :: lattice(3,3)
        integer(c_int), value, intent(in) :: natoms
        integer(c_int), intent(in) :: species(natoms)
        real(c_double), intent(inout) :: coordinates(3,natoms)
        logical(c_bool), value, intent(in) :: cartesian, pbc
        if (.not. is_loaded) error stop "accelnet_nbl_init requires loaded potentials"
        if (neighbor_list_initialized) call accelnet_nbl_final()
        neighbor_structure%natoms = natoms
        neighbor_structure%pbc = pbc
        neighbor_structure%lattice = lattice
        allocate(neighbor_structure%positions(3,natoms), neighbor_structure%species(natoms))
        neighbor_structure%species = species
        if (cartesian) then
            neighbor_structure%positions = coordinates
        else
            neighbor_structure%positions = matmul(lattice, coordinates)
            coordinates = neighbor_structure%positions
        end if
        call build_neighbor_list(neighbor_structure, real(accelnet_Rc_max,real64), global_neighbors, &
                                 real(accelnet_Rc_min,real64))
        neighbor_list_initialized = .true.
    end subroutine accelnet_nbl_init

    subroutine accelnet_nbl_final() bind(C)
        if (allocated(neighbor_structure%positions)) deallocate(neighbor_structure%positions)
        if (allocated(neighbor_structure%species)) deallocate(neighbor_structure%species)
        if (allocated(global_neighbors%offsets)) deallocate(global_neighbors%offsets)
        if (allocated(global_neighbors%atom_indices)) deallocate(global_neighbors%atom_indices)
        if (allocated(global_neighbors%image_shifts)) deallocate(global_neighbors%image_shifts)
        if (allocated(global_neighbors%positions)) deallocate(global_neighbors%positions)
        if (allocated(global_neighbors%displacements)) deallocate(global_neighbors%displacements)
        neighbor_list_initialized = .false.
    end subroutine accelnet_nbl_final

    subroutine accelnet_nbl_neighbors(iatom, nnb, nbcoo, nbdist, nblist, nbtype) bind(C)
        integer(c_int), value, intent(in) :: iatom
        integer(c_int), intent(inout) :: nnb
        real(c_double), intent(out) :: nbcoo(3,nnb), nbdist(nnb)
        integer(c_int), intent(out) :: nblist(nnb), nbtype(nnb)
        integer :: first, last, actual, entry
        if (.not. neighbor_list_initialized) error stop "accelnet neighbor list is not initialized"
        if (iatom < 1 .or. iatom > neighbor_structure%natoms) error stop "central atom index is out of range"
        first = global_neighbors%offsets(iatom)
        last = global_neighbors%offsets(iatom + 1) - 1
        actual = max(0, last - first + 1)
        if (actual > nnb) error stop "neighbor output arrays are too small"
        do entry = 1, actual
            nbcoo(:,entry) = global_neighbors%positions(:,first + entry - 1)
            nbdist(entry) = sqrt(sum(global_neighbors%displacements(:,first + entry - 1)**2))
            nblist(entry) = global_neighbors%atom_indices(first + entry - 1)
            nbtype(entry) = neighbor_structure%species(nblist(entry))
        end do
        nnb = actual
    end subroutine accelnet_nbl_neighbors

    subroutine accelnet_sfb_init(species, radial_order, angular_order, radial_cutoff, angular_cutoff, stat)
        character(len=*), intent(in) :: species(:)
        integer, intent(in) :: radial_order, angular_order
        real(real64), intent(in) :: radial_cutoff, angular_cutoff
        integer, intent(out) :: stat
        type(predictor_model) :: basis_model
        integer :: i,j,d,ns
        stat = ACCELNET_OK
        if (sfb_initialized) then
            stat = ACCELNET_ERR_INIT
            return
        end if
        if (size(species) < 1 .or. radial_order < 0 .or. angular_order < 0 .or. &
            radial_cutoff <= 0.0_real64 .or. angular_cutoff <= 0.0_real64) then
            stat = ACCELNET_ERR_ARGUMENT
            return
        end if
        call initialize_config(sfb_config, size(species), radial_cutoff, radial_order, &
                               angular_cutoff, angular_order, version=0)
        ns=size(species); d=sfb_config%num_descriptors()
        allocate(basis_model%setups(ns),basis_model%networks(ns))
        do i=1,ns
            basis_model%setups(i)%global_to_local=[(j,j=1,ns)]
            call add_chebyshev(basis_model%setups(i)%model,sfb_config)
            basis_model%networks(i)%nlayers=2
            basis_model%networks(i)%nodes=[d,1]
            basis_model%networks(i)%activation=[0]
            basis_model%networks(i)%weight_offsets=[0]
            allocate(basis_model%networks(i)%weights(d+1),basis_model%networks(i)%descriptor_shift(d), &
                basis_model%networks(i)%descriptor_scale(d),basis_model%networks(i)%atomic_references(ns))
            basis_model%networks(i)%weights=0; basis_model%networks(i)%descriptor_shift=0
            basis_model%networks(i)%descriptor_scale=1; basis_model%networks(i)%atomic_references=0
        end do
        call sfb_model%initialize(basis_model,use_host=.true.,status=stat)
        if (stat/=0) then
            stat=ACCELNET_ERR_ARGUMENT
            return
        end if
        sfb_radial_cutoff = radial_cutoff
        sfb_initialized = .true.
    end subroutine accelnet_sfb_init

    subroutine accelnet_sfb_init_c(ntypes, species_pointers, radial_order, angular_order, &
                                   radial_cutoff, angular_cutoff, stat) bind(C, name="accelnet_sfb_init")
        integer(c_int), value, intent(in) :: ntypes, radial_order, angular_order
        type(c_ptr), intent(in) :: species_pointers(ntypes)
        real(c_double), value, intent(in) :: radial_cutoff, angular_cutoff
        integer(c_int), intent(out) :: stat
        character(len=TYPE_LENGTH) :: species(ntypes)
        integer :: i, local_status
        do i = 1, ntypes
            call copy_c_pointer_string(species_pointers(i), species(i))
        end do
        call accelnet_sfb_init(species, radial_order, angular_order, radial_cutoff, angular_cutoff, local_status)
        stat = local_status
    end subroutine accelnet_sfb_init_c

    subroutine accelnet_sfb_final(stat) bind(C)
        integer(c_int), intent(out) :: stat
        stat = ACCELNET_OK
        if (.not. sfb_initialized) then
            stat = ACCELNET_ERR_INIT
            return
        end if
        call sfb_work%release()
        call sfb_model%release()
        sfb_initialized = .false.
        sfb_radial_cutoff = 0.0_real64
    end subroutine accelnet_sfb_final

    integer(c_int) function accelnet_sfb_nvalues() bind(C)
        if (sfb_initialized) then
            accelnet_sfb_nvalues = sfb_config%num_descriptors()
        else
            accelnet_sfb_nvalues = 0
        end if
    end function accelnet_sfb_nvalues

    subroutine accelnet_sfb_eval(type_i, coo_i, n_j, type_j, coo_j, nvalues, values, stat) bind(C)
        integer(c_int), value, intent(in) :: type_i, n_j, nvalues
        real(c_double), intent(in) :: coo_i(3), coo_j(3,n_j)
        integer(c_int), intent(in) :: type_j(n_j)
        real(c_double), intent(inout) :: values(nvalues)
        integer(c_int), intent(out) :: stat
        real(real64), allocatable :: displacements(:, :)
        real(real64), allocatable :: raw(:,:)
        real(real64) :: energy(1),force(3,n_j+1)
        integer :: j
        stat = ACCELNET_OK
        if (.not. sfb_initialized) then
            stat = ACCELNET_ERR_INIT
            return
        end if
        if (type_i < 1 .or. type_i > sfb_config%num_species .or. &
            any(type_j < 1) .or. any(type_j > sfb_config%num_species)) then
            stat = ACCELNET_ERR_TYPE
            return
        end if
        if (nvalues < sfb_config%num_descriptors()) then
            stat = ACCELNET_ERR_ARGUMENT
            return
        end if
        allocate(displacements(3,n_j),raw(1,sfb_config%num_descriptors()))
        displacements = coo_j - spread(coo_i,2,n_j)
        force=0
        call evaluate_batch_target(sfb_model,[type_i,type_j],[1],[1,n_j+1],[(j+1,j=1,n_j)], &
            displacements,energy,force,sfb_work,status=stat,descriptor_values=raw)
        if (stat/=0) then
            stat=ACCELNET_ERR_ARGUMENT
            return
        end if
        values(:size(raw,2))=raw(1,:)
    end subroutine accelnet_sfb_eval

    subroutine accelnet_sfb_reconstruct_radial(nvalues, values, nx, x, y, stat) bind(C)
        integer(c_int), value, intent(in) :: nvalues, nx
        real(c_double), intent(in) :: values(nvalues)
        real(c_double), intent(out) :: x(nx), y(nx)
        integer(c_int), intent(out) :: stat
        real(real64) :: basis(sfb_config%radial_order + 1), ratio, weight, dx
        integer :: point, coefficient, nr
        stat = ACCELNET_OK
        if (.not. sfb_initialized) then
            stat = ACCELNET_ERR_INIT
            return
        end if
        nr = sfb_config%radial_order + 1
        if (nvalues < nr .or. nx < 2) then
            stat = ACCELNET_ERR_ARGUMENT
            return
        end if
        dx = sfb_radial_cutoff/real(nx - 1,real64)
        do point = 1, nx
            x(point) = real(point - 1,real64)*dx
        end do
        do point = 1, nx - 1
            call chebyshev_values(x(point),0.0_real64,sfb_radial_cutoff,sfb_config%radial_order,basis)
            ratio = x(point)/sfb_radial_cutoff
            if (ratio <= 0.0_real64) then
                y(point) = 0.0_real64
                cycle
            end if
            weight = 0.5_real64/(PI_COMPAT*sqrt(ratio - ratio*ratio))
            basis = basis*weight
            basis(1) = 0.5_real64*basis(1)
            y(point) = 0.0_real64
            do coefficient = nr, 1, -1
                y(point) = y(point) + values(coefficient)*basis(coefficient)
            end do
        end do
        y(nx) = 0.0_real64
    end subroutine accelnet_sfb_reconstruct_radial

    integer function find_species(name, species) result(position)
        character(len=*), intent(in) :: name, species(:)
        integer :: i
        position = 0
        do i = 1, size(species)
            if (trim(name) == trim(species(i))) then
                position = i
                return
            end if
        end do
    end function find_species

    subroutine copy_c_pointer_string(pointer, output)
        type(c_ptr), intent(in) :: pointer
        character(len=*), intent(out) :: output
        character(c_char), pointer :: input(:)
        integer :: i
        call c_f_pointer(pointer,input,[len(output)+1])
        output = ""
        do i = 1, len(output)
            if (input(i) == c_null_char) exit
            output(i:i) = input(i)
        end do
    end subroutine copy_c_pointer_string

    subroutine copy_c_string(input, output)
        character(c_char), intent(in) :: input(*)
        character(len=*), intent(out) :: output
        integer :: i
        output = ""
        do i = 1, len(output)
            if (input(i) == c_null_char) exit
            output(i:i) = input(i)
        end do
    end subroutine copy_c_string

end module accelnet
