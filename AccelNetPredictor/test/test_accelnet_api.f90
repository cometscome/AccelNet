program test_accelnet_api
    use iso_c_binding, only: c_bool, c_double, c_int
    use iso_fortran_env, only: real64
    use accelnet
    use accelnet_descriptors, only: atomic_structure, read_xsf
    use accelnet_predictor, only: predictor_model, load_predictor_from_networks
    implicit none

    character(len=1024) :: network_files(2), structure_file
    character(len=16) :: species(2), input_species(2)
    type(atomic_structure) :: structure
    type(predictor_model) :: reference_model
    real(real64), allocatable :: reference_forces(:,:), api_forces(:,:)
    real(c_double), allocatable :: nbcoo(:,:), nbdist(:), values(:), x(:), y(:)
    integer(c_int), allocatable :: nblist(:), nbtype(:)
    integer(c_int) :: stat, nnb, capacity, atom, nvalues
    integer :: converted(2), input_ids(2), evaluation_mode
    logical :: loaded
    real(real64) :: reference_energy, api_energy, atomic_energy, force_energy

    call get_command_argument(1, network_files(1))
    call get_command_argument(2, network_files(2))
    call get_command_argument(3, structure_file)
    species = [character(len=16) :: "Ti", "O"]

    call load_predictor_from_networks(network_files, reference_model, chebyshev_version=0)
    call read_xsf(trim(structure_file), species, structure)
    allocate(reference_forces(3,structure%natoms))
    call reference_model%predict_energy_forces(structure, reference_energy, reference_forces)

    call accelnet_init(species, stat)
    call require(stat == ACCELNET_OK, "init")
    loaded = accelnet_all_loaded()
    call require(.not. loaded, "all_loaded before load")
    call accelnet_set_chebyshev_version(0_c_int, stat)
    call require(stat == ACCELNET_OK, "Chebyshev version")
    call require(accelnet_get_chebyshev_evaluation() == ACCELNET_CHEBYSHEV_AUTO, &
        "Chebyshev default mode")
    call accelnet_set_chebyshev_evaluation(ACCELNET_CHEBYSHEV_DIRECT, stat)
    call require(stat == ACCELNET_OK .and. &
        accelnet_get_chebyshev_evaluation() == ACCELNET_CHEBYSHEV_DIRECT, &
        "Chebyshev direct before load")
    call accelnet_set_chebyshev_evaluation(3_c_int, stat)
    call require(stat == ACCELNET_ERR_ARGUMENT, "invalid Chebyshev mode")
    call accelnet_set_chebyshev_evaluation(ACCELNET_CHEBYSHEV_MOMENT, stat)
    call require(stat == ACCELNET_OK, "Chebyshev moment before load")
    call accelnet_set_g5_evaluation(ACCELNET_G5_DIRECT, stat)
    call require(stat == ACCELNET_ERR_INIT, "G5 mode before load")
    call accelnet_load_potential(1, trim(network_files(1)), stat, is_ascii=.true.)
    call require(stat == ACCELNET_OK, "load Ti")
    call accelnet_load_potential(2, trim(network_files(2)), stat, is_ascii=.true.)
    call require(stat == ACCELNET_OK, "load O")
    loaded = accelnet_all_loaded()
    call require(loaded, "all_loaded after load")
    call require(accelnet_get_chebyshev_evaluation() == ACCELNET_CHEBYSHEV_MOMENT, &
        "Chebyshev moment after load")
    call accelnet_set_chebyshev_evaluation(ACCELNET_CHEBYSHEV_DIRECT, stat)
    call require(stat == ACCELNET_OK .and. &
        accelnet_get_chebyshev_evaluation() == ACCELNET_CHEBYSHEV_DIRECT, &
        "Chebyshev direct after load")
    call accelnet_set_chebyshev_evaluation(ACCELNET_CHEBYSHEV_AUTO, stat)
    call require(stat == ACCELNET_OK, "Chebyshev auto after load")
    call accelnet_set_g5_evaluation(ACCELNET_G5_DIRECT, stat)
    call require(stat == ACCELNET_OK, "G5 direct mode")
    call accelnet_set_g5_evaluation(ACCELNET_G5_MOMENT_FORCE, stat)
    call require(stat == ACCELNET_OK, "G5 forced-moment mode")
    call accelnet_set_g5_evaluation(4_c_int, stat)
    call require(stat == ACCELNET_ERR_ARGUMENT, "invalid G5 mode")
    call accelnet_set_g5_evaluation(ACCELNET_G5_AUTO, stat)
    call require(stat == ACCELNET_OK, "G5 auto mode")
    call require(accelnet_nsf_max > 0 .and. accelnet_nnb_max > 0, "published maxima")
    call require(accelnet_Rc_min > 0.0_c_double .and. accelnet_Rc_max > accelnet_Rc_min, "published cutoffs")
    call require(abs(accelnet_free_atom_energy(1_c_int) - &
        reference_model%networks(1)%atomic_references(1)) < 1.0e-12_real64, "free atom energy")

    input_species = [character(len=16) :: "O", "Ti"]
    input_ids = [1, 2]
    call accelnet_convert_atom_types(input_species, input_ids, converted, stat)
    call require(stat == ACCELNET_OK .and. all(converted == [2, 1]), "atom type conversion")

    call accelnet_nbl_init(structure%lattice, int(structure%natoms,c_int), &
        int(structure%species,c_int), structure%positions, .true._c_bool, .true._c_bool)
    capacity = max(accelnet_nnb_max, int(structure%natoms,c_int))
    allocate(nbcoo(3,capacity), nbdist(capacity), nblist(capacity), nbtype(capacity))
    allocate(api_forces(3,structure%natoms))
    do evaluation_mode = ACCELNET_CHEBYSHEV_AUTO, ACCELNET_CHEBYSHEV_MOMENT
        call accelnet_set_chebyshev_evaluation(int(evaluation_mode,c_int), stat)
        call require(stat == ACCELNET_OK .and. &
            accelnet_get_chebyshev_evaluation() == evaluation_mode, "select Chebyshev evaluation mode")
        api_forces = 0.0_real64
        api_energy = 0.0_real64
        force_energy = 0.0_real64
        do atom = 1, structure%natoms
            nnb = capacity
            call accelnet_nbl_neighbors(int(atom,c_int), nnb, nbcoo, nbdist, nblist, nbtype)
            call accelnet_atomic_energy(structure%positions(:,atom), int(structure%species(atom),c_int), &
                nnb, nbcoo(:,1:nnb), nbtype(1:nnb), atomic_energy, stat)
            call require(stat == ACCELNET_OK, "atomic energy")
            api_energy = api_energy + atomic_energy
            call accelnet_atomic_energy_and_forces(structure%positions(:,atom), &
                int(structure%species(atom),c_int), int(atom,c_int), nnb, nbcoo(:,1:nnb), &
                nbtype(1:nnb), nblist(1:nnb), int(structure%natoms,c_int), atomic_energy, api_forces, stat)
            call require(stat == ACCELNET_OK, "atomic energy and forces")
            force_energy = force_energy + atomic_energy
        end do
        call require(abs(api_energy-reference_energy) < 1.0e-9_real64, "energy parity by mode")
        call require(abs(force_energy-reference_energy) < 1.0e-9_real64, "force-call energy parity by mode")
        call require(maxval(abs(api_forces-reference_forces)) < 1.0e-9_real64, "force parity by mode")
    end do
    call accelnet_nbl_final()

    call accelnet_final(stat)
    loaded = accelnet_all_loaded()
    call require(stat == ACCELNET_OK .and. .not. loaded, "final")

    call accelnet_sfb_init(species, 3, 2, 5.0_real64, 4.0_real64, stat)
    call require(stat == ACCELNET_OK, "SFB init")
    nvalues = accelnet_sfb_nvalues()
    call require(nvalues > 0, "SFB size")
    allocate(values(nvalues), x(17), y(17)); values = 0.0_real64
    call accelnet_sfb_eval(1_c_int, structure%positions(:,1), 1_c_int, &
        int(structure%species(2:2),c_int), structure%positions(:,2:2), nvalues, values, stat)
    call require(stat == ACCELNET_OK .and. all(values == values), "SFB eval")
    call accelnet_sfb_reconstruct_radial(nvalues, values, 17_c_int, x, y, stat)
    call require(stat == ACCELNET_OK .and. all(y == y), "SFB reconstruction")
    call accelnet_sfb_final(stat)
    call require(stat == ACCELNET_OK .and. accelnet_sfb_nvalues() == 0, "SFB final")

    write(*,'(A,ES14.6,A,ES14.6)') "API energy error=", abs(api_energy-reference_energy), &
        " force max error=", maxval(abs(api_forces-reference_forces))

contains
    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,'(2A)') "FAILED: ", trim(label)
            error stop 1
        end if
    end subroutine require
end program test_accelnet_api
