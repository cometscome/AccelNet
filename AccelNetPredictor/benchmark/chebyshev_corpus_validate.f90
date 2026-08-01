program chebyshev_corpus_validate
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model, load_predictor
    use accelnet_descriptors, only: atomic_structure, read_xsf
    use aenet, only: aenet_init, aenet_final, aenet_load_potential, aenet_convert_atom_types, &
        aenet_atomic_energy_and_forces, aenet_nnb_max, aenet_Rc_min, aenet_Rc_max
    use geometry, only: geo_init, geo_final, pbc, latticeVec, nAtoms, atomType, atomTypeName, cooLatt
    use lclist, only: lcl_init, lcl_final, lcl_nbdist_cart
    implicit none

    real(real64), parameter :: ENERGY_TOLERANCE = 1.0e-8_real64
    real(real64), parameter :: FORCE_TOLERANCE = 1.0e-8_real64
    integer, parameter :: MAX_REPORTED_FAILURES = 20
    type(predictor_model) :: model
    type(atomic_structure) :: structure
    character(len=1024), allocatable :: setups(:), networks(:)
    character(len=1024) :: argument, corpus_dir, filename
    character(len=2), allocatable :: species_names(:)
    integer, allocatable :: original_species(:), converted_species(:)
    integer, allocatable :: neighbor_atoms(:), neighbor_species(:)
    real(real64), allocatable :: neighbor_coordinates(:, :), neighbor_distances(:)
    real(real64), allocatable :: predictor_forces(:, :), aenet_forces(:, :)
    real(real64) :: predictor_energy, aenet_energy, energy_error, force_error
    real(real64) :: maximum_energy_error, maximum_force_error
    real(real64) :: predictor_energy_sum, aenet_energy_sum
    real(real64) :: predictor_force_checksum, aenet_force_checksum
    real(real64) :: atomic_energy, center(3)
    integer :: count, nstructures, species, structure_index, atom, n, status, local_status
    integer :: energy_failure_count, force_failure_count, reported_failure_count
    integer :: maximum_energy_structure, maximum_force_structure, maximum_force_atom, maximum_force_component
    integer :: force_location(2)

    if (command_argument_count() < 7) then
        write(*, "(A)") &
            "usage: chebyshev-corpus-validate NSPECIES SETUP... NETWORK... DIRECTORY NSTRUCTURES"
        error stop 2
    end if
    call get_command_argument(1, argument); read(argument, *) count
    if (command_argument_count() /= 2*count + 3) error stop "invalid number of arguments"
    allocate(setups(count), networks(count), species_names(count))
    do species = 1, count
        call get_command_argument(1 + species, setups(species))
        call get_command_argument(1 + count + species, networks(species))
    end do
    call get_command_argument(2*count + 2, corpus_dir)
    call get_command_argument(2*count + 3, argument); read(argument, *) nstructures
    if (nstructures < 1 .or. nstructures > 9999) error stop "invalid NSTRUCTURES"

    call load_predictor(setups, networks, model)
    species_names = model%species_names
    call aenet_init(species_names, status); call check_status("aenet_init", status)
    do species = 1, count
        call aenet_load_potential(species, trim(networks(species)), status, is_ascii=.true.)
        call check_status("aenet_load_potential", status)
    end do
    allocate(neighbor_coordinates(3, aenet_nnb_max), neighbor_distances(aenet_nnb_max), &
             neighbor_atoms(aenet_nnb_max), neighbor_species(aenet_nnb_max))

    maximum_energy_error = 0.0_real64
    maximum_force_error = 0.0_real64
    maximum_energy_structure = 0
    maximum_force_structure = 0
    maximum_force_atom = 0
    maximum_force_component = 0
    energy_failure_count = 0
    force_failure_count = 0
    reported_failure_count = 0
    predictor_energy_sum = 0.0_real64
    aenet_energy_sum = 0.0_real64
    predictor_force_checksum = 0.0_real64
    aenet_force_checksum = 0.0_real64

    do structure_index = 1, nstructures
        call structure_filename(corpus_dir, structure_index, filename)
        call read_xsf(trim(filename), model%species_names, structure)
        allocate(predictor_forces(3, structure%natoms))
        call model%predict_energy_forces(structure, predictor_energy, predictor_forces)

        call geo_init(trim(filename), "xsf")
        if (nAtoms /= structure%natoms) error stop "XSF readers disagree on atom count"
        allocate(original_species(nAtoms), converted_species(nAtoms), aenet_forces(3, nAtoms))
        original_species = atomType
        call aenet_convert_atom_types(atomTypeName, original_species, converted_species, status)
        call check_status("aenet_convert_atom_types", status)
        aenet_energy = 0.0_real64
        aenet_forces = 0.0_real64
        call lcl_init(aenet_Rc_min, aenet_Rc_max, latticeVec, nAtoms, converted_species, cooLatt, pbc)
        do atom = 1, nAtoms
            center = matmul(latticeVec, cooLatt(:, atom))
            n = aenet_nnb_max
            call lcl_nbdist_cart(atom, n, neighbor_coordinates, neighbor_distances, aenet_Rc_max, &
                                 nblist=neighbor_atoms, nbtype=neighbor_species)
            call aenet_atomic_energy_and_forces(center, converted_species(atom), atom, n, &
                neighbor_coordinates, neighbor_species, neighbor_atoms, nAtoms, atomic_energy, &
                aenet_forces, local_status)
            call check_status("aenet_atomic_energy_and_forces", local_status)
            aenet_energy = aenet_energy + atomic_energy
        end do
        call lcl_final()
        call geo_final()

        energy_error = abs(predictor_energy - aenet_energy)
        force_error = maxval(abs(predictor_forces - aenet_forces))
        if (energy_error > maximum_energy_error) then
            maximum_energy_error = energy_error
            maximum_energy_structure = structure_index
        end if
        if (force_error > maximum_force_error) then
            maximum_force_error = force_error
            maximum_force_structure = structure_index
            force_location = maxloc(abs(predictor_forces - aenet_forces))
            maximum_force_component = force_location(1)
            maximum_force_atom = force_location(2)
        end if
        if (energy_error > ENERGY_TOLERANCE) then
            energy_failure_count = energy_failure_count + 1
            call report_failure("energy", structure_index, energy_error)
        end if
        if (force_error > FORCE_TOLERANCE) then
            force_failure_count = force_failure_count + 1
            call report_failure("force", structure_index, force_error)
        end if
        predictor_energy_sum = predictor_energy_sum + predictor_energy
        aenet_energy_sum = aenet_energy_sum + aenet_energy
        predictor_force_checksum = predictor_force_checksum + sum(abs(predictor_forces))
        aenet_force_checksum = aenet_force_checksum + sum(abs(aenet_forces))
        deallocate(predictor_forces, aenet_forces, original_species, converted_species)
    end do

    call aenet_final(status)
    write(*, "(A,1X,I0)") "STRUCTURES_CHECKED", nstructures
    write(*, "(A,1X,ES24.16,1X,A,1X,I0)") &
        "MAX_ENERGY_ABS_ERROR_EV", maximum_energy_error, "STRUCTURE", maximum_energy_structure
    write(*, "(A,1X,ES24.16,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,I0)") &
        "MAX_FORCE_COMPONENT_ABS_ERROR", maximum_force_error, "STRUCTURE", maximum_force_structure, &
        "ATOM", maximum_force_atom, "COMPONENT", maximum_force_component
    write(*, "(A,1X,I0)") "ENERGY_FAILURES", energy_failure_count
    write(*, "(A,1X,I0)") "FORCE_FAILURES", force_failure_count
    write(*, "(A,1X,ES24.16)") "PREDICTOR_ENERGY_SUM_EV", predictor_energy_sum
    write(*, "(A,1X,ES24.16)") "AENET_ENERGY_SUM_EV", aenet_energy_sum
    write(*, "(A,1X,ES24.16)") "PREDICTOR_FORCE_ABS_CHECKSUM", predictor_force_checksum
    write(*, "(A,1X,ES24.16)") "AENET_FORCE_ABS_CHECKSUM", aenet_force_checksum
    if (energy_failure_count > 0 .or. force_failure_count > 0) &
        error stop "Chebyshev corpus validation failed"

contains
    subroutine report_failure(kind, index, error_value)
        character(len=*), intent(in) :: kind
        integer, intent(in) :: index
        real(real64), intent(in) :: error_value
        if (reported_failure_count < MAX_REPORTED_FAILURES) then
            write(*, "(A,1X,A,1X,I0,1X,ES24.16)") "FAILURE", trim(kind), index, error_value
            reported_failure_count = reported_failure_count + 1
        end if
    end subroutine report_failure

    subroutine structure_filename(directory, index, path)
        character(len=*), intent(in) :: directory
        integer, intent(in) :: index
        character(len=*), intent(out) :: path
        character(len=4) :: number
        write(number, "(I4.4)") index
        path = trim(directory)//"/structure"//number//".xsf"
    end subroutine structure_filename

    subroutine check_status(operation, value)
        character(len=*), intent(in) :: operation
        integer, intent(in) :: value
        if (value /= 0) then
            write(*, "(A,1X,A,1X,I0)") "aenet operation failed:", trim(operation), value
            error stop "aenet validation failed"
        end if
    end subroutine check_status
end program chebyshev_corpus_validate
