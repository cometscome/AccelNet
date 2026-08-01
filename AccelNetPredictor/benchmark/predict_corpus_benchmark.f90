program predict_corpus_benchmark
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model, load_predictor
    use accelnet_descriptors, only: atomic_structure, read_xsf
    implicit none

    type(predictor_model) :: model
    type(atomic_structure), allocatable :: structures(:)
    character(len=1024), allocatable :: setups(:), networks(:)
    character(len=1024) :: argument, corpus_dir, filename
    real(real64), allocatable :: forces(:, :)
    real(real64) :: energy, energy_checksum, force_checksum
    real(real64) :: start_time, end_time, energy_time, force_time
    integer :: count, repeats, nstructures, species, structure_index, iteration, max_atoms

    if (command_argument_count() < 7) then
        write(*, "(A)") &
            "usage: accelnet-corpus-benchmark NSPECIES SETUP... NETWORK... DIRECTORY NSTRUCTURES REPEATS"
        error stop 2
    end if
    call get_command_argument(1, argument); read(argument, *) count
    if (command_argument_count() /= 2*count + 4) error stop "invalid number of arguments"
    allocate(setups(count), networks(count))
    do species = 1, count
        call get_command_argument(1 + species, setups(species))
        call get_command_argument(1 + count + species, networks(species))
    end do
    call get_command_argument(2*count + 2, corpus_dir)
    call get_command_argument(2*count + 3, argument); read(argument, *) nstructures
    call get_command_argument(2*count + 4, argument); read(argument, *) repeats
    if (nstructures < 1 .or. nstructures > 9999) error stop "invalid NSTRUCTURES"
    if (repeats < 1) error stop "REPEATS must be positive"

    call load_predictor(setups, networks, model)
    allocate(structures(nstructures))
    max_atoms = 0
    do structure_index = 1, nstructures
        call structure_filename(corpus_dir, structure_index, filename)
        call read_xsf(trim(filename), model%species_names, structures(structure_index))
        max_atoms = max(max_atoms, structures(structure_index)%natoms)
    end do
    allocate(forces(3, max_atoms))

    ! Untimed full-corpus pass both warms the implementation and records checksums.
    energy_checksum = 0.0_real64
    force_checksum = 0.0_real64
    do structure_index = 1, nstructures
        call model%predict_energy_forces(structures(structure_index), energy, &
            forces(:, 1:structures(structure_index)%natoms))
        energy_checksum = energy_checksum + energy
        force_checksum = force_checksum + sum(abs(forces(:, 1:structures(structure_index)%natoms)))
    end do

    call cpu_time(start_time)
    do iteration = 1, repeats
        do structure_index = 1, nstructures
            call model%predict_energy(structures(structure_index), energy)
        end do
    end do
    call cpu_time(end_time)
    energy_time = end_time - start_time

    call cpu_time(start_time)
    do iteration = 1, repeats
        do structure_index = 1, nstructures
            call model%predict_energy_forces(structures(structure_index), energy, &
                forces(:, 1:structures(structure_index)%natoms))
        end do
    end do
    call cpu_time(end_time)
    force_time = end_time - start_time

    write(*, "(A,1X,I0)") "STRUCTURES", nstructures
    write(*, "(A,1X,I0)") "REPEATS", repeats
    write(*, "(A,1X,ES24.16)") "ENERGY_CHECKSUM_EV", energy_checksum
    write(*, "(A,1X,ES24.16)") "FORCE_ABS_CHECKSUM", force_checksum
    write(*, "(A,1X,ES24.16)") "ENERGY_SECONDS_PER_STRUCTURE", &
        energy_time/real(nstructures*repeats, real64)
    write(*, "(A,1X,ES24.16)") "ENERGY_FORCE_SECONDS_PER_STRUCTURE", &
        force_time/real(nstructures*repeats, real64)

contains
    subroutine structure_filename(directory, index, path)
        character(len=*), intent(in) :: directory
        integer, intent(in) :: index
        character(len=*), intent(out) :: path
        character(len=4) :: number
        write(number, "(I4.4)") index
        path = trim(directory)//"/structure"//number//".xsf"
    end subroutine structure_filename
end program predict_corpus_benchmark
