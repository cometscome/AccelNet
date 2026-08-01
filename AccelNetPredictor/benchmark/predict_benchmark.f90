program predict_benchmark
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model, load_predictor, load_predictor_from_networks
    use accelnet_descriptors, only: atomic_structure, read_xsf
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: structure
    character(len=1024), allocatable :: setups(:), networks(:)
    character(len=1024) :: argument, xsf
    integer :: count, repeats, species, iteration
    real(real64), allocatable :: forces(:, :)
    real(real64) :: energy, start_time, end_time, energy_api_time, force_api_time, file_time
    if (command_argument_count() < 5) then
        write(*, "(A)") "usage: accelnet-predict-benchmark NSPECIES SETUP... NETWORK... XSF REPEATS"
        error stop 2
    end if
    call get_command_argument(1, argument); read(argument, *) count
    if (command_argument_count() /= 2*count + 3) error stop "invalid number of arguments"
    allocate(setups(count), networks(count))
    do species = 1, count
        call get_command_argument(1 + species, setups(species))
        call get_command_argument(1 + count + species, networks(species))
    end do
    call get_command_argument(2*count + 2, xsf)
    call get_command_argument(2*count + 3, argument); read(argument, *) repeats
    if (repeats < 1) error stop "REPEATS must be positive"
    if (index(trim(networks(1)), ".ascii") > 0) then
        call load_predictor(setups, networks, model)
    else
        call load_predictor_from_networks(networks, model)
    end if
    call read_xsf(trim(xsf), model%species_names, structure)
    allocate(forces(3, structure%natoms))
    call model%predict_energy(structure, energy) ! warm-up
    call cpu_time(start_time)
    do iteration = 1, repeats
        call model%predict_energy(structure, energy)
    end do
    call cpu_time(end_time)
    energy_api_time = end_time - start_time
    call model%predict_energy_forces(structure, energy, forces) ! warm-up
    call cpu_time(start_time)
    do iteration = 1, repeats
        call model%predict_energy_forces(structure, energy, forces)
    end do
    call cpu_time(end_time)
    force_api_time = end_time - start_time
    call model%predict_energy(xsf, energy) ! warm-up
    call cpu_time(start_time)
    do iteration = 1, repeats
        call model%predict_energy(xsf, energy)
    end do
    call cpu_time(end_time)
    file_time = end_time - start_time
    write(*, "(A,1X,ES24.16)") "TOTAL_ENERGY_EV", energy
    write(*, "(A,1X,ES24.16)") "FORCE_ABS_CHECKSUM", sum(abs(forces))
    write(*, "(A,1X,I0)") "REPEATS", repeats
    write(*, "(A,1X,ES24.16)") "ENERGY_API_SECONDS_PER_STRUCTURE", energy_api_time/repeats
    write(*, "(A,1X,ES24.16)") "ENERGY_FORCE_API_SECONDS_PER_STRUCTURE", force_api_time/repeats
    write(*, "(A,1X,ES24.16)") "ENERGY_FILE_SECONDS_PER_STRUCTURE", file_time/repeats
end program predict_benchmark
