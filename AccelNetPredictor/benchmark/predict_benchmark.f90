program predict_benchmark
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model, load_predictor, load_predictor_from_networks, &
        load_predictor_from_n2p2
    use accelnet_descriptors, only: atomic_structure, read_xsf
    use accelnet_behler, only: G5_EVALUATION_AUTO, G5_EVALUATION_DIRECT, G5_EVALUATION_MOMENT
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: structure
    character(len=1024), allocatable :: setups(:), networks(:)
    character(len=1024) :: argument, xsf, model_directory
    integer :: count, repeats, species, iteration, g5_mode
    real(real64), allocatable :: forces(:, :)
    real(real64) :: energy, start_time, end_time, energy_api_time, force_api_time, file_time
    g5_mode = G5_EVALUATION_AUTO
    if (command_argument_count() == 4 .or. command_argument_count() == 5) then
        call get_command_argument(1, argument)
        if (trim(argument) /= "--n2p2") call usage()
        call get_command_argument(2, model_directory)
        call get_command_argument(3, xsf)
        call get_command_argument(4, argument); read(argument, *) repeats
        if (repeats < 1) error stop "REPEATS must be positive"
        call load_predictor_from_n2p2(trim(model_directory), model)
        if (command_argument_count() == 5) then
            call get_command_argument(5, argument)
            select case(trim(argument))
            case("auto");   g5_mode = G5_EVALUATION_AUTO
            case("direct"); g5_mode = G5_EVALUATION_DIRECT
            case("moment"); g5_mode = G5_EVALUATION_MOMENT
            case default; error stop "G5_MODE must be auto, direct, or moment"
            end select
        end if
    else if (command_argument_count() >= 5) then
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
    else
        call usage()
    end if
    call model%set_g5_evaluation(g5_mode)
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
    select case(g5_mode)
    case(G5_EVALUATION_AUTO);   write(*, "(A)") "G5_MODE auto"
    case(G5_EVALUATION_DIRECT); write(*, "(A)") "G5_MODE direct"
    case(G5_EVALUATION_MOMENT); write(*, "(A)") "G5_MODE moment"
    end select
    write(*, "(A,1X,ES24.16)") "ENERGY_API_SECONDS_PER_STRUCTURE", energy_api_time/repeats
    write(*, "(A,1X,ES24.16)") "ENERGY_FORCE_API_SECONDS_PER_STRUCTURE", force_api_time/repeats
    write(*, "(A,1X,ES24.16)") "ENERGY_FILE_SECONDS_PER_STRUCTURE", file_time/repeats

contains

    subroutine usage()
        write(*, "(A)") "usage: accelnet-predict-benchmark NSPECIES SETUP... NETWORK... XSF REPEATS"
        write(*, "(A)") "   or: accelnet-predict-benchmark --n2p2 MODEL_DIR XSF REPEATS [G5_MODE]"
        error stop 2
    end subroutine usage
end program predict_benchmark
