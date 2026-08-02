program benchmark_g5_scaling
    use iso_fortran_env, only: real64
    use accelnet_behler, only: behler_config, initialize_behler_config, add_g5, &
        evaluate_behler_values, contract_behler_derivatives, set_behler_g5_evaluation, &
        G5_EVALUATION_DIRECT, G5_EVALUATION_MOMENT_FORCE
    implicit none

    integer, parameter :: neighbor_counts(*) = [4, 8, 12, 16, 24, 32, 48, 64, 96, 128]
    real(real64), parameter :: eta_values(*) = [0.000357_real64, 0.028569_real64, 0.089277_real64]
    integer, parameter :: zeta_values(*) = [1, 2, 4]
    type(behler_config) :: config
    real(real64), allocatable :: displacements(:, :), values(:), direct_values(:), coefficients(:)
    real(real64), allocatable :: contracted_neighbors(:, :), direct_neighbors(:, :)
    integer, allocatable :: species(:)
    real(real64) :: contracted_center(3), direct_center(3)
    real(real64) :: start_time, end_time, direct_value_time, moment_value_time
    real(real64) :: direct_force_time, moment_force_time, radius, cosine_polar, sine_polar, azimuth
    real(real64) :: value_error, force_error
    integer :: repeats, count_index, nneighbors, neighbor, eta_index, zeta_index, pair, iteration
    character(len=64) :: argument

    repeats = 500
    if (command_argument_count() > 1) error stop "usage: benchmark-g5-scaling [REPEATS]"
    if (command_argument_count() == 1) then
        call get_command_argument(1, argument)
        read(argument, *) repeats
        if (repeats < 1) error stop "REPEATS must be positive"
    end if

    call initialize_behler_config(config, 2)
    do eta_index = 1, size(eta_values)
        do zeta_index = 1, size(zeta_values)
            do pair = 1, 3
                select case(pair)
                case(1)
                    call add_pair(1, 1)
                case(2)
                    call add_pair(1, 2)
                case(3)
                    call add_pair(2, 2)
                end select
            end do
        end do
    end do

    allocate(displacements(3, maxval(neighbor_counts)), species(maxval(neighbor_counts)))
    allocate(values(config%num_descriptors()), direct_values(config%num_descriptors()), &
             coefficients(config%num_descriptors()))
    allocate(contracted_neighbors(3, maxval(neighbor_counts)), direct_neighbors(3, maxval(neighbor_counts)))
    coefficients = [(sin(0.37_real64*real(neighbor, real64)), neighbor=1, size(coefficients))]

    write(*, "(A,1X,I0)") "DESCRIPTORS", config%num_descriptors()
    write(*, "(A,1X,I0)") "REPEATS", repeats
    write(*, "(A)") "NEIGHBORS DIRECT_VALUE MOMENT_VALUE VALUE_SPEEDUP DIRECT_FORCE MOMENT_FORCE FORCE_SPEEDUP VALUE_ERROR FORCE_ERROR"
    do count_index = 1, size(neighbor_counts)
        nneighbors = neighbor_counts(count_index)
        call make_environment(nneighbors)

        call set_behler_g5_evaluation(config, G5_EVALUATION_DIRECT)
        call evaluate_behler_values(config, displacements(:, 1:nneighbors), species(1:nneighbors), values)
        call contract_behler_derivatives(config, displacements(:, 1:nneighbors), species(1:nneighbors), &
            coefficients, contracted_center, contracted_neighbors(:, 1:nneighbors))
        direct_values = values
        direct_center = contracted_center
        direct_neighbors(:, 1:nneighbors) = contracted_neighbors(:, 1:nneighbors)
        call cpu_time(start_time)
        do iteration = 1, repeats
            call evaluate_behler_values(config, displacements(:, 1:nneighbors), species(1:nneighbors), values)
        end do
        call cpu_time(end_time)
        direct_value_time = (end_time - start_time)/real(repeats, real64)
        call cpu_time(start_time)
        do iteration = 1, repeats
            call contract_behler_derivatives(config, displacements(:, 1:nneighbors), species(1:nneighbors), &
                coefficients, contracted_center, contracted_neighbors(:, 1:nneighbors))
        end do
        call cpu_time(end_time)
        direct_force_time = (end_time - start_time)/real(repeats, real64)

        call set_behler_g5_evaluation(config, G5_EVALUATION_MOMENT_FORCE)
        call evaluate_behler_values(config, displacements(:, 1:nneighbors), species(1:nneighbors), values)
        call contract_behler_derivatives(config, displacements(:, 1:nneighbors), species(1:nneighbors), &
            coefficients, contracted_center, contracted_neighbors(:, 1:nneighbors))
        value_error = maxval(abs(values - direct_values))
        force_error = max(maxval(abs(contracted_center - direct_center)), &
            maxval(abs(contracted_neighbors(:, 1:nneighbors) - direct_neighbors(:, 1:nneighbors))))
        call cpu_time(start_time)
        do iteration = 1, repeats
            call evaluate_behler_values(config, displacements(:, 1:nneighbors), species(1:nneighbors), values)
        end do
        call cpu_time(end_time)
        moment_value_time = (end_time - start_time)/real(repeats, real64)
        call cpu_time(start_time)
        do iteration = 1, repeats
            call contract_behler_derivatives(config, displacements(:, 1:nneighbors), species(1:nneighbors), &
                coefficients, contracted_center, contracted_neighbors(:, 1:nneighbors))
        end do
        call cpu_time(end_time)
        moment_force_time = (end_time - start_time)/real(repeats, real64)

        write(*, "(I9,8(1X,ES14.6))") nneighbors, direct_value_time, moment_value_time, &
            direct_value_time/moment_value_time, direct_force_time, moment_force_time, &
            direct_force_time/moment_force_time, value_error, force_error
    end do

contains

    subroutine add_pair(species1, species2)
        integer, intent(in) :: species1, species2
        call add_g5(config, species1, species2, 6.5_real64, -1.0_real64, &
            real(zeta_values(zeta_index), real64), eta_values(eta_index), 0.35_real64)
        call add_g5(config, species1, species2, 6.5_real64, 1.0_real64, &
            real(zeta_values(zeta_index), real64), eta_values(eta_index), 0.35_real64)
    end subroutine add_pair

    subroutine make_environment(count)
        integer, intent(in) :: count
        integer :: index
        real(real64), parameter :: golden_angle = 2.39996322972865332_real64
        do index = 1, count
            radius = 1.5_real64 + 3.5_real64*real(mod(7*index, 31), real64)/31.0_real64
            cosine_polar = 1.0_real64 - 2.0_real64*(real(index, real64) - 0.5_real64)/real(count, real64)
            sine_polar = sqrt(max(0.0_real64, 1.0_real64 - cosine_polar*cosine_polar))
            azimuth = golden_angle*real(index, real64)
            displacements(:, index) = radius*[sine_polar*cos(azimuth), sine_polar*sin(azimuth), cosine_polar]
            species(index) = 1 + mod(index, 2)
        end do
    end subroutine make_environment
end program benchmark_g5_scaling
