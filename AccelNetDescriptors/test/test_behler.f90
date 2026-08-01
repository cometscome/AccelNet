program test_behler
    use iso_fortran_env, only: real64
    use accelnet_behler
    use accelnet_descriptor_models
    implicit none

    type(behler_config) :: config, g4_config, g5_config
    type(descriptor_model) :: model
    real(real64) :: displacements(3, 3), shifted(3, 3)
    integer :: species(3)
    real(real64), allocatable :: values(:), plus_values(:), minus_values(:)
    real(real64), allocatable :: center_derivative(:, :), neighbor_derivative(:, :, :)
    real(real64), allocatable :: model_values(:)
    real(real64) :: finite_difference, error
    real(real64) :: g4_values(1), g4_plus(1), g4_minus(1)
    real(real64) :: g4_center(3, 1), g4_neighbors(3, 1, 3)
    real(real64) :: g5_values(2), g5_plus(2), g5_minus(2)
    real(real64) :: g5_center(3, 2), g5_neighbors(3, 2, 3)
    integer :: neighbor, component, coefficient
    real(real64), parameter :: step = 1.0e-6_real64

    call initialize_behler_config(config, 2)
    call add_g1(config, 1, 4.5_real64)
    call add_g2(config, 2, 4.5_real64, 0.3_real64, 0.7_real64)
    call add_g3(config, 1, 4.5_real64, 1.2_real64)
    call add_g4(config, 1, 2, 4.5_real64, 1.0_real64, 2.0_real64, 0.15_real64)
    call add_g5(config, 1, 2, 4.5_real64, -1.0_real64, 3.0_real64, 0.2_real64)
    displacements(:, 1) = [1.1_real64, 0.2_real64, -0.1_real64]
    displacements(:, 2) = [-0.4_real64, 1.3_real64, 0.3_real64]
    displacements(:, 3) = [0.2_real64, -0.5_real64, 1.5_real64]
    species = [1, 2, 1]

    allocate(values(config%num_descriptors()), plus_values(config%num_descriptors()))
    allocate(minus_values(config%num_descriptors()))
    allocate(center_derivative(3, config%num_descriptors()))
    allocate(neighbor_derivative(3, config%num_descriptors(), size(species)))
    call evaluate_behler_values_derivatives(config, displacements, species, values, &
                                            center_derivative, neighbor_derivative)
    if (.not. all(values == values)) error stop "non-finite Behler values"
    error = maxval(abs(center_derivative + sum(neighbor_derivative, dim=3)))
    if (error > 2.0e-12_real64) error stop "Behler translation derivative failed"

    do neighbor = 1, size(species)
        do component = 1, 3
            shifted = displacements
            shifted(component, neighbor) = shifted(component, neighbor) + step
            call evaluate_behler_values(config, shifted, species, plus_values)
            shifted(component, neighbor) = shifted(component, neighbor) - 2.0_real64*step
            call evaluate_behler_values(config, shifted, species, minus_values)
            do coefficient = 1, config%num_descriptors()
                finite_difference = (plus_values(coefficient) - minus_values(coefficient))/(2.0_real64*step)
                error = abs(finite_difference - neighbor_derivative(component, coefficient, neighbor))
                if (error > 3.0e-8_real64*max(1.0_real64, abs(finite_difference))) then
                    write(*, "(A,3(I0,1X),3(ES16.8,1X))") "Behler derivative mismatch: ", &
                        neighbor, component, coefficient, finite_difference, &
                        neighbor_derivative(component, coefficient, neighbor), error
                    error stop "Behler finite difference failed"
                end if
            end do
        end do
    end do

    call initialize_behler_config(g4_config, 2)
    call add_g4(g4_config, 1, 2, 4.5_real64, 1.0_real64, 2.0_real64, 0.15_real64)
    call evaluate_behler_values_derivatives(g4_config, displacements, species, g4_values, &
                                            g4_center, g4_neighbors)
    do neighbor = 1, size(species)
        do component = 1, 3
            shifted = displacements
            shifted(component, neighbor) = shifted(component, neighbor) + step
            call evaluate_behler_values(g4_config, shifted, species, g4_plus)
            shifted(component, neighbor) = shifted(component, neighbor) - 2.0_real64*step
            call evaluate_behler_values(g4_config, shifted, species, g4_minus)
            finite_difference = (g4_plus(1) - g4_minus(1))/(2.0_real64*step)
            error = abs(finite_difference - g4_neighbors(component, 1, neighbor))
            if (error > 3.0e-8_real64*max(1.0_real64, abs(finite_difference))) &
                error stop "G4-only finite difference failed"
        end do
    end do

    call initialize_behler_config(g5_config, 2)
    call add_g5(g5_config, 1, 2, 4.5_real64, -1.0_real64, 3.0_real64, 0.2_real64)
    call add_g5(g5_config, 1, 2, 4.5_real64, -1.0_real64, 3.0_real64, 0.4_real64)
    call evaluate_behler_values_derivatives(g5_config, displacements, species, g5_values, &
                                            g5_center, g5_neighbors)
    do neighbor = 1, size(species)
        do component = 1, 3
            shifted = displacements
            shifted(component, neighbor) = shifted(component, neighbor) + step
            call evaluate_behler_values(g5_config, shifted, species, g5_plus)
            shifted(component, neighbor) = shifted(component, neighbor) - 2.0_real64*step
            call evaluate_behler_values(g5_config, shifted, species, g5_minus)
            do coefficient = 1, 2
                finite_difference = (g5_plus(coefficient) - g5_minus(coefficient))/(2.0_real64*step)
                error = abs(finite_difference - g5_neighbors(component, coefficient, neighbor))
                if (error > 3.0e-8_real64*max(1.0_real64, abs(finite_difference))) &
                    error stop "G5-only finite difference failed"
            end do
        end do
    end do

    call add_behler(model, config)
    allocate(model_values(model%num_descriptors()))
    call evaluate_model_values(model, displacements, species, model_values)
    if (maxval(abs(model_values - values)) > 0.0_real64) error stop "Behler model dispatch differs"
    if (abs(model%maximum_cutoff - 4.5_real64) > 0.0_real64) error stop "Behler model cutoff"
    write(*, "(A)") "Behler G1-G5 tests passed"
end program test_behler
