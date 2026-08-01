program test_behler
    use iso_fortran_env, only: real64
    use accelnet_behler
    use accelnet_descriptor_models
    implicit none

    type(behler_config) :: config, g4_config, g5_config, moment_config
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
    real(real64) :: contracted_center(3), contracted_neighbors(3, 3)
    real(real64) :: contraction_coefficients(2)
    real(real64) :: moment_displacements(3, 20), moment_values(4), moment_direct_values(4)
    real(real64) :: moment_center(3, 4), moment_neighbors(3, 4, 20)
    real(real64) :: moment_contracted_center(3), moment_contracted_neighbors(3, 20)
    real(real64) :: moment_coefficients(4), radius, azimuth, polar
    integer :: moment_species(20)
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
    call add_g4(g4_config, 1, 2, 4.5_real64, 1.0_real64, 2.0_real64, 0.15_real64, 0.45_real64)
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
    call add_g5(g5_config, 1, 2, 4.5_real64, -1.0_real64, 3.0_real64, 0.2_real64, 0.35_real64)
    call add_g5(g5_config, 1, 2, 4.5_real64, -1.0_real64, 3.0_real64, 0.4_real64, 0.70_real64)
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
    contraction_coefficients = [0.7_real64, -0.3_real64]
    call contract_behler_derivatives(g5_config, displacements, species, contraction_coefficients, &
        contracted_center, contracted_neighbors)
    if (maxval(abs(contracted_center - matmul(g5_center, contraction_coefficients))) > 2.0e-13_real64) &
        error stop "G5 direct center contraction failed"
    do neighbor = 1, size(species)
        if (maxval(abs(contracted_neighbors(:, neighbor) - &
            matmul(g5_neighbors(:, :, neighbor), contraction_coefficients))) > 2.0e-13_real64) &
            error stop "G5 direct neighbor contraction failed"
    end do

    call initialize_behler_config(moment_config, 2)
    call add_g5(moment_config, 1, 1, 4.5_real64, 1.0_real64, 1.0_real64, 0.15_real64, 0.20_real64)
    call add_g5(moment_config, 1, 2, 4.5_real64, -1.0_real64, 2.0_real64, 0.25_real64, 0.55_real64)
    call add_g5(moment_config, 1, 2, 4.0_real64, 1.0_real64, 3.0_real64, 0.35_real64, 0.80_real64)
    call add_g5(moment_config, 2, 2, 4.0_real64, -1.0_real64, 4.0_real64, 0.45_real64, 1.10_real64)
    call set_behler_g5_evaluation(moment_config, G5_EVALUATION_MOMENT)
    do neighbor = 1, 20
        radius = 1.0_real64 + 0.11_real64*real(mod(neighbor, 17), real64)
        azimuth = 0.73_real64*real(neighbor, real64)
        polar = 0.35_real64 + 0.09_real64*real(mod(3*neighbor, 23), real64)
        moment_displacements(:, neighbor) = radius*[sin(polar)*cos(azimuth), &
            sin(polar)*sin(azimuth), cos(polar)]
        moment_species(neighbor) = 1 + mod(neighbor, 2)
    end do
    call evaluate_behler_values(moment_config, moment_displacements, moment_species, moment_values)
    call evaluate_behler_values_derivatives(moment_config, moment_displacements, moment_species, &
        moment_direct_values, moment_center, moment_neighbors)
    if (maxval(abs(moment_values - moment_direct_values)) > 2.0e-12_real64) then
        write(*, "(A,4(ES16.8,1X))") "G5 moment value error: ", moment_values - moment_direct_values
        error stop "G5 moment values differ from direct pairs"
    end if
    moment_coefficients = [0.7_real64, -0.3_real64, 0.2_real64, -0.5_real64]
    call contract_behler_derivatives(moment_config, moment_displacements, moment_species, &
        moment_coefficients, moment_contracted_center, moment_contracted_neighbors)
    if (maxval(abs(moment_contracted_center - matmul(moment_center, moment_coefficients))) > 2.0e-11_real64) &
        error stop "G5 moment center contraction differs from direct pairs"
    do neighbor = 1, 20
        if (maxval(abs(moment_contracted_neighbors(:, neighbor) - &
            matmul(moment_neighbors(:, :, neighbor), moment_coefficients))) > 2.0e-11_real64) &
            error stop "G5 moment neighbor contraction differs from direct pairs"
    end do

    call add_behler(model, config)
    allocate(model_values(model%num_descriptors()))
    call evaluate_model_values(model, displacements, species, model_values)
    if (maxval(abs(model_values - values)) > 0.0_real64) error stop "Behler model dispatch differs"
    if (abs(model%maximum_cutoff - 4.5_real64) > 0.0_real64) error stop "Behler model cutoff"
    write(*, "(A)") "Behler G1-G5 tests passed"
end program test_behler
