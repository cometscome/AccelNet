program test_lj_model
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: descriptor_config, initialize_config, evaluate_atom
    use accelnet_lj
    use accelnet_descriptor_models
    implicit none

    type(lj_config) :: lj
    type(descriptor_config) :: chebyshev
    type(descriptor_model) :: model
    real(real64) :: displacements(3, 2), lj_values(4), expected(4)
    real(real64) :: center_derivative(3, 4), neighbor_derivative(3, 4, 2)
    real(real64) :: plus_values(4), minus_values(4), shifted(3, 2), finite_difference
    real(real64), allocatable :: model_values(:), chebyshev_values(:)
    real(real64), allocatable :: model_center_derivative(:, :), model_neighbor_derivative(:, :, :)
    integer :: species(2), coefficient
    real(real64), parameter :: step = 1.0e-6_real64

    call initialize_lj_config(lj, 2, 3.0_real64)
    displacements = 0.0_real64
    displacements(1, 1) = 1.0_real64
    displacements(2, 2) = 2.0_real64
    species = [1, 2]
    expected = [1.0_real64, 1.0_real64, 1.0_real64/64.0_real64, 1.0_real64/4096.0_real64]

    call evaluate_lj_values(lj, displacements, species, lj_values)
    call assert_close(maxval(abs(lj_values - expected)), 0.0_real64, 1.0e-15_real64, "known LJ values")

    call evaluate_lj_values_derivatives(lj, displacements, species, lj_values, &
                                        center_derivative, neighbor_derivative)
    call assert_close(neighbor_derivative(1, 1, 1), -6.0_real64, 1.0e-14_real64, "r^-6 derivative")
    call assert_close(neighbor_derivative(1, 2, 1), -12.0_real64, 1.0e-14_real64, "r^-12 derivative")
    call assert_close(maxval(abs(center_derivative + sum(neighbor_derivative, dim=3))), &
                      0.0_real64, 1.0e-14_real64, "translation derivative")

    do coefficient = 1, lj%num_descriptors()
        shifted = displacements
        shifted(2, 2) = shifted(2, 2) + step
        call evaluate_lj_values(lj, shifted, species, plus_values)
        shifted(2, 2) = shifted(2, 2) - 2.0_real64*step
        call evaluate_lj_values(lj, shifted, species, minus_values)
        finite_difference = (plus_values(coefficient) - minus_values(coefficient))/(2.0_real64*step)
        call assert_close(neighbor_derivative(2, coefficient, 2), finite_difference, &
                          2.0e-9_real64, "LJ finite difference")
    end do

    call initialize_config(chebyshev, 2, 3.0_real64, 4, 3.0_real64, 2, version=0)
    call add_chebyshev(model, chebyshev)
    call add_lj(model, lj)
    allocate(model_values(model%num_descriptors()))
    allocate(chebyshev_values(chebyshev%num_descriptors()))
    allocate(model_center_derivative(3, model%num_descriptors()))
    allocate(model_neighbor_derivative(3, model%num_descriptors(), size(species)))
    call evaluate_model_values(model, displacements, species, model_values)
    call evaluate_atom(chebyshev, displacements, species, chebyshev_values)
    call assert_close(maxval(abs(model_values(1:size(chebyshev_values)) - chebyshev_values)), &
                      0.0_real64, 0.0_real64, "model Chebyshev fast path")
    call assert_close(maxval(abs(model_values(size(chebyshev_values) + 1:) - expected)), &
                      0.0_real64, 1.0e-15_real64, "model LJ offset")
    call evaluate_model_values_derivatives(model, displacements, species, model_values, &
                                           model_center_derivative, model_neighbor_derivative)
    call assert_close(maxval(abs(model_neighbor_derivative(:, size(chebyshev_values) + 1:, :) - &
                                  neighbor_derivative)), 0.0_real64, 0.0_real64, "model LJ derivatives")
    if (abs(model%maximum_cutoff - 3.0_real64) > 0.0_real64) error stop "model maximum cutoff"

    write(*, "(A)") "LJ and descriptor model tests passed"

contains

    subroutine assert_close(actual, reference, tolerance, label)
        real(real64), intent(in) :: actual, reference, tolerance
        character(len=*), intent(in) :: label
        if (abs(actual - reference) > tolerance) then
            write(*, "(A,2(1X,ES24.16))") trim(label), actual, reference
            error stop "assert_close failed"
        end if
    end subroutine assert_close

end program test_lj_model
