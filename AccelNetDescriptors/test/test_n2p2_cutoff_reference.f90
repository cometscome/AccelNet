program test_n2p2_cutoff_reference
    use iso_c_binding, only: c_double, c_int
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: cutoff_value, cutoff_derivative, &
        CUTOFF_HARD, CUTOFF_FRACTIONAL
    implicit none

    interface
        subroutine n2p2_cutoff_value_derivative(kind, radius, alpha, distance, &
                                                 value, derivative) bind(C)
            import :: c_double, c_int
            integer(c_int), value :: kind
            real(c_double), value :: radius, alpha, distance
            real(c_double), intent(out) :: value, derivative
        end subroutine n2p2_cutoff_value_derivative
    end interface

    real(real64), parameter :: radius = 5.0_real64
    real(real64), parameter :: alphas(2) = [0.0_real64, 0.25_real64]
    real(real64), parameter :: distances(7) = [0.0_real64, 0.75_real64, 1.25_real64, &
        2.3_real64, 4.999_real64, 5.0_real64, 5.1_real64]
    real(real64) :: reference_value, reference_derivative, actual_value, actual_derivative
    real(real64) :: x, width, fractional_value, fractional_derivative
    integer :: kind, ia, ir

    do kind = CUTOFF_HARD, 8
        do ia = 1, size(alphas)
            do ir = 1, size(distances)
                call n2p2_cutoff_value_derivative(int(kind, c_int), radius, alphas(ia), &
                    distances(ir), reference_value, reference_derivative)
                actual_value = cutoff_value(distances(ir), radius, kind, alphas(ia))
                actual_derivative = cutoff_derivative(distances(ir), radius, kind, alphas(ia))
                call assert_close("n2p2 cutoff value", kind, alphas(ia), distances(ir), &
                    actual_value, reference_value)
                call assert_close("n2p2 cutoff derivative", kind, alphas(ia), distances(ir), &
                    actual_derivative, reference_derivative)
            end do
        end do
    end do

    ! Type 9 is an AccelNet extension, so compare it with its published
    ! X^2/(1+X^2) definition instead of pretending n2p2 has a reference.
    width = 0.2_real64*radius
    do ir = 1, size(distances)
        if (distances(ir) >= radius) then
            fractional_value = 0.0_real64
            fractional_derivative = 0.0_real64
        else
            x = (distances(ir) - radius)/width
            fractional_value = x*x/(1.0_real64 + x*x)
            fractional_derivative = 2.0_real64*x/(width*(1.0_real64 + x*x)**2)
        end if
        call assert_close("fractional cutoff value", CUTOFF_FRACTIONAL, 0.2_real64, &
            distances(ir), cutoff_value(distances(ir), radius, CUTOFF_FRACTIONAL, 0.2_real64), &
            fractional_value)
        call assert_close("fractional cutoff derivative", CUTOFF_FRACTIONAL, 0.2_real64, &
            distances(ir), cutoff_derivative(distances(ir), radius, CUTOFF_FRACTIONAL, 0.2_real64), &
            fractional_derivative)
    end do

contains

    subroutine assert_close(label, kind, alpha, distance, actual, reference)
        character(len=*), intent(in) :: label
        integer, intent(in) :: kind
        real(real64), intent(in) :: alpha, distance, actual, reference
        real(real64) :: scale
        scale = max(1.0_real64, abs(actual), abs(reference))
        if (abs(actual - reference) > 2.0e-13_real64*scale) then
            write(*, "(A,1X,I0,4(1X,ES24.16))") trim(label), kind, alpha, distance, actual, reference
            error stop "cutoff reference mismatch"
        end if
    end subroutine assert_close
end program test_n2p2_cutoff_reference
