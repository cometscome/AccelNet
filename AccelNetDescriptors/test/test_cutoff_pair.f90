program test_cutoff_pair
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet_descriptors, only: cutoff_value, cutoff_derivative, cutoff_value_derivative
    implicit none
    integer :: kind, plateau, i
    real(real64) :: rc, alpha, r, value, derivative, numerical, points(9)
    real(real64), parameter :: step = 1.0e-6_real64

    rc = 4.5_real64
    do kind = 0, 9
        do plateau = 0, 1
            alpha = 0.25_real64*plateau
            if (kind == 9 .and. plateau == 0) alpha = 0.1_real64
            points = [0.0_real64, alpha*rc, 0.15_real64*rc, 0.43_real64*rc, 0.79_real64*rc, &
                rc*(1.0_real64 - 1.0e-5_real64), rc, rc + step, 2.0_real64*rc]
            do i = 1, size(points)
                r = points(i)
                call cutoff_value_derivative(r, rc, kind, alpha, value, derivative)
                if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(derivative)) &
                    error stop 'non-finite fused cutoff'
                if (abs(value - cutoff_value(r, rc, kind, alpha)) > 3.0e-14_real64) &
                    error stop 'fused cutoff value mismatch'
                if (abs(derivative - cutoff_derivative(r, rc, kind, alpha)) > 3.0e-14_real64) &
                    error stop 'fused cutoff derivative mismatch'
                if (kind == 0 .or. r <= step .or. abs(r - alpha*rc) <= step .or. abs(r - rc) <= step) cycle
                numerical = (cutoff_value(r + step, rc, kind, alpha) - &
                    cutoff_value(r - step, rc, kind, alpha))/(2.0_real64*step)
                if (abs(numerical - derivative) > 5.0e-8_real64) error stop 'cutoff finite difference mismatch'
            end do
        end do
    end do
    ! Preserve the optional argument defaults of the separate cutoff functions.
    do kind = 0, 9
        call cutoff_value_derivative(1.3_real64, rc, cutoff_type=kind, value=value, derivative=derivative)
        if (abs(value - cutoff_value(1.3_real64, rc, kind)) > 3.0e-14_real64) error stop 'default alpha value'
        if (abs(derivative - cutoff_derivative(1.3_real64, rc, kind)) > 3.0e-14_real64) &
            error stop 'default alpha derivative'
    end do
    call cutoff_value_derivative(1.3_real64, rc, value=value, derivative=derivative)
    if (abs(value - cutoff_value(1.3_real64, rc)) > 3.0e-14_real64) error stop 'default cutoff value'
    if (abs(derivative - cutoff_derivative(1.3_real64, rc)) > 3.0e-14_real64) error stop 'default cutoff derivative'
    print *, 'Fused cutoff values and derivatives passed for all ten cutoff types.'
end program test_cutoff_pair
