! Device scalar formulas mirror the CPU implementation; CPU objects are unchanged.
module accelnet_target_math
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: CUTOFF_HARD, CUTOFF_COS, CUTOFF_TANHU, CUTOFF_TANH, &
        CUTOFF_EXP, CUTOFF_POLY1, CUTOFF_POLY2, CUTOFF_POLY3, CUTOFF_POLY4, CUTOFF_FRACTIONAL
    use aenet_network, only: ACTIVATION_LINEAR, ACTIVATION_TANH, ACTIVATION_LOGISTIC, &
        ACTIVATION_AENET_MTANH, ACTIVATION_AENET_TWIST, ACTIVATION_SOFTPLUS, ACTIVATION_RELU, &
        ACTIVATION_GAUSSIAN, ACTIVATION_COSINE, ACTIVATION_REVERSE_LOGISTIC, ACTIVATION_EXPONENTIAL, ACTIVATION_HARMONIC
    implicit none
    private
    real(real64), parameter :: PI_ACCELNET = 3.14159265358979_real64
    public :: target_cutoff_value, target_cutoff_derivative, target_activate, target_activation_derivative
contains
    pure real(real64) function target_cutoff_value(distance, rc, cutoff_type, alpha) result(value)
        !$omp declare target
        real(real64), intent(in) :: distance, rc
        integer, intent(in), optional :: cutoff_type
        real(real64), intent(in), optional :: alpha
        integer :: kind
        real(real64) :: inner, inverse_width, x, t, width
        kind = CUTOFF_COS
        if (present(cutoff_type)) kind = cutoff_type
        if (kind == CUTOFF_HARD) then
            value = merge(1.0_real64, 0.0_real64, distance < rc)
            return
        end if
        if (distance >= rc) then
            value = 0.0_real64
            return
        end if
        select case(kind)
        case(CUTOFF_TANHU, CUTOFF_TANH)
            t = tanh(1.0_real64 - distance/rc)
            value = t*t*t
            if (kind == CUTOFF_TANH) value = value/tanh(1.0_real64)**3
        case(CUTOFF_FRACTIONAL)
            width = rc
            if (present(alpha)) width = alpha*rc
            x = (distance - rc)/width
            value = x*x/(1.0_real64 + x*x)
        case default
            inner = 0.0_real64
            if (present(alpha)) inner = alpha*rc
            if (distance < inner) then
                value = 1.0_real64
                return
            end if
            inverse_width = 1.0_real64/(rc - inner)
            x = (distance - inner)*inverse_width
            select case(kind)
            case(CUTOFF_COS)
                value = 0.5_real64*(cos(PI_ACCELNET*x) + 1.0_real64)
            case(CUTOFF_EXP)
                value = exp(1.0_real64 + 1.0_real64/(x*x - 1.0_real64))
            case(CUTOFF_POLY1)
                value = (2.0_real64*x - 3.0_real64)*x*x + 1.0_real64
            case(CUTOFF_POLY2)
                value = ((15.0_real64 - 6.0_real64*x)*x - 10.0_real64)*x*x*x + 1.0_real64
            case(CUTOFF_POLY3)
                value = (x*(x*(20.0_real64*x - 70.0_real64) + 84.0_real64) - 35.0_real64)*x**4 + 1.0_real64
            case(CUTOFF_POLY4)
                value = (x*(x*((315.0_real64 - 70.0_real64*x)*x - 540.0_real64) + &
                    420.0_real64) - 126.0_real64)*x**5 + 1.0_real64
            case default
                value = 0.0_real64
            end select
        end select
    end function target_cutoff_value

    pure real(real64) function target_cutoff_derivative(distance, rc, cutoff_type, alpha) result(value)
        !$omp declare target
        real(real64), intent(in) :: distance, rc
        integer, intent(in), optional :: cutoff_type
        real(real64), intent(in), optional :: alpha
        integer :: kind
        real(real64) :: inner, inverse_width, x, t, core_derivative, width
        kind = CUTOFF_COS
        if (present(cutoff_type)) kind = cutoff_type
        if (distance >= rc) then
            value = 0.0_real64
            return
        end if
        select case(kind)
        case(CUTOFF_HARD)
            value = 0.0_real64
        case(CUTOFF_TANHU, CUTOFF_TANH)
            t = tanh(1.0_real64 - distance/rc)
            value = 3.0_real64*t*t*(t*t - 1.0_real64)/rc
            if (kind == CUTOFF_TANH) value = value/tanh(1.0_real64)**3
        case(CUTOFF_FRACTIONAL)
            width = rc
            if (present(alpha)) width = alpha*rc
            x = (distance - rc)/width
            value = 2.0_real64*x/(width*(1.0_real64 + x*x)**2)
        case default
            inner = 0.0_real64
            if (present(alpha)) inner = alpha*rc
            if (distance < inner) then
                value = 0.0_real64
                return
            end if
            inverse_width = 1.0_real64/(rc - inner)
            x = (distance - inner)*inverse_width
            select case(kind)
            case(CUTOFF_COS)
                core_derivative = -0.5_real64*PI_ACCELNET*sin(PI_ACCELNET*x)
            case(CUTOFF_EXP)
                t = 1.0_real64/(x*x - 1.0_real64)
                core_derivative = -2.0_real64*x*t*t*exp(1.0_real64 + t)
            case(CUTOFF_POLY1)
                core_derivative = x*(6.0_real64*x - 6.0_real64)
            case(CUTOFF_POLY2)
                core_derivative = x*x*((60.0_real64 - 30.0_real64*x)*x - 30.0_real64)
            case(CUTOFF_POLY3)
                core_derivative = x**3*(x*(x*(140.0_real64*x - 420.0_real64) + &
                    420.0_real64) - 140.0_real64)
            case(CUTOFF_POLY4)
                core_derivative = x**4*(x*(x*((2520.0_real64 - 630.0_real64*x)*x - &
                    3780.0_real64) + 2520.0_real64) - 630.0_real64)
            case default
                core_derivative = 0.0_real64
            end select
            value = inverse_width*core_derivative
        end select
    end function target_cutoff_derivative

    pure real(real64) function target_activate(x, code) result(y)
        !$omp declare target
        real(real64), intent(in) :: x
        integer, intent(in) :: code
        real(real64), parameter :: a = 1.7159_real64
        real(real64), parameter :: b = 0.666666666666667_real64
        real(real64), parameter :: c = 0.1_real64
        select case(code)
        case(ACTIVATION_LINEAR); y = x
        case(ACTIVATION_TANH); y = tanh(x)
        case(ACTIVATION_LOGISTIC); y = 1.0_real64/(1.0_real64 + exp(-x))
        case(ACTIVATION_AENET_MTANH)
            y = a*tanh(b*x)
        case(ACTIVATION_AENET_TWIST)
            y = a*tanh(b*x) + c*x
        case(ACTIVATION_SOFTPLUS)
            if (x > 0.0_real64) then
                y = x + log(1.0_real64 + exp(-x))
            else
                y = log(1.0_real64 + exp(x))
            end if
        case(ACTIVATION_RELU); y = max(x, 0.0_real64)
        case(ACTIVATION_GAUSSIAN); y = exp(-0.5_real64*x*x)
        case(ACTIVATION_COSINE); y = cos(x)
        case(ACTIVATION_REVERSE_LOGISTIC); y = 1.0_real64 - 1.0_real64/(1.0_real64 + exp(-x))
        case(ACTIVATION_EXPONENTIAL); y = exp(-x)
        case(ACTIVATION_HARMONIC); y = x*x
        case default; y = 0.0_real64 ! rejected during model packing
        end select
    end function target_activate

    pure real(real64) function target_activation_derivative(x, y, code) result(value)
        !$omp declare target
        real(real64), intent(in) :: x, y
        integer, intent(in) :: code
        real(real64), parameter :: a = 1.7159_real64
        real(real64), parameter :: b = 0.666666666666667_real64
        real(real64), parameter :: c = 0.1_real64
        real(real64) :: tanhbx
        select case(code)
        case(ACTIVATION_LINEAR); value = 1.0_real64
        case(ACTIVATION_TANH); value = 1.0_real64 - y*y
        case(ACTIVATION_LOGISTIC); value = y*(1.0_real64 - y)
        case(ACTIVATION_AENET_MTANH)
            tanhbx = tanh(b*x)
            value = a*b*(1.0_real64 - tanhbx*tanhbx)
        case(ACTIVATION_AENET_TWIST)
            tanhbx = tanh(b*x)
            value = a*b*(1.0_real64 - tanhbx*tanhbx) + c
        case(ACTIVATION_SOFTPLUS); value = 1.0_real64 - exp(-y)
        case(ACTIVATION_RELU); value = merge(1.0_real64, 0.0_real64, y > 0.0_real64)
        case(ACTIVATION_GAUSSIAN); value = -x*y
        case(ACTIVATION_COSINE); value = -sin(x)
        case(ACTIVATION_REVERSE_LOGISTIC); value = -y*(1.0_real64 - y)
        case(ACTIVATION_EXPONENTIAL); value = -y
        case(ACTIVATION_HARMONIC); value = 2.0_real64*x
        case default; value = 0.0_real64 ! rejected during model packing
        end select
    end function target_activation_derivative
end module
