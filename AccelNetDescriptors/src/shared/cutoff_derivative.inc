! Shared CPU/device scalar formula.
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
