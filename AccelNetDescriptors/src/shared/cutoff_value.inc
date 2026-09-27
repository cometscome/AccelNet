! Shared CPU/device scalar formula.
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
