! Shared CPU/device scalar formula.
        real(real64), intent(in) :: cosine, lambda, zeta
        integer, intent(in) :: integer_zeta
        real(real64) :: value, base
        base = 0.5_real64*(1.0_real64 + lambda*cosine)
        if (integer_zeta > 0) then
            select case(integer_zeta)
            case(1)
                value = base
            case(2)
                value = base*base
            case(3)
                value = base*base*base
            case(4)
                value = (base*base)*(base*base)
            case default
                value = base**integer_zeta
            end select
        else
            value = base**zeta
        end if
