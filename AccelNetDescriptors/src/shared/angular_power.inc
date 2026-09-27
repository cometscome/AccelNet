! Shared CPU/device scalar formula.
        real(real64), intent(in) :: cosine, lambda, zeta, derivative_prefactor
        integer, intent(in) :: integer_zeta
        real(real64), intent(out) :: value, derivative
        real(real64) :: base
        base = 0.5_real64*(1.0_real64 + lambda*cosine)
        if (integer_zeta > 0) then
            select case(integer_zeta)
            case(1)
                value = base
                derivative = derivative_prefactor
            case(2)
                value = base*base
                derivative = derivative_prefactor*base
            case(3)
                value = base*base*base
                derivative = derivative_prefactor*base*base
            case(4)
                value = (base*base)*(base*base)
                derivative = derivative_prefactor*base*base*base
            case default
                value = base**integer_zeta
                derivative = derivative_prefactor*base**(integer_zeta - 1)
            end select
        else
            value = base**zeta
            derivative = derivative_prefactor*base**(zeta - 1.0_real64)
        end if
