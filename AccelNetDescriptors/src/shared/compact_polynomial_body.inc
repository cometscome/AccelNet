        if (subtype > 10) then
            chain=chain*2*(1-t); t=t*(2-t)
        end if
        if (t >= 1) return ! Avoid division by zero after endpoint rounding.
        select case(mod(subtype,10))
        case(1)
            value=1+t*t*(2*t-3); derivative=6*t*(t-1)
        case(2)
            value=1+t**3*(-10+t*(15-6*t)); derivative=-30*t*t*(1-t)**2
        case(3)
            value=1+t**4*(-35+t*(84+t*(-70+20*t))); derivative=-140*t**3*(1-t)**3
        case(4)
            value=1+t**5*(-126+t*(420+t*(-540+t*(315-70*t)))); derivative=-630*t**4*(1-t)**4
        case(5)
            value=exp(1+1/(t*t-1)); derivative=-2*t*value/(t*t-1)**2
        end select
        derivative=derivative*chain
