        integer, intent(in) :: kind
        real(real64), intent(in) :: p(7),c,theta,inv_sin
        real(real64), intent(out) :: a,da
        real(real64), intent(in), optional :: angular_cache(2)
        real(real64) :: t,power,y,chain,width
        integer :: zi,subtype
        a=0; da=0
        if (kind == 13) then
            ! Reuse the established G4/G5 small-integer multiplication paths.
            ! Never round a near-integer exponent into an integer one.
            zi=0
            if (p(5) >= 1 .and. p(5) <= 16) then
                if (p(5) == real(nint(p(5)),real64)) zi=nint(p(5))
            end if
            if (zi > 0) then
                call angular_power(c,p(4),p(5),zi,0.5_real64*p(4)*p(5),a,da)
                a=2*a; da=2*da
            else
                t=max(0.0_real64,0.5_real64*(1+p(4)*c))
                power=t**(p(5)-1)
                a=2*t*power; da=p(4)*p(5)*power
            end if
        else
            ! Match n2p2's exclusion of exactly collinear compact angular pairs.
            if (c <= -1 .or. c >= 1) return
            if (theta <= p(3) .or. theta >= p(4)) return
            if (present(angular_cache)) then
                y=(theta-angular_cache(2))*angular_cache(1)
                t=abs(y); chain=sign(angular_cache(1),y)
            else
                width=0.5_real64*(p(4)-p(3))
                y=(theta-0.5_real64*(p(4)+p(3)))/width
                t=abs(y); chain=sign(1.0_real64,y)/width
            end if
            subtype=mod(nint(p(5)),10)
            associate(value=>a,derivative=>da)
                include 'compact_polynomial_body.inc'
            end associate
            da=-da*inv_sin
        end if
