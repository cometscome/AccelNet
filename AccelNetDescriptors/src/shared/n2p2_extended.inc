! Scalar n2p2 extensions, compiled identically for host and device.
! Compact p = [Rc, Rlow, theta_left(rad), theta_right(rad), subtype, 0, 0].
! Weighted exponential p = [Rc, eta, Rs, lambda, zeta, 0, cutoff_alpha].
    pure subroutine compact_window(x,left,right,subtype,value,derivative)
        !$omp declare target
        real(real64), intent(in) :: x,left,right
        integer, intent(in) :: subtype
        real(real64), intent(out) :: value,derivative
        real(real64) :: y,t,chain,width
        value=0; derivative=0
        if (x <= left .or. x >= right) return
        width=0.5_real64*(right-left)
        y=(x-0.5_real64*(right+left))/width
        t=abs(y); chain=sign(1.0_real64,y)/width
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
    end subroutine

    pure subroutine extended_radial(kind,ct,p,r,v,dv)
        !$omp declare target
        integer, intent(in) :: kind,ct
        real(real64), intent(in) :: p(7),r
        real(real64), intent(out) :: v,dv
        real(real64) :: e,fc,dfc
        v=0; dv=0
        if (r <= 1e-12_real64 .or. r >= p(1)) return
        if (kind >= 20) then
            call compact_window(r,p(2),p(1),nint(p(5)),v,dv)
        else
            e=exp(-p(2)*(r-p(3))**2)
            fc=sf_cutoff_value(r,p(1),ct,p(7)); dfc=sf_cutoff_derivative(r,p(1),ct,p(7))
            v=e*fc; dv=e*(dfc-2*p(2)*(r-p(3))*fc)
        end if
    end subroutine

    pure subroutine extended_pair(kind,ct,p,uj,uk,rj,rk,qj,qk,dqj,dqk,v,gj,gk)
        !$omp declare target
        integer, intent(in) :: kind,ct
        real(real64), intent(in) :: p(7),uj(3),uk(3),rj,rk,qj,qk,dqj,dqk
        real(real64), intent(out) :: v,gj(3),gk(3)
        real(real64) :: c,theta,a,da,t,ujk(3),rjk,qjk,dqjk,prod
        v=0; gj=0; gk=0
        if (rj <= 1e-12_real64 .or. rk <= 1e-12_real64) return
        if (qj == 0 .or. qk == 0) return
        c=max(-1.0_real64,min(1.0_real64,sum(uj*uk)))
        qjk=1; dqjk=0; ujk=0
        if (kind == 13 .or. kind == 21 .or. kind == 24) then
            ujk=rk*uk-rj*uj; rjk=sqrt(sum(ujk**2))
            if (rjk <= 1e-12_real64) return
            ujk=ujk/rjk
            call extended_radial(kind,ct,p,rjk,qjk,dqjk)
            if (qjk == 0) return
        end if
        if (kind == 13) then
            t=max(0.0_real64,0.5_real64*(1+p(4)*c))
            a=2*t**p(5)
            if (p(5) == 1) then
                da=p(4)
            else
                da=p(4)*p(5)*t**(p(5)-1)
            end if
        else
            ! Match n2p2's exclusion of exactly collinear compact angular pairs.
            if (c <= -1 .or. c >= 1) return
            theta=acos(c)
            call compact_window(theta,p(3),p(4),mod(nint(p(5)),10),a,da)
            da=-da/sqrt(1-c*c)
        end if
        prod=qj*qk*qjk; v=a*prod
        gj=da*prod*(uk-c*uj)/rj+a*qk*(dqj*qjk*uj-qj*dqjk*ujk)
        gk=da*prod*(uj-c*uk)/rk+a*qj*(dqk*qjk*uk+qk*dqjk*ujk)
    end subroutine
