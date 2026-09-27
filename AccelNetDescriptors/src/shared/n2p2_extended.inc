    pure subroutine extended_angular(kind,p,c,theta,inv_sin,a,da,angular_cache)
        !$omp declare target
        include 'extended_angular_body.inc'
    end subroutine

! Scalar n2p2 extensions, compiled identically for host and device.
! Compact p = [Rc, Rlow, theta_left(rad), theta_right(rad), subtype, 0, 0].
! Weighted exponential p = [Rc, eta, Rs, lambda, zeta, 0, cutoff_alpha].
    pure subroutine compact_window(x,left,right,subtype,value,derivative,angular_cache)
        !$omp declare target
        include 'compact_window_body.inc'
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
        theta=0; t=0
        if (kind /= 13 .and. c > -1 .and. c < 1) then
            theta=acos(c); t=1/sqrt(1-c*c)
        end if
        call extended_angular(kind,p,c,theta,t,a,da)
        prod=qj*qk*qjk; v=a*prod
        gj=da*prod*(uk-c*uj)/rj+a*qk*(dqj*qjk*uj-qj*dqjk*ujk)
        gk=da*prod*(uj-c*uk)/rk+a*qj*(dqk*qjk*uk+qk*dqjk*ujk)
    end subroutine
