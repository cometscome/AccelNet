! Scalar formulas are shared with CPU sources; this module compiles their device variants.
module accelnet_target_math
    use iso_fortran_env, only: real64
    use accelnet_target_descriptors, only: PACKED_FEATURE_FIELDS
    use accelnet_descriptors, only: CUTOFF_HARD, CUTOFF_COS, CUTOFF_TANHU, CUTOFF_TANH, &
        CUTOFF_EXP, CUTOFF_POLY1, CUTOFF_POLY2, CUTOFF_POLY3, CUTOFF_POLY4, CUTOFF_FRACTIONAL
    use aenet_network, only: ACTIVATION_LINEAR, ACTIVATION_TANH, ACTIVATION_LOGISTIC, &
        ACTIVATION_AENET_MTANH, ACTIVATION_AENET_TWIST, ACTIVATION_SOFTPLUS, ACTIVATION_RELU, &
        ACTIVATION_GAUSSIAN, ACTIVATION_COSINE, ACTIVATION_REVERSE_LOGISTIC, ACTIVATION_EXPONENTIAL, ACTIVATION_HARMONIC
    implicit none
    private
    real(real64), parameter :: PI_ACCELNET = 3.14159265358979_real64
    public :: extended_radial, extended_pair, extended_angular
    public :: generic_pair_geometry, generic_pair_contracted, generic_pair_value
    public :: generic_radial, generic_radial_value, generic_lj, gaussian_radial
    public :: angular_power, angular_value, g5_moment_active
    public :: target_cutoff_value, target_cutoff_derivative, target_activate, target_activation_derivative
    public :: target_cutoff_pair
contains
    pure logical function g5_moment_active(mode,group,neighbors) result(active)
        !$omp declare target
        integer, intent(in) :: mode,group,neighbors
        ! Mode 2 keeps the established 16-neighbor threshold; 3 forces moments.
        active = group > 0 .and. (mode == 3 .or. ((mode == 0 .or. mode == 2) .and. neighbors >= 16))
    end function

    pure real(real64) function target_cutoff_value(distance, rc, cutoff_type, alpha) result(value)
        !$omp declare target
        include 'cutoff_value.inc'
    end function target_cutoff_value

    pure real(real64) function target_cutoff_derivative(distance, rc, cutoff_type, alpha) result(value)
        !$omp declare target
        include 'cutoff_derivative.inc'
    end function target_cutoff_derivative

    pure subroutine target_cutoff_pair(distance,rc,kind,alpha,value,derivative)
        !$omp declare target
        include 'cutoff_pair.inc'
    end subroutine

    pure real(real64) function target_activate(x, code) result(y)
        !$omp declare target
        real(real64), intent(in) :: x
        integer, intent(in) :: code
        real(real64), parameter :: a = 1.7159_real64
        real(real64), parameter :: b = 0.666666666666667_real64
        real(real64), parameter :: c = 0.1_real64
        select case(code)
        include 'activate.inc'
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
        include 'activation_derivative.inc'
        case default; value = 0.0_real64 ! rejected during model packing
        end select
    end function target_activation_derivative
    pure subroutine angular_power(cosine, lambda, zeta, integer_zeta, derivative_prefactor, value, derivative)
        !$omp declare target
        include 'angular_power.inc'
    end subroutine
    pure function angular_value(cosine, lambda, zeta, integer_zeta) result(value)
        !$omp declare target
        include 'angular_value.inc'
    end function angular_value
    ! Radial feature and derivative with respect to distance.
    pure subroutine gaussian_radial(r,eta,rs,fc,dfc,value,derivative)
        !$omp declare target
        real(real64), intent(in) :: r,eta,rs,fc,dfc
        real(real64), intent(out) :: value,derivative
        real(real64) :: h,dh
        h=exp(-eta*(r-rs)**2)
        dh=-2*eta*(r-rs)*h
        value=fc*h; derivative=dfc*h+fc*dh
    end subroutine

    pure subroutine generic_radial(kind, r, ct, p, value, derivative)
        !$omp declare target
        integer, intent(in) :: kind,ct
        real(real64), intent(in) :: r,p(7)
        real(real64), intent(out) :: value,derivative
        real(real64) :: fc,dfc,h,dh
        value = 0; derivative = 0
        if (r <= 1e-12_real64 .or. r > p(1)) return
        call target_cutoff_pair(r,p(1),ct,p(7),fc,dfc)
        h = 1; dh = 0
        select case(kind)
        case(2,4,5)
            call gaussian_radial(r,p(2),p(3),fc,dfc,value,derivative)
            return
        case(3)
            h = cos(p(6)*r); dh = -p(6)*sin(p(6)*r)
        case(6)
            h = r**(-6); dh = -6*h/r
        case(7)
            h = r**(-12); dh = -12*h/r
        end select
        value = fc*h; derivative = dfc*h+fc*dh
    end subroutine

    ! Value-only Behler radial terms; LJ uses the fused 6/12 helper.
    pure real(real64) function generic_radial_value(kind,r,ct,p) result(value)
        !$omp declare target
        integer, intent(in) :: kind,ct
        real(real64), intent(in) :: r,p(7)
        real(real64) :: h
        value = 0
        if (r <= 1e-12_real64 .or. r > p(1)) return
        h = 1
        select case(kind)
        case(2,4,5); h = exp(-p(2)*(r-p(3))**2)
        case(3); h = cos(p(6)*r)
        end select
        value = target_cutoff_value(r,p(1),ct,p(7))*h
    end function


    ! LJ columns are emitted in adjacent 6/12 pairs by the validated packer.
    pure subroutine generic_lj(r,ct,p,force,v,v12,dv,dv12)
        !$omp declare target
        integer, intent(in) :: ct
        real(real64), intent(in) :: r,p(7)
        logical, intent(in) :: force
        real(real64), intent(out) :: v,v12,dv,dv12
        real(real64) :: ir2,ir6,ir12,fc,dfc
        v=0; v12=0; dv=0; dv12=0
        if (r <= 1e-12_real64 .or. r > p(1)) return
        ir2=1.0_real64/(r*r); ir6=ir2*ir2*ir2; ir12=ir6*ir6
        fc=target_cutoff_value(r,p(1),ct,p(7))
        v=fc*ir6; v12=fc*ir12
        if (force) then
            dfc=target_cutoff_derivative(r,p(1),ct,p(7))
            dv=ir6*(dfc-6*fc/r); dv12=ir12*(dfc-12*fc/r)
        end if
    end subroutine

    pure subroutine generic_pair_geometry(f,p,uj,uk,rj,rk,cosine,qjk,dqjk,ujk,force)
        !$omp declare target
        integer, intent(in) :: f(PACKED_FEATURE_FIELDS)
        real(real64), intent(in) :: p(7),uj(3),uk(3),rj,rk
        logical, intent(in) :: force
        real(real64), intent(out) :: cosine,qjk,dqjk,ujk(3)
        real(real64) :: rjk,rjk2
        cosine=0; qjk=0; dqjk=0; ujk=0
        if (rj <= 1e-12_real64 .or. rk <= 1e-12_real64 .or. rj > p(1) .or. rk > p(1)) return
        qjk=1
        if (f(1) == 4) then
            ujk=rk*uk-rj*uj; rjk2=sum(ujk**2)
            if (rjk2 <= 1e-24_real64 .or. rjk2 >= p(1)*p(1)) then
                qjk=0; return
            end if
            rjk=sqrt(rjk2)
            if (force) then
                ujk=ujk/rjk
                call generic_radial(2,rjk,f(4),p,qjk,dqjk)
            else
                qjk=generic_radial_value(2,rjk,f(4),p)
            end if
        end if
        cosine=max(-1.0_real64,min(1.0_real64,sum(uj*uk)))
    end subroutine

    pure subroutine generic_pair_contracted(f,p,uj,uk,rj,rk,qj,qk,dqj,dqk,coeff,degree,both,gradient,gradient_k)
        !$omp declare target
        logical, intent(in) :: both
        integer, intent(in) :: f(PACKED_FEATURE_FIELDS),degree
        real(real64), intent(in) :: p(7),uj(3),uk(3),rj,rk,qj,qk,dqj,dqk,coeff(0:16)
        real(real64), intent(out) :: gradient(3),gradient_k(3)
        real(real64) :: cosine,qjk,dqjk,ujk(3),a,da,t,product,angular_j,angular_k,radial_j,radial_k,radial_jk
        integer :: d
        gradient=0; gradient_k=0
        call generic_pair_geometry(f,p,uj,uk,rj,rk,cosine,qjk,dqjk,ujk,.true.)
        if (qjk == 0) return
        if (degree > 0) then
            t=0.5_real64*(1+p(4)*cosine)
            a=coeff(degree); da=0
            do d=degree-1,0,-1
                da=da*t+a
                a=a*t+coeff(d)
            end do
            da=da*0.5_real64*p(4)
        else
            call angular_power(cosine,p(4),p(5),f(5),0.5_real64*p(5)*p(4),a,da)
            a=a*coeff(0); da=da*coeff(0)
        end if
        product=qj*qk*qjk
        ! Contract scalar factors before forming Cartesian components. This
        ! shares the j-k radial term and avoids three divisions per vector.
        angular_j=2*da*product/rj
        radial_j=2*a*qk*dqj*qjk-angular_j*cosine
        radial_jk=2*a*qj*qk*dqjk
        gradient=angular_j*uk+radial_j*uj-radial_jk*ujk
        if (both) then
            angular_k=2*da*product/rk
            radial_k=2*a*qj*dqk*qjk-angular_k*cosine
            gradient_k=angular_k*uj+radial_k*uk+radial_jk*ujk
        end if
    end subroutine
    pure real(real64) function generic_pair_value(f,p,uj,uk,rj,rk,qj,qk) result(value)
        !$omp declare target
        integer, intent(in) :: f(PACKED_FEATURE_FIELDS)
        real(real64), intent(in) :: p(7),uj(3),uk(3),rj,rk,qj,qk
        real(real64) :: qjk,rjk,rjk2,cosine
        value = 0
        if (rj <= 1e-12_real64 .or. rk <= 1e-12_real64 .or. rj > p(1) .or. rk > p(1)) return
        qjk = 1
        if (f(1) == 4) then
            rjk2 = sum((rk*uk-rj*uj)**2)
            if (rjk2 <= 1e-24_real64 .or. rjk2 >= p(1)*p(1)) return
            rjk = sqrt(rjk2)
            qjk = generic_radial_value(2,rjk,f(4),p)
        end if
        cosine = max(-1.0_real64,min(1.0_real64,sum(uj*uk)))
        value = 2*angular_value(cosine,p(4),p(5),f(5))*qj*qk*qjk
    end function
    pure real(real64) function sf_cutoff_value(distance,rc,cutoff_type,alpha) result(value)
        !$omp declare target
        include 'cutoff_value.inc'
    end function
    pure real(real64) function sf_cutoff_derivative(distance,rc,cutoff_type,alpha) result(value)
        !$omp declare target
        include 'cutoff_derivative.inc'
    end function
    include 'n2p2_extended.inc'
end module
