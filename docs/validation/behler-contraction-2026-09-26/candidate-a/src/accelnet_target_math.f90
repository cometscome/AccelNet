! Scalar formulas are shared with CPU sources; this module compiles their device variants.
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
    public :: generic_pair_geometry, generic_pair_contracted, generic_pair_value
    public :: generic_radial, generic_radial_value, generic_lj
    public :: angular_power, angular_value
    public :: target_cutoff_value, target_cutoff_derivative, target_activate, target_activation_derivative
contains
    pure real(real64) function target_cutoff_value(distance, rc, cutoff_type, alpha) result(value)
        !$omp declare target
        include 'cutoff_value.inc'
    end function target_cutoff_value

    pure real(real64) function target_cutoff_derivative(distance, rc, cutoff_type, alpha) result(value)
        !$omp declare target
        include 'cutoff_derivative.inc'
    end function target_cutoff_derivative

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
    pure subroutine generic_radial(kind, r, ct, p, value, derivative)
        !$omp declare target
        integer, intent(in) :: kind,ct
        real(real64), intent(in) :: r,p(7)
        real(real64), intent(out) :: value,derivative
        real(real64) :: fc,dfc,h,dh
        value = 0; derivative = 0
        if (r <= 1e-12_real64 .or. r > p(1)) return
        fc = target_cutoff_value(r,p(1),ct,p(7))
        dfc = target_cutoff_derivative(r,p(1),ct,p(7))
        h = 1; dh = 0
        select case(kind)
        case(2,4,5)
            h = exp(-p(2)*(r-p(3))**2)
            dh = -2*p(2)*(r-p(3))*h
        case(3)
            h = cos(p(6)*r); dh = -p(6)*sin(p(6)*r)
        case(6)
            h = r**(-6); dh = -6*h/r
        case(7)
            h = r**(-12); dh = -12*h/r
        end select
        value = fc*h; derivative = dfc*h+fc*dh
    end subroutine

    ! Unordered G4/G5 pair. derivative is with respect to the first neighbor.
    ! The force kernel visits the reverse pair to obtain the second derivative.

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
        integer, intent(in) :: f(:)
        real(real64), intent(in) :: p(7),uj(3),uk(3),rj,rk
        logical, intent(in) :: force
        real(real64), intent(out) :: cosine,qjk,dqjk,ujk(3)
        real(real64) :: rjk
        cosine=0; qjk=0; dqjk=0; ujk=0
        if (rj <= 1e-12_real64 .or. rk <= 1e-12_real64 .or. rj > p(1) .or. rk > p(1)) return
        qjk=1
        if (f(1) == 4) then
            ujk=rk*uk-rj*uj; rjk=sqrt(sum(ujk**2))
            if (rjk <= 1e-12_real64 .or. rjk > p(1)) then
                qjk=0; return
            end if
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
        integer, intent(in) :: f(:),degree
        real(real64), intent(in) :: p(7),uj(3),uk(3),rj,rk,qj,qk,dqj,dqk,coeff(0:16)
        real(real64), intent(out) :: gradient(3),gradient_k(3)
        real(real64) :: cosine,qjk,dqjk,ujk(3),a,da,t,product
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
        gradient=2*(da*product*(uk-cosine*uj)/rj+a*qk*(dqj*qjk*uj-qj*dqjk*ujk))
        if (both) gradient_k=2*(da*product*(uj-cosine*uk)/rk+a*qj*(dqk*qjk*uk+qk*dqjk*ujk))
    end subroutine
    pure real(real64) function generic_pair_value(f,p,uj,uk,rj,rk,qj,qk) result(value)
        !$omp declare target
        integer, intent(in) :: f(:)
        real(real64), intent(in) :: p(7),uj(3),uk(3),rj,rk,qj,qk
        real(real64) :: qjk,rjk,cosine
        value = 0
        if (rj <= 1e-12_real64 .or. rk <= 1e-12_real64 .or. rj > p(1) .or. rk > p(1)) return
        qjk = 1
        if (f(1) == 4) then
            rjk = sqrt(sum((rk*uk-rj*uj)**2))
            qjk = generic_radial_value(2,rjk,f(4),p)
        end if
        cosine = max(-1.0_real64,min(1.0_real64,sum(uj*uk)))
        value = 2*angular_value(cosine,p(4),p(5),f(5))*qj*qk*qjk
    end function
end module
