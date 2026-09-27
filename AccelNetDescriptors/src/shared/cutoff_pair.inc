! Shared fused value/derivative formula for serial and device kernels.
        use accelnet_descriptors, only: CUTOFF_TANHU,CUTOFF_TANH
        real(real64), intent(in) :: distance,rc,alpha
        integer, intent(in) :: kind
        real(real64), intent(out) :: value,derivative
        real(real64) :: t
        ! The water model uses tanh cutoffs. Reuse the expensive transcendental
        ! for value and derivative, just as compact_window reuses its polynomial.
        if (kind == CUTOFF_TANHU .or. kind == CUTOFF_TANH) then
            value=0; derivative=0
            if (distance >= rc) return
            t=tanh(1.0_real64-distance/rc)
            value=t*t*t
            derivative=3.0_real64*t*t*(t*t-1.0_real64)/rc
            if (kind == CUTOFF_TANH) then
                value=value/tanh(1.0_real64)**3
                derivative=derivative/tanh(1.0_real64)**3
            end if
        else
            value=target_cutoff_value(distance,rc,kind,alpha)
            derivative=target_cutoff_derivative(distance,rc,kind,alpha)
        end if
