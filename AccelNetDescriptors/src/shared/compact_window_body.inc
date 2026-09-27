        real(real64), intent(in) :: x,left,right
        integer, intent(in) :: subtype
        real(real64), intent(out) :: value,derivative
        real(real64), intent(in), optional :: angular_cache(2)
        real(real64) :: y,t,chain,width
        value=0; derivative=0
        if (x <= left .or. x >= right) return
        if (present(angular_cache)) then
            y=(x-angular_cache(2))*angular_cache(1)
            t=abs(y); chain=sign(angular_cache(1),y)
        else
            width=0.5_real64*(right-left)
            y=(x-0.5_real64*(right+left))/width
            t=abs(y); chain=sign(1.0_real64,y)/width
        end if
        include 'compact_polynomial_body.inc'
