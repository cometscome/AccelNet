program test_n2p2_activation_reference
    use iso_c_binding, only: c_double, c_int
    use iso_fortran_env, only: real64
    use aenet_network, only: atomic_network
    implicit none

    interface
        subroutine n2p2_activation_value_derivative(code, input, value, derivative) bind(C)
            import :: c_double, c_int
            integer(c_int), value :: code
            real(c_double), value :: input
            real(c_double), intent(out) :: value, derivative
        end subroutine n2p2_activation_value_derivative
    end interface

    integer, parameter :: codes(10) = [0, 1, 2, 5, 6, 7, 8, 9, 10, 11]
    real(real64), parameter :: inputs(7) = &
        [-20.0_real64, -2.0_real64, -0.7_real64, 0.0_real64, 0.7_real64, 2.0_real64, 20.0_real64]
    type(atomic_network) :: network
    real(real64) :: output, gradient(1), reference, reference_derivative, scale
    integer :: ic, i

    network%nlayers = 2
    network%maxnodes = 1
    network%nodes = [1, 1]
    network%weight_offsets = [0, 2]
    network%weights = [1.0_real64, 0.0_real64]

    do ic = 1, size(codes)
        network%activation = [codes(ic)]
        do i = 1, size(inputs)
            call n2p2_activation_value_derivative(int(codes(ic), c_int), inputs(i), &
                reference, reference_derivative)
            call network%input_gradient([inputs(i)], output, gradient)
            scale = max(1.0_real64, abs(reference), abs(output))
            if (abs(output - reference) > 4.0e-15_real64*scale) &
                error stop "n2p2 activation value mismatch"
            scale = max(1.0_real64, abs(reference_derivative), abs(gradient(1)))
            if (abs(gradient(1) - reference_derivative) > 4.0e-15_real64*scale) &
                error stop "n2p2 activation derivative mismatch"
        end do
    end do
end program test_n2p2_activation_reference
