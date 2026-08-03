program test_aenet_activation_reference
    use iso_fortran_env, only: real64
    use aenet_network, only: atomic_network
    use feedforward, only: ff_activate
    implicit none

    real(real64), parameter :: inputs(7) = &
        [-20.0_real64, -2.0_real64, -0.7_real64, 0.0_real64, 0.7_real64, 2.0_real64, 20.0_real64]
    type(atomic_network) :: network
    real(real64) :: output, gradient(1), reference, reference_derivative, scale
    integer :: code, i

    network%nlayers = 2
    network%maxnodes = 1
    network%nodes = [1, 1]
    network%weight_offsets = [0, 2]
    network%weights = [1.0_real64, 0.0_real64]

    do code = 0, 4
        network%activation = [code]
        do i = 1, size(inputs)
            call ff_activate(code, inputs(i), reference, reference_derivative)
            call network%input_gradient([inputs(i)], output, gradient)
            scale = max(1.0_real64, abs(reference), abs(output))
            if (abs(output - reference) > 3.0e-15_real64*scale) &
                error stop "aenet activation value mismatch"
            scale = max(1.0_real64, abs(reference_derivative), abs(gradient(1)))
            if (abs(gradient(1) - reference_derivative) > 3.0e-15_real64*scale) &
                error stop "aenet activation derivative mismatch"
        end do
    end do
end program test_aenet_activation_reference
