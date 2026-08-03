program test_network_activations
    use iso_fortran_env, only: real64
    use aenet_network, only: atomic_network, read_aenet_network_ascii, ACTIVATION_AENET_MTANH, &
        ACTIVATION_AENET_TWIST, ACTIVATION_SOFTPLUS
    implicit none
    type(atomic_network) :: network
    real(real64), parameter :: a = 1.7159_real64
    real(real64), parameter :: b = 0.666666666666667_real64
    real(real64), parameter :: c = 0.1_real64
    real(real64), parameter :: x = 0.7_real64
    real(real64) :: output, gradient(1), expected, expected_gradient, tanhbx
    character(len=1024) :: legacy_file

    network%nlayers = 2
    network%maxnodes = 1
    network%nodes = [1, 1]
    network%weight_offsets = [0, 2]
    network%weights = [1.0_real64, 0.0_real64]
    network%activation = [ACTIVATION_AENET_MTANH]

    tanhbx = tanh(b*x)
    expected = a*tanhbx
    expected_gradient = a*b*(1.0_real64 - tanhbx*tanhbx)
    call network%input_gradient([x], output, gradient)
    call assert_close(output, expected, "aenet mtanh value")
    call assert_close(gradient(1), expected_gradient, "aenet mtanh derivative")

    network%activation = [ACTIVATION_AENET_TWIST]
    call network%input_gradient([x], output, gradient)
    call assert_close(output, expected + c*x, "aenet twist value")
    call assert_close(gradient(1), expected_gradient + c, "aenet twist derivative")

    network%activation = [ACTIVATION_SOFTPLUS]
    expected = log(1.0_real64 + exp(x))
    expected_gradient = 1.0_real64/(1.0_real64 + exp(-x))
    call network%input_gradient([x], output, gradient)
    call assert_close(output, expected, "n2p2 softplus value")
    call assert_close(gradient(1), expected_gradient, "n2p2 softplus derivative")

    call get_command_argument(1, legacy_file)
    call read_aenet_network_ascii(trim(legacy_file), network)
    if (network%activation(1) /= ACTIVATION_SOFTPLUS) &
        error stop "legacy n2p2 softplus activation was not upgraded"

contains

    subroutine assert_close(actual, reference, label)
        real(real64), intent(in) :: actual, reference
        character(len=*), intent(in) :: label
        if (abs(actual - reference) > 2.0e-15_real64) then
            write(*, "(A,2(1X,ES24.16))") trim(label), actual, reference
            error stop "activation mismatch"
        end if
    end subroutine assert_close
end program test_network_activations
