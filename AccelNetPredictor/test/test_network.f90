program test_network
    use iso_fortran_env, only: real64
    use aenet_network, only: atomic_network, read_aenet_network_ascii
    implicit none
    type(atomic_network) :: network
    character(len=1024) :: filename
    real(real64), allocatable :: input(:), gradient(:)
    real(real64) :: output
    call get_command_argument(1, filename)
    call read_aenet_network_ascii(trim(filename), network)
    if (network%nodes(1) /= 56) error stop "unexpected input size"
    if (trim(network%atomtype) /= "Ti") error stop "unexpected species"
    allocate(input(network%nodes(1)), gradient(network%nodes(1))); input = 0.0_real64
    call network%input_gradient(input, output, gradient)
    if (.not. (output == output) .or. any(gradient /= gradient)) error stop "non-finite network output"
end program test_network
