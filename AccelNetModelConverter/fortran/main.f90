program accelnet_model_converter_fortran
    use model_conversion, only: convert_n2p2_to_accelnet, convert_accelnet_to_n2p2
    implicit none
    integer, parameter :: PATH_LENGTH = 4096
    integer :: count, i
    character(len=PATH_LENGTH) :: command, input, output
    character(len=PATH_LENGTH), allocatable :: networks(:)

    count = command_argument_count()
    if (count < 1) call usage()
    call get_command_argument(1, command)
    select case(trim(command))
    case("n2p2-to-accelnet")
        if (count /= 3) call usage()
        call get_command_argument(2, input)
        call get_command_argument(3, output)
        call convert_n2p2_to_accelnet(trim(input), trim(output))
    case("accelnet-to-n2p2")
        if (count < 3) call usage()
        call get_command_argument(2, output)
        allocate(networks(count - 2))
        do i = 1, size(networks)
            call get_command_argument(i + 2, networks(i))
        end do
        call convert_accelnet_to_n2p2(trim(output), networks)
    case default
        call usage()
    end select

contains
    subroutine usage()
        write(*, "(A)") "Usage:"
        write(*, "(A)") "  accelnet-model-converter-fortran n2p2-to-accelnet INPUT_DIR OUTPUT_DIR"
        write(*, "(A)") "  accelnet-model-converter-fortran accelnet-to-n2p2 OUTPUT_DIR NETWORK..."
        stop 2
    end subroutine usage
end program accelnet_model_converter_fortran
