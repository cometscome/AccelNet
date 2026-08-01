program descriptor_generate
    use accelnet_descriptors
    implicit none

    type(descriptor_config) :: config
    type(atomic_structure), allocatable :: structures(:)
    character(len=1024), allocatable :: labels(:)
    character(len=1024) :: output_file, input_file
    character(len=2), parameter :: species_names(2) = ["Ti", "O "]
    integer :: nargs, i

    nargs = command_argument_count()
    if (nargs < 2) then
        write(*, "(A)") "usage: accelnet-descriptor OUTPUT XSF [XSF ...]"
        error stop 2
    end if
    call get_command_argument(1, output_file)
    allocate(structures(nargs - 1), labels(nargs - 1))
    do i = 2, nargs
        call get_command_argument(i, input_file)
        labels(i - 1) = trim(input_file)
        call read_xsf(trim(input_file), species_names, structures(i - 1))
    end do

    call initialize_config(config, 2, 6.5d0, 20, 5.0d0, 6, version=0)
    call write_descriptor_file(trim(output_file), config, structures, labels)
end program descriptor_generate
