program compare_descriptor_files
    use iso_fortran_env, only: real64, int64
    implicit none

    character(len=1024) :: first_file, second_file, line1, line2
    integer :: first_unit, second_unit, ios, s, i, g
    integer :: nstructures1, nstructures2, nspecies1, nspecies2, ndescriptors1, ndescriptors2
    integer :: natoms1, natoms2, atom1, atom2, species1, species2
    integer(int64), allocatable :: bits1(:), bits2(:)
    real(real64) :: value1, value2, difference, scale, maximum_error

    if (command_argument_count() /= 2) error stop "usage: compare_descriptor_files FIRST SECOND"
    call get_command_argument(1, first_file)
    call get_command_argument(2, second_file)
    open(newunit=first_unit, file=trim(first_file), status="old", action="read")
    open(newunit=second_unit, file=trim(second_file), status="old", action="read")

    read(first_unit, "(A)") line1
    read(second_unit, "(A)") line2
    if (trim(line1) /= trim(line2)) error stop "descriptor format headers differ"
    read(first_unit, *) nstructures1, nspecies1, ndescriptors1
    read(second_unit, *) nstructures2, nspecies2, ndescriptors2
    if (nstructures1 /= nstructures2 .or. nspecies1 /= nspecies2 .or. &
        ndescriptors1 /= ndescriptors2) error stop "descriptor dimensions differ"
    allocate(bits1(ndescriptors1), bits2(ndescriptors1))

    maximum_error = 0.0_real64
    do s = 1, nstructures1
        read(first_unit, "(A)") line1
        read(second_unit, "(A)") line2
        if (trim(line1) /= trim(line2)) error stop "structure labels differ"
        read(first_unit, *) natoms1
        read(second_unit, *) natoms2
        if (natoms1 /= natoms2) error stop "atom counts differ"
        do i = 1, natoms1
            read(first_unit, *) atom1, species1
            read(second_unit, *) atom2, species2
            if (atom1 /= atom2 .or. species1 /= species2) error stop "atom metadata differ"
            read(first_unit, "(*(Z16,1X))", iostat=ios) bits1
            if (ios /= 0) error stop "cannot read first descriptor row"
            read(second_unit, "(*(Z16,1X))", iostat=ios) bits2
            if (ios /= 0) error stop "cannot read second descriptor row"
            do g = 1, ndescriptors1
                value1 = transfer(bits1(g), value1)
                value2 = transfer(bits2(g), value2)
                difference = abs(value1 - value2)
                scale = max(1.0_real64, abs(value1), abs(value2))
                maximum_error = max(maximum_error, difference/scale)
                if (difference > 5.0e-13_real64*scale) then
                    write(*, "(A,3(I0,1X),2(ES24.16,1X))") &
                        "descriptor mismatch at structure/atom/coefficient: ", s, i, g, value1, value2
                    error stop "descriptor tolerance exceeded"
                end if
            end do
        end do
    end do
    read(first_unit, "(A)", iostat=ios) line1
    if (ios == 0) error stop "extra data in first descriptor file"
    read(second_unit, "(A)", iostat=ios) line2
    if (ios == 0) error stop "extra data in second descriptor file"
    close(first_unit)
    close(second_unit)
    write(*, "(A,ES12.4)") "maximum scaled error: ", maximum_error
end program compare_descriptor_files
