program compare_n2p2_descriptors
    use iso_fortran_env, only: real64, int64
    implicit none

    character(len=1024) :: reference_file, accelnet_file, line
    integer :: reference_unit, accelnet_unit, ios, natoms, accelnet_atoms, descriptors
    integer :: atom, coefficient, atomic_number_reference, atomic_number_accelnet, atom_index
    integer(int64), allocatable :: descriptor_bits(:)
    real(real64), allocatable :: reference_values(:)
    real(real64) :: energy_row(4), reference, value, difference, scale, maximum_error

    if (command_argument_count() /= 2) &
        error stop "usage: compare_n2p2_descriptors FUNCTION_DATA ACCELNET_HEX"
    call get_command_argument(1, reference_file)
    call get_command_argument(2, accelnet_file)
    open(newunit=reference_unit, file=trim(reference_file), status="old", action="read")
    open(newunit=accelnet_unit, file=trim(accelnet_file), status="old", action="read")

    read(reference_unit, *) natoms
    read(accelnet_unit, "(A)") line
    if (trim(line) /= "ACCELNET_N2P2_DESCRIPTOR_HEX_V1") error stop "invalid AccelNet n2p2 descriptor header"
    read(accelnet_unit, *) accelnet_atoms, descriptors
    if (natoms /= accelnet_atoms .or. natoms <= 0 .or. descriptors <= 0) &
        error stop "n2p2/AccelNet descriptor dimensions differ"
    allocate(reference_values(descriptors), descriptor_bits(descriptors))
    maximum_error = 0.0_real64
    do atom = 1, natoms
        read(reference_unit, *, iostat=ios) atomic_number_reference, reference_values
        if (ios /= 0) error stop "cannot read n2p2 function.data atom row"
        read(accelnet_unit, *, iostat=ios) atom_index, atomic_number_accelnet
        if (ios /= 0 .or. atom_index /= atom .or. atomic_number_reference /= atomic_number_accelnet) &
            error stop "n2p2/AccelNet atom metadata differ"
        read(accelnet_unit, "(*(Z16,1X))", iostat=ios) descriptor_bits
        if (ios /= 0) error stop "cannot read AccelNet descriptor row"
        do coefficient = 1, descriptors
            reference = reference_values(coefficient)
            value = transfer(descriptor_bits(coefficient), value)
            difference = abs(reference - value)
            scale = max(1.0_real64, abs(reference), abs(value))
            maximum_error = max(maximum_error, difference/scale)
            ! nnp-scaling writes function.data with ten digits after the decimal point.
            if (difference > 7.0e-10_real64*scale) then
                write(*, "(A,2(I0,1X),2(ES24.16,1X))") &
                    "n2p2/AccelNet mismatch at atom/coefficient: ", atom, coefficient, reference, value
                error stop "n2p2 descriptor tolerance exceeded"
            end if
        end do
    end do
    read(reference_unit, *, iostat=ios) energy_row
    if (ios /= 0) error stop "missing n2p2 function.data energy row"
    read(reference_unit, "(A)", iostat=ios) line
    if (ios == 0) error stop "extra structures in n2p2 function.data"
    read(accelnet_unit, "(A)", iostat=ios) line
    if (ios == 0) error stop "extra data in AccelNet descriptor file"
    close(reference_unit)
    close(accelnet_unit)
    write(*, "(A,ES12.4)") "n2p2/AccelNet maximum scaled error: ", maximum_error
end program compare_n2p2_descriptors
