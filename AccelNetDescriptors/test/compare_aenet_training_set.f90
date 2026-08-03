program compare_aenet_training_set
    use iso_fortran_env, only: real64, int64
    implicit none

    character(len=1024) :: aenet_file, accelnet_file, line
    character(len=16), allocatable :: species_names(:)
    integer :: aenet_unit, accelnet_unit, ios, structure, atom, coefficient
    integer :: ntypes, total_atoms, nstructures, nstructures_accelnet
    integer :: ntypes_structure, nspecies_accelnet, ndescriptors
    integer :: natoms, natoms_accelnet, atom_type, atom_index, accelnet_species, nsf
    integer(int64), allocatable :: descriptor_bits(:)
    real(real64), allocatable :: aenet_values(:)
    real(real64) :: normalization_scale, normalization_shift, atomic_energies(2)
    real(real64) :: energy_statistics(3), cohesive_energy, coordinates(3), forces(3)
    real(real64) :: reference, value, difference, scale, maximum_error
    logical :: normalized

    if (command_argument_count() /= 2) &
        error stop "usage: compare_aenet_training_set AENET_ASCII ACCELNET_HEX"
    call get_command_argument(1, aenet_file)
    call get_command_argument(2, accelnet_file)
    open(newunit=aenet_unit, file=trim(aenet_file), status="old", action="read")
    open(newunit=accelnet_unit, file=trim(accelnet_file), status="old", action="read")

    ! Header emitted by aenet's official trnset2ASCII.x --raw.
    read(aenet_unit, "(A)") line
    read(aenet_unit, *) normalized
    read(aenet_unit, *) normalization_scale
    read(aenet_unit, *) normalization_shift
    if (normalized .or. normalization_scale /= 1.0_real64 .or. normalization_shift /= 0.0_real64) &
        error stop "aenet training set is not raw"
    read(aenet_unit, *) ntypes
    allocate(species_names(ntypes))
    do atom_type = 1, ntypes
        read(aenet_unit, "(A)") species_names(atom_type)
    end do
    if (ntypes /= 2) error stop "reference fixture must contain two species"
    read(aenet_unit, *) atomic_energies
    read(aenet_unit, *) total_atoms
    read(aenet_unit, *) nstructures
    read(aenet_unit, *) energy_statistics

    read(accelnet_unit, "(A)") line
    if (trim(line) /= "ACCELNET_DESCRIPTOR_HEX_V1") error stop "invalid AccelNet descriptor header"
    read(accelnet_unit, *) nstructures_accelnet, nspecies_accelnet, ndescriptors
    if (nstructures /= nstructures_accelnet .or. ntypes /= nspecies_accelnet) &
        error stop "reference dimensions differ"
    allocate(descriptor_bits(ndescriptors), aenet_values(ndescriptors))

    maximum_error = 0.0_real64
    do structure = 1, nstructures
        read(aenet_unit, "(A)") line
        read(aenet_unit, *) natoms, ntypes_structure
        read(aenet_unit, *) cohesive_energy
        if (ntypes_structure /= ntypes) error stop "aenet structure species count differs"
        read(accelnet_unit, "(A)") line
        read(accelnet_unit, *) natoms_accelnet
        if (natoms /= natoms_accelnet) error stop "structure atom counts differ"
        do atom = 1, natoms
            read(aenet_unit, *) atom_type
            read(aenet_unit, *) coordinates
            read(aenet_unit, *) forces
            read(aenet_unit, *) nsf
            if (nsf /= ndescriptors) error stop "descriptor counts differ"
            read(aenet_unit, *) aenet_values

            read(accelnet_unit, *) atom_index, accelnet_species
            if (atom_index /= atom .or. accelnet_species /= atom_type) &
                error stop "atom metadata differ"
            read(accelnet_unit, "(*(Z16,1X))", iostat=ios) descriptor_bits
            if (ios /= 0) error stop "cannot read AccelNet descriptor row"
            do coefficient = 1, ndescriptors
                reference = aenet_values(coefficient)
                value = transfer(descriptor_bits(coefficient), value)
                difference = abs(reference - value)
                scale = max(1.0_real64, abs(reference), abs(value))
                maximum_error = max(maximum_error, difference/scale)
                if (difference > 8.0e-12_real64*scale) then
                    write(*, "(A,3(I0,1X),2(ES24.16,1X))") &
                        "aenet/AccelNet mismatch at structure/atom/coefficient: ", &
                        structure, atom, coefficient, reference, value
                    error stop "aenet descriptor tolerance exceeded"
                end if
            end do
        end do
    end do
    read(aenet_unit, "(A)", iostat=ios) line
    if (ios == 0) error stop "extra data in aenet training set"
    read(accelnet_unit, "(A)", iostat=ios) line
    if (ios == 0) error stop "extra data in AccelNet descriptor file"
    if (total_atoms <= 0) error stop "empty aenet training set"
    close(aenet_unit)
    close(accelnet_unit)
    write(*, "(A,ES12.4)") "aenet/AccelNet maximum scaled error: ", maximum_error
end program compare_aenet_training_set
