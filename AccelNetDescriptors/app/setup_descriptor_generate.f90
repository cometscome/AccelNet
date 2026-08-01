program setup_descriptor_generate
    use iso_fortran_env, only: real64, int64
    use accelnet_descriptors
    use accelnet_descriptor_models
    use accelnet_setup
    implicit none

    type(descriptor_setup), allocatable :: setups(:)
    type(atomic_structure), allocatable :: structures(:)
    type(neighbor_data) :: neighbors
    character(len=16), allocatable :: global_species(:)
    character(len=1024), allocatable :: setup_files(:), labels(:)
    character(len=1024) :: output_file, argument
    real(real64), allocatable :: values(:)
    integer, allocatable :: global_neighbors(:), local_neighbors(:)
    integer :: nargs, number_of_setups, number_of_structures
    integer :: s, atom, first, last, n, coefficient, unit, descriptor_count
    real(real64) :: maximum_cutoff
    real(real64) :: minimum_distance

    nargs = command_argument_count()
    if (nargs < 4) then
        write(*, "(A)") "usage: accelnet-setup-descriptor OUTPUT NSETUP SETUP... XSF..."
        error stop 2
    end if
    call get_command_argument(1, output_file)
    call get_command_argument(2, argument)
    read(argument, *) number_of_setups
    if (number_of_setups < 1 .or. nargs < number_of_setups + 3) error stop "invalid setup count"
    allocate(setup_files(number_of_setups))
    do s = 1, number_of_setups
        call get_command_argument(s + 2, setup_files(s))
    end do
    call read_accelnet_setup_set(setup_files, setups, global_species)
    descriptor_count = setups(1)%model%num_descriptors()
    maximum_cutoff = setups(1)%model%maximum_cutoff
    minimum_distance = setups(1)%minimum_distance
    do s = 2, size(setups)
        if (setups(s)%model%num_descriptors() /= descriptor_count) &
            error stop "HEX regression output requires equal descriptor counts for all central species"
        maximum_cutoff = max(maximum_cutoff, setups(s)%model%maximum_cutoff)
        minimum_distance = min(minimum_distance, setups(s)%minimum_distance)
    end do

    number_of_structures = nargs - number_of_setups - 2
    allocate(structures(number_of_structures), labels(number_of_structures))
    do s = 1, number_of_structures
        call get_command_argument(number_of_setups + 2 + s, labels(s))
        call read_xsf(trim(labels(s)), global_species, structures(s))
    end do

    open(newunit=unit, file=trim(output_file), status="replace", action="write")
    write(unit, "(A)") "ACCELNET_DESCRIPTOR_HEX_V1"
    write(unit, "(3(I0,1X))") number_of_structures, number_of_setups, descriptor_count
    allocate(values(descriptor_count))
    do s = 1, number_of_structures
        write(unit, "(A)") trim(labels(s))
        write(unit, "(I0)") structures(s)%natoms
        call build_neighbor_list(structures(s), maximum_cutoff, neighbors, minimum_distance)
        do atom = 1, structures(s)%natoms
            first = neighbors%offsets(atom)
            last = neighbors%offsets(atom + 1) - 1
            n = last - first + 1
            allocate(global_neighbors(n), local_neighbors(n))
            global_neighbors = structures(s)%species(neighbors%atom_indices(first:last))
            call setups(structures(s)%species(atom))%map_species(global_neighbors, local_neighbors)
            call evaluate_model_values(setups(structures(s)%species(atom))%model, &
                neighbors%displacements(:, first:last), local_neighbors, values)
            write(unit, "(I0,1X,I0)") atom, structures(s)%species(atom)
            write(unit, "(*(Z16.16,1X))") &
                (transfer(values(coefficient), 0_int64), coefficient = 1, descriptor_count)
            deallocate(global_neighbors, local_neighbors)
        end do
    end do
    close(unit)
end program setup_descriptor_generate
