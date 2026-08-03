program dump_n2p2_descriptors
    use iso_fortran_env, only: real64, int64
    use accelnet_descriptors, only: atomic_structure, neighbor_data, read_xsf, build_neighbor_list
    use accelnet_descriptor_models, only: evaluate_model_values
    use accelnet_setup, only: descriptor_setup
    use n2p2_network, only: load_n2p2_setups, atomic_number
    implicit none

    type(descriptor_setup), allocatable :: setups(:)
    type(atomic_structure) :: structure
    type(neighbor_data) :: neighbors
    character(len=16), allocatable :: species_names(:)
    character(len=1024) :: model_directory, xsf_file, output_file
    real(real64), allocatable :: values(:)
    integer, allocatable :: global_neighbors(:), local_neighbors(:)
    integer :: atom, coefficient, descriptor_count, first, last, n, unit, species
    real(real64) :: maximum_cutoff, minimum_distance

    if (command_argument_count() /= 3) &
        error stop "usage: dump_n2p2_descriptors MODEL_DIRECTORY XSF OUTPUT"
    call get_command_argument(1, model_directory)
    call get_command_argument(2, xsf_file)
    call get_command_argument(3, output_file)

    call load_n2p2_setups(trim(model_directory), setups, species_names)
    call read_xsf(trim(xsf_file), species_names, structure)
    descriptor_count = setups(1)%model%num_descriptors()
    maximum_cutoff = setups(1)%model%maximum_cutoff
    minimum_distance = setups(1)%minimum_distance
    do species = 2, size(setups)
        if (setups(species)%model%num_descriptors() /= descriptor_count) &
            error stop "descriptor dump requires equal descriptor counts for all species"
        maximum_cutoff = max(maximum_cutoff, setups(species)%model%maximum_cutoff)
        minimum_distance = min(minimum_distance, setups(species)%minimum_distance)
    end do

    call build_neighbor_list(structure, maximum_cutoff, neighbors, minimum_distance)
    allocate(values(descriptor_count))
    open(newunit=unit, file=trim(output_file), status="replace", action="write")
    write(unit, "(A)") "ACCELNET_N2P2_DESCRIPTOR_HEX_V1"
    write(unit, "(2(I0,1X))") structure%natoms, descriptor_count
    do atom = 1, structure%natoms
        species = structure%species(atom)
        first = neighbors%offsets(atom)
        last = neighbors%offsets(atom + 1) - 1
        n = last - first + 1
        allocate(global_neighbors(n), local_neighbors(n))
        global_neighbors = structure%species(neighbors%atom_indices(first:last))
        call setups(species)%map_species(global_neighbors, local_neighbors)
        call evaluate_model_values(setups(species)%model, neighbors%displacements(:, first:last), &
            local_neighbors, values)
        write(unit, "(I0,1X,I0)") atom, atomic_number(species_names(species))
        write(unit, "(*(Z16.16,1X))") &
            (transfer(values(coefficient), 0_int64), coefficient = 1, descriptor_count)
        deallocate(global_neighbors, local_neighbors)
    end do
    close(unit)
end program dump_n2p2_descriptors
