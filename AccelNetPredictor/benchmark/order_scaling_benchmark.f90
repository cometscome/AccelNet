program order_scaling_benchmark
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: atomic_structure, neighbor_data, descriptor_config, &
        read_xsf, build_neighbor_list, initialize_config, evaluate_atom, contract_atom_derivatives
    implicit none
    integer, parameter :: orders(*) = [0, 2, 4, 6, 8, 10, 12, 14, 16]
    character(len=16), parameter :: species_names(2) = [character(len=16) :: "Ti", "O"]
    character(len=1024) :: xsf_file, argument
    type(atomic_structure) :: structure
    type(neighbor_data) :: neighbors
    type(descriptor_config) :: config
    real(real64), allocatable :: values(:), coefficients(:), contracted_neighbors(:,:)
    integer, allocatable :: neighbor_species(:)
    real(real64) :: contracted_center(3), checksum, start_time, end_time, value_seconds, force_seconds
    integer :: repeats, order_index, angular_order, moment_count, first, last, iteration

    if (command_argument_count() < 1 .or. command_argument_count() > 2) then
        write(*,'(A)') "usage: order-scaling-benchmark structure.xsf [REPEATS]"
        error stop 2
    end if
    call get_command_argument(1,xsf_file)
    repeats = 2000
    if (command_argument_count() == 2) then
        call get_command_argument(2,argument); read(argument,*) repeats
    end if
    if (repeats < 1) error stop "REPEATS must be positive"
    call read_xsf(trim(xsf_file),species_names,structure)
    call build_neighbor_list(structure,6.5_real64,neighbors,0.35_real64,preserve_legacy_order=.true.)
    first = neighbors%offsets(1); last = neighbors%offsets(2)-1
    allocate(neighbor_species(last-first+1))
    neighbor_species = structure%species(neighbors%atom_indices(first:last))

    write(*,'(A,1X,I0)') "NEIGHBORS",size(neighbor_species)
    write(*,'(A)') "ORDER MOMENTS VALUE_SECONDS_PER_ATOM FORCE_SECONDS_PER_ATOM CHECKSUM"
    do order_index = 1, size(orders)
        angular_order = orders(order_index)
        call initialize_config(config,2,6.5_real64,20,5.0_real64,angular_order,version=0,central_type_index=1)
        moment_count = config%number_of_angular_moments
        allocate(values(config%num_descriptors()),coefficients(config%num_descriptors()), &
                 contracted_neighbors(3,size(neighbor_species)))
        coefficients = 1.0_real64/real(config%num_descriptors(),real64)
        call evaluate_atom(config,neighbors%displacements(:,first:last),neighbor_species,values)
        call contract_atom_derivatives(config,neighbors%displacements(:,first:last),neighbor_species, &
            coefficients,contracted_center,contracted_neighbors)
        call cpu_time(start_time)
        do iteration = 1, repeats
            call evaluate_atom(config,neighbors%displacements(:,first:last),neighbor_species,values)
        end do
        call cpu_time(end_time); value_seconds = (end_time-start_time)/repeats
        call cpu_time(start_time)
        do iteration = 1, repeats
            call contract_atom_derivatives(config,neighbors%displacements(:,first:last),neighbor_species, &
                coefficients,contracted_center,contracted_neighbors)
        end do
        call cpu_time(end_time); force_seconds = (end_time-start_time)/repeats
        checksum = sum(values)+sum(contracted_center)+sum(contracted_neighbors)
        write(*,'(I5,1X,I7,3(1X,ES24.16))') angular_order,moment_count,value_seconds,force_seconds,checksum
        deallocate(values,coefficients,contracted_neighbors)
    end do
end program order_scaling_benchmark
