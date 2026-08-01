program test_lj_reference
    use iso_fortran_env, only: real64
    use accelnet_descriptors
    use accelnet_lj
    use LJbasis_fp, only: FP_LJBasis
    implicit none

    type(atomic_structure) :: structure
    type(neighbor_data) :: neighbors
    type(lj_config) :: config
    type(FP_LJBasis) :: original_basis
    character(len=1024) :: input_file
    character(len=2), parameter :: names(2) = ["Ti", "O "]
    real(real64), allocatable :: parameters(:, :), original_values(:), new_values(:), coordinates(:, :)
    integer, allocatable :: types(:)
    integer :: ltype(2), atom, j, first, last, n
    real(real64) :: maximum_scaled_error, difference, scale

    if (command_argument_count() /= 1) error stop "usage: test_lj_reference XSF"
    call get_command_argument(1, input_file)
    call read_xsf(trim(input_file), names, structure)
    call build_neighbor_list(structure, 6.5_real64, neighbors)
    call initialize_lj_config(config, 2, 6.5_real64)
    allocate(parameters(5, config%num_descriptors()), source=0.0_real64)
    parameters(1, 1) = 6.5_real64
    original_basis = FP_LJBasis(parameters, 2, names)
    allocate(original_values(config%num_descriptors()), new_values(config%num_descriptors()))
    ltype = [1, 2]
    maximum_scaled_error = 0.0_real64

    do atom = 1, structure%natoms
        first = neighbors%offsets(atom)
        last = neighbors%offsets(atom + 1) - 1
        n = last - first + 1
        allocate(coordinates(3, n), types(n))
        do j = 1, n
            coordinates(:, j) = neighbors%positions(:, first + j - 1)
            types(j) = structure%species(neighbors%atom_indices(first + j - 1))
        end do
        original_values = 0.0_real64
        call original_basis%evaluate(structure%species(atom), structure%positions(:, atom), n, &
                                     coordinates, types, ltype, original_values)
        call evaluate_lj_values(config, neighbors%displacements(:, first:last), types, new_values)
        do j = 1, config%num_descriptors()
            difference = abs(original_values(j) - new_values(j))
            scale = max(1.0_real64, abs(original_values(j)), abs(new_values(j)))
            maximum_scaled_error = max(maximum_scaled_error, difference/scale)
            if (difference > 5.0e-13_real64*scale) error stop "LJ reference values differ"
        end do
        deallocate(coordinates, types)
    end do
    write(*, "(A,ES12.4)") "LJ maximum scaled error: ", maximum_scaled_error
end program test_lj_reference
