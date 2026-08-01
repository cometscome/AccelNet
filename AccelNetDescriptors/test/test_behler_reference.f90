program test_behler_reference
    use iso_fortran_env, only: real64
    use accelnet_descriptors
    use accelnet_setup
    use accelnet_descriptor_models
    use behler2011basis_fp, only: FP_Behler2011basis
    implicit none

    integer, parameter :: number_of_functions = 5
    type(atomic_structure) :: structure
    type(neighbor_data) :: neighbors
    type(descriptor_setup) :: parsed_setup
    type(FP_Behler2011basis) :: original
    character(len=1024) :: input_file, setup_file
    character(len=2), parameter :: names(2) = ["Ti", "O "]
    real(real64) :: parameters(5, number_of_functions)
    integer :: kinds(number_of_functions), environments(2, number_of_functions)
    real(real64), allocatable :: coordinates(:, :), original_values(:), new_values(:)
    real(real64), allocatable :: original_center(:, :), original_neighbors(:, :, :)
    real(real64), allocatable :: new_center(:, :), new_neighbors(:, :, :)
    integer, allocatable :: types(:), local_types(:)
    integer :: ltype(2), atom, first, last, n, j, coefficient
    real(real64) :: difference, scale, maximum_scaled_error

    if (command_argument_count() /= 2) error stop "usage: test_behler_reference XSF SETUP"
    call get_command_argument(1, input_file)
    call get_command_argument(2, setup_file)
    call read_xsf(trim(input_file), names, structure)
    call build_neighbor_list(structure, 4.5_real64, neighbors)
    call read_accelnet_setup(trim(setup_file), names, parsed_setup)
    parameters = 0.0_real64
    parameters(:, 1) = [4.5_real64, -1.0_real64, 3.0_real64, 0.2_real64, 0.0_real64]
    parameters(:, 2) = [4.5_real64, 0.3_real64, 0.7_real64, 0.0_real64, 0.0_real64]
    parameters(:, 3) = [4.5_real64, 1.0_real64, 2.0_real64, 0.15_real64, 0.0_real64]
    parameters(:, 4) = [4.5_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    parameters(:, 5) = [4.5_real64, 1.2_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    kinds = [5, 2, 4, 1, 3]
    environments = 0
    environments(1, :) = [2, 2, 1, 1, 1]
    environments(2, 1) = 2
    environments(2, 3) = 2
    original = FP_Behler2011basis(number_of_functions, 2, number_of_functions, &
        parameters, 2, names, kinds, environments)
    ltype = [1, 2]
    allocate(original_values(number_of_functions), new_values(number_of_functions))
    maximum_scaled_error = 0.0_real64

    do atom = 1, structure%natoms
        first = neighbors%offsets(atom)
        last = neighbors%offsets(atom + 1) - 1
        n = last - first + 1
        allocate(coordinates(3, n), types(n), local_types(n))
        allocate(original_center(3, number_of_functions), original_neighbors(3, number_of_functions, n))
        allocate(new_center(3, number_of_functions), new_neighbors(3, number_of_functions, n))
        do j = 1, n
            coordinates(:, j) = neighbors%positions(:, first + j - 1)
            types(j) = structure%species(neighbors%atom_indices(first + j - 1))
        end do
        call parsed_setup%map_species(types, local_types)
        call original%evaluate(structure%species(atom), structure%positions(:, atom), n, &
            coordinates, types, ltype, original_values, original_center, original_neighbors)
        call evaluate_model_values_derivatives(parsed_setup%model, neighbors%displacements(:, first:last), &
            local_types, new_values, new_center, new_neighbors)
        do coefficient = 1, number_of_functions
            call compare_scalar(original_values(coefficient), new_values(coefficient))
            do j = 1, 3
                call compare_scalar(original_center(j, coefficient), new_center(j, coefficient))
            end do
            do j = 1, n
                call compare_vector(original_neighbors(:, coefficient, j), new_neighbors(:, coefficient, j))
            end do
        end do
        deallocate(coordinates, types, local_types, original_center, original_neighbors, new_center, new_neighbors)
    end do
    write(*, "(A,ES12.4)") "Behler maximum scaled error: ", maximum_scaled_error

contains
    subroutine compare_scalar(reference, value)
        real(real64), intent(in) :: reference, value
        difference = abs(reference - value)
        scale = max(1.0_real64, abs(reference), abs(value))
        maximum_scaled_error = max(maximum_scaled_error, difference/scale)
        if (difference > 2.0e-12_real64*scale) then
            write(*, "(A,2(ES24.16,1X))") "Behler reference mismatch: ", reference, value
            error stop "Behler reference mismatch"
        end if
    end subroutine
    subroutine compare_vector(reference, value)
        real(real64), intent(in) :: reference(3), value(3)
        integer :: component
        do component = 1, 3
            call compare_scalar(reference(component), value(component))
        end do
    end subroutine
end program test_behler_reference
