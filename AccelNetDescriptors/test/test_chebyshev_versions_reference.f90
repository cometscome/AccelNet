program test_chebyshev_versions_reference
    use iso_fortran_env, only: real64
    use accelnet_descriptors
    use chebyshevbasis_fp, only: FP_ChebyshevBasis
    implicit none

    integer, parameter :: number_of_neighbors = 4
    integer, parameter :: number_of_descriptors = 18
    character(len=2), parameter :: names(2) = ["Ti", "O "]
    integer, parameter :: neighbor_species(number_of_neighbors) = [1, 1, 2, 2]
    integer, parameter :: local_types(2) = [1, 2]
    real(real64), parameter :: center(3) = 0.0_real64
    real(real64) :: coordinates(3, number_of_neighbors)
    real(real64) :: parameters(5, number_of_descriptors)
    real(real64) :: original_values(number_of_descriptors), new_values(number_of_descriptors)
    real(real64) :: original_center(3, number_of_descriptors)
    real(real64) :: new_center(3, number_of_descriptors)
    real(real64) :: original_neighbors(3, number_of_descriptors, number_of_neighbors)
    real(real64) :: new_neighbors(3, number_of_descriptors, number_of_neighbors)
    type(FP_ChebyshevBasis) :: original_basis
    type(descriptor_config) :: config
    integer :: version

    coordinates(:, 1) = [1.1_real64, 0.2_real64, -0.1_real64]
    coordinates(:, 2) = [-0.4_real64, 1.3_real64, 0.3_real64]
    coordinates(:, 3) = [0.2_real64, -0.5_real64, 1.5_real64]
    coordinates(:, 4) = [-1.2_real64, -0.7_real64, 0.6_real64]

    do version = 1, 10, 9
        parameters = 0.0_real64
        parameters(:, 1) = [4.5_real64, 4.0_real64, 4.0_real64, 3.0_real64, real(version, real64)]
        original_basis = FP_ChebyshevBasis(parameters, 2, names)
        call original_basis%evaluate(2, center, number_of_neighbors, coordinates, neighbor_species, &
                                     local_types, original_values, original_center, original_neighbors)

        call initialize_config(config, 2, 4.5_real64, 4, 4.0_real64, 3, &
                               version=version, central_type_index=2)
        call evaluate_atom_with_derivatives(config, coordinates, neighbor_species, new_values, &
                                            new_center, new_neighbors)
        call assert_scaled_close("values", version, new_values, original_values)
        call assert_scaled_close_2d("center derivatives", version, new_center, original_center)
        call assert_scaled_close_3d("neighbor derivatives", version, new_neighbors, original_neighbors)
    end do
    write(*, "(A)") "Chebyshev version 1/10 reference tests passed"

contains

    subroutine assert_scaled_close(label, checked_version, actual, expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: checked_version
        real(real64), intent(in) :: actual(:), expected(:)
        real(real64) :: error
        error = maxval(abs(actual - expected)/(1.0_real64 + abs(expected)))
        if (error > 5.0e-14_real64) then
            write(*, "(A,I0,2A,ES24.16)") "version ", checked_version, " ", trim(label), error
            error stop "Chebyshev version reference mismatch"
        end if
    end subroutine assert_scaled_close

    subroutine assert_scaled_close_2d(label, checked_version, actual, expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: checked_version
        real(real64), intent(in) :: actual(:, :), expected(:, :)
        call assert_scaled_close(label, checked_version, reshape(actual, [size(actual)]), &
                                 reshape(expected, [size(expected)]))
    end subroutine assert_scaled_close_2d

    subroutine assert_scaled_close_3d(label, checked_version, actual, expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: checked_version
        real(real64), intent(in) :: actual(:, :, :), expected(:, :, :)
        call assert_scaled_close(label, checked_version, reshape(actual, [size(actual)]), &
                                 reshape(expected, [size(expected)]))
    end subroutine assert_scaled_close_3d

end program test_chebyshev_versions_reference
