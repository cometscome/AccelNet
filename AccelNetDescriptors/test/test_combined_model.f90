program test_combined_model
    use iso_fortran_env, only: real64
    use accelnet_descriptors
    use accelnet_lj
    use accelnet_behler
    use accelnet_descriptor_models
    implicit none

    type(descriptor_config) :: chebyshev
    type(behler_config) :: behler
    type(lj_config) :: lj
    type(descriptor_model) :: model
    real(real64) :: displacements(3, 3)
    integer :: species(3), nc, nb, nl, first
    real(real64), allocatable :: values(:), center(:, :), neighbors(:, :, :)
    real(real64), allocatable :: local_values(:), local_center(:, :), local_neighbors(:, :, :)

    displacements(:, 1) = [1.1_real64, 0.2_real64, -0.1_real64]
    displacements(:, 2) = [-0.4_real64, 1.3_real64, 0.3_real64]
    displacements(:, 3) = [0.2_real64, -0.5_real64, 1.5_real64]
    species = [1, 2, 1]
    call initialize_config(chebyshev, 2, 4.5_real64, 4, 4.5_real64, 2, version=0)
    call initialize_behler_config(behler, 2)
    call add_g1(behler, 1, 4.5_real64)
    call add_g2(behler, 2, 4.5_real64, 0.3_real64, 0.7_real64)
    call add_g4(behler, 1, 2, 4.5_real64, 1.0_real64, 2.0_real64, 0.15_real64)
    call initialize_lj_config(lj, 2, 4.5_real64)
    call add_chebyshev(model, chebyshev)
    call add_behler(model, behler)
    call add_lj(model, lj)
    nc = chebyshev%num_descriptors()
    nb = behler%num_descriptors()
    nl = lj%num_descriptors()
    if (model%num_descriptors() /= nc + nb + nl) error stop "combined model size"
    allocate(values(model%num_descriptors()), center(3, model%num_descriptors()))
    allocate(neighbors(3, model%num_descriptors(), size(species)))
    allocate(local_values(max(nc, nb, nl)), local_center(3, max(nc, nb, nl)))
    allocate(local_neighbors(3, max(nc, nb, nl), size(species)))
    call evaluate_model_values_derivatives(model, displacements, species, values, center, neighbors)

    call evaluate_atom_with_derivatives(chebyshev, displacements, species, local_values(1:nc), &
                                        local_center(:, 1:nc), local_neighbors(:, 1:nc, :))
    call compare_block(1, nc)
    first = nc + 1
    call evaluate_behler_values_derivatives(behler, displacements, species, local_values(1:nb), &
                                            local_center(:, 1:nb), local_neighbors(:, 1:nb, :))
    call compare_block(first, nb)
    first = nc + nb + 1
    call evaluate_lj_values_derivatives(lj, displacements, species, local_values(1:nl), &
                                        local_center(:, 1:nl), local_neighbors(:, 1:nl, :))
    call compare_block(first, nl)
    write(*, "(A)") "combined Chebyshev + Behler + LJ model test passed"

contains
    subroutine compare_block(offset, count)
        integer, intent(in) :: offset, count
        if (maxval(abs(values(offset:offset + count - 1) - local_values(1:count))) > 0.0_real64) &
            error stop "combined model values differ"
        if (maxval(abs(center(:, offset:offset + count - 1) - local_center(:, 1:count))) > 0.0_real64) &
            error stop "combined model center derivatives differ"
        if (maxval(abs(neighbors(:, offset:offset + count - 1, :) - &
                       local_neighbors(:, 1:count, :))) > 0.0_real64) &
            error stop "combined model neighbor derivatives differ"
    end subroutine
end program test_combined_model
