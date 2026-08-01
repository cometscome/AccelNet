program dump_descriptor_reference
    use iso_fortran_env, only: real64
    use accelnet_lj
    use accelnet_behler
    implicit none

    type(lj_config) :: lj
    type(behler_config) :: behler
    real(real64) :: displacements(3, 3)
    integer :: species(3)
    real(real64), allocatable :: values(:), center(:, :), neighbors(:, :, :)

    displacements(:, 1) = [1.1_real64, 0.2_real64, -0.1_real64]
    displacements(:, 2) = [-0.4_real64, 1.3_real64, 0.3_real64]
    displacements(:, 3) = [0.2_real64, -0.5_real64, 1.5_real64]
    species = [2, 1, 2]

    call initialize_lj_config(lj, 2, 4.5_real64)
    allocate(values(lj%num_descriptors()), center(3, lj%num_descriptors()))
    allocate(neighbors(3, lj%num_descriptors(), size(species)))
    call evaluate_lj_values_derivatives(lj, displacements, species, values, center, neighbors)
    call dump_arrays(values, center, neighbors)
    deallocate(values, center, neighbors)

    call initialize_behler_config(behler, 2)
    call add_g1(behler, 1, 4.5_real64)
    call add_g2(behler, 2, 4.5_real64, 0.3_real64, 0.7_real64)
    call add_g3(behler, 1, 4.5_real64, 1.2_real64)
    call add_g4(behler, 1, 2, 4.5_real64, 1.0_real64, 2.0_real64, 0.15_real64)
    call add_g5(behler, 1, 2, 4.5_real64, -1.0_real64, 3.0_real64, 0.2_real64)
    allocate(values(behler%num_descriptors()), center(3, behler%num_descriptors()))
    allocate(neighbors(3, behler%num_descriptors(), size(species)))
    call evaluate_behler_values_derivatives(behler, displacements, species, values, center, neighbors)
    call dump_arrays(values, center, neighbors)

contains

    subroutine dump_arrays(v, dc, dn)
        real(real64), intent(in) :: v(:), dc(:, :), dn(:, :, :)
        integer :: i, j, k
        write(*, "(I0,1X,I0)") size(v), size(dn, 3)
        do i = 1, size(v)
            write(*, "(ES25.17E3)") v(i)
        end do
        do j = 1, size(dc, 2)
            do i = 1, size(dc, 1)
                write(*, "(ES25.17E3)") dc(i, j)
            end do
        end do
        do k = 1, size(dn, 3)
            do j = 1, size(dn, 2)
                do i = 1, size(dn, 1)
                    write(*, "(ES25.17E3)") dn(i, j, k)
                end do
            end do
        end do
    end subroutine dump_arrays

end program dump_descriptor_reference
