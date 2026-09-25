program test_g4_groups
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: ieee_exceptions, only: ieee_set_flag, ieee_get_flag, ieee_invalid, ieee_divide_by_zero
    use accelnet_behler
    use accelnet_descriptors, only: cutoff_value
    implicit none
    integer, parameter :: ng = 8, nn = 6
    integer, parameter :: first(ng) = [1, 2, 1, 1, 2, 3, 1, 1], second(ng) = [2, 2, 2, 1, 2, 3, 2, 2]
    real(real64), parameter :: rc(ng) = &
        [4.5_real64, 3.7_real64, 4.5_real64, 2.9_real64, 3.7_real64, 5.2_real64, 3.7_real64, 4.5_real64]
    real(real64), parameter :: eta(ng) = &
        [0.15_real64, 0.31_real64, 0.15_real64, 0.0_real64, 0.31_real64, 0.42_real64, 0.31_real64, 0.31_real64]
    real(real64), parameter :: rs(ng) = &
        [0.0_real64, 0.4_real64, 0.0_real64, 0.0_real64, 0.4_real64, -0.2_real64, 0.0_real64, 0.4_real64]
    real(real64), parameter :: zeta(ng) = &
        [1.0_real64, 4.0_real64, 2.5_real64, 1.5_real64, 1.0_real64, 3.0_real64, 3.0_real64, 2.0_real64]
    real(real64), parameter :: lambda(ng) = &
        [1.0_real64, -1.0_real64, -1.0_real64, 1.0_real64, 1.0_real64, -1.0_real64, 1.0_real64, -1.0_real64]
    real(real64), parameter :: step = 2.0e-6_real64
    integer, parameter :: grouped_order(ng) = [1, 3, 7, 8, 2, 5, 4, 6]
    type(behler_config) :: config
    integer :: species(nn), output(ng), kind, layout, p, index, j, axis, scenario, permutation(nn), shift_case
    real(real64) :: shifts(ng)
    logical :: invalid, divided_by_zero
    real(real64) :: xyz(3, nn), moved(3, nn), alpha, reference, numerical, wide_difference, error
    real(real64), allocatable :: values(:), only_values(:), center(:, :), neighbors(:, :, :), reordered(:)

    ! Test mixed zero/nonzero shifts, all-fast and all-general configurations.
    ! The mixed case has two zero-shift cutoffs for pair (1,2), as well as a
    ! shifted descriptor at the same cutoff and a negative-shift species pair.
    do shift_case = 0, 2
        shifts = rs
        if (shift_case == 1) shifts = 0.0_real64
        if (shift_case == 2) then
            where (shifts == 0.0_real64) shifts = -0.15_real64
        end if
        ! Layout 0 has contiguous outputs within each pair. Layout 1 interleaves
        ! other descriptor families and extra radial/cutoff groups unused by G4.
        do layout = 0, 1
            do kind = 0, 9
                alpha = 0.25_real64
                call initialize_behler_config(config, 3, kind, alpha)
                do index = 1, ng
                    p = index
                    if (layout == 0) p = grouped_order(index)
                    if (layout == 1) then
                        call add_g2(config, 1, 6.1_real64, 0.3_real64, 0.27_real64)
                        call add_g5(config, 1, 2, 5.8_real64, 1.0_real64, 2.0_real64, 0.73_real64)
                    end if
                    call add_g4(config, first(p), second(p), rc(p), lambda(p), zeta(p), eta(p), shifts(p))
                    output(p) = config%num_descriptors()
                end do
                allocate(values(config%num_descriptors()), only_values(config%num_descriptors()), &
                    reordered(config%num_descriptors()), center(3, config%num_descriptors()), &
                    neighbors(3, config%num_descriptors(), nn))
                do scenario = 1, 4
                    xyz(:, 1) = [1.1_real64, 0.2_real64, -0.1_real64]
                    xyz(:, 2) = [-0.4_real64, 1.3_real64, 0.3_real64]
                    xyz(:, 3) = [0.2_real64, -0.5_real64, 1.5_real64]
                    xyz(:, 4) = [2.7_real64, 0.1_real64, 0.2_real64]
                    xyz(:, 5) = [-1.3_real64, -0.7_real64, 0.4_real64]
                    xyz(:, 6) = [0.5_real64, 0.3_real64, -1.2_real64]
                    species = [1, 2, 1, 2, 2, 3]
                    if (scenario == 2) then
                        ! Almost collinear bonds, including a zero angular base.
                        xyz(:, 1) = [1.0_real64, 0.0_real64, 0.0_real64]
                        xyz(:, 2) = [1.4_real64, 1.0e-7_real64, 0.0_real64]
                        xyz(:, 3) = [-1.2_real64, 0.0_real64, 0.0_real64]
                    else if (scenario == 3) then
                        ! Empty species groups plus a bond exactly at a cutoff.
                        species = [1, 1, 1, 1, 1, 1]
                        xyz(:, 4) = [2.9_real64, 0.0_real64, 0.0_real64]
                    else if (scenario == 4) then
                        ! Activate the negative-shift (3,3) pair.
                        species(3) = 3
                    end if
                    call evaluate_behler_values_derivatives(config, xyz, species, values, center, neighbors)
                    call evaluate_behler_values(config, xyz, species, only_values)
                    if (.not. all(ieee_is_finite(values)) .or. .not. all(ieee_is_finite(neighbors))) &
                        error stop 'non-finite grouped G4'
                    if (maxval(abs(values - only_values)) > 1.0e-12_real64) error stop 'G4 value-only mismatch'
                    if (any(values(output) /= only_values(output))) error stop 'G4 value paths differ in rounding'
                    if (maxval(abs(center + sum(neighbors, dim=3))) > 2.0e-12_real64) &
                        error stop 'G4 translation invariance'
                    permutation = [6, 4, 2, 5, 1, 3]
                    call evaluate_behler_values(config, xyz(:, permutation), species(permutation), reordered)
                    if (maxval(abs(reordered - values)) > 2.0e-12_real64) error stop 'G4 permutation invariance'
                    do p = 1, ng
                        reference = scalar_g4(xyz, species, p, kind, alpha)
                        if (abs(values(output(p)) - reference) > 2.0e-12_real64) error stop 'G4 scalar definition mismatch'
                        ! A hard cutoff has no derivative at its discontinuity.
                        if (kind == 0 .and. scenario == 3) cycle
                        do j = 1, nn
                            do axis = 1, 3
                                moved = xyz
                                moved(axis, j) = moved(axis, j) + step
                                numerical = scalar_g4(moved, species, p, kind, alpha)
                                moved(axis, j) = moved(axis, j) - 2.0_real64*step
                                numerical = (numerical - scalar_g4(moved, species, p, kind, alpha))/(2.0_real64*step)
                                if (scenario == 3) then
                                    ! A smooth cutoff can have a discontinuous second
                                    ! derivative at Rc. Cancel the O(h) boundary error.
                                    moved = xyz
                                    moved(axis, j) = moved(axis, j) + 2.0_real64*step
                                    wide_difference = scalar_g4(moved, species, p, kind, alpha)
                                    moved(axis, j) = moved(axis, j) - 4.0_real64*step
                                    wide_difference = (wide_difference - &
                                        scalar_g4(moved, species, p, kind, alpha))/(4.0_real64*step)
                                    numerical = 2.0_real64*numerical - wide_difference
                                end if
                                error = abs(numerical - neighbors(axis, output(p), j))
                                if (error > 2.0e-7_real64) then
                                    print *, 'G4 derivative error:', layout, kind, scenario, p, j, axis, error
                                    error stop 'G4 scalar finite difference mismatch'
                                end if
                            end do
                        end do
                    end do
                end do
                deallocate(values, only_values, reordered, center, neighbors)
            end do
        end do
    end do
    ! Energy-only calls must not evaluate an unused singular angular derivative.
    call initialize_behler_config(config, 1)
    call add_g4(config, 1, 1, 4.5_real64, -1.0_real64, 0.5_real64, 0.15_real64)
    call add_g4(config, 1, 1, 4.5_real64, -1.0_real64, 0.0_real64, 0.15_real64)
    allocate(values(2))
    xyz(:, 1) = [1.0_real64, 0.0_real64, 0.0_real64]
    xyz(:, 2) = [2.0_real64, 0.0_real64, 0.0_real64]
    species = 1
    call ieee_set_flag(ieee_invalid, .false.)
    call ieee_set_flag(ieee_divide_by_zero, .false.)
    call evaluate_behler_values(config, xyz(:, 1:2), species(1:2), values)
    call ieee_get_flag(ieee_invalid, invalid)
    call ieee_get_flag(ieee_divide_by_zero, divided_by_zero)
    if (invalid .or. divided_by_zero .or. .not. all(ieee_is_finite(values))) &
        error stop 'G4 energy-only call evaluated a singular derivative'
    if (values(1) /= 0.0_real64 .or. values(2) <= 0.0_real64) error stop 'G4 energy-only collinear value'
    print *, 'Grouped G4 values, finite differences and invariants passed.'
contains
    function scalar_g4(r, species, p, kind, alpha) result(value)
        real(real64), intent(in) :: r(:, :), alpha
        integer, intent(in) :: species(:), p, kind
        real(real64) :: value, a, b, c, cosine, angular
        integer :: j, k
        value = 0.0_real64
        ! Deliberately evaluate ordered pairs and the defining scalar formula,
        ! without the production kernel's group metadata or derivative algebra.
        do j = 1, size(species)
            do k = 1, size(species)
                if (j == k) cycle
                if (.not. ((species(j) == first(p) .and. species(k) == second(p)) .or. &
                    (species(k) == first(p) .and. species(j) == second(p)))) cycle
                a = norm2(r(:, j))
                b = norm2(r(:, k))
                c = norm2(r(:, k) - r(:, j))
                if (min(a, b, c) <= 1.0e-12_real64 .or. max(a, b, c) >= rc(p)) cycle
                cosine = max(-1.0_real64, min(1.0_real64, dot_product(r(:, j), r(:, k))/(a*b)))
                angular = 2.0_real64**(-zeta(p))*max(0.0_real64, 1.0_real64 + lambda(p)*cosine)**zeta(p)
                value = value + angular* &
                    exp(-eta(p)*((a - shifts(p))**2 + (b - shifts(p))**2 + (c - shifts(p))**2))* &
                    cutoff_value(a, rc(p), kind, alpha)*cutoff_value(b, rc(p), kind, alpha)* &
                    cutoff_value(c, rc(p), kind, alpha)
            end do
        end do
    end function scalar_g4
end program test_g4_groups
