program compare_performance
    use iso_fortran_env, only: real64
    use accelnet_descriptors
    use chebyshevbasis_fp, only: FP_ChebyshevBasis
    use lclist, only: lcl_init, lcl_final, lcl_nmax_nbdist, lcl_nbdist_cart
    implicit none

    type(descriptor_config) :: config
    type(descriptor_config) :: config_by_type(2)
    type(atomic_structure) :: structure
    type(neighbor_data) :: neighbors
    type(neighbor_data) :: timed_neighbors
    type(FP_ChebyshevBasis) :: original_basis
    character(len=1024) :: input_file, argument
    character(len=2), parameter :: names(2) = ["Ti", "O "]
    real(real64), allocatable :: parameters(:, :), original_values(:), new_values(:)
    real(real64), allocatable :: absolute_neighbor_positions(:, :)
    real(real64), allocatable, target :: fractional(:, :)
    real(real64), allocatable :: lcl_coordinates(:, :), lcl_distances(:), structure_values(:, :)
    integer, allocatable :: neighbor_species(:)
    integer, allocatable :: lcl_indices(:), lcl_types(:)
    integer :: ltype(2), repetitions, repetition, i, entry, first, last, n, nmax, version
    real(real64) :: t0, t1, original_seconds, new_seconds, original_full_seconds, new_full_seconds
    real(real64) :: original_checksum, new_checksum, ratio
    real(real64) :: center(3)

    if (command_argument_count() < 1) then
        write(*, "(A)") "usage: compare_performance XSF [REPETITIONS] [VERSION]"
        error stop 2
    end if
    call get_command_argument(1, input_file)
    repetitions = 20
    if (command_argument_count() >= 2) then
        call get_command_argument(2, argument)
        read(argument, *) repetitions
    end if
    version = 0
    if (command_argument_count() >= 3) then
        call get_command_argument(3, argument)
        read(argument, *) version
    end if

    call initialize_config(config_by_type(1), 2, 6.5_real64, 20, 5.0_real64, 6, &
                           version=version, central_type_index=1)
    call initialize_config(config_by_type(2), 2, 6.5_real64, 20, 5.0_real64, 6, &
                           version=version, central_type_index=2)
    config = config_by_type(1)
    call read_xsf(trim(input_file), names, structure)
    call build_neighbor_list(structure, 6.5_real64, neighbors, preserve_legacy_order=version == 10)

    allocate(parameters(5, config%num_descriptors()), source=0.0_real64)
    parameters(:, 1) = [6.5_real64, 20.0_real64, 5.0_real64, 6.0_real64, real(version, real64)]
    original_basis = FP_ChebyshevBasis(parameters, 2, names)
    ltype = [1, 2]
    allocate(original_values(config%num_descriptors()), new_values(config%num_descriptors()))
    allocate(absolute_neighbor_positions(3, size(neighbors%atom_indices)))
    allocate(neighbor_species(size(neighbors%atom_indices)))
    allocate(fractional(3, structure%natoms))
    fractional = matmul(inverse_matrix(structure%lattice), structure%positions)
    nmax = lcl_nmax_nbdist(1.0_real64, 6.5_real64)
    allocate(lcl_coordinates(3, nmax), lcl_distances(nmax), lcl_indices(nmax), lcl_types(nmax))
    allocate(structure_values(config%num_descriptors(), structure%natoms))
    do i = 1, structure%natoms
        do entry = neighbors%offsets(i), neighbors%offsets(i + 1) - 1
            absolute_neighbor_positions(:, entry) = neighbors%positions(:, entry)
            neighbor_species(entry) = structure%species(neighbors%atom_indices(entry))
        end do
    end do

    ! Warm both implementations and force their results to be observed.
    call run_original(original_checksum)
    call run_new(new_checksum)
    if (abs(original_checksum - new_checksum) > 5.0e-12_real64*max(1.0_real64, abs(original_checksum))) &
        error stop "warm-up checksums differ"

    call cpu_time(t0)
    do repetition = 1, repetitions
        call run_original(original_checksum)
    end do
    call cpu_time(t1)
    original_seconds = t1 - t0

    call cpu_time(t0)
    do repetition = 1, repetitions
        call run_new(new_checksum)
    end do
    call cpu_time(t1)
    new_seconds = t1 - t0

    if (abs(original_checksum - new_checksum) > 5.0e-12_real64*max(1.0_real64, abs(original_checksum))) &
        error stop "benchmark checksums differ"
    ratio = new_seconds/original_seconds
    write(*, "(A,I0)") "atoms: ", structure%natoms
    write(*, "(A,I0)") "Chebyshev version: ", version
    write(*, "(A,I0)") "neighbor entries: ", size(neighbors%atom_indices)
    write(*, "(A,F12.6)") "original AccelNet seconds: ", original_seconds
    write(*, "(A,F12.6)") "AccelNetDescriptors seconds: ", new_seconds
    write(*, "(A,F10.4)") "new/original ratio: ", ratio
    write(*, "(A,ES24.16)") "checksum: ", new_checksum
    if (ratio > 1.05_real64) error stop "performance regression exceeds 5 percent"

    ! Include neighbor-list construction and the high-level structure API.
    call run_original_full(original_checksum)
    call run_new_full(new_checksum)
    if (abs(original_checksum - new_checksum) > 5.0e-13_real64*max(1.0_real64, abs(original_checksum))) &
        error stop "end-to-end warm-up checksums differ"
    call cpu_time(t0)
    do repetition = 1, repetitions
        call run_original_full(original_checksum)
    end do
    call cpu_time(t1)
    original_full_seconds = t1 - t0
    call cpu_time(t0)
    do repetition = 1, repetitions
        call run_new_full(new_checksum)
    end do
    call cpu_time(t1)
    new_full_seconds = t1 - t0
    ratio = new_full_seconds/original_full_seconds
    write(*, "(A,F12.6)") "original end-to-end seconds: ", original_full_seconds
    write(*, "(A,F12.6)") "new end-to-end seconds: ", new_full_seconds
    write(*, "(A,F10.4)") "end-to-end new/original ratio: ", ratio
    if (ratio > 1.05_real64) error stop "end-to-end performance regression exceeds 5 percent"

contains

    subroutine run_original(checksum)
        real(real64), intent(out) :: checksum
        checksum = 0.0_real64
        do i = 1, structure%natoms
            first = neighbors%offsets(i)
            last = neighbors%offsets(i + 1) - 1
            n = last - first + 1
            call original_basis%evaluate(structure%species(i), structure%positions(:, i), n, &
                absolute_neighbor_positions(:, first:last), neighbor_species(first:last), &
                ltype, original_values)
            checksum = checksum + sum(original_values)
        end do
    end subroutine run_original

    subroutine run_new(checksum)
        real(real64), intent(out) :: checksum
        checksum = 0.0_real64
        do i = 1, structure%natoms
            first = neighbors%offsets(i)
            last = neighbors%offsets(i + 1) - 1
            call evaluate_atom(config_by_type(structure%species(i)), neighbors%displacements(:, first:last), &
                               neighbor_species(first:last), new_values)
            checksum = checksum + sum(new_values)
        end do
    end subroutine run_new

    subroutine run_original_full(checksum)
        real(real64), intent(out) :: checksum
        checksum = 0.0_real64
        call lcl_init(1.0_real64, 6.5_real64, structure%lattice, structure%natoms, &
                      structure%species, fractional, structure%pbc)
        do i = 1, structure%natoms
            n = nmax
            call lcl_nbdist_cart(i, n, lcl_coordinates, lcl_distances, r_cut=6.5_real64, &
                                 nblist=lcl_indices, nbtype=lcl_types)
            center = matmul(structure%lattice, fractional(:, i))
            call original_basis%evaluate(structure%species(i), center, n, &
                lcl_coordinates(:, 1:n), lcl_types(1:n), ltype, original_values)
            checksum = checksum + sum(original_values)
        end do
        call lcl_final()
    end subroutine run_original_full

    subroutine run_new_full(checksum)
        real(real64), intent(out) :: checksum
        call build_neighbor_list(structure, 6.5_real64, timed_neighbors, preserve_legacy_order=version == 10)
        call evaluate_structure(config, structure, timed_neighbors, structure_values)
        checksum = sum(structure_values)
    end subroutine run_new_full

    pure function inverse_matrix(matrix) result(inverse)
        real(real64), intent(in) :: matrix(3, 3)
        real(real64) :: inverse(3, 3), determinant

        determinant = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) &
                    - matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) &
                    + matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
        inverse(1,1) =  (matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))/determinant
        inverse(1,2) = -(matrix(1,2)*matrix(3,3)-matrix(1,3)*matrix(3,2))/determinant
        inverse(1,3) =  (matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2))/determinant
        inverse(2,1) = -(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))/determinant
        inverse(2,2) =  (matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1))/determinant
        inverse(2,3) = -(matrix(1,1)*matrix(2,3)-matrix(1,3)*matrix(2,1))/determinant
        inverse(3,1) =  (matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))/determinant
        inverse(3,2) = -(matrix(1,1)*matrix(3,2)-matrix(1,2)*matrix(3,1))/determinant
        inverse(3,3) =  (matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1))/determinant
    end function inverse_matrix

end program compare_performance
