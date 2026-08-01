program compare_extension_performance
    use iso_fortran_env, only: real64
    use accelnet_descriptors
    use accelnet_lj
    use accelnet_descriptor_models
    use accelnet_behler
    use LJbasis_fp, only: FP_LJBasis
    use behler2011basis_fp, only: FP_Behler2011basis
    implicit none

    type(descriptor_config) :: chebyshev
    type(lj_config) :: lj
    type(descriptor_model) :: chebyshev_model
    type(atomic_structure) :: structure
    type(neighbor_data) :: neighbors
    type(FP_LJBasis) :: original_lj
    type(behler_config) :: behler
    type(FP_Behler2011basis) :: original_behler
    character(len=1024) :: input_file, argument
    character(len=2), parameter :: names(2) = ["Ti", "O "]
    real(real64), allocatable :: parameters(:, :), direct_values(:), model_values(:)
    real(real64), allocatable :: original_lj_values(:), new_lj_values(:)
    real(real64), allocatable :: original_behler_values(:), new_behler_values(:)
    real(real64) :: behler_parameters(5, 5)
    integer :: behler_kinds(5), behler_environments(2, 5)
    integer, allocatable :: all_species(:)
    integer :: repetitions, repetition, atom, entry, first, last, n, ltype(2)
    real(real64) :: t0, t1, direct_seconds, model_seconds, original_lj_seconds, new_lj_seconds
    real(real64) :: original_behler_seconds, new_behler_seconds
    real(real64) :: direct_checksum, model_checksum, original_lj_checksum, new_lj_checksum, ratio

    if (command_argument_count() < 1) error stop "usage: compare_extension_performance XSF [REPETITIONS]"
    call get_command_argument(1, input_file)
    repetitions = 200
    if (command_argument_count() >= 2) then
        call get_command_argument(2, argument)
        read(argument, *) repetitions
    end if

    call read_xsf(trim(input_file), names, structure)
    call build_neighbor_list(structure, 6.5_real64, neighbors)
    call initialize_config(chebyshev, 2, 6.5_real64, 20, 5.0_real64, 6, version=0)
    call initialize_lj_config(lj, 2, 6.5_real64)
    call initialize_behler_config(behler, 2)
    call add_g1(behler, 1, 4.5_real64)
    call add_g2(behler, 2, 4.5_real64, 0.3_real64, 0.7_real64)
    call add_g3(behler, 1, 4.5_real64, 1.2_real64)
    call add_g4(behler, 1, 2, 4.5_real64, 1.0_real64, 2.0_real64, 0.15_real64)
    call add_g5(behler, 2, 2, 4.5_real64, -1.0_real64, 3.0_real64, 0.2_real64)
    call add_chebyshev(chebyshev_model, chebyshev)
    allocate(direct_values(chebyshev%num_descriptors()))
    allocate(model_values(chebyshev_model%num_descriptors()))
    allocate(original_lj_values(lj%num_descriptors()), new_lj_values(lj%num_descriptors()))
    allocate(original_behler_values(behler%num_descriptors()), new_behler_values(behler%num_descriptors()))
    allocate(all_species(size(neighbors%atom_indices)))
    do entry = 1, size(neighbors%atom_indices)
        all_species(entry) = structure%species(neighbors%atom_indices(entry))
    end do

    allocate(parameters(5, lj%num_descriptors()), source=0.0_real64)
    parameters(1, 1) = 6.5_real64
    original_lj = FP_LJBasis(parameters, 2, names)
    behler_parameters = 0.0_real64
    behler_parameters(:, 1) = [4.5_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    behler_parameters(:, 2) = [4.5_real64, 0.3_real64, 0.7_real64, 0.0_real64, 0.0_real64]
    behler_parameters(:, 3) = [4.5_real64, 1.2_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    behler_parameters(:, 4) = [4.5_real64, 1.0_real64, 2.0_real64, 0.15_real64, 0.0_real64]
    behler_parameters(:, 5) = [4.5_real64, -1.0_real64, 3.0_real64, 0.2_real64, 0.0_real64]
    behler_kinds = [1, 2, 3, 4, 5]
    behler_environments = 0
    behler_environments(1, :) = [1, 2, 1, 1, 2]
    behler_environments(2, 4:5) = [2, 2]
    original_behler = FP_Behler2011basis(5, 2, 5, behler_parameters, 2, names, &
                                         behler_kinds, behler_environments)
    ltype = [1, 2]

    call run_chebyshev_direct(direct_checksum)
    call run_chebyshev_model(model_checksum)
    if (direct_checksum /= model_checksum) error stop "model changes Chebyshev values"
    call cpu_time(t0)
    do repetition = 1, repetitions
        call run_chebyshev_direct(direct_checksum)
    end do
    call cpu_time(t1)
    direct_seconds = t1 - t0
    call cpu_time(t0)
    do repetition = 1, repetitions
        call run_chebyshev_model(model_checksum)
    end do
    call cpu_time(t1)
    model_seconds = t1 - t0
    ratio = model_seconds/direct_seconds
    write(*, "(A,F12.6)") "direct Chebyshev seconds: ", direct_seconds
    write(*, "(A,F12.6)") "model Chebyshev seconds: ", model_seconds
    write(*, "(A,F10.4)") "model/direct Chebyshev ratio: ", ratio
    if (ratio > 1.03_real64) error stop "common model regresses Chebyshev by more than 3 percent"

    call run_original_lj(original_lj_checksum)
    call run_new_lj(new_lj_checksum)
    if (abs(original_lj_checksum - new_lj_checksum) > &
        5.0e-13_real64*max(1.0_real64, abs(original_lj_checksum))) error stop "LJ checksums differ"
    call cpu_time(t0)
    do repetition = 1, repetitions
        call run_original_lj(original_lj_checksum)
    end do
    call cpu_time(t1)
    original_lj_seconds = t1 - t0
    call cpu_time(t0)
    do repetition = 1, repetitions
        call run_new_lj(new_lj_checksum)
    end do
    call cpu_time(t1)
    new_lj_seconds = t1 - t0
    ratio = new_lj_seconds/original_lj_seconds
    write(*, "(A,F12.6)") "original LJ seconds: ", original_lj_seconds
    write(*, "(A,F12.6)") "new LJ seconds: ", new_lj_seconds
    write(*, "(A,F10.4)") "new/original LJ ratio: ", ratio
    if (ratio > 1.05_real64) error stop "LJ performance regression exceeds 5 percent"

    call run_original_behler(original_lj_checksum)
    call run_new_behler(new_lj_checksum)
    if (abs(original_lj_checksum - new_lj_checksum) > &
        2.0e-12_real64*max(1.0_real64, abs(original_lj_checksum))) error stop "Behler checksums differ"
    call cpu_time(t0)
    do repetition = 1, repetitions
        call run_original_behler(original_lj_checksum)
    end do
    call cpu_time(t1)
    original_behler_seconds = t1 - t0
    call cpu_time(t0)
    do repetition = 1, repetitions
        call run_new_behler(new_lj_checksum)
    end do
    call cpu_time(t1)
    new_behler_seconds = t1 - t0
    ratio = new_behler_seconds/original_behler_seconds
    write(*, "(A,F12.6)") "original Behler seconds: ", original_behler_seconds
    write(*, "(A,F12.6)") "new Behler seconds: ", new_behler_seconds
    write(*, "(A,F10.4)") "new/original Behler ratio: ", ratio
    if (ratio > 1.05_real64) error stop "Behler performance regression exceeds 5 percent"

contains

    subroutine run_chebyshev_direct(checksum)
        real(real64), intent(out) :: checksum
        checksum = 0.0_real64
        do atom = 1, structure%natoms
            first = neighbors%offsets(atom)
            last = neighbors%offsets(atom + 1) - 1
            call evaluate_atom(chebyshev, neighbors%displacements(:, first:last), &
                               all_species(first:last), direct_values)
            checksum = checksum + sum(direct_values)
        end do
    end subroutine run_chebyshev_direct

    subroutine run_chebyshev_model(checksum)
        real(real64), intent(out) :: checksum
        checksum = 0.0_real64
        do atom = 1, structure%natoms
            first = neighbors%offsets(atom)
            last = neighbors%offsets(atom + 1) - 1
            call evaluate_model_values(chebyshev_model, neighbors%displacements(:, first:last), &
                                       all_species(first:last), model_values)
            checksum = checksum + sum(model_values)
        end do
    end subroutine run_chebyshev_model

    subroutine run_original_lj(checksum)
        real(real64), intent(out) :: checksum
        checksum = 0.0_real64
        do atom = 1, structure%natoms
            first = neighbors%offsets(atom)
            last = neighbors%offsets(atom + 1) - 1
            n = last - first + 1
            original_lj_values = 0.0_real64
            call original_lj%evaluate(structure%species(atom), structure%positions(:, atom), n, &
                                      neighbors%positions(:, first:last), all_species(first:last), &
                                      ltype, original_lj_values)
            checksum = checksum + sum(original_lj_values)
        end do
    end subroutine run_original_lj

    subroutine run_new_lj(checksum)
        real(real64), intent(out) :: checksum
        checksum = 0.0_real64
        do atom = 1, structure%natoms
            first = neighbors%offsets(atom)
            last = neighbors%offsets(atom + 1) - 1
            call evaluate_lj_values(lj, neighbors%displacements(:, first:last), &
                                    all_species(first:last), new_lj_values)
            checksum = checksum + sum(new_lj_values)
        end do
    end subroutine run_new_lj

    subroutine run_original_behler(checksum)
        real(real64), intent(out) :: checksum
        checksum = 0.0_real64
        do atom = 1, structure%natoms
            first = neighbors%offsets(atom)
            last = neighbors%offsets(atom + 1) - 1
            n = last - first + 1
            call original_behler%evaluate(structure%species(atom), structure%positions(:, atom), n, &
                neighbors%positions(:, first:last), all_species(first:last), ltype, original_behler_values)
            checksum = checksum + sum(original_behler_values)
        end do
    end subroutine run_original_behler

    subroutine run_new_behler(checksum)
        real(real64), intent(out) :: checksum
        checksum = 0.0_real64
        do atom = 1, structure%natoms
            first = neighbors%offsets(atom)
            last = neighbors%offsets(atom + 1) - 1
            call evaluate_behler_values(behler, neighbors%displacements(:, first:last), &
                                        all_species(first:last), new_behler_values)
            checksum = checksum + sum(new_behler_values)
        end do
    end subroutine run_new_behler

end program compare_extension_performance
