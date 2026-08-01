program test_tio2
    use iso_fortran_env, only: real64
    use accelnet_descriptors
    implicit none

    character(len=1024) :: filename
    character(len=2), parameter :: names(2) = ["Ti", "O "]
    type(descriptor_config) :: config
    type(atomic_structure) :: structure
    type(neighbor_data) :: neighbors
    real(real64), allocatable :: values(:, :), values2(:, :)
    integer :: i

    if (command_argument_count() /= 1) error stop "test_tio2 requires one XSF"
    call get_command_argument(1, filename)
    call read_xsf(trim(filename), names, structure)
    if (structure%natoms /= 24) error stop "unexpected fixture atom count"
    if (count(structure%species == 1) /= 8) error stop "unexpected Ti count"
    if (count(structure%species == 2) /= 16) error stop "unexpected O count"

    call initialize_config(config, 2, 6.5_real64, 20, 5.0_real64, 6, version=0)
    if (config%num_descriptors() /= 56) error stop "unexpected descriptor count"
    call build_neighbor_list(structure, 6.5_real64, neighbors)
    if (any([(neighbors%count_for_atom(i), i=1,structure%natoms)] <= 0)) &
        error stop "empty periodic neighbor list"

    allocate(values(56, structure%natoms), values2(56, structure%natoms))
    call evaluate_structure(config, structure, neighbors, values)
    call evaluate_structure(config, structure, neighbors, values2)
    if (any(values /= values2)) error stop "descriptor is not deterministic"
    if (any(.not. (abs(values) < huge(1.0_real64)))) error stop "non-finite descriptor"
    if (sum(abs(values)) <= 0.0_real64) error stop "zero descriptor"

    write(*, "(A,I0,A,I0)") "TiO2 regression fixture passed: ", structure%natoms, &
                            " atoms, ", size(neighbors%atom_indices)
end program test_tio2
