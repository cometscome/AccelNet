program test_setup_parser
    use iso_fortran_env, only: real64
    use accelnet_descriptors
    use accelnet_behler
    use accelnet_descriptor_models
    use accelnet_setup
    implicit none

    character(len=1024) :: chebyshev_file, lj_file, behler_file, multi_file
    character(len=1024) :: version1_file, version10_file
    character(len=2), parameter :: global_species(2) = ["Ti", "O "]
    type(descriptor_setup) :: chebyshev_setup, lj_setup, behler_setup, multi_setup
    type(descriptor_setup) :: version1_setup, version10_setup
    type(descriptor_config) :: manual_chebyshev
    type(behler_config) :: manual_behler
    type(descriptor_model) :: manual_model
    real(real64) :: displacements(3, 3)
    integer :: global_neighbors(3), local_neighbors(3)
    real(real64), allocatable :: parsed_values(:), manual_values(:)
    real(real64) :: parsed_chebyshev(56), manual_chebyshev_values(56)

    if (command_argument_count() /= 6) error stop "usage: test_setup_parser CHEB LJ BEHLER MULTI V1 V10"
    call get_command_argument(1, chebyshev_file)
    call get_command_argument(2, lj_file)
    call get_command_argument(3, behler_file)
    call get_command_argument(4, multi_file)
    call get_command_argument(5, version1_file)
    call get_command_argument(6, version10_file)
    call read_accelnet_setup(trim(chebyshev_file), global_species, chebyshev_setup)
    call read_accelnet_setup(trim(lj_file), global_species, lj_setup)
    call read_accelnet_setup(trim(behler_file), global_species, behler_setup)
    call read_accelnet_setup(trim(multi_file), global_species, multi_setup)
    call read_accelnet_setup(trim(version1_file), global_species, version1_setup)
    call read_accelnet_setup(trim(version10_file), global_species, version10_setup)

    if (chebyshev_setup%model%num_descriptors() /= 56) error stop "parsed Chebyshev size"
    if (lj_setup%model%num_descriptors() /= 4) error stop "parsed LJ size"
    if (behler_setup%model%num_descriptors() /= 5) error stop "parsed Behler size"
    if (multi_setup%model%num_descriptors() /= 60) error stop "parsed multi size"
    if (version1_setup%model%num_descriptors() /= 56) error stop "parsed version 1 size"
    if (version10_setup%model%num_descriptors() /= 56) error stop "parsed version 10 size"
    if (trim(chebyshev_setup%central_species) /= "Ti") error stop "parsed central species"
    if (chebyshev_setup%central_global_species /= 1) error stop "central global mapping"
    if (abs(chebyshev_setup%minimum_distance - 0.35_real64) > 0.0_real64) error stop "parsed RMIN"
    if (index(chebyshev_setup%description, "Descriptor for Ti") == 0) error stop "parsed description"

    global_neighbors = [2, 1, 2]
    call behler_setup%map_species(global_neighbors, local_neighbors)
    if (any(local_neighbors /= [2, 1, 2])) error stop "parsed species mapping"
    displacements(:, 1) = [1.1_real64, 0.2_real64, -0.1_real64]
    displacements(:, 2) = [-0.4_real64, 1.3_real64, 0.3_real64]
    displacements(:, 3) = [0.2_real64, -0.5_real64, 1.5_real64]

    ! AccelNet orders Behler coefficients by species, then G kind, not by file order.
    call initialize_behler_config(manual_behler, 2)
    call add_g1(manual_behler, 1, 4.5_real64)
    call add_g3(manual_behler, 1, 4.5_real64, 1.2_real64)
    call add_g4(manual_behler, 1, 2, 4.5_real64, 1.0_real64, 2.0_real64, 0.15_real64)
    call add_g2(manual_behler, 2, 4.5_real64, 0.3_real64, 0.7_real64)
    call add_g5(manual_behler, 2, 2, 4.5_real64, -1.0_real64, 3.0_real64, 0.2_real64)
    call add_behler(manual_model, manual_behler)
    allocate(parsed_values(5), manual_values(5))
    call evaluate_model_values(behler_setup%model, displacements, local_neighbors, parsed_values)
    call evaluate_model_values(manual_model, displacements, local_neighbors, manual_values)
    if (maxval(abs(parsed_values - manual_values)) > 0.0_real64) then
        write(*, "(A,*(ES16.8,1X))") "parsed: ", parsed_values
        write(*, "(A,*(ES16.8,1X))") "manual: ", manual_values
        error stop "parsed Behler model differs"
    end if

    call initialize_config(manual_chebyshev, 2, 6.5_real64, 20, 5.0_real64, 6, &
                           version=1, central_type_index=1)
    call evaluate_model_values(version1_setup%model, displacements, local_neighbors, parsed_chebyshev)
    call evaluate_atom(manual_chebyshev, displacements, local_neighbors, manual_chebyshev_values)
    if (any(parsed_chebyshev /= manual_chebyshev_values)) error stop "parsed Chebyshev version 1 differs"

    call initialize_config(manual_chebyshev, 2, 6.5_real64, 20, 5.0_real64, 6, &
                           version=10, central_type_index=1)
    call evaluate_model_values(version10_setup%model, displacements, local_neighbors, parsed_chebyshev)
    call evaluate_atom(manual_chebyshev, displacements, local_neighbors, manual_chebyshev_values)
    if (any(parsed_chebyshev /= manual_chebyshev_values)) error stop "parsed Chebyshev version 10 differs"
    write(*, "(A)") "AccelNet setup parser tests passed"
end program test_setup_parser
