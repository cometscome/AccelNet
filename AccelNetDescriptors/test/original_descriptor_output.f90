program original_descriptor_output
    use iso_fortran_env, only: real64, int64
    use accelnet_descriptors
    use chebyshevbasis_fp, only: FP_ChebyshevBasis
    use lclist, only: lcl_init, lcl_final, lcl_nmax_nbdist, lcl_nbdist_cart
    implicit none

    type(descriptor_config) :: config
    type(atomic_structure), allocatable :: structures(:)
    type(neighbor_data) :: neighbors
    type(FP_ChebyshevBasis) :: original_basis
    character(len=1024) :: output_file, input_file
    character(len=2), parameter :: names(2) = ["Ti", "O "]
    real(real64), allocatable :: parameters(:, :), values(:), coordinates(:, :)
    real(real64), allocatable, target :: fractional(:, :)
    real(real64), allocatable :: distances(:)
    integer, allocatable :: types(:)
    integer, allocatable :: atom_indices(:)
    integer :: unit, nargs, s, i, n, nmax, g, version, input_start, number_of_structures
    integer :: ltype(2)
    real(real64) :: center(3)

    nargs = command_argument_count()
    if (nargs < 2) then
        write(*, "(A)") "usage: original_descriptor_output OUTPUT [--version=N] XSF [XSF ...]"
        error stop 2
    end if
    call get_command_argument(1, output_file)
    version = 0
    input_start = 2
    call get_command_argument(2, input_file)
    if (index(input_file, "--version=") == 1) then
        read(input_file(11:), *) version
        input_start = 3
    end if
    number_of_structures = nargs - input_start + 1
    if (number_of_structures < 1) error stop "at least one XSF input is required"

    call initialize_config(config, 2, 6.5_real64, 20, 5.0_real64, 6, &
                           version=version, central_type_index=1)
    allocate(structures(number_of_structures))
    do s = 1, number_of_structures
        call get_command_argument(input_start + s - 1, input_file)
        call read_xsf(trim(input_file), names, structures(s))
    end do
    allocate(parameters(5, config%num_descriptors()), source=0.0_real64)
    parameters(:, 1) = [6.5_real64, 20.0_real64, 5.0_real64, 6.0_real64, real(version, real64)]
    original_basis = FP_ChebyshevBasis(parameters, 2, names)
    ltype = [1, 2]
    allocate(values(config%num_descriptors()))

    open(newunit=unit, file=trim(output_file), status="replace", action="write")
    write(unit, "(A)") "ACCELNET_DESCRIPTOR_HEX_V1"
    write(unit, "(3(I0,1X))") size(structures), 2, config%num_descriptors()
    do s = 1, size(structures)
        call get_command_argument(input_start + s - 1, input_file)
        write(unit, "(A)") trim(input_file)
        write(unit, "(I0)") structures(s)%natoms
        allocate(fractional(3, structures(s)%natoms))
        fractional = matmul(inverse_matrix(structures(s)%lattice), structures(s)%positions)
        call lcl_init(1.0_real64, 6.5_real64, structures(s)%lattice, structures(s)%natoms, &
                      structures(s)%species, fractional, structures(s)%pbc)
        nmax = lcl_nmax_nbdist(1.0_real64, 6.5_real64)
        allocate(coordinates(3, nmax), distances(nmax), types(nmax), atom_indices(nmax))
        do i = 1, structures(s)%natoms
            n = nmax
            call lcl_nbdist_cart(i, n, coordinates, distances, r_cut=6.5_real64, &
                                 nblist=atom_indices, nbtype=types)
            center = matmul(structures(s)%lattice, fractional(:, i))
            call original_basis%evaluate(structures(s)%species(i), center, n, &
                                         coordinates(:, 1:n), types(1:n), ltype, values)
            write(unit, "(I0,1X,I0)") i, structures(s)%species(i)
            write(unit, "(*(Z16.16,1X))") (transfer(values(g), 0_int64), g = 1, config%num_descriptors())
        end do
        deallocate(coordinates, distances, types, atom_indices)
        call lcl_final()
        deallocate(fractional)
    end do
    close(unit)

contains

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
end program original_descriptor_output
