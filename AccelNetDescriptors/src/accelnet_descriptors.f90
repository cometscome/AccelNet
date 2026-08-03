module accelnet_descriptors
    use iso_fortran_env, only: real64, int64
    use accelnet_legacy_lcl, only: lcl_init, lcl_final, lcl_nmax_nbdist, lcl_nbdist_cart
    implicit none
    private

    real(real64), parameter :: PI_ACCELNET = 3.14159265358979_real64
    real(real64), parameter :: EPS_DISTANCE = 1.0e-12_real64
    real(real64), parameter :: NEIGHBOR_SKIN = 1.0e-3_real64
    integer, parameter :: MOMENT_MIN_ANGULAR_NEIGHBORS = 16
    integer, parameter, public :: CHEBYSHEV_EVALUATION_AUTO = 0
    integer, parameter, public :: CHEBYSHEV_EVALUATION_DIRECT = 1
    integer, parameter, public :: CHEBYSHEV_EVALUATION_MOMENT = 2

    integer, parameter, public :: CUTOFF_HARD = 0
    integer, parameter, public :: CUTOFF_COS = 1
    integer, parameter, public :: CUTOFF_TANHU = 2
    integer, parameter, public :: CUTOFF_TANH = 3
    integer, parameter, public :: CUTOFF_EXP = 4
    integer, parameter, public :: CUTOFF_POLY1 = 5
    integer, parameter, public :: CUTOFF_POLY2 = 6
    integer, parameter, public :: CUTOFF_POLY3 = 7
    integer, parameter, public :: CUTOFF_POLY4 = 8
    integer, parameter, public :: CUTOFF_FRACTIONAL = 9

    type, public :: descriptor_config
        real(real64) :: radial_rc = 0.0_real64
        real(real64) :: angular_rc = 0.0_real64
        integer :: cutoff_type = CUTOFF_COS
        real(real64) :: cutoff_alpha = 0.0_real64
        integer :: radial_order = 0
        integer :: angular_order = 0
        integer :: num_species = 0
        integer :: version = 0
        integer :: evaluation_mode = CHEBYSHEV_EVALUATION_AUTO
        integer :: central_type_index = 0
        real(real64), allocatable :: species_weights(:)
        integer :: number_of_angular_moments = 0
        integer, allocatable :: moment_x_power(:), moment_y_power(:), moment_z_power(:)
        integer, allocatable :: moment_first(:), moment_last(:)
        real(real64), allocatable :: moment_multinomial(:)
        real(real64), allocatable :: angular_power_coefficients(:, :)
    contains
        procedure :: num_descriptors
    end type descriptor_config

    type, public :: atomic_structure
        integer :: natoms = 0
        logical :: pbc = .false.
        real(real64) :: lattice(3, 3) = 0.0_real64
        real(real64), allocatable :: positions(:, :)
        integer, allocatable :: species(:)
    end type atomic_structure

    type, public :: neighbor_data
        integer, allocatable :: offsets(:)
        integer, allocatable :: atom_indices(:)
        integer, allocatable :: image_shifts(:, :)
        real(real64), allocatable :: positions(:, :)
        real(real64), allocatable :: displacements(:, :)
    contains
        procedure :: count_for_atom
    end type neighbor_data

    public :: initialize_config
    public :: read_xsf
    public :: build_neighbor_list
    public :: evaluate_atom
    public :: evaluate_atom_with_derivatives
    public :: contract_atom_derivatives
    public :: evaluate_structure
    public :: write_descriptor_file
    public :: cutoff_value, cutoff_derivative, validate_cutoff_parameters
    public :: chebyshev_values, chebyshev_values_derivatives
    public :: set_chebyshev_evaluation, chebyshev_uses_moments

contains

    subroutine set_chebyshev_evaluation(config, mode)
        type(descriptor_config), intent(inout) :: config
        integer, intent(in) :: mode
        if (mode < CHEBYSHEV_EVALUATION_AUTO .or. mode > CHEBYSHEV_EVALUATION_MOMENT) &
            error stop "Chebyshev evaluation mode must be auto, direct, or moment"
        config%evaluation_mode = mode
    end subroutine set_chebyshev_evaluation

    pure logical function chebyshev_uses_moments(config, angular_neighbors) result(use_moments)
        type(descriptor_config), intent(in) :: config
        integer, intent(in) :: angular_neighbors
        select case (config%evaluation_mode)
        case (CHEBYSHEV_EVALUATION_DIRECT)
            use_moments = .false.
        case (CHEBYSHEV_EVALUATION_MOMENT)
            use_moments = .true.
        case default
            use_moments = angular_neighbors >= MOMENT_MIN_ANGULAR_NEIGHBORS
        end select
    end function chebyshev_uses_moments

    subroutine initialize_config(config, num_species, radial_rc, radial_order, &
                                 angular_rc, angular_order, version, central_type_index, &
                                 cutoff_type, cutoff_alpha)
        type(descriptor_config), intent(out) :: config
        integer, intent(in) :: num_species, radial_order, angular_order
        real(real64), intent(in) :: radial_rc, angular_rc
        integer, intent(in), optional :: version, central_type_index, cutoff_type
        real(real64), intent(in), optional :: cutoff_alpha
        integer :: i, spin

        if (num_species < 1) error stop "num_species must be positive"
        if (radial_rc <= 0.0_real64 .or. angular_rc <= 0.0_real64) &
            error stop "cutoffs must be positive"
        if (radial_order < 0 .or. angular_order < 0) &
            error stop "Chebyshev orders must be non-negative"

        config%num_species = num_species
        config%radial_rc = radial_rc
        config%radial_order = radial_order
        config%angular_rc = angular_rc
        config%angular_order = angular_order
        if (present(cutoff_type)) config%cutoff_type = cutoff_type
        if (present(cutoff_alpha)) config%cutoff_alpha = cutoff_alpha
        call validate_cutoff_parameters(config%cutoff_type, config%cutoff_alpha)
        config%version = 0
        if (present(version)) config%version = version
        if (config%version /= 0 .and. config%version /= 1 .and. config%version /= 10) &
            error stop "Chebyshev version must be 0, 1, or 10"
        config%central_type_index = 0
        if (present(central_type_index)) config%central_type_index = central_type_index
        if (config%version == 10 .and. num_species > 1 .and. config%central_type_index < 1) &
            error stop "Chebyshev version=10 requires central_type_index for multiple species"

        allocate(config%species_weights(num_species))
        spin = -num_species/2
        do i = 1, num_species
            if (spin == 0 .and. mod(num_species, 2) == 0) spin = spin + 1
            config%species_weights(i) = real(spin, real64)
            spin = spin + 1
        end do
        call initialize_angular_moment_basis(config)
    end subroutine initialize_config

    subroutine initialize_angular_moment_basis(config)
        type(descriptor_config), intent(inout) :: config
        integer :: n, q, a, b, c, entry, maximum_moments
        real(real64) :: alpha, beta

        maximum_moments = (config%angular_order + 1)*(config%angular_order + 2)* &
                          (config%angular_order + 3)/6
        allocate(config%moment_x_power(maximum_moments), config%moment_y_power(maximum_moments), &
                 config%moment_z_power(maximum_moments), config%moment_multinomial(maximum_moments))
        allocate(config%moment_first(config%angular_order + 1), &
                 config%moment_last(config%angular_order + 1))
        allocate(config%angular_power_coefficients(config%angular_order + 1, &
                                                   config%angular_order + 1), source=0.0_real64)
        entry = 0
        do q = 0, config%angular_order
            config%moment_first(q + 1) = entry + 1
            do a = 0, q
                do b = 0, q - a
                    c = q - a - b
                    entry = entry + 1
                    config%moment_x_power(entry) = a
                    config%moment_y_power(entry) = b
                    config%moment_z_power(entry) = c
                    config%moment_multinomial(entry) = factorial_real(q)/ &
                        (factorial_real(a)*factorial_real(b)*factorial_real(c))
                end do
            end do
            config%moment_last(q + 1) = entry
        end do
        config%number_of_angular_moments = entry

        alpha = 1.0_real64
        beta = 0.0_real64
        if (config%version == 1) then
            alpha = 2.0_real64/PI_ACCELNET
            beta = -1.0_real64
        end if
        config%angular_power_coefficients(1, 1) = 1.0_real64
        if (config%angular_order >= 1) then
            config%angular_power_coefficients(2, 1) = beta
            config%angular_power_coefficients(2, 2) = alpha
        end if
        do n = 2, config%angular_order
            do q = 0, n - 1
                config%angular_power_coefficients(n + 1, q + 1) = &
                    config%angular_power_coefficients(n + 1, q + 1) + &
                    2.0_real64*beta*config%angular_power_coefficients(n, q + 1)
                config%angular_power_coefficients(n + 1, q + 2) = &
                    config%angular_power_coefficients(n + 1, q + 2) + &
                    2.0_real64*alpha*config%angular_power_coefficients(n, q + 1)
            end do
            do q = 0, n - 2
                config%angular_power_coefficients(n + 1, q + 1) = &
                    config%angular_power_coefficients(n + 1, q + 1) - &
                    config%angular_power_coefficients(n - 1, q + 1)
            end do
        end do
    end subroutine initialize_angular_moment_basis

    pure real(real64) function factorial_real(n) result(value)
        integer, intent(in) :: n
        integer :: i
        value = 1.0_real64
        do i = 2, n
            value = value*real(i, real64)
        end do
    end function factorial_real

    integer function num_descriptors(self) result(n)
        class(descriptor_config), intent(in) :: self
        n = self%radial_order + self%angular_order + 2
        if (self%num_species > 1) n = 2*n
    end function num_descriptors

    integer function count_for_atom(self, iatom) result(n)
        class(neighbor_data), intent(in) :: self
        integer, intent(in) :: iatom
        n = self%offsets(iatom + 1) - self%offsets(iatom)
    end function count_for_atom

    subroutine validate_cutoff_parameters(cutoff_type, alpha)
        integer, intent(in) :: cutoff_type
        real(real64), intent(in) :: alpha
        if (cutoff_type < CUTOFF_HARD .or. cutoff_type > CUTOFF_FRACTIONAL) &
            error stop "cutoff type must be between 0 and 9"
        if (alpha < 0.0_real64 .or. alpha >= 1.0_real64) &
            error stop "cutoff alpha must satisfy 0 <= alpha < 1"
        if (cutoff_type == CUTOFF_FRACTIONAL .and. alpha <= 0.0_real64) &
            error stop "fractional cutoff alpha=h/Rc must be positive"
    end subroutine validate_cutoff_parameters

    pure real(real64) function cutoff_value(distance, rc, cutoff_type, alpha) result(value)
        real(real64), intent(in) :: distance, rc
        integer, intent(in), optional :: cutoff_type
        real(real64), intent(in), optional :: alpha
        integer :: kind
        real(real64) :: inner, inverse_width, x, t, width
        kind = CUTOFF_COS
        if (present(cutoff_type)) kind = cutoff_type
        if (kind == CUTOFF_HARD) then
            value = merge(1.0_real64, 0.0_real64, distance < rc)
            return
        end if
        if (distance >= rc) then
            value = 0.0_real64
            return
        end if
        select case(kind)
        case(CUTOFF_TANHU, CUTOFF_TANH)
            t = tanh(1.0_real64 - distance/rc)
            value = t*t*t
            if (kind == CUTOFF_TANH) value = value/tanh(1.0_real64)**3
        case(CUTOFF_FRACTIONAL)
            width = rc
            if (present(alpha)) width = alpha*rc
            x = (distance - rc)/width
            value = x*x/(1.0_real64 + x*x)
        case default
            inner = 0.0_real64
            if (present(alpha)) inner = alpha*rc
            if (distance < inner) then
                value = 1.0_real64
                return
            end if
            inverse_width = 1.0_real64/(rc - inner)
            x = (distance - inner)*inverse_width
            select case(kind)
            case(CUTOFF_COS)
                value = 0.5_real64*(cos(PI_ACCELNET*x) + 1.0_real64)
            case(CUTOFF_EXP)
                value = exp(1.0_real64 + 1.0_real64/(x*x - 1.0_real64))
            case(CUTOFF_POLY1)
                value = (2.0_real64*x - 3.0_real64)*x*x + 1.0_real64
            case(CUTOFF_POLY2)
                value = ((15.0_real64 - 6.0_real64*x)*x - 10.0_real64)*x*x*x + 1.0_real64
            case(CUTOFF_POLY3)
                value = (x*(x*(20.0_real64*x - 70.0_real64) + 84.0_real64) - 35.0_real64)*x**4 + 1.0_real64
            case(CUTOFF_POLY4)
                value = (x*(x*((315.0_real64 - 70.0_real64*x)*x - 540.0_real64) + &
                    420.0_real64) - 126.0_real64)*x**5 + 1.0_real64
            case default
                value = 0.0_real64
            end select
        end select
    end function cutoff_value

    pure real(real64) function cutoff_derivative(distance, rc, cutoff_type, alpha) result(value)
        real(real64), intent(in) :: distance, rc
        integer, intent(in), optional :: cutoff_type
        real(real64), intent(in), optional :: alpha
        integer :: kind
        real(real64) :: inner, inverse_width, x, t, core_derivative, width
        kind = CUTOFF_COS
        if (present(cutoff_type)) kind = cutoff_type
        if (distance >= rc) then
            value = 0.0_real64
            return
        end if
        select case(kind)
        case(CUTOFF_HARD)
            value = 0.0_real64
        case(CUTOFF_TANHU, CUTOFF_TANH)
            t = tanh(1.0_real64 - distance/rc)
            value = 3.0_real64*t*t*(t*t - 1.0_real64)/rc
            if (kind == CUTOFF_TANH) value = value/tanh(1.0_real64)**3
        case(CUTOFF_FRACTIONAL)
            width = rc
            if (present(alpha)) width = alpha*rc
            x = (distance - rc)/width
            value = 2.0_real64*x/(width*(1.0_real64 + x*x)**2)
        case default
            inner = 0.0_real64
            if (present(alpha)) inner = alpha*rc
            if (distance < inner) then
                value = 0.0_real64
                return
            end if
            inverse_width = 1.0_real64/(rc - inner)
            x = (distance - inner)*inverse_width
            select case(kind)
            case(CUTOFF_COS)
                core_derivative = -0.5_real64*PI_ACCELNET*sin(PI_ACCELNET*x)
            case(CUTOFF_EXP)
                t = 1.0_real64/(x*x - 1.0_real64)
                core_derivative = -2.0_real64*x*t*t*exp(1.0_real64 + t)
            case(CUTOFF_POLY1)
                core_derivative = x*(6.0_real64*x - 6.0_real64)
            case(CUTOFF_POLY2)
                core_derivative = x*x*((60.0_real64 - 30.0_real64*x)*x - 30.0_real64)
            case(CUTOFF_POLY3)
                core_derivative = x**3*(x*(x*(140.0_real64*x - 420.0_real64) + &
                    420.0_real64) - 140.0_real64)
            case(CUTOFF_POLY4)
                core_derivative = x**4*(x*(x*((2520.0_real64 - 630.0_real64*x)*x - &
                    3780.0_real64) + 2520.0_real64) - 630.0_real64)
            case default
                core_derivative = 0.0_real64
            end select
            value = inverse_width*core_derivative
        end select
    end function cutoff_derivative

    pure subroutine chebyshev_values(r, r0, r1, order, values)
        real(real64), intent(in) :: r, r0, r1
        integer, intent(in) :: order
        real(real64), intent(out) :: values(order + 1)
        real(real64) :: x
        integer :: i

        x = (2.0_real64*r - r0 - r1)/(r1 - r0)
        values(1) = 1.0_real64
        if (order > 0) then
            values(2) = x
            do i = 3, order + 1
                values(i) = 2.0_real64*x*values(i - 1) - values(i - 2)
            end do
        end if
    end subroutine chebyshev_values

    pure subroutine chebyshev_values_derivatives(r, r0, r1, order, values, derivatives)
        real(real64), intent(in) :: r, r0, r1
        integer, intent(in) :: order
        real(real64), intent(out) :: values(order + 1), derivatives(order + 1)
        real(real64) :: x, u1, u2, u3, scale
        integer :: i

        call chebyshev_values(r, r0, r1, order, values)
        x = (2.0_real64*r - r0 - r1)/(r1 - r0)
        scale = 2.0_real64/(r1 - r0)
        derivatives(1) = 0.0_real64
        if (order > 0) then
            u1 = 1.0_real64
            derivatives(2) = u1
            u2 = 2.0_real64*x
            do i = 3, order + 1
                derivatives(i) = u2*real(i - 1, real64)
                u3 = 2.0_real64*x*u2 - u1
                u1 = u2
                u2 = u3
            end do
        end if
        derivatives = derivatives*scale
    end subroutine chebyshev_values_derivatives

    subroutine read_xsf(filename, species_names, structure)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: species_names(:)
        type(atomic_structure), intent(out) :: structure
        character(len=2048) :: line
        character(len=16) :: atom_name
        integer :: unit, ios, i, itype, ignored
        real(real64) :: x, y, z
        logical :: found_coordinates

        structure%pbc = .false.
        structure%lattice = 0.0_real64
        found_coordinates = .false.
        open(newunit=unit, file=filename, status="old", action="read")
        do
            read(unit, "(A)", iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (index(line, "PRIMVEC") == 1) then
                structure%pbc = .true.
                read(unit, *) structure%lattice(:, 1)
                read(unit, *) structure%lattice(:, 2)
                read(unit, *) structure%lattice(:, 3)
            else if (index(line, "PRIMCOORD") == 1) then
                read(unit, *) structure%natoms, ignored
                allocate(structure%positions(3, structure%natoms))
                allocate(structure%species(structure%natoms))
                do i = 1, structure%natoms
                    read(unit, "(A)") line
                    read(line, *) atom_name, x, y, z
                    itype = find_species(atom_name, species_names)
                    if (itype == 0) error stop "unknown species in XSF"
                    structure%positions(:, i) = [x, y, z]
                    structure%species(i) = itype
                end do
                found_coordinates = .true.
            end if
        end do
        close(unit)
        if (.not. found_coordinates) error stop "PRIMCOORD not found in XSF"
    end subroutine read_xsf

    integer function find_species(name, species_names) result(index_found)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: species_names(:)
        integer :: i
        index_found = 0
        do i = 1, size(species_names)
            if (trim(name) == trim(species_names(i))) then
                index_found = i
                return
            end if
        end do
    end function find_species

    subroutine build_neighbor_list(structure, cutoff, neighbors, minimum_distance, preserve_legacy_order)
        type(atomic_structure), intent(in) :: structure
        real(real64), intent(in) :: cutoff
        type(neighbor_data), intent(out) :: neighbors
        real(real64), intent(in), optional :: minimum_distance
        logical, intent(in), optional :: preserve_legacy_order
        real(real64), allocatable, target :: fractional(:, :)
        integer, allocatable, target :: species_copy(:)
        real(real64), allocatable :: coordinates(:, :), distances(:)
        integer, allocatable :: atom_indices(:), species(:)
        integer, allocatable :: temporary_atom_indices(:), temporary_image_shifts(:, :)
        real(real64), allocatable :: temporary_positions(:, :), temporary_displacements(:, :)
        real(real64) :: inverse_lattice(3, 3), rmin
        integer :: atom, nmax, n, total, entry, local_entry, maximum_entries, initial_capacity
        logical :: use_legacy_order

        ! The legacy linked-cell path can omit periodic images in skew cells.
        ! Keep its fast orthogonal-cell path and use the exact image search for
        ! non-orthogonal lattices until a triclinic linked-cell implementation
        ! with equivalent coverage is available.
        use_legacy_order = .false.
        if (present(preserve_legacy_order)) use_legacy_order = preserve_legacy_order
        if (structure%pbc .and. .not. use_legacy_order .and. &
            .not. lattice_vectors_are_orthogonal(structure%lattice)) then
            call build_neighbor_list_bruteforce(structure, cutoff, neighbors)
            return
        end if

        if (.not. structure%pbc) then
            call build_neighbor_list_bruteforce(structure, cutoff, neighbors)
            return
        end if

        rmin = 1.0_real64
        if (present(minimum_distance)) rmin = minimum_distance
        inverse_lattice = inverse3(structure%lattice)
        allocate(fractional(3, structure%natoms), species_copy(structure%natoms))
        fractional = matmul(inverse_lattice, structure%positions)
        species_copy = structure%species
        call lcl_init(rmin, cutoff, structure%lattice, structure%natoms, &
                      species_copy, fractional, structure%pbc)
        nmax = lcl_nmax_nbdist(rmin, cutoff)
        allocate(coordinates(3, nmax), distances(nmax), atom_indices(nmax), species(nmax))
        if (nmax > 0 .and. structure%natoms > huge(maximum_entries)/nmax) &
            error stop "neighbor-list capacity exceeds integer range"
        maximum_entries = structure%natoms*nmax
        if (structure%natoms > maximum_entries/64) then
            initial_capacity = maximum_entries
        else
            initial_capacity = min(maximum_entries, max(nmax, 64*structure%natoms))
        end if
        allocate(temporary_atom_indices(initial_capacity), temporary_image_shifts(3, initial_capacity), &
                 temporary_positions(3, initial_capacity), temporary_displacements(3, initial_capacity))

        allocate(neighbors%offsets(structure%natoms + 1))
        neighbors%offsets(1) = 1
        entry = 0
        do atom = 1, structure%natoms
            n = nmax
            call lcl_nbdist_cart(atom, n, coordinates, distances, r_cut=cutoff, &
                                 nblist=atom_indices, nbtype=species)
            if (entry + n > size(temporary_atom_indices)) &
                call grow_neighbor_buffers(temporary_atom_indices, temporary_image_shifts, &
                    temporary_positions, temporary_displacements, entry, entry + n, maximum_entries)
            do local_entry = 1, n
                entry = entry + 1
                temporary_atom_indices(entry) = atom_indices(local_entry)
                temporary_positions(:, entry) = coordinates(:, local_entry)
                temporary_displacements(:, entry) = coordinates(:, local_entry) - structure%positions(:, atom)
                temporary_image_shifts(:, entry) = nint(matmul(inverse_lattice, &
                    coordinates(:, local_entry) - structure%positions(:, atom_indices(local_entry))))
            end do
            neighbors%offsets(atom + 1) = entry + 1
        end do
        call lcl_final()
        total = entry
        allocate(neighbors%atom_indices(total), neighbors%image_shifts(3, total), &
                 neighbors%positions(3, total), neighbors%displacements(3, total))
        neighbors%atom_indices = temporary_atom_indices(1:total)
        neighbors%image_shifts = temporary_image_shifts(:, 1:total)
        neighbors%positions = temporary_positions(:, 1:total)
        neighbors%displacements = temporary_displacements(:, 1:total)
    end subroutine build_neighbor_list

    pure logical function lattice_vectors_are_orthogonal(lattice) result(orthogonal)
        real(real64), intent(in) :: lattice(3, 3)
        real(real64) :: scale12, scale13, scale23, tolerance
        tolerance = 1.0e-12_real64
        scale12 = sqrt(sum(lattice(:, 1)**2)*sum(lattice(:, 2)**2))
        scale13 = sqrt(sum(lattice(:, 1)**2)*sum(lattice(:, 3)**2))
        scale23 = sqrt(sum(lattice(:, 2)**2)*sum(lattice(:, 3)**2))
        orthogonal = abs(dot_product(lattice(:, 1), lattice(:, 2))) <= tolerance*scale12 .and. &
                     abs(dot_product(lattice(:, 1), lattice(:, 3))) <= tolerance*scale13 .and. &
                     abs(dot_product(lattice(:, 2), lattice(:, 3))) <= tolerance*scale23
    end function lattice_vectors_are_orthogonal

    subroutine grow_neighbor_buffers(atom_indices, image_shifts, positions, displacements, used, required, limit)
        integer, allocatable, intent(inout) :: atom_indices(:), image_shifts(:, :)
        real(real64), allocatable, intent(inout) :: positions(:, :), displacements(:, :)
        integer, intent(in) :: used, required, limit
        integer, allocatable :: new_atom_indices(:), new_image_shifts(:, :)
        real(real64), allocatable :: new_positions(:, :), new_displacements(:, :)
        integer :: new_capacity
        if (size(atom_indices) > limit/2) then
            new_capacity = limit
        else
            new_capacity = min(limit, max(required, 2*max(1, size(atom_indices))))
        end if
        if (new_capacity < required) error stop "neighbor-list buffer cannot grow"
        allocate(new_atom_indices(new_capacity), new_image_shifts(3, new_capacity), &
                 new_positions(3, new_capacity), new_displacements(3, new_capacity))
        if (used > 0) then
            new_atom_indices(1:used) = atom_indices(1:used)
            new_image_shifts(:, 1:used) = image_shifts(:, 1:used)
            new_positions(:, 1:used) = positions(:, 1:used)
            new_displacements(:, 1:used) = displacements(:, 1:used)
        end if
        call move_alloc(new_atom_indices, atom_indices)
        call move_alloc(new_image_shifts, image_shifts)
        call move_alloc(new_positions, positions)
        call move_alloc(new_displacements, displacements)
    end subroutine grow_neighbor_buffers

    subroutine build_neighbor_list_bruteforce(structure, cutoff, neighbors)
        type(atomic_structure), intent(in) :: structure
        real(real64), intent(in) :: cutoff
        type(neighbor_data), intent(out) :: neighbors
        integer :: i, j, nx, ny, nz, n1, n2, n3, total, entry
        integer :: mins(3), maxs(3)
        real(real64) :: inverse_lattice(3, 3), neighbor_position(3), displacement(3), cutoff2

        cutoff2 = (cutoff + NEIGHBOR_SKIN)**2
        if (structure%pbc) then
            inverse_lattice = inverse3(structure%lattice)
            n1 = ceiling(cutoff*sqrt(sum(inverse_lattice(1, :)**2))) + 1
            n2 = ceiling(cutoff*sqrt(sum(inverse_lattice(2, :)**2))) + 1
            n3 = ceiling(cutoff*sqrt(sum(inverse_lattice(3, :)**2))) + 1
            mins = [-n1, -n2, -n3]
            maxs = [ n1,  n2,  n3]
        else
            mins = 0
            maxs = 0
        end if

        allocate(neighbors%offsets(structure%natoms + 1))
        neighbors%offsets(1) = 1
        total = 0
        do i = 1, structure%natoms
            do j = 1, structure%natoms
                do nx = mins(1), maxs(1)
                    do ny = mins(2), maxs(2)
                        do nz = mins(3), maxs(3)
                            if (i == j .and. nx == 0 .and. ny == 0 .and. nz == 0) cycle
                            neighbor_position = structure%positions(:, j)
                            if (structure%pbc) neighbor_position = neighbor_position + &
                                matmul(structure%lattice, real([nx, ny, nz], real64))
                            displacement = neighbor_position - structure%positions(:, i)
                            if (sum(displacement*displacement) <= cutoff2) total = total + 1
                        end do
                    end do
                end do
            end do
            neighbors%offsets(i + 1) = total + 1
        end do

        allocate(neighbors%atom_indices(total))
        allocate(neighbors%image_shifts(3, total))
        allocate(neighbors%positions(3, total))
        allocate(neighbors%displacements(3, total))
        entry = 0
        do i = 1, structure%natoms
            do j = 1, structure%natoms
                do nx = mins(1), maxs(1)
                    do ny = mins(2), maxs(2)
                        do nz = mins(3), maxs(3)
                            if (i == j .and. nx == 0 .and. ny == 0 .and. nz == 0) cycle
                            neighbor_position = structure%positions(:, j)
                            if (structure%pbc) neighbor_position = neighbor_position + &
                                matmul(structure%lattice, real([nx, ny, nz], real64))
                            displacement = neighbor_position - structure%positions(:, i)
                            if (sum(displacement*displacement) <= cutoff2) then
                                entry = entry + 1
                                neighbors%atom_indices(entry) = j
                                neighbors%image_shifts(:, entry) = [nx, ny, nz]
                                neighbors%positions(:, entry) = neighbor_position
                                neighbors%displacements(:, entry) = displacement
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end subroutine build_neighbor_list_bruteforce

    pure function inverse3(matrix) result(inverse)
        real(real64), intent(in) :: matrix(3, 3)
        real(real64) :: inverse(3, 3), determinant

        determinant = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) &
                    - matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) &
                    + matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
        if (abs(determinant) <= tiny(1.0_real64)) error stop "singular lattice"
        inverse(1,1) =  (matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))/determinant
        inverse(1,2) = -(matrix(1,2)*matrix(3,3)-matrix(1,3)*matrix(3,2))/determinant
        inverse(1,3) =  (matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2))/determinant
        inverse(2,1) = -(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))/determinant
        inverse(2,2) =  (matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1))/determinant
        inverse(2,3) = -(matrix(1,1)*matrix(2,3)-matrix(1,3)*matrix(2,1))/determinant
        inverse(3,1) =  (matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))/determinant
        inverse(3,2) = -(matrix(1,1)*matrix(3,2)-matrix(1,2)*matrix(3,1))/determinant
        inverse(3,3) =  (matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1))/determinant
    end function inverse3

    subroutine evaluate_structure(config, structure, neighbors, values)
        type(descriptor_config), intent(in) :: config
        type(atomic_structure), intent(in) :: structure
        type(neighbor_data), intent(in) :: neighbors
        real(real64), intent(out) :: values(:, :)
        integer :: i, first, last, n, k
        integer, allocatable :: neighbor_species(:)

        if (size(values, 1) < config%num_descriptors()) error stop "values first dimension too small"
        if (size(values, 2) < structure%natoms) error stop "values second dimension too small"
        do i = 1, structure%natoms
            first = neighbors%offsets(i)
            last = neighbors%offsets(i + 1) - 1
            n = last - first + 1
            allocate(neighbor_species(n))
            do k = 1, n
                neighbor_species(k) = structure%species(neighbors%atom_indices(first + k - 1))
            end do
            call evaluate_atom(config, neighbors%displacements(:, first:last), &
                               neighbor_species, values(:, i))
            deallocate(neighbor_species)
        end do
    end subroutine evaluate_structure

    pure subroutine monomial_powers(value, order, powers)
        real(real64), intent(in) :: value
        integer, intent(in) :: order
        real(real64), intent(out) :: powers(0:order)
        integer :: q
        powers(0) = 1.0_real64
        do q = 1, order
            powers(q) = powers(q - 1)*value
        end do
    end subroutine monomial_powers

    subroutine evaluate_angular_moment_values(config, displacements, neighbor_species, distances, &
                                              angular_cutoffs, unweighted_values, weighted_values)
        type(descriptor_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :), distances(:), angular_cutoffs(:)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(inout) :: unweighted_values(:)
        real(real64), intent(inout), optional :: weighted_values(:)
        integer :: j, n, q, entry, a, b, c
        real(real64) :: rj, ux, uy, uz, fc, species_weight, weighted_cutoff, monomial
        real(real64) :: self0, self1, sum0, sum1, coefficient
        real(real64) :: moments0(config%number_of_angular_moments)
        real(real64) :: moments1(config%number_of_angular_moments)
        real(real64) :: pair0(config%angular_order + 1), pair1(config%angular_order + 1)
        real(real64) :: px(0:config%angular_order), py(0:config%angular_order), pz(0:config%angular_order)

        moments0 = 0.0_real64
        moments1 = 0.0_real64
        self0 = 0.0_real64
        self1 = 0.0_real64
        do j = 1, size(neighbor_species)
            rj = distances(j)
            if (rj > config%angular_rc .or. rj < EPS_DISTANCE) cycle
            ux = displacements(1, j)/rj
            uy = displacements(2, j)/rj
            uz = displacements(3, j)/rj
            fc = angular_cutoffs(j)
            species_weight = config%species_weights(neighbor_species(j))
            weighted_cutoff = species_weight*fc
            self0 = self0 + fc*fc
            self1 = self1 + weighted_cutoff*weighted_cutoff
            call monomial_powers(ux, config%angular_order, px)
            call monomial_powers(uy, config%angular_order, py)
            call monomial_powers(uz, config%angular_order, pz)
            do entry = 1, config%number_of_angular_moments
                a = config%moment_x_power(entry)
                b = config%moment_y_power(entry)
                c = config%moment_z_power(entry)
                monomial = px(a)*py(b)*pz(c)
                moments0(entry) = moments0(entry) + fc*monomial
                moments1(entry) = moments1(entry) + weighted_cutoff*monomial
            end do
        end do
        do q = 0, config%angular_order
            sum0 = 0.0_real64
            sum1 = 0.0_real64
            do entry = config%moment_first(q + 1), config%moment_last(q + 1)
                coefficient = config%moment_multinomial(entry)
                sum0 = sum0 + coefficient*moments0(entry)*moments0(entry)
                sum1 = sum1 + coefficient*moments1(entry)*moments1(entry)
            end do
            pair0(q + 1) = 0.5_real64*(sum0 - self0)
            pair1(q + 1) = 0.5_real64*(sum1 - self1)
        end do
        do n = 0, config%angular_order
            unweighted_values(n + 1) = unweighted_values(n + 1) + &
                dot_product(config%angular_power_coefficients(n + 1, 1:n + 1), pair0(1:n + 1))
            if (present(weighted_values)) weighted_values(n + 1) = weighted_values(n + 1) + &
                dot_product(config%angular_power_coefficients(n + 1, 1:n + 1), pair1(1:n + 1))
        end do
    end subroutine evaluate_angular_moment_values

    subroutine evaluate_atom(config, displacements, neighbor_species, values)
        type(descriptor_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: values(:)
        integer :: nr, na, j, k, angular_neighbors
        integer :: radial1, angular1, radial2, angular2
        real(real64) :: rj, rk, cosine, fcj, fck, sj, sk, weight
        real(real64) :: angular_r0, angular_r1
        real(real64) :: tr(config%radial_order + 1)
        real(real64) :: ta(config%angular_order + 1)
        real(real64) :: temp_radial(config%radial_order + 1)
        real(real64) :: temp_angular(config%angular_order + 1)
        real(real64) :: distances(size(neighbor_species))
        real(real64) :: radial_cutoffs(size(neighbor_species)), angular_cutoffs(size(neighbor_species))

        nr = config%radial_order + 1
        na = config%angular_order + 1
        radial1 = 1
        angular1 = radial1 + nr
        radial2 = angular1 + na
        angular2 = radial2 + nr
        angular_r0 = -1.0_real64
        angular_r1 = 1.0_real64
        if (config%version == 1) then
            angular_r0 = 0.0_real64
            angular_r1 = PI_ACCELNET
        end if
        values(1:config%num_descriptors()) = 0.0_real64
        call add_version10_center(config, neighbor_species, values)

        do j = 1, size(neighbor_species)
            distances(j) = sqrt(dot_product(displacements(:, j), displacements(:, j)))
            radial_cutoffs(j) = cutoff_value(distances(j), config%radial_rc, &
                config%cutoff_type, config%cutoff_alpha)
            angular_cutoffs(j) = cutoff_value(distances(j), config%angular_rc, &
                config%cutoff_type, config%cutoff_alpha)
        end do
        angular_neighbors = count(distances <= config%angular_rc .and. distances > EPS_DISTANCE)

        do j = 1, size(neighbor_species)
            rj = distances(j)
            if (rj <= config%radial_rc .and. rj > EPS_DISTANCE) then
                fcj = radial_cutoffs(j)
                call chebyshev_values(rj, 0.0_real64, config%radial_rc, config%radial_order, tr)
                temp_radial = fcj*tr
                values(radial1:radial1 + nr - 1) = values(radial1:radial1 + nr - 1) + temp_radial
                if (config%num_species > 1) then
                    sj = config%species_weights(neighbor_species(j))
                    values(radial2:radial2 + nr - 1) = values(radial2:radial2 + nr - 1) + sj*temp_radial
                end if
            end if
            if (chebyshev_uses_moments(config, angular_neighbors)) cycle
            if (rj > config%angular_rc .or. rj < EPS_DISTANCE) cycle
            if (config%num_species > 1) sj = config%species_weights(neighbor_species(j))
            do k = j + 1, size(neighbor_species)
                rk = distances(k)
                if (rk > config%angular_rc .or. rk < EPS_DISTANCE) cycle
                cosine = dot_product(displacements(:, j), displacements(:, k))/(rj*rk)
                fcj = angular_cutoffs(j)
                fck = angular_cutoffs(k)
                weight = fcj*fck
                call chebyshev_values(cosine, angular_r0, angular_r1, config%angular_order, ta)
                temp_angular = weight*ta
                values(angular1:angular1 + na - 1) = values(angular1:angular1 + na - 1) + temp_angular
                if (config%num_species > 1) then
                    sk = config%species_weights(neighbor_species(k))
                    values(angular2:angular2 + na - 1) = values(angular2:angular2 + na - 1) + &
                        sj*sk*temp_angular
                end if
            end do
        end do
        if (chebyshev_uses_moments(config, angular_neighbors)) then
            if (config%num_species > 1) then
                call evaluate_angular_moment_values(config, displacements, neighbor_species, distances, &
                    angular_cutoffs, values(angular1:angular1 + na - 1), &
                    values(angular2:angular2 + na - 1))
            else
                call evaluate_angular_moment_values(config, displacements, neighbor_species, distances, &
                    angular_cutoffs, values(angular1:angular1 + na - 1))
            end if
        end if
    end subroutine evaluate_atom

    subroutine evaluate_atom_with_derivatives(config, displacements, neighbor_species, values, &
                                              derivative_center, derivative_neighbors)
        type(descriptor_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: values(:)
        real(real64), intent(out) :: derivative_center(:, :)
        real(real64), intent(out) :: derivative_neighbors(:, :, :)
        integer :: nr, na, j, k, component, basis
        integer :: radial1, angular1, radial2, angular2, index
        real(real64) :: rj, rk, cosine, fcj, fck, dfcj, dfck, sj, sk, weight
        real(real64) :: angular_r0, angular_r1
        real(real64) :: inv_rj2, inv_rk2, inv_rjrk
        real(real64) :: tr(config%radial_order + 1), dtr(config%radial_order + 1)
        real(real64) :: ta(config%angular_order + 1), dta(config%angular_order + 1)
        real(real64) :: dcj(3), dck(3), dwj(3), dwk(3), termj, termk
        real(real64) :: distances(size(neighbor_species)), unit_vectors(3, size(neighbor_species))
        real(real64) :: radial_cutoffs(size(neighbor_species)), angular_cutoffs(size(neighbor_species))
        real(real64) :: radial_cutoff_derivatives(size(neighbor_species))
        real(real64) :: angular_cutoff_derivatives(size(neighbor_species))

        nr = config%radial_order + 1
        na = config%angular_order + 1
        radial1 = 1
        angular1 = radial1 + nr
        radial2 = angular1 + na
        angular2 = radial2 + nr
        angular_r0 = -1.0_real64
        angular_r1 = 1.0_real64
        if (config%version == 1) then
            angular_r0 = 0.0_real64
            angular_r1 = PI_ACCELNET
        end if
        values(1:config%num_descriptors()) = 0.0_real64
        derivative_center(:, 1:config%num_descriptors()) = 0.0_real64
        derivative_neighbors(:, 1:config%num_descriptors(), :) = 0.0_real64
        call add_version10_center(config, neighbor_species, values)

        unit_vectors = 0.0_real64
        do j = 1, size(neighbor_species)
            distances(j) = sqrt(dot_product(displacements(:, j), displacements(:, j)))
            if (distances(j) > EPS_DISTANCE) unit_vectors(:, j) = displacements(:, j)/distances(j)
            radial_cutoffs(j) = cutoff_value(distances(j), config%radial_rc, &
                config%cutoff_type, config%cutoff_alpha)
            angular_cutoffs(j) = cutoff_value(distances(j), config%angular_rc, &
                config%cutoff_type, config%cutoff_alpha)
            radial_cutoff_derivatives(j) = cutoff_derivative(distances(j), config%radial_rc, &
                config%cutoff_type, config%cutoff_alpha)
            angular_cutoff_derivatives(j) = cutoff_derivative(distances(j), config%angular_rc, &
                config%cutoff_type, config%cutoff_alpha)
        end do

        do j = 1, size(neighbor_species)
            rj = distances(j)
            if (rj <= config%radial_rc .and. rj > EPS_DISTANCE) then
                fcj = radial_cutoffs(j)
                dfcj = radial_cutoff_derivatives(j)
                call chebyshev_values_derivatives(rj, 0.0_real64, config%radial_rc, &
                                                  config%radial_order, tr, dtr)
                sj = config%species_weights(neighbor_species(j))
                do basis = 1, nr
                    values(radial1 + basis - 1) = values(radial1 + basis - 1) + fcj*tr(basis)
                    do component = 1, 3
                        termj = unit_vectors(component, j)*(dfcj*tr(basis) + fcj*dtr(basis))
                        derivative_neighbors(component, radial1 + basis - 1, j) = &
                            derivative_neighbors(component, radial1 + basis - 1, j) + termj
                        derivative_center(component, radial1 + basis - 1) = &
                            derivative_center(component, radial1 + basis - 1) - termj
                    end do
                    if (config%num_species > 1) then
                        index = radial2 + basis - 1
                        values(index) = values(index) + sj*fcj*tr(basis)
                        derivative_neighbors(:, index, j) = derivative_neighbors(:, index, j) + &
                            sj*unit_vectors(:, j)*(dfcj*tr(basis) + fcj*dtr(basis))
                        derivative_center(:, index) = derivative_center(:, index) - &
                            sj*unit_vectors(:, j)*(dfcj*tr(basis) + fcj*dtr(basis))
                    end if
                end do
            end if

            if (rj > config%angular_rc) cycle
            sj = config%species_weights(neighbor_species(j))
            do k = j + 1, size(neighbor_species)
                rk = distances(k)
                if (rk > config%angular_rc .or. rk < EPS_DISTANCE) cycle
                cosine = dot_product(unit_vectors(:, j), unit_vectors(:, k))
                fcj = angular_cutoffs(j)
                fck = angular_cutoffs(k)
                dfcj = angular_cutoff_derivatives(j)
                dfck = angular_cutoff_derivatives(k)
                weight = fcj*fck
                sk = config%species_weights(neighbor_species(k))
                call chebyshev_values_derivatives(cosine, angular_r0, angular_r1, &
                                                  config%angular_order, ta, dta)
                inv_rj2 = 1.0_real64/(rj*rj)
                inv_rk2 = 1.0_real64/(rk*rk)
                inv_rjrk = 1.0_real64/(rj*rk)
                dcj = -cosine*displacements(:, j)*inv_rj2 + displacements(:, k)*inv_rjrk
                dck = -cosine*displacements(:, k)*inv_rk2 + displacements(:, j)*inv_rjrk
                dwj = dfcj*fck*unit_vectors(:, j)
                dwk = fcj*dfck*unit_vectors(:, k)
                do basis = 1, na
                    index = angular1 + basis - 1
                    values(index) = values(index) + weight*ta(basis)
                    do component = 1, 3
                        termj = dwj(component)*ta(basis) + weight*dta(basis)*dcj(component)
                        termk = dwk(component)*ta(basis) + weight*dta(basis)*dck(component)
                        derivative_neighbors(component, index, j) = &
                            derivative_neighbors(component, index, j) + termj
                        derivative_neighbors(component, index, k) = &
                            derivative_neighbors(component, index, k) + termk
                        derivative_center(component, index) = &
                            derivative_center(component, index) - termj - termk
                    end do
                    if (config%num_species > 1) then
                        index = angular2 + basis - 1
                        values(index) = values(index) + sj*sk*weight*ta(basis)
                        derivative_neighbors(:, index, j) = derivative_neighbors(:, index, j) + &
                            sj*sk*(dwj*ta(basis) + weight*dta(basis)*dcj)
                        derivative_neighbors(:, index, k) = derivative_neighbors(:, index, k) + &
                            sj*sk*(dwk*ta(basis) + weight*dta(basis)*dck)
                        derivative_center(:, index) = derivative_center(:, index) - &
                            sj*sk*(dwj*ta(basis) + weight*dta(basis)*dcj) - &
                            sj*sk*(dwk*ta(basis) + weight*dta(basis)*dck)
                    end if
                end do
            end do
        end do
    end subroutine evaluate_atom_with_derivatives

    subroutine contract_angular_moment_derivatives(config, displacements, neighbor_species, distances, &
                                                   angular_cutoffs, angular_cutoff_derivatives, &
                                                   coefficients0, contracted_center, contracted_neighbors, &
                                                   coefficients1)
        type(descriptor_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :), distances(:), angular_cutoffs(:)
        real(real64), intent(in) :: angular_cutoff_derivatives(:), coefficients0(:)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(inout) :: contracted_center(3), contracted_neighbors(:, :)
        real(real64), intent(in), optional :: coefficients1(:)
        integer :: j, n, q, entry, a, b, c
        real(real64) :: rj, fc, dfc, species_weight, weighted_cutoff, monomial
        real(real64) :: ux, uy, uz, gx, gy, gz, udot, coefficient, scale, correction
        real(real64) :: dmonomial_x, dmonomial_y, dmonomial_z
        real(real64) :: derivative_x, derivative_y, derivative_z, dm0x, dm0y, dm0z
        real(real64) :: moments0(config%number_of_angular_moments)
        real(real64) :: moments1(config%number_of_angular_moments)
        real(real64) :: power_coefficients0(config%angular_order + 1)
        real(real64) :: power_coefficients1(config%angular_order + 1)
        real(real64) :: px(0:config%angular_order), py(0:config%angular_order), pz(0:config%angular_order)

        moments0 = 0.0_real64
        moments1 = 0.0_real64
        do j = 1, size(neighbor_species)
            rj = distances(j)
            if (rj > config%angular_rc .or. rj < EPS_DISTANCE) cycle
            ux = displacements(1, j)/rj
            uy = displacements(2, j)/rj
            uz = displacements(3, j)/rj
            fc = angular_cutoffs(j)
            species_weight = config%species_weights(neighbor_species(j))
            call monomial_powers(ux, config%angular_order, px)
            call monomial_powers(uy, config%angular_order, py)
            call monomial_powers(uz, config%angular_order, pz)
            do entry = 1, config%number_of_angular_moments
                monomial = px(config%moment_x_power(entry))*py(config%moment_y_power(entry))* &
                           pz(config%moment_z_power(entry))
                moments0(entry) = moments0(entry) + fc*monomial
                moments1(entry) = moments1(entry) + species_weight*fc*monomial
            end do
        end do

        power_coefficients0 = 0.0_real64
        power_coefficients1 = 0.0_real64
        do q = 0, config%angular_order
            do n = q, config%angular_order
                power_coefficients0(q + 1) = power_coefficients0(q + 1) + &
                    coefficients0(n + 1)*config%angular_power_coefficients(n + 1, q + 1)
                if (present(coefficients1)) power_coefficients1(q + 1) = power_coefficients1(q + 1) + &
                    coefficients1(n + 1)*config%angular_power_coefficients(n + 1, q + 1)
            end do
        end do

        do j = 1, size(neighbor_species)
            rj = distances(j)
            if (rj > config%angular_rc .or. rj < EPS_DISTANCE) cycle
            ux = displacements(1, j)/rj
            uy = displacements(2, j)/rj
            uz = displacements(3, j)/rj
            fc = angular_cutoffs(j)
            dfc = angular_cutoff_derivatives(j)
            species_weight = config%species_weights(neighbor_species(j))
            weighted_cutoff = species_weight*fc
            call monomial_powers(ux, config%angular_order, px)
            call monomial_powers(uy, config%angular_order, py)
            call monomial_powers(uz, config%angular_order, pz)
            derivative_x = 0.0_real64
            derivative_y = 0.0_real64
            derivative_z = 0.0_real64
            do q = 0, config%angular_order
                do entry = config%moment_first(q + 1), config%moment_last(q + 1)
                    a = config%moment_x_power(entry)
                    b = config%moment_y_power(entry)
                    c = config%moment_z_power(entry)
                    monomial = px(a)*py(b)*pz(c)
                    gx = 0.0_real64
                    gy = 0.0_real64
                    gz = 0.0_real64
                    if (a > 0) gx = real(a, real64)*px(a - 1)*py(b)*pz(c)
                    if (b > 0) gy = real(b, real64)*px(a)*py(b - 1)*pz(c)
                    if (c > 0) gz = real(c, real64)*px(a)*py(b)*pz(c - 1)
                    udot = ux*gx + uy*gy + uz*gz
                    dmonomial_x = (gx - ux*udot)/rj
                    dmonomial_y = (gy - uy*udot)/rj
                    dmonomial_z = (gz - uz*udot)/rj
                    dm0x = dfc*ux*monomial + fc*dmonomial_x
                    dm0y = dfc*uy*monomial + fc*dmonomial_y
                    dm0z = dfc*uz*monomial + fc*dmonomial_z
                    coefficient = config%moment_multinomial(entry)
                    scale = power_coefficients0(q + 1)*coefficient*moments0(entry)
                    if (present(coefficients1)) scale = scale + power_coefficients1(q + 1)* &
                        coefficient*moments1(entry)*species_weight
                    derivative_x = derivative_x + scale*dm0x
                    derivative_y = derivative_y + scale*dm0y
                    derivative_z = derivative_z + scale*dm0z
                end do
                correction = power_coefficients0(q + 1)*fc*dfc
                if (present(coefficients1)) correction = correction + power_coefficients1(q + 1)* &
                    weighted_cutoff*species_weight*dfc
                derivative_x = derivative_x - correction*ux
                derivative_y = derivative_y - correction*uy
                derivative_z = derivative_z - correction*uz
            end do
            contracted_neighbors(1, j) = contracted_neighbors(1, j) + derivative_x
            contracted_neighbors(2, j) = contracted_neighbors(2, j) + derivative_y
            contracted_neighbors(3, j) = contracted_neighbors(3, j) + derivative_z
            contracted_center(1) = contracted_center(1) - derivative_x
            contracted_center(2) = contracted_center(2) - derivative_y
            contracted_center(3) = contracted_center(3) - derivative_z
        end do
    end subroutine contract_angular_moment_derivatives

    subroutine contract_atom_derivatives(config, displacements, neighbor_species, coefficients, &
                                         contracted_center, contracted_neighbors)
        type(descriptor_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :), coefficients(:)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: contracted_center(3), contracted_neighbors(:, :)
        integer :: nr, na, j, k, basis, radial1, angular1, radial2, angular2, angular_neighbors
        real(real64) :: rj, rk, cosine, fcj, fck, dfcj, dfck, sj, sk, weight, coefficient
        real(real64) :: angular_r0, angular_r1, radial_contraction
        real(real64) :: value_contraction, derivative_contraction
        real(real64) :: inv_rj2, inv_rk2, inv_rjrk
        real(real64) :: tr(config%radial_order + 1), dtr(config%radial_order + 1)
        real(real64) :: ta(config%angular_order + 1), dta(config%angular_order + 1)
        real(real64) :: dcjx, dcjy, dcjz, dckx, dcky, dckz
        real(real64) :: dwjx, dwjy, dwjz, dwkx, dwky, dwkz
        real(real64) :: derivative_jx, derivative_jy, derivative_jz
        real(real64) :: derivative_kx, derivative_ky, derivative_kz
        real(real64) :: distances(size(neighbor_species)), unit_vectors(3, size(neighbor_species))
        real(real64) :: radial_cutoffs(size(neighbor_species)), angular_cutoffs(size(neighbor_species))
        real(real64) :: radial_cutoff_derivatives(size(neighbor_species))
        real(real64) :: angular_cutoff_derivatives(size(neighbor_species))

        if (size(coefficients) < config%num_descriptors()) error stop "Chebyshev coefficient array too small"
        if (size(contracted_neighbors, 1) /= 3 .or. &
            size(contracted_neighbors, 2) < size(neighbor_species)) error stop "Chebyshev contracted array too small"
        nr = config%radial_order + 1
        na = config%angular_order + 1
        radial1 = 1
        angular1 = radial1 + nr
        radial2 = angular1 + na
        angular2 = radial2 + nr
        angular_r0 = -1.0_real64
        angular_r1 = 1.0_real64
        if (config%version == 1) then
            angular_r0 = 0.0_real64
            angular_r1 = PI_ACCELNET
        end if
        contracted_center = 0.0_real64
        contracted_neighbors(:, 1:size(neighbor_species)) = 0.0_real64

        unit_vectors = 0.0_real64
        do j = 1, size(neighbor_species)
            distances(j) = sqrt(dot_product(displacements(:, j), displacements(:, j)))
            if (distances(j) > EPS_DISTANCE) unit_vectors(:, j) = displacements(:, j)/distances(j)
            radial_cutoffs(j) = cutoff_value(distances(j), config%radial_rc, &
                config%cutoff_type, config%cutoff_alpha)
            angular_cutoffs(j) = cutoff_value(distances(j), config%angular_rc, &
                config%cutoff_type, config%cutoff_alpha)
            radial_cutoff_derivatives(j) = cutoff_derivative(distances(j), config%radial_rc, &
                config%cutoff_type, config%cutoff_alpha)
            angular_cutoff_derivatives(j) = cutoff_derivative(distances(j), config%angular_rc, &
                config%cutoff_type, config%cutoff_alpha)
        end do
        angular_neighbors = count(distances <= config%angular_rc .and. distances > EPS_DISTANCE)

        do j = 1, size(neighbor_species)
            rj = distances(j)
            if (rj <= config%radial_rc .and. rj > EPS_DISTANCE) then
                fcj = radial_cutoffs(j)
                dfcj = radial_cutoff_derivatives(j)
                call chebyshev_values_derivatives(rj, 0.0_real64, config%radial_rc, &
                                                  config%radial_order, tr, dtr)
                sj = config%species_weights(neighbor_species(j))
                radial_contraction = 0.0_real64
                do basis = 1, nr
                    coefficient = coefficients(radial1 + basis - 1)
                    if (config%num_species > 1) coefficient = coefficient + &
                        sj*coefficients(radial2 + basis - 1)
                    radial_contraction = radial_contraction + coefficient*(dfcj*tr(basis) + fcj*dtr(basis))
                end do
                derivative_jx = radial_contraction*unit_vectors(1, j)
                derivative_jy = radial_contraction*unit_vectors(2, j)
                derivative_jz = radial_contraction*unit_vectors(3, j)
                contracted_neighbors(1, j) = contracted_neighbors(1, j) + derivative_jx
                contracted_neighbors(2, j) = contracted_neighbors(2, j) + derivative_jy
                contracted_neighbors(3, j) = contracted_neighbors(3, j) + derivative_jz
                contracted_center(1) = contracted_center(1) - derivative_jx
                contracted_center(2) = contracted_center(2) - derivative_jy
                contracted_center(3) = contracted_center(3) - derivative_jz
            end if

            if (chebyshev_uses_moments(config, angular_neighbors)) cycle
            if (rj > config%angular_rc .or. rj < EPS_DISTANCE) cycle
            sj = config%species_weights(neighbor_species(j))
            do k = j + 1, size(neighbor_species)
                rk = distances(k)
                if (rk > config%angular_rc .or. rk < EPS_DISTANCE) cycle
                cosine = dot_product(unit_vectors(:, j), unit_vectors(:, k))
                fcj = angular_cutoffs(j)
                fck = angular_cutoffs(k)
                dfcj = angular_cutoff_derivatives(j)
                dfck = angular_cutoff_derivatives(k)
                weight = fcj*fck
                sk = config%species_weights(neighbor_species(k))
                call chebyshev_values_derivatives(cosine, angular_r0, angular_r1, &
                                                  config%angular_order, ta, dta)
                inv_rj2 = 1.0_real64/(rj*rj)
                inv_rk2 = 1.0_real64/(rk*rk)
                inv_rjrk = 1.0_real64/(rj*rk)
                dcjx = -cosine*displacements(1, j)*inv_rj2 + displacements(1, k)*inv_rjrk
                dcjy = -cosine*displacements(2, j)*inv_rj2 + displacements(2, k)*inv_rjrk
                dcjz = -cosine*displacements(3, j)*inv_rj2 + displacements(3, k)*inv_rjrk
                dckx = -cosine*displacements(1, k)*inv_rk2 + displacements(1, j)*inv_rjrk
                dcky = -cosine*displacements(2, k)*inv_rk2 + displacements(2, j)*inv_rjrk
                dckz = -cosine*displacements(3, k)*inv_rk2 + displacements(3, j)*inv_rjrk
                dwjx = dfcj*fck*unit_vectors(1, j)
                dwjy = dfcj*fck*unit_vectors(2, j)
                dwjz = dfcj*fck*unit_vectors(3, j)
                dwkx = fcj*dfck*unit_vectors(1, k)
                dwky = fcj*dfck*unit_vectors(2, k)
                dwkz = fcj*dfck*unit_vectors(3, k)
                value_contraction = 0.0_real64
                derivative_contraction = 0.0_real64
                do basis = 1, na
                    coefficient = coefficients(angular1 + basis - 1)
                    if (config%num_species > 1) coefficient = coefficient + &
                        sj*sk*coefficients(angular2 + basis - 1)
                    value_contraction = value_contraction + coefficient*ta(basis)
                    derivative_contraction = derivative_contraction + coefficient*dta(basis)
                end do
                derivative_contraction = weight*derivative_contraction
                derivative_jx = dwjx*value_contraction + dcjx*derivative_contraction
                derivative_jy = dwjy*value_contraction + dcjy*derivative_contraction
                derivative_jz = dwjz*value_contraction + dcjz*derivative_contraction
                derivative_kx = dwkx*value_contraction + dckx*derivative_contraction
                derivative_ky = dwky*value_contraction + dcky*derivative_contraction
                derivative_kz = dwkz*value_contraction + dckz*derivative_contraction
                contracted_neighbors(1, j) = contracted_neighbors(1, j) + derivative_jx
                contracted_neighbors(2, j) = contracted_neighbors(2, j) + derivative_jy
                contracted_neighbors(3, j) = contracted_neighbors(3, j) + derivative_jz
                contracted_neighbors(1, k) = contracted_neighbors(1, k) + derivative_kx
                contracted_neighbors(2, k) = contracted_neighbors(2, k) + derivative_ky
                contracted_neighbors(3, k) = contracted_neighbors(3, k) + derivative_kz
                contracted_center(1) = contracted_center(1) - derivative_jx - derivative_kx
                contracted_center(2) = contracted_center(2) - derivative_jy - derivative_ky
                contracted_center(3) = contracted_center(3) - derivative_jz - derivative_kz
            end do
        end do
        if (chebyshev_uses_moments(config, angular_neighbors)) then
            if (config%num_species > 1) then
                call contract_angular_moment_derivatives(config, displacements, neighbor_species, distances, &
                    angular_cutoffs, angular_cutoff_derivatives, coefficients(angular1:angular1 + na - 1), &
                    contracted_center, contracted_neighbors, coefficients(angular2:angular2 + na - 1))
            else
                call contract_angular_moment_derivatives(config, displacements, neighbor_species, distances, &
                    angular_cutoffs, angular_cutoff_derivatives, coefficients(angular1:angular1 + na - 1), &
                    contracted_center, contracted_neighbors)
            end if
        end if
    end subroutine contract_atom_derivatives

    subroutine add_version10_center(config, neighbor_species, values)
        type(descriptor_config), intent(in) :: config
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(inout) :: values(:)
        integer :: basis, nr, radial2, center_species
        real(real64) :: center_value

        if (config%version /= 10) return
        nr = config%radial_order + 1
        radial2 = nr + config%angular_order + 2
        center_species = 1
        if (config%num_species > 1) then
            ! AccelNet version 10 indexes the mapped neighbor-type array with
            ! the global central type ID. Preserve that observable behavior.
            if (config%central_type_index > size(neighbor_species)) &
                error stop "too few neighbors for AccelNet version=10 central type lookup"
            center_species = neighbor_species(config%central_type_index)
            if (center_species < 1 .or. center_species > config%num_species) &
                error stop "version=10 center species index is out of range"
        end if
        do basis = 1, nr
            if (mod(basis - 1, 2) == 0) then
                center_value = 1.0_real64
            else
                center_value = -1.0_real64
            end if
            values(basis) = values(basis) + center_value
            if (config%num_species > 1) values(radial2 + basis - 1) = &
                values(radial2 + basis - 1) + config%species_weights(center_species)*center_value
        end do
    end subroutine add_version10_center

    subroutine write_descriptor_file(filename, config, structures, labels)
        character(len=*), intent(in) :: filename
        type(descriptor_config), intent(in) :: config
        type(atomic_structure), intent(in) :: structures(:)
        character(len=*), intent(in), optional :: labels(:)
        type(neighbor_data) :: neighbors
        real(real64), allocatable :: values(:, :)
        integer :: unit, s, i, g
        character(len=1024) :: label

        open(newunit=unit, file=filename, status="replace", action="write")
        write(unit, "(A)") "ACCELNET_DESCRIPTOR_HEX_V1"
        write(unit, "(3(I0,1X))") size(structures), config%num_species, config%num_descriptors()
        do s = 1, size(structures)
            label = ""
            if (present(labels)) label = trim(labels(s))
            write(unit, "(A)") trim(label)
            write(unit, "(I0)") structures(s)%natoms
            call build_neighbor_list(structures(s), max(config%radial_rc, config%angular_rc), neighbors)
            allocate(values(config%num_descriptors(), structures(s)%natoms))
            call evaluate_structure(config, structures(s), neighbors, values)
            do i = 1, structures(s)%natoms
                write(unit, "(I0,1X,I0)") i, structures(s)%species(i)
                write(unit, "(*(Z16.16,1X))") (transfer(values(g, i), 0_int64), g = 1, config%num_descriptors())
            end do
            deallocate(values)
        end do
        close(unit)
    end subroutine write_descriptor_file

end module accelnet_descriptors
