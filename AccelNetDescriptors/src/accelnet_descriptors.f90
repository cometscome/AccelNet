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
    public :: read_n2p2_data
    public :: build_neighbor_list
    public :: evaluate_atom
    public :: evaluate_atom_with_derivatives
    public :: contract_atom_derivatives
    public :: evaluate_structure
    public :: write_descriptor_file
    public :: cutoff_value, cutoff_derivative, validate_cutoff_parameters
    public :: cutoff_value_derivative
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
        include 'cutoff_value.inc'
    end function cutoff_value

    pure real(real64) function cutoff_derivative(distance, rc, cutoff_type, alpha) result(value)
        include 'cutoff_derivative.inc'
    end function cutoff_derivative

    pure subroutine cutoff_value_derivative(distance, rc, cutoff_type, alpha, value, derivative)
        real(real64), intent(in) :: distance, rc
        integer, intent(in), optional :: cutoff_type
        real(real64), intent(in), optional :: alpha
        real(real64), intent(out) :: value, derivative
        integer :: kind
        real(real64) :: inner, inverse_width, x, t, width, core_derivative

        ! Share the expensive elementary functions when both outputs are needed.
        kind = CUTOFF_COS
        if (present(cutoff_type)) kind = cutoff_type
        value = 0.0_real64
        derivative = 0.0_real64
        if (distance >= rc) return
        select case(kind)
        case(CUTOFF_HARD)
            value = 1.0_real64
        case(CUTOFF_TANHU, CUTOFF_TANH)
            t = tanh(1.0_real64 - distance/rc)
            value = t*t*t
            derivative = 3.0_real64*t*t*(t*t - 1.0_real64)/rc
            if (kind == CUTOFF_TANH) then
                value = value/tanh(1.0_real64)**3
                derivative = derivative/tanh(1.0_real64)**3
            end if
        case(CUTOFF_FRACTIONAL)
            width = rc
            if (present(alpha)) width = alpha*rc
            x = (distance - rc)/width
            value = x*x/(1.0_real64 + x*x)
            derivative = 2.0_real64*x/(width*(1.0_real64 + x*x)**2)
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
                core_derivative = -0.5_real64*PI_ACCELNET*sin(PI_ACCELNET*x)
            case(CUTOFF_EXP)
                t = 1.0_real64/(x*x - 1.0_real64)
                value = exp(1.0_real64 + t)
                core_derivative = -2.0_real64*x*t*t*value
            case(CUTOFF_POLY1)
                value = (2.0_real64*x - 3.0_real64)*x*x + 1.0_real64
                core_derivative = x*(6.0_real64*x - 6.0_real64)
            case(CUTOFF_POLY2)
                value = ((15.0_real64 - 6.0_real64*x)*x - 10.0_real64)*x*x*x + 1.0_real64
                core_derivative = x*x*((60.0_real64 - 30.0_real64*x)*x - 30.0_real64)
            case(CUTOFF_POLY3)
                value = (x*(x*(20.0_real64*x - 70.0_real64) + 84.0_real64) - 35.0_real64)*x**4 + 1.0_real64
                core_derivative = x**3*(x*(x*(140.0_real64*x - 420.0_real64) + &
                    420.0_real64) - 140.0_real64)
            case(CUTOFF_POLY4)
                value = (x*(x*((315.0_real64 - 70.0_real64*x)*x - 540.0_real64) + &
                    420.0_real64) - 126.0_real64)*x**5 + 1.0_real64
                core_derivative = x**4*(x*(x*((2520.0_real64 - 630.0_real64*x)*x - &
                    3780.0_real64) + 2520.0_real64) - 630.0_real64)
            case default
                core_derivative = 0.0_real64
            end select
            derivative = inverse_width*core_derivative
        end select
    end subroutine cutoff_value_derivative

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

    subroutine read_n2p2_data(filename, species_names, structures)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: species_names(:)
        type(atomic_structure), allocatable, intent(out) :: structures(:)
        character(len=2048) :: line
        character(len=32) :: key
        character(len=16) :: atom_name
        integer, allocatable :: atom_counts(:), lattice_counts(:), atom_indices(:)
        integer :: unit, ios, number_of_structures, current, i, itype
        real(real64) :: x, y, z
        logical :: inside

        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "cannot open n2p2 input.data"
        number_of_structures = 0
        do
            read(unit, "(A)", iostat=ios) line
            if (ios /= 0) exit
            call n2p2_line_key(line, key, ios)
            if (ios == 0 .and. trim(key) == "begin") number_of_structures = number_of_structures + 1
        end do
        if (number_of_structures < 1) error stop "n2p2 input.data contains no structures"
        allocate(atom_counts(number_of_structures)); atom_counts = 0
        rewind(unit); current = 0; inside = .false.
        do
            read(unit, "(A)", iostat=ios) line
            if (ios /= 0) exit
            call n2p2_line_key(line, key, ios)
            if (ios /= 0) cycle
            select case(trim(key))
            case("begin")
                if (inside) error stop "nested begin in n2p2 input.data"
                current = current + 1; inside = .true.
            case("atom")
                if (.not. inside) error stop "atom outside structure in n2p2 input.data"
                atom_counts(current) = atom_counts(current) + 1
            case("end")
                if (.not. inside) error stop "end outside structure in n2p2 input.data"
                inside = .false.
            end select
        end do
        if (inside) error stop "unterminated structure in n2p2 input.data"
        if (any(atom_counts < 1)) error stop "empty structure in n2p2 input.data"

        allocate(structures(number_of_structures), lattice_counts(number_of_structures), &
                 atom_indices(number_of_structures))
        lattice_counts = 0; atom_indices = 0
        do i = 1, number_of_structures
            structures(i)%natoms = atom_counts(i)
            allocate(structures(i)%positions(3, atom_counts(i)), structures(i)%species(atom_counts(i)))
        end do
        rewind(unit); current = 0; inside = .false.
        do
            read(unit, "(A)", iostat=ios) line
            if (ios /= 0) exit
            call n2p2_line_key(line, key, ios)
            if (ios /= 0) cycle
            select case(trim(key))
            case("begin")
                current = current + 1; inside = .true.
            case("lattice")
                if (.not. inside) error stop "lattice outside structure in n2p2 input.data"
                lattice_counts(current) = lattice_counts(current) + 1
                if (lattice_counts(current) > 3) error stop "too many lattice vectors in n2p2 input.data"
                read(line, *, iostat=ios) key, x, y, z
                if (ios /= 0) error stop "invalid lattice line in n2p2 input.data"
                structures(current)%lattice(:, lattice_counts(current)) = [x, y, z]
            case("atom")
                atom_indices(current) = atom_indices(current) + 1
                read(line, *, iostat=ios) key, x, y, z, atom_name
                if (ios /= 0) error stop "invalid atom line in n2p2 input.data"
                itype = find_species(atom_name, species_names)
                if (itype == 0) error stop "unknown species in n2p2 input.data"
                structures(current)%positions(:, atom_indices(current)) = [x, y, z]
                structures(current)%species(atom_indices(current)) = itype
            case("end")
                if (lattice_counts(current) /= 0 .and. lattice_counts(current) /= 3) &
                    error stop "n2p2 structure must contain zero or three lattice vectors"
                structures(current)%pbc = lattice_counts(current) == 3
                inside = .false.
            end select
        end do
        close(unit)
    end subroutine read_n2p2_data

    subroutine n2p2_line_key(line, key, ios)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: key
        integer, intent(out) :: ios
        integer :: i, code
        key = ""
        read(line, *, iostat=ios) key
        if (ios /= 0) return
        do i = 1, len_trim(key)
            code = iachar(key(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) key(i:i) = achar(code + 32)
        end do
    end subroutine n2p2_line_key

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
        real(real64) :: center_translation(3), image_bound
        integer :: atom, nmax, n, total, entry, local_entry, maximum_entries, initial_capacity, neighbor_status
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
        ! A small model rmin makes the packing bound enormous even for a few
        ! hundred atoms. Bound each atom's periodic copies geometrically too:
        ! an interval of length 2*Rc*|inverse_lattice_row| contains at most
        ! ceil(length)+1 integers, independently of its fractional origin.
        image_bound = real(structure%natoms,real64)* &
            product(ceiling(2*(cutoff+NEIGHBOR_SKIN)*sqrt(sum(inverse_lattice**2,dim=2)))+1.0_real64)
        nmax = min(nmax,int(min(image_bound,real(huge(nmax),real64))))
        allocate(coordinates(3, nmax), distances(nmax), atom_indices(nmax), species(nmax))
        ! CSR offsets need room for the terminal entry; only actual growth,
        ! not the loose N*nmax upper bound, should exhaust the index range.
        maximum_entries = int(min(real(structure%natoms,real64)*real(nmax,real64), &
            real(huge(maximum_entries)-1,real64)))
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
            ! lcl_init wraps its fractional coordinates in place. Restore
            ! image positions relative to the caller's (possibly unwrapped)
            ! central atom before forming force and virial displacements.
            center_translation = matmul(structure%lattice, &
                real(nint(matmul(inverse_lattice, structure%positions(:,atom)) - fractional(:,atom)), real64))
            n = nmax
            call lcl_nbdist_cart(atom, n, coordinates, distances, r_cut=cutoff, &
                                 nblist=atom_indices, nbtype=species, stat=neighbor_status)
            if (neighbor_status /= 0) error stop "neighbor-list geometric capacity exceeded"
            coordinates(:,1:n) = coordinates(:,1:n) + spread(center_translation, 2, n)
            if (entry > maximum_entries-n) error stop "neighbor-list entries exceed integer range"
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
        integer :: mins(3), maxs(3), base_shift(3)
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
                ! The bounded image search must be centered on the central
                ! atom even when atoms carry independent unwrapped coordinates.
                ! Truncation leaves primary-cell pair ordering unchanged and
                ! brings each fractional separation into (-1,1).
                base_shift = 0
                if (structure%pbc) base_shift = -int(matmul(inverse_lattice, &
                    structure%positions(:,j) - structure%positions(:,i)))
                do nx = mins(1), maxs(1)
                    do ny = mins(2), maxs(2)
                        do nz = mins(3), maxs(3)
                            if (i == j .and. nx == 0 .and. ny == 0 .and. nz == 0) cycle
                            neighbor_position = structure%positions(:, j)
                            if (structure%pbc) neighbor_position = neighbor_position + &
                                matmul(structure%lattice, real(base_shift + [nx, ny, nz], real64))
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
                base_shift = 0
                if (structure%pbc) base_shift = -int(matmul(inverse_lattice, &
                    structure%positions(:,j) - structure%positions(:,i)))
                do nx = mins(1), maxs(1)
                    do ny = mins(2), maxs(2)
                        do nz = mins(3), maxs(3)
                            if (i == j .and. nx == 0 .and. ny == 0 .and. nz == 0) cycle
                            neighbor_position = structure%positions(:, j)
                            if (structure%pbc) neighbor_position = neighbor_position + &
                                matmul(structure%lattice, real(base_shift + [nx, ny, nz], real64))
                            displacement = neighbor_position - structure%positions(:, i)
                            if (sum(displacement*displacement) <= cutoff2) then
                                entry = entry + 1
                                neighbors%atom_indices(entry) = j
                                neighbors%image_shifts(:, entry) = base_shift + [nx, ny, nz]
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

    include 'legacy_chebyshev_evaluation.inc'

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
