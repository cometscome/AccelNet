module accelnet_behler
    use iso_fortran_env, only: error_unit, real64
    use accelnet_descriptors, only: cutoff_value, cutoff_derivative, cutoff_value_derivative, &
        validate_cutoff_parameters, CUTOFF_COS
    use accelnet_descriptors, only: sf_cutoff_value => cutoff_value, sf_cutoff_derivative => cutoff_derivative
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    real(real64), parameter :: EPS_DISTANCE = 1.0e-12_real64
    integer, parameter, public :: MIN_G5_MOMENT_NEIGHBORS = 16
    integer, parameter, public :: MAX_G5_MOMENT_ORDER = 10
    integer, parameter, public :: G5_EVALUATION_AUTO = 0
    integer, parameter, public :: G5_EVALUATION_DIRECT = 1
    integer, parameter, public :: G5_EVALUATION_MOMENT = 2
    integer, parameter, public :: G5_EVALUATION_MOMENT_FORCE = 3

    type :: g1_parameter
        integer :: species = 0, output = 0, cutoff_group = 0
        real(real64) :: rc = 0.0_real64
    end type
    type :: g2_parameter
        integer :: species = 0, output = 0, cutoff_group = 0, exponential_group = 0
        real(real64) :: rc = 0.0_real64, rs = 0.0_real64, eta = 0.0_real64
    end type
    type :: g3_parameter
        integer :: species = 0, output = 0, cutoff_group = 0
        real(real64) :: rc = 0.0_real64, kappa = 0.0_real64
    end type
    type :: angular_parameter
        integer :: species1 = 0, species2 = 0, output = 0
        integer :: cutoff_group = 0, exponential_group = 0
        integer :: next_in_pair = 0
        integer :: integer_zeta = 0
        real(real64) :: rc = 0.0_real64, lambda = 0.0_real64
        real(real64) :: zeta = 0.0_real64, eta = 0.0_real64, rs = 0.0_real64
        real(real64) :: derivative_prefactor = 0.0_real64
    end type
    type :: g4_parameter_group
        integer :: count = 0, angular_count = 0
        integer :: first_output = 0
        logical :: outputs_contiguous = .true.
        integer, allocatable :: output(:), cutoff_group(:), exponential_group(:), angular_group(:)
        ! Only combinations used by this species pair are evaluated for G4.
        integer, allocatable :: active_cutoffs(:), active_exponentials(:)
        integer, allocatable :: product_group(:), product_cutoff(:), product_exponential(:)
        real(real64), allocatable :: rc(:)
        integer, allocatable :: angular_integer_zeta(:)
        real(real64), allocatable :: angular_lambda(:), angular_zeta(:), angular_derivative_prefactor(:)
    end type
    type :: g1_parameter_group
        integer :: count = 0
        integer, allocatable :: output(:), cutoff_group(:)
        real(real64), allocatable :: rc(:)
    end type
    type :: g4_cutoff_group
        integer :: cutoff_group = 0
        type(g4_parameter_group) :: parameters
    end type
    type :: g4_zero_shift_pair
        type(g4_cutoff_group), allocatable :: cutoffs(:)
    end type
    type :: g2_parameter_group
        integer :: count = 0
        integer, allocatable :: output(:), cutoff_group(:), exponential_group(:)
        real(real64), allocatable :: rc(:), rs(:), eta(:)
    end type
    type :: g3_parameter_group
        integer :: count = 0
        integer, allocatable :: output(:), cutoff_group(:)
        real(real64), allocatable :: rc(:), kappa(:)
    end type

    type :: extended_parameter
        integer :: kind=0, species1=0, species2=0, output=0
        real(real64) :: p(7)=0
    end type

    type, public :: behler_config
        integer :: num_species = 0
        integer :: cutoff_type = CUTOFF_COS
        real(real64) :: cutoff_alpha = 0.0_real64
        integer :: number_of_descriptors = 0
        integer :: maximum_g5_integer_zeta = 0
        integer :: number_of_g5_moments = 0
        integer :: g5_evaluation_mode = G5_EVALUATION_AUTO
        logical :: g5_high_order_warning_emitted = .false.
        real(real64) :: maximum_cutoff = 0.0_real64
        real(real64) :: maximum_angular_cutoff = 0.0_real64
        real(real64) :: maximum_g4_cutoff = 0.0_real64
        real(real64), allocatable :: cutoff_radii(:)
        real(real64), allocatable :: radial_eta(:), radial_shift(:)
        real(real64), allocatable :: angular_eta(:), angular_shift(:)
        integer, allocatable :: g5_moment_x_power(:), g5_moment_y_power(:), g5_moment_z_power(:)
        integer, allocatable :: g5_moment_first(:), g5_moment_last(:)
        real(real64), allocatable :: g5_moment_multinomial(:)
        integer, allocatable :: g4_pair_first(:), g4_pair_last(:)
        integer, allocatable :: g5_pair_first(:), g5_pair_last(:)
        type(g4_parameter_group), allocatable :: g4_groups(:)
        type(g4_zero_shift_pair), allocatable :: g4_zero_shift_groups(:)
        logical :: has_zero_shift_g4 = .false., has_shifted_g4 = .false.
        type(g4_parameter_group), allocatable :: g5_groups(:)
        type(g1_parameter_group), allocatable :: g1_groups(:)
        type(g2_parameter_group), allocatable :: g2_groups(:)
        type(g3_parameter_group), allocatable :: g3_groups(:)
        type(g1_parameter), allocatable :: g1(:)
        type(g2_parameter), allocatable :: g2(:)
        type(g3_parameter), allocatable :: g3(:)
        type(angular_parameter), allocatable :: g4(:)
        type(angular_parameter), allocatable :: g5(:)
        type(extended_parameter), allocatable :: extended(:)
        real(real64), allocatable :: species_weights(:)
    contains
        procedure :: num_descriptors => behler_num_descriptors
    end type behler_config

    public :: initialize_behler_config, add_g1, add_g2, add_g3, add_g4, add_g5
    public :: evaluate_behler_values, evaluate_behler_values_derivatives
    public :: contract_behler_derivatives, behler_supports_direct_contraction
    public :: add_extended, compact_subtype, compact_subtype_name
    public :: set_behler_g5_evaluation

contains

    subroutine initialize_behler_config(config, num_species, cutoff_type, cutoff_alpha)
        type(behler_config), intent(out) :: config
        integer, intent(in) :: num_species
        integer, intent(in), optional :: cutoff_type
        real(real64), intent(in), optional :: cutoff_alpha
        if (num_species < 1) error stop "Behler num_species must be positive"
        config%num_species = num_species
        allocate(config%species_weights(num_species)); config%species_weights=1
        if (present(cutoff_type)) config%cutoff_type = cutoff_type
        if (present(cutoff_alpha)) config%cutoff_alpha = cutoff_alpha
        call validate_cutoff_parameters(config%cutoff_type, config%cutoff_alpha)
        allocate(config%g4_groups(num_species*(num_species + 1)/2))
        allocate(config%g4_zero_shift_groups(num_species*(num_species + 1)/2))
        allocate(config%g5_groups(num_species*(num_species + 1)/2))
        allocate(config%g1_groups(num_species), config%g2_groups(num_species), config%g3_groups(num_species))
    end subroutine

    integer function behler_num_descriptors(self) result(n)
        class(behler_config), intent(in) :: self
        n = self%number_of_descriptors
    end function

    subroutine set_behler_g5_evaluation(config, mode)
        type(behler_config), intent(inout) :: config
        integer, intent(in) :: mode
        if (mode < G5_EVALUATION_AUTO .or. mode > G5_EVALUATION_MOMENT_FORCE) &
            error stop "invalid G5 evaluation mode"
        config%g5_evaluation_mode = mode
    end subroutine set_behler_g5_evaluation

    pure logical function use_g5_moments(config, distances, for_force) result(use_moments)
        type(behler_config), intent(in) :: config
        real(real64), intent(in) :: distances(:)
        logical, intent(in) :: for_force
        use_moments = config%maximum_g5_integer_zeta > 0 .and. &
            config%g5_evaluation_mode /= G5_EVALUATION_DIRECT
        if (use_moments .and. config%g5_evaluation_mode /= G5_EVALUATION_MOMENT_FORCE) &
            use_moments = count(distances > EPS_DISTANCE .and. &
                distances <= config%maximum_angular_cutoff) >= MIN_G5_MOMENT_NEIGHBORS
    end function use_g5_moments

    subroutine validate_radial(config, species, rc)
        type(behler_config), intent(in) :: config
        integer, intent(in) :: species
        real(real64), intent(in) :: rc
        if (species < 1 .or. species > config%num_species) error stop "Behler species out of range"
        if (rc <= 0.0_real64) error stop "Behler cutoff must be positive"
    end subroutine

    subroutine validate_angular(config, species1, species2, rc)
        type(behler_config), intent(in) :: config
        integer, intent(in) :: species1, species2
        real(real64), intent(in) :: rc
        call validate_radial(config, species1, rc)
        call validate_radial(config, species2, rc)
    end subroutine

    subroutine add_g1(config, species, rc)
        type(behler_config), intent(inout) :: config
        integer, intent(in) :: species
        real(real64), intent(in) :: rc
        type(g1_parameter) :: parameter
        call validate_radial(config, species, rc)
        config%number_of_descriptors = config%number_of_descriptors + 1
        parameter%species = species
        parameter%output = config%number_of_descriptors
        parameter%rc = rc
        call find_or_add_value(config%cutoff_radii, rc, parameter%cutoff_group)
        if (allocated(config%g1)) then
            config%g1 = [config%g1, parameter]
        else
            config%g1 = [parameter]
        end if
        call append_g1_group(config%g1_groups(species), parameter)
        config%maximum_cutoff = max(config%maximum_cutoff, rc)
    end subroutine

    subroutine add_g2(config, species, rc, rs, eta)
        type(behler_config), intent(inout) :: config
        integer, intent(in) :: species
        real(real64), intent(in) :: rc, rs, eta
        type(g2_parameter) :: parameter
        call validate_radial(config, species, rc)
        config%number_of_descriptors = config%number_of_descriptors + 1
        parameter%species = species
        parameter%output = config%number_of_descriptors
        parameter%rc = rc
        parameter%rs = rs
        parameter%eta = eta
        call find_or_add_value(config%cutoff_radii, rc, parameter%cutoff_group)
        call find_or_add_pair(config%radial_eta, config%radial_shift, eta, rs, parameter%exponential_group)
        if (allocated(config%g2)) then
            config%g2 = [config%g2, parameter]
        else
            config%g2 = [parameter]
        end if
        call append_g2_group(config%g2_groups(species), parameter)
        config%maximum_cutoff = max(config%maximum_cutoff, rc)
    end subroutine

    subroutine add_g3(config, species, rc, kappa)
        type(behler_config), intent(inout) :: config
        integer, intent(in) :: species
        real(real64), intent(in) :: rc, kappa
        type(g3_parameter) :: parameter
        call validate_radial(config, species, rc)
        config%number_of_descriptors = config%number_of_descriptors + 1
        parameter%species = species
        parameter%output = config%number_of_descriptors
        parameter%rc = rc
        parameter%kappa = kappa
        call find_or_add_value(config%cutoff_radii, rc, parameter%cutoff_group)
        if (allocated(config%g3)) then
            config%g3 = [config%g3, parameter]
        else
            config%g3 = [parameter]
        end if
        call append_g3_group(config%g3_groups(species), parameter)
        config%maximum_cutoff = max(config%maximum_cutoff, rc)
    end subroutine

    subroutine add_g4(config, species1, species2, rc, lambda, zeta, eta, rs)
        type(behler_config), intent(inout) :: config
        integer, intent(in) :: species1, species2
        real(real64), intent(in) :: rc, lambda, zeta, eta
        real(real64), intent(in), optional :: rs
        type(angular_parameter) :: parameter
        call validate_angular(config, species1, species2, rc)
        config%number_of_descriptors = config%number_of_descriptors + 1
        parameter%species1 = species1
        parameter%species2 = species2
        parameter%output = config%number_of_descriptors
        parameter%rc = rc
        parameter%lambda = lambda
        parameter%zeta = zeta
        parameter%eta = eta
        if (present(rs)) parameter%rs = rs
        parameter%integer_zeta = integer_zeta_kind(zeta)
        parameter%derivative_prefactor = 0.5_real64*zeta*lambda
        call find_or_add_value(config%cutoff_radii, rc, parameter%cutoff_group)
        call find_or_add_pair(config%angular_eta, config%angular_shift, eta, parameter%rs, &
                              parameter%exponential_group)
        if (allocated(config%g4)) then
            config%g4 = [config%g4, parameter]
        else
            config%g4 = [parameter]
        end if
        call link_angular_pair(config%g4, config%g4_pair_first, config%g4_pair_last, &
                               config%num_species, species1, species2)
        ! Partition exactly at Rs=0. Nonzero shifts, including negative and
        ! very small shifts, retain the general kernel. Different cutoffs can
        ! still use the fast path in separate groups of the same species pair.
        if (parameter%rs == 0.0_real64) then
            call append_zero_shift_g4(config%g4_zero_shift_groups( &
                unordered_pair_index(config%num_species, species1, species2)), parameter)
            config%has_zero_shift_g4 = .true.
        else
            call append_g4_group(config%g4_groups(unordered_pair_index(config%num_species, species1, species2)), &
                                 parameter)
            config%has_shifted_g4 = .true.
        end if
        config%maximum_cutoff = max(config%maximum_cutoff, rc)
        config%maximum_angular_cutoff = max(config%maximum_angular_cutoff, rc)
        config%maximum_g4_cutoff = max(config%maximum_g4_cutoff, rc)
    end subroutine

    subroutine add_g5(config, species1, species2, rc, lambda, zeta, eta, rs)
        type(behler_config), intent(inout) :: config
        integer, intent(in) :: species1, species2
        real(real64), intent(in) :: rc, lambda, zeta, eta
        real(real64), intent(in), optional :: rs
        type(angular_parameter) :: parameter
        call validate_angular(config, species1, species2, rc)
        config%number_of_descriptors = config%number_of_descriptors + 1
        parameter%species1 = species1
        parameter%species2 = species2
        parameter%output = config%number_of_descriptors
        parameter%rc = rc
        parameter%lambda = lambda
        parameter%zeta = zeta
        parameter%eta = eta
        if (present(rs)) parameter%rs = rs
        parameter%integer_zeta = integer_zeta_kind(zeta)
        parameter%derivative_prefactor = 0.5_real64*zeta*lambda
        call find_or_add_value(config%cutoff_radii, rc, parameter%cutoff_group)
        call find_or_add_pair(config%angular_eta, config%angular_shift, eta, parameter%rs, &
                              parameter%exponential_group)
        if (allocated(config%g5)) then
            config%g5 = [config%g5, parameter]
        else
            config%g5 = [parameter]
        end if
        call link_angular_pair(config%g5, config%g5_pair_first, config%g5_pair_last, &
                               config%num_species, species1, species2)
        call append_g4_group(config%g5_groups(unordered_pair_index(config%num_species, species1, species2)), &
                             parameter)
        if (parameter%integer_zeta > 0 .and. parameter%integer_zeta <= MAX_G5_MOMENT_ORDER .and. &
            parameter%integer_zeta > config%maximum_g5_integer_zeta) &
            call initialize_g5_moment_basis(config, parameter%integer_zeta)
        if (parameter%integer_zeta > MAX_G5_MOMENT_ORDER .and. .not. config%g5_high_order_warning_emitted) then
            write(error_unit, "(A,I0,A,I0,A)") "WARNING: G5 zeta=", parameter%integer_zeta, &
                " exceeds the legacy atomic moment maximum order ", MAX_G5_MOMENT_ORDER, &
                "; the legacy atomic evaluator uses direct evaluation (shared batch policy is separate)."
            config%g5_high_order_warning_emitted = .true.
        end if
        config%maximum_cutoff = max(config%maximum_cutoff, rc)
        config%maximum_angular_cutoff = max(config%maximum_angular_cutoff, rc)
    end subroutine

    subroutine initialize_g5_moment_basis(config, maximum_zeta)
        type(behler_config), intent(inout) :: config
        integer, intent(in) :: maximum_zeta
        integer :: q, a, b, c, entry, maximum_moments
        if (maximum_zeta <= config%maximum_g5_integer_zeta) return
        maximum_moments = (maximum_zeta + 1)*(maximum_zeta + 2)*(maximum_zeta + 3)/6
        if (allocated(config%g5_moment_x_power)) then
            deallocate(config%g5_moment_x_power, config%g5_moment_y_power, config%g5_moment_z_power, &
                       config%g5_moment_first, config%g5_moment_last, config%g5_moment_multinomial)
        end if
        allocate(config%g5_moment_x_power(maximum_moments), &
                 config%g5_moment_y_power(maximum_moments), &
                 config%g5_moment_z_power(maximum_moments), &
                 config%g5_moment_multinomial(maximum_moments), &
                 config%g5_moment_first(maximum_zeta + 1), &
                 config%g5_moment_last(maximum_zeta + 1))
        entry = 0
        do q = 0, maximum_zeta
            config%g5_moment_first(q + 1) = entry + 1
            do a = 0, q
                do b = 0, q - a
                    c = q - a - b
                    entry = entry + 1
                    config%g5_moment_x_power(entry) = a
                    config%g5_moment_y_power(entry) = b
                    config%g5_moment_z_power(entry) = c
                    config%g5_moment_multinomial(entry) = factorial_real(q)/ &
                        (factorial_real(a)*factorial_real(b)*factorial_real(c))
                end do
            end do
            config%g5_moment_last(q + 1) = entry
        end do
        config%maximum_g5_integer_zeta = maximum_zeta
        config%number_of_g5_moments = entry
    end subroutine initialize_g5_moment_basis

    pure real(real64) function factorial_real(n) result(value)
        integer, intent(in) :: n
        integer :: i
        value = 1.0_real64
        do i = 2, n
            value = value*real(i, real64)
        end do
    end function factorial_real

    pure real(real64) function binomial_real(n, k) result(value)
        integer, intent(in) :: n, k
        value = factorial_real(n)/(factorial_real(k)*factorial_real(n - k))
    end function binomial_real

    pure real(real64) function moment_monomial(unit_vector, x_power, y_power, z_power) result(value)
        real(real64), intent(in) :: unit_vector(3)
        integer, intent(in) :: x_power, y_power, z_power
        value = unit_vector(1)**x_power*unit_vector(2)**y_power*unit_vector(3)**z_power
    end function moment_monomial

    pure subroutine moment_monomial_gradient(unit_vector, x_power, y_power, z_power, gradient)
        real(real64), intent(in) :: unit_vector(3)
        integer, intent(in) :: x_power, y_power, z_power
        real(real64), intent(out) :: gradient(3)
        gradient = 0.0_real64
        if (x_power > 0) gradient(1) = real(x_power, real64)*unit_vector(1)**(x_power - 1)* &
            unit_vector(2)**y_power*unit_vector(3)**z_power
        if (y_power > 0) gradient(2) = real(y_power, real64)*unit_vector(1)**x_power* &
            unit_vector(2)**(y_power - 1)*unit_vector(3)**z_power
        if (z_power > 0) gradient(3) = real(z_power, real64)*unit_vector(1)**x_power* &
            unit_vector(2)**y_power*unit_vector(3)**(z_power - 1)
    end subroutine moment_monomial_gradient

    subroutine find_or_add_value(values, candidate, index)
        real(real64), allocatable, intent(inout) :: values(:)
        real(real64), intent(in) :: candidate
        integer, intent(out) :: index
        integer :: i
        if (allocated(values)) then
            do i = 1, size(values)
                if (values(i) == candidate) then
                    index = i
                    return
                end if
            end do
            values = [values, candidate]
        else
            values = [candidate]
        end if
        index = size(values)
    end subroutine find_or_add_value

    subroutine find_or_add_pair(first, second, candidate_first, candidate_second, index)
        real(real64), allocatable, intent(inout) :: first(:), second(:)
        real(real64), intent(in) :: candidate_first, candidate_second
        integer, intent(out) :: index
        integer :: i
        if (allocated(first)) then
            do i = 1, size(first)
                if (first(i) == candidate_first .and. second(i) == candidate_second) then
                    index = i
                    return
                end if
            end do
            first = [first, candidate_first]
            second = [second, candidate_second]
        else
            first = [candidate_first]
            second = [candidate_second]
        end if
        index = size(first)
    end subroutine find_or_add_pair

    pure integer function unordered_pair_index(num_species, species1, species2) result(index)
        integer, intent(in) :: num_species, species1, species2
        integer :: low_species, high_species, first
        low_species = min(species1, species2)
        high_species = max(species1, species2)
        index = 0
        do first = 1, low_species - 1
            index = index + num_species - first + 1
        end do
        index = index + high_species - low_species + 1
    end function unordered_pair_index

    subroutine link_angular_pair(parameters, first_for_pair, last_for_pair, num_species, species1, species2)
        type(angular_parameter), intent(inout) :: parameters(:)
        integer, allocatable, intent(inout) :: first_for_pair(:), last_for_pair(:)
        integer, intent(in) :: num_species, species1, species2
        integer :: pair, current, number_of_pairs
        number_of_pairs = num_species*(num_species + 1)/2
        if (.not. allocated(first_for_pair)) then
            allocate(first_for_pair(number_of_pairs), source=0)
            allocate(last_for_pair(number_of_pairs), source=0)
        end if
        pair = unordered_pair_index(num_species, species1, species2)
        current = size(parameters)
        if (last_for_pair(pair) == 0) then
            first_for_pair(pair) = current
        else
            parameters(last_for_pair(pair))%next_in_pair = current
        end if
        last_for_pair(pair) = current
    end subroutine link_angular_pair

    subroutine append_zero_shift_g4(pair, parameter)
        type(g4_zero_shift_pair), intent(inout) :: pair
        type(angular_parameter), intent(in) :: parameter
        type(g4_cutoff_group) :: new_group
        integer :: i
        if (allocated(pair%cutoffs)) then
            do i = 1, size(pair%cutoffs)
                if (pair%cutoffs(i)%cutoff_group /= parameter%cutoff_group) cycle
                call append_g4_group(pair%cutoffs(i)%parameters, parameter)
                return
            end do
        end if
        new_group%cutoff_group = parameter%cutoff_group
        call append_g4_group(new_group%parameters, parameter)
        if (allocated(pair%cutoffs)) then
            pair%cutoffs = [pair%cutoffs, new_group]
        else
            pair%cutoffs = [new_group]
        end if
    end subroutine append_zero_shift_g4

    subroutine append_g4_group(group, parameter)
        type(g4_parameter_group), intent(inout) :: group
        type(angular_parameter), intent(in) :: parameter
        integer :: i, angular_index, product_index
        call append_unique_index(group%active_cutoffs, parameter%cutoff_group)
        call append_unique_index(group%active_exponentials, parameter%exponential_group)
        product_index = 0
        if (allocated(group%product_cutoff)) then
            do i = 1, size(group%product_cutoff)
                if (group%product_cutoff(i) == parameter%cutoff_group .and. &
                    group%product_exponential(i) == parameter%exponential_group) then
                    product_index = i
                    exit
                end if
            end do
        end if
        if (product_index == 0) then
            if (allocated(group%product_cutoff)) then
                group%product_cutoff = [group%product_cutoff, parameter%cutoff_group]
                group%product_exponential = [group%product_exponential, parameter%exponential_group]
            else
                group%product_cutoff = [parameter%cutoff_group]
                group%product_exponential = [parameter%exponential_group]
            end if
            product_index = size(group%product_cutoff)
        end if
        if (allocated(group%product_group)) then
            group%product_group = [group%product_group, product_index]
        else
            group%product_group = [product_index]
        end if
        angular_index = 0
        do i = 1, group%angular_count
            if (group%angular_lambda(i) == parameter%lambda .and. group%angular_zeta(i) == parameter%zeta) then
                angular_index = i
                exit
            end if
        end do
        if (angular_index == 0) then
            angular_index = group%angular_count + 1
            if (group%angular_count == 0) then
                group%angular_lambda = [parameter%lambda]
                group%angular_zeta = [parameter%zeta]
                group%angular_integer_zeta = [parameter%integer_zeta]
                group%angular_derivative_prefactor = [parameter%derivative_prefactor]
            else
                group%angular_lambda = [group%angular_lambda, parameter%lambda]
                group%angular_zeta = [group%angular_zeta, parameter%zeta]
                group%angular_integer_zeta = [group%angular_integer_zeta, parameter%integer_zeta]
                group%angular_derivative_prefactor = &
                    [group%angular_derivative_prefactor, parameter%derivative_prefactor]
            end if
            group%angular_count = angular_index
        end if
        if (group%count == 0) then
            group%first_output = parameter%output
            group%outputs_contiguous = .true.
            group%output = [parameter%output]
            group%cutoff_group = [parameter%cutoff_group]
            group%exponential_group = [parameter%exponential_group]
            group%angular_group = [angular_index]
            group%rc = [parameter%rc]
        else
            if (parameter%output /= group%first_output + group%count) group%outputs_contiguous = .false.
            group%output = [group%output, parameter%output]
            group%cutoff_group = [group%cutoff_group, parameter%cutoff_group]
            group%exponential_group = [group%exponential_group, parameter%exponential_group]
            group%angular_group = [group%angular_group, angular_index]
            group%rc = [group%rc, parameter%rc]
        end if
        group%count = group%count + 1
    end subroutine append_g4_group

    subroutine append_unique_index(indices, index)
        integer, allocatable, intent(inout) :: indices(:)
        integer, intent(in) :: index
        if (allocated(indices)) then
            if (any(indices == index)) return
            indices = [indices, index]
        else
            indices = [index]
        end if
    end subroutine append_unique_index

    subroutine append_g1_group(group, parameter)
        type(g1_parameter_group), intent(inout) :: group
        type(g1_parameter), intent(in) :: parameter
        if (group%count == 0) then
            group%output = [parameter%output]
            group%cutoff_group = [parameter%cutoff_group]
            group%rc = [parameter%rc]
        else
            group%output = [group%output, parameter%output]
            group%cutoff_group = [group%cutoff_group, parameter%cutoff_group]
            group%rc = [group%rc, parameter%rc]
        end if
        group%count = group%count + 1
    end subroutine append_g1_group

    subroutine append_g2_group(group, parameter)
        type(g2_parameter_group), intent(inout) :: group
        type(g2_parameter), intent(in) :: parameter
        if (group%count == 0) then
            group%output = [parameter%output]
            group%cutoff_group = [parameter%cutoff_group]
            group%exponential_group = [parameter%exponential_group]
            group%rc = [parameter%rc]
            group%rs = [parameter%rs]
            group%eta = [parameter%eta]
        else
            group%output = [group%output, parameter%output]
            group%cutoff_group = [group%cutoff_group, parameter%cutoff_group]
            group%exponential_group = [group%exponential_group, parameter%exponential_group]
            group%rc = [group%rc, parameter%rc]
            group%rs = [group%rs, parameter%rs]
            group%eta = [group%eta, parameter%eta]
        end if
        group%count = group%count + 1
    end subroutine append_g2_group

    subroutine append_g3_group(group, parameter)
        type(g3_parameter_group), intent(inout) :: group
        type(g3_parameter), intent(in) :: parameter
        if (group%count == 0) then
            group%output = [parameter%output]
            group%cutoff_group = [parameter%cutoff_group]
            group%rc = [parameter%rc]
            group%kappa = [parameter%kappa]
        else
            group%output = [group%output, parameter%output]
            group%cutoff_group = [group%cutoff_group, parameter%cutoff_group]
            group%rc = [group%rc, parameter%rc]
            group%kappa = [group%kappa, parameter%kappa]
        end if
        group%count = group%count + 1
    end subroutine append_g3_group

    pure logical function species_pair_matches(first, second, expected1, expected2)
        integer, intent(in) :: first, second, expected1, expected2
        species_pair_matches = (first == expected1 .and. second == expected2) .or. &
                               (first == expected2 .and. second == expected1)
    end function

    pure integer function allocated_size(values) result(n)
        real(real64), allocatable, intent(in) :: values(:)
        n = 0
        if (allocated(values)) n = size(values)
    end function allocated_size

    pure integer function integer_zeta_kind(zeta) result(value)
        real(real64), intent(in) :: zeta
        integer :: candidate
        candidate = nint(zeta)
        value = 0
        if (candidate >= 1 .and. abs(zeta - real(candidate, real64)) < 1.0e-10_real64) value = candidate
    end function integer_zeta_kind

    pure subroutine angular_power(cosine, lambda, zeta, integer_zeta, derivative_prefactor, value, derivative)
        include 'angular_power.inc'
    end subroutine

    pure function angular_value(cosine, lambda, zeta, integer_zeta) result(value)
        include 'angular_value.inc'
    end function angular_value

    include 'legacy_behler_evaluation.inc'

    integer function compact_subtype(name) result(code)
        character(len=*), intent(in) :: name
        integer :: ios
        code=0
        if (trim(name) == 'e') then
            code=5
        else if (len_trim(name) == 2 .or. len_trim(name) == 3) then
            if (name(1:1) /= 'p') error stop 'invalid compact subtype'
            read(name(2:2),*,iostat=ios) code
            if (ios /= 0) error stop 'invalid compact subtype'
            if (code < 1 .or. code > 4) error stop 'invalid compact polynomial'
            if (len_trim(name) == 3) then
                if (name(3:3) /= 'a') error stop 'invalid compact asymmetry'
                code=code+10
            end if
        end if
        if (code == 0) error stop 'invalid compact subtype'
    end function

    function compact_subtype_name(code) result(name)
        integer, intent(in) :: code
        character(len=3) :: name
        name='e'
        if (code /= 5) then
            write(name,'(A,I1)') 'p',mod(code,10)
            if (code > 10) name(3:3)='a'
        end if
    end function

    subroutine add_extended(config,kind,species1,species2,parameters)
        type(behler_config), intent(inout) :: config
        integer, intent(in) :: kind,species1,species2
        real(real64), intent(in) :: parameters(7)
        type(extended_parameter) :: ep
        integer :: subtype
        real(real64) :: middle
        if (.not. (kind == 12 .or. kind == 13 .or. (kind >= 20 .and. kind <= 25))) &
            error stop 'invalid extended symmetry function'
        if (.not. all(ieee_is_finite(parameters))) error stop 'nonfinite extended parameters'
        if (parameters(1) <= 0) error stop 'invalid extended cutoff'
        if (kind >= 20 .and. kind <= 22) then
            call validate_radial(config,species1,parameters(1))
            if (kind /= 20) call validate_radial(config,species2,parameters(1))
        else
            if (species1 /= 0 .or. species2 /= 0) error stop 'weighted SF must include all species'
        end if
        ep%kind=kind; ep%species1=species1; ep%species2=species2; ep%p=parameters
        if (kind == 13) then
            if (abs(parameters(4)) > 1 .or. parameters(5) < 1) error stop 'invalid weighted angular power'
        end if
        if (kind >= 20) then
            if (parameters(2) >= parameters(1)) error stop 'empty compact radial interval'
            subtype=nint(parameters(5))
            if (parameters(5) /= real(subtype,real64)) error stop 'invalid compact subtype code'
            if (.not. (subtype >= 1 .and. subtype <= 5) .and. &
                .not. (subtype >= 11 .and. subtype <= 14)) error stop 'invalid compact subtype code'
            if (kind /= 20 .and. kind /= 23) then
                middle=0.5_real64*(parameters(3)+parameters(4))
                if (parameters(3) >= parameters(4) .or. parameters(4)-parameters(3) > 360) &
                    error stop 'invalid compact angle interval'
                if ((parameters(3) < 0 .and. middle /= 0) .or. &
                    (parameters(4) > 180 .and. middle /= 180)) error stop 'invalid compact angle center'
                ep%p(3:4)=parameters(3:4)*(acos(-1.0_real64)/180)
            end if
        end if
        config%number_of_descriptors=config%number_of_descriptors+1
        ep%output=config%number_of_descriptors
        if (allocated(config%extended)) then
            config%extended=[config%extended,ep]
        else
            config%extended=[ep]
        end if
        config%maximum_cutoff=max(config%maximum_cutoff,parameters(1))
    end subroutine

    include 'n2p2_extended.inc'

end module accelnet_behler
