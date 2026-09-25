module accelnet_behler
    use iso_fortran_env, only: error_unit, real64
    use accelnet_descriptors, only: cutoff_value, cutoff_derivative, cutoff_value_derivative, &
        validate_cutoff_parameters, CUTOFF_COS
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
    contains
        procedure :: num_descriptors => behler_num_descriptors
    end type behler_config

    public :: initialize_behler_config, add_g1, add_g2, add_g3, add_g4, add_g5
    public :: evaluate_behler_values, evaluate_behler_values_derivatives
    public :: contract_behler_derivatives, behler_supports_direct_contraction
    public :: set_behler_g5_evaluation

contains

    subroutine initialize_behler_config(config, num_species, cutoff_type, cutoff_alpha)
        type(behler_config), intent(out) :: config
        integer, intent(in) :: num_species
        integer, intent(in), optional :: cutoff_type
        real(real64), intent(in), optional :: cutoff_alpha
        if (num_species < 1) error stop "Behler num_species must be positive"
        config%num_species = num_species
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
                " exceeds the moment maximum order ", MAX_G5_MOMENT_ORDER, &
                "; high-order descriptors will use direct evaluation."
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
        real(real64), intent(in) :: cosine, lambda, zeta, derivative_prefactor
        integer, intent(in) :: integer_zeta
        real(real64), intent(out) :: value, derivative
        real(real64) :: base
        base = 0.5_real64*(1.0_real64 + lambda*cosine)
        if (integer_zeta > 0) then
            select case(integer_zeta)
            case(1)
                value = base
                derivative = derivative_prefactor
            case(2)
                value = base*base
                derivative = derivative_prefactor*base
            case(3)
                value = base*base*base
                derivative = derivative_prefactor*base*base
            case(4)
                value = (base*base)*(base*base)
                derivative = derivative_prefactor*base*base*base
            case default
                value = base**integer_zeta
                derivative = derivative_prefactor*base**(integer_zeta - 1)
            end select
        else
            value = base**zeta
            derivative = derivative_prefactor*base**(zeta - 1.0_real64)
        end if
    end subroutine

    pure function angular_value(cosine, lambda, zeta, integer_zeta) result(value)
        real(real64), intent(in) :: cosine, lambda, zeta
        integer, intent(in) :: integer_zeta
        real(real64) :: value, base
        base = 0.5_real64*(1.0_real64 + lambda*cosine)
        if (integer_zeta > 0) then
            select case(integer_zeta)
            case(1)
                value = base
            case(2)
                value = base*base
            case(3)
                value = base*base*base
            case(4)
                value = (base*base)*(base*base)
            case default
                value = base**integer_zeta
            end select
        else
            value = base**zeta
        end if
    end function angular_value

    subroutine evaluate_behler_values(config, displacements, neighbor_species, values)
        type(behler_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: values(:)
        real(real64) :: center_derivative(0, 0), neighbor_derivatives(0, 0, 0)

        call evaluate_core(config, displacements, neighbor_species, values, .false., &
                           center_derivative, neighbor_derivatives)
    end subroutine

    subroutine evaluate_behler_values_derivatives(config, displacements, neighbor_species, values, &
                                                  derivative_center, derivative_neighbors)
        type(behler_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: values(:)
        real(real64), intent(out) :: derivative_center(:, :)
        real(real64), intent(out) :: derivative_neighbors(:, :, :)
        call evaluate_core(config, displacements, neighbor_species, values, .true., &
                           derivative_center, derivative_neighbors)
    end subroutine

    pure logical function behler_supports_direct_contraction(config) result(supported)
        type(behler_config), intent(in) :: config
        supported = .not. allocated(config%g4)
    end function behler_supports_direct_contraction

    subroutine contract_behler_derivatives(config, displacements, neighbor_species, coefficients, &
                                           contracted_center, contracted_neighbors)
        type(behler_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :), coefficients(:)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: contracted_center(3), contracted_neighbors(:, :)
        integer :: nneighbors, ncutoff, nradial_exp, nangular_exp
        integer :: j, k, p, pair, group, cutoff_index, exponential_index, descriptor
        real(real64) :: distance_j, distance_k, cutoff_j, cutoff_k, exponential
        real(real64) :: radial_value, radial_derivative, cosine_jk, product
        real(real64) :: angular_term, angular_derivative, weight
        real(real64) :: vector_j(3), vector_k(3)
        real(real64) :: derivative(3)
        real(real64) :: dcos_j1, dcos_j2, dcos_j3, dcos_k1, dcos_k2, dcos_k3
        real(real64) :: angular_coefficient, radial_coefficient_j, radial_coefficient_k
        real(real64) :: pair_j1, pair_j2, pair_j3, pair_k1, pair_k2, pair_k3
        real(real64) :: distances(size(neighbor_species)), unit_vectors(3, size(neighbor_species))
        real(real64) :: cutoffs(size(neighbor_species), max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: cutoff_derivatives(size(neighbor_species), max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: radial_exponentials(size(neighbor_species), max(1, allocated_size(config%radial_eta)))
        real(real64) :: angular_exponentials(size(neighbor_species), max(1, allocated_size(config%angular_eta)))
        real(real64) :: products(max(1, allocated_size(config%angular_eta)), &
                                 max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: radial_j(max(1, allocated_size(config%angular_eta)), &
                                max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: radial_k(max(1, allocated_size(config%angular_eta)), &
                                max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: angular_values(max(1, allocated_size_angular(config%g5)))
        real(real64) :: angular_derivatives(max(1, allocated_size_angular(config%g5)))
        logical :: use_integer_moments

        if (.not. behler_supports_direct_contraction(config)) &
            error stop "direct Behler contraction does not support G4"
        if (size(coefficients) < config%number_of_descriptors) &
            error stop "Behler coefficient array too small"
        if (size(contracted_neighbors, 1) /= 3 .or. &
            size(contracted_neighbors, 2) < size(neighbor_species)) &
            error stop "Behler contracted array too small"

        nneighbors = size(neighbor_species)
        ncutoff = allocated_size(config%cutoff_radii)
        nradial_exp = allocated_size(config%radial_eta)
        nangular_exp = allocated_size(config%angular_eta)
        contracted_center = 0.0_real64
        contracted_neighbors(:, 1:nneighbors) = 0.0_real64
        unit_vectors = 0.0_real64

        do j = 1, nneighbors
            distances(j) = sqrt(dot_product(displacements(:, j), displacements(:, j)))
            if (distances(j) > EPS_DISTANCE) unit_vectors(:, j) = displacements(:, j)/distances(j)
            do group = 1, ncutoff
                cutoffs(j, group) = cutoff_value(distances(j), config%cutoff_radii(group), &
                    config%cutoff_type, config%cutoff_alpha)
                cutoff_derivatives(j, group) = cutoff_derivative(distances(j), config%cutoff_radii(group), &
                    config%cutoff_type, config%cutoff_alpha)
            end do
            do group = 1, nradial_exp
                radial_exponentials(j, group) = exp(-config%radial_eta(group)* &
                    (distances(j) - config%radial_shift(group))**2)
            end do
            do group = 1, nangular_exp
                angular_exponentials(j, group) = exp(-config%angular_eta(group)* &
                    (distances(j) - config%angular_shift(group))**2)
            end do
        end do

        do j = 1, nneighbors
            distance_j = distances(j)
            if (distance_j <= EPS_DISTANCE) cycle
            vector_j = unit_vectors(:, j)
            associate(parameters => config%g1_groups(neighbor_species(j)))
            do p = 1, parameters%count
                if (distance_j > parameters%rc(p)) cycle
                derivative = coefficients(parameters%output(p))* &
                    cutoff_derivatives(j, parameters%cutoff_group(p))*vector_j
                contracted_neighbors(:, j) = contracted_neighbors(:, j) + derivative
                contracted_center = contracted_center - derivative
            end do
            end associate
            associate(parameters => config%g2_groups(neighbor_species(j)))
            do p = 1, parameters%count
                if (distance_j > parameters%rc(p)) cycle
                cutoff_j = cutoffs(j, parameters%cutoff_group(p))
                exponential = radial_exponentials(j, parameters%exponential_group(p))
                radial_derivative = exponential*(cutoff_derivatives(j, parameters%cutoff_group(p)) - &
                    2.0_real64*parameters%eta(p)*(distance_j - parameters%rs(p))*cutoff_j)
                derivative = coefficients(parameters%output(p))*radial_derivative*vector_j
                contracted_neighbors(:, j) = contracted_neighbors(:, j) + derivative
                contracted_center = contracted_center - derivative
            end do
            end associate
            associate(parameters => config%g3_groups(neighbor_species(j)))
            do p = 1, parameters%count
                if (distance_j > parameters%rc(p)) cycle
                cutoff_j = cutoffs(j, parameters%cutoff_group(p))
                radial_value = cos(parameters%kappa(p)*distance_j)
                radial_derivative = cutoff_derivatives(j, parameters%cutoff_group(p))*radial_value - &
                    parameters%kappa(p)*cutoff_j*sin(parameters%kappa(p)*distance_j)
                derivative = coefficients(parameters%output(p))*radial_derivative*vector_j
                contracted_neighbors(:, j) = contracted_neighbors(:, j) + derivative
                contracted_center = contracted_center - derivative
            end do
            end associate
        end do

        if (.not. allocated(config%g5)) return
        use_integer_moments = use_g5_moments(config, distances, .true.)
        if (use_integer_moments) call contract_g5_integer_moments(config, distances, unit_vectors, &
            neighbor_species, cutoffs, cutoff_derivatives, angular_exponentials, coefficients, &
            contracted_center, contracted_neighbors)
        do j = 1, nneighbors
            distance_j = distances(j)
            if (distance_j <= EPS_DISTANCE .or. distance_j > config%maximum_angular_cutoff) cycle
            vector_j = unit_vectors(:, j)
            do k = j + 1, nneighbors
                distance_k = distances(k)
                if (distance_k <= EPS_DISTANCE .or. distance_k > config%maximum_angular_cutoff) cycle
                pair = unordered_pair_index(config%num_species, neighbor_species(j), neighbor_species(k))
                associate(parameters => config%g5_groups(pair))
                if (parameters%count == 0) cycle
                vector_k = unit_vectors(:, k)
                cosine_jk = max(-1.0_real64, min(1.0_real64, dot_product(vector_j, vector_k)))
                dcos_j1 = (vector_k(1) - cosine_jk*vector_j(1))/distance_j
                dcos_j2 = (vector_k(2) - cosine_jk*vector_j(2))/distance_j
                dcos_j3 = (vector_k(3) - cosine_jk*vector_j(3))/distance_j
                dcos_k1 = (vector_j(1) - cosine_jk*vector_k(1))/distance_k
                dcos_k2 = (vector_j(2) - cosine_jk*vector_k(2))/distance_k
                dcos_k3 = (vector_j(3) - cosine_jk*vector_k(3))/distance_k
                do exponential_index = 1, nangular_exp
                    do cutoff_index = 1, ncutoff
                        cutoff_j = cutoffs(j, cutoff_index)
                        cutoff_k = cutoffs(k, cutoff_index)
                        products(exponential_index, cutoff_index) = &
                            angular_exponentials(j, exponential_index)*cutoff_j* &
                            angular_exponentials(k, exponential_index)*cutoff_k
                        radial_derivative = angular_exponentials(j, exponential_index)* &
                            (cutoff_derivatives(j, cutoff_index) - &
                             2.0_real64*config%angular_eta(exponential_index)* &
                             (distance_j - config%angular_shift(exponential_index))*cutoff_j)
                        radial_j(exponential_index, cutoff_index) = radial_derivative* &
                            angular_exponentials(k, exponential_index)*cutoff_k
                        radial_derivative = angular_exponentials(k, exponential_index)* &
                            (cutoff_derivatives(k, cutoff_index) - &
                             2.0_real64*config%angular_eta(exponential_index)* &
                             (distance_k - config%angular_shift(exponential_index))*cutoff_k)
                        radial_k(exponential_index, cutoff_index) = &
                            angular_exponentials(j, exponential_index)*cutoff_j*radial_derivative
                    end do
                end do
                do group = 1, parameters%angular_count
                    if (use_integer_moments .and. parameters%angular_integer_zeta(group) > 0 .and. &
                        parameters%angular_integer_zeta(group) <= MAX_G5_MOMENT_ORDER) cycle
                    call angular_power(cosine_jk, parameters%angular_lambda(group), &
                        parameters%angular_zeta(group), parameters%angular_integer_zeta(group), &
                        parameters%angular_derivative_prefactor(group), angular_values(group), &
                        angular_derivatives(group))
                end do
                pair_j1 = 0.0_real64; pair_j2 = 0.0_real64; pair_j3 = 0.0_real64
                pair_k1 = 0.0_real64; pair_k2 = 0.0_real64; pair_k3 = 0.0_real64
                do p = 1, parameters%count
                    group = parameters%angular_group(p)
                    if (use_integer_moments .and. parameters%angular_integer_zeta(group) > 0 .and. &
                        parameters%angular_integer_zeta(group) <= MAX_G5_MOMENT_ORDER) cycle
                    if (distance_j > parameters%rc(p) .or. distance_k > parameters%rc(p)) cycle
                    cutoff_index = parameters%cutoff_group(p)
                    exponential_index = parameters%exponential_group(p)
                    descriptor = parameters%output(p)
                    product = products(exponential_index, cutoff_index)
                    angular_term = angular_values(group)
                    angular_derivative = angular_derivatives(group)
                    weight = 2.0_real64*coefficients(descriptor)
                    angular_coefficient = weight*angular_derivative*product
                    radial_coefficient_j = weight*angular_term*radial_j(exponential_index, cutoff_index)
                    radial_coefficient_k = weight*angular_term*radial_k(exponential_index, cutoff_index)
                    pair_j1 = pair_j1 + angular_coefficient*dcos_j1 + radial_coefficient_j*vector_j(1)
                    pair_j2 = pair_j2 + angular_coefficient*dcos_j2 + radial_coefficient_j*vector_j(2)
                    pair_j3 = pair_j3 + angular_coefficient*dcos_j3 + radial_coefficient_j*vector_j(3)
                    pair_k1 = pair_k1 + angular_coefficient*dcos_k1 + radial_coefficient_k*vector_k(1)
                    pair_k2 = pair_k2 + angular_coefficient*dcos_k2 + radial_coefficient_k*vector_k(2)
                    pair_k3 = pair_k3 + angular_coefficient*dcos_k3 + radial_coefficient_k*vector_k(3)
                end do
                contracted_neighbors(1, j) = contracted_neighbors(1, j) + pair_j1
                contracted_neighbors(2, j) = contracted_neighbors(2, j) + pair_j2
                contracted_neighbors(3, j) = contracted_neighbors(3, j) + pair_j3
                contracted_neighbors(1, k) = contracted_neighbors(1, k) + pair_k1
                contracted_neighbors(2, k) = contracted_neighbors(2, k) + pair_k2
                contracted_neighbors(3, k) = contracted_neighbors(3, k) + pair_k3
                contracted_center(1) = contracted_center(1) - pair_j1 - pair_k1
                contracted_center(2) = contracted_center(2) - pair_j2 - pair_k2
                contracted_center(3) = contracted_center(3) - pair_j3 - pair_k3
                end associate
            end do
        end do
    end subroutine contract_behler_derivatives

    subroutine contract_g5_integer_moments(config, distances, unit_vectors, neighbor_species, cutoffs, &
                                           cutoff_derivatives, angular_exponentials, coefficients, &
                                           contracted_center, contracted_neighbors)
        type(behler_config), intent(in) :: config
        real(real64), intent(in) :: distances(:), unit_vectors(:, :), cutoffs(:, :)
        real(real64), intent(in) :: cutoff_derivatives(:, :), angular_exponentials(:, :), coefficients(:)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(inout) :: contracted_center(3), contracted_neighbors(:, :)
        integer :: nneighbors, ncutoff, nangular_exp
        integer :: j, p, entry, q, species, species1, species2, cutoff_index, exponential_index
        integer :: x_power, y_power, z_power
        real(real64) :: h, polynomial_coefficient, radial_derivative, adjoint, self_scale
        real(real64) :: ux, uy, uz, gradient_x, gradient_y, gradient_z, gradient_dot_u
        real(real64) :: dhx, dhy, dhz, dphix, dphiy, dphiz
        real(real64) :: derivative_x, derivative_y, derivative_z
        real(real64) :: moments(max(1, config%number_of_g5_moments), &
                                max(1, allocated_size(config%cutoff_radii)), &
                                max(1, allocated_size(config%angular_eta)), config%num_species)
        real(real64) :: moment_adjoints(max(1, config%number_of_g5_moments), &
                                        max(1, allocated_size(config%cutoff_radii)), &
                                        max(1, allocated_size(config%angular_eta)), config%num_species)
        real(real64) :: self_adjoints(max(1, allocated_size(config%cutoff_radii)), &
                                      max(1, allocated_size(config%angular_eta)), config%num_species)
        real(real64) :: monomials(max(1, config%number_of_g5_moments))
        real(real64) :: monomial_gradient_x(max(1, config%number_of_g5_moments))
        real(real64) :: monomial_gradient_y(max(1, config%number_of_g5_moments))
        real(real64) :: monomial_gradient_z(max(1, config%number_of_g5_moments))

        nneighbors = size(neighbor_species)
        ncutoff = allocated_size(config%cutoff_radii)
        nangular_exp = allocated_size(config%angular_eta)
        moments = 0.0_real64
        do j = 1, nneighbors
            if (distances(j) <= EPS_DISTANCE .or. distances(j) > config%maximum_angular_cutoff) cycle
            do entry = 1, config%number_of_g5_moments
                monomials(entry) = moment_monomial(unit_vectors(:, j), config%g5_moment_x_power(entry), &
                    config%g5_moment_y_power(entry), config%g5_moment_z_power(entry))
            end do
            do exponential_index = 1, nangular_exp
                do cutoff_index = 1, ncutoff
                    h = angular_exponentials(j, exponential_index)*cutoffs(j, cutoff_index)
                    if (h == 0.0_real64) cycle
                    do entry = 1, config%number_of_g5_moments
                        moments(entry, cutoff_index, exponential_index, neighbor_species(j)) = &
                            moments(entry, cutoff_index, exponential_index, neighbor_species(j)) + &
                            h*monomials(entry)
                    end do
                end do
            end do
        end do

        moment_adjoints = 0.0_real64
        self_adjoints = 0.0_real64
        do p = 1, size(config%g5)
            if (config%g5(p)%integer_zeta <= 0 .or. &
                config%g5(p)%integer_zeta > MAX_G5_MOMENT_ORDER) cycle
            species1 = config%g5(p)%species1
            species2 = config%g5(p)%species2
            cutoff_index = config%g5(p)%cutoff_group
            exponential_index = config%g5(p)%exponential_group
            do q = 0, config%g5(p)%integer_zeta
                polynomial_coefficient = coefficients(config%g5(p)%output)* &
                    2.0_real64**(1 - config%g5(p)%integer_zeta)* &
                    binomial_real(config%g5(p)%integer_zeta, q)*config%g5(p)%lambda**q
                do entry = config%g5_moment_first(q + 1), config%g5_moment_last(q + 1)
                    if (species1 == species2) then
                        moment_adjoints(entry, cutoff_index, exponential_index, species1) = &
                            moment_adjoints(entry, cutoff_index, exponential_index, species1) + &
                            polynomial_coefficient*config%g5_moment_multinomial(entry)* &
                            moments(entry, cutoff_index, exponential_index, species1)
                    else
                        moment_adjoints(entry, cutoff_index, exponential_index, species1) = &
                            moment_adjoints(entry, cutoff_index, exponential_index, species1) + &
                            polynomial_coefficient*config%g5_moment_multinomial(entry)* &
                            moments(entry, cutoff_index, exponential_index, species2)
                        moment_adjoints(entry, cutoff_index, exponential_index, species2) = &
                            moment_adjoints(entry, cutoff_index, exponential_index, species2) + &
                            polynomial_coefficient*config%g5_moment_multinomial(entry)* &
                            moments(entry, cutoff_index, exponential_index, species1)
                    end if
                end do
                if (species1 == species2) self_adjoints(cutoff_index, exponential_index, species1) = &
                    self_adjoints(cutoff_index, exponential_index, species1) - 0.5_real64*polynomial_coefficient
            end do
        end do

        do j = 1, nneighbors
            if (distances(j) <= EPS_DISTANCE .or. distances(j) > config%maximum_angular_cutoff) cycle
            species = neighbor_species(j)
            ux = unit_vectors(1, j)
            uy = unit_vectors(2, j)
            uz = unit_vectors(3, j)
            do entry = 1, config%number_of_g5_moments
                x_power = config%g5_moment_x_power(entry)
                y_power = config%g5_moment_y_power(entry)
                z_power = config%g5_moment_z_power(entry)
                monomials(entry) = moment_monomial(unit_vectors(:, j), x_power, y_power, z_power)
                gradient_x = 0.0_real64
                gradient_y = 0.0_real64
                gradient_z = 0.0_real64
                if (x_power > 0) gradient_x = real(x_power, real64)*ux**(x_power - 1)*uy**y_power*uz**z_power
                if (y_power > 0) gradient_y = real(y_power, real64)*ux**x_power*uy**(y_power - 1)*uz**z_power
                if (z_power > 0) gradient_z = real(z_power, real64)*ux**x_power*uy**y_power*uz**(z_power - 1)
                gradient_dot_u = ux*gradient_x + uy*gradient_y + uz*gradient_z
                monomial_gradient_x(entry) = (gradient_x - gradient_dot_u*ux)/distances(j)
                monomial_gradient_y(entry) = (gradient_y - gradient_dot_u*uy)/distances(j)
                monomial_gradient_z(entry) = (gradient_z - gradient_dot_u*uz)/distances(j)
            end do
            derivative_x = 0.0_real64
            derivative_y = 0.0_real64
            derivative_z = 0.0_real64
            do exponential_index = 1, nangular_exp
                do cutoff_index = 1, ncutoff
                    h = angular_exponentials(j, exponential_index)*cutoffs(j, cutoff_index)
                    if (h == 0.0_real64) cycle
                    radial_derivative = angular_exponentials(j, exponential_index)* &
                        (cutoff_derivatives(j, cutoff_index) - &
                         2.0_real64*config%angular_eta(exponential_index)* &
                         (distances(j) - config%angular_shift(exponential_index))*cutoffs(j, cutoff_index))
                    dhx = radial_derivative*ux
                    dhy = radial_derivative*uy
                    dhz = radial_derivative*uz
                    self_scale = 2.0_real64*self_adjoints(cutoff_index, exponential_index, species)*h
                    derivative_x = derivative_x + self_scale*dhx
                    derivative_y = derivative_y + self_scale*dhy
                    derivative_z = derivative_z + self_scale*dhz
                    do entry = 1, config%number_of_g5_moments
                        dphix = dhx*monomials(entry) + h*monomial_gradient_x(entry)
                        dphiy = dhy*monomials(entry) + h*monomial_gradient_y(entry)
                        dphiz = dhz*monomials(entry) + h*monomial_gradient_z(entry)
                        adjoint = moment_adjoints(entry, cutoff_index, exponential_index, species)
                        derivative_x = derivative_x + adjoint*dphix
                        derivative_y = derivative_y + adjoint*dphiy
                        derivative_z = derivative_z + adjoint*dphiz
                    end do
                end do
            end do
            contracted_neighbors(1, j) = contracted_neighbors(1, j) + derivative_x
            contracted_neighbors(2, j) = contracted_neighbors(2, j) + derivative_y
            contracted_neighbors(3, j) = contracted_neighbors(3, j) + derivative_z
            contracted_center(1) = contracted_center(1) - derivative_x
            contracted_center(2) = contracted_center(2) - derivative_y
            contracted_center(3) = contracted_center(3) - derivative_z
        end do
    end subroutine contract_g5_integer_moments

    pure integer function allocated_size_angular(values) result(n)
        type(angular_parameter), allocatable, intent(in) :: values(:)
        n = 0
        if (allocated(values)) n = size(values)
    end function allocated_size_angular

    subroutine evaluate_core(config, displacements, neighbor_species, values, do_derivatives, &
                             derivative_center, derivative_neighbors)
        type(behler_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: values(:)
        logical, intent(in) :: do_derivatives
        real(real64), intent(inout) :: derivative_center(:, :)
        real(real64), intent(inout) :: derivative_neighbors(:, :, :)
        integer :: j, k, p, group, pair, output, nneighbors, ncutoff, nradial_exp, nangular_exp
        real(real64) :: rj, rk, rjk, rjk_squared, cosine, fcj, fck, fcjk, dfcj, dfck, dfcjk
        real(real64) :: exponential, radial_value, radial_derivative
        real(real64) :: angular, dangular, qj, qk, qjk, dqj, dqk, dqjk, factor
        real(real64) :: uj(3), uk(3), ujk(3), dcj(3), dck(3), gradient_j(3), gradient_k(3)
        real(real64) :: displacement_jk(3)
        logical :: g4_pair_in_range
        real(real64) :: distances(size(neighbor_species)), unit_vectors(3, size(neighbor_species))
        real(real64) :: cutoffs(size(neighbor_species), max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: cutoff_derivatives(size(neighbor_species), max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: radial_exponentials(size(neighbor_species), max(1, allocated_size(config%radial_eta)))
        real(real64) :: angular_exponentials(size(neighbor_species), max(1, allocated_size(config%angular_eta)))
        real(real64) :: pair_cutoffs(max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: pair_cutoff_derivatives(max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: pair_exponentials(max(1, allocated_size(config%angular_eta)))
        real(real64) :: g4_products(max(1, allocated_size(config%angular_eta)), &
                                    max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: g5_products(max(1, allocated_size(config%angular_eta)), &
                                    max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: g4_radial_j(3, max(1, allocated_size(config%angular_eta)), &
                                    max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: g4_radial_k(3, max(1, allocated_size(config%angular_eta)), &
                                    max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: g5_radial_j(3, max(1, allocated_size(config%angular_eta)), &
                                    max(1, allocated_size(config%cutoff_radii)))
        real(real64) :: g5_radial_k(3, max(1, allocated_size(config%angular_eta)), &
                                    max(1, allocated_size(config%cutoff_radii)))
        nneighbors = size(neighbor_species)
        ncutoff = 0
        nradial_exp = 0
        nangular_exp = 0
        if (allocated(config%cutoff_radii)) ncutoff = size(config%cutoff_radii)
        if (allocated(config%radial_eta)) nradial_exp = size(config%radial_eta)
        if (allocated(config%angular_eta)) nangular_exp = size(config%angular_eta)
        values(1:config%number_of_descriptors) = 0.0_real64
        if (do_derivatives) then
            derivative_center(:, 1:config%number_of_descriptors) = 0.0_real64
            derivative_neighbors(:, 1:config%number_of_descriptors, :) = 0.0_real64
        end if

        unit_vectors(:, 1:nneighbors) = 0.0_real64
        do j = 1, nneighbors
            distances(j) = sqrt(dot_product(displacements(:, j), displacements(:, j)))
            if (distances(j) > EPS_DISTANCE) unit_vectors(:, j) = displacements(:, j)/distances(j)
            do group = 1, ncutoff
                if (do_derivatives) then
                    call cutoff_value_derivative(distances(j), config%cutoff_radii(group), &
                        config%cutoff_type, config%cutoff_alpha, cutoffs(j, group), cutoff_derivatives(j, group))
                else
                    cutoffs(j, group) = cutoff_value(distances(j), config%cutoff_radii(group), &
                        config%cutoff_type, config%cutoff_alpha)
                end if
            end do
            do group = 1, nradial_exp
                radial_exponentials(j, group) = exp(-config%radial_eta(group)* &
                    (distances(j) - config%radial_shift(group))**2)
            end do
            do group = 1, nangular_exp
                angular_exponentials(j, group) = exp(-config%angular_eta(group)* &
                    (distances(j) - config%angular_shift(group))**2)
            end do
        end do

        do j = 1, nneighbors
            rj = distances(j)
            if (rj <= EPS_DISTANCE) cycle
            uj = unit_vectors(:, j)
            associate(parameters => config%g1_groups(neighbor_species(j)))
            if (parameters%count > 0) then
                do p = 1, parameters%count
                    if (rj > parameters%rc(p)) cycle
                    output = parameters%output(p)
                    fcj = cutoffs(j, parameters%cutoff_group(p))
                    values(output) = values(output) + fcj
                    if (do_derivatives) then
                        gradient_j = cutoff_derivatives(j, parameters%cutoff_group(p))*uj
                        call add_gradient(output, j, gradient_j)
                    end if
                end do
            end if
            end associate
            associate(parameters => config%g2_groups(neighbor_species(j)))
            if (parameters%count > 0) then
                do p = 1, parameters%count
                    if (rj > parameters%rc(p)) cycle
                    output = parameters%output(p)
                    fcj = cutoffs(j, parameters%cutoff_group(p))
                    exponential = radial_exponentials(j, parameters%exponential_group(p))
                    values(output) = values(output) + fcj*exponential
                    if (do_derivatives) then
                        radial_derivative = exponential*(cutoff_derivatives(j, parameters%cutoff_group(p)) - &
                            2.0_real64*parameters%eta(p)*(rj - parameters%rs(p))*fcj)
                        call add_gradient(output, j, radial_derivative*uj)
                    end if
                end do
            end if
            end associate
            associate(parameters => config%g3_groups(neighbor_species(j)))
            if (parameters%count > 0) then
                do p = 1, parameters%count
                    if (rj > parameters%rc(p)) cycle
                    output = parameters%output(p)
                    fcj = cutoffs(j, parameters%cutoff_group(p))
                    radial_value = cos(parameters%kappa(p)*rj)
                    values(output) = values(output) + fcj*radial_value
                    if (do_derivatives) then
                        radial_derivative = cutoff_derivatives(j, parameters%cutoff_group(p))*radial_value - &
                            parameters%kappa(p)*fcj*sin(parameters%kappa(p)*rj)
                        call add_gradient(output, j, radial_derivative*uj)
                    end if
                end do
            end if
            end associate
            cycle
            if (.not. allocated(config%g4) .and. .not. allocated(config%g5)) cycle
            if (rj > config%maximum_angular_cutoff) cycle
            do k = j + 1, nneighbors
                rk = distances(k)
                if (rk <= EPS_DISTANCE .or. rk > config%maximum_angular_cutoff) cycle
                g4_pair_in_range = .false.
                if (allocated(config%g4)) then
                    displacement_jk = displacements(:, k) - displacements(:, j)
                    rjk_squared = dot_product(displacement_jk, displacement_jk)
                    g4_pair_in_range = rjk_squared > EPS_DISTANCE**2 .and. &
                        rjk_squared < config%maximum_g4_cutoff**2
                    if (.not. g4_pair_in_range .and. .not. allocated(config%g5)) cycle
                end if
                uk = unit_vectors(:, k)
                cosine = max(-1.0_real64, min(1.0_real64, dot_product(uj, uk)))
                if (do_derivatives) then
                    dcj = (uk - cosine*uj)/rj
                    dck = (uj - cosine*uk)/rk
                end if
                pair = unordered_pair_index(config%num_species, neighbor_species(j), neighbor_species(k))
                if (allocated(config%g5)) then
                    do group = 1, nangular_exp
                        do output = 1, ncutoff
                            fcj = cutoffs(j, output); fck = cutoffs(k, output)
                            qj = angular_exponentials(j, group)*fcj
                            qk = angular_exponentials(k, group)*fck
                            g5_products(group, output) = qj*qk
                            if (do_derivatives) then
                                dqj = angular_exponentials(j, group)*(cutoff_derivatives(j, output) - &
                                    2.0_real64*config%angular_eta(group)* &
                                    (rj - config%angular_shift(group))*fcj)
                                dqk = angular_exponentials(k, group)*(cutoff_derivatives(k, output) - &
                                    2.0_real64*config%angular_eta(group)* &
                                    (rk - config%angular_shift(group))*fck)
                                g5_radial_j(:, group, output) = dqj*uj*qk
                                g5_radial_k(:, group, output) = qj*dqk*uk
                            end if
                        end do
                    end do
                end if
                if (g4_pair_in_range) then
                    rjk = sqrt(rjk_squared)
                    if (do_derivatives) ujk = displacement_jk/rjk
                    do group = 1, ncutoff
                        pair_cutoffs(group) = cutoff_value(rjk, config%cutoff_radii(group), &
                            config%cutoff_type, config%cutoff_alpha)
                        if (do_derivatives) pair_cutoff_derivatives(group) = &
                            cutoff_derivative(rjk, config%cutoff_radii(group), &
                            config%cutoff_type, config%cutoff_alpha)
                    end do
                    do group = 1, nangular_exp
                        pair_exponentials(group) = exp(-config%angular_eta(group)* &
                            (rjk - config%angular_shift(group))**2)
                    end do
                    do group = 1, nangular_exp
                        do output = 1, ncutoff
                            fcj = cutoffs(j, output); fck = cutoffs(k, output); fcjk = pair_cutoffs(output)
                            qj = angular_exponentials(j, group)*fcj
                            qk = angular_exponentials(k, group)*fck
                            qjk = pair_exponentials(group)*fcjk
                            g4_products(group, output) = qj*qk*qjk
                            if (do_derivatives) then
                                dqj = angular_exponentials(j, group)*(cutoff_derivatives(j, output) - &
                                    2.0_real64*config%angular_eta(group)* &
                                    (rj - config%angular_shift(group))*fcj)
                                dqk = angular_exponentials(k, group)*(cutoff_derivatives(k, output) - &
                                    2.0_real64*config%angular_eta(group)* &
                                    (rk - config%angular_shift(group))*fck)
                                dqjk = pair_exponentials(group)*(pair_cutoff_derivatives(output) - &
                                    2.0_real64*config%angular_eta(group)* &
                                    (rjk - config%angular_shift(group))*fcjk)
                                g4_radial_j(:, group, output) = dqj*uj*qk*qjk - qj*qk*dqjk*ujk
                                g4_radial_k(:, group, output) = qj*dqk*uk*qjk + qj*qk*dqjk*ujk
                            end if
                        end do
                    end do
                    p = config%g4_pair_first(pair)
                    do while (p > 0)
                        if (rj > config%g4(p)%rc .or. rk > config%g4(p)%rc .or. rjk > config%g4(p)%rc) then
                            p = config%g4(p)%next_in_pair
                            cycle
                        end if
                        call angular_power(cosine, config%g4(p)%lambda, config%g4(p)%zeta, &
                            config%g4(p)%integer_zeta, config%g4(p)%derivative_prefactor, angular, dangular)
                        group = config%g4(p)%exponential_group
                        output = config%g4(p)%output
                        factor = 2.0_real64
                        values(output) = values(output) + factor*angular* &
                            g4_products(group, config%g4(p)%cutoff_group)
                        if (do_derivatives) then
                            gradient_j = factor*(dangular*dcj*g4_products(group, config%g4(p)%cutoff_group) + &
                                angular*g4_radial_j(:, group, config%g4(p)%cutoff_group))
                            gradient_k = factor*(dangular*dck*g4_products(group, config%g4(p)%cutoff_group) + &
                                angular*g4_radial_k(:, group, config%g4(p)%cutoff_group))
                            call add_pair_gradients(output, j, k, gradient_j, gradient_k)
                        end if
                        p = config%g4(p)%next_in_pair
                    end do
                end if
                if (allocated(config%g5)) then
                    p = config%g5_pair_first(pair)
                    do while (p > 0)
                        if (rj > config%g5(p)%rc .or. rk > config%g5(p)%rc) then
                            p = config%g5(p)%next_in_pair
                            cycle
                        end if
                        call angular_power(cosine, config%g5(p)%lambda, config%g5(p)%zeta, &
                            config%g5(p)%integer_zeta, config%g5(p)%derivative_prefactor, angular, dangular)
                        group = config%g5(p)%exponential_group
                        output = config%g5(p)%output
                        factor = 2.0_real64
                        values(output) = values(output) + factor*angular* &
                            g5_products(group, config%g5(p)%cutoff_group)
                        if (do_derivatives) then
                            gradient_j = factor*(dangular*dcj*g5_products(group, config%g5(p)%cutoff_group) + &
                                angular*g5_radial_j(:, group, config%g5(p)%cutoff_group))
                            gradient_k = factor*(dangular*dck*g5_products(group, config%g5(p)%cutoff_group) + &
                                angular*g5_radial_k(:, group, config%g5(p)%cutoff_group))
                            call add_pair_gradients(output, j, k, gradient_j, gradient_k)
                        end if
                        p = config%g5(p)%next_in_pair
                    end do
                end if
            end do
        end do
        if (config%has_zero_shift_g4) call evaluate_g4_zero_shift()
        if (config%has_shifted_g4) call evaluate_g4_only()
        if (allocated(config%g5)) call evaluate_g5_only()

    contains
        subroutine evaluate_g4_zero_shift()
            integer :: jj, kk, pp, pair_index, cutoff_group, cutoff_index, angular_index, exponential_index, output, active
            real(real64) :: distance_j, distance_k, distance_jk, distance_jk_squared, radius
            real(real64) :: inverse_j, inverse_k, inverse_jk, inverse_product, cosine_jk, squared_sum
            real(real64) :: cutoff_j, cutoff_k, cutoff_jk, dcutoff_jk, cutoff_product
            real(real64) :: radial_j, radial_k, radial_jk, eta_product
            real(real64) :: energy_coefficient, angular_coefficient, geometry_j, geometry_k, p1, p2, p3
            real(real64) :: pair_displacement(3), edge_j(3), edge_k(3), edge_jk(3)
            real(real64) :: exponentials(max(1, nangular_exp))
            real(real64) :: angular_values(max(1, size(config%g4)))
            real(real64) :: angular_derivatives(max(1, size(config%g4)))

            ! Model setup partitions by species pair, cutoff and exact Rs=0.
            ! Value-only and derivative calls share the same value arithmetic.
            do jj = 1, nneighbors
                distance_j = distances(jj)
                if (distance_j <= EPS_DISTANCE .or. distance_j >= config%maximum_g4_cutoff) cycle
                inverse_j = 1.0_real64/distance_j
                do kk = jj + 1, nneighbors
                    distance_k = distances(kk)
                    if (distance_k <= EPS_DISTANCE .or. distance_k >= config%maximum_g4_cutoff) cycle
                    pair_index = unordered_pair_index(config%num_species, neighbor_species(jj), neighbor_species(kk))
                    if (.not. allocated(config%g4_zero_shift_groups(pair_index)%cutoffs)) cycle
                    pair_displacement = displacements(:, kk) - displacements(:, jj)
                    distance_jk_squared = dot_product(pair_displacement, pair_displacement)
                    if (distance_jk_squared <= EPS_DISTANCE**2 .or. &
                        distance_jk_squared >= config%maximum_g4_cutoff**2) cycle
                    distance_jk = sqrt(distance_jk_squared)
                    inverse_k = 1.0_real64/distance_k
                    inverse_jk = 1.0_real64/distance_jk
                    inverse_product = inverse_j*inverse_k
                    cosine_jk = max(-1.0_real64, min(1.0_real64, &
                        dot_product(displacements(:, jj), displacements(:, kk))*inverse_product))
                    geometry_j = inverse_product - cosine_jk*inverse_j*inverse_j
                    geometry_k = inverse_product - cosine_jk*inverse_k*inverse_k
                    squared_sum = distance_j*distance_j + distance_k*distance_k + distance_jk_squared

                    do cutoff_group = 1, size(config%g4_zero_shift_groups(pair_index)%cutoffs)
                        associate(group_data => config%g4_zero_shift_groups(pair_index)%cutoffs(cutoff_group))
                        cutoff_index = group_data%cutoff_group
                        radius = config%cutoff_radii(cutoff_index)
                        if (max(distance_j, distance_k, distance_jk) >= radius) cycle
                        cutoff_j = cutoffs(jj, cutoff_index)
                        cutoff_k = cutoffs(kk, cutoff_index)
                        call cutoff_value_derivative(distance_jk, radius, config%cutoff_type, &
                            config%cutoff_alpha, cutoff_jk, dcutoff_jk)
                        cutoff_product = cutoff_j*cutoff_k*cutoff_jk
                        if (do_derivatives) then
                            radial_j = cutoff_derivatives(jj, cutoff_index)*cutoff_k*cutoff_jk*inverse_j
                            radial_k = cutoff_j*cutoff_derivatives(kk, cutoff_index)*cutoff_jk*inverse_k
                            radial_jk = cutoff_j*cutoff_k*dcutoff_jk*inverse_jk
                        end if
                        associate(parameters => group_data%parameters)
                        do active = 1, size(parameters%active_exponentials)
                            exponential_index = parameters%active_exponentials(active)
                            exponentials(exponential_index) = exp(-config%angular_eta(exponential_index)*squared_sum)
                        end do
                        do angular_index = 1, parameters%angular_count
                            if (do_derivatives) then
                                call angular_power(cosine_jk, parameters%angular_lambda(angular_index), &
                                    parameters%angular_zeta(angular_index), parameters%angular_integer_zeta(angular_index), &
                                    parameters%angular_derivative_prefactor(angular_index), &
                                    angular_values(angular_index), angular_derivatives(angular_index))
                            else
                                angular_values(angular_index) = angular_value(cosine_jk, &
                                    parameters%angular_lambda(angular_index), parameters%angular_zeta(angular_index), &
                                    parameters%angular_integer_zeta(angular_index))
                            end if
                        end do
                        do pp = 1, parameters%count
                            output = parameters%output(pp)
                            angular_index = parameters%angular_group(pp)
                            exponential_index = parameters%exponential_group(pp)
                            energy_coefficient = 2.0_real64*angular_values(angular_index)*exponentials(exponential_index)
                            values(output) = values(output) + energy_coefficient*cutoff_product
                            if (.not. do_derivatives) cycle
                            angular_coefficient = 2.0_real64*angular_derivatives(angular_index)* &
                                exponentials(exponential_index)*cutoff_product
                            eta_product = 2.0_real64*config%angular_eta(exponential_index)*cutoff_product
                            p1 = angular_coefficient*geometry_j + energy_coefficient*(radial_j - eta_product)
                            p2 = angular_coefficient*geometry_k + energy_coefficient*(radial_k - eta_product)
                            p3 = angular_coefficient*inverse_product - energy_coefficient*(radial_jk - eta_product)
                            ! Three shared edge products supply both neighbor
                            ! gradients and the center gradient, without dividing
                            ! by cutoff values or the angular factor.
                            edge_j = p1*displacements(:, jj)
                            edge_k = p2*displacements(:, kk)
                            edge_jk = p3*pair_displacement
                            derivative_neighbors(:, output, jj) = derivative_neighbors(:, output, jj) + edge_j + edge_jk
                            derivative_neighbors(:, output, kk) = derivative_neighbors(:, output, kk) + edge_k - edge_jk
                            derivative_center(:, output) = derivative_center(:, output) - edge_j - edge_k
                        end do
                        end associate
                        end associate
                    end do
                end do
            end do
        end subroutine evaluate_g4_zero_shift

        subroutine evaluate_g4_only()
            integer :: jj, kk, pp, cutoff_index, exponential_index, descriptor, active_index, product_index
            real(real64) :: distance_j, distance_k, distance_jk, distance_jk_squared
            real(real64) :: cosine_jk, cutoff_j, cutoff_k, cutoff_jk
            real(real64) :: exp_j, exp_k, exp_jk, product, angular_term, angular_derivative
            real(real64) :: dcutoff_j, dcutoff_k, dcutoff_jk
            real(real64) :: pair_displacement(3)
            real(real64) :: vector_j(3), vector_k(3), vector_jk(3)
            real(real64) :: angular_coefficient, radial_coefficient, exponential_product
            real(real64) :: inverse_j, inverse_k, cross_j, cross_k, along_j, along_k, along_jk
            ! One scalar radial derivative per triangle edge and active product.
            ! No division by cutoff values: they may vanish at the boundary.
            real(real64) :: products(max(1, size(config%g4)))
            real(real64) :: radial_j(max(1, size(config%g4)))
            real(real64) :: radial_k(max(1, size(config%g4)))
            real(real64) :: radial_jk(max(1, size(config%g4)))
            real(real64) :: dj1, dj2, dj3, dk1, dk2, dk3
            real(real64) :: angular_values(max(1, size(config%g4)))
            real(real64) :: angular_derivatives(max(1, size(config%g4)))

            do jj = 1, nneighbors
                distance_j = distances(jj)
                if (distance_j <= EPS_DISTANCE .or. distance_j > config%maximum_g4_cutoff) cycle
                vector_j = unit_vectors(:, jj)
                do kk = jj + 1, nneighbors
                    distance_k = distances(kk)
                    if (distance_k <= EPS_DISTANCE .or. distance_k > config%maximum_g4_cutoff) cycle
                    pair = unordered_pair_index(config%num_species, neighbor_species(jj), neighbor_species(kk))
                    if (config%g4_groups(pair)%count == 0) cycle
                    pair_displacement = displacements(:, kk) - displacements(:, jj)
                    distance_jk_squared = dot_product(pair_displacement, pair_displacement)
                    if (distance_jk_squared <= EPS_DISTANCE**2 .or. &
                        distance_jk_squared >= config%maximum_g4_cutoff**2) cycle

                    vector_k = unit_vectors(:, kk)
                    cosine_jk = max(-1.0_real64, min(1.0_real64, dot_product(vector_j, vector_k)))
                    associate(parameters => config%g4_groups(pair))

                    distance_jk = sqrt(distance_jk_squared)
                    if (do_derivatives) then
                        vector_jk = pair_displacement/distance_jk
                        inverse_j = 1.0_real64/distance_j
                        inverse_k = 1.0_real64/distance_k
                    end if

                    do active_index = 1, size(parameters%active_cutoffs)
                        group = parameters%active_cutoffs(active_index)
                        if (do_derivatives) then
                            call cutoff_value_derivative(distance_jk, config%cutoff_radii(group), &
                                config%cutoff_type, config%cutoff_alpha, pair_cutoffs(group), pair_cutoff_derivatives(group))
                        else
                            pair_cutoffs(group) = cutoff_value(distance_jk, config%cutoff_radii(group), &
                                config%cutoff_type, config%cutoff_alpha)
                        end if
                    end do
                    do active_index = 1, size(parameters%active_exponentials)
                        group = parameters%active_exponentials(active_index)
                        pair_exponentials(group) = exp(-config%angular_eta(group)* &
                            (distance_jk - config%angular_shift(group))**2)
                    end do
                    do product_index = 1, size(parameters%product_cutoff)
                        exponential_index = parameters%product_exponential(product_index)
                        cutoff_index = parameters%product_cutoff(product_index)
                        cutoff_j = cutoffs(jj, cutoff_index)
                        cutoff_k = cutoffs(kk, cutoff_index)
                        cutoff_jk = pair_cutoffs(cutoff_index)
                        exp_j = angular_exponentials(jj, exponential_index)
                        exp_k = angular_exponentials(kk, exponential_index)
                        exp_jk = pair_exponentials(exponential_index)
                        products(product_index) = &
                            exp_j*cutoff_j*exp_k*cutoff_k*exp_jk*cutoff_jk
                        if (do_derivatives) then
                            dcutoff_j = cutoff_derivatives(jj, cutoff_index) - &
                                2.0_real64*config%angular_eta(exponential_index)* &
                                (distance_j - config%angular_shift(exponential_index))*cutoff_j
                            dcutoff_k = cutoff_derivatives(kk, cutoff_index) - &
                                2.0_real64*config%angular_eta(exponential_index)* &
                                (distance_k - config%angular_shift(exponential_index))*cutoff_k
                            dcutoff_jk = pair_cutoff_derivatives(cutoff_index) - &
                                2.0_real64*config%angular_eta(exponential_index)* &
                                (distance_jk - config%angular_shift(exponential_index))*cutoff_jk
                            exponential_product = exp_j*exp_k*exp_jk
                            radial_j(product_index) = exponential_product*dcutoff_j*cutoff_k*cutoff_jk
                            radial_k(product_index) = exponential_product*cutoff_j*dcutoff_k*cutoff_jk
                            radial_jk(product_index) = exponential_product*cutoff_j*cutoff_k*dcutoff_jk
                        end if
                    end do

                    do group = 1, parameters%angular_count
                        if (do_derivatives) then
                            call angular_power(cosine_jk, parameters%angular_lambda(group), &
                                parameters%angular_zeta(group), parameters%angular_integer_zeta(group), &
                                parameters%angular_derivative_prefactor(group), angular_values(group), &
                                angular_derivatives(group))
                        else
                            angular_values(group) = angular_value(cosine_jk, parameters%angular_lambda(group), &
                                parameters%angular_zeta(group), parameters%angular_integer_zeta(group))
                        end if
                    end do

                    if (do_derivatives .and. parameters%outputs_contiguous) then
                    do pp = 1, parameters%count
                        if (distance_j > parameters%rc(pp) .or. distance_k > parameters%rc(pp) .or. &
                            distance_jk > parameters%rc(pp)) cycle
                        product_index = parameters%product_group(pp)
                        descriptor = parameters%first_output + pp - 1
                        product = products(product_index)
                        group = parameters%angular_group(pp)
                        angular_term = angular_values(group)
                        angular_derivative = angular_derivatives(group)
                        values(descriptor) = values(descriptor) + 2.0_real64*angular_term*product
                        angular_coefficient = 2.0_real64*angular_derivative*product
                        radial_coefficient = 2.0_real64*angular_term
                        cross_j = angular_coefficient*inverse_j
                        cross_k = angular_coefficient*inverse_k
                        along_j = radial_coefficient*radial_j(product_index) - cross_j*cosine_jk
                        along_k = radial_coefficient*radial_k(product_index) - cross_k*cosine_jk
                        along_jk = radial_coefficient*radial_jk(product_index)
                        dj1 = along_j*vector_j(1) + cross_j*vector_k(1) - along_jk*vector_jk(1)
                        dj2 = along_j*vector_j(2) + cross_j*vector_k(2) - along_jk*vector_jk(2)
                        dj3 = along_j*vector_j(3) + cross_j*vector_k(3) - along_jk*vector_jk(3)
                        dk1 = along_k*vector_k(1) + cross_k*vector_j(1) + along_jk*vector_jk(1)
                        dk2 = along_k*vector_k(2) + cross_k*vector_j(2) + along_jk*vector_jk(2)
                        dk3 = along_k*vector_k(3) + cross_k*vector_j(3) + along_jk*vector_jk(3)
                        derivative_neighbors(1, descriptor, jj) = derivative_neighbors(1, descriptor, jj) + dj1
                        derivative_neighbors(2, descriptor, jj) = derivative_neighbors(2, descriptor, jj) + dj2
                        derivative_neighbors(3, descriptor, jj) = derivative_neighbors(3, descriptor, jj) + dj3
                        derivative_neighbors(1, descriptor, kk) = derivative_neighbors(1, descriptor, kk) + dk1
                        derivative_neighbors(2, descriptor, kk) = derivative_neighbors(2, descriptor, kk) + dk2
                        derivative_neighbors(3, descriptor, kk) = derivative_neighbors(3, descriptor, kk) + dk3
                        derivative_center(1, descriptor) = derivative_center(1, descriptor) - dj1 - dk1
                        derivative_center(2, descriptor) = derivative_center(2, descriptor) - dj2 - dk2
                        derivative_center(3, descriptor) = derivative_center(3, descriptor) - dj3 - dk3
                    end do
                    else
                    do pp = 1, parameters%count
                        if (distance_j > parameters%rc(pp) .or. distance_k > parameters%rc(pp) .or. &
                            distance_jk > parameters%rc(pp)) cycle
                        product_index = parameters%product_group(pp)
                        descriptor = parameters%output(pp)
                        product = products(product_index)
                        group = parameters%angular_group(pp)
                        angular_term = angular_values(group)
                        values(descriptor) = values(descriptor) + 2.0_real64*angular_term*product
                        if (do_derivatives) then
                            angular_derivative = angular_derivatives(group)
                            angular_coefficient = 2.0_real64*angular_derivative*product
                            radial_coefficient = 2.0_real64*angular_term
                            cross_j = angular_coefficient*inverse_j
                            cross_k = angular_coefficient*inverse_k
                            along_j = radial_coefficient*radial_j(product_index) - cross_j*cosine_jk
                            along_k = radial_coefficient*radial_k(product_index) - cross_k*cosine_jk
                            along_jk = radial_coefficient*radial_jk(product_index)
                            dj1 = along_j*vector_j(1) + cross_j*vector_k(1) - along_jk*vector_jk(1)
                            dj2 = along_j*vector_j(2) + cross_j*vector_k(2) - along_jk*vector_jk(2)
                            dj3 = along_j*vector_j(3) + cross_j*vector_k(3) - along_jk*vector_jk(3)
                            dk1 = along_k*vector_k(1) + cross_k*vector_j(1) + along_jk*vector_jk(1)
                            dk2 = along_k*vector_k(2) + cross_k*vector_j(2) + along_jk*vector_jk(2)
                            dk3 = along_k*vector_k(3) + cross_k*vector_j(3) + along_jk*vector_jk(3)
                            derivative_neighbors(1, descriptor, jj) = derivative_neighbors(1, descriptor, jj) + dj1
                            derivative_neighbors(2, descriptor, jj) = derivative_neighbors(2, descriptor, jj) + dj2
                            derivative_neighbors(3, descriptor, jj) = derivative_neighbors(3, descriptor, jj) + dj3
                            derivative_neighbors(1, descriptor, kk) = derivative_neighbors(1, descriptor, kk) + dk1
                            derivative_neighbors(2, descriptor, kk) = derivative_neighbors(2, descriptor, kk) + dk2
                            derivative_neighbors(3, descriptor, kk) = derivative_neighbors(3, descriptor, kk) + dk3
                            derivative_center(1, descriptor) = derivative_center(1, descriptor) - dj1 - dk1
                            derivative_center(2, descriptor) = derivative_center(2, descriptor) - dj2 - dk2
                            derivative_center(3, descriptor) = derivative_center(3, descriptor) - dj3 - dk3
                        end if
                    end do
                    end if
                    end associate
                end do
            end do
        end subroutine evaluate_g4_only

        subroutine evaluate_g5_only()
            integer :: jj, kk, pp, cutoff_index, exponential_index, descriptor
            real(real64) :: distance_j, distance_k, cosine_jk, cutoff_j, cutoff_k
            real(real64) :: exp_j, exp_k, product, angular_term, angular_derivative
            real(real64) :: dcutoff_j, dcutoff_k
            real(real64) :: vector_j(3), vector_k(3), dcosine_j(3), dcosine_k(3)
            real(real64) :: derivative_j(3), derivative_k(3)
            real(real64) :: angular_values(max(1, size(config%g5)))
            real(real64) :: angular_derivatives(max(1, size(config%g5)))
            logical :: use_integer_moments

            use_integer_moments = .not. do_derivatives .and. use_g5_moments(config, distances, .false.)
            if (use_integer_moments) call evaluate_g5_integer_moments()

            do jj = 1, nneighbors
                distance_j = distances(jj)
                if (distance_j <= EPS_DISTANCE .or. distance_j > config%maximum_angular_cutoff) cycle
                vector_j = unit_vectors(:, jj)
                do kk = jj + 1, nneighbors
                    distance_k = distances(kk)
                    if (distance_k <= EPS_DISTANCE .or. distance_k > config%maximum_angular_cutoff) cycle
                    pair = unordered_pair_index(config%num_species, neighbor_species(jj), neighbor_species(kk))
                    associate(parameters => config%g5_groups(pair))
                    if (parameters%count == 0) cycle

                    vector_k = unit_vectors(:, kk)
                    cosine_jk = max(-1.0_real64, min(1.0_real64, dot_product(vector_j, vector_k)))
                    if (do_derivatives) then
                        dcosine_j = (vector_k - cosine_jk*vector_j)/distance_j
                        dcosine_k = (vector_j - cosine_jk*vector_k)/distance_k
                    end if

                    do exponential_index = 1, nangular_exp
                        do cutoff_index = 1, ncutoff
                            cutoff_j = cutoffs(jj, cutoff_index)
                            cutoff_k = cutoffs(kk, cutoff_index)
                            exp_j = angular_exponentials(jj, exponential_index)
                            exp_k = angular_exponentials(kk, exponential_index)
                            g5_products(exponential_index, cutoff_index) = exp_j*cutoff_j*exp_k*cutoff_k
                            if (do_derivatives) then
                                dcutoff_j = cutoff_derivatives(jj, cutoff_index) - &
                                    2.0_real64*config%angular_eta(exponential_index)* &
                                    (distance_j - config%angular_shift(exponential_index))*cutoff_j
                                dcutoff_k = cutoff_derivatives(kk, cutoff_index) - &
                                    2.0_real64*config%angular_eta(exponential_index)* &
                                    (distance_k - config%angular_shift(exponential_index))*cutoff_k
                                g5_radial_j(:, exponential_index, cutoff_index) = &
                                    exp_j*dcutoff_j*vector_j*exp_k*cutoff_k
                                g5_radial_k(:, exponential_index, cutoff_index) = &
                                    exp_j*cutoff_j*exp_k*dcutoff_k*vector_k
                            end if
                        end do
                    end do

                    do group = 1, parameters%angular_count
                        if (use_integer_moments .and. parameters%angular_integer_zeta(group) > 0 .and. &
                            parameters%angular_integer_zeta(group) <= MAX_G5_MOMENT_ORDER) cycle
                        if (do_derivatives) then
                            call angular_power(cosine_jk, parameters%angular_lambda(group), &
                                parameters%angular_zeta(group), parameters%angular_integer_zeta(group), &
                                parameters%angular_derivative_prefactor(group), angular_values(group), &
                                angular_derivatives(group))
                        else
                            angular_values(group) = angular_value(cosine_jk, parameters%angular_lambda(group), &
                                parameters%angular_zeta(group), parameters%angular_integer_zeta(group))
                        end if
                    end do

                    do pp = 1, parameters%count
                        group = parameters%angular_group(pp)
                        if (use_integer_moments .and. parameters%angular_integer_zeta(group) > 0 .and. &
                            parameters%angular_integer_zeta(group) <= MAX_G5_MOMENT_ORDER) cycle
                        if (distance_j > parameters%rc(pp) .or. distance_k > parameters%rc(pp)) cycle
                        cutoff_index = parameters%cutoff_group(pp)
                        exponential_index = parameters%exponential_group(pp)
                        descriptor = parameters%output(pp)
                        product = g5_products(exponential_index, cutoff_index)
                        angular_term = angular_values(group)
                        values(descriptor) = values(descriptor) + 2.0_real64*angular_term*product
                        if (do_derivatives) then
                            angular_derivative = angular_derivatives(group)
                            derivative_j = 2.0_real64*(angular_derivative*dcosine_j*product + &
                                angular_term*g5_radial_j(:, exponential_index, cutoff_index))
                            derivative_k = 2.0_real64*(angular_derivative*dcosine_k*product + &
                                angular_term*g5_radial_k(:, exponential_index, cutoff_index))
                            derivative_neighbors(:, descriptor, jj) = &
                                derivative_neighbors(:, descriptor, jj) + derivative_j
                            derivative_neighbors(:, descriptor, kk) = &
                                derivative_neighbors(:, descriptor, kk) + derivative_k
                            derivative_center(:, descriptor) = derivative_center(:, descriptor) - &
                                derivative_j - derivative_k
                        end if
                    end do
                    end associate
                end do
            end do
        end subroutine evaluate_g5_only

        subroutine evaluate_g5_integer_moments()
            integer :: jj, pp, entry, q, species1, species2, cutoff_index, exponential_index
            real(real64) :: h, pair_moment, polynomial_coefficient
            real(real64) :: moments(max(1, config%number_of_g5_moments), max(1, ncutoff), &
                                    max(1, nangular_exp), config%num_species)
            real(real64) :: self_terms(max(1, ncutoff), max(1, nangular_exp), config%num_species)

            moments = 0.0_real64
            self_terms = 0.0_real64
            do jj = 1, nneighbors
                if (distances(jj) <= EPS_DISTANCE .or. distances(jj) > config%maximum_angular_cutoff) cycle
                do exponential_index = 1, nangular_exp
                    do cutoff_index = 1, ncutoff
                        h = angular_exponentials(jj, exponential_index)*cutoffs(jj, cutoff_index)
                        if (h == 0.0_real64) cycle
                        self_terms(cutoff_index, exponential_index, neighbor_species(jj)) = &
                            self_terms(cutoff_index, exponential_index, neighbor_species(jj)) + h*h
                        do entry = 1, config%number_of_g5_moments
                            moments(entry, cutoff_index, exponential_index, neighbor_species(jj)) = &
                                moments(entry, cutoff_index, exponential_index, neighbor_species(jj)) + h* &
                                moment_monomial(unit_vectors(:, jj), config%g5_moment_x_power(entry), &
                                                config%g5_moment_y_power(entry), &
                                                config%g5_moment_z_power(entry))
                        end do
                    end do
                end do
            end do

            do pp = 1, size(config%g5)
                if (config%g5(pp)%integer_zeta <= 0 .or. &
                    config%g5(pp)%integer_zeta > MAX_G5_MOMENT_ORDER) cycle
                species1 = config%g5(pp)%species1
                species2 = config%g5(pp)%species2
                cutoff_index = config%g5(pp)%cutoff_group
                exponential_index = config%g5(pp)%exponential_group
                do q = 0, config%g5(pp)%integer_zeta
                    pair_moment = 0.0_real64
                    do entry = config%g5_moment_first(q + 1), config%g5_moment_last(q + 1)
                        pair_moment = pair_moment + config%g5_moment_multinomial(entry)* &
                            moments(entry, cutoff_index, exponential_index, species1)* &
                            moments(entry, cutoff_index, exponential_index, species2)
                    end do
                    if (species1 == species2) pair_moment = 0.5_real64*(pair_moment - &
                        self_terms(cutoff_index, exponential_index, species1))
                    polynomial_coefficient = 2.0_real64**(1 - config%g5(pp)%integer_zeta)* &
                        binomial_real(config%g5(pp)%integer_zeta, q)*config%g5(pp)%lambda**q
                    values(config%g5(pp)%output) = values(config%g5(pp)%output) + &
                        polynomial_coefficient*pair_moment
                end do
            end do
        end subroutine evaluate_g5_integer_moments

        subroutine add_gradient(coefficient, neighbor, gradient)
            integer, intent(in) :: coefficient, neighbor
            real(real64), intent(in) :: gradient(3)
            derivative_neighbors(:, coefficient, neighbor) = &
                derivative_neighbors(:, coefficient, neighbor) + gradient
            derivative_center(:, coefficient) = derivative_center(:, coefficient) - gradient
        end subroutine
        subroutine add_pair_gradients(coefficient, neighbor1, neighbor2, gradient1, gradient2)
            integer, intent(in) :: coefficient, neighbor1, neighbor2
            real(real64), intent(in) :: gradient1(3), gradient2(3)
            derivative_neighbors(:, coefficient, neighbor1) = &
                derivative_neighbors(:, coefficient, neighbor1) + gradient1
            derivative_neighbors(:, coefficient, neighbor2) = &
                derivative_neighbors(:, coefficient, neighbor2) + gradient2
            derivative_center(:, coefficient) = derivative_center(:, coefficient) - gradient1 - gradient2
        end subroutine
    end subroutine evaluate_core

end module accelnet_behler
