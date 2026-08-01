module accelnet_behler
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: cutoff_value, cutoff_derivative
    implicit none
    private

    real(real64), parameter :: EPS_DISTANCE = 1.0e-12_real64

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
        real(real64) :: zeta = 0.0_real64, eta = 0.0_real64
        real(real64) :: derivative_prefactor = 0.0_real64
    end type
    type :: g4_parameter_group
        integer :: count = 0, angular_count = 0
        integer :: first_output = 0
        logical :: outputs_contiguous = .true.
        integer, allocatable :: output(:), cutoff_group(:), exponential_group(:), angular_group(:)
        real(real64), allocatable :: rc(:)
        integer, allocatable :: angular_integer_zeta(:)
        real(real64), allocatable :: angular_lambda(:), angular_zeta(:), angular_derivative_prefactor(:)
    end type
    type :: g1_parameter_group
        integer :: count = 0
        integer, allocatable :: output(:), cutoff_group(:)
        real(real64), allocatable :: rc(:)
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
        integer :: number_of_descriptors = 0
        real(real64) :: maximum_cutoff = 0.0_real64
        real(real64) :: maximum_angular_cutoff = 0.0_real64
        real(real64) :: maximum_g4_cutoff = 0.0_real64
        real(real64), allocatable :: cutoff_radii(:)
        real(real64), allocatable :: radial_eta(:), radial_shift(:)
        real(real64), allocatable :: angular_eta(:)
        integer, allocatable :: g4_pair_first(:), g4_pair_last(:)
        integer, allocatable :: g5_pair_first(:), g5_pair_last(:)
        type(g4_parameter_group), allocatable :: g4_groups(:)
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

contains

    subroutine initialize_behler_config(config, num_species)
        type(behler_config), intent(out) :: config
        integer, intent(in) :: num_species
        if (num_species < 1) error stop "Behler num_species must be positive"
        config%num_species = num_species
        allocate(config%g4_groups(num_species*(num_species + 1)/2))
        allocate(config%g5_groups(num_species*(num_species + 1)/2))
        allocate(config%g1_groups(num_species), config%g2_groups(num_species), config%g3_groups(num_species))
    end subroutine

    integer function behler_num_descriptors(self) result(n)
        class(behler_config), intent(in) :: self
        n = self%number_of_descriptors
    end function

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

    subroutine add_g4(config, species1, species2, rc, lambda, zeta, eta)
        type(behler_config), intent(inout) :: config
        integer, intent(in) :: species1, species2
        real(real64), intent(in) :: rc, lambda, zeta, eta
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
        parameter%integer_zeta = integer_zeta_kind(zeta)
        parameter%derivative_prefactor = 0.5_real64*zeta*lambda
        call find_or_add_value(config%cutoff_radii, rc, parameter%cutoff_group)
        call find_or_add_value(config%angular_eta, eta, parameter%exponential_group)
        if (allocated(config%g4)) then
            config%g4 = [config%g4, parameter]
        else
            config%g4 = [parameter]
        end if
        call link_angular_pair(config%g4, config%g4_pair_first, config%g4_pair_last, &
                               config%num_species, species1, species2)
        call append_g4_group(config%g4_groups(unordered_pair_index(config%num_species, species1, species2)), &
                             parameter)
        config%maximum_cutoff = max(config%maximum_cutoff, rc)
        config%maximum_angular_cutoff = max(config%maximum_angular_cutoff, rc)
        config%maximum_g4_cutoff = max(config%maximum_g4_cutoff, rc)
    end subroutine

    subroutine add_g5(config, species1, species2, rc, lambda, zeta, eta)
        type(behler_config), intent(inout) :: config
        integer, intent(in) :: species1, species2
        real(real64), intent(in) :: rc, lambda, zeta, eta
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
        parameter%integer_zeta = integer_zeta_kind(zeta)
        parameter%derivative_prefactor = 0.5_real64*zeta*lambda
        call find_or_add_value(config%cutoff_radii, rc, parameter%cutoff_group)
        call find_or_add_value(config%angular_eta, eta, parameter%exponential_group)
        if (allocated(config%g5)) then
            config%g5 = [config%g5, parameter]
        else
            config%g5 = [parameter]
        end if
        call link_angular_pair(config%g5, config%g5_pair_first, config%g5_pair_last, &
                               config%num_species, species1, species2)
        call append_g4_group(config%g5_groups(unordered_pair_index(config%num_species, species1, species2)), &
                             parameter)
        config%maximum_cutoff = max(config%maximum_cutoff, rc)
        config%maximum_angular_cutoff = max(config%maximum_angular_cutoff, rc)
    end subroutine

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

    subroutine append_g4_group(group, parameter)
        type(g4_parameter_group), intent(inout) :: group
        type(angular_parameter), intent(in) :: parameter
        integer :: i, angular_index
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
                cutoffs(j, group) = cutoff_value(distances(j), config%cutoff_radii(group))
                if (do_derivatives) cutoff_derivatives(j, group) = &
                    cutoff_derivative(distances(j), config%cutoff_radii(group))
            end do
            do group = 1, nradial_exp
                radial_exponentials(j, group) = exp(-config%radial_eta(group)* &
                    (distances(j) - config%radial_shift(group))**2)
            end do
            do group = 1, nangular_exp
                angular_exponentials(j, group) = exp(-config%angular_eta(group)*distances(j)**2)
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
                                    2.0_real64*config%angular_eta(group)*rj*fcj)
                                dqk = angular_exponentials(k, group)*(cutoff_derivatives(k, output) - &
                                    2.0_real64*config%angular_eta(group)*rk*fck)
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
                        pair_cutoffs(group) = cutoff_value(rjk, config%cutoff_radii(group))
                        if (do_derivatives) pair_cutoff_derivatives(group) = &
                            cutoff_derivative(rjk, config%cutoff_radii(group))
                    end do
                    do group = 1, nangular_exp
                        pair_exponentials(group) = exp(-config%angular_eta(group)*rjk*rjk)
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
                                    2.0_real64*config%angular_eta(group)*rj*fcj)
                                dqk = angular_exponentials(k, group)*(cutoff_derivatives(k, output) - &
                                    2.0_real64*config%angular_eta(group)*rk*fck)
                                dqjk = pair_exponentials(group)*(pair_cutoff_derivatives(output) - &
                                    2.0_real64*config%angular_eta(group)*rjk*fcjk)
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
        if (allocated(config%g4)) call evaluate_g4_only()
        if (allocated(config%g5)) call evaluate_g5_only()

    contains
        subroutine evaluate_g4_only()
            integer :: jj, kk, pp, cutoff_index, exponential_index, descriptor
            real(real64) :: distance_j, distance_k, distance_jk, distance_jk_squared
            real(real64) :: cosine_jk, cutoff_j, cutoff_k, cutoff_jk
            real(real64) :: exp_j, exp_k, exp_jk, product, angular_term, angular_derivative
            real(real64) :: dcutoff_j, dcutoff_k, dcutoff_jk
            real(real64) :: pair_displacement(3)
            real(real64) :: vector_j(3), vector_k(3), vector_jk(3)
            real(real64) :: dcosine_j(3), dcosine_k(3), derivative_j(3), derivative_k(3)
            real(real64) :: angular_coefficient, radial_coefficient
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
                    pair_displacement = displacements(:, kk) - displacements(:, jj)
                    distance_jk_squared = dot_product(pair_displacement, pair_displacement)
                    if (distance_jk_squared <= EPS_DISTANCE**2 .or. &
                        distance_jk_squared >= config%maximum_g4_cutoff**2) cycle

                    vector_k = unit_vectors(:, kk)
                    cosine_jk = max(-1.0_real64, min(1.0_real64, dot_product(vector_j, vector_k)))
                    pair = unordered_pair_index(config%num_species, neighbor_species(jj), neighbor_species(kk))
                    associate(parameters => config%g4_groups(pair))
                    if (parameters%count == 0) cycle

                    distance_jk = sqrt(distance_jk_squared)
                    if (do_derivatives) then
                        vector_jk = pair_displacement/distance_jk
                        dcosine_j = (vector_k - cosine_jk*vector_j)/distance_j
                        dcosine_k = (vector_j - cosine_jk*vector_k)/distance_k
                    end if

                    do group = 1, ncutoff
                        pair_cutoffs(group) = cutoff_value(distance_jk, config%cutoff_radii(group))
                        if (do_derivatives) pair_cutoff_derivatives(group) = &
                            cutoff_derivative(distance_jk, config%cutoff_radii(group))
                    end do
                    do group = 1, nangular_exp
                        pair_exponentials(group) = exp(-config%angular_eta(group)*distance_jk_squared)
                    end do
                    do exponential_index = 1, nangular_exp
                        do cutoff_index = 1, ncutoff
                            cutoff_j = cutoffs(jj, cutoff_index)
                            cutoff_k = cutoffs(kk, cutoff_index)
                            cutoff_jk = pair_cutoffs(cutoff_index)
                            exp_j = angular_exponentials(jj, exponential_index)
                            exp_k = angular_exponentials(kk, exponential_index)
                            exp_jk = pair_exponentials(exponential_index)
                            g4_products(exponential_index, cutoff_index) = &
                                exp_j*cutoff_j*exp_k*cutoff_k*exp_jk*cutoff_jk
                            if (do_derivatives) then
                                dcutoff_j = cutoff_derivatives(jj, cutoff_index) - &
                                    2.0_real64*config%angular_eta(exponential_index)*distance_j*cutoff_j
                                dcutoff_k = cutoff_derivatives(kk, cutoff_index) - &
                                    2.0_real64*config%angular_eta(exponential_index)*distance_k*cutoff_k
                                dcutoff_jk = pair_cutoff_derivatives(cutoff_index) - &
                                    2.0_real64*config%angular_eta(exponential_index)*distance_jk*cutoff_jk
                                g4_radial_j(:, exponential_index, cutoff_index) = &
                                    exp_j*dcutoff_j*vector_j*exp_k*cutoff_k*exp_jk*cutoff_jk - &
                                    exp_j*cutoff_j*exp_k*cutoff_k*exp_jk*dcutoff_jk*vector_jk
                                g4_radial_k(:, exponential_index, cutoff_index) = &
                                    exp_j*cutoff_j*exp_k*dcutoff_k*vector_k*exp_jk*cutoff_jk + &
                                    exp_j*cutoff_j*exp_k*cutoff_k*exp_jk*dcutoff_jk*vector_jk
                            end if
                        end do
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
                        cutoff_index = parameters%cutoff_group(pp)
                        exponential_index = parameters%exponential_group(pp)
                        descriptor = parameters%first_output + pp - 1
                        product = g4_products(exponential_index, cutoff_index)
                        group = parameters%angular_group(pp)
                        angular_term = angular_values(group)
                        angular_derivative = angular_derivatives(group)
                        values(descriptor) = values(descriptor) + 2.0_real64*angular_term*product
                        angular_coefficient = 2.0_real64*angular_derivative*product
                        radial_coefficient = 2.0_real64*angular_term
                        dj1 = angular_coefficient*dcosine_j(1) + &
                            radial_coefficient*g4_radial_j(1, exponential_index, cutoff_index)
                        dj2 = angular_coefficient*dcosine_j(2) + &
                            radial_coefficient*g4_radial_j(2, exponential_index, cutoff_index)
                        dj3 = angular_coefficient*dcosine_j(3) + &
                            radial_coefficient*g4_radial_j(3, exponential_index, cutoff_index)
                        dk1 = angular_coefficient*dcosine_k(1) + &
                            radial_coefficient*g4_radial_k(1, exponential_index, cutoff_index)
                        dk2 = angular_coefficient*dcosine_k(2) + &
                            radial_coefficient*g4_radial_k(2, exponential_index, cutoff_index)
                        dk3 = angular_coefficient*dcosine_k(3) + &
                            radial_coefficient*g4_radial_k(3, exponential_index, cutoff_index)
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
                        cutoff_index = parameters%cutoff_group(pp)
                        exponential_index = parameters%exponential_group(pp)
                        descriptor = parameters%output(pp)
                        product = g4_products(exponential_index, cutoff_index)
                        group = parameters%angular_group(pp)
                        angular_term = angular_values(group)
                        values(descriptor) = values(descriptor) + 2.0_real64*angular_term*product
                        if (do_derivatives) then
                            angular_derivative = angular_derivatives(group)
                            derivative_j = 2.0_real64*(angular_derivative*dcosine_j*product + &
                                angular_term*g4_radial_j(:, exponential_index, cutoff_index))
                            derivative_k = 2.0_real64*(angular_derivative*dcosine_k*product + &
                                angular_term*g4_radial_k(:, exponential_index, cutoff_index))
                            derivative_neighbors(:, descriptor, jj) = &
                                derivative_neighbors(:, descriptor, jj) + derivative_j
                            derivative_neighbors(:, descriptor, kk) = &
                                derivative_neighbors(:, descriptor, kk) + derivative_k
                            derivative_center(:, descriptor) = derivative_center(:, descriptor) - derivative_j - derivative_k
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
                                    2.0_real64*config%angular_eta(exponential_index)*distance_j*cutoff_j
                                dcutoff_k = cutoff_derivatives(kk, cutoff_index) - &
                                    2.0_real64*config%angular_eta(exponential_index)*distance_k*cutoff_k
                                g5_radial_j(:, exponential_index, cutoff_index) = &
                                    exp_j*dcutoff_j*vector_j*exp_k*cutoff_k
                                g5_radial_k(:, exponential_index, cutoff_index) = &
                                    exp_j*cutoff_j*exp_k*dcutoff_k*vector_k
                            end if
                        end do
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

                    do pp = 1, parameters%count
                        if (distance_j > parameters%rc(pp) .or. distance_k > parameters%rc(pp)) cycle
                        cutoff_index = parameters%cutoff_group(pp)
                        exponential_index = parameters%exponential_group(pp)
                        descriptor = parameters%output(pp)
                        product = g5_products(exponential_index, cutoff_index)
                        group = parameters%angular_group(pp)
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
