! Retained CPU evaluator from gpu checkpoint 543b180; reference/descriptor-only compatibility.
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
        supported = .not. allocated(config%g4) .and. .not. allocated(config%extended)
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
                cutoffs(j, group) = cutoff_value(distances(j), config%cutoff_radii(group), &
                    config%cutoff_type, config%cutoff_alpha)
                if (do_derivatives) cutoff_derivatives(j, group) = &
                    cutoff_derivative(distances(j), config%cutoff_radii(group), &
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
        if (allocated(config%extended)) then
            do p=1,size(config%extended)
                associate(ep => config%extended(p))
                output=ep%output
                do j=1,nneighbors
                    if (ep%species1 > 0 .and. ep%species2 == 0) then
                        if (neighbor_species(j) /= ep%species1) cycle
                    end if
                    call extended_radial(ep%kind,config%cutoff_type,ep%p,distances(j),qj,dqj)
                    if (qj == 0 .and. dqj == 0) cycle
                    if (ep%kind == 12 .or. ep%kind == 20 .or. ep%kind == 23) then
                        factor=1
                        if (ep%species1 == 0) factor=config%species_weights(neighbor_species(j))
                        values(output)=values(output)+factor*qj
                        if (do_derivatives) call add_gradient(output,j,factor*dqj*unit_vectors(:,j))
                    else
                        do k=j+1,nneighbors
                            if (ep%species1 > 0) then
                                if (.not. species_pair_matches(neighbor_species(j),neighbor_species(k), &
                                    ep%species1,ep%species2)) cycle
                            end if
                            call extended_radial(ep%kind,config%cutoff_type,ep%p,distances(k),qk,dqk)
                            call extended_pair(ep%kind,config%cutoff_type,ep%p,unit_vectors(:,j),unit_vectors(:,k), &
                                distances(j),distances(k),qj,qk,dqj,dqk,radial_value,gradient_j,gradient_k)
                            factor=1
                            if (ep%species1 == 0) factor=config%species_weights(neighbor_species(j))* &
                                config%species_weights(neighbor_species(k))
                            values(output)=values(output)+factor*radial_value
                            if (do_derivatives) call add_pair_gradients(output,j,k,factor*gradient_j,factor*gradient_k)
                        end do
                    end if
                end do
                end associate
            end do
        end if
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
                        pair_cutoffs(group) = cutoff_value(distance_jk, config%cutoff_radii(group), &
                            config%cutoff_type, config%cutoff_alpha)
                        if (do_derivatives) pair_cutoff_derivatives(group) = &
                            cutoff_derivative(distance_jk, config%cutoff_radii(group), &
                            config%cutoff_type, config%cutoff_alpha)
                    end do
                    do group = 1, nangular_exp
                        pair_exponentials(group) = exp(-config%angular_eta(group)* &
                            (distance_jk - config%angular_shift(group))**2)
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
                                    2.0_real64*config%angular_eta(exponential_index)* &
                                    (distance_j - config%angular_shift(exponential_index))*cutoff_j
                                dcutoff_k = cutoff_derivatives(kk, cutoff_index) - &
                                    2.0_real64*config%angular_eta(exponential_index)* &
                                    (distance_k - config%angular_shift(exponential_index))*cutoff_k
                                dcutoff_jk = pair_cutoff_derivatives(cutoff_index) - &
                                    2.0_real64*config%angular_eta(exponential_index)* &
                                    (distance_jk - config%angular_shift(exponential_index))*cutoff_jk
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


