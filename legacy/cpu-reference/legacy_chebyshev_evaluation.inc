! Retained CPU evaluator from gpu checkpoint 543b180; reference/descriptor-only compatibility.
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

