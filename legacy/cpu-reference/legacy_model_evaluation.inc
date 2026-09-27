! Retained CPU evaluator from gpu checkpoint 543b180; reference/descriptor-only compatibility.
    subroutine evaluate_model_values(model, displacements, neighbor_species, values)
        type(descriptor_model), intent(in) :: model
        real(real64), intent(in) :: displacements(:, :)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: values(:)
        integer :: component, first, last

        if (size(values) < model%num_outputs) error stop "model values array too small"
        values(1:model%num_outputs) = 0.0_real64
        if (allocated(model%chebyshev)) then
            do component = 1, size(model%chebyshev)
                first = model%chebyshev(component)%output_offset + 1
                last = first + model%chebyshev(component)%config%num_descriptors() - 1
                call evaluate_atom(model%chebyshev(component)%config, displacements, &
                                   neighbor_species, values(first:last))
            end do
        end if
        if (allocated(model%lj)) then
            do component = 1, size(model%lj)
                first = model%lj(component)%output_offset + 1
                last = first + model%lj(component)%config%num_descriptors() - 1
                call evaluate_lj_values(model%lj(component)%config, displacements, &
                                        neighbor_species, values(first:last))
            end do
        end if
        if (allocated(model%behler)) then
            do component = 1, size(model%behler)
                first = model%behler(component)%output_offset + 1
                last = first + model%behler(component)%config%num_descriptors() - 1
                call evaluate_behler_values(model%behler(component)%config, displacements, &
                                            neighbor_species, values(first:last))
            end do
        end if
    end subroutine evaluate_model_values

    subroutine evaluate_model_values_derivatives(model, displacements, neighbor_species, values, &
                                                 derivative_center, derivative_neighbors)
        type(descriptor_model), intent(in) :: model
        real(real64), intent(in) :: displacements(:, :)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: values(:)
        real(real64), intent(out) :: derivative_center(:, :)
        real(real64), intent(out) :: derivative_neighbors(:, :, :)
        integer :: component, first, last

        if (size(values) < model%num_outputs .or. size(derivative_center, 2) < model%num_outputs .or. &
            size(derivative_neighbors, 2) < model%num_outputs) error stop "model output arrays too small"
        ! Each component evaluator initializes its own disjoint output slice.
        ! Clearing the complete model arrays here as well made the force path
        ! write every descriptor Jacobian twice before doing useful work.
        if (allocated(model%chebyshev)) then
            do component = 1, size(model%chebyshev)
                first = model%chebyshev(component)%output_offset + 1
                last = first + model%chebyshev(component)%config%num_descriptors() - 1
                call evaluate_atom_with_derivatives(model%chebyshev(component)%config, displacements, &
                    neighbor_species, values(first:last), derivative_center(:, first:last), &
                    derivative_neighbors(:, first:last, :))
            end do
        end if
        if (allocated(model%lj)) then
            do component = 1, size(model%lj)
                first = model%lj(component)%output_offset + 1
                last = first + model%lj(component)%config%num_descriptors() - 1
                call evaluate_lj_values_derivatives(model%lj(component)%config, displacements, &
                    neighbor_species, values(first:last), derivative_center(:, first:last), &
                    derivative_neighbors(:, first:last, :))
            end do
        end if
        if (allocated(model%behler)) then
            do component = 1, size(model%behler)
                first = model%behler(component)%output_offset + 1
                last = first + model%behler(component)%config%num_descriptors() - 1
                call evaluate_behler_values_derivatives(model%behler(component)%config, displacements, &
                    neighbor_species, values(first:last), derivative_center(:, first:last), &
                    derivative_neighbors(:, first:last, :))
            end do
        end if
    end subroutine evaluate_model_values_derivatives

    pure logical function model_supports_direct_contraction(model) result(supported)
        type(descriptor_model), intent(in) :: model
        integer :: component
        supported = .not. allocated(model%lj) .and. &
                    (allocated(model%chebyshev) .or. allocated(model%behler))
        if (allocated(model%behler)) then
            do component = 1, size(model%behler)
                supported = supported .and. &
                    behler_supports_direct_contraction(model%behler(component)%config)
            end do
        end if
    end function model_supports_direct_contraction

    subroutine contract_model_derivatives(model, displacements, neighbor_species, coefficients, &
                                          contracted_center, contracted_neighbors)
        type(descriptor_model), intent(in) :: model
        real(real64), intent(in) :: displacements(:, :), coefficients(:)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: contracted_center(3), contracted_neighbors(:, :)
        real(real64) :: component_center(3)
        real(real64) :: component_neighbors(3, size(neighbor_species))
        integer :: component, first, last

        if (.not. model_supports_direct_contraction(model)) &
            error stop "direct derivative contraction is not implemented for this model"
        if (size(coefficients) < model%num_outputs) error stop "model coefficient array too small"
        if (size(contracted_neighbors, 1) /= 3 .or. &
            size(contracted_neighbors, 2) < size(neighbor_species)) error stop "model contracted array too small"
        contracted_center = 0.0_real64
        contracted_neighbors(:, 1:size(neighbor_species)) = 0.0_real64
        if (allocated(model%chebyshev)) then
            do component = 1, size(model%chebyshev)
                first = model%chebyshev(component)%output_offset + 1
                last = first + model%chebyshev(component)%config%num_descriptors() - 1
                call contract_atom_derivatives(model%chebyshev(component)%config, displacements, &
                    neighbor_species, coefficients(first:last), component_center, component_neighbors)
                contracted_center = contracted_center + component_center
                contracted_neighbors(:, 1:size(neighbor_species)) = &
                    contracted_neighbors(:, 1:size(neighbor_species)) + component_neighbors
            end do
        end if
        if (allocated(model%behler)) then
            do component = 1, size(model%behler)
                first = model%behler(component)%output_offset + 1
                last = first + model%behler(component)%config%num_descriptors() - 1
                call contract_behler_derivatives(model%behler(component)%config, displacements, &
                    neighbor_species, coefficients(first:last), component_center, component_neighbors)
                contracted_center = contracted_center + component_center
                contracted_neighbors(:, 1:size(neighbor_species)) = &
                    contracted_neighbors(:, 1:size(neighbor_species)) + component_neighbors
            end do
        end if
    end subroutine contract_model_derivatives

