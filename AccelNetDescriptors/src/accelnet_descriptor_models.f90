module accelnet_descriptor_models
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: descriptor_config, evaluate_atom, evaluate_atom_with_derivatives, &
                                    contract_atom_derivatives, set_chebyshev_evaluation
    use accelnet_lj, only: lj_config, evaluate_lj_values, evaluate_lj_values_derivatives
    use accelnet_behler, only: behler_config, evaluate_behler_values, evaluate_behler_values_derivatives, &
                              contract_behler_derivatives, behler_supports_direct_contraction, &
                              set_behler_g5_evaluation
    implicit none
    private

    type :: chebyshev_component
        type(descriptor_config) :: config
        integer :: output_offset = 0
    end type chebyshev_component

    type :: lj_component
        type(lj_config) :: config
        integer :: output_offset = 0
    end type lj_component

    type :: behler_component
        type(behler_config) :: config
        integer :: output_offset = 0
    end type behler_component

    type, public :: descriptor_model
        type(chebyshev_component), allocatable :: chebyshev(:)
        type(lj_component), allocatable :: lj(:)
        type(behler_component), allocatable :: behler(:)
        integer :: num_outputs = 0
        real(real64) :: maximum_cutoff = 0.0_real64
    contains
        procedure :: num_descriptors => model_num_descriptors
    end type descriptor_model

    public :: add_chebyshev
    public :: add_lj
    public :: add_behler
    public :: evaluate_model_values
    public :: evaluate_model_values_derivatives
    public :: contract_model_derivatives, model_supports_direct_contraction
    public :: set_model_g5_evaluation
    public :: set_model_chebyshev_evaluation

contains

    subroutine set_model_chebyshev_evaluation(model, mode)
        type(descriptor_model), intent(inout) :: model
        integer, intent(in) :: mode
        integer :: component
        if (.not. allocated(model%chebyshev)) return
        do component = 1, size(model%chebyshev)
            call set_chebyshev_evaluation(model%chebyshev(component)%config, mode)
        end do
    end subroutine set_model_chebyshev_evaluation

    subroutine set_model_g5_evaluation(model, mode)
        type(descriptor_model), intent(inout) :: model
        integer, intent(in) :: mode
        integer :: component
        if (.not. allocated(model%behler)) return
        do component = 1, size(model%behler)
            call set_behler_g5_evaluation(model%behler(component)%config, mode)
        end do
    end subroutine set_model_g5_evaluation

    integer function model_num_descriptors(self) result(n)
        class(descriptor_model), intent(in) :: self
        n = self%num_outputs
    end function model_num_descriptors

    subroutine add_chebyshev(model, config)
        type(descriptor_model), intent(inout) :: model
        type(descriptor_config), intent(in) :: config
        type(chebyshev_component), allocatable :: expanded(:)
        integer :: old_size

        old_size = 0
        if (allocated(model%chebyshev)) old_size = size(model%chebyshev)
        allocate(expanded(old_size + 1))
        if (old_size > 0) expanded(1:old_size) = model%chebyshev
        expanded(old_size + 1)%config = config
        expanded(old_size + 1)%output_offset = model%num_outputs
        call move_alloc(expanded, model%chebyshev)
        model%num_outputs = model%num_outputs + config%num_descriptors()
        model%maximum_cutoff = max(model%maximum_cutoff, config%radial_rc, config%angular_rc)
    end subroutine add_chebyshev

    subroutine add_lj(model, config)
        type(descriptor_model), intent(inout) :: model
        type(lj_config), intent(in) :: config
        type(lj_component), allocatable :: expanded(:)
        integer :: old_size

        old_size = 0
        if (allocated(model%lj)) old_size = size(model%lj)
        allocate(expanded(old_size + 1))
        if (old_size > 0) expanded(1:old_size) = model%lj
        expanded(old_size + 1)%config = config
        expanded(old_size + 1)%output_offset = model%num_outputs
        call move_alloc(expanded, model%lj)
        model%num_outputs = model%num_outputs + config%num_descriptors()
        model%maximum_cutoff = max(model%maximum_cutoff, config%radial_rc)
    end subroutine add_lj

    subroutine add_behler(model, config)
        type(descriptor_model), intent(inout) :: model
        type(behler_config), intent(in) :: config
        type(behler_component), allocatable :: expanded(:)
        integer :: old_size

        old_size = 0
        if (allocated(model%behler)) old_size = size(model%behler)
        allocate(expanded(old_size + 1))
        if (old_size > 0) expanded(1:old_size) = model%behler
        expanded(old_size + 1)%config = config
        expanded(old_size + 1)%output_offset = model%num_outputs
        call move_alloc(expanded, model%behler)
        model%num_outputs = model%num_outputs + config%num_descriptors()
        model%maximum_cutoff = max(model%maximum_cutoff, config%maximum_cutoff)
    end subroutine add_behler

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

end module accelnet_descriptor_models
