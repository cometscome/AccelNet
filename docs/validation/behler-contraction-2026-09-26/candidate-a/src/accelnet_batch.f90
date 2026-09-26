! Batch evaluation is deliberately separate from the established CPU APIs.
! CSR rows describe central atoms; targets may include ghosts or repeated
! periodic images. Displacements must already include the image translation.
module accelnet_batch
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model
    use accelnet_batch_target_serial, only: target_model, target_workspace, evaluate_batch_target
    use accelnet_descriptor_models, only: evaluate_model_values, evaluate_model_values_derivatives, &
        contract_model_derivatives, model_supports_direct_contraction
    implicit none
    private
    public :: batch_workspace, evaluate_batch, evaluate_batch_reference

    type :: batch_workspace
        private
        real(real64), allocatable :: descriptor(:), normalized(:), gradient(:), coefficients(:)
        real(real64), allocatable :: center_derivatives(:,:), neighbor_derivatives(:,:,:)
        real(real64), allocatable :: neighbor_forces(:,:)
        integer, allocatable :: local_species(:), global_species(:)
        logical :: direct = .true.
        integer :: dimension_capacity = 0, neighbor_capacity = 0, growth_count = 0
        type(target_workspace) :: shared_work
    contains
        procedure :: reserve => reserve_workspace
        procedure :: release => release_workspace
        procedure :: allocations => workspace_allocations
    end type
contains
    ! One growth event can allocate several arrays. No event occurs when the
    ! existing capacities suffice, including after a model reload of equal size.
    integer function workspace_allocations(self) result(n)
        class(batch_workspace), intent(in) :: self
        n = self%growth_count + self%shared_work%allocations()
    end function

    subroutine release_workspace(self)
        class(batch_workspace), intent(inout) :: self
        type(batch_workspace) :: empty
        select type (self)
        type is (batch_workspace)
            self = empty
        end select
    end subroutine

    subroutine reserve_workspace(self, model, maximum_neighbors)
        class(batch_workspace), intent(inout) :: self
        type(predictor_model), intent(in) :: model
        integer, intent(in) :: maximum_neighbors
        integer :: dimension, species, d, n
        logical :: grew, need_jacobian
        if (.not. allocated(model%networks) .or. .not. allocated(model%setups)) &
            error stop 'batch: model is not initialized'
        if (size(model%networks) == 0 .or. size(model%networks) /= size(model%setups)) &
            error stop 'batch: invalid model species'
        if (maximum_neighbors < 0) error stop 'batch: negative neighbor capacity'
        dimension = 0
        grew = .false.
        ! Refresh metadata on every call: the caller may reload/change a model.
        self%direct = .true.
        do species = 1, size(model%networks)
            dimension = max(dimension, model%networks(species)%nodes(1))
            self%direct = self%direct .and. model_supports_direct_contraction(model%setups(species)%model)
        end do
        ! Match the established structure CPU policy. Mixing per-species paths
        ! made small G4/G5 models slower (G5 values/geometry evaluated twice).
        ! GPU species grouping will be a separate, measured backend decision.
        d = max(1, dimension, self%dimension_capacity)
        n = max(1, maximum_neighbors, self%neighbor_capacity)
        if (.not. allocated(self%descriptor) .or. d > self%dimension_capacity) then
            if (allocated(self%descriptor)) &
                deallocate(self%descriptor, self%normalized, self%gradient, self%coefficients)
            allocate(self%descriptor(d), self%normalized(d), self%gradient(d), self%coefficients(d))
            grew = .true.
        end if
        if (.not. allocated(self%local_species) .or. n > self%neighbor_capacity) then
            if (allocated(self%local_species)) &
                deallocate(self%local_species, self%global_species, self%neighbor_forces)
            allocate(self%local_species(n), self%global_species(n), self%neighbor_forces(3,n))
            grew = .true.
        end if
        need_jacobian = .not. self%direct
        if (allocated(self%center_derivatives)) then
            if (size(self%center_derivatives,2) < d .or. size(self%neighbor_derivatives,3) < n) &
                deallocate(self%center_derivatives, self%neighbor_derivatives)
        end if
        if (need_jacobian .and. .not. allocated(self%center_derivatives)) then
            allocate(self%center_derivatives(3,d), self%neighbor_derivatives(3,d,n))
            grew = .true.
        end if
        self%dimension_capacity = d
        self%neighbor_capacity = n
        if (grew) self%growth_count = self%growth_count + 1
    end subroutine

    subroutine evaluate_batch(model, species, centers, offsets, indices, displacements, energies, forces, work, virial)
        type(predictor_model), intent(in) :: model
        integer, intent(in) :: species(:), centers(:), offsets(:), indices(:)
        real(real64), intent(in) :: displacements(:,:)
        real(real64), intent(out) :: energies(:)
        real(real64), intent(inout) :: forces(:,:)
        type(batch_workspace), intent(inout) :: work
        real(real64), intent(inout), optional :: virial(3,3)
        type(target_model) :: packed
        integer :: s, status
        logical :: chebyshev

        chebyshev = allocated(model%setups) .and. allocated(model%networks)
        if (chebyshev) then
            chebyshev = size(model%setups) > 0
            do s = 1, size(model%setups)
                if (.not. allocated(model%setups(s)%model%chebyshev)) then
                    chebyshev = .false.
                    exit
                end if
                if (size(model%setups(s)%model%chebyshev) /= 1 .or. &
                    allocated(model%setups(s)%model%lj) .or. allocated(model%setups(s)%model%behler)) then
                    chebyshev = .false.
                    exit
                end if
            end do
        end if
        if (chebyshev) then
            ! Repack to honor direct edits/reloads of the public model fields.
            ! The persistent workspace reuses buffers and unchanged model data.
            call packed%initialize(model,use_host=.true.,status=status)
            if (status == 0) then
                call evaluate_batch_target(packed,species,centers,offsets,indices,displacements, &
                    energies,forces,work%shared_work,virial)
                return
            end if
        end if
        call evaluate_batch_reference(model,species,centers,offsets,indices,displacements,energies,forces,work,virial)
    end subroutine

    ! Retained as an independent numerical/performance reference and as the
    ! fallback for descriptor families not yet selected for shared CPU kernels.
    subroutine evaluate_batch_reference(model, species, centers, offsets, indices, displacements, energies, forces, work, virial)
        type(predictor_model), intent(in) :: model
        integer, intent(in) :: species(:), centers(:), offsets(:), indices(:)
        real(real64), intent(in) :: displacements(:,:)
        real(real64), intent(out) :: energies(:)
        real(real64), intent(inout) :: forces(:,:)
        type(batch_workspace), intent(inout) :: work
        real(real64), intent(inout), optional :: virial(3,3)
        integer :: row, atom, kind, first, last, n, dimension, j, target, component, coefficient, max_neighbors
        real(real64) :: network_energy, center_force(3), neighbor_force(3)

        ! Validate the CSR boundary before indexing it. Offsets may refer to a
        ! subrange of the edge arrays, allowing zero-copy partitioning of rows.
        if (size(offsets) /= size(centers)+1) error stop 'batch: wrong CSR offset count'
        if (size(energies) /= size(centers)) error stop 'batch: wrong energy count'
        if (size(forces,1) /= 3 .or. size(forces,2) /= size(species)) error stop 'batch: wrong force shape'
        if (size(displacements,1) /= 3 .or. size(displacements,2) /= size(indices)) &
            error stop 'batch: wrong displacement shape'
        if (any(offsets < 1) .or. any(offsets > size(indices)+1)) error stop 'batch: CSR offset out of range'
        if (any(offsets(2:) < offsets(:size(centers)))) error stop 'batch: CSR offsets must be nondecreasing'
        if (any(centers < 1) .or. any(centers > size(species))) error stop 'batch: center index out of range'
        first = offsets(1)
        last = offsets(size(offsets))-1
        if (any(indices(first:last) < 1) .or. any(indices(first:last) > size(species))) &
            error stop 'batch: neighbor index out of range'
        max_neighbors = 0
        if (size(centers) > 0) max_neighbors = maxval(offsets(2:) - offsets(:size(centers)))
        call work%reserve(model, max_neighbors)
        if (any(species < 1) .or. any(species > size(model%networks))) error stop 'batch: species out of range'

        do row = 1, size(centers)
            atom = centers(row)
            kind = species(atom)
            dimension = model%networks(kind)%nodes(1)
            first = offsets(row)
            last = offsets(row+1)-1
            n = last-first+1
            work%global_species(1:n) = species(indices(first:last))
            call model%setups(kind)%map_species(work%global_species(1:n), work%local_species(1:n))
            if (work%direct) then
                call evaluate_model_values(model%setups(kind)%model, displacements(:,first:last), &
                    work%local_species(1:n), work%descriptor(1:dimension))
            else
                call evaluate_model_values_derivatives(model%setups(kind)%model, displacements(:,first:last), &
                    work%local_species(1:n), work%descriptor(1:dimension), &
                    work%center_derivatives(:,1:dimension), work%neighbor_derivatives(:,1:dimension,1:n))
            end if
            work%normalized(1:dimension) = (work%descriptor(1:dimension) - model%networks(kind)%descriptor_shift) * &
                model%networks(kind)%descriptor_scale
            call model%networks(kind)%input_gradient(work%normalized(1:dimension), network_energy, &
                work%gradient(1:dimension))
            energies(row) = network_energy/model%networks(kind)%energy_scale + model%networks(kind)%energy_shift + &
                model%networks(kind)%atomic_references(kind)
            work%coefficients(1:dimension) = -work%gradient(1:dimension) * &
                model%networks(kind)%descriptor_scale/model%networks(kind)%energy_scale
            if (work%direct) then
                call contract_model_derivatives(model%setups(kind)%model, displacements(:,first:last), &
                    work%local_species(1:n), work%coefficients(1:dimension), center_force, work%neighbor_forces(:,1:n))
                forces(:,atom) = forces(:,atom) + center_force
            else
                do coefficient = 1, dimension
                    forces(:,atom) = forces(:,atom) + &
                        work%coefficients(coefficient)*work%center_derivatives(:,coefficient)
                end do
            end if
            do j = 1, n
                target = indices(first+j-1)
                if (work%direct) then
                    neighbor_force = work%neighbor_forces(:,j)
                    forces(:,target) = forces(:,target) + neighbor_force
                else
                    do coefficient = 1, dimension
                        forces(:,target) = forces(:,target) + &
                            work%coefficients(coefficient)*work%neighbor_derivatives(:,coefficient,j)
                    end do
                    if (present(virial)) neighbor_force = matmul(work%neighbor_derivatives(:,1:dimension,j), &
                        work%coefficients(1:dimension))
                end if
                ! Compute virial using image displacements BEFORE folding forces
                ! onto target IDs. No half factor and no volume normalization.
                if (present(virial)) then
                    do component = 1, 3
                        virial(:,component) = virial(:,component) + displacements(:,first+j-1)*neighbor_force(component)
                    end do
                end if
            end do
        end do
    end subroutine
end module accelnet_batch
