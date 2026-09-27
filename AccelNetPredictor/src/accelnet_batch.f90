! All public batch inference uses the common numerical implementation.
module accelnet_batch
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model
    use accelnet_batch_target_serial, only: target_model,target_workspace,evaluate_batch_target
    use accelnet_cpu_reference, only: reference_workspace=>batch_workspace, &
        legacy_evaluate=>evaluate_batch_reference
    implicit none
    private
    public :: batch_workspace,evaluate_batch,evaluate_batch_reference
    type :: batch_workspace
        private
        type(target_workspace) :: shared_work
        type(reference_workspace) :: reference_work
    contains
        procedure :: reserve => reserve_workspace
        procedure :: release => release_workspace
        procedure :: allocations => workspace_allocations
    end type
contains
    subroutine reserve_workspace(self,model,maximum_neighbors)
        class(batch_workspace), intent(inout) :: self
        class(predictor_model), intent(in) :: model
        integer, intent(in) :: maximum_neighbors
        ! Preserve explicit reference preallocation; common scratch uses actual CSR sizes.
        call self%reference_work%reserve(model,maximum_neighbors)
    end subroutine
    subroutine release_workspace(self)
        class(batch_workspace), intent(inout) :: self
        call self%shared_work%release()
        call self%reference_work%release()
    end subroutine
    integer function workspace_allocations(self) result(n)
        class(batch_workspace), intent(in) :: self
        n=self%shared_work%allocations()+self%reference_work%allocations()
    end function
    subroutine evaluate_batch(model,species,centers,offsets,indices,displacements,energies,forces,work,virial,energy_only)
        class(predictor_model), intent(in) :: model
        integer, intent(in) :: species(:),centers(:),offsets(:),indices(:)
        real(real64), intent(in) :: displacements(:,:)
        real(real64), intent(out) :: energies(:)
        real(real64), intent(inout) :: forces(:,:)
        type(batch_workspace), intent(inout) :: work
        real(real64), optional, intent(inout) :: virial(3,3)
        logical, optional, intent(in) :: energy_only
        type(target_model) :: packed
        call packed%initialize(model,use_host=.true.)
        call evaluate_batch_target(packed,species,centers,offsets,indices,displacements,energies,forces, &
            work%shared_work,virial,energy_only=energy_only)
    end subroutine
    ! Compatibility spelling, never a production fallback.
    subroutine evaluate_batch_reference(model,species,centers,offsets,indices,displacements,energies,forces,work,virial)
        type(predictor_model), intent(in) :: model
        integer, intent(in) :: species(:),centers(:),offsets(:),indices(:)
        real(real64), intent(in) :: displacements(:,:)
        real(real64), intent(out) :: energies(:)
        real(real64), intent(inout) :: forces(:,:)
        type(batch_workspace), intent(inout) :: work
        real(real64), optional, intent(inout) :: virial(3,3)
        call legacy_evaluate(model,species,centers,offsets,indices,displacements,energies,forces,work%reference_work,virial)
    end subroutine
end module
