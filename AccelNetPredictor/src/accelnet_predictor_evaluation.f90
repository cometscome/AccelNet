! Public structure/file inference delegates to the common serial batch kernel.
submodule (accelnet_predictor) accelnet_predictor_evaluation
    use accelnet_batch, only: batch_workspace, evaluate_batch
contains
    module procedure predict_energy_structure
        real(real64), allocatable :: forces(:,:)
        allocate(forces(3,structure%natoms))
        call evaluate_common_structure(self,structure,total_energy,forces,energy_only=.true.)
    end procedure

    module procedure predict_energy_forces_structure
        call evaluate_common_structure(self,structure,total_energy,forces,virial)
    end procedure

    subroutine evaluate_common_structure(self,structure,total_energy,forces,virial,energy_only)
        class(predictor_model), intent(in) :: self
        type(atomic_structure), intent(in) :: structure
        real(real64), intent(out) :: total_energy,forces(:,:)
        real(real64), optional, intent(out) :: virial(3,3)
        logical, optional, intent(in) :: energy_only
        type(neighbor_data) :: neighbors
        type(batch_workspace) :: work
        real(real64), allocatable :: energies(:)
        integer, allocatable :: centers(:)
        integer :: i
        if (size(forces,1)/=3 .or. size(forces,2)/=structure%natoms) &
            error stop 'forces must have shape (3, structure%natoms)'
        call build_neighbor_list(structure,self%maximum_cutoff,neighbors,self%minimum_distance)
        allocate(energies(structure%natoms),centers(structure%natoms))
        centers=[(i,i=1,structure%natoms)]
        forces=0
        if (present(virial)) virial=0
        call evaluate_batch(self,structure%species,centers,neighbors%offsets,neighbors%atom_indices, &
            neighbors%displacements,energies,forces,work,virial,energy_only)
        total_energy=sum(energies)
    end subroutine
end submodule
