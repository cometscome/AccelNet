program test_virial_invariants
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet
    use accelnet_descriptors, only: atomic_structure, neighbor_data, build_neighbor_list
    use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2
    use accelnet_descriptor_models, only: model_supports_direct_contraction
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: s, transformed
    character(len=1024) :: directory, mixed_paths
    real(real64) :: e, et, f(3,4), ft(3,4), w(3,3), wt(3,3), rotation(3,3), angle
    integer, parameter :: permutation(4)=[3,1,4,2]
    integer :: stat, geometry, mode
    call get_command_argument(1,directory)
    call get_command_argument(2,mixed_paths)
    call load_predictor_from_n2p2(trim(directory),model)
    call accelnet_init_n2p2(trim(directory),stat)
    call require(stat == ACCELNET_OK,'init')
    if (trim(mixed_paths) == 'mixed-paths') then
        call require(model_supports_direct_contraction(model%setups(1)%model), 'H contraction path')
        call require(.not. model_supports_direct_contraction(model%setups(2)%model), 'O G4 Jacobian path')
    end if
    s%natoms=4
    allocate(s%positions(3,4),s%species(4))
    s%species=[1,2,1,2]
    angle=0.37_real64
    rotation=reshape([cos(angle),sin(angle),0.0_real64,-sin(angle),cos(angle),0.0_real64, &
                      0.0_real64,0.0_real64,1.0_real64],[3,3])
    do mode=ACCELNET_G5_DIRECT,ACCELNET_G5_MOMENT_FORCE
        call model%set_g5_evaluation(mode)
        call accelnet_set_g5_evaluation(mode,stat)
        call require(stat == ACCELNET_OK,'G5 mode')
        do geometry=1,3
            s%pbc=geometry /= 1
            s%positions=reshape([0.2_real64,0.3_real64,0.4_real64, 1.4_real64,0.6_real64,0.9_real64, &
                0.7_real64,1.8_real64,1.1_real64, 2.2_real64,1.3_real64,2.0_real64],[3,4])
            s%lattice=0.0_real64
            s%lattice(1,1)=3.2_real64
            s%lattice(2,2)=3.4_real64
            s%lattice(3,3)=3.6_real64
            if (geometry == 3) then
                s%lattice(1,2)=0.7_real64
                s%lattice(1,3)=-0.4_real64
                s%lattice(2,3)=0.5_real64
            end if
            call evaluate(s,e,f,w)
            call require(maxval(abs(w)) > 1e-6_real64,'nontrivial virial')
            call require(maxval(abs(sum(f,dim=2))) < 1e-10_real64,'zero net force')
            call require(maxval(abs(w-transpose(w))) < 1e-10_real64,'symmetric virial')
            call check_strain(s,w)

            transformed=s
            transformed%species=s%species(permutation)
            transformed%positions=s%positions(:,permutation)
            call evaluate(transformed,et,ft,wt)
            call compare(e,f(:,permutation),w,et,ft,wt,'atom/species permutation')

            transformed=s
            transformed%positions=matmul(rotation,s%positions)
            transformed%lattice=matmul(rotation,s%lattice)
            call evaluate(transformed,et,ft,wt)
            call compare(e,matmul(rotation,f),matmul(rotation,matmul(w,transpose(rotation))), &
                         et,ft,wt,'rigid rotation: W transforms as R W R^T')

            if (s%pbc) then
                ! Independent shifts deliberately exceed the cutoff search range.
                transformed=s
                transformed%positions(:,1)=s%positions(:,1)+matmul(s%lattice,[5.0_real64,-4.0_real64,3.0_real64])
                transformed%positions(:,3)=s%positions(:,3)+matmul(s%lattice,[-6.0_real64,2.0_real64,-5.0_real64])
                call evaluate(transformed,et,ft,wt)
                call compare(e,f,w,et,ft,wt,'independent periodic image shifts')
            end if
        end do
    end do
    call accelnet_final(stat)
    call require(stat == ACCELNET_OK,'final')
    print *, 'Multispecies virial permutations, rotations and periodic images passed'
contains
    subroutine check_strain(s,w)
        type(atomic_structure), intent(in) :: s
        real(real64), intent(in) :: w(3,3)
        type(atomic_structure) :: perturbed
        real(real64), parameter :: h=1e-5_real64
        real(real64) :: deformation(3,3),ep,em,fd
        integer :: a,b,i
        do a=1,3
            do b=1,3
                deformation=0.0_real64
                do i=1,3
                    deformation(i,i)=1.0_real64
                end do
                deformation(b,a)=deformation(b,a)+h
                perturbed=s
                perturbed%positions=matmul(deformation,s%positions)
                perturbed%lattice=matmul(deformation,s%lattice)
                call model%predict_energy(perturbed,ep)
                deformation(b,a)=deformation(b,a)-2*h
                perturbed%positions=matmul(deformation,s%positions)
                perturbed%lattice=matmul(deformation,s%lattice)
                call model%predict_energy(perturbed,em)
                fd=-(ep-em)/(2*h)
                call require(abs(w(a,b)-fd) <= 2e-6_real64+2e-7_real64*abs(w(a,b)), &
                    'G5 evaluation mode strain derivative')
            end do
        end do
    end subroutine

    subroutine evaluate(s,e,f,w)
        type(atomic_structure), intent(in) :: s
        real(real64), intent(out) :: e,f(3,4),w(3,3)
        type(neighbor_data) :: neighbors
        real(real64) :: ea,ei,fa(3,4),wa(3,3),image_position(3)
        integer :: atom,first,last,n,status,j
        call model%predict_energy_forces(s,e,f,w)
        call require(ieee_is_finite(e) .and. all(ieee_is_finite(f)) .and. all(ieee_is_finite(w)), 'finite outputs')
        call build_neighbor_list(s,model%maximum_cutoff,neighbors,model%minimum_distance)
        ea=0.0_real64; fa=0.0_real64; wa=0.0_real64
        do atom=1,s%natoms
            first=neighbors%offsets(atom); last=neighbors%offsets(atom+1)-1; n=last-first+1
            do j=first,last
                image_position=s%positions(:,neighbors%atom_indices(j))
                if (s%pbc) image_position=image_position+matmul(s%lattice,real(neighbors%image_shifts(:,j),real64))
                call require(maxval(abs(image_position-neighbors%positions(:,j))) < 1e-11_real64, 'image metadata')
                call require(maxval(abs(neighbors%displacements(:,j) - &
                    (image_position-s%positions(:,atom)))) < 1e-11_real64,'image displacement metadata')
            end do
            call accelnet_atomic_energy_and_forces_virial(s%positions(:,atom),s%species(atom),atom,n, &
                neighbors%positions(:,first:last),s%species(neighbors%atom_indices(first:last)), &
                neighbors%atom_indices(first:last),s%natoms,ei,fa,wa,status)
            call require(status == ACCELNET_OK,'atomic API')
            ea=ea+ei
        end do
        call compare(e,f,w,ea,fa,wa,'atomic/structure API (mixed derivative paths)')
    end subroutine

    subroutine compare(e,f,w,et,ft,wt,label)
        real(real64), intent(in) :: e,f(3,4),w(3,3),et,ft(3,4),wt(3,3)
        character(len=*), intent(in) :: label
        if (max(abs(e-et),maxval(abs(f-ft)),maxval(abs(w-wt))) > 2e-9_real64) then
            print *, 'Errors E/F/W:',abs(e-et),maxval(abs(f-ft)),maxval(abs(w-wt))
            call require(.false.,label)
        end if
    end subroutine

    subroutine require(condition,label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            print *, 'FAILED: ',label
            error stop 1
        end if
    end subroutine
end program
