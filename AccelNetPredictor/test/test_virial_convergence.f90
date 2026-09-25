! Independent energy-only central differences of all nine strain components.
! The CSV records every component and step, including the roundoff-dominated tail.
program test_virial_convergence
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet
    use accelnet_descriptors, only: atomic_structure, neighbor_data, build_neighbor_list, read_xsf
    use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2, load_predictor_from_networks
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: structure, sheared
    character(len=1024) :: format, directory, networks(2), xsf, csv_path
    real(real64) :: deformation(3,3)
    integer :: csv, stat, mode, geometry, i

    call get_command_argument(1, format)
    select case(trim(format))
    case('--n2p2')
        if (command_argument_count() /= 3) error stop 'usage: --n2p2 model_directory output.csv'
        call get_command_argument(2, directory)
        call get_command_argument(3, csv_path)
        call load_predictor_from_n2p2(trim(directory), model)
        call accelnet_init_n2p2(trim(directory), stat)
    case('--aenet')
        if (command_argument_count() /= 5) error stop 'usage: --aenet Ti.nn O.nn structure.xsf output.csv'
        call get_command_argument(2, networks(1))
        call get_command_argument(3, networks(2))
        call get_command_argument(4, xsf)
        call get_command_argument(5, csv_path)
        call load_predictor_from_networks(networks, model)
        call read_xsf(trim(xsf), model%species_names, structure)
        call accelnet_init(model%species_names, stat)
        call require(stat == ACCELNET_OK, 'aenet init')
        do i=1,2
            call accelnet_load_potential(i, trim(networks(i)), stat)
            call require(stat == ACCELNET_OK, 'aenet load')
        end do
    case default
        error stop 'expected --n2p2 or --aenet'
    end select
    call require(stat == ACCELNET_OK, 'model initialization')
    open(newunit=csv, file=trim(csv_path), status='replace', action='write')
    write(csv,'(A)') 'case,mode,natoms,h,a,b,analytic,numerical,abs_error,tolerance'
    if (trim(format) == '--n2p2') then
        structure%natoms=3
        allocate(structure%positions(3,3), structure%species(3))
        structure%species=1
        structure%species(2)=min(2,size(model%species_names))
        structure%positions(:,1)=[0.05_real64,0.06_real64,0.07_real64]*model%maximum_cutoff
        structure%positions(:,2)=[0.45_real64,0.49_real64,0.51_real64]*model%maximum_cutoff
        structure%positions(:,3)=[0.62_real64,0.16_real64,0.31_real64]*model%maximum_cutoff
        do geometry=1,3
            structure%pbc=geometry /= 1
            structure%lattice=0.0_real64
            structure%lattice(1,1)=0.85_real64
            structure%lattice(2,2)=0.92_real64
            structure%lattice(3,3)=0.88_real64
            if (geometry == 3) then
                structure%lattice(1,2)=0.17_real64
                structure%lattice(1,3)=0.08_real64
                structure%lattice(2,3)=0.14_real64
            end if
            structure%lattice=structure%lattice*model%maximum_cutoff
            select case(geometry)
            case(1)
                call sweep(structure, 'molecule', 0)
            case(2)
                call sweep(structure, 'orthogonal', 0)
            case(3)
                call sweep(structure, 'triclinic', 0)
            end select
        end do
        structure%natoms=1
        structure%positions=structure%positions(:,:1)
        structure%species=[1]
        call sweep(structure, 'self_images', 0)
    else
        deformation=0.0_real64
        do i=1,3
            deformation(i,i)=1.0_real64
        end do
        deformation(1,2)=0.07_real64
        deformation(2,3)=-0.04_real64
        sheared=structure
        sheared%lattice=matmul(deformation,structure%lattice)
        sheared%positions=matmul(deformation,structure%positions)
        do mode=ACCELNET_CHEBYSHEV_AUTO,ACCELNET_CHEBYSHEV_MOMENT
            call model%set_chebyshev_evaluation(mode)
            call accelnet_set_chebyshev_evaluation(mode,stat)
            call require(stat == ACCELNET_OK, 'evaluation mode')
            call sweep(structure, 'original', mode)
            call sweep(sheared, 'sheared', mode)
        end do
    end if
    close(csv)
    call accelnet_final(stat)
    call require(stat == ACCELNET_OK, 'finalize')
contains
    subroutine sweep(s,label,mode)
        type(atomic_structure), intent(in) :: s
        character(len=*), intent(in) :: label
        integer, intent(in) :: mode
        real(real64), parameter :: steps(*)=[1e-2_real64,3e-3_real64,1e-3_real64,3e-4_real64, &
            1e-4_real64,3e-5_real64,1e-5_real64,3e-6_real64,1e-6_real64,3e-7_real64, &
            1e-7_real64,3e-8_real64,1e-8_real64]
        ! Per-component tolerance does not scale with the large reference energy.
        real(real64), parameter :: absolute_tolerance=2e-6_real64, relative_tolerance=2e-7_real64
        type(atomic_structure) :: perturbed
        type(neighbor_data) :: neighbors
        real(real64) :: energy, eplus, eminus, ei, atomic_energy, h, transform(3,3)
        real(real64) :: force(3,s%natoms), atomic_force(3,s%natoms), w(3,3), atomic_w(3,3), fd(3,3)
        real(real64) :: tolerances(3,3), errors(size(steps)), scaled_errors(size(steps)), api_error
        integer :: k,a,b,j,atom,first,last,n,status,best
        call model%predict_energy_forces(s,energy,force,w)
        call require(all(ieee_is_finite(w)), 'finite analytical virial')
        atomic_force=0.0_real64; atomic_w=0.0_real64; atomic_energy=0.0_real64
        call build_neighbor_list(s,model%maximum_cutoff,neighbors,model%minimum_distance)
        do atom=1,s%natoms
            first=neighbors%offsets(atom); last=neighbors%offsets(atom+1)-1; n=last-first+1
            call accelnet_atomic_energy_and_forces_virial(s%positions(:,atom),s%species(atom),atom,n, &
                neighbors%positions(:,first:last),s%species(neighbors%atom_indices(first:last)), &
                neighbors%atom_indices(first:last),s%natoms,ei,atomic_force,atomic_w,status)
            call require(status == ACCELNET_OK, 'atomic virial')
            atomic_energy=atomic_energy+ei
        end do
        api_error=maxval(abs(w-atomic_w))
        call require(abs(energy-atomic_energy) < 1e-8_real64, 'atomic energy parity')
        call require(maxval(abs(force-atomic_force)) < 1e-9_real64, 'atomic force parity')
        call require(api_error < 1e-9_real64, 'atomic virial parity')
        tolerances=absolute_tolerance+relative_tolerance*abs(w)
        do k=1,size(steps)
            h=steps(k)
            do a=1,3
                do b=1,3
                    transform=0.0_real64
                    do j=1,3
                        transform(j,j)=1.0_real64
                    end do
                    ! Coordinates AND cell change at fixed fractional coordinates.
                    ! W(a,b)=-dE/d epsilon(b,a); only energy is used for the reference.
                    transform(b,a)=transform(b,a)+h
                    perturbed=s
                    perturbed%positions=matmul(transform,s%positions)
                    perturbed%lattice=matmul(transform,s%lattice)
                    call model%predict_energy(perturbed,eplus)
                    transform(b,a)=transform(b,a)-2*h
                    perturbed%positions=matmul(transform,s%positions)
                    perturbed%lattice=matmul(transform,s%lattice)
                    call model%predict_energy(perturbed,eminus)
                    fd(a,b)=-(eplus-eminus)/(2*h)
                    write(csv,'(A,",",I0,",",I0,",",ES24.16,",",I0,",",I0,4(",",ES24.16))') &
                        label,mode,s%natoms,h,a,b,w(a,b),fd(a,b),abs(w(a,b)-fd(a,b)),tolerances(a,b)
                end do
            end do
            call require(all(ieee_is_finite(fd)), 'finite numerical virial')
            errors(k)=maxval(abs(w-fd))
            scaled_errors(k)=maxval(abs(w-fd)/tolerances)
        end do
        best=minloc(errors,dim=1)
        write(*,'(A,1X,A,1X,I0,1X,A,ES10.3,1X,A,ES10.3,1X,A,ES10.3)') &
            'VIRIAL',label,mode,'best_h=',steps(best),'max_error=',errors(best),'atomic_error=',api_error
        ! Require multiple neighboring useful steps, not an accidentally good
        ! component-specific step. Very small h is allowed to amplify roundoff.
        call require(any(scaled_errors(:size(steps)-1) <= 1 .and. scaled_errors(2:) <= 1), &
            'two adjacent step sizes agree in all nine components')
        if (errors(1) > 10*maxval(tolerances)) &
            call require(minval(errors) < errors(1)/100, 'convergence from coarse step')
    end subroutine

    subroutine require(condition,label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*,'(2A)') 'FAILED: ',label
            error stop 1
        end if
    end subroutine
end program
