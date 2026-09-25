program test_virial
    use iso_fortran_env, only: real64
    use accelnet
    use accelnet_descriptors, only: atomic_structure, neighbor_data, build_neighbor_list
    use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2, load_predictor_from_networks
    use aenet_network, only: atomic_network, write_aenet_network_ascii
    use accelnet_descriptor_models, only: model_supports_direct_contraction
    implicit none
    type(predictor_model) :: model
    character(len=1024) :: directory
    integer :: stat, mode, sweep_unit, case_number = 0

    open(newunit=sweep_unit, file='virial-convergence-synthetic.csv', status='replace', action='write')
    write(sweep_unit,'(A)') 'case,model,mode,natoms,periodic,h,max_abs_error,max_scaled_error'
    call get_command_argument(1, directory)
    call load_predictor_from_n2p2(trim(directory), model)
    call accelnet_init_n2p2(trim(directory), stat)
    call require(stat == ACCELNET_OK, 'n2p2 init')
    call check_geometries()
    call check_errors()
    call accelnet_final(stat)
    call require(stat == ACCELNET_OK, 'n2p2 final')

    ! Synthetic ASCII networks make both derivative paths testable without
    ! an external model corpus. Chebyshev also exercises angular forces.
    call synthetic_model('Chebyshev', 6)
    call require(model_supports_direct_contraction(model%setups(1)%model), 'contraction path')
    do mode = ACCELNET_CHEBYSHEV_AUTO, ACCELNET_CHEBYSHEV_MOMENT
        call accelnet_set_chebyshev_evaluation(mode, stat)
        call require(stat == ACCELNET_OK, 'Chebyshev mode')
        call model%set_chebyshev_evaluation(mode)
        call check_geometries()
    end do
    call accelnet_final(stat)

    call synthetic_model('LJ', 2)
    call require(.not. model_supports_direct_contraction(model%setups(1)%model), 'Jacobian path')
    call check_geometries()
    call accelnet_final(stat)
    call check_errors(uninitialized=.true.)
    close(sweep_unit)
    print *, 'Atomic and structure virial finite differences passed'
contains
    subroutine synthetic_model(name, dimension)
        character(len=*), intent(in) :: name
        integer, intent(in) :: dimension
        type(atomic_network) :: network
        integer :: i
        network%atomtype = 'H'
        network%description = 'Synthetic virial regression model'
        network%descriptor_name = name
        network%minimum_radius = 0.35_real64
        network%maximum_radius = 3.4_real64
        network%environment_names = [character(len=16) :: 'H']
        network%species_names = network%environment_names
        network%atomic_references = [0.3_real64]
        network%energy_scale = 2.5_real64
        network%energy_shift = 0.7_real64
        network%nlayers = 2
        network%maxnodes = dimension
        network%nodes = [dimension, 1]
        network%activation = [0]
        network%weight_offsets = [0, dimension + 1]
        network%weights = [(0.1_real64*i, i=1,dimension+1)]
        allocate(network%descriptor_kinds(dimension), network%descriptor_environments(2,dimension), &
                 network%descriptor_parameters(6,dimension), network%descriptor_shift(dimension), &
                 network%descriptor_scale(dimension))
        network%descriptor_kinds = 1
        network%descriptor_environments = 1
        network%descriptor_parameters = 0.0_real64
        network%descriptor_parameters(1,:) = network%maximum_radius
        network%descriptor_parameters(2,:) = 2.0_real64
        network%descriptor_parameters(3,:) = network%maximum_radius
        network%descriptor_parameters(4,:) = 2.0_real64
        network%descriptor_parameters(5,:) = 1.0_real64
        network%descriptor_shift = 0.1_real64
        network%descriptor_scale = 0.8_real64
        call write_aenet_network_ascii('virial-test.nn.ascii', network)
        call load_predictor_from_networks(['virial-test.nn.ascii'], model)
        call accelnet_init(network%species_names, stat)
        call require(stat == ACCELNET_OK, 'synthetic init')
        call accelnet_load_potential(1, 'virial-test.nn.ascii', stat, is_ascii=.true.)
        call require(stat == ACCELNET_OK, 'synthetic load')
    end subroutine

    subroutine check_geometries()
        type(atomic_structure) :: s
        integer :: geometry
        s%natoms = 2
        allocate(s%positions(3,2), s%species(2))
        s%species = 1
        s%positions(:,1) = [0.05_real64, 0.06_real64, 0.07_real64]*model%maximum_cutoff
        s%positions(:,2) = [0.45_real64, 0.49_real64, 0.51_real64]*model%maximum_cutoff
        do geometry = 1, 3
            s%pbc = geometry /= 1
            s%lattice = 0.0_real64
            s%lattice(1,1) = 0.85_real64
            s%lattice(2,2) = 0.92_real64
            s%lattice(3,3) = 0.88_real64
            if (geometry == 3) then
                s%lattice(1,2) = 0.17_real64
                s%lattice(1,3) = 0.08_real64
                s%lattice(2,3) = 0.14_real64
            end if
            s%lattice = s%lattice*model%maximum_cutoff
            call check_structure(s)
        end do
        ! All periodic neighbors now fold onto the central atom: the total
        ! force vanishes, while cell derivatives (and virial) do not.
        s%natoms = 1
        s%positions = s%positions(:,:1)
        s%species = [1]
        call check_structure(s)
        ! An isolated atom has no neighbors and zero virial.
        s%pbc = .false.
        call check_structure(s)
    end subroutine

    subroutine check_structure(s)
        type(atomic_structure), intent(in) :: s
        type(atomic_structure) :: shifted
        real(real64) :: e, ea, old_e, ep, em, w(3,3), wa(3,3), ws(3,3), fd(3,3), transform(3,3)
        real(real64) :: f(3,s%natoms), fa(3,s%natoms), old_f(3,s%natoms), fs(3,s%natoms)
        real(real64), allocatable :: file_forces(:,:)
        real(real64), parameter :: steps(*)=[1e-2_real64,3e-3_real64,1e-3_real64,3e-4_real64, &
            1e-4_real64,3e-5_real64,1e-5_real64,3e-6_real64,1e-6_real64,3e-7_real64, &
            1e-7_real64,3e-8_real64,1e-8_real64]
        real(real64) :: h, force_fd, scale, scaled_errors(size(steps)), tolerance(3,3)
        integer :: a, b, i, unit, k
        h = 1.0e-6_real64
        call model%predict_energy_forces(s, e, f, w)
        call model%predict_energy_forces(s, old_e, old_f)
        call require(abs(e-old_e) < 1.0e-12_real64, 'object energy compatibility')
        call require(maxval(abs(f-old_f)) < 1.0e-12_real64, 'object force compatibility')
        fa = 0.0_real64; wa = 0.0_real64
        call atomic_sum(s, ea, fa, wa)
        scale = max(1.0_real64, abs(e), maxval(abs(f)), maxval(abs(w)))
        call require(abs(e-ea) < 1.0e-11_real64*scale, 'atomic/object energy')
        call require(maxval(abs(f-fa)) < 1.0e-11_real64*scale, 'atomic/object forces')
        call require(maxval(abs(w-wa)) < 1.0e-11_real64*scale, 'atomic/object virial')
        call require(maxval(abs(w-transpose(w))) < 1.0e-10_real64*scale, 'virial symmetry')
        fa = fa + 0.25_real64; wa = wa + 0.5_real64
        call atomic_sum(s, ea, fa, wa)
        call require(maxval(abs(fa-2*f-0.25_real64)) < 1.0e-10_real64*scale, 'additive forces')
        call require(maxval(abs(wa-2*w-0.5_real64)) < 1.0e-10_real64*scale, 'additive virial')
        case_number=case_number+1
        tolerance=2e-6_real64+2e-7_real64*abs(w)
        do k=1,size(steps)
            h=steps(k)
            do a=1,3
                do b=1,3
                    transform=0.0_real64
                    do i=1,3
                        transform(i,i)=1.0_real64
                    end do
                    ! W(a,b) = -dE/d epsilon(b,a), r'=(I+epsilon) r.
                    transform(b,a)=transform(b,a)+h
                    shifted=s
                    shifted%positions=matmul(transform,s%positions)
                    shifted%lattice=matmul(transform,s%lattice)
                    call model%predict_energy(shifted,ep)
                    transform(b,a)=transform(b,a)-2*h
                    shifted%positions=matmul(transform,s%positions)
                    shifted%lattice=matmul(transform,s%lattice)
                    call model%predict_energy(shifted,em)
                    fd(a,b)=-(ep-em)/(2*h)
                end do
            end do
            scaled_errors(k)=maxval(abs(w-fd)/tolerance)
            write(sweep_unit,'(I0,",",A,",",I0,",",I0,",",L1,3(",",ES24.16))') &
                case_number,trim(model%networks(1)%descriptor_name),accelnet_get_chebyshev_evaluation(), &
                s%natoms,s%pbc,h,maxval(abs(w-fd)),scaled_errors(k)
            ! Retain the original single-step regression bound too.
            if (h == 1e-6_real64) &
                call require(maxval(abs(w-fd)) < 2.0e-7_real64*scale, 'virial strain finite difference')
        end do
        call require(any(scaled_errors(:size(steps)-1) <= 1 .and. scaled_errors(2:) <= 1), &
            'virial step-size convergence in all nine components')
        h=1e-6_real64
        do a=1,3
            do i=1,s%natoms
                shifted=s
                shifted%positions(a,i)=s%positions(a,i)+h
                call model%predict_energy(shifted,ep)
                shifted%positions(a,i)=s%positions(a,i)-h
                call model%predict_energy(shifted,em)
                force_fd=-(ep-em)/(2*h)
                call require(abs(f(a,i)-force_fd) < 2.0e-7_real64*scale, 'force finite difference')
            end do
        end do
        shifted=s
        shifted%positions=s%positions + spread([1.7_real64,-0.6_real64,2.3_real64],2,s%natoms)
        ws=100.0_real64
        call model%predict_energy_forces(shifted,ea,fs,ws)
        if (abs(e-ea) >= 1.0e-10_real64*scale) then
            print *, 'Translation case: ', trim(model%networks(1)%descriptor_name), s%pbc, s%natoms
            print *, 'E before/after, scale: ', e, ea, scale
            print *, 'Cell: ', s%lattice
        end if
        call require(abs(e-ea) < 1.0e-10_real64*scale, 'translated energy')
        call require(maxval(abs(f-fs)) < 1.0e-10_real64*scale, 'translated forces')
        call require(maxval(abs(w-ws)) < 1.0e-10_real64*scale, 'translated virial / output reset')
        if (s%pbc) then
            open(newunit=unit, file='virial-test.xsf', status='replace', action='write')
            write(unit,'(A)') 'CRYSTAL', 'PRIMVEC'
            do i=1,3
                write(unit,*) s%lattice(:,i)
            end do
            write(unit,'(A)') 'PRIMCOORD'
            write(unit,*) s%natoms, 1
            do i=1,s%natoms
                write(unit,*) 'H ', s%positions(:,i)
            end do
            close(unit)
            call model%predict_energy_forces('virial-test.xsf', ea, file_forces, ws)
            call require(abs(e-ea) < 1.0e-10_real64*scale, 'file overload energy')
            call require(maxval(abs(f-file_forces)) < 1.0e-10_real64*scale, 'file overload forces')
            call require(maxval(abs(w-ws)) < 1.0e-10_real64*scale, 'file overload virial')
        end if
        if (s%natoms == 1) then
            call require(maxval(abs(f)) < 1.0e-10_real64*scale, 'one-atom force')
            if (s%pbc) then
                call require(maxval(abs(w)) > 1.0e-6_real64, 'self-image virial is nonzero')
            else
                call require(maxval(abs(w)) == 0.0_real64, 'empty neighborhood virial')
            end if
        end if
    end subroutine

    subroutine atomic_sum(s,e,f,w)
        type(atomic_structure), intent(in) :: s
        real(real64), intent(out) :: e
        real(real64), intent(inout) :: f(3,s%natoms), w(3,3)
        type(neighbor_data) :: neighbors
        real(real64) :: ei, old_e, old_f(3,s%natoms)
        integer :: i, first, last, n, status
        call build_neighbor_list(s, model%maximum_cutoff, neighbors, model%minimum_distance)
        old_f=f
        e=0.0_real64
        do i=1,s%natoms
            first=neighbors%offsets(i); last=neighbors%offsets(i+1)-1; n=last-first+1
            call accelnet_atomic_energy_and_forces(s%positions(:,i),s%species(i),i,n, &
                neighbors%positions(:,first:last),s%species(neighbors%atom_indices(first:last)), &
                neighbors%atom_indices(first:last),s%natoms,old_e,old_f,status)
            call require(status == ACCELNET_OK, 'old atomic force call')
            call accelnet_atomic_energy_and_forces_virial(s%positions(:,i),s%species(i),i,n, &
                neighbors%positions(:,first:last),s%species(neighbors%atom_indices(first:last)), &
                neighbors%atom_indices(first:last),s%natoms,ei,f,w,status)
            call require(status == ACCELNET_OK, 'atomic virial call')
            call require(abs(ei-old_e) < 1.0e-12_real64, 'atomic energy compatibility')
            e=e+ei
        end do
        call require(maxval(abs(f-old_f)) < 1.0e-12_real64, 'atomic force compatibility')
    end subroutine

    subroutine check_errors(uninitialized)
        logical, optional, intent(in) :: uninitialized
        real(real64) :: f(3,1), w(3,3), energy, center(3), coo(3,0)
        integer :: ids(0), status, expected, index
        expected=ACCELNET_ERR_ARGUMENT
        index=0
        if (present(uninitialized)) then
            expected=ACCELNET_ERR_INIT
            index=1
        end if
        f=2.0_real64; w=3.0_real64; center=0.0_real64
        call accelnet_atomic_energy_and_forces_virial(center,1,index,0,coo,ids,ids,1,energy,f,w,status)
        call require(status == expected, 'invalid/uninitialized status')
        call require(all(f == 2.0_real64) .and. all(w == 3.0_real64), 'error preserves accumulators')
        call require(energy == 0.0_real64, 'error energy initialized')
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
