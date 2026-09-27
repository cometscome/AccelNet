program test_batch_target
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet_batch, only: batch_workspace, evaluate_batch => evaluate_batch_reference
    use accelnet_batch, only: evaluate_batch_default => evaluate_batch
    use accelnet_batch_target, only: target_model, target_workspace, evaluate_batch_target
    use accelnet_predictor, only: predictor_model, load_predictor_from_networks, load_predictor_from_n2p2
    use accelnet_descriptors, only: atomic_structure, neighbor_data, build_neighbor_list, read_xsf, &
        descriptor_config, initialize_config
    use accelnet_descriptor_models, only: descriptor_model, add_chebyshev, add_behler
    use accelnet_behler, only: behler_config, initialize_behler_config, add_g1, add_g2, add_g3, add_g4, add_g5
    use batch_test_support
    implicit none
    type(predictor_model) :: model, source_model
    type(target_model) :: packed, snapshot
    type(target_workspace) :: work, fresh_work
    type(batch_workspace) :: cpuwork, defaultwork
    type(atomic_structure) :: s
    type(descriptor_config) :: config
    type(behler_config) :: grouped_config
    character(len=1024) :: argument, directory, filenames(2)
    integer :: v, version, mode, geometry, cutoff, activation, i, n, l, nw, checks, gpu_mode
    logical :: host
    character(len=16), parameter :: families(5) = [character(len=16) :: 'lj','behler','g4','g5','lj-behler']
    real(real64) :: max_error, degree_test, lambda_test
    checks = 0; max_error = 0
    gpu_mode = 0; host = .false.
    call get_command_argument(1,argument)
    if (argument == '--host') then
        host = .true.
        call get_command_argument(2,argument)
    end if
    if (argument == '--direct' .or. argument == '--moment') then
        gpu_mode = merge(1,2,argument == '--direct')
        argument = ''
    end if
    if (argument == '--switch-device') then
        call make_model('chebyshev',model,order=4)
        call make_structure(8,2,s)
        call packed%initialize(model,device=0,mode=2,use_host=host)
        call check(s)
        call packed%initialize(model,device=1,mode=2,use_host=host)
        call check(s)
        call packed%initialize(model,device=0,mode=2,use_host=host)
        call check(s)
    else if (argument == '--quick') then
        ! Small, nontrivial case for device memory/race checking tools.
        call make_model('chebyshev',model,order=4)
        call packed%initialize(model,mode=gpu_mode,use_host=host)
        call make_structure(8,2,s)
        call check(s)
        call finite_differences(s)
    else if (argument == '--g5-moments') then
        ! Compare shared moments with the independent CPU evaluator across
        ! cutoffs, mixed integer/fractional powers, group reuse and sparse rows.
        do v = 1,4
            select case(v)
            case(1); argument = 'g5'
            case(2); argument = 'g5-series'
            case(3); argument = 'behler'
            case(4); argument = 'lj-behler'
            end select
            do cutoff = 0,9
                call make_model(trim(argument),model,order=8)
                do i = 1,2
                    model%setups(i)%model%behler(1)%config%cutoff_type = cutoff
                    model%setups(i)%model%behler(1)%config%cutoff_alpha = 0.2_real64
                end do
                do mode = 0,3
                    call model%set_g5_evaluation(mode)
                    call packed%initialize(model,use_host=host)
                    call make_structure(8,2,s)
                    s%positions = 1.03_real64*s%positions; s%lattice = 1.03_real64*s%lattice
                    s%lattice(1,2) = 0.4_real64
                    call check(s)
                    if (mode == 3 .and. cutoff /= 0) call finite_differences(s)
                end do
            end do
            call make_structure(1,2,s)
            call check(s) ! periodic self images
            s%pbc = .false.
            call check(s) ! no neighbors
        end do
        do v = 1,5
            call make_model('g5',model)
            do i = 1,2
                model%setups(i)%model = descriptor_model()
                ! Separate components may request different modes while sharing
                ! a radial cache; only eligible exact powers enter moments.
                do n = 1,4
                    call initialize_behler_config(grouped_config,2,cutoff_type=1)
                    degree_test = real(n,real64)
                    if (n == 2) degree_test = 10
                    if (n == 3) degree_test = 1.5_real64
                    if (n == 4) degree_test = 11
                    if (v == 2 .and. n == 2) degree_test = 2+4e-13_real64
                    if (v == 3 .and. n == 4) degree_test = 16
                    lambda_test = -1
                    if (v == 4) lambda_test = 0.7_real64
                    if (v == 5) lambda_test = 0
                    call add_g5(grouped_config,1,merge(1,2,mod(n,2) == 0), &
                        3.4_real64,lambda_test,degree_test,merge(0.2_real64,0.3_real64,n < 3),0.4_real64)
                    grouped_config%g5_evaluation_mode = 3
                    if (n == 1 .and. i == 1) grouped_config%g5_evaluation_mode = 1
                    call add_behler(model%setups(i)%model,grouped_config)
                end do
            end do
            call packed%initialize(model,use_host=host)
            call make_structure(8,2,s)
            call check(s)
            call finite_differences(s)
            call make_structure(4,2,s)
            s%pbc = .false.; s%positions = 0
            s%positions(1,:) = [0.0_real64,0.9_real64,-1.2_real64,1.8_real64]
            call check(s) ! zero direction components and both angular endpoints
            call finite_differences(s)
        end do
        ! Mixed cutoffs/components: auto counts neighbors inside each component's
        ! angular cutoff, not the full CSR list or a different component's Rc.
        call make_model('g5',model)
        do i = 1,2
            model%setups(i)%model = descriptor_model()
            do n = 1,2
                call initialize_behler_config(grouped_config,2)
                call add_g5(grouped_config,1,1,real(n+1,real64),1.0_real64,2.0_real64,0.2_real64)
                call add_g5(grouped_config,1,2,real(n+1,real64),-1.0_real64,4.0_real64,0.2_real64)
                call add_behler(model%setups(i)%model,grouped_config)
            end do
        end do
        call make_structure(64,2,s)
        do mode = 0,3
            call model%set_g5_evaluation(mode)
            call packed%initialize(model,use_host=host)
            call check(s)
        end do
        ! Mixed central-element algorithms and reload of group/basis capacities.
        call make_model('g5-series',model,order=4)
        do mode = 1,3
            model%setups(1)%model%behler(1)%config%g5_evaluation_mode = mode
            model%setups(2)%model%behler(1)%config%g5_evaluation_mode = 4-mode
            call packed%initialize(model,use_host=host)
            call make_structure(8,2,s)
            call check(s)
        end do
    else if (argument == '--descriptors') then
        do v = 1,size(families)
            do cutoff = 0,9
                print *, "DESCRIPTOR CASE ",trim(families(v)),cutoff
                call make_model(trim(families(v)),model)
                do i = 1,2
                    if (allocated(model%setups(i)%model%lj)) then
                        model%setups(i)%model%lj(1)%config%cutoff_type = cutoff
                        model%setups(i)%model%lj(1)%config%cutoff_alpha = 0.2_real64
                    end if
                    if (allocated(model%setups(i)%model%behler)) then
                        model%setups(i)%model%behler(1)%config%cutoff_type = cutoff
                        model%setups(i)%model%behler(1)%config%cutoff_alpha = 0.2_real64
                    end if
                end do
                call packed%initialize(model,use_host=host)
                call make_structure(8,2,s)
                s%positions = 1.03_real64*s%positions; s%lattice = 1.03_real64*s%lattice
                s%lattice(1,2) = 0.4_real64
                call check(s)
                if (cutoff /= 0) call finite_differences(s)
                if (trim(families(v)) == 'g4') then
                    ! Exercise the squared j-k cutoff test just inside, exactly
                    ! on, and just outside Rc. Do not finite-difference a jump.
                    do geometry = -1,1
                        call make_structure(3,2,s)
                        s%pbc = .false.; s%positions = 0; s%species = 1
                        s%positions(1,:) = [0.0_real64,1.2_real64,-2.2_real64+geometry*1e-10_real64]
                        call check(s)
                    end do
                end if
            end do
            call make_structure(64,2,s)
            call check(s)
            call make_structure(1,2,s)
            call check(s)
            s%pbc = .false.
            call check(s)
        end do
        ! Change group sharing without changing model dimensions: all G4/G5
        ! share one radial group, split one feature, then merge it again. This
        ! exercises cache growth/reuse and independence from lambda/species.
        call work%release()
        call make_model('behler',model)
        call make_structure(8,2,s)
        do v = 1,3
            do i = 1,2
                call initialize_behler_config(grouped_config,2,cutoff_type=1)
                do n = 1,2
                    call add_g1(grouped_config,n,3.2_real64)
                    call add_g2(grouped_config,n,3.4_real64,0.3_real64,0.6_real64)
                    call add_g3(grouped_config,n,3.1_real64,1.3_real64)
                    call add_g4(grouped_config,n,n,3.4_real64,-1.0_real64,2.0_real64,0.2_real64,0.4_real64)
                    call add_g4(grouped_config,1,2,3.4_real64,1.0_real64,1.5_real64,0.2_real64,0.4_real64)
                    call add_g5(grouped_config,n,n,3.4_real64,1.0_real64,real(n,real64),0.2_real64,0.4_real64)
                    call add_g5(grouped_config,1,2,3.4_real64,-1.0_real64,2.5_real64, &
                        merge(0.17_real64,0.2_real64,v == 2 .and. n == 2),0.4_real64)
                end do
                model%setups(i)%model%behler(1)%config = grouped_config
            end do
            call packed%initialize(model,use_host=host)
            nw = work%allocations()
            call check(s)
            if (v == 2) call require(work%allocations() == nw+1,'radial cache group growth')
            if (v == 3) call require(work%allocations() == nw,'radial cache group reuse')
            call finite_differences(s)
        end do
        ! Sparse polynomial groups, changed coefficient capacity, fractional
        ! and high-degree fallbacks, reversed species pairs, and lambda changes.
        do geometry = 1,2
            call make_model(trim(families(geometry+2)),model)
            call work%release()
            do v = 1,6
                do i = 1,2
                    call initialize_behler_config(grouped_config,2,cutoff_type=1)
                    do n = 1,4
                        degree_test = real(2**(n-1),real64)
                        if (v >= 2 .and. n == 4) degree_test = 16
                        if (v == 4 .and. mod(n,2) == 0) degree_test = 1.5_real64
                        if (v == 5 .and. mod(n,2) == 0) degree_test = 17
                        ! Near-integer powers retain the original derivative prefactor.
                        if (v == 6 .and. mod(n,2) == 0) degree_test = 2+real(n,real64)*1e-13_real64
                        lambda_test = -1
                        if (v == 3 .and. n == 4) lambda_test = 0.7_real64
                        if (geometry == 1) then
                            call add_g4(grouped_config,1+mod(n,2),2-mod(n,2),3.4_real64, &
                                lambda_test,degree_test,0.2_real64,0.4_real64)
                        else
                            call add_g5(grouped_config,1+mod(n,2),2-mod(n,2),3.4_real64, &
                                lambda_test,degree_test,0.2_real64,0.4_real64)
                        end if
                    end do
                    model%setups(i)%model%behler(1)%config = grouped_config
                end do
                call packed%initialize(model,use_host=host)
                call make_structure(8,2,s)
                nw = work%allocations()
                call check(s)
                if (v == 2) call require(work%allocations() == nw,'angular coefficient storage reuse')
                call finite_differences(s)
                ! Collinear neighbors exercise t=0 and t=1, including zeta=1.
                call make_structure(4,2,s)
                s%pbc = .false.; s%positions = 0
                s%positions(1,:) = [0.0_real64,0.9_real64,-1.2_real64,1.8_real64]
                call check(s)
                call finite_differences(s)
            end do
        end do
        ! Four distinct G4 inputs: the second mixed-species term changes eta,
        ! so it must not be broadcast from the first term's angular group.
        call make_model('g4-distinct',model)
        call packed%initialize(model,use_host=host)
        call make_structure(8,2,s)
        call check(s)
        call finite_differences(s)
        ! G4 Jacobian storage must grow when G4 replaces a same-sized G5 model,
        ! then remain reusable when switching away from and back to G4.
        call work%release()
        call make_structure(8,2,s)
        do v = 1,4
            if (mod(v,2) == 1) then
                call make_model('g5',model)
            else
                call make_model('g4',model)
            end if
            call packed%initialize(model,use_host=host)
            nw = work%allocations()
            call check(s)
            if (v == 2) call require(work%allocations() == nw+1,'G4 Jacobian growth')
            if (v >= 3) call require(work%allocations() == nw,'G4 Jacobian reuse')
        end do
        ! Independent cutoff/exponential/angular reuse across G4 components.
        ! An angular representative with a different Rc must still initialize
        ! its follower's power; cutoff type/alpha must never be conflated.
        call make_model('g4-distinct',model)
        do i = 1,2
            model%setups(i)%model = descriptor_model()
            call initialize_behler_config(grouped_config,2,cutoff_type=1,cutoff_alpha=0.2_real64)
            call add_g4(grouped_config,1,2,2.8_real64,1.0_real64,2.0_real64,0.2_real64,0.4_real64)
            call add_g4(grouped_config,2,1,3.4_real64,1.0_real64,2.0_real64,0.2_real64,0.4_real64)
            call add_behler(model%setups(i)%model,grouped_config)
            call initialize_behler_config(grouped_config,2,cutoff_type=4,cutoff_alpha=0.1_real64)
            call add_g4(grouped_config,1,2,3.4_real64,1.0_real64,1.5_real64,0.3_real64,0.4_real64)
            call add_g4(grouped_config,1,1,3.4_real64,1.0_real64,2.0_real64,0.3_real64,0.2_real64)
            call add_behler(model%setups(i)%model,grouped_config)
        end do
        call packed%initialize(model,use_host=host)
        call check(s)
        call finite_differences(s)
        ! Interleaved radial and angular groups must refresh scalar pair caches
        ! when a previously visited group is encountered again.
        do i = 1,2
            model%setups(i)%model = descriptor_model()
            call initialize_behler_config(grouped_config,2)
            call add_g4(grouped_config,1,2,3.4_real64,1.0_real64,2.0_real64,0.2_real64,0.4_real64)
            call add_g4(grouped_config,2,1,3.4_real64,-1.0_real64,3.0_real64,0.3_real64,0.1_real64)
            call add_g4(grouped_config,1,2,3.4_real64,1.0_real64,1.5_real64,0.2_real64,0.4_real64)
            call add_g4(grouped_config,2,1,3.4_real64,1.0_real64,2.0_real64,0.3_real64,0.1_real64)
            call add_behler(model%setups(i)%model,grouped_config)
        end do
        call packed%initialize(model,use_host=host)
        call check(s)
        call finite_differences(s)
        ! More G4 columns than the GPU owner count: each owner must retain all
        ! its columns across every pair without races or skipped derivatives.
        call make_model('g4-series',model,order=12)
        call packed%initialize(model,use_host=host)
        call check(s)
        call finite_differences(s)
        ! Mixed species families share NN and scatter kernels, but descriptors
        ! take different paths. Reload metadata and the local species map.
        call make_model('chebyshev',source_model,order=4)
        call make_model('behler',model)
        model%setups(1) = source_model%setups(1); model%networks(1) = source_model%networks(1)
        call packed%initialize(model,use_host=host)
        call make_structure(8,2,s)
        call check(s)
        call finite_differences(s)
        call make_model('behler',model)
        call packed%initialize(model,use_host=host)
        call check(s)
        model%setups(1)%global_to_local = [2,1]
        call packed%initialize(model,use_host=host)
        i = work%uploads()
        call check(s)
        call require(work%uploads() == i+1,'local species map reload')
        model%setups(1)%model%behler(1)%config%cutoff_alpha = 0.15_real64
        call packed%initialize(model,use_host=host)
        i = work%uploads()
        call check(s)
        call require(work%uploads() == i+1,'descriptor parameter reload')
        ! An exactly representable hard-cutoff boundary and points on both
        ! sides; finite differences of a discontinuous cutoff are undefined.
        call make_model('lj',model)
        do i = 1,2
            model%setups(i)%model%lj(1)%config%cutoff_type = 0
            model%setups(i)%model%lj(1)%config%radial_rc = 2.0_real64
        end do
        call packed%initialize(model,use_host=host)
        call make_structure(2,2,s)
        s%pbc = .false.; s%positions = 0
        do i = -1,1
            s%positions(1,2) = 2.0_real64+real(i,real64)*1e-6_real64
            call check(s)
        end do
    else if (argument == '--n2p2') then
        call get_command_argument(2,directory)
        if (host) call get_command_argument(3,directory)
        call load_predictor_from_n2p2(trim(directory),model)
        call packed%initialize(model,use_host=host)
        do geometry = 1,3
            call make_structure(8,size(model%networks),s)
            s%pbc = geometry /= 1
            if (geometry == 3) s%lattice(1,2) = 0.4_real64
            call check(s)
            call finite_differences(s)
        end do
    else if (argument == '--embedded') then
        call get_command_argument(2,filenames(1)); call get_command_argument(3,filenames(2))
        call load_predictor_from_networks(filenames,model)
        call packed%initialize(model,use_host=host)
        call make_structure(8,2,s)
        call check(s)
        call finite_differences(s)
    else if (argument == '--aenet') then
        call get_command_argument(2,directory)
        call get_command_argument(3,argument)
        if (len_trim(argument) > 0) read(argument,*) gpu_mode
        filenames(1) = trim(directory)//'/Ti.nn.ascii'; filenames(2) = trim(directory)//'/O.nn.ascii'
        do v = 0, 2
            version = merge(10,v,v == 2)
            call load_predictor_from_networks(filenames,model,chebyshev_version=version)
            call packed%initialize(model,mode=gpu_mode,use_host=host)
            do geometry = 1, 2
                if (geometry == 1) then
                    call read_xsf(trim(directory)//'/structure0001.xsf',model%species_names,s)
                else
                    call read_xsf(trim(directory)//'/structure2935.xsf',model%species_names,s)
                    s%lattice(1,:) = s%lattice(1,:)+0.04_real64*s%lattice(2,:)
                    s%positions(1,:) = s%positions(1,:)+0.04_real64*s%positions(2,:)
                end if
                do mode = 0, 2
                    call model%set_chebyshev_evaluation(mode)
                    call check(s)
                end do
                if (v == 0) call finite_differences(s)
            end do
        end do
    else if (len_trim(argument) > 0) then
        call reject(trim(argument))
        stop 0
    else
        do v = 0, 2
            version = merge(10,v,v == 2)
            call make_model('chebyshev',model,order=3+v,version=version)
            call packed%initialize(model,mode=gpu_mode,use_host=host)
            do mode = 0, 2
                call model%set_chebyshev_evaluation(mode)
                do geometry = 1, 3
                    call make_structure(8,2,s)
                    s%pbc = geometry /= 1
                    if (geometry == 3) s%lattice(1,2) = 0.4_real64
                    call check(s)
                end do
            end do
            call finite_differences(s)
        end do
        do cutoff = 0, 9
            call make_model('chebyshev',model)
            do i = 1, 2
                model%setups(i)%model%chebyshev(1)%config%cutoff_type = cutoff
                model%setups(i)%model%chebyshev(1)%config%cutoff_alpha = 0.2_real64
                model%setups(i)%model%chebyshev(1)%config%angular_rc = 2.9_real64
            end do
            call packed%initialize(model,mode=gpu_mode,use_host=host)
            call make_structure(8,2,s)
            s%pbc = .false.
            call check(s)
            if (cutoff /= 0) call finite_differences(s)
        end do
        ! Different element depths/widths, plus every supported activation in
        ! hidden AND output layers. Weight layout is independently constructed.
        do activation = 0, 11
            call make_model('chebyshev',model,order=2)
            model%networks(2)%nlayers = 3
            model%networks(2)%nodes = [12,5,1]
            model%networks(2)%activation = [activation,activation]
            model%networks(2)%weight_offsets = [0,65,71]
            model%networks(2)%weights = [(0.04_real64*sin(real(i,real64)),i=1,71)]
            model%networks(1)%activation = activation
            call packed%initialize(model,mode=gpu_mode,use_host=host)
            call check(s)
            call finite_differences(s)
        end do
        call make_model('chebyshev',model,order=8)
        call packed%initialize(model,mode=gpu_mode,use_host=host)
        call make_structure(64,2,s)
        call check(s)
        call make_structure(512,2,s)
        call check(s)
        call make_structure(1,2,s)
        call check(s) ! periodic self images, force cancellation, nonzero virial
        s%pbc = .false.
        call check(s) ! isolated atom, zero edges
        ! Clenshaw/Horner endpoints and their empty/one-step recurrences. Collinear
        ! neighbors exercise cos(theta)=+1 and -1, including version 1's
        ! affine angular argument outside [-1,1].
        do v = 0, 1
            do n = 0, 12
                if (n /= 0 .and. n /= 1 .and. n /= 6 .and. n /= 12) cycle
                call make_model('chebyshev',model,order=n,version=v)
                call make_structure(4,2,s)
                s%pbc = .false.; s%positions = 0
                s%positions(1,:) = [-1.4_real64,0.0_real64,1.2_real64,2.5_real64]
                do mode = 1, 2
                    call packed%initialize(model,mode=mode,use_host=host)
                    call check(s)
                    call finite_differences(s)
                end do
            end do
        end do
        ! Different polynomial degrees share one padded moment buffer. Horner
        ! must traverse each species' own compact coefficient extent.
        call make_model('chebyshev',model,order=1)
        call make_model('chebyshev',source_model,order=7)
        model%setups(2) = source_model%setups(2)
        model%networks(2) = source_model%networks(2)
        call packed%initialize(model,mode=2,use_host=host)
        call make_structure(8,2,s)
        call check(s)
        call finite_differences(s)
        ! Force coefficients reuse moment scratch only for moment rows. Exercise
        ! direct/moment species together, then swap them in the same workspace.
        call make_model('chebyshev',model,order=4)
        call make_structure(64,2,s)
        do mode = 1, 2
            model%setups(1)%model%chebyshev(1)%config%evaluation_mode = mode
            model%setups(2)%model%chebyshev(1)%config%evaluation_mode = 3-mode
            call packed%initialize(model,use_host=host)
            call check(s)
            call finite_differences(s)
        end do
        call work%release()
        call require(work%allocations() == 0,'workspace release')
        call check(s)
        ! Single-species order-zero model; no weighted descriptor channels.
        call make_model('chebyshev',model,order=0)
        ! Avoid self-slicing allocatable derived-type components: NVHPC 25.3
        ! retains their old outer extent under -fast in this fixture.
        source_model = model
        deallocate(model%networks,model%setups,model%species_names)
        allocate(model%networks(1),source=source_model%networks(:1))
        allocate(model%setups(1),source=source_model%setups(:1))
        allocate(model%species_names(1),source=source_model%species_names(:1))
        deallocate(model%setups(1)%model%chebyshev)
        model%setups(1)%model%num_outputs = 0
        call initialize_config(config,1,3.4_real64,0,3.4_real64,0)
        call add_chebyshev(model%setups(1)%model,config)
        model%setups(1)%global_to_local = [1]
        model%networks(1)%nodes = [2,8,4,1]
        model%networks(1)%descriptor_shift = [0.03_real64,0.06_real64]
        model%networks(1)%descriptor_scale = [0.21_real64,0.22_real64]
        model%networks(1)%atomic_references = [0.3_real64]
        nw = 0
        do l = 1, 3
            model%networks(1)%weight_offsets(l) = nw
            nw = nw+(model%networks(1)%nodes(l)+1)*model%networks(1)%nodes(l+1)
        end do
        model%networks(1)%weight_offsets(4) = nw
        model%networks(1)%weights = [(0.04_real64*sin(real(i,real64)),i=1,nw)]
        call packed%initialize(model,mode=gpu_mode,use_host=host)
        call make_structure(8,1,s)
        call check(s)
        call finite_differences(s)
        ! Packed models own a snapshot: copies and release must not alias data.
        snapshot = packed
        call packed%release()
        packed = snapshot
        call snapshot%release()
        call check(s)
        ! A copied scratch workspace must never share ownership of mappings.
        fresh_work = work
        call require(fresh_work%allocations() == 0 .and. fresh_work%uploads() == 0,'workspace copy is a fresh cache')
        call fresh_work%release()
        call check(s)
        ! Same-shape reload: changed weights/scaling must trigger a fresh upload.
        model%networks(1)%weights = model%networks(1)%weights+0.007_real64
        model%networks(1)%descriptor_scale = model%networks(1)%descriptor_scale*0.9_real64
        call packed%initialize(model,mode=gpu_mode,use_host=host)
        i = work%uploads()
        call check(s)
        call require(work%uploads() == i+1,'model reload must replace resident parameters')
        block
            type(target_workspace) :: scoped
            type(neighbor_data) :: nb
            real(real64) :: e(8), f(3,8)
            call build_neighbor_list(s,model%maximum_cutoff,nb,model%minimum_distance)
            f = 0
            call evaluate_batch_target(packed,s%species,[(i,i=1,8)],nb%offsets,nb%atom_indices,nb%displacements,e,f,scoped)
        end block ! finalizer must delete mappings before host allocations vanish
        call check(s)
    end if
    call work%release()
    print *, 'GPU checks:',checks,' maximum absolute E/F/W error:',max_error
    print *, 'GPU batch FP64 equivalence and finite differences passed'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not. ok) then
            print *, message
            error stop 1
        end if
    end subroutine

    subroutine close_array(actual,expected,message)
        real(real64), intent(in) :: actual(:), expected(:)
        character(len=*), intent(in) :: message
        call require(all(ieee_is_finite(actual)),message//' finite')
        max_error = max(max_error,maxval(abs(actual-expected)))
        if (.not. all(abs(actual-expected) <= 2e-10_real64+2e-10_real64*abs(expected))) then
            print *, 'error:',maxval(abs(actual-expected)), 'version/cutoff/activation:',version,cutoff,activation
            call require(.false.,message)
        end if
    end subroutine

    subroutine check(s)
        type(atomic_structure), intent(in) :: s
        type(neighbor_data) :: nb
        real(real64) :: e(s%natoms), eg(s%natoms), split_e(s%natoms), f(3,s%natoms), fg(3,s%natoms), w(3,3), wg(3,3)
        integer :: centers(s%natoms), rows(s%natoms), split, count, j, uploads
        centers = [(j,j=1,s%natoms)]
        call build_neighbor_list(s,model%maximum_cutoff,nb,model%minimum_distance)
        f = 0; w = 0
        call evaluate_batch(model,s%species,centers,nb%offsets,nb%atom_indices,nb%displacements,e,f,cpuwork,w)
        ! Exercise the public CPU dispatch against the independent reference
        ! for every descriptor/cutoff/network case, even in GPU-enabled builds.
        fg = 0; wg = 0
        call evaluate_batch_default(model,s%species,centers,nb%offsets,nb%atom_indices, &
            nb%displacements,eg,fg,defaultwork,wg)
        call close_array(eg,e,'default CPU atomic energies')
        call close_array(reshape(fg,[3*s%natoms]),reshape(f,[3*s%natoms]),'default CPU forces')
        call close_array(reshape(wg,[9]),reshape(w,[9]),'default CPU virial')
        fg = 0; wg = 0
        call evaluate_batch_target(packed,s%species,centers,nb%offsets,nb%atom_indices,nb%displacements,eg,fg,work,wg)
        call close_array(eg,e,'atomic energies')
        call close_array(reshape(fg,[3*s%natoms]),reshape(f,[3*s%natoms]),'forces')
        call close_array(reshape(wg,[9]),reshape(w,[9]),'virial')
        count = work%allocations(); uploads = work%uploads(); split = s%natoms/2
        fg = 0.25_real64; wg = 0.5_real64
        call evaluate_batch_target(packed,s%species,centers(:split),nb%offsets(:split+1),nb%atom_indices, &
            nb%displacements,split_e(:split),fg,work,wg)
        call evaluate_batch_target(packed,s%species,centers(split+1:),nb%offsets(split+1:),nb%atom_indices, &
            nb%displacements,split_e(split+1:),fg,work,wg)
        call close_array(split_e,e,'partition atomic energies')
        call close_array(reshape(fg-0.25_real64,[3*s%natoms]),reshape(f,[3*s%natoms]),'partition additive forces')
        call close_array(reshape(wg-0.5_real64,[9]),reshape(w,[9]),'partition additive virial')
        call require(work%allocations() == count,'workspace reuse')
        call require(work%uploads() == uploads,'unchanged model must remain resident')
        ! Force-only call and reversed, noncontiguous central rows. Each CSR row
        ! keeps its own neighbors, including ghost and periodic image targets.
        fg = 0
        do j = s%natoms, 1, -1
            call evaluate_batch_target(packed,s%species,centers(j:j),nb%offsets(j:j+1),nb%atom_indices, &
                nb%displacements,split_e(j:j),fg,work)
        end do
        call close_array(split_e,e,'reordered energies')
        call close_array(reshape(fg,[3*s%natoms]),reshape(f,[3*s%natoms]),'reordered forces')
        call evaluate_batch_target(packed,[integer::],[integer::],[1],[integer::], &
            nb%displacements(:,:0),eg(:0),fg(:,:0),work,wg)
        checks = checks+1
    end subroutine

    subroutine gpu_eval(s,e,f,w)
        type(atomic_structure), intent(in) :: s
        real(real64), intent(out) :: e, f(3,s%natoms), w(3,3)
        type(neighbor_data) :: nb
        real(real64) :: energies(s%natoms)
        integer :: j
        call build_neighbor_list(s,model%maximum_cutoff,nb,model%minimum_distance)
        f = 0; w = 0
        call evaluate_batch_target(packed,s%species,[(j,j=1,s%natoms)],nb%offsets,nb%atom_indices, &
            nb%displacements,energies,f,work,w)
        e = sum(energies)
    end subroutine

    subroutine finite_differences(s)
        type(atomic_structure), intent(in) :: s
        type(atomic_structure) :: shifted
        real(real64) :: e, ep, em, f(3,s%natoms), dummy_f(3,s%natoms), w(3,3), dummy_w(3,3), fd
        real(real64), parameter :: h = 2e-5_real64
        integer :: a, b
        call gpu_eval(s,e,f,w)
        do a = 1, 3
            shifted = s; shifted%positions(a,1) = s%positions(a,1)+h
            call gpu_eval(shifted,ep,dummy_f,dummy_w)
            shifted = s; shifted%positions(a,1) = s%positions(a,1)-h
            call gpu_eval(shifted,em,dummy_f,dummy_w)
            fd = -(ep-em)/(2*h)
            call require(abs(fd-f(a,1)) < 3e-6_real64*max(1.0_real64,abs(f(a,1))),'GPU force finite difference')
        end do
        ! dr(a)*F(b) convention: differentiate deformation of coordinate b by a.
        do a = 1, 3
            do b = 1, 3
                shifted = s
                shifted%positions(b,:) = s%positions(b,:)+h*s%positions(a,:)
                shifted%lattice(b,:) = s%lattice(b,:)+h*s%lattice(a,:)
                call gpu_eval(shifted,ep,dummy_f,dummy_w)
                shifted = s
                shifted%positions(b,:) = s%positions(b,:)-h*s%positions(a,:)
                shifted%lattice(b,:) = s%lattice(b,:)-h*s%lattice(a,:)
                call gpu_eval(shifted,em,dummy_f,dummy_w)
                fd = -(ep-em)/(2*h)
                call require(abs(fd-w(a,b)) < 1e-5_real64*max(1.0_real64,abs(w(a,b))),'GPU virial finite difference')
            end do
        end do
    end subroutine

    subroutine reject(which)
        character(len=*), intent(in) :: which
        integer :: species(2), centers(2), offsets(3), indices(2)
        real(real64) :: dr(3,2), e(2), f(3,2)
        call make_model('chebyshev',model)
        if (which == 'unsupported') call make_model('combined',model)
        if (which == 'activation') model%networks(1)%activation(1) = 12
        if (which == 'version10') then
            do i = 1, 2
                model%setups(i)%model%chebyshev(1)%config%version = 10
            end do
        end if
        if (which /= 'uninitialized') call packed%initialize(model,mode=gpu_mode,use_host=host)
        species = [1,2]; centers = [1,2]; offsets = [1,2,3]; indices = [2,1]
        dr(:,1) = [1.0_real64,0.0_real64,0.0_real64]; dr(:,2) = -dr(:,1); f = 0
        select case(which)
        case('offset'); offsets(1) = 0
        case('monotonic'); offsets = [1,3,2]
        case('target'); indices(1) = 3
        case('center'); centers(1) = 0
        case('species'); species(2) = 3
        case('shape')
            call evaluate_batch_target(packed,species,centers,offsets(:2),indices,dr,e,f,work)
            return
        end select
        call evaluate_batch_target(packed,species,centers,offsets,indices,dr,e,f,work)
    end subroutine
end program
