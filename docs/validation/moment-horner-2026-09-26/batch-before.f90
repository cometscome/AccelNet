! Shared packed batch backend. CMake also compiles a directive-free serial
! module instance for the CPU batch API, without an OpenMP runtime requirement.
module accelnet_batch_target
    use iso_fortran_env, only: real64
    use iso_c_binding, only: c_ptr, c_f_pointer
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet_target_runtime, only: omp_get_initial_device, omp_get_num_devices, omp_get_default_device, &
        omp_is_initial_device, omp_get_wtime
    use accelnet_predictor, only: predictor_model
    use accelnet_descriptors, only: descriptor_config, validate_cutoff_parameters
    use aenet_network, only: atomic_network
    use accelnet_target_descriptors, only: pack_generic_descriptors
    use accelnet_target_kernels, only: run_target_batch
    implicit none
    private
    public :: target_model, target_workspace, target_profile, evaluate_batch_target, evaluate_lammps_target

    type :: target_profile
        real(real64) :: prepare = 0, upload = 0, descriptors = 0, network = 0, forces = 0, download = 0, total = 0
    end type

    type :: target_model
        private
        integer :: device = -1, maxnodes = 0, maxlayers = 0
        integer, allocatable :: meta(:,:), nodes(:,:), acts(:,:), woffset(:,:), local_species(:,:), mp(:,:,:)
        integer, allocatable :: features(:,:,:)
        real(real64), allocatable :: feature_params(:,:,:)
        real(real64), allocatable :: multiplicity(:,:), polynomial(:,:,:)
        real(real64), allocatable :: weights(:,:), params(:,:), shift(:,:), scale(:,:), spin(:,:)
    contains
        procedure :: initialize => initialize_target_model
        procedure :: release => release_target_model
    end type

    type :: target_workspace
        private
        type(target_model) :: cached
        real(real64), allocatable :: g(:,:), values(:,:,:), deriv(:,:,:), delta(:,:,:), moments(:,:,:), powers(:,:,:)
        integer, allocatable :: species(:), centers(:), offsets(:), indices(:), use_moment(:), edge_row(:)
        real(real64), allocatable :: dr(:,:), energies(:), forces(:,:), virial(:,:), geom(:,:), edge_force(:,:)
        integer :: growth_count = 0, upload_count = 0
    contains
        procedure :: release => release_target_workspace
        procedure :: allocations => target_allocations
        procedure :: uploads => target_uploads
        procedure, private :: copy_workspace
        generic :: assignment(=) => copy_workspace
        final :: finalize_workspace, finalize_workspace_array
    end type
contains
    subroutine release_target_model(self)
        class(target_model), intent(inout) :: self
        if (allocated(self%meta)) deallocate(self%meta, self%nodes, self%acts, self%woffset, self%local_species, &
            self%features, self%feature_params, self%weights, self%params, self%shift, self%scale, &
            self%spin, self%mp, self%multiplicity, self%polynomial)
        self%device = -1; self%maxnodes = 0; self%maxlayers = 0
    end subroutine

    subroutine initialize_target_model(self, model, device, mode, status, message, use_host)
        class(target_model), intent(inout) :: self
        type(predictor_model), intent(in) :: model
        integer, intent(in), optional :: device, mode
        integer, optional, intent(out) :: status
        character(len=*), optional, intent(out) :: message
        logical, optional, intent(in) :: use_host
        type(descriptor_config) :: config
        type(atomic_network) :: net
        integer :: ns, s, d, ml, mn, mw, l, count_weights, j, local, mm, ma, nm, na
        integer, allocatable :: fi(:,:)
        real(real64), allocatable :: fr(:,:)
        integer :: descriptor_species, pack_status
        character(len=256) :: pack_message
        logical :: on_host, chebyshev, host_requested
        if (present(status)) status = 0
        if (present(message)) message = ''
        call self%release()
        if (.not. allocated(model%networks) .or. .not. allocated(model%setups)) &
            then
            call target_failure('OpenMP target: uninitialized model',status,message)
            return
        end if
        ns = size(model%networks)
        if (ns == 0 .or. size(model%setups) /= ns) then
            call target_failure('OpenMP target: invalid species count',status,message)
            return
        end if
        mn = 0; ml = 0; mw = 0; mm = 0; ma = 0
        do s = 1, ns
            chebyshev = allocated(model%setups(s)%model%chebyshev)
            net = model%networks(s)
            if (chebyshev) then
                if (size(model%setups(s)%model%chebyshev) /= 1 .or. &
                    allocated(model%setups(s)%model%lj) .or. allocated(model%setups(s)%model%behler)) then
                    call target_failure('OpenMP target: mixed or multiple Chebyshev components are not supported',status,message)
                    return
                end if
            config = model%setups(s)%model%chebyshev(1)%config
            net = model%networks(s)
            if (config%cutoff_type < 0 .or. config%cutoff_type > 9 .or. &
                config%cutoff_alpha < 0 .or. config%cutoff_alpha >= 1 .or. &
                (config%cutoff_type == 9 .and. config%cutoff_alpha <= 0)) &
                then
                call target_failure('OpenMP target: invalid cutoff parameters',status,message)
                return
            end if
            if (config%radial_rc <= 0 .or. config%angular_rc <= 0 .or. &
                config%radial_order < 0 .or. config%angular_order < 0) then
                call target_failure('OpenMP target: invalid Chebyshev config',status,message)
                return
            end if
            if (config%version /= 0 .and. config%version /= 1 .and. config%version /= 10) &
                then
                call target_failure('OpenMP target: unsupported Chebyshev version',status,message)
                return
            end if
            if (config%version == 10 .and. config%num_species > 1 .and. config%central_type_index < 1) &
                then
                call target_failure('OpenMP target: invalid version 10 center lookup',status,message)
                return
            end if
            if (.not. allocated(config%species_weights)) then
                call target_failure('OpenMP target: missing species weights',status,message)
                return
            end if
            if (size(config%species_weights) /= config%num_species) then
                call target_failure('OpenMP target: invalid species weights',status,message)
                return
            end if
                descriptor_species = config%num_species
            else
                call pack_generic_descriptors(model%setups(s)%model,fi,fr,pack_status,pack_message)
                if (pack_status /= 0) then
                    call target_failure(trim(pack_message),status,message)
                    return
                end if
                descriptor_species = minval(fi(6,:))
            end if
            if (.not. allocated(model%setups(s)%global_to_local)) then
                call target_failure('OpenMP target: missing species map',status,message)
                return
            end if
            if (size(model%setups(s)%global_to_local) /= ns) then
                call target_failure('OpenMP target: invalid species map size',status,message)
                return
            end if
            if (any(model%setups(s)%global_to_local < 0) .or. &
                any(model%setups(s)%global_to_local > descriptor_species)) then
                call target_failure('OpenMP target: invalid species map values',status,message)
                return
            end if
            if (.not. allocated(net%nodes) .or. .not. allocated(net%activation) .or. &
                .not. allocated(net%weight_offsets) .or. .not. allocated(net%weights)) &
                then
                call target_failure('OpenMP target: uninitialized network',status,message)
                return
            end if
            if (net%nlayers < 2 .or. size(net%nodes) /= net%nlayers) then
                call target_failure('OpenMP target: invalid network depth',status,message)
                return
            end if
            if (any(net%nodes < 1) .or. net%nodes(net%nlayers) /= 1) then
                call target_failure('OpenMP target: invalid network width',status,message)
                return
            end if
            d = net%nodes(1)
            if (d /= model%setups(s)%model%num_outputs) then
                call target_failure('OpenMP target: descriptor/network size mismatch',status,message)
                return
            end if
            if (size(net%activation) < net%nlayers-1 .or. size(net%weight_offsets) < net%nlayers-1) &
                then
                call target_failure('OpenMP target: invalid network metadata',status,message)
                return
            end if
            if (any(net%activation(:net%nlayers-1) < 0) .or. any(net%activation(:net%nlayers-1) > 11)) &
                then
                call target_failure('OpenMP target: unsupported network activation',status,message)
                return
            end if
            count_weights = 0
            do l = 1, net%nlayers-1
                if (net%weight_offsets(l) /= count_weights) then
                    call target_failure('OpenMP target: invalid weight offsets',status,message)
                    return
                end if
                if (net%nodes(l) >= huge(j)/net%nodes(l+1)-1) then
                    call target_failure('OpenMP target: network too large',status,message)
                    return
                end if
                j = (net%nodes(l)+1)*net%nodes(l+1)
                if (j > huge(j)-count_weights) then
                    call target_failure('OpenMP target: network too large',status,message)
                    return
                end if
                count_weights = count_weights+j
            end do
            if (size(net%weights) /= count_weights) then
                call target_failure('OpenMP target: invalid weights',status,message)
                return
            end if
            if (.not. allocated(net%descriptor_shift) .or. .not. allocated(net%descriptor_scale) .or. &
                .not. allocated(net%atomic_references)) then
                call target_failure('OpenMP target: missing normalization',status,message)
                return
            end if
            if (size(net%descriptor_shift) /= d .or. size(net%descriptor_scale) /= d .or. &
                size(net%atomic_references) /= ns) then
                call target_failure('OpenMP target: invalid normalization sizes',status,message)
                return
            end if
            if (.not. ieee_is_finite(net%energy_scale) .or. net%energy_scale == 0) &
                then
                call target_failure('OpenMP target: invalid energy scale',status,message)
                return
            end if
            if (chebyshev) then
            if (config%evaluation_mode < 0 .or. config%evaluation_mode > 2) then
                call target_failure('OpenMP target: invalid evaluation mode',status,message)
                return
            end if
            end if
            if (present(mode)) then
                if (mode < 0 .or. mode > 2) then
                    call target_failure('OpenMP target: mode must be auto=0, direct=1 or moment=2',status,message)
                    return
                end if
            end if
            if (chebyshev) then
            if (.not. allocated(config%moment_x_power) .or. .not. allocated(config%angular_power_coefficients)) &
                then
                call target_failure('OpenMP target: uninitialized moment basis',status,message)
                return
            end if
            mm = max(mm,config%number_of_angular_moments); ma = max(ma,config%angular_order+1)
            end if
            mn = max(mn,maxval(net%nodes)); ml = max(ml,net%nlayers); mw = max(mw,count_weights)
        end do
        ! Explicit host execution permits correctness/performance comparison of
        ! exactly the same kernels. Accidental GPU fallback remains an error.
        host_requested = .false.
        if (present(use_host)) host_requested = use_host
        self%device = omp_get_default_device()
        if (present(device)) self%device = device
        if (host_requested) then
            self%device = omp_get_initial_device()
        else
        if (self%device < 0 .or. self%device >= omp_get_num_devices()) then
            call target_failure('OpenMP target: GPU device unavailable',status,message)
            return
        end if
        end if
        ! Probe once per model initialization, not once per evaluation.
        on_host = .true.
        !$omp target device(self%device) if(.not.host_requested) map(from:on_host)
        on_host = omp_is_initial_device()
        !$omp end target
        if (on_host .neqv. host_requested) then
            call target_failure('OpenMP target: actual execution device differs from requested backend',status,message)
            return
        end if
        mm = max(1,mm); ma = max(1,ma)
        allocate(self%features(6,mn,ns),self%feature_params(7,mn,ns))
        self%features = 0; self%feature_params = 0
        allocate(self%mp(4,mm,ns),self%multiplicity(mm,ns),self%polynomial(ma,ma,ns))
        self%mp = 0; self%multiplicity = 0; self%polynomial = 0
        allocate(self%meta(10,ns), self%nodes(ml,ns), self%acts(ml,ns), self%woffset(ml,ns), &
            self%local_species(ns,ns), self%weights(mw,ns), self%params(5,ns), &
            self%shift(mn,ns), self%scale(mn,ns), self%spin(ns,ns))
        self%meta = 0; self%params = 0
        self%nodes = 0; self%acts = 0; self%woffset = 0; self%weights = 0
        self%shift = 0; self%scale = 0; self%spin = 0
        do s = 1, ns
            net = model%networks(s); d = net%nodes(1); l = net%nlayers
            if (allocated(model%setups(s)%model%chebyshev)) then
            config = model%setups(s)%model%chebyshev(1)%config
            self%meta(:9,s) = [config%radial_order+1, config%angular_order+1, merge(1,0,config%num_species > 1), &
                config%version, config%cutoff_type, l, config%central_type_index, &
                config%evaluation_mode, config%number_of_angular_moments]
            if (present(mode)) self%meta(8,s) = mode
            nm = config%number_of_angular_moments; na = config%angular_order+1
            self%mp(1,:nm,s) = config%moment_x_power
            self%mp(2,:nm,s) = config%moment_y_power
            self%mp(3,:nm,s) = config%moment_z_power
            self%mp(4,:nm,s) = config%moment_x_power+config%moment_y_power+config%moment_z_power+1
            self%multiplicity(:nm,s) = config%moment_multinomial
            self%polynomial(:na,:na,s) = config%angular_power_coefficients
            self%params(:,s) = [config%radial_rc, config%angular_rc, config%cutoff_alpha, &
                net%energy_scale, net%energy_shift+net%atomic_references(s)]
            self%local_species(:,s) = model%setups(s)%global_to_local
            do j = 1, ns
                local = self%local_species(j,s)
                if (local > 0) self%spin(j,s) = config%species_weights(local)
            end do
            else
                call pack_generic_descriptors(model%setups(s)%model,fi,fr,pack_status,pack_message)
                self%features(:,:d,s) = fi; self%feature_params(:,:d,s) = fr
                self%meta(6,s) = l; self%meta(8,s) = 1; self%meta(10,s) = 1
                self%params(:,s) = [model%maximum_cutoff,model%maximum_cutoff,0.0_real64, &
                    net%energy_scale,net%energy_shift+net%atomic_references(s)]
                self%local_species(:,s) = model%setups(s)%global_to_local
            end if
            self%nodes(:l,s) = net%nodes
            self%acts(:l-1,s) = net%activation(:l-1)
            self%woffset(:l-1,s) = net%weight_offsets(:l-1)
            self%weights(:size(net%weights),s) = net%weights
            self%shift(:d,s) = net%descriptor_shift; self%scale(:d,s) = net%descriptor_scale
        end do
        self%maxnodes = mn; self%maxlayers = ml
    end subroutine

    ! Workspace assignment starts a fresh scratch cache. Device ownership is
    ! never copied; model snapshots retain normal independent value semantics.
    impure elemental subroutine copy_workspace(self, other)
        class(target_workspace), intent(out) :: self
        type(target_workspace), intent(in) :: other
        self%growth_count = 0
        self%upload_count = 0
    end subroutine

    ! NVHPC 25.3 dispatches an elemental finalizer incorrectly when a workspace
    ! is deallocated through the C handle. Explicit rank finalizers avoid it.
    subroutine finalize_workspace(self)
        type(target_workspace), intent(inout) :: self
        call self%release()
    end subroutine

    subroutine finalize_workspace_array(self)
        type(target_workspace), intent(inout) :: self(:)
        integer :: i
        do i = 1,size(self)
            call self(i)%release()
        end do
    end subroutine

    subroutine release_buffers(self)
        class(target_workspace), intent(inout) :: self
        if (.not. allocated(self%g)) return
        call map_buffers(self%cached%device,.false.,self%g,self%values,self%deriv,self%delta,self%moments,self%powers, &
            self%species,self%centers,self%offsets,self%indices,self%use_moment,self%dr,self%energies,self%forces,self%virial, &
            self%geom,self%edge_row,self%edge_force)
        deallocate(self%g,self%values,self%deriv,self%delta,self%moments,self%powers,self%species,self%centers, &
            self%offsets,self%indices,self%use_moment,self%dr,self%energies,self%forces,self%virial, &
            self%geom,self%edge_row,self%edge_force)
    end subroutine

    subroutine release_target_workspace(self)
        class(target_workspace), intent(inout) :: self
        call release_buffers(self)
        if (allocated(self%cached%meta)) call map_model(self%cached,.false.)
        call self%cached%release()
        self%growth_count = 0; self%upload_count = 0
    end subroutine

    integer function target_uploads(self) result(n)
        class(target_workspace), intent(in) :: self
        n = self%upload_count
    end function

    logical function same_model(a,b) result(same)
        type(target_model), intent(in) :: a,b
        same = .false.
        if (.not. allocated(b%meta)) return
        if (a%device /= b%device .or. a%maxnodes /= b%maxnodes .or. a%maxlayers /= b%maxlayers) return
        if (size(a%meta,2) /= size(b%meta,2)) return
        if (any(shape(a%mp) /= shape(b%mp))) return
        if (any(shape(a%polynomial) /= shape(b%polynomial))) return
        if (any(a%meta /= b%meta)) return
        if (any(a%features /= b%features) .or. any(a%feature_params /= b%feature_params)) return
        if (any(a%local_species /= b%local_species)) return
        if (any(a%nodes /= b%nodes)) return
        if (any(a%acts /= b%acts) .or. any(a%woffset /= b%woffset)) return
        if (any(a%weights /= b%weights) .or. any(a%params /= b%params)) return
        if (any(a%shift /= b%shift) .or. any(a%scale /= b%scale) .or. any(a%spin /= b%spin)) return
        if (any(a%mp /= b%mp) .or. any(a%multiplicity /= b%multiplicity) .or. &
            any(a%polynomial /= b%polynomial)) return
        same = .true.
    end function

    integer function target_allocations(self) result(n)
        class(target_workspace), intent(in) :: self
        n = self%growth_count
    end function

    subroutine evaluate_batch_target(model, species, centers, offsets, indices, displacements, energies, forces, work, &
                                     virial, profile, status, message)
        type(target_model), intent(in) :: model
        integer, intent(in) :: species(:), centers(:), offsets(:), indices(:)
        real(real64), intent(in) :: displacements(:,:)
        real(real64), intent(out) :: energies(:)
        real(real64), intent(inout) :: forces(:,:)
        type(target_workspace), intent(inout) :: work
        real(real64), intent(inout), optional :: virial(3,3)
        type(target_profile), intent(out), optional :: profile
        integer, optional, intent(out) :: status
        character(len=*), optional, intent(out) :: message
        type(target_profile) :: timing
        real(real64) :: started, mark
        integer :: first, last, row, s, na, nrw, nedges
        if (present(status)) status = 0
        if (present(message)) message = ""
        started = omp_get_wtime()
        if (present(profile)) profile = target_profile()
        if (.not. allocated(model%meta)) then
            call target_failure('OpenMP target: initialize target_model first',status,message)
            return
        end if
        if (size(offsets) /= size(centers)+1) then
            call target_failure('OpenMP target: wrong CSR offset count',status,message)
            return
        end if
        if (size(energies) /= size(centers)) then
            call target_failure('OpenMP target: wrong energy count',status,message)
            return
        end if
        if (size(forces,1) /= 3 .or. size(forces,2) /= size(species)) then
            call target_failure('OpenMP target: wrong force shape',status,message)
            return
        end if
        if (size(displacements,1) /= 3 .or. size(displacements,2) /= size(indices)) &
            then
            call target_failure('OpenMP target: wrong displacement shape',status,message)
            return
        end if
        if (any(offsets < 1) .or. any(offsets > size(indices)+1)) then
            call target_failure('OpenMP target: CSR offset out of range',status,message)
            return
        end if
        if (any(offsets(2:) < offsets(:size(centers)))) then
            call target_failure('OpenMP target: CSR offsets must be nondecreasing',status,message)
            return
        end if
        if (any(centers < 1) .or. any(centers > size(species))) then
            call target_failure('OpenMP target: center index out of range',status,message)
            return
        end if
        first = offsets(1); last = offsets(size(offsets))-1
        if (any(indices(first:last) < 1) .or. any(indices(first:last) > size(species))) &
            then
            call target_failure('OpenMP target: neighbor index out of range',status,message)
            return
        end if
        if (any(species < 1) .or. any(species > size(model%meta,2))) then
            call target_failure('OpenMP target: species out of range',status,message)
            return
        end if
        do row = 1, size(centers)
            s = species(centers(row)); first = offsets(row); last = offsets(row+1)-1
            if (any(model%local_species(species(indices(first:last)),s) == 0)) &
                then
                call target_failure('OpenMP target: neighbor species absent from model environment',status,message)
                return
            end if
            if (model%meta(4,s) == 10 .and. model%meta(3,s) == 1 .and. last-first+1 < model%meta(7,s)) &
                then
                call target_failure('OpenMP target: too few neighbors for version 10 center lookup',status,message)
                return
            end if
        end do
        if (size(centers) == 0) return
        nrw = size(centers); na = size(species)
        first = offsets(1); last = offsets(nrw+1)-1; nedges = last-first+1
        call prepare_workspace(model,work,nrw,na,nedges)
        timing%prepare = omp_get_wtime()-started
        mark = omp_get_wtime()
        work%species(:size(species)) = species; work%centers(:nrw) = centers
        work%offsets(:nrw+1) = offsets-first+1
        work%indices(:nedges) = indices(first:last); work%dr(:,:nedges) = displacements(:,first:last)
        call upload_inputs(model%device,size(species),nrw,nedges,work%species,work%centers,work%offsets,work%indices,work%dr)
        timing%upload = omp_get_wtime()-mark
        call execute_workspace(model,work,nrw,size(species),nedges,energies,forces,virial,timing)
        timing%total = omp_get_wtime()-started
        if (present(profile)) profile = timing
    end subroutine

    subroutine target_failure(text,status,message)
        character(len=*), intent(in) :: text
        integer, optional, intent(out) :: status
        character(len=*), optional, intent(out) :: message
        if (.not. present(status)) then
            print *, text
            error stop 1
        end if
        status = 1
        if (present(message)) message = text
    end subroutine

    subroutine prepare_workspace(model,work,nrw,natoms,nedges)
        type(target_model), intent(in) :: model
        type(target_workspace), intent(inout) :: work
        integer, intent(in) :: nrw,natoms,nedges
        integer :: n,mn,ml,ne,na,mm,mpower
        logical :: grow
        na = natoms
        ! Each workspace owns its model cache; compare values, not addresses, so
        ! reloads and independently copied target_model objects cannot go stale.
        if (allocated(work%cached%meta)) then
            if (work%cached%device /= model%device) call work%release()
        end if
        if (.not. same_model(model,work%cached)) then
            if (allocated(work%cached%meta)) call map_model(work%cached,.false.)
            work%cached = model
            call map_model(work%cached,.true.)
            work%upload_count = work%upload_count+1
        end if
        n = nrw; mn = model%maxnodes; ml = model%maxlayers
        ne = max(1,nedges); mm = size(model%mp,2); mpower = size(model%polynomial,1)
        grow = .true.
        if (allocated(work%g)) then
            n = max(n,size(work%g,1)); mn = max(mn,size(work%g,2)); ml = max(ml,size(work%values,3))
            ne = max(ne,size(work%indices)); na = max(na,size(work%species)); mm = max(mm,size(work%moments,2))
            mpower = max(mpower,size(work%powers,2))
            grow = n > size(work%g,1) .or. mn > size(work%g,2) .or. ml > size(work%values,3) .or. &
                ne > size(work%indices) .or. na > size(work%species) .or. mm > size(work%moments,2) .or. &
                mpower > size(work%powers,2)
        end if
        if (grow) then
            call release_buffers(work)
            allocate(work%g(n,mn),work%values(n,mn,ml),work%deriv(n,mn,ml),work%delta(n,mn,2), &
                work%moments(n,mm,2),work%powers(ne,mpower,3),work%species(na),work%centers(n),work%offsets(n+1), &
                work%indices(ne),work%use_moment(n),work%dr(3,ne),work%geom(7,ne),work%edge_row(ne),work%edge_force(3,ne), &
                work%energies(n),work%forces(3,na),work%virial(3,3))
            call map_buffers(model%device,.true.,work%g,work%values,work%deriv,work%delta,work%moments,work%powers, &
                work%species,work%centers,work%offsets,work%indices,work%use_moment,work%dr, &
                work%energies,work%forces,work%virial,work%geom,work%edge_row,work%edge_force)
            work%growth_count = work%growth_count+1
        end if
    end subroutine

    subroutine execute_workspace(model,work,nrw,natoms,nedges,energies,forces,virial,timing)
        type(target_model), intent(in) :: model
        type(target_workspace), intent(inout) :: work
        integer, intent(in) :: nrw,natoms,nedges
        real(real64), intent(out) :: energies(:)
        real(real64), intent(inout) :: forces(:,:)
        real(real64), optional, intent(inout) :: virial(3,3)
        type(target_profile), intent(inout) :: timing
        real(real64) :: stages(3),mark
        call run_target_batch(model%device,work%cached%meta,work%cached%nodes,work%cached%acts,work%cached%woffset, &
            work%cached%weights,work%cached%params,work%cached%shift,work%cached%scale,work%cached%spin, &
            work%cached%features,work%cached%feature_params,work%cached%local_species, &
            work%cached%mp,work%cached%multiplicity,work%cached%polynomial,work%species,work%centers,work%offsets, &
            work%indices,work%dr,work%energies,work%forces,work%virial,work%g,work%values,work%deriv,work%delta, &
            work%moments,work%powers,work%use_moment,work%geom,work%edge_row,work%edge_force,nrw,natoms,nedges,stages)
        timing%descriptors = stages(1); timing%network = stages(2); timing%forces = stages(3)
        mark = omp_get_wtime()
        call download_outputs(model%device,natoms,nrw,work%energies,work%forces,work%virial)
        energies = work%energies(:nrw)
        forces = forces+work%forces(:,:natoms)
        if (present(virial)) virial = virial+work%virial
        timing%download = omp_get_wtime()-mark
    end subroutine

    ! Borrow CUDA device arrays owned by lib/gpu. No host neighbor round trip.
    ! Layout is the GPU package's unpacked matrix with threads_per_atom=1.
    subroutine evaluate_lammps_target(model,work,natoms,nrows,maxnb,pitch,xptr,nptr,energies,forces,virial)
        type(target_model), intent(in) :: model
        type(target_workspace), intent(inout) :: work
        integer, intent(in) :: natoms,nrows,maxnb,pitch
        type(c_ptr), intent(in) :: xptr,nptr
        real(real64), intent(out) :: energies(nrows)
        real(real64), intent(inout) :: forces(3,natoms),virial(3,3)
        real(real64), pointer :: x(:,:)
        integer, pointer :: nb(:,:)
        integer :: ne
        type(target_profile) :: timing
        if (nrows == 0) return
        call prepare_workspace(model,work,nrows,natoms,nrows*maxnb)
        call c_f_pointer(xptr,x,[4,natoms])
        call c_f_pointer(nptr,nb,[pitch,maxnb+2])
        call import_lammps_inputs(model%device,natoms,nrows,pitch,maxnb,x,nb,work%species,work%centers, &
            work%offsets,work%indices,work%dr,ne)
        call execute_workspace(model,work,nrows,natoms,ne,energies,forces,virial,timing)
    end subroutine

    subroutine import_lammps_inputs(device,natoms,nrows,pitch,maxnb,x,nb,species,centers,offsets,indices,dr,ne)
        integer, intent(in) :: device,natoms,nrows,pitch,maxnb
        real(real64), intent(in) :: x(4,natoms)
        integer, intent(in) :: nb(pitch,maxnb+2)
        integer, contiguous, intent(inout) :: species(:),centers(:),offsets(:),indices(:)
        real(real64), contiguous, intent(inout) :: dr(:,:)
        integer, intent(out) :: ne
        integer :: row,j,k,i,e,total
        !$omp target device(device) if(device /= omp_get_initial_device()) is_device_ptr(nb) map(alloc:offsets) map(from:total)
        offsets(1) = 1
        do row = 1,nrows
            offsets(row+1) = offsets(row)+nb(row,2)
        end do
        total = offsets(nrows+1)-1
        !$omp end target
        ne = total
        !$omp target teams distribute parallel do device(device) if(device /= omp_get_initial_device()) is_device_ptr(x,nb) &
        !$omp& map(alloc:species,centers,offsets,indices,dr) private(i,j,k,e)
        do row = 1,max(natoms,nrows)
            if (row <= natoms) species(row) = int(x(4,row))
            if (row > nrows) cycle
            i = nb(row,1)+1
            centers(row) = i
            do k = 1,nb(row,2)
                j = iand(nb(row,k+2),int(z'3FFFFFFF'))+1
                e = offsets(row)+k-1
                indices(e) = j
                dr(:,e) = x(1:3,j)-x(1:3,i)
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine

    subroutine map_model(model,enter)
        type(target_model), intent(inout) :: model
        logical, intent(in) :: enter
        call map_model_arrays(model%device,enter,model%meta,model%nodes,model%acts,model%woffset,model%weights, &
            model%params,model%shift,model%scale,model%spin,model%mp,model%multiplicity,model%polynomial, &
            model%features,model%feature_params,model%local_species)
    end subroutine

    subroutine map_model_arrays(device,enter,meta,nodes,acts,woffset,weights,params,shift,scale,spin,mp,multiplicity,polynomial, &
                                features,feature_params,local_species)
        integer, intent(in) :: device
        logical, intent(in) :: enter
        integer, contiguous, intent(inout) :: meta(:,:),nodes(:,:),acts(:,:),woffset(:,:),mp(:,:,:)
        real(real64), contiguous, intent(inout) :: weights(:,:),params(:,:),shift(:,:),scale(:,:),spin(:,:)
        real(real64), contiguous, intent(inout) :: multiplicity(:,:),polynomial(:,:,:)
        integer, contiguous, intent(inout) :: features(:,:,:),local_species(:,:)
        real(real64), contiguous, intent(inout) :: feature_params(:,:,:)
        if (enter) then
            !$omp target enter data device(device) if(device /= omp_get_initial_device()) &
            !$omp& map(to:features,feature_params,local_species) &
            !$omp& map(to:meta,nodes,acts,woffset,weights,params,shift,scale,spin,mp,multiplicity,polynomial)
        else
            !$omp target exit data device(device) if(device /= omp_get_initial_device()) &
            !$omp& map(delete:features,feature_params,local_species) &
            !$omp& map(delete:meta,nodes,acts,woffset,weights,params,shift,scale,spin,mp,multiplicity,polynomial)
        end if
    end subroutine

    subroutine map_buffers(device,enter,g,values,deriv,delta,moments,powers,species,centers,offsets,indices, &
                           use_moment,dr,energies,forces,virial,geom,edge_row,edge_force)
        integer, intent(in) :: device
        logical, intent(in) :: enter
        real(real64), contiguous, intent(inout) :: g(:,:),values(:,:,:),deriv(:,:,:),delta(:,:,:),moments(:,:,:),powers(:,:,:)
        integer, contiguous, intent(inout) :: species(:),centers(:),offsets(:),indices(:),use_moment(:),edge_row(:)
        real(real64), contiguous, intent(inout) :: dr(:,:),energies(:),forces(:,:),virial(:,:),geom(:,:),edge_force(:,:)
        if (enter) then
            !$omp target enter data device(device) if(device /= omp_get_initial_device()) &
            !$omp& map(alloc:g,values,deriv,delta,moments,powers,species,centers,offsets,indices, &
            !$omp& use_moment,dr,energies,forces,virial,geom,edge_row,edge_force)
        else
            !$omp target exit data device(device) if(device /= omp_get_initial_device()) &
            !$omp& map(delete:g,values,deriv,delta,moments,powers,species,centers,offsets,indices, &
            !$omp& use_moment,dr,energies,forces,virial,geom,edge_row,edge_force)
        end if
    end subroutine

    subroutine upload_inputs(device,natoms,nrows,nedges,species,centers,offsets,indices,dr)
        integer, intent(in) :: device,natoms,nrows,nedges
        integer, contiguous, intent(in) :: species(:),centers(:),offsets(:),indices(:)
        real(real64), contiguous, intent(in) :: dr(:,:)
        !$omp target update device(device) &
        !$omp& if(device /= omp_get_initial_device()) to(species(:natoms),centers(:nrows),offsets(:nrows+1))
        if (nedges > 0) then
            !$omp target update device(device) if(device /= omp_get_initial_device()) to(indices(:nedges),dr(:,:nedges))
        end if
    end subroutine

    subroutine download_outputs(device,natoms,nrows,energies,forces,virial)
        integer, intent(in) :: device,natoms,nrows
        real(real64), contiguous, intent(inout) :: energies(:),forces(:,:),virial(:,:)
        !$omp target update device(device) if(device /= omp_get_initial_device()) from(energies(:nrows),forces(:,:natoms),virial)
    end subroutine
end module
