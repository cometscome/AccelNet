! Optional OpenMP target backend. The established CPU library does not depend
! on this module or on an OpenMP runtime.
module accelnet_batch_target
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use omp_lib, only: omp_get_num_devices, omp_get_default_device
    use accelnet_predictor, only: predictor_model
    use accelnet_descriptors, only: descriptor_config, validate_cutoff_parameters
    use aenet_network, only: atomic_network
    use accelnet_target_kernels, only: run_target_batch
    implicit none
    private
    public :: target_model, target_workspace, evaluate_batch_target

    type :: target_model
        private
        integer :: device = -1, maxnodes = 0, maxlayers = 0
        integer, allocatable :: meta(:,:), nodes(:,:), acts(:,:), woffset(:,:), local_species(:,:)
        real(real64), allocatable :: weights(:,:), params(:,:), shift(:,:), scale(:,:), spin(:,:)
    contains
        procedure :: initialize => initialize_target_model
        procedure :: release => release_target_model
    end type

    type :: target_workspace
        private
        real(real64), allocatable :: g(:,:), values(:,:,:), deriv(:,:,:), delta(:,:,:)
        integer :: growth_count = 0
    contains
        procedure :: release => release_target_workspace
        procedure :: allocations => target_allocations
    end type
contains
    subroutine release_target_model(self)
        class(target_model), intent(inout) :: self
        if (allocated(self%meta)) deallocate(self%meta, self%nodes, self%acts, self%woffset, self%local_species, &
            self%weights, self%params, self%shift, self%scale, self%spin)
        self%device = -1; self%maxnodes = 0; self%maxlayers = 0
    end subroutine

    subroutine initialize_target_model(self, model, device)
        class(target_model), intent(inout) :: self
        type(predictor_model), intent(in) :: model
        integer, intent(in), optional :: device
        type(descriptor_config) :: config
        type(atomic_network) :: net
        integer :: ns, s, d, ml, mn, mw, l, count_weights, j, local
        call self%release()
        if (.not. allocated(model%networks) .or. .not. allocated(model%setups)) &
            error stop 'OpenMP target: uninitialized model'
        ns = size(model%networks)
        if (ns == 0 .or. size(model%setups) /= ns) error stop 'OpenMP target: invalid species count'
        mn = 0; ml = 0; mw = 0
        do s = 1, ns
            if (allocated(model%setups(s)%model%lj) .or. allocated(model%setups(s)%model%behler)) &
                error stop 'OpenMP target: only a single Chebyshev component per species is supported'
            if (.not. allocated(model%setups(s)%model%chebyshev)) &
                error stop 'OpenMP target: model has no Chebyshev component'
            if (size(model%setups(s)%model%chebyshev) /= 1) &
                error stop 'OpenMP target: composite Chebyshev models are not supported yet'
            config = model%setups(s)%model%chebyshev(1)%config
            net = model%networks(s)
            call validate_cutoff_parameters(config%cutoff_type, config%cutoff_alpha)
            if (config%radial_rc <= 0 .or. config%angular_rc <= 0 .or. &
                config%radial_order < 0 .or. config%angular_order < 0) error stop 'OpenMP target: invalid Chebyshev config'
            if (config%version /= 0 .and. config%version /= 1 .and. config%version /= 10) &
                error stop 'OpenMP target: unsupported Chebyshev version'
            if (config%version == 10 .and. config%num_species > 1 .and. config%central_type_index < 1) &
                error stop 'OpenMP target: invalid version 10 center lookup'
            if (.not. allocated(config%species_weights)) error stop 'OpenMP target: missing species weights'
            if (size(config%species_weights) /= config%num_species) error stop 'OpenMP target: invalid species weights'
            if (.not. allocated(model%setups(s)%global_to_local)) error stop 'OpenMP target: missing species map'
            if (size(model%setups(s)%global_to_local) /= ns) error stop 'OpenMP target: invalid species map size'
            if (any(model%setups(s)%global_to_local < 0) .or. &
                any(model%setups(s)%global_to_local > config%num_species)) error stop 'OpenMP target: invalid species map values'
            if (.not. allocated(net%nodes) .or. .not. allocated(net%activation) .or. &
                .not. allocated(net%weight_offsets) .or. .not. allocated(net%weights)) &
                error stop 'OpenMP target: uninitialized network'
            if (net%nlayers < 2 .or. size(net%nodes) /= net%nlayers) error stop 'OpenMP target: invalid network depth'
            if (any(net%nodes < 1) .or. net%nodes(net%nlayers) /= 1) error stop 'OpenMP target: invalid network width'
            d = net%nodes(1)
            if (d /= config%num_descriptors()) error stop 'OpenMP target: descriptor/network size mismatch'
            if (size(net%activation) < net%nlayers-1 .or. size(net%weight_offsets) < net%nlayers-1) &
                error stop 'OpenMP target: invalid network metadata'
            if (any(net%activation(:net%nlayers-1) < 0) .or. any(net%activation(:net%nlayers-1) > 11)) &
                error stop 'OpenMP target: unsupported network activation'
            count_weights = 0
            do l = 1, net%nlayers-1
                if (net%weight_offsets(l) /= count_weights) error stop 'OpenMP target: invalid weight offsets'
                if (net%nodes(l) >= huge(j)/net%nodes(l+1)-1) error stop 'OpenMP target: network too large'
                j = (net%nodes(l)+1)*net%nodes(l+1)
                if (j > huge(j)-count_weights) error stop 'OpenMP target: network too large'
                count_weights = count_weights+j
            end do
            if (size(net%weights) /= count_weights) error stop 'OpenMP target: invalid weights'
            if (.not. allocated(net%descriptor_shift) .or. .not. allocated(net%descriptor_scale) .or. &
                .not. allocated(net%atomic_references)) error stop 'OpenMP target: missing normalization'
            if (size(net%descriptor_shift) /= d .or. size(net%descriptor_scale) /= d .or. &
                size(net%atomic_references) /= ns) error stop 'OpenMP target: invalid normalization sizes'
            if (.not. ieee_is_finite(net%energy_scale) .or. net%energy_scale == 0) &
                error stop 'OpenMP target: invalid energy scale'
            mn = max(mn,maxval(net%nodes)); ml = max(ml,net%nlayers); mw = max(mw,count_weights)
        end do
        self%device = omp_get_default_device()
        if (present(device)) self%device = device
        if (self%device < 0 .or. self%device >= omp_get_num_devices()) error stop 'OpenMP target: GPU device unavailable'
        allocate(self%meta(7,ns), self%nodes(ml,ns), self%acts(ml,ns), self%woffset(ml,ns), &
            self%local_species(ns,ns), self%weights(mw,ns), self%params(5,ns), &
            self%shift(mn,ns), self%scale(mn,ns), self%spin(ns,ns))
        self%nodes = 0; self%acts = 0; self%woffset = 0; self%weights = 0
        self%shift = 0; self%scale = 0; self%spin = 0
        do s = 1, ns
            config = model%setups(s)%model%chebyshev(1)%config
            net = model%networks(s); d = net%nodes(1); l = net%nlayers
            self%meta(:,s) = [config%radial_order+1, config%angular_order+1, merge(1,0,config%num_species > 1), &
                config%version, config%cutoff_type, l, config%central_type_index]
            self%params(:,s) = [config%radial_rc, config%angular_rc, config%cutoff_alpha, &
                net%energy_scale, net%energy_shift+net%atomic_references(s)]
            self%local_species(:,s) = model%setups(s)%global_to_local
            do j = 1, ns
                local = self%local_species(j,s)
                if (local > 0) self%spin(j,s) = config%species_weights(local)
            end do
            self%nodes(:l,s) = net%nodes
            self%acts(:l-1,s) = net%activation(:l-1)
            self%woffset(:l-1,s) = net%weight_offsets(:l-1)
            self%weights(:size(net%weights),s) = net%weights
            self%shift(:d,s) = net%descriptor_shift; self%scale(:d,s) = net%descriptor_scale
        end do
        self%maxnodes = mn; self%maxlayers = ml
    end subroutine

    subroutine release_target_workspace(self)
        class(target_workspace), intent(inout) :: self
        if (allocated(self%g)) deallocate(self%g, self%values, self%deriv, self%delta)
        self%growth_count = 0
    end subroutine

    integer function target_allocations(self) result(n)
        class(target_workspace), intent(in) :: self
        n = self%growth_count
    end function

    subroutine evaluate_batch_target(model, species, centers, offsets, indices, displacements, energies, forces, work, virial)
        type(target_model), intent(in) :: model
        integer, intent(in) :: species(:), centers(:), offsets(:), indices(:)
        real(real64), intent(in) :: displacements(:,:)
        real(real64), intent(out) :: energies(:)
        real(real64), intent(inout) :: forces(:,:)
        type(target_workspace), intent(inout) :: work
        real(real64), intent(inout), optional :: virial(3,3)
        real(real64) :: w(3,3)
        integer :: first, last, row, s, n, mn, ml, old_growth
        logical :: grow
        if (.not. allocated(model%meta)) error stop 'OpenMP target: initialize target_model first'
        if (size(offsets) /= size(centers)+1) error stop 'OpenMP target: wrong CSR offset count'
        if (size(energies) /= size(centers)) error stop 'OpenMP target: wrong energy count'
        if (size(forces,1) /= 3 .or. size(forces,2) /= size(species)) error stop 'OpenMP target: wrong force shape'
        if (size(displacements,1) /= 3 .or. size(displacements,2) /= size(indices)) &
            error stop 'OpenMP target: wrong displacement shape'
        if (any(offsets < 1) .or. any(offsets > size(indices)+1)) error stop 'OpenMP target: CSR offset out of range'
        if (any(offsets(2:) < offsets(:size(centers)))) error stop 'OpenMP target: CSR offsets must be nondecreasing'
        if (any(centers < 1) .or. any(centers > size(species))) error stop 'OpenMP target: center index out of range'
        first = offsets(1); last = offsets(size(offsets))-1
        if (any(indices(first:last) < 1) .or. any(indices(first:last) > size(species))) &
            error stop 'OpenMP target: neighbor index out of range'
        if (any(species < 1) .or. any(species > size(model%meta,2))) error stop 'OpenMP target: species out of range'
        do row = 1, size(centers)
            s = species(centers(row)); first = offsets(row); last = offsets(row+1)-1
            if (any(model%local_species(species(indices(first:last)),s) == 0)) &
                error stop 'OpenMP target: neighbor species absent from model environment'
            if (model%meta(4,s) == 10 .and. model%meta(3,s) == 1 .and. last-first+1 < model%meta(7,s)) &
                error stop 'OpenMP target: too few neighbors for version 10 center lookup'
        end do
        if (size(centers) == 0) return
        n = size(centers); mn = model%maxnodes; ml = model%maxlayers
        grow = .true.
        if (allocated(work%g)) then
            n = max(n,size(work%g,1)); mn = max(mn,size(work%g,2)); ml = max(ml,size(work%values,3))
            grow = n > size(work%g,1) .or. mn > size(work%g,2) .or. ml > size(work%values,3)
        end if
        if (grow) then
            old_growth = work%growth_count
            call work%release()
            allocate(work%g(n,mn),work%values(n,mn,ml),work%deriv(n,mn,ml),work%delta(n,mn,2))
            work%growth_count = old_growth+1
        end if
        w = 0.0_real64
        call run_target_batch(model%device,model%meta,model%nodes,model%acts,model%woffset,model%weights, &
            model%params,model%shift,model%scale,model%spin,species,centers,offsets,indices,displacements, &
            energies,forces,w,work%g,work%values,work%deriv,work%delta)
        if (present(virial)) virial = virial+w
    end subroutine
end module
