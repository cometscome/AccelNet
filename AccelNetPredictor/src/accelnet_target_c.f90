! Instance-owned C interface for the optional GPU backend.
module accelnet_target_c
    use iso_c_binding
    use accelnet_predictor, only: predictor_model, load_predictor_from_network_data, load_predictor_from_n2p2
    use aenet_network, only: atomic_network, read_aenet_network
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet_batch_target
    implicit none
    private
    ! Keep C entry points public, including with GNU Fortran 16.2 (GCC PR126872).
    public :: accelnet_target_create, accelnet_target_create_modes
    public :: accelnet_target_create_versioned, accelnet_target_create_n2p2
    public :: accelnet_target_get_species, accelnet_target_destroy
    public :: accelnet_target_compute, accelnet_target_compute_device

    type :: context
        type(target_model), pointer :: model => null()
        type(target_workspace), pointer :: work => null()
        character(len=16), allocatable :: species(:)
    end type
contains
    ! Preserve the existing ABI and its Chebyshev mode argument.
    function accelnet_target_create(nspecies,paths,device,mode,cutoff,handle,message) result(status) bind(C)
        integer(c_int), value :: nspecies,device,mode
        type(c_ptr), intent(in) :: paths(nspecies)
        real(c_double), intent(out) :: cutoff
        type(c_ptr), intent(out) :: handle
        character(c_char), intent(out) :: message(512)
        integer(c_int) :: status
        status = accelnet_target_create_modes(nspecies,paths,device,mode,0_c_int,cutoff,handle,message)
    end function

    function accelnet_target_create_modes(nspecies,paths,device,mode,g5_mode,cutoff,handle,message) result(status) bind(C)
        integer(c_int), value :: nspecies,device,mode,g5_mode
        type(c_ptr), intent(in) :: paths(nspecies)
        real(c_double), intent(out) :: cutoff
        type(c_ptr), intent(out) :: handle
        character(c_char), intent(out) :: message(512)
        integer(c_int) :: status
        status = accelnet_target_create_versioned(nspecies,paths,device,0_c_int,mode,g5_mode,cutoff,handle,message)
    end function

    function accelnet_target_create_versioned(nspecies,paths,device,version,mode,g5_mode,cutoff,handle,message) &
        result(status) bind(C)
        integer(c_int), value :: nspecies,device,mode,g5_mode
        integer(c_int), value :: version
        type(c_ptr), intent(in) :: paths(nspecies)
        real(c_double), intent(out) :: cutoff
        type(c_ptr), intent(out) :: handle
        character(c_char), intent(out) :: message(512)
        integer(c_int) :: status
        type(predictor_model) :: source
        type(atomic_network), allocatable :: nets(:)
        integer :: k, nr, na, kind, env_count, expected_dim
        character(len=100) :: family
        real(c_double) :: rc, ac, alpha, middle
        character(c_char), pointer :: path(:)
        character(len=4096), allocatable :: filenames(:)
        integer, allocatable :: order(:)
        integer :: i,j,unit,ios
        handle = c_null_ptr; cutoff = 0; message = c_null_char; status = 1
        call c_message('Invalid or incomplete embedded model metadata',message)
        if (version /= 0 .and. version /= 1 .and. version /= 10) then
            call c_message('Chebyshev version must be 0, 1, or 10',message)
            return
        end if
        if (nspecies < 1 .or. mode < 0 .or. mode > 2 .or. g5_mode < 0 .or. g5_mode > 3) then
            call c_message('Invalid species count or evaluation mode',message)
            return
        end if
        allocate(filenames(nspecies))
        do i = 1,nspecies
            if (.not. c_associated(paths(i))) return
            call c_f_pointer(paths(i),path,[4096])
            filenames(i) = ''
            do j = 1,4096
                if (path(j) == c_null_char) exit
                filenames(i)(j:j) = path(j)
            end do
            if (j > 4096) return
            open(newunit=unit,file=trim(filenames(i)),status='old',action='read',iostat=ios)
            if (ios /= 0) then
                call c_message('Cannot open network: '//trim(filenames(i)),message)
                return
            end if
            close(unit)
        end do
        allocate(nets(nspecies))
        allocate(order(nspecies),source=0)
        do i = 1,nspecies
            call read_aenet_network(trim(filenames(i)),nets(i),ios)
            if (ios /= 0) then
                call c_message('Invalid or incomplete network: '//trim(filenames(i)),message)
                return
            end if
            family = nets(i)%descriptor_name
            do j = 1,len_trim(family)
                k = iachar(family(j:j))
                if (k >= iachar('A') .and. k <= iachar('Z')) family(j:j) = achar(k+32)
            end do
            if (size(nets(i)%species_names) /= nspecies .or. size(nets(i)%descriptor_parameters,1) < 1 .or. &
                size(nets(i)%descriptor_parameters,2) < 1) then
                call c_message('Network species count or descriptor metadata is invalid',message)
                return
            end if
            if (any(nets(i)%species_names /= nets(1)%species_names)) then
                call c_message('Networks disagree on embedded global species ordering',message)
                return
            end if
            do j=1,nspecies
                if (count(nets(i)%species_names == nets(i)%species_names(j)) /= 1) then
                    call c_message('Duplicate embedded model species',message)
                    return
                end if
            end do
            k=0
            do j=1,nspecies
                if (nets(i)%atomtype == nets(1)%species_names(j)) k=j
            end do
            if (k == 0) then
                call c_message('Network element is missing from embedded species',message)
                return
            end if
            if (order(k) /= 0) then
                call c_message('Duplicate network element',message)
                return
            end if
            order(k)=i
            do k = 1,nspecies
                if (.not. any(nets(i)%environment_names == nets(i)%species_names(k))) then
                    call c_message('Embedded environment species are incomplete',message)
                    return
                end if
            end do
            if (.not. all(ieee_is_finite(nets(i)%descriptor_parameters))) return
            alpha=nets(i)%descriptor_cutoff_alpha
            if (nets(i)%descriptor_cutoff_type < 0 .or. nets(i)%descriptor_cutoff_type > 9 .or. &
                alpha < 0 .or. alpha >= 1 .or. (nets(i)%descriptor_cutoff_type == 9 .and. alpha <= 0)) then
                call c_message('Invalid cutoff parameters',message)
                return
            end if
            env_count = size(nets(i)%environment_names)
            select case(trim(family))
            case('chebyshev')
                if (size(nets(i)%descriptor_parameters,1) < 4) return
                rc=nets(i)%descriptor_parameters(1,1); ac=nets(i)%descriptor_parameters(3,1)
                nr=nint(nets(i)%descriptor_parameters(2,1)); na=nint(nets(i)%descriptor_parameters(4,1))
                if (rc <= 0 .or. ac <= 0 .or. nr < 0 .or. na < 0 .or. nr > 10000 .or. na > 100) then
                    call c_message('Invalid or oversized Chebyshev configuration',message)
                    return
                end if
                expected_dim = (nr+na+2)*merge(2,1,env_count>1)
            case('lj')
                if (nets(i)%descriptor_parameters(1,1) <= 0) return
                expected_dim = 2*env_count
            case('n2p2_extended')
                expected_dim=size(nets(i)%descriptor_kinds)
                if (size(nets(i)%descriptor_parameters,1) /= 9 .or. &
                    size(nets(i)%descriptor_parameters,2) /= expected_dim .or. &
                    size(nets(i)%descriptor_environments,1) /= 2 .or. &
                    size(nets(i)%descriptor_environments,2) /= expected_dim) return
                do j=1,expected_dim
                    kind=nets(i)%descriptor_kinds(j)
                    if (.not. (kind == 2 .or. kind == 4 .or. kind == 5 .or. kind == 12 .or. &
                        kind == 13 .or. (kind >= 20 .and. kind <= 25))) return
                    if (nets(i)%descriptor_parameters(1,j) <= 0) return
                    k=nets(i)%descriptor_environments(1,j)
                    if (kind < 12 .or. (kind >= 20 .and. kind <= 22)) then
                        if (k < 1 .or. k > env_count) return
                        if (kind /= 2 .and. kind /= 20) then
                            k=nets(i)%descriptor_environments(2,j)
                            if (k < 1 .or. k > env_count) return
                        end if
                    else
                        if (any(nets(i)%descriptor_environments(:,j) /= 0)) return
                    end if
                    if (kind == 4 .or. kind == 5) then
                        if (abs(nets(i)%descriptor_parameters(2,j)) > 1 .or. &
                            nets(i)%descriptor_parameters(3,j) < 1 .or. &
                            nets(i)%descriptor_parameters(3,j) > real(huge(1)-1,c_double)) return
                    else if (kind == 13) then
                        if (abs(nets(i)%descriptor_parameters(4,j)) > 1 .or. &
                            nets(i)%descriptor_parameters(5,j) < 1) return
                    end if
                    if (kind >= 20) then
                        if (nets(i)%descriptor_parameters(2,j) >= nets(i)%descriptor_parameters(1,j)) return
                        k=nint(nets(i)%descriptor_parameters(5,j))
                        if (real(k,c_double) /= nets(i)%descriptor_parameters(5,j)) return
                        if (.not. (k >= 1 .and. k <= 5) .and. .not. (k >= 11 .and. k <= 14)) return
                        if (kind /= 20 .and. kind /= 23) then
                            rc=nets(i)%descriptor_parameters(3,j); ac=nets(i)%descriptor_parameters(4,j)
                            middle=0.5_c_double*(rc+ac)
                            if (rc >= ac .or. ac-rc > 360) return
                            if ((rc < 0 .and. middle /= 0) .or. (ac > 180 .and. middle /= 180)) return
                        end if
                    end if
                end do
            case('behler2011')
                expected_dim = size(nets(i)%descriptor_kinds)
                if (size(nets(i)%descriptor_parameters,1) < 4 .or. &
                    size(nets(i)%descriptor_parameters,2) /= expected_dim .or. &
                    size(nets(i)%descriptor_environments,1) < 2 .or. &
                    size(nets(i)%descriptor_environments,2) /= expected_dim) return
                do j = 1,expected_dim
                    kind = nets(i)%descriptor_kinds(j)
                    if (kind < 1 .or. kind > 5) then
                        call c_message('Unsupported embedded Behler descriptor kind',message)
                        return
                    end if
                    if (nets(i)%descriptor_parameters(1,j) <= 0) return
                    k = nets(i)%descriptor_environments(1,j)
                    if (k < 1 .or. k > env_count) return
                    if (kind >= 4) then
                        k = nets(i)%descriptor_environments(2,j)
                        if (k < 1 .or. k > env_count) return
                        if (abs(nets(i)%descriptor_parameters(2,j)) > 1 .or. &
                            nets(i)%descriptor_parameters(3,j) < 1 .or. &
                            nets(i)%descriptor_parameters(3,j) > real(huge(1)-1,c_double)) return
                    end if
                end do
            case default
                call c_message('GPU supports embedded Chebyshev, LJ, Behler2011 and n2p2_extended descriptors',message)
                return
            end select
            if (expected_dim /= nets(i)%nodes(1)) then
                call c_message('Embedded descriptor/network dimension mismatch',message)
                return
            end if
        end do
        call load_predictor_from_network_data(nets(order),source,chebyshev_version=version)
        status = create_context(source,device,mode,g5_mode,cutoff,handle,message)
    end function

    function accelnet_target_create_n2p2(directory,device,g5_mode,nspecies,cutoff,handle,message) result(status) bind(C)
        type(c_ptr), value :: directory
        integer(c_int), value :: device,g5_mode
        integer(c_int), intent(out) :: nspecies
        real(c_double), intent(out) :: cutoff
        type(c_ptr), intent(out) :: handle
        character(c_char), intent(out) :: message(512)
        integer(c_int) :: status
        type(predictor_model) :: source
        character(c_char), pointer :: path(:)
        character(len=4096) :: filename
        character(len=512) :: detail
        integer :: i,unit,ios
        status = 1; handle = c_null_ptr; cutoff = 0; nspecies = 0; message = c_null_char
        if (g5_mode < 0 .or. g5_mode > 3 .or. device < -1) then
            call c_message('Invalid device or G5 evaluation mode',message)
            return
        end if
        if (.not. c_associated(directory)) then
            call c_message('Missing n2p2 model directory',message)
            return
        end if
        call c_f_pointer(directory,path,[4096])
        filename = ''
        do i = 1,4086 ! leave room for /input.nn and a terminator
            if (path(i) == c_null_char) exit
            filename(i:i) = path(i)
        end do
        if (i == 1 .or. i > 4086) then
            call c_message('Empty or overlong n2p2 model directory',message)
            return
        end if
        open(newunit=unit,file=trim(filename)//'/input.nn',status='old',action='read',iostat=ios)
        if (ios /= 0) then
            call c_message('Cannot open n2p2 input.nn',message)
            return
        end if
        close(unit)
        call load_predictor_from_n2p2(trim(filename),source,status,detail)
        if (status /= 0) then
            call c_message(detail,message)
            return
        end if
        status = create_context(source,device,0_c_int,g5_mode,cutoff,handle,message)
        if (status == 0) nspecies = size(source%species_names)
    end function

    function create_context(source,device,mode,g5_mode,cutoff,handle,message) result(status)
        type(predictor_model), intent(in) :: source
        integer(c_int), intent(in) :: device,mode,g5_mode
        real(c_double), intent(out) :: cutoff
        type(c_ptr), intent(out) :: handle
        character(c_char), intent(out) :: message(512)
        integer(c_int) :: status
        type(context), pointer :: ctx
        character(len=511) :: detail
        status = 1; handle = c_null_ptr; cutoff = 0; message = c_null_char
        if (device < -1) then
            call c_message('Invalid device; use a GPU index or ACCELNET_TARGET_HOST',message)
            return
        end if
        allocate(ctx)
        allocate(ctx%model,ctx%work)
        call ctx%model%initialize(source,device,mode,status,detail,use_host=(device == -1),g5_mode=g5_mode)
        if (status /= 0) then
            call c_message(detail,message)
            deallocate(ctx%model,ctx%work)
            deallocate(ctx)
            return
        end if
        ctx%species = source%species_names
        cutoff = source%maximum_cutoff
        handle = c_loc(ctx)
    end function

    function accelnet_target_get_species(handle,species,capacity,symbol,message) result(status) bind(C)
        type(c_ptr), value :: handle
        integer(c_int), value :: species,capacity
        character(c_char), intent(out) :: symbol(capacity),message(512)
        integer(c_int) :: status
        type(context), pointer :: ctx
        integer :: i,n
        status = 1; message = c_null_char
        if (capacity > 0) symbol(1) = c_null_char
        if (.not. c_associated(handle)) then
            call c_message('Missing target model handle',message)
            return
        end if
        call c_f_pointer(handle,ctx)
        if (species < 1 .or. species > size(ctx%species)) then
            call c_message('Species index is outside the model (one-based)',message)
            return
        end if
        n = len_trim(ctx%species(species))
        if (capacity <= n) then
            call c_message('Species buffer is too small (include the null terminator)',message)
            return
        end if
        do i = 1,n
            symbol(i) = ctx%species(species)(i:i)
        end do
        symbol(n+1) = c_null_char
        status = 0
    end function

    subroutine accelnet_target_destroy(handle) bind(C)
        type(c_ptr), value :: handle
        type(context), pointer :: ctx
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle,ctx)
        call ctx%work%release()
        call ctx%model%release()
        deallocate(ctx%model,ctx%work)
        deallocate(ctx)
    end subroutine

    function accelnet_target_compute(handle,natoms,nrows,nedges,species,centers,offsets,indices,dr,e,f,w,message) &
        result(status) bind(C)
        type(c_ptr), value :: handle
        integer(c_int), value :: natoms,nrows,nedges
        integer(c_int), intent(in) :: species(natoms),centers(nrows),offsets(nrows+1),indices(nedges)
        real(c_double), intent(in) :: dr(3,nedges)
        real(c_double), intent(out) :: e(nrows)
        real(c_double), intent(inout) :: f(3,natoms),w(3,3)
        character(c_char), intent(out) :: message(512)
        integer(c_int) :: status
        character(len=511) :: detail
        type(context), pointer :: ctx
        status = 1; message = c_null_char
        if (.not. c_associated(handle)) then
            call c_message('Missing target model handle',message)
            return
        end if
        call c_f_pointer(handle,ctx)
        call evaluate_batch_target(ctx%model,species,centers,offsets,indices,dr,e,f,ctx%work,w, &
            status=status,message=detail)
        if (status /= 0) call c_message(detail,message)
    end function

    function accelnet_target_compute_device(handle,natoms,nrows,maxnb,pitch,x,nb,e,f,w) result(status) bind(C)
        type(c_ptr), value :: handle,x,nb
        integer(c_int), value :: natoms,nrows,maxnb,pitch
        real(c_double), intent(out) :: e(nrows)
        real(c_double), intent(inout) :: f(3,natoms),w(3,3)
        integer(c_int) :: status
        type(context), pointer :: ctx
        status = 1
        if (.not. c_associated(handle)) return
        if (natoms < nrows .or. nrows < 0 .or. maxnb < 0 .or. pitch < nrows) return
        if (nrows > 0) then
            if (maxnb > huge(1)/nrows) return
            if (.not. c_associated(x) .or. .not. c_associated(nb)) return
        end if
        call c_f_pointer(handle,ctx)
        call evaluate_lammps_target(ctx%model,ctx%work,natoms,nrows,maxnb,pitch,x,nb,e,f,w)
        status = 0
    end function

    subroutine c_message(text,message)
        character(len=*), intent(in) :: text
        character(c_char), intent(out) :: message(512)
        integer :: i
        message = c_null_char
        do i = 1,min(len_trim(text),511)
            message(i) = text(i:i)
        end do
    end subroutine
end module
