! Instance-owned C interface for the optional GPU backend.
module accelnet_target_c
    use iso_c_binding
    use accelnet_predictor, only: predictor_model, load_predictor_from_network_data
    use aenet_network, only: atomic_network, read_aenet_network
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet_batch_target
    implicit none
    private
    type :: context
        type(target_model), pointer :: model => null()
        type(target_workspace), pointer :: work => null()
    end type
contains
    function accelnet_target_create(nspecies,paths,device,mode,cutoff,handle,message) result(status) bind(C)
        integer(c_int), value :: nspecies,device,mode
        type(c_ptr), intent(in) :: paths(nspecies)
        real(c_double), intent(out) :: cutoff
        type(c_ptr), intent(out) :: handle
        character(c_char), intent(out) :: message(512)
        integer(c_int) :: status
        type(context), pointer :: ctx
        type(predictor_model) :: source
        type(atomic_network), allocatable :: nets(:)
        integer :: k, nr, na, kind, env_count, expected_dim
        character(len=100) :: family
        real(c_double) :: rc, ac, alpha
        character(c_char), pointer :: path(:)
        character(len=4096), allocatable :: filenames(:)
        character(len=511) :: detail
        integer :: i,j,unit,ios
        handle = c_null_ptr; cutoff = 0; message = c_null_char; status = 1
        if (nspecies < 1 .or. mode < 0 .or. mode > 2) then
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
            if (any(nets(i)%species_names /= nets(1)%species_names) .or. &
                nets(i)%atomtype /= nets(i)%species_names(i)) then
                call c_message('Network paths must follow embedded global species ordering',message)
                return
            end if
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
                call c_message('GPU supports embedded Chebyshev, LJ and Behler2011 descriptors',message)
                return
            end select
            if (expected_dim /= nets(i)%nodes(1)) then
                call c_message('Embedded descriptor/network dimension mismatch',message)
                return
            end if
        end do
        call load_predictor_from_network_data(nets,source)
        allocate(ctx)
        allocate(ctx%model,ctx%work)
        call ctx%model%initialize(source,device,mode,status,detail)
        if (status /= 0) then
            call c_message(detail,message)
            deallocate(ctx%model,ctx%work)
            deallocate(ctx)
            return
        end if
        cutoff = source%maximum_cutoff
        handle = c_loc(ctx)
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
        if (.not. c_associated(handle)) return
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
