! Synchronous wall-clock measurements including transfers and optional neighbor
! construction. Persistent device allocation/model upload are warmed up first.
program target_benchmark
    use iso_fortran_env, only: real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet_batch, only: batch_workspace, evaluate_batch => evaluate_batch_reference
    use accelnet_batch, only: evaluate_batch_shared => evaluate_batch
    use accelnet_batch_target, only: target_model, target_workspace, target_profile, evaluate_batch_target
    use accelnet_target_runtime, only: omp_get_wtime
    use accelnet_predictor, only: predictor_model
    use accelnet_descriptors, only: atomic_structure, neighbor_data, build_neighbor_list
    use batch_test_support
    implicit none
    type(predictor_model) :: model
    type(target_model) :: packed
    type(target_workspace) :: gpuwork
    type(target_profile) :: profile
    type(batch_workspace) :: cpuwork, sharedwork
    type(atomic_structure) :: s
    type(neighbor_data) :: nb
    real(real64), allocatable :: ec(:), eg(:), fc(:,:), fg(:,:)
    integer, allocatable :: centers(:)
    real(real64) :: wc(3,3), wg(3,3), seconds, elapsed, times(4), error, phase(8), phase_sum(8), mark, neighbor_seconds, spacing
    integer(int64) :: start, finish, rate, repeats
    integer :: n, i, sample, order, method, j, mode, allocations, uploads, methods, fixed_neighbors
    character(len=128) :: arg, family, backend
    logical :: host
    if (command_argument_count() < 3 .or. command_argument_count() > 9) &
        error stop 'usage: accelnet-target-benchmark N ORDER SECONDS [MODE] [SPACING] '// &
            '[gpu|host|cpu-shared] [FAMILY] [no-neighbors] [ENV_NEIGHBORS]'
    call get_command_argument(1,arg); read(arg,*) n
    call get_command_argument(2,arg); read(arg,*) order
    call get_command_argument(3,arg); read(arg,*) seconds
    methods = 4
    if (command_argument_count() >= 8) then
        call get_command_argument(8,arg)
        if (arg /= 'no-neighbors') error stop 'expected no-neighbors'
        methods = 2
    end if
    fixed_neighbors = 0
    if (command_argument_count() == 9) then
        call get_command_argument(9,arg); read(arg,*) fixed_neighbors
        if (fixed_neighbors < 1) error stop 'positive environment neighbor count required'
    end if
    mode = 0; spacing = 1.7_real64
    if (command_argument_count() >= 4) then
        call get_command_argument(4,arg); read(arg,*) mode
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5,arg); read(arg,*) spacing
    end if
    if (n < 1 .or. order < 0 .or. seconds <= 0) error stop 'invalid benchmark parameters'
    family = 'chebyshev'; backend = 'gpu'
    if (command_argument_count() >= 6) call get_command_argument(6,backend)
    if (command_argument_count() >= 7) call get_command_argument(7,family)
    if (backend /= 'host' .and. backend /= 'gpu' .and. backend /= 'cpu-shared') error stop 'invalid backend'
    host = backend /= 'gpu'
    call make_model(trim(family),model,order=order)
    call model%set_chebyshev_evaluation(mode)
    ! Generic mode 0 also measures the established CPU's automatic G5 policy;
    ! mode 1 compares direct algorithms. Chebyshev uses its own mode above.
    if (family == 'chebyshev') then
        call model%set_g5_evaluation(1)
    else
        call model%set_g5_evaluation(mode)
    end if
    call make_structure(n,2,s)
    if (spacing <= 0) error stop 'spacing must be positive'
    s%positions = s%positions*(spacing/1.7_real64); s%lattice = s%lattice*(spacing/1.7_real64)
    if (family == 'chebyshev') then
        call packed%initialize(model,mode=mode,use_host=host)
    else
        call packed%initialize(model,g5_mode=mode,use_host=host)
    end if
    allocate(ec(n),eg(n),fc(3,n),fg(3,n),centers(n))
    centers = [(i,i=1,n)]
    if (fixed_neighbors > 0) then
        call make_g5_scaling_neighbors(n,fixed_neighbors,nb)
    else
        call build_neighbor_list(s,model%maximum_cutoff,nb,model%minimum_distance)
    end if
    call evaluate(1); call evaluate(2)
    error = max(maxval(abs(ec-eg)),maxval(abs(fc-fg)),maxval(abs(wc-wg)))
    call check()
    allocations = gpuwork%allocations(); uploads = gpuwork%uploads()
    write(*,*) 'BACKEND ',trim(backend),' FAMILY ',trim(family),' MODE ',mode
    write(*,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'CASE',n,order,size(nb%atom_indices),error
    times = 0
    do sample = 1, 5
        do j = 1, methods
            method = j
            if (mod(sample,2) == 0) method = methods+1-j
            call evaluate(method) ! warmup also excludes one-time allocations
            call system_clock(start,rate)
            if (rate <= 0) error stop 'wall clock unavailable'
            repeats = 0
            phase_sum = 0
            do
                call evaluate(method)
                phase_sum = phase_sum+phase
                repeats = repeats+1
                call system_clock(finish)
                elapsed = real(finish-start,real64)/real(rate,real64)
                if (elapsed >= seconds) exit
            end do
            times(method) = elapsed/real(repeats,real64)
            call check()
            if (mod(method,2) == 0) write(*,'(A,2(1X,I0),8(1X,ES24.16))') &
                'PROFILE',sample,method,phase_sum/real(repeats,real64)
        end do
        write(*,'(A,1X,I0,4(1X,ES24.16))') 'TIMING',sample,times
    end do
    if (gpuwork%allocations() /= allocations .or. gpuwork%uploads() /= uploads) &
        error stop 'resident buffers or model were reallocated during steady-state timing'
    write(*,'(A,2(1X,I0))') 'RESIDENCY',gpuwork%allocations(),gpuwork%uploads()
    call gpuwork%release()
contains
    subroutine evaluate(method)
        integer, intent(in) :: method
        phase = 0; neighbor_seconds = 0
        if (method > 2) then
            mark = omp_get_wtime()
            call build_neighbor_list(s,model%maximum_cutoff,nb,model%minimum_distance)
            neighbor_seconds = omp_get_wtime()-mark
        end if
        if (mod(method,2) == 1) then
            fc = 0; wc = 0
            call evaluate_batch(model,s%species,centers,nb%offsets,nb%atom_indices,nb%displacements,ec,fc,cpuwork,wc)
        else
            fg = 0; wg = 0
            if (backend == 'cpu-shared') then
                call evaluate_batch_shared(model,s%species,centers,nb%offsets,nb%atom_indices, &
                    nb%displacements,eg,fg,sharedwork,wg)
                profile = target_profile()
            else
                call evaluate_batch_target(packed,s%species,centers,nb%offsets,nb%atom_indices, &
                    nb%displacements,eg,fg,gpuwork,wg,profile)
            end if
            phase = [neighbor_seconds,profile%prepare,profile%upload,profile%descriptors,profile%network, &
                profile%forces,profile%download,profile%total]
        end if
    end subroutine

    subroutine check()
        if (.not. all(ieee_is_finite(eg)) .or. .not. all(ieee_is_finite(fg)) .or. &
            .not. all(ieee_is_finite(wg))) error stop 'nonfinite GPU result'
        if (any(abs(ec-eg) > 2e-10_real64+2e-10_real64*abs(ec))) error stop 'GPU energy mismatch'
        if (any(abs(fc-fg) > 2e-10_real64+2e-10_real64*abs(fc))) error stop 'GPU force mismatch'
        if (any(abs(wc-wg) > 2e-10_real64+2e-10_real64*abs(wc))) error stop 'GPU virial mismatch'
    end subroutine
end program
