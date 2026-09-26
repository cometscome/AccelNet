! Synchronous wall-clock measurements including target allocation and transfers.
program target_benchmark
    use iso_fortran_env, only: real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet_batch, only: batch_workspace, evaluate_batch
    use accelnet_batch_target, only: target_model, target_workspace, evaluate_batch_target
    use accelnet_predictor, only: predictor_model
    use accelnet_descriptors, only: atomic_structure, neighbor_data, build_neighbor_list
    use batch_test_support
    implicit none
    type(predictor_model) :: model
    type(target_model) :: packed
    type(target_workspace) :: gpuwork
    type(batch_workspace) :: cpuwork
    type(atomic_structure) :: s
    type(neighbor_data) :: nb
    real(real64), allocatable :: ec(:), eg(:), fc(:,:), fg(:,:)
    integer, allocatable :: centers(:)
    real(real64) :: wc(3,3), wg(3,3), seconds, elapsed, times(4), error
    integer(int64) :: start, finish, rate, repeats
    integer :: n, i, sample, order, method, j
    character(len=128) :: arg
    if (command_argument_count() /= 3) error stop 'usage: accelnet-target-benchmark NATOMS ORDER SECONDS'
    call get_command_argument(1,arg); read(arg,*) n
    call get_command_argument(2,arg); read(arg,*) order
    call get_command_argument(3,arg); read(arg,*) seconds
    if (n < 1 .or. order < 0 .or. seconds <= 0) error stop 'invalid benchmark parameters'
    call make_model('chebyshev',model,order=order)
    call make_structure(n,2,s)
    call packed%initialize(model)
    allocate(ec(n),eg(n),fc(3,n),fg(3,n),centers(n))
    centers = [(i,i=1,n)]
    call build_neighbor_list(s,model%maximum_cutoff,nb,model%minimum_distance)
    call evaluate(1); call evaluate(2)
    error = max(maxval(abs(ec-eg)),maxval(abs(fc-fg)),maxval(abs(wc-wg)))
    call check()
    write(*,'(A,1X,I0,1X,I0,1X,I0,1X,ES24.16)') 'CASE',n,order,size(nb%atom_indices),error
    do sample = 1, 5
        do j = 1, 4
            method = j
            if (mod(sample,2) == 0) method = 5-j
            call evaluate(method) ! warmup also excludes one-time allocations
            call system_clock(start,rate)
            if (rate <= 0) error stop 'wall clock unavailable'
            repeats = 0
            do
                call evaluate(method)
                repeats = repeats+1
                call system_clock(finish)
                elapsed = real(finish-start,real64)/real(rate,real64)
                if (elapsed >= seconds) exit
            end do
            times(method) = elapsed/real(repeats,real64)
            call check()
        end do
        write(*,'(A,1X,I0,4(1X,ES24.16))') 'TIMING',sample,times
    end do
contains
    subroutine evaluate(method)
        integer, intent(in) :: method
        if (method > 2) call build_neighbor_list(s,model%maximum_cutoff,nb,model%minimum_distance)
        if (mod(method,2) == 1) then
            fc = 0; wc = 0
            call evaluate_batch(model,s%species,centers,nb%offsets,nb%atom_indices,nb%displacements,ec,fc,cpuwork,wc)
        else
            fg = 0; wg = 0
            call evaluate_batch_target(packed,s%species,centers,nb%offsets,nb%atom_indices,nb%displacements,eg,fg,gpuwork,wg)
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
