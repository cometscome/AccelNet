! Same driver for archived/current target libraries; neighbor construction is excluded.
program target_energy_benchmark
    use iso_fortran_env, only: real64,int64
    use accelnet_predictor, only: predictor_model,load_predictor_from_networks
    use accelnet_descriptors, only: atomic_structure,neighbor_data,build_neighbor_list
    use accelnet_batch, only: batch_workspace,evaluate_batch_reference
    use accelnet_batch_target, only: target_model,target_workspace,target_profile,evaluate_batch_target
    use batch_test_support, only: make_model,make_structure,make_g5_scaling_neighbors
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: structure
    type(neighbor_data) :: neighbors
    type(batch_workspace) :: reference_work
    type(target_model) :: packed
    type(target_workspace) :: work
    type(target_profile) :: profile
    character(len=1024) :: arg,directory,files(2)
    integer :: n,order,mode,i
    integer, allocatable :: centers(:)
    integer(int64) :: started,finished,rate,repeats
    real(real64) :: seconds,elapsed,w(3,3),reference_w(3,3)
    real(real64), allocatable :: e(:),reference_e(:),f(:,:),reference_f(:,:)
    if (command_argument_count()/=5) error stop 'N ORDER SECONDS MODE MODEL_DIR_OR_DASH'
    call get_command_argument(1,arg); read(arg,*) n
    call get_command_argument(2,arg); read(arg,*) order
    call get_command_argument(3,arg); read(arg,*) seconds
    call get_command_argument(4,arg); read(arg,*) mode
    call get_command_argument(5,directory)
    if (n<1.or.order<0.or.seconds<=0.or.mode<0.or.mode>2) error stop 'invalid arguments'
    if (directory=='-') then
        call make_model('chebyshev',model,order=order)
    else
        files=[trim(directory)//'/Ti.nn.ascii',trim(directory)//'/O.nn.ascii']
        call load_predictor_from_networks(files,model)
    end if
    call model%set_chebyshev_evaluation(mode)
    call packed%initialize(model)
    call make_structure(n,size(model%networks),structure)
    if (directory=='-') then
        call make_g5_scaling_neighbors(n,64,neighbors)
    else
        call build_neighbor_list(structure,model%maximum_cutoff,neighbors,model%minimum_distance)
    end if
    allocate(centers(n),e(n),reference_e(n),f(3,n),reference_f(3,n))
    centers=[(i,i=1,n)]
    reference_f=0; reference_w=0
    call evaluate_batch_reference(model,structure%species,centers,neighbors%offsets, &
        neighbors%atom_indices,neighbors%displacements,reference_e,reference_f,reference_work,reference_w)
    f=0; w=0
    call evaluate_batch_target(packed,structure%species,centers,neighbors%offsets, &
        neighbors%atom_indices,neighbors%displacements,e,f,work,w)
    call check(e,reference_e)
    call check(reshape(f,[3*n]),reshape(reference_f,[3*n]))
    call check(reshape(w,[9]),reshape(reference_w,[9]))
    f=0.25_real64; w=0.5_real64
    call evaluate()
    call system_clock(started,rate); repeats=0
    do
        call evaluate()
        repeats=repeats+1
        call system_clock(finished)
        elapsed=real(finished-started,real64)/real(rate,real64)
        if (elapsed>=seconds) exit
    end do
    call check(e,reference_e)
    if (any(f/=0.25_real64).or.any(w/=0.5_real64)) error stop 'energy-only modified accumulators'
    write(*,'(A,ES24.16,1X,I0)') 'TIMING ',elapsed/repeats,repeats
    write(*,'(A,7ES24.16)') 'PROFILE ',profile%prepare,profile%upload,profile%descriptors, &
        profile%network,profile%forces,profile%download,profile%total
    write(*,'(A,ES24.16)') 'MAX_ENERGY_ERROR ',maxval(abs(e-reference_e))
contains
    subroutine evaluate()
        call evaluate_batch_target(packed,structure%species,centers,neighbors%offsets, &
            neighbors%atom_indices,neighbors%displacements,e,f,work,w,profile,energy_only=.true.)
    end subroutine
    subroutine check(actual,expected)
        real(real64), intent(in) :: actual(:),expected(:)
        if (.not.all(abs(actual-expected)<=2e-10_real64*(1+abs(expected)))) error stop 'reference mismatch'
    end subroutine
end program
