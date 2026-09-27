! Without ACCELNET_HAVE_BATCH this driver can be linked against the pre-change
! CPU library. Keep fixtures and timing code identical for both builds.
program cpu_regression_benchmark
    use iso_fortran_env, only: real64, int64
    use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2
    use accelnet_descriptors, only: atomic_structure, neighbor_data, build_neighbor_list
    use batch_test_support
#ifdef ACCELNET_HAVE_BATCH
    use accelnet_batch, only: batch_workspace, evaluate_batch
    use accelnet_cpu_reference, only: reference_predict_energy_forces
#endif
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: s
    type(neighbor_data) :: neighbors
#ifdef ACCELNET_HAVE_BATCH
    type(batch_workspace) :: work
#endif
    real(real64), allocatable :: f(:,:), energies(:)
    integer, allocatable :: centers(:)
    real(real64) :: e, w(3,3), duration, seconds
    integer(int64) :: start, finish, rate, repeats
    integer :: natoms, i
    character(len=1024) :: family, directory, argument, mode
    if (command_argument_count() /= 5) &
        error stop 'usage: cpu-regression-benchmark FAMILY NATOMS SECONDS DATA_DIR structure|batch'
    call get_command_argument(1, family)
    call get_command_argument(2, argument); read(argument,*) natoms
    call get_command_argument(3, argument); read(argument,*) duration
    call get_command_argument(4, directory)
    call get_command_argument(5, mode)
    if (natoms < 1 .or. duration <= 0) error stop 'positive atom count and duration required'
    if (mode /= 'structure' .and. mode /= 'batch' .and. mode /= 'reference') error stop 'invalid evaluation mode'
    select case(trim(family))
    case('chebyshev', 'lj')
        call make_model(trim(family), model)
    case('n2p2-g5')
        call load_predictor_from_n2p2(trim(directory)//'/n2p2-per-element', model)
    case('n2p2-g4')
        call load_predictor_from_n2p2(trim(directory)//'/n2p2-virial-angular', model)
    case default
        error stop 'unknown benchmark family'
    end select
    call make_structure(natoms, size(model%networks), s)
    allocate(f(3,natoms), energies(natoms), centers(natoms))
    centers = [(i,i=1,natoms)]
    call evaluate() ! warm up code and allocations
    call system_clock(start, rate)
    if (rate <= 0) error stop 'wall clock unavailable'
    repeats = 0
    do
        call evaluate()
        repeats = repeats+1
        call system_clock(finish)
        seconds = real(finish-start,real64)/real(rate,real64)
        if (seconds >= duration) exit
    end do
    write(*,'(A,1X,ES24.16,1X,I0)') 'TIMING', seconds/real(repeats,real64), repeats
    write(*,'(A,1X,ES24.16)') 'ENERGY', e
    do i = 1, natoms
        write(*,'(A,3(1X,ES24.16))') 'FORCE', f(:,i)
    end do
    do i = 1, 3
        write(*,'(A,3(1X,ES24.16))') 'VIRIAL', w(:,i)
    end do
contains
    subroutine evaluate()
        if (mode == 'structure') then
            call model%predict_energy_forces(s, e, f, w)
#ifdef ACCELNET_HAVE_BATCH
        else if (mode == 'reference') then
            call reference_predict_energy_forces(model,s,e,f,w)
#endif
        else
#ifdef ACCELNET_HAVE_BATCH
            ! Include neighbor construction in both modes for a fair end-to-end
            ! CPU comparison. The reusable workspace persists across iterations.
            call build_neighbor_list(s, model%maximum_cutoff, neighbors, model%minimum_distance)
            f = 0; w = 0
            call evaluate_batch(model, s%species, centers, neighbors%offsets, neighbors%atom_indices, &
                neighbors%displacements, energies, f, work, w)
            e = sum(energies)
#else
            error stop 'baseline driver has no batch implementation'
#endif
        end if
    end subroutine
end program
