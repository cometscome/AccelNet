! Same model/structure inputs as the independent n2p2 benchmark.
program n2p2_candidate_benchmark
    use iso_fortran_env, only: real64,int64
    use accelnet_predictor, only: predictor_model,load_predictor_from_n2p2,load_predictor_from_networks
    use accelnet_descriptors, only: atomic_structure,neighbor_data,read_xsf,build_neighbor_list
    use accelnet_batch, only: batch_workspace,evaluate_batch,evaluate_batch_reference
    use accelnet_batch_target, only: target_model,target_workspace,evaluate_batch_target
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: s
    type(neighbor_data) :: nb
    type(batch_workspace) :: work,refwork
    type(target_model) :: packed
    type(target_workspace) :: gpuwork
    real(real64), allocatable :: e(:),f(:,:),er(:),fr(:,:)
    integer, allocatable :: centers(:)
    real(real64) :: w(3,3),wr(3,3),duration,elapsed
    integer(int64) :: start,finish,rate,count
    integer :: i,sample,unit,ios,nfiles,g5_mode
    character(len=4096), allocatable :: filenames(:)
    character(len=4096) :: element,filename
    character(len=4096) :: dir,path,arg,backend,scope,format
    if (command_argument_count() < 6 .or. command_argument_count() > 7) &
        error stop 'usage: candidate MODEL INPUT_XSF SECONDS cpu|gpu|host fixed|full n2p2|native [G5_MODE]'
    call get_command_argument(1,dir);call get_command_argument(2,path)
    call get_command_argument(3,arg);read(arg,*) duration
    call get_command_argument(4,backend);call get_command_argument(5,scope);call get_command_argument(6,format)
    if (duration < 0 .or. (scope /= 'fixed' .and. scope /= 'full')) error stop 'invalid timing arguments'
    if (backend /= 'cpu' .and. backend /= 'gpu' .and. backend /= 'host') error stop 'invalid backend'
    if (format == 'n2p2') then
        call load_predictor_from_n2p2(trim(dir),model)
    else if (format == 'native') then
        open(newunit=unit,file=trim(dir)//'/networks.list',status='old',action='read')
        nfiles=0
        do
            read(unit,*,iostat=ios) element,filename
            if (ios /= 0) exit
            nfiles=nfiles+1
        end do
        rewind(unit);allocate(filenames(nfiles))
        do i=1,nfiles
            read(unit,*) element,filename
            filenames(i)=trim(dir)//'/'//trim(filename)
        end do
        close(unit)
        call load_predictor_from_networks(filenames,model)
    else
        error stop 'invalid model format'
    end if
    g5_mode=0
    if (command_argument_count() == 7) then
        call get_command_argument(7,arg);read(arg,*) g5_mode
        call model%set_g5_evaluation(g5_mode)
    end if
    call read_xsf(trim(path),model%species_names,s)
    call packed%initialize(model,use_host=backend /= 'gpu',g5_mode=g5_mode)
    allocate(e(s%natoms),f(3,s%natoms),er(s%natoms),fr(3,s%natoms),centers(s%natoms))
    centers=[(i,i=1,s%natoms)]
    call build_neighbor_list(s,model%maximum_cutoff,nb,model%minimum_distance)
    fr=0;wr=0
    call evaluate_batch_reference(model,s%species,centers,nb%offsets,nb%atom_indices,nb%displacements,er,fr,refwork,wr)
    call evaluate()
    if (maxval(abs(e-er)) > 2e-9_real64 .or. maxval(abs(f-fr)) > 2e-9_real64 .or. &
        maxval(abs(w-wr)) > 2e-8_real64) error stop 'batch/reference mismatch'
    do sample=1,merge(5,0,duration > 0)
        call system_clock(start,rate);count=0
        do
            call evaluate();count=count+1
            call system_clock(finish);elapsed=real(finish-start,real64)/real(rate,real64)
            if (elapsed >= duration) exit
        end do
        write(*,'(A,1X,ES24.16,1X,I0)') 'TIMING',elapsed/real(count,real64),count
    end do
    write(*,'(A,1X,ES24.16)') 'ENERGY',sum(e)
    do i=1,s%natoms
        write(*,'(A,3(1X,ES24.16))') 'FORCE',f(:,i)
    end do
    do i=1,3
        write(*,'(A,3(1X,ES24.16))') 'VIRIAL',w(:,i)
    end do
contains
    subroutine evaluate()
        if (scope == 'full') call build_neighbor_list(s,model%maximum_cutoff,nb,model%minimum_distance)
        f=0;w=0
        if (backend == 'cpu') then
            call evaluate_batch(model,s%species,centers,nb%offsets,nb%atom_indices,nb%displacements,e,f,work,w)
        else
            call evaluate_batch_target(packed,s%species,centers,nb%offsets,nb%atom_indices,nb%displacements,e,f,gpuwork,w)
        end if
    end subroutine
end program
