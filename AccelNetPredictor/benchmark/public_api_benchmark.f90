! The same driver is compiled against archived and current libraries.
program public_api_benchmark
    use iso_fortran_env, only: real64,int64
    use accelnet_predictor, only: predictor_model,load_predictor_from_networks,load_predictor_from_n2p2
    use accelnet_descriptors, only: atomic_structure,neighbor_data,build_neighbor_list
    use accelnet_batch, only: batch_workspace,evaluate_batch
    use accelnet, only: accelnet_init,accelnet_init_n2p2,accelnet_load_potential,accelnet_final, &
        accelnet_atomic_energy,accelnet_atomic_energy_and_forces_virial
    use batch_test_support, only: make_structure,make_model
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: structure
    type(neighbor_data) :: neighbors
    type(batch_workspace) :: work
    character(len=1024) :: family,directory,path,arg,files(2)
    real(real64), allocatable :: forces(:,:),energies(:)
    integer, allocatable :: centers(:)
    real(real64) :: total,virial(3,3),duration,seconds
    integer(int64) :: start,finish,rate,repeats
    integer :: natoms,i,stat
    if (command_argument_count()/=5) error stop 'FAMILY MODEL_DIR NATOMS SECONDS PATH'
    call get_command_argument(1,family); call get_command_argument(2,directory)
    call get_command_argument(3,arg); read(arg,*) natoms
    call get_command_argument(4,arg); read(arg,*) duration
    call get_command_argument(5,path)
    if (natoms<1.or.duration<=0) error stop 'invalid dimensions or duration'
    select case(trim(family))
    case('aenet')
        files=[trim(directory)//'/Ti.nn.ascii',trim(directory)//'/O.nn.ascii']
        call load_predictor_from_networks(files,model)
        call accelnet_init(model%species_names,stat); call checked()
        do i=1,2
            call accelnet_load_potential(i,trim(files(i)),stat,is_ascii=.true.); call checked()
        end do
    case('n2p2')
        call load_predictor_from_n2p2(trim(directory),model)
        call accelnet_init_n2p2(trim(directory),stat); call checked()
    case('combined','multi-chebyshev','mixed-components')
        call make_model(trim(family),model)
        if (index(path,'atomic')==1) error stop 'synthetic composite needs object/batch API'
    case default
        error stop 'invalid family'
    end select
    call make_structure(natoms,size(model%networks),structure)
    call build_neighbor_list(structure,model%maximum_cutoff,neighbors,model%minimum_distance)
    allocate(forces(3,natoms),energies(natoms),centers(natoms))
    centers=[(i,i=1,natoms)]
    call evaluate()
    call system_clock(start,rate); repeats=0
    do
        call evaluate()
        repeats=repeats+1
        call system_clock(finish)
        seconds=real(finish-start,real64)/real(rate,real64)
        if (seconds>=duration) exit
    end do
    write(*,'(A,1X,ES24.16,1X,I0)') 'TIMING',seconds/repeats,repeats
    write(*,'(A,1X,ES24.16)') 'ENERGY',total
    do i=1,natoms
        write(*,'(A,3(1X,ES24.16))') 'FORCE',forces(:,i)
    end do
    do i=1,3
        write(*,'(A,3(1X,ES24.16))') 'VIRIAL',virial(i,:)
    end do
    if (family=='aenet'.or.family=='n2p2') then
        call accelnet_final(stat); call checked()
    end if
contains
    subroutine checked()
        if (stat/=0) error stop 'public API error'
    end subroutine
    subroutine evaluate()
        integer :: row,first,last,n
        real(real64) :: center(3)
        forces=0; virial=0; center=0; total=0
        select case(trim(path))
        case('structure')
            call model%predict_energy_forces(structure,total,forces,virial)
        case('energy')
            call model%predict_energy(structure,total)
        case('batch')
            call evaluate_batch(model,structure%species,centers,neighbors%offsets,neighbors%atom_indices, &
                neighbors%displacements,energies,forces,work,virial)
            total=sum(energies)
        case('atomic','atomic-energy')
            do row=1,natoms
                first=neighbors%offsets(row); last=neighbors%offsets(row+1)-1; n=last-first+1
                if (path=='atomic') then
                    call accelnet_atomic_energy_and_forces_virial(center,structure%species(row),row,n, &
                        neighbors%displacements(:,first:last),structure%species(neighbors%atom_indices(first:last)), &
                        neighbors%atom_indices(first:last),natoms,energies(row),forces,virial,stat)
                else
                    call accelnet_atomic_energy(center,structure%species(row),n,neighbors%displacements(:,first:last), &
                        structure%species(neighbors%atom_indices(first:last)),energies(row),stat)
                end if
                call checked()
            end do
            total=sum(energies)
        case default
            error stop 'invalid API path'
        end select
    end subroutine
end program
