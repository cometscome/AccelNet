! Serialized models and independent legacy-CPU reference values for C callers.
program write_target_c_loading_fixtures
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model, load_predictor_from_networks, load_predictor_from_n2p2
    use accelnet_batch, only: batch_workspace, evaluate_batch_reference
    use aenet_network, only: write_aenet_network_ascii
    use batch_test_support, only: make_model
    implicit none
    type(predictor_model) :: model, loaded
    type(batch_workspace) :: work
    character(len=4096) :: directory, data, paths(2), filename
    character(len=40), parameter :: fixtures(4) = [character(len=40) :: &
        'n2p2', 'n2p2-per-element', 'n2p2-per-element-depth', 'n2p2-virial-angular']
    integer :: s,v,j,k,edge,unit,nspecies,species(4),centers(4),offsets(5),indices(12)
    real(real64) :: positions(3,4),dr(3,12),energies(4),forces(3,4),virial(3,3)
    call get_command_argument(1,directory)
    call get_command_argument(2,data)
    call make_model('chebyshev',model,order=3)
    do s=1,2
        associate(net => model%networks(s))
        net%species_names=model%species_names
        net%environment_names=model%species_names
        if (s==2) net%environment_names=model%species_names([2,1])
        net%minimum_radius=model%minimum_distance
        net%maximum_radius=model%maximum_cutoff
        net%description='Synthetic C API version test (not a physical potential)'
        net%descriptor_name='Chebyshev'
        allocate(net%descriptor_kinds(net%nodes(1)),net%descriptor_environments(2,net%nodes(1)))
        net%descriptor_kinds=1; net%descriptor_environments=1
        net%descriptor_parameters=spread([3.4_real64,3.0_real64,3.4_real64,3.0_real64],2,net%nodes(1))
        paths(s)=trim(directory)//'/'//trim(net%atomtype)//'.cheb.nn.ascii'
        call write_aenet_network_ascii(trim(paths(s)),net)
        end associate
    end do
    positions=reshape([0.0_real64,0.0_real64,0.0_real64,1.1_real64,0.2_real64,0.1_real64, &
        0.3_real64,1.3_real64,-0.2_real64,1.5_real64,1.1_real64,0.4_real64],[3,4])
    centers=[1,2,3,4]; offsets=[1,4,7,10,13]
    edge=0
    do j=1,4
        do k=1,4
            if (j==k) cycle
            edge=edge+1
            indices(edge)=k
            dr(:,edge)=positions(:,k)-positions(:,j)
        end do
    end do
    do v=0,2
        call load_predictor_from_networks(paths,loaded,chebyshev_version=merge(10,v,v==2))
        write(filename,'(a,"/cheb-",i0,".ref")') trim(directory),merge(10,v,v==2)
        call write_reference()
    end do
    do j=1,size(fixtures)
        call load_predictor_from_n2p2(trim(data)//'/'//trim(fixtures(j)),loaded)
        filename=trim(directory)//'/'//trim(fixtures(j))//'.ref'
        call write_reference()
    end do
contains
    subroutine write_reference()
        integer :: atom
        nspecies=size(loaded%species_names)
        species=[(1+mod(atom-1,nspecies),atom=1,4)]
        forces=0; virial=0
        call evaluate_batch_reference(loaded,species,centers,offsets,indices,dr,energies,forces,work,virial)
        open(newunit=unit,file=trim(filename),status='replace')
        write(unit,*) nspecies
        write(unit,'(a)') loaded%species_names
        write(unit,'(es26.17)') loaded%maximum_cutoff,energies,forces,virial
        close(unit)
    end subroutine
end program
