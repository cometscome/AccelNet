program write_target_fixtures
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model
    use aenet_network, only: write_aenet_network_ascii
    use batch_test_support, only: make_model
    implicit none
    type(predictor_model) :: model
    integer :: family,s,d,j,b
    character(len=1024) :: directory,filename
    character(len=8) :: name
    call get_command_argument(1,directory)
    do family = 1,2
        name = 'lj'
        if (family == 2) name = 'behler'
        call make_model(trim(name),model)
        do s = 1,2
            associate(net => model%networks(s))
            d = net%nodes(1)
            net%species_names = model%species_names
            net%environment_names = model%species_names
            if (s == 2) net%environment_names = model%species_names([2,1])
            net%minimum_radius = model%minimum_distance; net%maximum_radius = model%maximum_cutoff
            net%description = 'Synthetic GPU descriptor regression (not a physical potential)'
            allocate(net%descriptor_kinds(d),net%descriptor_parameters(7,d),net%descriptor_environments(2,d))
            net%descriptor_kinds = 1; net%descriptor_parameters = 0; net%descriptor_environments = 1
            net%descriptor_parameters(5,:) = 1.0_real64 ! serialized cosine cutoff
            if (family == 1) then
                net%descriptor_name = 'LJ'
                net%descriptor_parameters(1,:) = 3.4_real64
            else
                net%descriptor_name = 'Behler2011'
                associate(cfg => model%setups(s)%model%behler(1)%config)
                do j = 1,size(cfg%g1)
                    associate(p => cfg%g1(j))
                    b = p%output
                    net%descriptor_kinds(b) = 1
                    net%descriptor_parameters(1,b) = p%rc
                    net%descriptor_environments(1,b) = p%species
                    end associate
                end do
                do j = 1,size(cfg%g2)
                    associate(p => cfg%g2(j))
                    b = p%output
                    net%descriptor_kinds(b) = 2
                    net%descriptor_parameters(1,b) = p%rc
                    net%descriptor_environments(1,b) = p%species
                    net%descriptor_parameters(2:3,b) = [p%rs,p%eta]
                    end associate
                end do
                do j = 1,size(cfg%g3)
                    associate(p => cfg%g3(j))
                    b = p%output
                    net%descriptor_kinds(b) = 3
                    net%descriptor_parameters(1,b) = p%rc
                    net%descriptor_environments(1,b) = p%species
                    net%descriptor_parameters(2,b) = p%kappa
                    end associate
                end do
                do j = 1,size(cfg%g4)
                    associate(p => cfg%g4(j))
                    b = p%output
                    net%descriptor_kinds(b) = 4
                    net%descriptor_parameters(1,b) = p%rc
                    net%descriptor_environments(:,b) = [p%species1,p%species2]
                    net%descriptor_parameters(2:4,b) = [p%lambda,p%zeta,p%eta]
                    net%descriptor_parameters(7,b) = p%rs
                    end associate
                end do
                do j = 1,size(cfg%g5)
                    associate(p => cfg%g5(j))
                    b = p%output
                    net%descriptor_kinds(b) = 5
                    net%descriptor_parameters(1,b) = p%rc
                    net%descriptor_environments(:,b) = [p%species1,p%species2]
                    net%descriptor_parameters(2:4,b) = [p%lambda,p%zeta,p%eta]
                    net%descriptor_parameters(7,b) = p%rs
                    end associate
                end do
                end associate
            end if
            filename = trim(directory)//'/'//trim(net%atomtype)//'.'//trim(name)//'.nn.ascii'
            call write_aenet_network_ascii(trim(filename),net)
            end associate
        end do
    end do
end program
