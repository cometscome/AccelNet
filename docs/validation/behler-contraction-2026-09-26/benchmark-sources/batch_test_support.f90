! Deterministic, self-contained models shared by correctness and timing drivers.
! Uses only APIs available before the batch implementation so the timing driver
! can also be linked against an archived CPU library.
module batch_test_support
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model
    use accelnet_descriptors, only: atomic_structure, descriptor_config, initialize_config
    use accelnet_lj, only: lj_config, initialize_lj_config
    use accelnet_descriptor_models, only: add_chebyshev, add_lj, add_behler
    use accelnet_behler, only: behler_config, initialize_behler_config, add_g1, add_g2, add_g3, add_g4, add_g5
    implicit none
contains
    subroutine make_model(family, model, order, version)
        character(len=*), intent(in) :: family
        type(predictor_model), intent(out) :: model
        integer, intent(in), optional :: order, version
        type(descriptor_config) :: chebyshev
        type(lj_config) :: lj
        type(behler_config) :: behler
        integer :: s, d, l, nw, j, degree, v, pair, sign, t1, t2
        degree = 3
        v = 0
        if (present(order)) degree = order
        if (present(version)) v = version
        model%species_names = [character(len=16) :: 'H', 'O']
        allocate(model%setups(2), model%networks(2))
        model%maximum_cutoff = 3.4_real64
        model%minimum_distance = 0.35_real64
        do s = 1, 2
            model%setups(s)%global_to_local = [1,2]
            if (s == 2) model%setups(s)%global_to_local = [2,1]
            select case (family)
            case ('chebyshev', 'combined')
                call initialize_config(chebyshev, 2, 3.4_real64, degree, 3.4_real64, degree, &
                    version=v, central_type_index=s)
                call add_chebyshev(model%setups(s)%model, chebyshev)
            case ('g4-series', 'g5-series')
                call initialize_behler_config(behler,2)
                do pair = 1,3
                    t1 = merge(2,1,pair == 3); t2 = merge(1,2,pair == 1)
                    do sign = 1,2
                        do j = 1,max(2,degree)
                            if (family == 'g4-series') then
                                call add_g4(behler,t1,t2,3.4_real64,real(2*sign-3,real64), &
                                    real(j,real64),0.2_real64,0.4_real64)
                            else
                                call add_g5(behler,t1,t2,3.4_real64,real(2*sign-3,real64), &
                                    real(j,real64),0.2_real64,0.4_real64)
                            end if
                        end do
                    end do
                end do
                call add_behler(model%setups(s)%model,behler)
            case ('behler', 'g4', 'g5', 'lj-behler')
                call initialize_behler_config(behler,2)
                do j = 1,2
                    if (family == 'behler' .or. family == 'lj-behler') then
                        call add_g1(behler,j,3.2_real64)
                        call add_g2(behler,j,3.4_real64,0.3_real64,0.6_real64)
                        call add_g3(behler,j,3.1_real64,1.3_real64)
                    end if
                    if (family /= 'g5') then
                        call add_g4(behler,j,j,3.4_real64,-1.0_real64,2.0_real64,0.15_real64,0.2_real64)
                        call add_g4(behler,1,2,3.3_real64,1.0_real64,1.5_real64,0.12_real64,0.1_real64)
                    end if
                    if (family /= 'g4') then
                        call add_g5(behler,j,j,3.4_real64,1.0_real64,real(j,real64),0.2_real64,0.4_real64)
                        call add_g5(behler,1,2,3.0_real64,-1.0_real64,2.5_real64,0.1_real64,0.2_real64)
                    end if
                end do
                call add_behler(model%setups(s)%model,behler)
            case ('lj')
            case default
                error stop 'unknown synthetic model'
            end select
            if (family == 'lj' .or. family == 'combined' .or. family == 'lj-behler') then
                ! Smooth cutoff: finite differences must not cross a hard jump.
                call initialize_lj_config(lj, 2, 3.4_real64, cutoff_type=1)
                call add_lj(model%setups(s)%model, lj)
            end if
            d = model%setups(s)%model%num_descriptors()
            model%networks(s)%atomtype = model%species_names(s)
            model%networks(s)%nlayers = 4
            model%networks(s)%nodes = [d,8,4,1]
            model%networks(s)%maxnodes = max(d,8)
            model%networks(s)%activation = [1,1,0]
            allocate(model%networks(s)%weight_offsets(4))
            nw = 0
            do l = 1, 3
                model%networks(s)%weight_offsets(l) = nw
                nw = nw + (model%networks(s)%nodes(l)+1)*model%networks(s)%nodes(l+1)
            end do
            model%networks(s)%weight_offsets(4) = nw
            model%networks(s)%weights = [(0.08_real64*sin(real(j+3*s,real64)), j=1,nw)]
            model%networks(s)%descriptor_shift = [(0.03_real64*j, j=1,d)]
            model%networks(s)%descriptor_scale = [(0.2_real64+0.01_real64*j, j=1,d)]
            model%networks(s)%energy_scale = 2.5_real64
            model%networks(s)%energy_shift = 0.7_real64
            model%networks(s)%atomic_references = [0.3_real64,-0.2_real64]
        end do
    end subroutine

    subroutine make_structure(natoms, ntypes, s)
        integer, intent(in) :: natoms, ntypes
        type(atomic_structure), intent(out) :: s
        integer :: i, side
        s%natoms = natoms
        side = ceiling(real(max(1,natoms),real64)**(1.0_real64/3.0_real64))
        do while (side**3 < natoms)
            side = side+1
        end do
        s%pbc = .true.
        s%lattice = 0.0_real64
        do i = 1, 3
            s%lattice(i,i) = 1.7_real64*side
        end do
        allocate(s%positions(3,natoms), s%species(natoms))
        do i = 1, natoms
            s%positions(:,i) = 1.7_real64*real([mod(i-1,side),mod((i-1)/side,side),(i-1)/side**2],real64) &
                + 0.07_real64*sin(real([i,2*i,3*i],real64))
            s%species(i) = 1+mod(i-1,ntypes)
        end do
    end subroutine
end module
