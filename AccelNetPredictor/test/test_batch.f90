program test_batch
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use accelnet_batch, only: batch_workspace, evaluate_batch, evaluate_batch_reference
    use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2, load_predictor_from_networks
    use accelnet_descriptors, only: atomic_structure, neighbor_data, build_neighbor_list, read_xsf
    use batch_test_support
    implicit none
    type(predictor_model) :: model
    type(batch_workspace) :: work, reference_work
    type(atomic_structure) :: s
    character(len=1024) :: directory, argument, network_files(2)
    integer :: family, version, mode, geometry
    character(len=16), parameter :: generic_families(7) = [character(len=16) :: &
        'g4-distinct', 'g4-series', 'g5', 'g5-series', 'behler', 'lj-behler', 'lj']
    call get_command_argument(1, directory)
    call get_command_argument(2, argument)
    if (argument == '--aenet') then
        call get_command_argument(3, directory)
        network_files(1) = trim(directory)//'/Ti.nn.ascii'
        network_files(2) = trim(directory)//'/O.nn.ascii'
        family = 0
        geometry = 0
        do version = 0, 2
            call load_predictor_from_networks(network_files, model, &
                chebyshev_version=merge(10,version,version == 2))
            call read_xsf(trim(directory)//'/structure2935.xsf', model%species_names, s)
            do mode = 0, 2
                call model%set_chebyshev_evaluation(mode)
                call check(s)
            end do
        end do
        print *, 'Batch real aenet model equivalence passed'
        stop
    end if
    if (len_trim(argument) > 0) call invalid_input(trim(argument))

    do family = 1, 3
        do version = 0, 2
            select case(family)
            case(1)
                call make_model('chebyshev', model, order=3+version, version=merge(10,version,version == 2))
            case(2)
                call make_model('lj', model)
            case(3)
                call make_model('combined', model)
            end select
            do mode = 0, 2
                call model%set_chebyshev_evaluation(mode)
                do geometry = 1, 3
                    call make_structure(8, 2, s)
                    s%pbc = geometry /= 1
                    if (geometry == 3) s%lattice(1,2) = 0.4_real64
                    call check(s)
                end do
            end do
        end do
    end do
    ! Check the default common path and forced-G5 fallback with one workspace,
    ! including partitioned/reordered rows and additive forces/virials.
    do family = 1, size(generic_families)
        call make_model(trim(generic_families(family)), model, order=8)
        do mode = 0, 3
            call model%set_g5_evaluation(mode)
            do geometry = 1, 3
                call make_structure(8, 2, s)
                s%pbc = geometry /= 1
                if (geometry == 3) s%lattice(1,2) = 0.4_real64
                call check(s)
            end do
        end do
    end do
    ! Reuse the same workspace across model reloads, species counts, topology,
    ! descriptor sizes, G4 Jacobians, and G5 contraction/moment modes.
    do family = 1, 4
        select case(family)
        case(1)
            call load_predictor_from_n2p2(trim(directory)//'/n2p2', model)
        case(2)
            call load_predictor_from_n2p2(trim(directory)//'/n2p2-per-element', model)
        case(3)
            call load_predictor_from_n2p2(trim(directory)//'/n2p2-per-element-depth', model)
        case(4)
            call load_predictor_from_n2p2(trim(directory)//'/n2p2-virial-angular', model)
        end select
        do mode = 1, 3
            call model%set_g5_evaluation(mode)
            do geometry = 1, 3
                call make_structure(8, size(model%networks), s)
                s%pbc = geometry /= 1
                if (geometry == 3) s%lattice(1,2) = 0.4_real64
                call check(s)
            end do
        end do
    end do
    call make_model('chebyshev', model, order=7)
    call make_structure(64, 2, s)
    call check(s)
    call make_structure(1, 2, s)
    call check(s) ! repeated periodic self images: nonzero virial, zero force
    s%pbc = .false.
    call check(s) ! zero neighbors
    call work%release()
    call require(work%allocations() == 0, 'release resets workspace')
    call check(s)
    print *, 'Batch energies, forces, virials, partitions and workspace reuse passed'
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not. ok) then
            print *, message
            error stop 1
        end if
    end subroutine

    subroutine close_array(a, b, message)
        real(real64), intent(in) :: a(:), b(:)
        character(len=*), intent(in) :: message
        call require(all(ieee_is_finite(a)), message//' finite')
        call require(all(abs(a-b) <= 2e-10_real64 + 2e-10_real64*abs(b)), message)
    end subroutine

    subroutine check(s)
        type(atomic_structure), intent(in) :: s
        type(atomic_structure) :: shifted
        type(neighbor_data) :: neighbors
        real(real64) :: e, ep, em, f(3,s%natoms), fb(3,s%natoms), w(3,3), wb(3,3)
        real(real64) :: energies(s%natoms), split_energies(s%natoms), h
        integer :: centers(s%natoms), k, allocations, split
        centers = [(k,k=1,s%natoms)]
        call model%predict_energy_forces(s, e, f, w)
        call build_neighbor_list(s, model%maximum_cutoff, neighbors, model%minimum_distance)
        fb=0; wb=0
        call evaluate_batch_reference(model,s%species,centers,neighbors%offsets,neighbors%atom_indices, &
            neighbors%displacements,energies,fb,reference_work,wb)
        call close_array([sum(energies)],[e],'structure/reference energy')
        call close_array(reshape(fb,[3*s%natoms]),reshape(f,[3*s%natoms]),'structure/reference force')
        call close_array(reshape(wb,[9]),reshape(w,[9]),'structure/reference virial')
        fb = 0; wb = 0
        call evaluate_batch(model, s%species, centers, neighbors%offsets, neighbors%atom_indices, &
            neighbors%displacements, energies, fb, work, wb)
        call close_array([sum(energies)], [e], 'batch energy')
        call close_array(reshape(fb,[3*s%natoms]), reshape(f,[3*s%natoms]), 'batch force')
        call close_array(reshape(wb,[9]), reshape(w,[9]), 'batch virial')
        allocations = work%allocations()

        ! Partition CSR rows without repacking edges. Contributions to atoms
        ! outside each central subset (ghosts in MD) must still be accumulated.
        split = s%natoms/2
        fb = 0.25_real64; wb = 0.5_real64
        call evaluate_batch(model, s%species, centers(:split), neighbors%offsets(:split+1), neighbors%atom_indices, &
            neighbors%displacements, split_energies(:split), fb, work, wb)
        call evaluate_batch(model, s%species, centers(split+1:), neighbors%offsets(split+1:), neighbors%atom_indices, &
            neighbors%displacements, split_energies(split+1:), fb, work, wb)
        call close_array(split_energies, energies, 'partition atomic energy')
        call close_array(reshape(fb-0.25_real64,[3*s%natoms]), reshape(f,[3*s%natoms]), 'partition additive force')
        call close_array(reshape(wb-0.5_real64,[9]), reshape(w,[9]), 'partition additive virial')
        call require(work%allocations() == allocations, 'workspace must not grow on reuse')
        ! Reversed central order tests the row-to-atom mapping independently.
        fb = 0
        do k = s%natoms, 1, -1
            call evaluate_batch(model, s%species, centers(k:k), neighbors%offsets(k:k+1), neighbors%atom_indices, &
                neighbors%displacements, split_energies(k:k), fb, work)
        end do
        call close_array(split_energies, energies, 'reordered atomic energy')
        call close_array(reshape(fb,[3*s%natoms]), reshape(f,[3*s%natoms]), 'reordered force without virial')
        ! Empty batch is a no-op for accumulators, including no atoms at all.
        call evaluate_batch(model, [integer::], [integer::], [1], [integer::], &
            neighbors%displacements(:,:0), energies(:0), fb(:,:0), work, wb)

        h = 1e-5_real64
        shifted = s
        shifted%positions(1,1) = s%positions(1,1)+h
        call model%predict_energy(shifted, ep)
        shifted%positions(1,1) = s%positions(1,1)-h
        call model%predict_energy(shifted, em)
        if (abs(f(1,1)+(ep-em)/(2*h)) >= 2e-6_real64*max(1.0_real64,abs(f(1,1)))) &
            print *, 'FD case family/version/mode/geometry/N/f/FD:', family,version,mode,geometry,s%natoms, &
                f(1,1), -(ep-em)/(2*h)
        call require(abs(f(1,1)+(ep-em)/(2*h)) < 2e-6_real64*max(1.0_real64,abs(f(1,1))), &
            'finite difference force')
    end subroutine

    subroutine invalid_input(which)
        character(len=*), intent(in) :: which
        integer :: species(2), centers(2), offsets(3), indices(2)
        real(real64) :: dr(3,2), e(2), f(3,2)
        call make_model('chebyshev', model)
        species = [1,2]; centers = [1,2]; offsets = [1,2,3]; indices = [2,1]
        dr(:,1) = [1.0_real64,0.0_real64,0.0_real64]; dr(:,2) = -dr(:,1); f = 0
        select case(which)
        case('offset'); offsets(1) = 0
        case('monotonic'); offsets = [1,3,2]
        case('target'); indices(1) = 3
        case('center'); centers(1) = 0
        case('species'); species(2) = 3
        case('shape')
            call evaluate_batch(model, species, centers, offsets(:2), indices, dr, e, f, work)
            stop 0
        end select
        call evaluate_batch(model, species, centers, offsets, indices, dr, e, f, work)
        stop 0 ! CTest WILL_FAIL must reject a silently accepted input.
    end subroutine
end program
