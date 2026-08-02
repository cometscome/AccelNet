module model_conversion
    use iso_fortran_env, only: real64
    use aenet_network, only: atomic_network, read_aenet_network, write_aenet_network_ascii
    use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2
    use n2p2_network, only: atomic_number
    implicit none
    private
    integer, parameter :: PATH_LENGTH = 4096
    public :: convert_n2p2_to_accelnet, convert_accelnet_to_n2p2

contains

    subroutine convert_n2p2_to_accelnet(input_directory, output_directory)
        character(len=*), intent(in) :: input_directory, output_directory
        type(predictor_model) :: model
        integer :: i, unit
        character(len=PATH_LENGTH) :: path

        call make_directory(output_directory)
        call load_predictor_from_n2p2(input_directory, model)
        do i = 1, size(model%networks)
            call reorder_descriptors(model%networks(i), .true.)
            path = join_path(output_directory, trim(model%networks(i)%atomtype)//".nn.ascii")
            call write_aenet_network_ascii(trim(path), model%networks(i))
            write(*, "(A)") "WROTE "//trim(path)
        end do
        path = join_path(output_directory, "networks.list")
        open(newunit=unit, file=trim(path), status="replace", action="write")
        do i = 1, size(model%networks)
            write(unit, "(A,1X,A)") trim(model%networks(i)%atomtype), &
                trim(model%networks(i)%atomtype)//".nn.ascii"
        end do
        close(unit)
        write(*, "(A)") "WROTE "//trim(path)
    end subroutine convert_n2p2_to_accelnet

    subroutine convert_accelnet_to_n2p2(output_directory, network_files)
        character(len=*), intent(in) :: output_directory
        character(len=*), intent(in) :: network_files(:)
        type(atomic_network), allocatable :: networks(:)
        integer, allocatable :: order(:)
        integer :: i

        if (size(network_files) == 0) error stop "at least one AccelNet network is required"
        allocate(networks(size(network_files)))
        do i = 1, size(network_files)
            call read_aenet_network(trim(network_files(i)), networks(i))
        end do
        call validate_networks(networks, order)
        do i = 1, size(networks)
            call reorder_descriptors(networks(i), .false.)
        end do
        call make_directory(output_directory)
        call write_n2p2_settings(output_directory, networks, order)
        call write_n2p2_scaling(output_directory, networks, order)
        do i = 1, size(order)
            call write_n2p2_weights(output_directory, networks(order(i)))
        end do
    end subroutine convert_accelnet_to_n2p2

    subroutine validate_networks(networks, order)
        type(atomic_network), intent(in) :: networks(:)
        integer, allocatable, intent(out) :: order(:)
        integer :: i, j, candidate

        if (size(networks) /= size(networks(1)%species_names)) &
            error stop "one network is required for every embedded species"
        allocate(order(size(networks)))
        do i = 1, size(networks)
            order(i) = i
            if (trim(lowercase(networks(i)%descriptor_name)) /= "behler2011") &
                error stop "only Behler2011 AccelNet networks can be converted to n2p2"
            if (any(networks(i)%descriptor_kinds /= 2 .and. &
                    networks(i)%descriptor_kinds /= 4 .and. networks(i)%descriptor_kinds /= 5)) &
                error stop "only Behler G2/G4/G5 descriptors can be converted to n2p2"
            if (any(networks(i)%nodes /= networks(1)%nodes) .or. &
                any(networks(i)%activation /= networks(1)%activation)) &
                error stop "n2p2 conversion requires a common network topology"
            if (any(networks(i)%species_names /= networks(1)%species_names) .or. &
                any(networks(i)%atomic_references /= networks(1)%atomic_references) .or. &
                networks(i)%energy_scale /= networks(1)%energy_scale .or. &
                networks(i)%energy_shift /= networks(1)%energy_shift) &
                error stop "AccelNet networks have inconsistent global metadata"
            if (networks(i)%descriptor_cutoff_type /= networks(1)%descriptor_cutoff_type .or. &
                networks(i)%descriptor_cutoff_alpha /= networks(1)%descriptor_cutoff_alpha) &
                error stop "AccelNet networks use inconsistent descriptor cutoffs"
        end do
        do i = 2, size(order)
            candidate = order(i); j = i - 1
            do while (j >= 1)
                if (atomic_number(trim(networks(order(j))%atomtype)) < &
                    atomic_number(trim(networks(candidate)%atomtype))) exit
                order(j + 1) = order(j); j = j - 1
            end do
            order(j + 1) = candidate
        end do
        do i = 1, size(order)
            if (atomic_number(trim(networks(order(i))%atomtype)) == 0) &
                error stop "unknown element in AccelNet network"
            if (i > 1) then
                if (trim(networks(order(i))%atomtype) == trim(networks(order(i - 1))%atomtype)) &
                    error stop "duplicate element network"
            end if
        end do
    end subroutine validate_networks

    subroutine reorder_descriptors(network, aenet_order)
        type(atomic_network), intent(inout) :: network
        logical, intent(in) :: aenet_order
        integer, allocatable :: permutation(:), kinds(:), environments(:, :)
        real(real64), allocatable :: parameters(:, :), shifts(:), scales(:), weights(:)
        integer :: i, j, candidate, nout, old_index, new_index

        allocate(permutation(size(network%descriptor_kinds)))
        permutation = [(i, i=1, size(permutation))]
        do i = 2, size(permutation)
            candidate = permutation(i); j = i - 1
            do while (j >= 1)
                if (.not. descriptor_less(network, candidate, permutation(j), aenet_order)) exit
                permutation(j + 1) = permutation(j); j = j - 1
            end do
            permutation(j + 1) = candidate
        end do
        if (all(permutation == [(i, i=1, size(permutation))])) return

        kinds = network%descriptor_kinds(permutation)
        parameters = network%descriptor_parameters(:, permutation)
        environments = network%descriptor_environments(:, permutation)
        shifts = network%descriptor_shift(permutation)
        scales = network%descriptor_scale(permutation)
        weights = network%weights
        nout = network%nodes(2)
        do new_index = 1, size(permutation)
            old_index = permutation(new_index)
            weights((new_index - 1)*nout + 1:new_index*nout) = &
                network%weights((old_index - 1)*nout + 1:old_index*nout)
        end do
        network%descriptor_kinds = kinds
        network%descriptor_parameters = parameters
        network%descriptor_environments = environments
        network%descriptor_shift = shifts
        network%descriptor_scale = scales
        network%weights = weights
    end subroutine reorder_descriptors

    logical function descriptor_less(network, ia, ib, aenet_order) result(less)
        type(atomic_network), intent(in) :: network
        integer, intent(in) :: ia, ib
        logical, intent(in) :: aenet_order
        integer :: ka(8), kb(8)
        real(real64) :: ra(6), rb(6)

        call descriptor_key(network, ia, aenet_order, ka, ra)
        call descriptor_key(network, ib, aenet_order, kb, rb)
        if (aenet_order) then
            if (any(ka(1:4) /= kb(1:4))) then
                less = integer_key_less(ka(1:4), kb(1:4)); return
            end if
            if (ka(5) /= kb(5)) then
                less = ka(5) < kb(5); return
            end if
            if (any(ra /= rb)) then
                less = real_key_less(ra, rb); return
            end if
            less = integer_key_less(ka(6:7), kb(6:7))
        else
            if (ka(1) /= kb(1)) then
                less = ka(1) < kb(1); return
            end if
            if (any(ra /= rb)) then
                less = real_key_less(ra, rb); return
            end if
            less = integer_key_less(ka(2:3), kb(2:3))
        end if
    end function descriptor_less

    subroutine descriptor_key(network, i, aenet_order, integers, reals)
        type(atomic_network), intent(in) :: network
        integer, intent(in) :: i
        logical, intent(in) :: aenet_order
        integer, intent(out) :: integers(8)
        real(real64), intent(out) :: reals(6)
        integer :: kind, first, second, n2kind

        kind = network%descriptor_kinds(i)
        first = network%descriptor_environments(1, i)
        second = network%descriptor_environments(2, i)
        select case(kind)
        case(2); n2kind = 2
        case(4); n2kind = 3
        case(5); n2kind = 9
        case default; error stop "unsupported Behler descriptor kind"
        end select
        integers = 0
        if (aenet_order) then
            if (kind == 2) then
                integers(1:4) = [first, 0, 0, kind]
            else
                integers(1:4) = [min(first, second), 1, max(first, second), kind]
            end if
            integers(5:7) = [n2kind, first, second]
        else
            integers(1:3) = [n2kind, first, second]
        end if
        if (kind == 2) then
            reals = [network%descriptor_parameters(1, i), network%descriptor_parameters(3, i), &
                network%descriptor_parameters(2, i), 0.0_real64, 0.0_real64, 0.0_real64]
        else
            reals = [network%descriptor_parameters(1, i), network%descriptor_parameters(4, i), &
                0.0_real64, network%descriptor_parameters(3, i), &
                network%descriptor_parameters(2, i), 0.0_real64]
        end if
    end subroutine descriptor_key

    pure logical function integer_key_less(a, b) result(less)
        integer, intent(in) :: a(:), b(:)
        integer :: i
        less = .false.
        do i = 1, size(a)
            if (a(i) < b(i)) then; less = .true.; return
            elseif (a(i) > b(i)) then; return
            end if
        end do
    end function integer_key_less

    pure logical function real_key_less(a, b) result(less)
        real(real64), intent(in) :: a(:), b(:)
        integer :: i
        less = .false.
        do i = 1, size(a)
            if (a(i) < b(i)) then; less = .true.; return
            elseif (a(i) > b(i)) then; return
            end if
        end do
    end function real_key_less

    subroutine write_n2p2_settings(directory, networks, order)
        character(len=*), intent(in) :: directory
        type(atomic_network), intent(in) :: networks(:)
        integer, intent(in) :: order(:)
        integer :: unit, i, j, kind
        character(len=PATH_LENGTH) :: path
        character(len=1) :: activation

        path = join_path(directory, "input.nn")
        open(newunit=unit, file=trim(path), status="replace", action="write")
        write(unit, "(A)") "# n2p2 2G-HDNNP model converted by the Fortran converter"
        write(unit, "(A,1X,I0)") "number_of_elements", size(order)
        write(unit, "(A)", advance="no") "elements"
        do i = 1, size(order); write(unit, "(1X,A)", advance="no") trim(networks(order(i))%atomtype); end do
        write(unit, *)
        do i = 1, size(order)
            j = species_position(networks(1), trim(networks(order(i))%atomtype))
            write(unit, "(A,1X,A,1X,ES25.17E3)") "atom_energy", trim(networks(order(i))%atomtype), &
                networks(1)%atomic_references(j)
        end do
        if (networks(1)%energy_scale /= 1.0_real64 .or. networks(1)%energy_shift /= 0.0_real64) then
            write(unit, "(A,1X,ES25.17E3)") "mean_energy", networks(1)%energy_shift
            write(unit, "(A,1X,ES25.17E3)") "conv_energy", networks(1)%energy_scale
            write(unit, "(A)") "conv_length 1.0"
        end if
        write(unit, "(A,1X,I0,1X,ES25.17E3)") "cutoff_type", &
            networks(1)%descriptor_cutoff_type, networks(1)%descriptor_cutoff_alpha
        write(unit, "(A)") "scale_symmetry_functions_sigma"
        write(unit, "(A)") "scale_min_short 0.0"
        write(unit, "(A)") "scale_max_short 1.0"
        write(unit, "(A,1X,I0)") "global_hidden_layers_short", networks(1)%nlayers - 2
        write(unit, "(A)", advance="no") "global_nodes_short"
        do i = 2, networks(1)%nlayers - 1
            write(unit, "(1X,I0)", advance="no") networks(1)%nodes(i)
        end do
        write(unit, *)
        write(unit, "(A)", advance="no") "global_activation_short"
        do i = 1, networks(1)%nlayers - 1
            activation = activation_name(networks(1)%activation(i))
            write(unit, "(1X,A)", advance="no") activation
        end do
        write(unit, *)
        do i = 1, size(order)
            do j = 1, size(networks(order(i))%descriptor_kinds)
                kind = networks(order(i))%descriptor_kinds(j)
                if (kind == 2) then
                    write(unit, "(A,1X,A,1X,I0,1X,A,3(1X,ES25.17E3))") "symfunction_short", &
                        trim(networks(order(i))%atomtype), 2, &
                        trim(environment_name(networks(order(i)), 1, j)), &
                        networks(order(i))%descriptor_parameters(3, j), &
                        networks(order(i))%descriptor_parameters(2, j), &
                        networks(order(i))%descriptor_parameters(1, j)
                else
                    write(unit, "(A,1X,A,1X,I0,2(1X,A),5(1X,ES25.17E3))") "symfunction_short", &
                        trim(networks(order(i))%atomtype), merge(3, 9, kind == 4), &
                        trim(environment_name(networks(order(i)), 1, j)), &
                        trim(environment_name(networks(order(i)), 2, j)), &
                        networks(order(i))%descriptor_parameters(4, j), &
                        networks(order(i))%descriptor_parameters(2, j), &
                        networks(order(i))%descriptor_parameters(3, j), &
                        networks(order(i))%descriptor_parameters(1, j), 0.0_real64
                end if
            end do
        end do
        close(unit)
        write(*, "(A)") "WROTE "//trim(path)
    end subroutine write_n2p2_settings

    subroutine write_n2p2_scaling(directory, networks, order)
        character(len=*), intent(in) :: directory
        type(atomic_network), intent(in) :: networks(:)
        integer, intent(in) :: order(:)
        integer :: unit, i, j
        real(real64) :: sigma, minimum, maximum
        character(len=PATH_LENGTH) :: path

        path = join_path(directory, "scaling.data")
        open(newunit=unit, file=trim(path), status="replace", action="write")
        write(unit, "(A)") "# element sf minimum maximum mean sigma"
        do i = 1, size(order)
            do j = 1, size(networks(order(i))%descriptor_shift)
                if (networks(order(i))%descriptor_scale(j) == 0.0_real64) &
                    error stop "zero descriptor scale cannot be represented by n2p2"
                sigma = 1.0_real64/networks(order(i))%descriptor_scale(j)
                minimum = networks(order(i))%descriptor_shift(j) - 1.0e6_real64*abs(sigma)
                maximum = networks(order(i))%descriptor_shift(j) + 1.0e6_real64*abs(sigma)
                write(unit, "(I0,1X,I0,4(1X,ES25.17E3))") i, j, minimum, maximum, &
                    networks(order(i))%descriptor_shift(j), sigma
            end do
        end do
        close(unit)
        write(*, "(A)") "WROTE "//trim(path)
    end subroutine write_n2p2_scaling

    subroutine write_n2p2_weights(directory, network)
        character(len=*), intent(in) :: directory
        type(atomic_network), intent(in) :: network
        integer :: unit, layer, source, target, index, connection
        character(len=PATH_LENGTH) :: path, filename

        write(filename, '("weights.",I3.3,".data")') atomic_number(trim(network%atomtype))
        path = join_path(directory, trim(filename))
        open(newunit=unit, file=trim(path), status="replace", action="write")
        write(unit, "(A)") "# Neural network connections converted from AccelNet."
        index = 1; connection = 1
        do layer = 1, network%nlayers - 1
            do source = 1, network%nodes(layer)
                do target = 1, network%nodes(layer + 1)
                    write(unit, "(ES25.17E3,1X,A,5(1X,I0))") network%weights(index), "a", &
                        connection, layer - 1, source, layer, target
                    index = index + 1; connection = connection + 1
                end do
            end do
            do target = 1, network%nodes(layer + 1)
                write(unit, "(ES25.17E3,1X,A,3(1X,I0))") network%weights(index), "b", &
                    connection, layer, target
                index = index + 1; connection = connection + 1
            end do
        end do
        close(unit)
        write(*, "(A)") "WROTE "//trim(path)
    end subroutine write_n2p2_weights

    integer function species_position(network, symbol) result(position)
        type(atomic_network), intent(in) :: network
        character(len=*), intent(in) :: symbol
        integer :: i
        position = 0
        do i = 1, size(network%species_names)
            if (trim(network%species_names(i)) == trim(symbol)) then; position = i; return; end if
        end do
        error stop "network atom type is absent from embedded species"
    end function species_position

    function environment_name(network, row, descriptor) result(name)
        type(atomic_network), intent(in) :: network
        integer, intent(in) :: row, descriptor
        character(len=16) :: name
        integer :: index
        index = network%descriptor_environments(row, descriptor)
        if (index < 1 .or. index > size(network%environment_names)) &
            error stop "descriptor environment index is invalid"
        name = network%environment_names(index)
    end function environment_name

    character(len=1) function activation_name(code) result(name)
        integer, intent(in) :: code
        select case(code)
        case(0); name = "l"
        case(1); name = "t"
        case(2); name = "s"
        case(3); name = "p"
        case(5); name = "r"
        case(6); name = "g"
        case(7); name = "c"
        case(8); name = "S"
        case(9); name = "e"
        case(10); name = "h"
        case default; error stop "activation has no n2p2 equivalent"
        end select
    end function activation_name

    subroutine make_directory(path)
        character(len=*), intent(in) :: path
        integer :: status
        if (index(path, "'") /= 0) error stop "apostrophes in output paths are not supported"
        call execute_command_line("cmake -E make_directory '"//trim(path)//"'", exitstat=status)
        if (status /= 0) error stop "cannot create output directory"
    end subroutine make_directory

    function join_path(directory, name) result(path)
        character(len=*), intent(in) :: directory, name
        character(len=PATH_LENGTH) :: path
        integer :: n
        n = len_trim(directory)
        if (n > 0 .and. directory(n:n) == "/") then
            path = trim(directory)//trim(name)
        else
            path = trim(directory)//"/"//trim(name)
        end if
    end function join_path

    pure function lowercase(input) result(output)
        character(len=*), intent(in) :: input
        character(len=len(input)) :: output
        integer :: i, code
        output = input
        do i = 1, len(input)
            code = iachar(input(i:i))
            if (code >= iachar("A") .and. code <= iachar("Z")) output(i:i) = achar(code + 32)
        end do
    end function lowercase
end module model_conversion
