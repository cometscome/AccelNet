program test_binary_network
    use iso_fortran_env, only: real64
    use aenet_network, only: atomic_network, read_aenet_network
    use accelnet_predictor, only: predictor_model, load_predictor_from_networks
    implicit none
    type(atomic_network) :: network
    type(predictor_model) :: ascii_model, binary_model
    character(len=1024) :: ascii_files(2), binary_files(2), structure
    real(real64) :: ascii_energy, binary_energy
    integer :: i

    do i = 1, 2
        call get_command_argument(i, ascii_files(i))
        call get_command_argument(i + 2, binary_files(i))
        call convert_like_accelnet(trim(ascii_files(i)), trim(binary_files(i)))
    end do
    call get_command_argument(5, structure)

    call read_aenet_network(trim(binary_files(1)), network)
    if (trim(network%atomtype) /= "Ti") error stop "binary atom type differs"
    if (network%nodes(1) /= 56) error stop "binary network dimensions differ"

    call load_predictor_from_networks(ascii_files, ascii_model)
    call load_predictor_from_networks(binary_files, binary_model)
    call ascii_model%predict_energy(trim(structure), ascii_energy)
    call binary_model%predict_energy(trim(structure), binary_energy)
    if (abs(binary_energy - ascii_energy) > 1.0e-12_real64) then
        write(*, *) "ASCII/binary energy difference:", binary_energy - ascii_energy
        error stop "binary aenet network changes prediction"
    end if
    call ascii_model%reload(binary_files)
    call ascii_model%predict_energy(trim(structure), binary_energy)
    if (abs(binary_energy - ascii_energy) > 1.0e-12_real64) &
        error stop "reloaded binary network changes prediction"

contains

    ! This writes the native sequential-unformatted record layout used by
    ! AccelNet's nnASCII2bin utility.  Keeping the fixture conversion in the
    ! test makes the production predictor independent of the original tree.
    subroutine convert_like_accelnet(ascii_file, binary_file)
        character(len=*), intent(in) :: ascii_file, binary_file
        integer :: input, output, nlayers, maxnodes, nweights, nvalues
        integer :: nenv, nsf, nparam, neval, ntypes, natoms, nstructures
        integer, allocatable :: nodes(:), activation(:), weight_offsets(:), value_offsets(:)
        integer, allocatable :: kinds(:), environments(:, :)
        real(real64), allocatable :: weights(:), parameters(:, :), minima(:), maxima(:), averages(:), moments(:)
        real(real64), allocatable :: atomic_references(:)
        real(real64) :: minimum_radius, maximum_radius, scale, shift, emin, emax, eavg
        character(len=1024) :: description, training_file
        character(len=100) :: descriptor_name
        character(len=2) :: atomtype
        character(len=2), allocatable :: environment_names(:), species_names(:)
        logical :: normalized

        open(newunit=input, file=ascii_file, status="old", action="read")
        open(newunit=output, file=binary_file, status="replace", action="write", form="unformatted")

        read(input, *) nlayers; read(input, *) maxnodes
        read(input, *) nweights; read(input, *) nvalues
        allocate(nodes(nlayers), activation(nlayers - 1), weight_offsets(nlayers), &
                 value_offsets(nlayers), weights(nweights))
        read(input, *) nodes; read(input, *) activation; read(input, *) weight_offsets
        read(input, *) value_offsets; read(input, *) weights
        write(output) nlayers; write(output) maxnodes; write(output) nweights; write(output) nvalues
        write(output) nodes; write(output) activation; write(output) weight_offsets
        write(output) value_offsets; write(output) weights

        read(input, *) description; read(input, *) atomtype; read(input, *) nenv
        allocate(environment_names(nenv)); read(input, *) environment_names
        read(input, *) minimum_radius; read(input, *) maximum_radius
        read(input, *) descriptor_name; read(input, *) nsf; read(input, *) nparam
        allocate(kinds(nsf), parameters(nparam, nsf), environments(2, nsf), &
                 minima(nsf), maxima(nsf), averages(nsf), moments(nsf))
        read(input, *) kinds; read(input, *) parameters; read(input, *) environments
        read(input, *) neval; read(input, *) minima; read(input, *) maxima
        read(input, *) averages; read(input, *) moments
        write(output) description; write(output) atomtype; write(output) nenv
        write(output) environment_names; write(output) minimum_radius; write(output) maximum_radius
        write(output) descriptor_name; write(output) nsf; write(output) nparam
        write(output) kinds; write(output) parameters; write(output) environments
        write(output) neval; write(output) minima; write(output) maxima
        write(output) averages; write(output) moments

        read(input, *) training_file; read(input, *) normalized
        read(input, *) scale; read(input, *) shift; read(input, *) ntypes
        allocate(species_names(ntypes), atomic_references(ntypes))
        read(input, *) species_names; read(input, *) atomic_references
        read(input, *) natoms; read(input, *) nstructures; read(input, *) emin, emax, eavg
        write(output) training_file; write(output) normalized
        write(output) scale; write(output) shift; write(output) ntypes
        write(output) species_names; write(output) atomic_references
        write(output) natoms; write(output) nstructures; write(output) emin, emax, eavg
        close(input); close(output)
    end subroutine convert_like_accelnet
end program test_binary_network
