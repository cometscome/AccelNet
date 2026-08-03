module aenet_network
    use iso_fortran_env, only: real64
    implicit none
    private

    integer, parameter, public :: ACTIVATION_LINEAR = 0
    integer, parameter, public :: ACTIVATION_TANH = 1
    integer, parameter, public :: ACTIVATION_LOGISTIC = 2
    integer, parameter, public :: ACTIVATION_AENET_MTANH = 3
    integer, parameter, public :: ACTIVATION_AENET_TWIST = 4
    integer, parameter, public :: ACTIVATION_RELU = 5
    integer, parameter, public :: ACTIVATION_GAUSSIAN = 6
    integer, parameter, public :: ACTIVATION_COSINE = 7
    integer, parameter, public :: ACTIVATION_REVERSE_LOGISTIC = 8
    integer, parameter, public :: ACTIVATION_EXPONENTIAL = 9
    integer, parameter, public :: ACTIVATION_HARMONIC = 10
    ! Keep native aenet codes 3/4 available for mtanh/twist.  Softplus is an
    ! AccelNet extension used by imported n2p2 models and therefore gets a
    ! non-conflicting serialized code.
    integer, parameter, public :: ACTIVATION_SOFTPLUS = 11

    type, public :: atomic_network
        character(len=16) :: atomtype = ""
        character(len=1024) :: description = ""
        character(len=100) :: descriptor_name = ""
        real(real64) :: minimum_radius = 1.0_real64
        real(real64) :: maximum_radius = 0.0_real64
        character(len=16), allocatable :: environment_names(:)
        integer, allocatable :: descriptor_kinds(:), descriptor_environments(:, :)
        real(real64), allocatable :: descriptor_parameters(:, :)
        integer :: descriptor_cutoff_type = 1
        real(real64) :: descriptor_cutoff_alpha = 0.0_real64
        integer :: nlayers = 0
        integer :: maxnodes = 0
        integer, allocatable :: nodes(:), activation(:), weight_offsets(:)
        real(real64), allocatable :: weights(:)
        real(real64), allocatable :: descriptor_shift(:), descriptor_scale(:)
        character(len=16), allocatable :: species_names(:)
        real(real64), allocatable :: atomic_references(:)
        real(real64) :: energy_scale = 1.0_real64
        real(real64) :: energy_shift = 0.0_real64
    contains
        procedure :: evaluate => evaluate_network
        procedure :: input_gradient => network_input_gradient
    end type atomic_network

    public :: read_aenet_network, read_aenet_network_ascii, read_aenet_network_binary
    public :: write_aenet_network_ascii

contains

    subroutine write_aenet_network_ascii(filename, network)
        character(len=*), intent(in) :: filename
        type(atomic_network), intent(in) :: network
        integer :: unit, ios, nsf, nparam, nvalues, i, layer
        integer, allocatable :: value_offsets(:)
        real(real64), allocatable :: minima(:), maxima(:), moments(:)

        nsf = size(network%descriptor_kinds)
        nparam = size(network%descriptor_parameters, 1)
        allocate(value_offsets(network%nlayers), minima(nsf), maxima(nsf), moments(nsf))
        value_offsets(1) = 0
        do layer = 1, network%nlayers - 1
            value_offsets(layer + 1) = value_offsets(layer) + network%nodes(layer) + 1
        end do
        nvalues = value_offsets(network%nlayers) + network%nodes(network%nlayers)
        do i = 1, nsf
            if (network%descriptor_scale(i) == 0.0_real64) &
                error stop "cannot write a network with zero descriptor scale"
            minima(i) = network%descriptor_shift(i) - abs(1.0_real64/network%descriptor_scale(i))
            maxima(i) = network%descriptor_shift(i) + abs(1.0_real64/network%descriptor_scale(i))
            moments(i) = network%descriptor_shift(i)**2 + &
                1.0_real64/network%descriptor_scale(i)**2
        end do

        open(newunit=unit, file=trim(filename), status="replace", action="write", iostat=ios)
        if (ios /= 0) error stop "cannot create aenet ASCII network"
        write(unit, *) network%nlayers
        write(unit, *) network%maxnodes
        write(unit, *) size(network%weights)
        write(unit, *) nvalues
        write(unit, *) network%nodes
        write(unit, *) network%activation
        write(unit, *) network%weight_offsets
        write(unit, *) value_offsets
        write(unit, *) network%weights
        write(unit, "(A)") trim(network%description)
        write(unit, "(A)") trim(network%atomtype)
        write(unit, *) size(network%environment_names)
        write(unit, *) network%environment_names
        write(unit, *) network%minimum_radius
        write(unit, *) network%maximum_radius
        write(unit, "(A)") trim(network%descriptor_name)
        write(unit, *) nsf
        write(unit, *) nparam
        write(unit, *) network%descriptor_kinds
        write(unit, *) network%descriptor_parameters
        write(unit, *) network%descriptor_environments
        write(unit, *) 0
        write(unit, *) minima
        write(unit, *) maxima
        write(unit, *) network%descriptor_shift
        write(unit, *) moments
        write(unit, "(A)") "converted model"
        write(unit, *) .true.
        write(unit, *) network%energy_scale
        write(unit, *) network%energy_shift
        write(unit, *) size(network%species_names)
        write(unit, *) network%species_names
        write(unit, *) network%atomic_references
        write(unit, *) 0
        write(unit, *) 0
        write(unit, *) 0.0_real64, 0.0_real64, 0.0_real64
        close(unit)
    end subroutine write_aenet_network_ascii

    subroutine read_aenet_network_ascii(filename, network)
        character(len=*), intent(in) :: filename
        type(atomic_network), intent(out) :: network
        integer :: unit, ios, nweights, nvalues, nenv, nsf, nparam, neval
        integer :: ntypes, natoms, nstructures
        integer, allocatable :: integers(:)
        real(real64), allocatable :: reals(:), averages(:), moments(:)
        real(real64) :: rcmin, rcmax, emin, emax, eavg, variance
        character(len=1024) :: line, training_file
        logical :: normalized

        open(newunit=unit, file=trim(filename), status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "cannot open aenet ASCII network"
        read(unit, *) network%nlayers
        read(unit, *) network%maxnodes
        read(unit, *) nweights
        read(unit, *) nvalues
        allocate(network%nodes(network%nlayers), network%activation(network%nlayers - 1))
        allocate(network%weight_offsets(network%nlayers), network%weights(nweights))
        read(unit, *) network%nodes
        read(unit, *) network%activation
        read(unit, *) network%weight_offsets
        read(unit, *) ! value offsets
        read(unit, *) network%weights

        read(unit, "(A)") network%description
        read(unit, "(A)") network%atomtype
        network%atomtype = adjustl(network%atomtype)
        read(unit, *) nenv
        allocate(network%environment_names(nenv)); read(unit, *) network%environment_names
        read(unit, *) network%minimum_radius
        read(unit, *) network%maximum_radius
        read(unit, "(A)") network%descriptor_name
        read(unit, *) nsf
        read(unit, *) nparam
        if (nsf /= network%nodes(1)) error stop "network and descriptor dimensions differ"
        allocate(network%descriptor_kinds(nsf), network%descriptor_parameters(nparam, nsf), &
                 network%descriptor_environments(2, nsf))
        allocate(integers(max(nsf, 2*nsf)), reals(max(1, nsf, nparam*nsf)))
        read(unit, *) network%descriptor_kinds
        read(unit, *) network%descriptor_parameters
        call recover_cutoff_metadata(network)
        read(unit, *) network%descriptor_environments
        read(unit, *) neval
        allocate(averages(nsf), moments(nsf), network%descriptor_shift(nsf), &
                 network%descriptor_scale(nsf))
        read(unit, *) reals(1:nsf) ! minima
        read(unit, *) reals(1:nsf) ! maxima
        read(unit, *) averages
        read(unit, *) moments
        network%descriptor_shift = averages
        block
            integer :: i
            do i = 1, nsf
                variance = max(moments(i) - averages(i)*averages(i), 0.0_real64)
                if (variance > 0.0_real64) then
                    network%descriptor_scale(i) = 1.0_real64/sqrt(variance)
                else
                    network%descriptor_scale(i) = 1.0_real64
                end if
            end do
        end block
        read(unit, "(A)") training_file
        read(unit, *) normalized
        read(unit, *) network%energy_scale
        read(unit, *) network%energy_shift
        read(unit, *) ntypes
        allocate(network%species_names(ntypes), network%atomic_references(ntypes))
        read(unit, *) network%species_names
        read(unit, *) network%atomic_references
        read(unit, *) natoms
        read(unit, *) nstructures
        read(unit, *) emin, emax, eavg
        close(unit)
        call upgrade_legacy_n2p2_activations(network)
    end subroutine read_aenet_network_ascii

    subroutine read_aenet_network(filename, network)
        character(len=*), intent(in) :: filename
        type(atomic_network), intent(out) :: network
        integer :: name_length
        name_length = len_trim(filename)
        if (name_length >= 6 .and. filename(name_length - 5:name_length) == ".ascii") then
            call read_aenet_network_ascii(filename, network)
        else
            call read_aenet_network_binary(filename, network)
        end if
    end subroutine read_aenet_network

    subroutine read_aenet_network_binary(filename, network)
        character(len=*), intent(in) :: filename
        type(atomic_network), intent(out) :: network
        integer :: unit, ios, nweights, nvalues, nenv, nsf, nparam, neval
        integer :: ntypes, natoms, nstructures, i
        integer, allocatable :: value_offsets(:)
        real(real64), allocatable :: minima(:), maxima(:), moments(:)
        character(len=1024) :: training_file
        character(len=2) :: atomtype
        character(len=2), allocatable :: environment_names(:), species_names(:)
        logical :: normalized
        real(real64) :: emin, emax, eavg, variance
        open(newunit=unit, file=trim(filename), status="old", action="read", &
             form="unformatted", iostat=ios)
        if (ios /= 0) then
            write(*, "(A,1X,A)") "cannot open aenet binary network:", trim(filename)
            error stop "cannot open aenet binary network"
        end if
        read(unit) network%nlayers; read(unit) network%maxnodes
        read(unit) nweights; read(unit) nvalues
        allocate(network%nodes(network%nlayers), network%activation(network%nlayers - 1), &
                 network%weight_offsets(network%nlayers), value_offsets(network%nlayers), &
                 network%weights(nweights))
        read(unit) network%nodes; read(unit) network%activation
        read(unit) network%weight_offsets; read(unit) value_offsets; read(unit) network%weights
        read(unit) network%description; read(unit) atomtype; read(unit) nenv
        network%atomtype = atomtype
        allocate(environment_names(nenv), network%environment_names(nenv)); read(unit) environment_names
        network%environment_names = environment_names
        read(unit) network%minimum_radius; read(unit) network%maximum_radius
        read(unit) network%descriptor_name; read(unit) nsf; read(unit) nparam
        if (nsf /= network%nodes(1)) error stop "network and descriptor dimensions differ"
        allocate(network%descriptor_kinds(nsf), network%descriptor_parameters(nparam, nsf), &
                 network%descriptor_environments(2, nsf), minima(nsf), maxima(nsf), &
                 network%descriptor_shift(nsf), network%descriptor_scale(nsf), moments(nsf))
        read(unit) network%descriptor_kinds; read(unit) network%descriptor_parameters
        call recover_cutoff_metadata(network)
        read(unit) network%descriptor_environments; read(unit) neval
        read(unit) minima; read(unit) maxima; read(unit) network%descriptor_shift; read(unit) moments
        read(unit) training_file; read(unit) normalized; read(unit) network%energy_scale
        read(unit) network%energy_shift; read(unit) ntypes
        allocate(species_names(ntypes), network%species_names(ntypes), network%atomic_references(ntypes))
        read(unit) species_names; network%species_names = species_names
        read(unit) network%atomic_references
        read(unit) natoms; read(unit) nstructures; read(unit) emin, emax, eavg
        close(unit)
        call upgrade_legacy_n2p2_activations(network)
        do i = 1, nsf
            variance = max(moments(i) - network%descriptor_shift(i)**2, 0.0_real64)
            if (variance > 0.0_real64) then
                network%descriptor_scale(i) = 1.0_real64/sqrt(variance)
            else
                network%descriptor_scale(i) = 1.0_real64
            end if
        end do
    end subroutine read_aenet_network_binary

    subroutine recover_cutoff_metadata(network)
        type(atomic_network), intent(inout) :: network
        integer :: nparam
        if (trim(lowercase(network%descriptor_name)) == "lj") then
            network%descriptor_cutoff_type = 0
        else
            network%descriptor_cutoff_type = 1
        end if
        network%descriptor_cutoff_alpha = 0.0_real64
        nparam = size(network%descriptor_parameters, 1)
        if (nparam < 6) return
        network%descriptor_cutoff_type = nint(network%descriptor_parameters(5, 1))
        network%descriptor_cutoff_alpha = network%descriptor_parameters(6, 1)
        if (any(nint(network%descriptor_parameters(5, :)) /= network%descriptor_cutoff_type) .or. &
            any(network%descriptor_parameters(6, :) /= network%descriptor_cutoff_alpha)) &
            error stop "embedded descriptors use inconsistent cutoff metadata"
    end subroutine recover_cutoff_metadata

    subroutine upgrade_legacy_n2p2_activations(network)
        type(atomic_network), intent(inout) :: network
        ! AccelNet 0.1.0 wrote n2p2 softplus with integer code 3, which is
        ! native aenet's mtanh code.  Preserve those generated files while all
        ! newly written files use ACTIVATION_SOFTPLUS=11.
        if (index(trim(network%description), "Imported n2p2 2G-HDNNP model") == 1) then
            where (network%activation == ACTIVATION_AENET_MTANH)
                network%activation = ACTIVATION_SOFTPLUS
            end where
        end if
    end subroutine upgrade_legacy_n2p2_activations

    subroutine evaluate_network(self, input, output)
        class(atomic_network), intent(in) :: self
        real(real64), intent(in) :: input(:)
        real(real64), intent(out) :: output
        real(real64) :: values(self%maxnodes, self%nlayers)
        integer :: layer, i, j, nin, nout, offset
        values = 0.0_real64
        values(1:self%nodes(1), 1) = input
        do layer = 1, self%nlayers - 1
            nin = self%nodes(layer); nout = self%nodes(layer + 1)
            offset = self%weight_offsets(layer) + 1
            do j = 1, nout
                values(j, layer + 1) = self%weights(offset + nin*nout + j - 1)
                do i = 1, nin
                    values(j, layer + 1) = values(j, layer + 1) + &
                        self%weights(offset + (i - 1)*nout + j - 1)*values(i, layer)
                end do
                values(j, layer + 1) = activate(values(j, layer + 1), self%activation(layer))
            end do
        end do
        output = values(1, self%nlayers)
    end subroutine evaluate_network

    subroutine network_input_gradient(self, input, output, gradient)
        class(atomic_network), intent(in) :: self
        real(real64), intent(in) :: input(:)
        real(real64), intent(out) :: output, gradient(:)
        real(real64) :: values(self%maxnodes, self%nlayers), preactivation(self%maxnodes, self%nlayers)
        real(real64) :: delta(self%maxnodes), previous(self%maxnodes)
        integer :: layer, i, j, nin, nout, offset
        values = 0.0_real64; preactivation = 0.0_real64; values(1:self%nodes(1), 1) = input
        do layer = 1, self%nlayers - 1
            nin = self%nodes(layer); nout = self%nodes(layer + 1); offset = self%weight_offsets(layer) + 1
            do j = 1, nout
                values(j, layer + 1) = self%weights(offset + nin*nout + j - 1)
                do i = 1, nin
                    values(j, layer + 1) = values(j, layer + 1) + self%weights(offset + (i - 1)*nout + j - 1)*values(i, layer)
                end do
                preactivation(j, layer + 1) = values(j, layer + 1)
                values(j, layer + 1) = activate(preactivation(j, layer + 1), self%activation(layer))
            end do
        end do
        output = values(1, self%nlayers); delta = 0.0_real64; delta(1) = 1.0_real64
        do layer = self%nlayers - 1, 1, -1
            nin = self%nodes(layer); nout = self%nodes(layer + 1); offset = self%weight_offsets(layer) + 1
            do j = 1, nout
                delta(j) = delta(j)*activation_derivative(preactivation(j, layer + 1), &
                    values(j, layer + 1), self%activation(layer))
            end do
            previous = 0.0_real64
            do i = 1, nin
                do j = 1, nout
                    previous(i) = previous(i) + self%weights(offset + (i - 1)*nout + j - 1)*delta(j)
                end do
            end do
            delta = previous
        end do
        gradient = delta(1:self%nodes(1))
    end subroutine network_input_gradient

    pure real(real64) function activate(x, code) result(y)
        real(real64), intent(in) :: x
        integer, intent(in) :: code
        real(real64), parameter :: a = 1.7159_real64
        real(real64), parameter :: b = 0.666666666666667_real64
        real(real64), parameter :: c = 0.1_real64
        select case(code)
        case(ACTIVATION_LINEAR); y = x
        case(ACTIVATION_TANH); y = tanh(x)
        case(ACTIVATION_LOGISTIC); y = 1.0_real64/(1.0_real64 + exp(-x))
        case(ACTIVATION_AENET_MTANH)
            y = a*tanh(b*x)
        case(ACTIVATION_AENET_TWIST)
            y = a*tanh(b*x) + c*x
        case(ACTIVATION_SOFTPLUS)
            if (x > 0.0_real64) then
                y = x + log(1.0_real64 + exp(-x))
            else
                y = log(1.0_real64 + exp(x))
            end if
        case(ACTIVATION_RELU); y = max(x, 0.0_real64)
        case(ACTIVATION_GAUSSIAN); y = exp(-0.5_real64*x*x)
        case(ACTIVATION_COSINE); y = cos(x)
        case(ACTIVATION_REVERSE_LOGISTIC); y = 1.0_real64 - 1.0_real64/(1.0_real64 + exp(-x))
        case(ACTIVATION_EXPONENTIAL); y = exp(-x)
        case(ACTIVATION_HARMONIC); y = x*x
        case default; error stop "unsupported aenet activation"
        end select
    end function activate

    pure real(real64) function activation_derivative(x, y, code) result(value)
        real(real64), intent(in) :: x, y
        integer, intent(in) :: code
        real(real64), parameter :: a = 1.7159_real64
        real(real64), parameter :: b = 0.666666666666667_real64
        real(real64), parameter :: c = 0.1_real64
        real(real64) :: tanhbx
        select case(code)
        case(ACTIVATION_LINEAR); value = 1.0_real64
        case(ACTIVATION_TANH); value = 1.0_real64 - y*y
        case(ACTIVATION_LOGISTIC); value = y*(1.0_real64 - y)
        case(ACTIVATION_AENET_MTANH)
            tanhbx = tanh(b*x)
            value = a*b*(1.0_real64 - tanhbx*tanhbx)
        case(ACTIVATION_AENET_TWIST)
            tanhbx = tanh(b*x)
            value = a*b*(1.0_real64 - tanhbx*tanhbx) + c
        case(ACTIVATION_SOFTPLUS); value = 1.0_real64 - exp(-y)
        case(ACTIVATION_RELU); value = merge(1.0_real64, 0.0_real64, y > 0.0_real64)
        case(ACTIVATION_GAUSSIAN); value = -x*y
        case(ACTIVATION_COSINE); value = -sin(x)
        case(ACTIVATION_REVERSE_LOGISTIC); value = -y*(1.0_real64 - y)
        case(ACTIVATION_EXPONENTIAL); value = -y
        case(ACTIVATION_HARMONIC); value = 2.0_real64*x
        case default; error stop "unsupported aenet activation"
        end select
    end function activation_derivative

    pure function lowercase(input) result(output)
        character(len=*), intent(in) :: input
        character(len=len(input)) :: output
        integer :: i, code
        output = input
        do i = 1, len(input)
            code = iachar(input(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) output(i:i) = achar(code + 32)
        end do
    end function lowercase
end module aenet_network
