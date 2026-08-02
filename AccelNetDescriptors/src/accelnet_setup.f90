module accelnet_setup
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: descriptor_config, initialize_config, validate_cutoff_parameters
    use accelnet_lj, only: lj_config, initialize_lj_config
    use accelnet_behler
    use accelnet_descriptor_models
    implicit none
    private

    integer, parameter :: LINE_LENGTH = 2048
    integer, parameter :: NAME_LENGTH = 16

    type :: raw_behler_function
        integer :: kind = 0, species1 = 0, species2 = 0
        real(real64) :: rc = 0.0_real64
        real(real64) :: parameter2 = 0.0_real64
        real(real64) :: parameter3 = 0.0_real64
        real(real64) :: parameter4 = 0.0_real64
    end type raw_behler_function

    type, public :: descriptor_setup
        character(len=NAME_LENGTH) :: central_species = ""
        character(len=1024) :: description = ""
        real(real64) :: minimum_distance = 1.0_real64
        character(len=NAME_LENGTH), allocatable :: environment_species(:)
        integer, allocatable :: global_to_local(:)
        integer :: central_global_species = 0
        type(descriptor_model) :: model
    contains
        procedure :: map_species => map_setup_species
    end type descriptor_setup

    public :: read_accelnet_setup
    public :: read_accelnet_setup_set

contains

    subroutine read_accelnet_setup_set(filenames, setups, global_species)
        character(len=*), intent(in) :: filenames(:)
        type(descriptor_setup), allocatable, intent(out) :: setups(:)
        character(len=NAME_LENGTH), allocatable, intent(out) :: global_species(:)
        integer :: i, j

        if (size(filenames) < 1) error stop "at least one setup file is required"
        allocate(setups(size(filenames)), global_species(size(filenames)))
        do i = 1, size(filenames)
            call scan_central_species(filenames(i), global_species(i))
            do j = 1, i - 1
                if (trim(global_species(i)) == trim(global_species(j))) &
                    error stop "duplicate ATOM species in setup set"
            end do
        end do
        do i = 1, size(filenames)
            call read_accelnet_setup(trim(filenames(i)), global_species, setups(i))
            if (setups(i)%central_global_species /= i) &
                error stop "setup set central-species ordering is inconsistent"
        end do
    end subroutine read_accelnet_setup_set

    subroutine scan_central_species(filename, species)
        character(len=*), intent(in) :: filename
        character(len=*), intent(out) :: species
        character(len=LINE_LENGTH) :: line, keyword
        integer :: unit, ios
        logical :: eof
        species = ""
        open(newunit=unit, file=trim(filename), status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "cannot open AccelNet setup file"
        do
            call read_valid_line(unit, line, eof)
            if (eof) exit
            read(line, *, iostat=ios) keyword
            if (ios == 0 .and. lowercase(trim(keyword)) == "atom") then
                read(line(len_trim(keyword) + 1:), *) species
                exit
            end if
        end do
        close(unit)
        if (len_trim(species) == 0) error stop "setup file has no ATOM keyword"
    end subroutine scan_central_species

    subroutine read_accelnet_setup(filename, global_species, setup)
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: global_species(:)
        type(descriptor_setup), intent(out) :: setup
        character(len=LINE_LENGTH) :: line, keyword
        integer :: unit, ios, i, number_of_kinds, number_of_functions, cutoff_type
        real(real64) :: cutoff_alpha
        logical :: eof

        open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "cannot open AccelNet setup file"
        do
            call read_valid_line(unit, line, eof)
            if (eof) exit
            read(line, *, iostat=ios) keyword
            if (ios /= 0) cycle
            select case(lowercase(trim(keyword)))
            case ("descr")
                call read_description(unit, setup%description)
            case ("atom")
                read(line(len_trim(keyword) + 1:), *) setup%central_species
                setup%central_global_species = species_position(setup%central_species, global_species)
                if (setup%central_global_species == 0) &
                    error stop "setup central species is not global"
            case ("env")
                read(line(len_trim(keyword) + 1:), *) number_of_kinds
                if (number_of_kinds < 1) error stop "ENV count must be positive"
                allocate(setup%environment_species(number_of_kinds))
                do i = 1, number_of_kinds
                    call read_valid_line(unit, line, eof)
                    if (eof) error stop "unexpected EOF in ENV block"
                    read(line, *) setup%environment_species(i)
                end do
            case ("rmin")
                read(line(len_trim(keyword) + 1:), *) setup%minimum_distance
            case ("basis")
                call require_environment(setup)
                call parse_basis_header(unit, line, setup)
            case ("symmfunc", "functions")
                call require_environment(setup)
                call parse_behler_header(line, cutoff_type, cutoff_alpha)
                call read_valid_line(unit, line, eof)
                if (eof) error stop "missing Behler function count"
                read(line, *) number_of_functions
                call parse_behler_block(unit, line, number_of_functions, cutoff_type, cutoff_alpha, setup)
            case default
                error stop "unknown keyword in AccelNet setup file"
            end select
        end do
        close(unit)
        call finalize_species_mapping(setup, global_species)
        if (setup%model%num_descriptors() < 1) error stop "setup contains no supported descriptors"
    end subroutine read_accelnet_setup

    subroutine parse_behler_header(line, cutoff_type, cutoff_alpha)
        character(len=*), intent(in) :: line
        integer, intent(out) :: cutoff_type
        real(real64), intent(out) :: cutoff_alpha
        character(len=LINE_LENGTH) :: function_type
        call get_string_value(line, "type", function_type)
        if (lowercase(trim(function_type)) /= "behler2011") &
            error stop "unsupported SYMMFUNC type in AccelNet setup"
        call get_integer_value(line, "cutoff_type", cutoff_type, 1)
        call get_real_value(line, "cutoff_alpha", cutoff_alpha, 0.0_real64)
        call validate_cutoff_parameters(cutoff_type, cutoff_alpha)
    end subroutine parse_behler_header

    subroutine parse_basis_header(unit, header, setup)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: header
        type(descriptor_setup), intent(inout) :: setup
        character(len=LINE_LENGTH) :: basis_type, line
        integer :: number_of_kinds, i, cutoff_type
        real(real64) :: cutoff_alpha
        logical :: eof

        call get_string_value(header, "type", basis_type)
        call get_integer_value(header, "cutoff_type", cutoff_type, -1)
        call get_real_value(header, "cutoff_alpha", cutoff_alpha, 0.0_real64)
        select case(lowercase(trim(basis_type)))
        case ("chebyshev")
            if (cutoff_type < 0) cutoff_type = 1
            call validate_cutoff_parameters(cutoff_type, cutoff_alpha)
            call read_valid_line(unit, line, eof)
            if (eof) error stop "missing Chebyshev parameter line"
            call add_chebyshev_from_line(line, setup, cutoff_type, cutoff_alpha)
        case ("lj")
            if (cutoff_type < 0) cutoff_type = 0
            call validate_cutoff_parameters(cutoff_type, cutoff_alpha)
            call read_valid_line(unit, line, eof)
            if (eof) error stop "missing LJ parameter line"
            call add_lj_from_line(line, setup, cutoff_type, cutoff_alpha)
        case ("multi")
            call read_valid_line(unit, line, eof)
            if (eof) error stop "missing multi basis count"
            read(line, *) number_of_kinds
            do i = 1, number_of_kinds
                call read_valid_line(unit, basis_type, eof)
                if (eof) error stop "missing multi basis type"
                call read_valid_line(unit, line, eof)
                if (eof) error stop "missing multi basis parameters"
                select case(lowercase(trim(adjustl(basis_type))))
                case ("chebyshev")
                    if (cutoff_type < 0) then
                        call add_chebyshev_from_line(line, setup, 1, cutoff_alpha)
                    else
                        call add_chebyshev_from_line(line, setup, cutoff_type, cutoff_alpha)
                    end if
                case ("lj")
                    if (cutoff_type < 0) then
                        call add_lj_from_line(line, setup, 0, cutoff_alpha)
                    else
                        call add_lj_from_line(line, setup, cutoff_type, cutoff_alpha)
                    end if
                case default
                    error stop "unsupported basis inside AccelNet multi setup"
                end select
            end do
        case default
            error stop "unsupported BASIS type in AccelNet setup"
        end select
    end subroutine parse_basis_header

    subroutine add_chebyshev_from_line(line, setup, cutoff_type, cutoff_alpha)
        character(len=*), intent(in) :: line
        type(descriptor_setup), intent(inout) :: setup
        integer, intent(in) :: cutoff_type
        real(real64), intent(in) :: cutoff_alpha
        type(descriptor_config) :: config
        real(real64) :: radial_rc, angular_rc
        integer :: radial_order, angular_order, version
        call get_real_value(line, "radial_Rc", radial_rc, 8.0_real64)
        call get_integer_value(line, "radial_N", radial_order, 10)
        call get_real_value(line, "angular_Rc", angular_rc, 6.5_real64)
        call get_integer_value(line, "angular_N", angular_order, 4)
        call get_integer_value(line, "version", version, 0)
        call initialize_config(config, size(setup%environment_species), radial_rc, radial_order, &
                               angular_rc, angular_order, version, setup%central_global_species, &
                               cutoff_type, cutoff_alpha)
        call add_chebyshev(setup%model, config)
    end subroutine add_chebyshev_from_line

    subroutine add_lj_from_line(line, setup, cutoff_type, cutoff_alpha)
        character(len=*), intent(in) :: line
        type(descriptor_setup), intent(inout) :: setup
        integer, intent(in) :: cutoff_type
        real(real64), intent(in) :: cutoff_alpha
        type(lj_config) :: config
        real(real64) :: radial_rc
        call get_real_value(line, "radial_Rc", radial_rc, 0.0_real64)
        if (radial_rc <= 0.0_real64) error stop "LJ setup requires radial_Rc"
        call initialize_lj_config(config, size(setup%environment_species), radial_rc, &
            cutoff_type, cutoff_alpha)
        call add_lj(setup%model, config)
    end subroutine add_lj_from_line

    subroutine parse_behler_block(unit, line, number_of_functions, cutoff_type, cutoff_alpha, setup)
        integer, intent(in) :: unit, number_of_functions, cutoff_type
        real(real64), intent(in) :: cutoff_alpha
        character(len=*), intent(inout) :: line
        type(descriptor_setup), intent(inout) :: setup
        type(raw_behler_function), allocatable :: functions(:)
        type(behler_config) :: config
        integer :: i, kind, species1, species2
        logical :: eof

        if (number_of_functions < 1) error stop "Behler function count must be positive"
        allocate(functions(number_of_functions))
        do i = 1, number_of_functions
            call read_valid_line(unit, line, eof)
            if (eof) error stop "unexpected EOF in Behler function block"
            call get_integer_value(line, "G", kind, 0)
            functions(i)%kind = kind
            call get_string_species(line, "type2", setup%environment_species, species1)
            functions(i)%species1 = species1
            call get_real_value(line, "Rc", functions(i)%rc, 0.0_real64)
            select case(kind)
            case (1)
            case (2)
                call get_real_value(line, "Rs", functions(i)%parameter2, 0.0_real64)
                call get_real_value(line, "eta", functions(i)%parameter3, 0.0_real64)
            case (3)
                call get_real_value(line, "kappa", functions(i)%parameter2, 0.0_real64)
            case (4, 5)
                call get_string_species(line, "type3", setup%environment_species, species2)
                functions(i)%species2 = species2
                call get_real_value(line, "lambda", functions(i)%parameter2, 0.0_real64)
                call get_real_value(line, "zeta", functions(i)%parameter3, 0.0_real64)
                call get_real_value(line, "eta", functions(i)%parameter4, 0.0_real64)
            case default
                error stop "unsupported Behler G function"
            end select
            if (functions(i)%rc <= 0.0_real64) error stop "Behler function requires positive Rc"
        end do
        call initialize_behler_config(config, size(setup%environment_species), cutoff_type, cutoff_alpha)
        call add_sorted_behler_functions(config, functions)
        call add_behler(setup%model, config)
    end subroutine parse_behler_block

    subroutine add_sorted_behler_functions(config, functions)
        type(behler_config), intent(inout) :: config
        type(raw_behler_function), intent(in) :: functions(:)
        integer :: species1, species2, kind, i, low_species, high_species
        do species1 = 1, config%num_species
            do kind = 1, 3
                do i = 1, size(functions)
                    if (functions(i)%kind /= kind .or. functions(i)%species1 /= species1) cycle
                    select case(kind)
                    case (1)
                        call add_g1(config, species1, functions(i)%rc)
                    case (2)
                        call add_g2(config, species1, functions(i)%rc, functions(i)%parameter2, &
                                    functions(i)%parameter3)
                    case (3)
                        call add_g3(config, species1, functions(i)%rc, functions(i)%parameter2)
                    end select
                end do
            end do
            do species2 = species1, config%num_species
                do kind = 4, 5
                    do i = 1, size(functions)
                        low_species = min(functions(i)%species1, functions(i)%species2)
                        high_species = max(functions(i)%species1, functions(i)%species2)
                        if (functions(i)%kind /= kind .or. low_species /= species1 .or. &
                            high_species /= species2) cycle
                        if (kind == 4) then
                            call add_g4(config, species1, species2, functions(i)%rc, &
                                functions(i)%parameter2, functions(i)%parameter3, functions(i)%parameter4)
                        else
                            call add_g5(config, species1, species2, functions(i)%rc, &
                                functions(i)%parameter2, functions(i)%parameter3, functions(i)%parameter4)
                        end if
                    end do
                end do
            end do
        end do
    end subroutine add_sorted_behler_functions

    subroutine finalize_species_mapping(setup, global_species)
        type(descriptor_setup), intent(inout) :: setup
        character(len=*), intent(in) :: global_species(:)
        integer :: global, local
        allocate(setup%global_to_local(size(global_species)), source=0)
        do global = 1, size(global_species)
            if (trim(global_species(global)) == trim(setup%central_species)) &
                setup%central_global_species = global
            do local = 1, size(setup%environment_species)
                if (trim(global_species(global)) == trim(setup%environment_species(local))) &
                    setup%global_to_local(global) = local
            end do
        end do
        if (setup%central_global_species == 0) error stop "setup central species is not global"
        if (any(setup%global_to_local == 0)) error stop "global species is missing from setup ENV"
    end subroutine finalize_species_mapping

    integer function species_position(name, species_names) result(position)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: species_names(:)
        integer :: i
        position = 0
        do i = 1, size(species_names)
            if (trim(name) == trim(species_names(i))) then
                position = i
                return
            end if
        end do
    end function species_position

    subroutine map_setup_species(self, global_species, local_species)
        class(descriptor_setup), intent(in) :: self
        integer, intent(in) :: global_species(:)
        integer, intent(out) :: local_species(:)
        integer :: i
        if (size(local_species) /= size(global_species)) error stop "species mapping arrays differ"
        do i = 1, size(global_species)
            if (global_species(i) < 1 .or. global_species(i) > size(self%global_to_local)) &
                error stop "global species index out of range"
            local_species(i) = self%global_to_local(global_species(i))
        end do
    end subroutine map_setup_species

    subroutine require_environment(setup)
        type(descriptor_setup), intent(in) :: setup
        if (.not. allocated(setup%environment_species)) error stop "BASIS appears before ENV"
    end subroutine require_environment

    subroutine read_description(unit, description)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: description
        character(len=LINE_LENGTH) :: line
        integer :: ios
        description = ""
        do
            read(unit, "(A)", iostat=ios) line
            if (ios /= 0) error stop "unterminated DESCR block"
            if (index(lowercase(trim(adjustl(line))), "end descr") == 1) exit
            if (len_trim(description) > 0) description = trim(description)//" "
            if (len_trim(description) + len_trim(line) <= len(description)) &
                description = trim(description)//trim(adjustl(line))
        end do
    end subroutine read_description

    subroutine read_valid_line(unit, line, eof)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: line
        logical, intent(out) :: eof
        integer :: ios
        eof = .false.
        do
            read(unit, "(A)", iostat=ios) line
            if (ios < 0) then
                eof = .true.
                return
            else if (ios > 0) then
                error stop "error reading AccelNet setup"
            end if
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (index("#!%", line(1:1)) > 0) cycle
            return
        end do
    end subroutine read_valid_line

    subroutine get_string_species(line, key, species_names, species_index)
        character(len=*), intent(in) :: line, key
        character(len=*), intent(in) :: species_names(:)
        integer, intent(out) :: species_index
        character(len=LINE_LENGTH) :: value
        integer :: i
        call get_string_value(line, key, value)
        species_index = 0
        do i = 1, size(species_names)
            if (trim(value) == trim(species_names(i))) species_index = i
        end do
        if (species_index == 0) error stop "Behler function species is not in ENV"
    end subroutine get_string_species

    subroutine get_string_value(line, key, value)
        character(len=*), intent(in) :: line, key
        character(len=*), intent(out) :: value
        integer :: first, last
        call value_bounds(line, key, first, last)
        value = line(first:last)
    end subroutine get_string_value

    subroutine get_real_value(line, key, value, default_value)
        character(len=*), intent(in) :: line, key
        real(real64), intent(out) :: value
        real(real64), intent(in) :: default_value
        integer :: first, last, ios
        logical :: found
        call value_bounds_optional(line, key, first, last, found)
        if (.not. found) then
            value = default_value
            return
        end if
        read(line(first:last), *, iostat=ios) value
        if (ios /= 0) error stop "invalid real value in setup"
    end subroutine get_real_value

    subroutine get_integer_value(line, key, value, default_value)
        character(len=*), intent(in) :: line, key
        integer, intent(out) :: value
        integer, intent(in) :: default_value
        integer :: first, last, ios
        logical :: found
        call value_bounds_optional(line, key, first, last, found)
        if (.not. found) then
            value = default_value
            return
        end if
        read(line(first:last), *, iostat=ios) value
        if (ios /= 0) error stop "invalid integer value in setup"
    end subroutine get_integer_value

    subroutine value_bounds(line, key, first, last)
        character(len=*), intent(in) :: line, key
        integer, intent(out) :: first, last
        logical :: found
        call value_bounds_optional(line, key, first, last, found)
        if (.not. found) error stop "required key is missing from setup line"
    end subroutine value_bounds

    subroutine value_bounds_optional(line, key, first, last, found)
        character(len=*), intent(in) :: line, key
        integer, intent(out) :: first, last
        logical, intent(out) :: found
        character(len=LINE_LENGTH) :: lower_line
        integer :: position, equals, length, candidate, after_key, key_length
        lower_line = lowercase(line)
        length = len_trim(line)
        key_length = len_trim(key)
        position = 0
        do candidate = 1, length - key_length + 1
            if (lower_line(candidate:candidate + key_length - 1) /= lowercase(trim(key))) cycle
            if (candidate > 1) then
                if (line(candidate - 1:candidate - 1) /= " " .and. &
                    line(candidate - 1:candidate - 1) /= char(9)) cycle
            end if
            after_key = candidate + key_length
            do while (after_key <= length .and. &
                      (line(after_key:after_key) == " " .or. line(after_key:after_key) == char(9)))
                after_key = after_key + 1
            end do
            if (after_key <= length .and. line(after_key:after_key) == "=") then
                position = candidate
                exit
            end if
        end do
        found = position > 0
        if (.not. found) then
            first = 0
            last = 0
            return
        end if
        equals = index(line(position + key_length:length), "=")
        if (equals == 0) error stop "key without equals sign in setup"
        first = position + key_length + equals
        do while (first <= length .and. line(first:first) == " ")
            first = first + 1
        end do
        last = first
        do while (last <= length .and. line(last:last) /= " " .and. line(last:last) /= char(9))
            last = last + 1
        end do
        last = last - 1
    end subroutine value_bounds_optional

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

end module accelnet_setup
