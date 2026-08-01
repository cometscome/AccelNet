module predict_input
    implicit none
    private
    integer, parameter :: PATH_LENGTH = 1024
    type, public :: predict_input_data
        character(len=16), allocatable :: species(:)
        character(len=PATH_LENGTH), allocatable :: networks(:)
        character(len=PATH_LENGTH), allocatable :: structures(:)
        integer :: chebyshev_version = 0
    end type
    public :: read_predict_input
contains
    subroutine read_predict_input(filename, input)
        character(len=*), intent(in) :: filename
        type(predict_input_data), intent(out) :: input
        integer :: unit, ios, count, i, position
        character(len=PATH_LENGTH) :: line, species_name, path
        logical :: found
        open(newunit=unit, file=trim(filename), status="old", action="read", iostat=ios)
        if (ios /= 0) error stop "cannot open predict.in"
        call read_optional_version(unit, input%chebyshev_version)
        rewind(unit)
        call find_section(unit, "TYPES", found)
        if (.not. found) error stop "predict.in has no TYPES section"
        call next_data_line(unit, line, ios); read(line, *) count
        allocate(input%species(count))
        do i = 1, count
            call next_data_line(unit, line, ios); read(line, *) input%species(i)
        end do
        rewind(unit); call find_section(unit, "NETWORKS", found)
        if (.not. found) error stop "predict.in has no NETWORKS section"
        allocate(input%networks(count))
        do i = 1, count
            call next_data_line(unit, line, ios)
            position = scan(trim(line), " "//achar(9))
            if (position == 0) error stop "invalid NETWORKS entry in predict.in"
            species_name = trim(line(:position - 1))
            path = adjustl(line(position + 1:))
            if (len_trim(path) >= 2) then
                if ((path(1:1) == '"' .and. path(len_trim(path):len_trim(path)) == '"') .or. &
                    (path(1:1) == "'" .and. path(len_trim(path):len_trim(path)) == "'")) then
                    path = path(2:len_trim(path) - 1)
                end if
            end if
            input%networks(i) = trim(path)
            if (trim(species_name) /= trim(input%species(i))) error stop "TYPES and NETWORKS ordering differs"
        end do
        rewind(unit); call find_section(unit, "FILES", found)
        if (.not. found) error stop "predict.in has no FILES section"
        call next_data_line(unit, line, ios); read(line, *) count
        allocate(input%structures(count))
        do i = 1, count
            call next_data_line(unit, line, ios); input%structures(i) = trim(line)
        end do
        close(unit)
    end subroutine read_predict_input

    subroutine read_optional_version(unit, version)
        integer, intent(in) :: unit
        integer, intent(out) :: version
        character(len=PATH_LENGTH) :: line, key, value
        integer :: ios, separator
        version = 0 ! ænet 2.04 and later Chebyshev convention.
        rewind(unit)
        do
            call next_data_line(unit, line, ios)
            if (ios /= 0) exit
            separator = index(line, "=")
            if (separator > 0) then
                key = uppercase(trim(adjustl(line(:separator - 1))))
                value = adjustl(line(separator + 1:))
            else
                separator = scan(trim(line), " "//achar(9))
                if (separator > 0) then
                    key = uppercase(trim(line(:separator - 1)))
                    value = adjustl(line(separator + 1:))
                else
                    key = uppercase(trim(line))
                    value = ""
                end if
            end if
            if (trim(key) == "VERSION" .or. trim(key) == "CHEBYSHEV_VERSION") then
                if (len_trim(value) == 0) call next_data_line(unit, value, ios)
                if (ios /= 0 .or. len_trim(value) == 0) error stop "VERSION in predict.in has no value"
                read(value, *, iostat=ios) version
                if (ios /= 0) error stop "invalid VERSION value in predict.in"
                if (version /= 0 .and. version /= 1 .and. version /= 10) &
                    error stop "predict.in VERSION must be 0, 1, or 10"
                return
            end if
        end do
    end subroutine read_optional_version

    subroutine find_section(unit, requested, found)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: requested
        logical, intent(out) :: found
        character(len=PATH_LENGTH) :: line
        integer :: ios
        found = .false.
        do
            call next_data_line(unit, line, ios)
            if (ios /= 0) return
            if (uppercase(trim(line)) == uppercase(trim(requested))) then
                found = .true.; return
            end if
        end do
    end subroutine find_section

    subroutine next_data_line(unit, line, ios)
        integer, intent(in) :: unit
        character(len=*), intent(out) :: line
        integer, intent(out) :: ios
        integer :: comment1, comment2, cut
        do
            read(unit, "(A)", iostat=ios) line
            if (ios /= 0) return
            comment1 = index(line, "!"); comment2 = index(line, "#")
            cut = len(line) + 1
            if (comment1 > 0) cut = min(cut, comment1)
            if (comment2 > 0) cut = min(cut, comment2)
            if (cut <= len(line)) line(cut:) = " "
            line = adjustl(line)
            if (len_trim(line) > 0) return
        end do
    end subroutine next_data_line

    pure function uppercase(input) result(output)
        character(len=*), intent(in) :: input
        character(len=len(input)) :: output
        integer :: i, code
        output = input
        do i = 1, len(input)
            code = iachar(input(i:i))
            if (code >= iachar('a') .and. code <= iachar('z')) output(i:i) = achar(code - 32)
        end do
    end function uppercase
end module predict_input
