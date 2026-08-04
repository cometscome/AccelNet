program accelnet_predict
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model, load_predictor, load_predictor_from_networks, &
        load_predictor_from_n2p2
    use predict_input, only: predict_input_data, read_predict_input
    use accelnet_descriptors, only: atomic_structure, read_n2p2_data
    implicit none
    type(predictor_model) :: model
    type(predict_input_data) :: input
    type(atomic_structure), allocatable :: data_structures(:)
    character(len=1024), allocatable :: setups(:), networks(:)
    character(len=1024) :: argument, xsf
    integer :: count, species
    real(real64) :: energy
    real(real64), allocatable :: forces(:, :)
    logical :: with_forces
    if (command_argument_count() == 2) then
        call get_command_argument(1, argument)
        if (trim(argument) == "--n2p2" .or. trim(argument) == "--n2p2-data") then
            call get_command_argument(2, xsf)
            call load_predictor_from_n2p2(".", model)
            if (trim(argument) == "--n2p2-data") then
                call predict_data_file(trim(xsf))
            else
                call model%predict_energy_forces(trim(xsf), energy, forces)
                call print_prediction(trim(xsf), energy, forces)
            end if
            stop
        end if
    elseif (command_argument_count() == 3) then
        call get_command_argument(1, argument)
        if (trim(argument) == "--n2p2" .or. trim(argument) == "--n2p2-data") then
            call get_command_argument(2, xsf)
            call load_predictor_from_n2p2(trim(xsf), model)
            call get_command_argument(3, xsf)
            if (trim(argument) == "--n2p2-data") then
                call predict_data_file(trim(xsf))
            else
                call model%predict_energy_forces(trim(xsf), energy, forces)
                call print_prediction(trim(xsf), energy, forces)
            end if
            stop
        end if
    end if
    if (command_argument_count() == 1) then
        call get_command_argument(1, argument)
        call read_predict_input(trim(argument), input)
        call load_predictor_from_networks(input%networks, model, input%chebyshev_version)
        do species = 1, size(input%structures)
            call model%predict_energy_forces(trim(input%structures(species)), energy, forces)
            call print_prediction(trim(input%structures(species)), energy, forces)
        end do
        stop
    elseif (command_argument_count() < 4) then
        write(*, "(A)") "usage: accelnet-predict NSPECIES SETUP... NETWORK... XSF"
        write(*, "(A)") "   or: accelnet-predict predict.in"
        write(*, "(A)") "   or: accelnet-predict --n2p2 [MODEL_DIR] XSF"
        write(*, "(A)") "   or: accelnet-predict --n2p2-data [MODEL_DIR] input.data"
        error stop 2
    end if
    call get_command_argument(1, argument); read(argument, *) count
    with_forces = command_argument_count() == 2*count + 3
    if (.not. with_forces .and. command_argument_count() /= 2*count + 2) error stop "invalid number of arguments"
    allocate(setups(count), networks(count))
    do species = 1, count
        call get_command_argument(1 + species, setups(species))
        call get_command_argument(1 + count + species, networks(species))
    end do
    call get_command_argument(2*count + 2, xsf)
    call load_predictor(setups, networks, model)
    if (with_forces) then
        call get_command_argument(2*count + 3, argument)
        if (trim(argument) /= "--forces") error stop "unknown option"
        call model%predict_energy_forces(xsf, energy, forces)
    else
        call model%predict_energy(xsf, energy)
    end if
    write(*, "(A,1X,ES24.16)") "TOTAL_ENERGY_EV", energy
    if (with_forces) then
        do species = 1, size(forces, 2)
            write(*, "(A,1X,I0,3(1X,ES24.16))") "FORCE_EV_PER_ANGSTROM", species, forces(:, species)
        end do
    end if
contains
    subroutine predict_data_file(filename)
        character(len=*), intent(in) :: filename
        integer :: structure_index
        call read_n2p2_data(filename, model%species_names, data_structures)
        do structure_index = 1, size(data_structures)
            allocate(forces(3, data_structures(structure_index)%natoms))
            call model%predict_energy_forces(data_structures(structure_index), energy, forces)
            call print_prediction(trim(filename)//":"//integer_string(structure_index), energy, forces)
            deallocate(forces)
        end do
    end subroutine predict_data_file

    function integer_string(value) result(text)
        integer, intent(in) :: value
        character(len=32) :: text
        write(text, "(I0)") value
    end function integer_string

    subroutine print_prediction(structure_file, predicted_energy, predicted_forces)
        character(len=*), intent(in) :: structure_file
        real(real64), intent(in) :: predicted_energy
        real(real64), intent(in) :: predicted_forces(:, :)
        integer :: atom
        write(*, "(A,1X,A)") "STRUCTURE", trim(structure_file)
        write(*, "(A,1X,ES24.16)") "TOTAL_ENERGY_EV", predicted_energy
        do atom = 1, size(predicted_forces, 2)
            write(*, "(A,1X,I0,3(1X,ES24.16))") &
                "FORCE_EV_PER_ANGSTROM", atom, predicted_forces(:, atom)
        end do
    end subroutine print_prediction
end program accelnet_predict
