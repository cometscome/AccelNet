program test_n2p2_atomic_api
    use iso_c_binding, only: c_int
    use iso_fortran_env, only: real64
    use accelnet
    use accelnet_descriptors, only: atomic_structure, read_n2p2_data
    use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2
    implicit none
    character(len=1024) :: directory, data_file
    character(len=16) :: species(1)
    type(atomic_structure), allocatable :: structures(:)
    type(predictor_model) :: reference_model
    real(real64) :: reference_energy, atomic_energy, api_energy
    real(real64) :: neighbor(3, 1)
    integer(c_int) :: stat
    integer :: atom, other

    call get_command_argument(1, directory)
    call get_command_argument(2, data_file)
    species = [character(len=16) :: "H"]
    call read_n2p2_data(trim(data_file), species, structures)
    call require(size(structures) == 2, "structure count")
    call require(structures(1)%natoms == 2 .and. .not. structures(1)%pbc, "molecular structure")
    call require(structures(2)%natoms == 1 .and. structures(2)%pbc, "periodic structure")
    call require(maxval(abs(structures(2)%lattice - reshape([10.0_real64,0.0_real64,0.0_real64, &
        0.0_real64,9.0_real64,0.0_real64, 0.0_real64,0.0_real64,8.0_real64], [3,3]))) < 1.0e-14_real64, &
        "lattice")

    call load_predictor_from_n2p2(trim(directory), reference_model)
    call reference_model%predict_energy(structures(1), reference_energy)
    call accelnet_init(species, stat)
    call require(stat == ACCELNET_OK, "atomic API init")
    call accelnet_load_n2p2(trim(directory), stat)
    call require(stat == ACCELNET_OK .and. accelnet_all_loaded(), "atomic API n2p2 load")
    api_energy = 0.0_real64
    do atom = 1, 2
        other = 3 - atom
        neighbor(:, 1) = structures(1)%positions(:, other)
        call accelnet_atomic_energy(structures(1)%positions(:, atom), 1_c_int, 1_c_int, &
            neighbor, [1_c_int], atomic_energy, stat)
        call require(stat == ACCELNET_OK, "atomic API evaluation")
        api_energy = api_energy + atomic_energy
    end do
    call require(abs(api_energy - reference_energy) < 1.0e-13_real64, "atomic API energy parity")
    call accelnet_final(stat)
    call require(stat == ACCELNET_OK, "atomic API final")
    call accelnet_init_n2p2(trim(directory), stat)
    call require(stat == ACCELNET_OK .and. accelnet_all_loaded(), "atomic API one-shot n2p2 init")
    call accelnet_final(stat)
    call require(stat == ACCELNET_OK, "atomic API one-shot final")
contains
    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label
        if (.not. condition) then
            write(*, "(2A)") "FAILED: ", trim(label)
            error stop 1
        end if
    end subroutine require
end program test_n2p2_atomic_api
