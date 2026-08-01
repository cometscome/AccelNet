program test_n2p2_network
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: atomic_structure
    use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: structure
    character(len=1024) :: directory
    real(real64) :: energy, forces(3, 1)

    call get_command_argument(1, directory)
    call load_predictor_from_n2p2(trim(directory), model)
    if (size(model%species_names) /= 1 .or. trim(model%species_names(1)) /= "H") &
        error stop "n2p2 species import failed"
    structure%natoms = 1
    allocate(structure%positions(3, 1), source=0.0_real64)
    allocate(structure%species(1), source=1)
    call model%predict_energy_forces(structure, energy, forces)
    if (abs(energy - 2.5_real64) > 1.0e-12_real64) error stop "n2p2 energy import failed"
    if (maxval(abs(forces)) > 1.0e-12_real64) error stop "n2p2 force import failed"
end program test_n2p2_network
