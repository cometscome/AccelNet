program test_n2p2_network
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: atomic_structure
    use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: structure
    character(len=1024) :: directory
    real(real64) :: energy, expected, q, dq, fc, dfc, x, forces(3, 2)

    call get_command_argument(1, directory)
    call load_predictor_from_n2p2(trim(directory), model)
    if (size(model%species_names) /= 1 .or. trim(model%species_names(1)) /= "H") &
        error stop "n2p2 species import failed"
    structure%natoms = 1
    allocate(structure%positions(3, 1), source=0.0_real64)
    allocate(structure%species(1), source=1)
    call model%predict_energy_forces(structure, energy, forces(:, 1:1))
    if (abs(energy - 1.75_real64) > 1.0e-12_real64) error stop "n2p2 normalized energy import failed"
    if (maxval(abs(forces(:, 1))) > 1.0e-12_real64) error stop "n2p2 force import failed"

    structure%natoms = 2
    deallocate(structure%positions, structure%species)
    allocate(structure%positions(3, 2), source=0.0_real64)
    allocate(structure%species(2), source=1)
    structure%positions(1, 2) = 2.0_real64
    call model%predict_energy_forces(structure, energy, forces)
    x = (2.0_real64 - 0.2_real64*3.0_real64)/(3.0_real64 - 0.2_real64*3.0_real64)
    fc = (x*(x*(20.0_real64*x - 70.0_real64) + 84.0_real64) - 35.0_real64)*x**4 + 1.0_real64
    dfc = x**3*(x*(x*(140.0_real64*x - 420.0_real64) + 420.0_real64) - 140.0_real64)/2.4_real64
    q = exp(-4.0_real64)*fc
    dq = exp(-4.0_real64)*(dfc - 4.0_real64*fc)
    expected = 3.5_real64 + q
    if (abs(energy - expected) > 2.0e-13_real64) error stop "n2p2 polynomial cutoff energy failed"
    if (abs(forces(1, 1) - dq) > 2.0e-12_real64 .or. &
        abs(forces(1, 2) + dq) > 2.0e-12_real64 .or. &
        maxval(abs(forces(2:3, :))) > 2.0e-12_real64) error stop "n2p2 polynomial cutoff force failed"
end program test_n2p2_network
