program test_predictor
    use iso_fortran_env, only: real64
    use accelnet_predictor, only: predictor_model, load_predictor, load_predictor_from_networks
    use accelnet_descriptors, only: atomic_structure, read_xsf
    implicit none
    type(predictor_model) :: model
    type(atomic_structure) :: structure
    character(len=1024) :: setup_files(2), network_files(2), energy_xsf, force_xsf, force_reference
    real(real64), allocatable :: forces(:, :)
    real(real64) :: energy, expected_total, expected_force(3, 23)
    integer :: unit, atom
    call get_command_argument(1, setup_files(1)); call get_command_argument(2, setup_files(2))
    call get_command_argument(3, network_files(1)); call get_command_argument(4, network_files(2))
    call get_command_argument(5, energy_xsf); call get_command_argument(6, force_xsf)
    call get_command_argument(7, force_reference)
    call load_predictor(setup_files, network_files, model)
    call model%predict_energy(energy_xsf, energy)
    expected_total = -203.874928_real64 + 8*(-1604.604515075_real64) + 16*(-432.503149303_real64)
    if (abs(energy - expected_total) > 5.0e-6_real64) error stop "energy golden mismatch"
    call read_xsf(trim(energy_xsf), model%species_names, structure)
    call model%predict_energy(structure, energy)
    if (abs(energy - expected_total) > 5.0e-6_real64) error stop "in-memory structure energy mismatch"
    call model%predict_energy_forces(force_xsf, energy, forces)
    if (abs((energy - 8*(-1604.604515075_real64) - 15*(-432.503149303_real64)) - &
            (-195.474317_real64)) > 5.0e-6_real64) error stop "force-structure energy mismatch"
    open(newunit=unit, file=trim(force_reference), status="old", action="read")
    do atom = 1, 23
        read(unit, *) expected_force(:, atom)
    end do
    close(unit)
    if (maxval(abs(forces - expected_force)) > 1.0e-6_real64) error stop "force golden mismatch"
    if (maxval(abs(sum(forces, dim=2))) > 2.0e-12_real64) error stop "net force is nonzero"
    call read_xsf(trim(force_xsf), model%species_names, structure)
    deallocate(forces); allocate(forces(3, structure%natoms))
    call model%predict_energy_forces(structure, energy, forces)
    if (maxval(abs(forces - expected_force)) > 1.0e-6_real64) &
        error stop "in-memory structure force mismatch"
    call load_predictor_from_networks(network_files, model)
    call model%predict_energy(energy_xsf, energy)
    if (abs(energy - expected_total) > 5.0e-6_real64) error stop "embedded setup energy mismatch"
    call model%predict_energy_forces(force_xsf, energy, forces)
    if (abs((energy - 8*(-1604.604515075_real64) - 15*(-432.503149303_real64)) - &
            (-195.474317_real64)) > 5.0e-6_real64) error stop "embedded setup force-structure energy mismatch"
    if (maxval(abs(forces - expected_force)) > 1.0e-6_real64) &
        error stop "embedded setup force mismatch"
    call model%reload(setup_files, network_files)
    call model%predict_energy(energy_xsf, energy)
    if (abs(energy - expected_total) > 5.0e-6_real64) error stop "setup reload energy mismatch"
end program test_predictor
