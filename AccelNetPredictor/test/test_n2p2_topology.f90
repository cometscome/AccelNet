program test_n2p2_topology
    use iso_fortran_env, only: real64
    use aenet_network, only: atomic_network, ACTIVATION_LINEAR, ACTIVATION_TANH, ACTIVATION_SOFTPLUS
    use accelnet_descriptors, only: atomic_structure
    use accelnet_setup, only: descriptor_setup
    use n2p2_network, only: load_n2p2_model
    use accelnet_predictor, only: predictor_model, load_predictor_from_n2p2
    implicit none

    character(len=1024) :: directory, depth_directory
    character(len=16), allocatable :: species_names(:)
    type(atomic_network), allocatable :: networks(:)
    type(descriptor_setup), allocatable :: setups(:)
    type(predictor_model) :: model
    type(atomic_structure) :: structure
    real(real64) :: energy, forces(3, 3)
    real(real64), parameter :: expected_energy = 4.6166056085325979e-2_real64
    real(real64), parameter :: expected_forces(3, 3) = reshape([ &
         4.7281397979032514e-3_real64,  7.5454842542942921e-3_real64,  1.1515651384941376e-3_real64, &
        -4.3488915459515285e-3_real64,  1.4308173250058413e-3_real64, -6.2640936984803124e-5_real64, &
        -3.7924825195172348e-4_real64, -8.9763015793001347e-3_real64, -1.0889242015093344e-3_real64], [3, 3])

    call get_command_argument(1, directory)
    call load_n2p2_model(trim(directory), networks, setups, species_names)

    if (size(networks) /= 2 .or. trim(species_names(1)) /= "H" .or. trim(species_names(2)) /= "O") &
        error stop "per-element topology species import failed"
    if (any(networks(1)%nodes /= [2, 2, 2, 1])) error stop "global H topology import failed"
    if (any(networks(2)%nodes /= [2, 3, 2, 1])) error stop "per-element O node topology import failed"
    if (any(networks(1)%activation /= [ACTIVATION_TANH, ACTIVATION_TANH, ACTIVATION_LINEAR])) &
        error stop "global H activation topology import failed"
    if (any(networks(2)%activation /= [ACTIVATION_SOFTPLUS, ACTIVATION_TANH, ACTIVATION_LINEAR])) &
        error stop "per-element O activation topology import failed"

    ! n2p2 normalize_nodes divides every weight and bias feeding a layer by
    ! the number of nodes in the preceding layer. The importer folds those
    ! factors into the stored connections.
    if (abs(networks(1)%weights(1) - 0.30_real64/2.0_real64) > 1.0e-15_real64 .or. &
        abs(networks(1)%weights(7) - 0.70_real64/2.0_real64) > 1.0e-15_real64 .or. &
        abs(networks(1)%weights(13) - (-0.15_real64)/2.0_real64) > 1.0e-15_real64 .or. &
        abs(networks(2)%weights(1) - 0.25_real64/2.0_real64) > 1.0e-15_real64 .or. &
        abs(networks(2)%weights(10) - 0.40_real64/3.0_real64) > 1.0e-15_real64 .or. &
        abs(networks(2)%weights(18) - 0.55_real64/2.0_real64) > 1.0e-15_real64) &
        error stop "normalize_nodes connection import failed"

    call load_predictor_from_n2p2(trim(directory), model)
    structure%natoms = 3
    allocate(structure%positions(3, 3), structure%species(3))
    structure%positions = reshape([ &
        0.0_real64, 0.0_real64, 0.0_real64, &
        1.2_real64, 0.3_real64, 0.1_real64, &
        0.4_real64, 1.5_real64, 0.2_real64], [3, 3])
    structure%species = [1, 2, 1]
    call model%predict_energy_forces(structure, energy, forces)
    if (abs(energy - expected_energy) > 2.0e-14_real64) &
        error stop "per-element topology/normalize_nodes energy differs from n2p2"
    if (maxval(abs(forces - expected_forces)) > 2.0e-14_real64) &
        error stop "per-element topology/normalize_nodes forces differ from n2p2"

    ! Exercise a genuine per-element hidden-layer-count override separately.
    ! n2p2 v2.3.0 has an upstream setup bug for this particular combination,
    ! so the end-to-end numerical reference above uses equal depth and
    ! different widths/activations, while this fixture covers parser depth.
    call get_command_argument(2, depth_directory)
    call load_n2p2_model(trim(depth_directory), networks, setups, species_names)
    if (any(networks(1)%nodes /= [1, 1])) error stop "zero-hidden-layer global topology failed"
    if (any(networks(2)%nodes /= [1, 2, 1])) error stop "per-element hidden-layer override failed"
end program test_n2p2_topology
