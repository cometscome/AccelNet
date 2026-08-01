program original_accelnet_benchmark
    use iso_fortran_env, only: real64
    use bfio, only: FILEPATHLEN
    use structureinfo, only: Structuredata
    use fromaenet, only: geo_type_conv
    use lclist, only: lcl_init, lcl_final, lcl_nbdist_cart
    use shared_potentials, only: initialize_potentials, AccelNet_atomic_energy_and_forces_novirial, &
        num_of_maxatoms_in_sphere, Rc_min, Rc_max
    implicit none
    type(Structuredata) :: structure
    character(len=2) :: species_names(2)
    character(len=FILEPATHLEN) :: network_files(2)
    character(len=1024) :: xsf_file, argument
    real(real64), allocatable :: neighbor_coordinates(:, :), neighbor_distances(:), forces(:, :)
    integer, allocatable :: neighbor_atoms(:), neighbor_species(:)
    integer :: repeats, iteration
    real(real64) :: energy, start_time, end_time

    if (command_argument_count() /= 4) then
        write(*, "(A)") "usage: original_accelnet_benchmark Ti.nn O.nn structure.xsf REPEATS"
        error stop 2
    end if
    call get_command_argument(1, network_files(1)); call get_command_argument(2, network_files(2))
    call get_command_argument(3, xsf_file); call get_command_argument(4, argument); read(argument, *) repeats
    species_names = [character(len=2) :: "Ti", "O"]
    call initialize_potentials(2, species_names, network_files)
    structure = Structuredata(); call structure%load(trim(xsf_file))
    allocate(neighbor_coordinates(3, num_of_maxatoms_in_sphere), &
             neighbor_distances(num_of_maxatoms_in_sphere), neighbor_atoms(num_of_maxatoms_in_sphere), &
             neighbor_species(num_of_maxatoms_in_sphere), forces(3, structure%num_of_atoms))
    call evaluate_structure(energy, forces)
    call cpu_time(start_time)
    do iteration = 1, repeats
        call evaluate_structure(energy, forces)
    end do
    call cpu_time(end_time)
    write(*, "(A,1X,ES24.16)") "TOTAL_ENERGY_EV", energy
    write(*, "(A,1X,ES24.16)") "FORCE_ABS_CHECKSUM", sum(abs(forces))
    write(*, "(A,1X,I0)") "REPEATS", repeats
    write(*, "(A,1X,ES24.16)") "ENERGY_FORCE_API_SECONDS_PER_STRUCTURE", &
        (end_time - start_time)/repeats

contains

    subroutine evaluate_structure(total_energy, output_forces)
        real(real64), intent(out) :: total_energy
        real(real64), intent(out) :: output_forces(:, :)
        integer :: atom, local, n, status, central_species
        real(real64) :: atomic_energy
        output_forces = 0.0_real64; total_energy = 0.0_real64
        call lcl_init(Rc_min, Rc_max, structure%latticeVec, structure%num_of_atoms, &
                      structure%atomType, structure%abc_ith_atomic_coordinates, structure%pbc)
        do atom = 1, structure%num_of_atoms
            central_species = geo_type_conv(structure%atomType(atom), structure%num_of_types, &
                structure%atomTypeName, 2, species_names)
            n = num_of_maxatoms_in_sphere
            call lcl_nbdist_cart(atom, n, neighbor_coordinates, neighbor_distances, r_cut=Rc_max, &
                                 nblist=neighbor_atoms, nbtype=neighbor_species)
            do local = 1, n
                neighbor_species(local) = geo_type_conv(neighbor_species(local), structure%num_of_types, &
                    structure%atomTypeName, 2, species_names)
            end do
            call AccelNet_atomic_energy_and_forces_novirial( &
                structure%calc_Cartesian_coordinates(atom), central_species, atom, n, &
                neighbor_coordinates, neighbor_species, neighbor_atoms, structure%num_of_atoms, &
                atomic_energy, output_forces, status)
            total_energy = total_energy + atomic_energy
        end do
        call lcl_final()
    end subroutine evaluate_structure
end program original_accelnet_benchmark
