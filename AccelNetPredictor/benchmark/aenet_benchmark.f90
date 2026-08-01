program aenet_benchmark
    use iso_fortran_env, only: real64
    use aenet, only: aenet_init, aenet_final, aenet_load_potential, aenet_convert_atom_types, &
        aenet_atomic_energy, aenet_atomic_energy_and_forces, aenet_nnb_max, aenet_Rc_min, aenet_Rc_max
    use geometry, only: geo_init, geo_final, pbc, latticeVec, nAtoms, nTypes, &
        atomType, atomTypeName, cooLatt
    use lclist, only: lcl_init, lcl_final, lcl_nbdist_cart
    implicit none
    character(len=2) :: species_names(2)
    character(len=1024) :: network_files(2), xsf_file, argument
    integer, allocatable :: original_species(:), neighbor_atoms(:), neighbor_species(:)
    real(real64), allocatable :: neighbor_coordinates(:, :), neighbor_distances(:), forces(:, :)
    integer :: repeats, iteration, status, species
    logical :: networks_are_ascii
    real(real64) :: energy, start_time, end_time, energy_time, force_time

    if (command_argument_count() /= 4) then
        write(*, "(A)") "usage: aenet-benchmark Ti.nn.ascii O.nn.ascii structure.xsf REPEATS"
        error stop 2
    end if
    call get_command_argument(1, network_files(1)); call get_command_argument(2, network_files(2))
    call get_command_argument(3, xsf_file); call get_command_argument(4, argument); read(argument, *) repeats
    networks_are_ascii = index(trim(network_files(1)), ".ascii") > 0
    species_names = [character(len=2) :: "Ti", "O"]
    call aenet_init(species_names, status); call check_status("aenet_init", status)
    do species = 1, 2
        call aenet_load_potential(species, trim(network_files(species)), status, is_ascii=networks_are_ascii)
        call check_status("aenet_load_potential", status)
    end do
    call geo_init(trim(xsf_file), "xsf")
    allocate(original_species(nAtoms)); original_species = atomType
    call aenet_convert_atom_types(atomTypeName, original_species, atomType, status)
    call check_status("aenet_convert_atom_types", status)
    allocate(neighbor_coordinates(3, aenet_nnb_max), neighbor_distances(aenet_nnb_max), &
             neighbor_atoms(aenet_nnb_max), neighbor_species(aenet_nnb_max), forces(3, nAtoms))

    call evaluate_energy(energy)
    call cpu_time(start_time)
    do iteration = 1, repeats
        call evaluate_energy(energy)
    end do
    call cpu_time(end_time); energy_time = end_time - start_time
    call evaluate_energy_forces(energy, forces)
    call cpu_time(start_time)
    do iteration = 1, repeats
        call evaluate_energy_forces(energy, forces)
    end do
    call cpu_time(end_time); force_time = end_time - start_time
    write(*, "(A,1X,ES24.16)") "TOTAL_ENERGY_EV", energy
    write(*, "(A,1X,ES24.16)") "FORCE_ABS_CHECKSUM", sum(abs(forces))
    write(*, "(A,1X,I0)") "REPEATS", repeats
    write(*, "(A,1X,ES24.16)") "ENERGY_API_SECONDS_PER_STRUCTURE", energy_time/repeats
    write(*, "(A,1X,ES24.16)") "ENERGY_FORCE_API_SECONDS_PER_STRUCTURE", force_time/repeats
    call geo_final(); call aenet_final(status)

contains

    subroutine evaluate_energy(total_energy)
        real(real64), intent(out) :: total_energy
        integer :: atom, n, local_status
        real(real64) :: atomic_energy, center(3)
        total_energy = 0.0_real64
        call lcl_init(aenet_Rc_min, aenet_Rc_max, latticeVec, nAtoms, atomType, cooLatt, pbc)
        do atom = 1, nAtoms
            center = matmul(latticeVec, cooLatt(:, atom)); n = aenet_nnb_max
            call lcl_nbdist_cart(atom, n, neighbor_coordinates, neighbor_distances, aenet_Rc_max, &
                                 nblist=neighbor_atoms, nbtype=neighbor_species)
            call aenet_atomic_energy(center, atomType(atom), n, neighbor_coordinates, &
                                     neighbor_species, atomic_energy, local_status)
            call check_status("aenet_atomic_energy", local_status)
            total_energy = total_energy + atomic_energy
        end do
        call lcl_final()
    end subroutine evaluate_energy

    subroutine evaluate_energy_forces(total_energy, output_forces)
        real(real64), intent(out) :: total_energy
        real(real64), intent(out) :: output_forces(:, :)
        integer :: atom, n, local_status
        real(real64) :: atomic_energy, center(3)
        total_energy = 0.0_real64; output_forces = 0.0_real64
        call lcl_init(aenet_Rc_min, aenet_Rc_max, latticeVec, nAtoms, atomType, cooLatt, pbc)
        do atom = 1, nAtoms
            center = matmul(latticeVec, cooLatt(:, atom)); n = aenet_nnb_max
            call lcl_nbdist_cart(atom, n, neighbor_coordinates, neighbor_distances, aenet_Rc_max, &
                                 nblist=neighbor_atoms, nbtype=neighbor_species)
            call aenet_atomic_energy_and_forces(center, atomType(atom), atom, n, neighbor_coordinates, &
                neighbor_species, neighbor_atoms, nAtoms, atomic_energy, output_forces, local_status)
            call check_status("aenet_atomic_energy_and_forces", local_status)
            total_energy = total_energy + atomic_energy
        end do
        call lcl_final()
    end subroutine evaluate_energy_forces

    subroutine check_status(operation, value)
        character(len=*), intent(in) :: operation
        integer, intent(in) :: value
        if (value /= 0) then
            write(*, "(A,1X,A,1X,I0)") "aenet operation failed:", trim(operation), value
            error stop "aenet benchmark failed"
        end if
    end subroutine check_status
end program aenet_benchmark
