program aenet_corpus_benchmark
    use iso_fortran_env, only: real64
    use aenet, only: aenet_init, aenet_final, aenet_load_potential, aenet_convert_atom_types, &
        aenet_atomic_energy, aenet_atomic_energy_and_forces, aenet_nnb_max, aenet_Rc_min, aenet_Rc_max
    use geometry, only: geo_init, geo_final, pbc, latticeVec, nAtoms, atomType, atomTypeName, cooLatt
    use lclist, only: lcl_init, lcl_final, lcl_nbdist_cart
    implicit none

    type :: corpus_structure
        integer :: natoms = 0
        logical :: periodic = .false.
        real(real64) :: lattice(3, 3) = 0.0_real64
        real(real64), allocatable :: fractional(:, :)
        integer, allocatable :: species(:)
    end type corpus_structure

    type(corpus_structure), allocatable :: structures(:)
    character(len=2) :: species_names(2)
    character(len=1024) :: network_files(2), corpus_dir, argument, filename
    integer, allocatable :: neighbor_atoms(:), neighbor_species(:)
    real(real64), allocatable :: neighbor_coordinates(:, :), neighbor_distances(:), forces(:, :)
    integer :: repeats, nstructures, iteration, structure_index, status, species, max_atoms
    real(real64) :: energy, energy_checksum, force_checksum
    real(real64) :: start_time, end_time, energy_time, force_time

    if (command_argument_count() /= 5) then
        write(*, "(A)") &
            "usage: aenet-corpus-benchmark Ti.nn.ascii O.nn.ascii DIRECTORY NSTRUCTURES REPEATS"
        error stop 2
    end if
    call get_command_argument(1, network_files(1))
    call get_command_argument(2, network_files(2))
    call get_command_argument(3, corpus_dir)
    call get_command_argument(4, argument); read(argument, *) nstructures
    call get_command_argument(5, argument); read(argument, *) repeats
    if (nstructures < 1 .or. nstructures > 9999) error stop "invalid NSTRUCTURES"
    if (repeats < 1) error stop "REPEATS must be positive"

    species_names = [character(len=2) :: "Ti", "O"]
    call aenet_init(species_names, status); call check_status("aenet_init", status)
    do species = 1, 2
        call aenet_load_potential(species, trim(network_files(species)), status, is_ascii=.true.)
        call check_status("aenet_load_potential", status)
    end do

    allocate(structures(nstructures))
    max_atoms = 0
    do structure_index = 1, nstructures
        call structure_filename(corpus_dir, structure_index, filename)
        call load_structure(trim(filename), structures(structure_index))
        max_atoms = max(max_atoms, structures(structure_index)%natoms)
    end do
    allocate(neighbor_coordinates(3, aenet_nnb_max), neighbor_distances(aenet_nnb_max), &
             neighbor_atoms(aenet_nnb_max), neighbor_species(aenet_nnb_max), forces(3, max_atoms))

    energy_checksum = 0.0_real64
    force_checksum = 0.0_real64
    do structure_index = 1, nstructures
        call evaluate_energy_forces(structures(structure_index), energy, &
            forces(:, 1:structures(structure_index)%natoms))
        energy_checksum = energy_checksum + energy
        force_checksum = force_checksum + sum(abs(forces(:, 1:structures(structure_index)%natoms)))
    end do

    call cpu_time(start_time)
    do iteration = 1, repeats
        do structure_index = 1, nstructures
            call evaluate_energy(structures(structure_index), energy)
        end do
    end do
    call cpu_time(end_time)
    energy_time = end_time - start_time

    call cpu_time(start_time)
    do iteration = 1, repeats
        do structure_index = 1, nstructures
            call evaluate_energy_forces(structures(structure_index), energy, &
                forces(:, 1:structures(structure_index)%natoms))
        end do
    end do
    call cpu_time(end_time)
    force_time = end_time - start_time

    write(*, "(A,1X,I0)") "STRUCTURES", nstructures
    write(*, "(A,1X,I0)") "REPEATS", repeats
    write(*, "(A,1X,ES24.16)") "ENERGY_CHECKSUM_EV", energy_checksum
    write(*, "(A,1X,ES24.16)") "FORCE_ABS_CHECKSUM", force_checksum
    write(*, "(A,1X,ES24.16)") "ENERGY_SECONDS_PER_STRUCTURE", &
        energy_time/real(nstructures*repeats, real64)
    write(*, "(A,1X,ES24.16)") "ENERGY_FORCE_SECONDS_PER_STRUCTURE", &
        force_time/real(nstructures*repeats, real64)
    call aenet_final(status)

contains
    subroutine load_structure(path, structure)
        character(len=*), intent(in) :: path
        type(corpus_structure), intent(out) :: structure
        integer, allocatable :: original_species(:)
        call geo_init(path, "xsf")
        structure%natoms = nAtoms
        structure%periodic = pbc
        structure%lattice = latticeVec
        allocate(structure%fractional(3, nAtoms), structure%species(nAtoms), original_species(nAtoms))
        structure%fractional = cooLatt
        original_species = atomType
        call aenet_convert_atom_types(atomTypeName, original_species, structure%species, status)
        call check_status("aenet_convert_atom_types", status)
        call geo_final()
    end subroutine load_structure

    subroutine evaluate_energy(structure, total_energy)
        type(corpus_structure), intent(in) :: structure
        real(real64), intent(out) :: total_energy
        integer :: atom, n, local_status
        real(real64) :: atomic_energy, center(3)
        total_energy = 0.0_real64
        call lcl_init(aenet_Rc_min, aenet_Rc_max, structure%lattice, structure%natoms, &
            structure%species, structure%fractional, structure%periodic)
        do atom = 1, structure%natoms
            center = matmul(structure%lattice, structure%fractional(:, atom))
            n = aenet_nnb_max
            call lcl_nbdist_cart(atom, n, neighbor_coordinates, neighbor_distances, aenet_Rc_max, &
                                 nblist=neighbor_atoms, nbtype=neighbor_species)
            call aenet_atomic_energy(center, structure%species(atom), n, neighbor_coordinates, &
                                     neighbor_species, atomic_energy, local_status)
            call check_status("aenet_atomic_energy", local_status)
            total_energy = total_energy + atomic_energy
        end do
        call lcl_final()
    end subroutine evaluate_energy

    subroutine evaluate_energy_forces(structure, total_energy, output_forces)
        type(corpus_structure), intent(in) :: structure
        real(real64), intent(out) :: total_energy
        real(real64), intent(out) :: output_forces(:, :)
        integer :: atom, n, local_status
        real(real64) :: atomic_energy, center(3)
        total_energy = 0.0_real64
        output_forces = 0.0_real64
        call lcl_init(aenet_Rc_min, aenet_Rc_max, structure%lattice, structure%natoms, &
            structure%species, structure%fractional, structure%periodic)
        do atom = 1, structure%natoms
            center = matmul(structure%lattice, structure%fractional(:, atom))
            n = aenet_nnb_max
            call lcl_nbdist_cart(atom, n, neighbor_coordinates, neighbor_distances, aenet_Rc_max, &
                                 nblist=neighbor_atoms, nbtype=neighbor_species)
            call aenet_atomic_energy_and_forces(center, structure%species(atom), atom, n, &
                neighbor_coordinates, neighbor_species, neighbor_atoms, structure%natoms, &
                atomic_energy, output_forces, local_status)
            call check_status("aenet_atomic_energy_and_forces", local_status)
            total_energy = total_energy + atomic_energy
        end do
        call lcl_final()
    end subroutine evaluate_energy_forces

    subroutine structure_filename(directory, index, path)
        character(len=*), intent(in) :: directory
        integer, intent(in) :: index
        character(len=*), intent(out) :: path
        character(len=4) :: number
        write(number, "(I4.4)") index
        path = trim(directory)//"/structure"//number//".xsf"
    end subroutine structure_filename

    subroutine check_status(operation, value)
        character(len=*), intent(in) :: operation
        integer, intent(in) :: value
        if (value /= 0) then
            write(*, "(A,1X,A,1X,I0)") "aenet operation failed:", trim(operation), value
            error stop "aenet corpus benchmark failed"
        end if
    end subroutine check_status
end program aenet_corpus_benchmark
