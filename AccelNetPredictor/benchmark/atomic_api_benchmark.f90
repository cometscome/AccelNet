program atomic_api_benchmark
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: atomic_structure, neighbor_data, read_xsf, build_neighbor_list
    use accelnet, only: accel_init => accelnet_init, accel_final => accelnet_final, &
        accel_load => accelnet_load_potential, accel_energy => accelnet_atomic_energy, &
        accel_energy_forces => accelnet_atomic_energy_and_forces, &
        accel_rc_min => accelnet_Rc_min, accel_rc_max => accelnet_Rc_max
    use aenet, only: ref_init => aenet_init, ref_final => aenet_final, &
        ref_load => aenet_load_potential, ref_energy => aenet_atomic_energy, &
        ref_energy_forces => aenet_atomic_energy_and_forces, &
        ref_rc_min => aenet_Rc_min, ref_rc_max => aenet_Rc_max
    implicit none

    type :: local_environment
        real(real64), allocatable :: coordinates(:,:)
        integer, allocatable :: species(:), indices(:)
    end type local_environment

    character(len=16) :: species_names(2)
    character(len=1024) :: networks(2), xsf_file, argument
    type(atomic_structure) :: structure, primitive_structure
    type(neighbor_data) :: neighbors
    type(local_environment), allocatable :: environments(:)
    real(real64), allocatable :: accel_forces(:,:), ref_forces(:,:)
    real(real64) :: accel_energy_value, ref_energy_value
    real(real64) :: accel_energy_seconds, ref_energy_seconds
    real(real64) :: accel_force_seconds, ref_force_seconds
    real(real64) :: start_time, end_time, energy_error, force_error
    integer :: repeats, iteration, status, atom, first, last, n, total_neighbors
    integer :: replicate_x, replicate_y, replicate_z
    logical :: networks_are_ascii

    if (command_argument_count() /= 4 .and. command_argument_count() /= 7) then
        write(*,'(A)') "usage: atomic-api-benchmark Ti.nn O.nn structure.xsf REPEATS [NX NY NZ]"
        error stop 2
    end if
    call get_command_argument(1, networks(1))
    call get_command_argument(2, networks(2))
    call get_command_argument(3, xsf_file)
    call get_command_argument(4, argument); read(argument,*) repeats
    if (repeats < 1) error stop "REPEATS must be positive"
    replicate_x = 1; replicate_y = 1; replicate_z = 1
    if (command_argument_count() == 7) then
        call get_command_argument(5, argument); read(argument,*) replicate_x
        call get_command_argument(6, argument); read(argument,*) replicate_y
        call get_command_argument(7, argument); read(argument,*) replicate_z
        if (min(replicate_x,replicate_y,replicate_z) < 1) error stop "supercell sizes must be positive"
    end if
    networks_are_ascii = index(trim(networks(1)), ".ascii") > 0
    species_names = [character(len=16) :: "Ti", "O"]

    call accel_init(species_names, status); call check_status("accelnet_init", status)
    call ref_init(species_names, status); call check_status("aenet_init", status)
    do atom = 1, 2
        call accel_load(atom, trim(networks(atom)), status, is_ascii=networks_are_ascii)
        call check_status("accelnet_load_potential", status)
        call ref_load(atom, trim(networks(atom)), status, is_ascii=networks_are_ascii)
        call check_status("aenet_load_potential", status)
    end do
    if (abs(accel_rc_max-ref_rc_max) > 1.0e-12_real64 .or. &
        abs(accel_rc_min-ref_rc_min) > 1.0e-12_real64) error stop "model cutoff mismatch"

    call read_xsf(trim(xsf_file), species_names, structure)
    if (replicate_x*replicate_y*replicate_z > 1) then
        primitive_structure = structure
        call make_supercell(primitive_structure,replicate_x,replicate_y,replicate_z,structure)
    end if
    call build_neighbor_list(structure, real(accel_rc_max,real64), neighbors, &
                             real(accel_rc_min,real64), preserve_legacy_order=.true.)
    allocate(environments(structure%natoms))
    total_neighbors = 0
    do atom = 1, structure%natoms
        first = neighbors%offsets(atom)
        last = neighbors%offsets(atom+1)-1
        n = max(0,last-first+1)
        allocate(environments(atom)%coordinates(3,n), environments(atom)%species(n), &
                 environments(atom)%indices(n))
        if (n > 0) then
            environments(atom)%coordinates = neighbors%positions(:,first:last)
            environments(atom)%indices = neighbors%atom_indices(first:last)
            environments(atom)%species = structure%species(environments(atom)%indices)
        end if
        total_neighbors = total_neighbors+n
    end do
    allocate(accel_forces(3,structure%natoms), ref_forces(3,structure%natoms))

    ! Warm up every path before measuring. Neighbor construction and model loading
    ! are deliberately outside all timed regions.
    call evaluate_accel_energy(accel_energy_value)
    call evaluate_ref_energy(ref_energy_value)
    call evaluate_accel_forces(accel_energy_value,accel_forces)
    call evaluate_ref_forces(ref_energy_value,ref_forces)
    energy_error = abs(accel_energy_value-ref_energy_value)
    force_error = maxval(abs(accel_forces-ref_forces))
    if (energy_error > 1.0e-8_real64) error stop "atomic API energy mismatch"
    if (force_error > 1.0e-8_real64) error stop "atomic API force mismatch"

    call cpu_time(start_time)
    do iteration = 1, repeats
        call evaluate_accel_energy(accel_energy_value)
    end do
    call cpu_time(end_time); accel_energy_seconds = (end_time-start_time)/repeats
    call cpu_time(start_time)
    do iteration = 1, repeats
        call evaluate_ref_energy(ref_energy_value)
    end do
    call cpu_time(end_time); ref_energy_seconds = (end_time-start_time)/repeats
    call cpu_time(start_time)
    do iteration = 1, repeats
        call evaluate_accel_forces(accel_energy_value,accel_forces)
    end do
    call cpu_time(end_time); accel_force_seconds = (end_time-start_time)/repeats
    call cpu_time(start_time)
    do iteration = 1, repeats
        call evaluate_ref_forces(ref_energy_value,ref_forces)
    end do
    call cpu_time(end_time); ref_force_seconds = (end_time-start_time)/repeats

    write(*,'(A,1X,I0)') "NATOMS", structure%natoms
    write(*,'(A,3(1X,I0))') "SUPERCELL", replicate_x, replicate_y, replicate_z
    write(*,'(A,1X,I0)') "TOTAL_NEIGHBORS", total_neighbors
    write(*,'(A,1X,I0)') "REPEATS", repeats
    write(*,'(A,1X,ES24.16)') "ENERGY_ABS_ERROR", energy_error
    write(*,'(A,1X,ES24.16)') "FORCE_MAX_ABS_ERROR", force_error
    write(*,'(A,1X,ES24.16)') "ACCELNET_ENERGY_SECONDS_PER_STRUCTURE", accel_energy_seconds
    write(*,'(A,1X,ES24.16)') "AENET_ENERGY_SECONDS_PER_STRUCTURE", ref_energy_seconds
    write(*,'(A,1X,F12.6)') "ENERGY_SPEEDUP_VS_AENET", ref_energy_seconds/accel_energy_seconds
    write(*,'(A,1X,ES24.16)') "ACCELNET_ENERGY_FORCE_SECONDS_PER_STRUCTURE", accel_force_seconds
    write(*,'(A,1X,ES24.16)') "AENET_ENERGY_FORCE_SECONDS_PER_STRUCTURE", ref_force_seconds
    write(*,'(A,1X,F12.6)') "ENERGY_FORCE_SPEEDUP_VS_AENET", ref_force_seconds/accel_force_seconds
    write(*,'(A,1X,ES24.16)') "ACCELNET_FORCE_ABS_CHECKSUM", sum(abs(accel_forces))
    write(*,'(A,1X,ES24.16)') "AENET_FORCE_ABS_CHECKSUM", sum(abs(ref_forces))

    call accel_final(status); call check_status("accelnet_final", status)
    call ref_final(status); call check_status("aenet_final", status)

contains

    subroutine make_supercell(primitive,nx,ny,nz,supercell)
        type(atomic_structure), intent(in) :: primitive
        integer, intent(in) :: nx, ny, nz
        type(atomic_structure), intent(out) :: supercell
        integer :: ix, iy, iz, source, target
        real(real64) :: translation(3)
        supercell%natoms = primitive%natoms*nx*ny*nz
        supercell%pbc = primitive%pbc
        supercell%lattice(:,1) = real(nx,real64)*primitive%lattice(:,1)
        supercell%lattice(:,2) = real(ny,real64)*primitive%lattice(:,2)
        supercell%lattice(:,3) = real(nz,real64)*primitive%lattice(:,3)
        allocate(supercell%positions(3,supercell%natoms),supercell%species(supercell%natoms))
        target = 0
        do iz = 0, nz-1
            do iy = 0, ny-1
                do ix = 0, nx-1
                    translation = real(ix,real64)*primitive%lattice(:,1) + &
                                  real(iy,real64)*primitive%lattice(:,2) + &
                                  real(iz,real64)*primitive%lattice(:,3)
                    do source = 1, primitive%natoms
                        target = target+1
                        supercell%positions(:,target) = primitive%positions(:,source)+translation
                        supercell%species(target) = primitive%species(source)
                    end do
                end do
            end do
        end do
    end subroutine make_supercell

    subroutine evaluate_accel_energy(total_energy)
        real(real64), intent(out) :: total_energy
        real(real64) :: atomic_energy
        integer :: central, local_status
        total_energy = 0.0_real64
        do central = 1, structure%natoms
            call accel_energy(structure%positions(:,central), structure%species(central), &
                size(environments(central)%species), environments(central)%coordinates, &
                environments(central)%species, atomic_energy, local_status)
            call check_status("accelnet_atomic_energy", local_status)
            total_energy = total_energy+atomic_energy
        end do
    end subroutine evaluate_accel_energy

    subroutine evaluate_ref_energy(total_energy)
        real(real64), intent(out) :: total_energy
        real(real64) :: atomic_energy
        integer :: central, local_status
        total_energy = 0.0_real64
        do central = 1, structure%natoms
            call ref_energy(structure%positions(:,central), structure%species(central), &
                size(environments(central)%species), environments(central)%coordinates, &
                environments(central)%species, atomic_energy, local_status)
            call check_status("aenet_atomic_energy", local_status)
            total_energy = total_energy+atomic_energy
        end do
    end subroutine evaluate_ref_energy

    subroutine evaluate_accel_forces(total_energy, forces)
        real(real64), intent(out) :: total_energy, forces(:,:)
        real(real64) :: atomic_energy
        integer :: central, local_status
        total_energy = 0.0_real64; forces = 0.0_real64
        do central = 1, structure%natoms
            call accel_energy_forces(structure%positions(:,central), structure%species(central), central, &
                size(environments(central)%species), environments(central)%coordinates, &
                environments(central)%species, environments(central)%indices, structure%natoms, &
                atomic_energy, forces, local_status)
            call check_status("accelnet_atomic_energy_and_forces", local_status)
            total_energy = total_energy+atomic_energy
        end do
    end subroutine evaluate_accel_forces

    subroutine evaluate_ref_forces(total_energy, forces)
        real(real64), intent(out) :: total_energy, forces(:,:)
        real(real64) :: atomic_energy
        integer :: central, local_status
        total_energy = 0.0_real64; forces = 0.0_real64
        do central = 1, structure%natoms
            call ref_energy_forces(structure%positions(:,central), structure%species(central), central, &
                size(environments(central)%species), environments(central)%coordinates, &
                environments(central)%species, environments(central)%indices, structure%natoms, &
                atomic_energy, forces, local_status)
            call check_status("aenet_atomic_energy_and_forces", local_status)
            total_energy = total_energy+atomic_energy
        end do
    end subroutine evaluate_ref_forces

    subroutine check_status(operation,value)
        character(len=*), intent(in) :: operation
        integer, intent(in) :: value
        if (value /= 0) then
            write(*,'(2A,I0)') trim(operation), " failed with status ", value
            error stop 1
        end if
    end subroutine check_status
end program atomic_api_benchmark
