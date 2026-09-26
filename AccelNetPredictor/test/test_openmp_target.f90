! Opt-in environment check, not an AccelNet GPU evaluator. Fail explicitly if
! the runtime silently executes a target region on the host.
program test_openmp_target
    use iso_fortran_env, only: real64
    use omp_lib, only: omp_get_num_devices, omp_is_initial_device
    implicit none
    integer, parameter :: edges = 4096, atoms = 127
    integer :: i, target, on_host, indices(edges)
    real(real64) :: forces(atoms), expected(atoms)
    if (omp_get_num_devices() < 1) error stop 'no OpenMP target device'
    on_host = 1
    forces = 0
    expected = 0
    do i = 1, edges
        indices(i) = 1+mod(17*i,atoms)
        expected(indices(i)) = expected(indices(i)) + real(i,real64)
    end do
    !$omp target map(from:on_host)
    on_host = merge(1,0,omp_is_initial_device())
    !$omp end target
    !$omp target teams distribute parallel do map(to:indices) map(tofrom:forces) private(target)
    do i = 1, edges
        target = indices(i)
        !$omp atomic update
        forces(target) = forces(target) + real(i,real64)
    end do
    !$omp end target teams distribute parallel do
    if (on_host /= 0) error stop 'OpenMP target fell back to CPU'
    if (any(forces /= expected)) error stop 'device map/atomic accumulation mismatch'
    print *, 'OpenMP target device execution and FP64 atomic accumulation passed'
end program
