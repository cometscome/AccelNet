! Runtime boundary for the optional batch kernels. ACCELNET_TARGET_SERIAL
! compiles the identical numerical source with OpenMP directives ignored and
! without OpenMP compiler/link flags or parallel/offload runtime calls.
! NVHPC may still link libnvomp as part of its standard Fortran runtime.
module accelnet_target_runtime
#ifndef ACCELNET_TARGET_SERIAL
    use omp_lib, only: omp_get_num_devices, omp_get_default_device, omp_get_initial_device, &
        omp_is_initial_device, omp_get_wtime
    implicit none
#else
    use iso_fortran_env, only: real64, int64
    implicit none
contains
    integer function omp_get_num_devices() result(n)
        n = 0
    end function
    integer function omp_get_default_device() result(n)
        n = -1
    end function
    integer function omp_get_initial_device() result(n)
        n = -1
    end function
    logical function omp_is_initial_device() result(host)
        host = .true.
    end function
    real(real64) function omp_get_wtime() result(t)
        integer(int64) :: count,rate
        call system_clock(count,rate)
        t = real(count,real64)/real(rate,real64)
    end function
#endif
end module
