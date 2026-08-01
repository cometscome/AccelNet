program benchmark_g4_derivative
    use iso_fortran_env, only: real64
    implicit none
    integer, parameter :: nm = 54, np = 4096, nr = 400
    real(real64), allocatable :: av(:), ad(:), prod(:, :), dcj(:, :), dck(:, :)
    real(real64), allocatable :: rj(:, :, :), rk(:, :, :), oj(:, :, :), ok(:, :, :)
    real(real64), allocatable :: rjs(:, :, :), rks(:, :, :), ojs(:, :, :), oks(:, :, :)
    real(real64), allocatable :: reference_j(:, :, :), reference_k(:, :, :)
    integer, allocatable :: output(:)
    real(real64) :: ta, tf, ts, ca, cf, cs
    integer :: m, p, c

    allocate(av(nm), ad(nm), prod(nm, np), dcj(3, np), dck(3, np), &
             rj(3, nm, np), rk(3, nm, np), oj(3, nm, np), ok(3, nm, np), &
             rjs(nm, 3, np), rks(nm, 3, np), ojs(nm, 3, np), oks(nm, 3, np), &
             reference_j(3, nm, np), reference_k(3, nm, np), output(nm))
    do m = 1, nm
        av(m) = 0.1_real64 + 0.003_real64*m
        ad(m) = -0.2_real64 + 0.002_real64*m
        output(m) = m
    end do
    do p = 1, np
        dcj(:, p) = [sin(0.013_real64*p), cos(0.011_real64*p), sin(0.007_real64*p)]
        dck(:, p) = [cos(0.017_real64*p), sin(0.019_real64*p), cos(0.023_real64*p)]
        do m = 1, nm
            prod(m, p) = 0.5_real64 + 1.0e-5_real64*(m + p)
            do c = 1, 3
                rj(c, m, p) = sin(0.001_real64*(c + 3*m + p))
                rk(c, m, p) = cos(0.0013_real64*(2*c + m + p))
                rjs(m, c, p) = rj(c, m, p)
                rks(m, c, p) = rk(c, m, p)
            end do
        end do
    end do

    call array_kernel(ta, ca)
    reference_j = oj; reference_k = ok
    call fused_kernel(tf, cf)
    call check_result("fused", oj, ok)
    call soa_kernel(ts, cs)
    do p = 1, np
        do m = 1, nm
            oj(:, m, p) = ojs(m, :, p)
            ok(:, m, p) = oks(m, :, p)
        end do
    end do
    call check_result("SoA", oj, ok)

    write(*, "(A,I0,A,I0,A,I0)") "members=", nm, " pairs=", np, " repeats=", nr
    call report("array expression", ta, ca, ta)
    call report("fused contiguous", tf, cf, ta)
    call report("member-first SoA", ts, cs, ta)

contains
    subroutine array_kernel(seconds, checksum)
        real(real64), intent(out) :: seconds, checksum
        real(real64) :: t0, t1, scale
        integer :: rep, pair, member, descriptor
        checksum = 0.0_real64
        call cpu_time(t0)
        do rep = 1, nr
            scale = 2.0_real64 + 1.0e-12_real64*rep
            do pair = 1, np
                do member = 1, nm
                    descriptor = output(member)
                    oj(:, descriptor, pair) = scale*(ad(member)*dcj(:, pair)*prod(member, pair) + &
                                                      av(member)*rj(:, member, pair))
                    ok(:, descriptor, pair) = scale*(ad(member)*dck(:, pair)*prod(member, pair) + &
                                                      av(member)*rk(:, member, pair))
                end do
            end do
            checksum = checksum + oj(1, 1, modulo(rep - 1, np) + 1)
        end do
        call cpu_time(t1); seconds = t1 - t0
    end subroutine array_kernel

    subroutine fused_kernel(seconds, checksum)
        real(real64), intent(out) :: seconds, checksum
        real(real64) :: t0, t1, scale, ac, rc
        integer :: rep, pair, member, descriptor, component
        checksum = 0.0_real64
        call cpu_time(t0)
        do rep = 1, nr
            scale = 2.0_real64 + 1.0e-12_real64*rep
            do pair = 1, np
                do member = 1, nm
                    descriptor = member
                    ac = scale*ad(member)*prod(member, pair)
                    rc = scale*av(member)
                    do component = 1, 3
                        oj(component, descriptor, pair) = ac*dcj(component, pair) + &
                            rc*rj(component, member, pair)
                        ok(component, descriptor, pair) = ac*dck(component, pair) + &
                            rc*rk(component, member, pair)
                    end do
                end do
            end do
            checksum = checksum + oj(1, 1, modulo(rep - 1, np) + 1)
        end do
        call cpu_time(t1); seconds = t1 - t0
    end subroutine fused_kernel

    subroutine soa_kernel(seconds, checksum)
        real(real64), intent(out) :: seconds, checksum
        real(real64) :: t0, t1, scale, ac, rc
        integer :: rep, pair, member, component
        checksum = 0.0_real64
        call cpu_time(t0)
        do rep = 1, nr
            scale = 2.0_real64 + 1.0e-12_real64*rep
            do pair = 1, np
                do component = 1, 3
                    !GCC$ ivdep
                    do concurrent (member = 1:nm)
                        ojs(member, component, pair) = scale*ad(member)*prod(member, pair)* &
                            dcj(component, pair) + scale*av(member)*rjs(member, component, pair)
                        oks(member, component, pair) = scale*ad(member)*prod(member, pair)* &
                            dck(component, pair) + scale*av(member)*rks(member, component, pair)
                    end do
                end do
            end do
            checksum = checksum + ojs(1, 1, modulo(rep - 1, np) + 1)
        end do
        call cpu_time(t1); seconds = t1 - t0
    end subroutine soa_kernel

    subroutine check_result(name, actual_j, actual_k)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: actual_j(:, :, :), actual_k(:, :, :)
        real(real64) :: error
        error = max(maxval(abs(reference_j - actual_j)), maxval(abs(reference_k - actual_k)))
        if (error > 2.0e-15_real64) then
            write(*, "(A,1X,A,1X,ES12.4)") "kernel mismatch", trim(name), error
            error stop "G4 derivative microbenchmark mismatch"
        end if
    end subroutine check_result

    subroutine report(name, seconds, checksum, baseline)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: seconds, checksum, baseline
        real(real64) :: ns
        ns = 1.0e9_real64*seconds/real(nr*np*nm, real64)
        write(*, "(A22,2X,F8.4,A,2X,F7.3,A,2X,F6.3,A,2X,ES12.4)") trim(name), seconds, " s", &
            ns, " ns/member", seconds/baseline, "x", checksum
    end subroutine report
end program benchmark_g4_derivative
