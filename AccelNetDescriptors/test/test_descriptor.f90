program test_descriptor
    use iso_fortran_env, only: real64
    use accelnet_descriptors
    implicit none

    type(descriptor_config) :: config
    real(real64) :: t(4), dt(4), values(12), vp(12), vm(12), coefficients(12)
    real(real64) :: displacements(3, 2), dc(3, 12), dn(3, 12, 2)
    real(real64) :: contracted_center(3), contracted_neighbors(3, 2)
    integer :: species(2), failures, component, basis
    real(real64), parameter :: h = 1.0e-6_real64

    failures = 0

    call assert_close("cutoff(0)", cutoff_value(0.0_real64, 5.0_real64), 1.0_real64, 1.0e-15_real64)
    call assert_close("cutoff(Rc)", cutoff_value(5.0_real64, 5.0_real64), 0.0_real64, 0.0_real64)
    call assert_close("cutoff derivative at zero", cutoff_derivative(0.0_real64, 5.0_real64), &
                      0.0_real64, 1.0e-15_real64)

    call chebyshev_values_derivatives(0.25_real64, 0.0_real64, 1.0_real64, 3, t, dt)
    call assert_close("T0", t(1), 1.0_real64, 1.0e-15_real64)
    call assert_close("T1", t(2), -0.5_real64, 1.0e-15_real64)
    call assert_close("T2", t(3), -0.5_real64, 1.0e-15_real64)
    call assert_close("T3", t(4), 1.0_real64, 1.0e-15_real64)
    call assert_close("dT1/dr", dt(2), 2.0_real64, 1.0e-15_real64)

    call initialize_config(config, 2, 3.0_real64, 2, 3.0_real64, 2, version=0)
    if (config%num_descriptors() /= 12) call fail("descriptor count")
    call assert_close("species weight Ti", config%species_weights(1), -1.0_real64, 0.0_real64)
    call assert_close("species weight O", config%species_weights(2), 1.0_real64, 0.0_real64)

    displacements(:, 1) = [1.1_real64, 0.2_real64, -0.1_real64]
    displacements(:, 2) = [-0.3_real64, 1.2_real64, 0.4_real64]
    species = [1, 2]
    call evaluate_atom_with_derivatives(config, displacements, species, values, dc, dn)
    coefficients = [(0.05_real64*real(basis, real64), basis=1, 12)]
    call contract_atom_derivatives(config, displacements, species, coefficients, &
                                   contracted_center, contracted_neighbors)
    do component = 1, 3
        call assert_close("direct center contraction", contracted_center(component), &
                          dot_product(dc(component, :), coefficients), 2.0e-13_real64)
        call assert_close("direct neighbor 1 contraction", contracted_neighbors(component, 1), &
                          dot_product(dn(component, :, 1), coefficients), 2.0e-13_real64)
        call assert_close("direct neighbor 2 contraction", contracted_neighbors(component, 2), &
                          dot_product(dn(component, :, 2), coefficients), 2.0e-13_real64)
    end do

    do basis = 1, config%num_descriptors()
        do component = 1, 3
            call assert_close("translation derivative", dc(component, basis) + &
                              sum(dn(component, basis, :)), 0.0_real64, 2.0e-13_real64)
        end do
    end do

    call test_moment_path()

    do component = 1, 3
        displacements(component, 1) = displacements(component, 1) + h
        call evaluate_atom(config, displacements, species, vp)
        displacements(component, 1) = displacements(component, 1) - 2.0_real64*h
        call evaluate_atom(config, displacements, species, vm)
        displacements(component, 1) = displacements(component, 1) + h
        do basis = 1, config%num_descriptors()
            call assert_close("finite-difference neighbor derivative", &
                dn(component, basis, 1), (vp(basis) - vm(basis))/(2.0_real64*h), 2.0e-8_real64)
        end do
    end do

    if (failures /= 0) then
        write(*, "(I0,A)") failures, " descriptor unit tests failed"
        error stop 1
    end if
    write(*, "(A)") "descriptor unit tests passed"

contains

    subroutine test_moment_path()
        integer, parameter :: nneighbors = 20, ndescriptor = 24
        type(descriptor_config) :: moment_config
        real(real64) :: r(3, nneighbors), norm, phase, radius
        real(real64) :: direct_values(ndescriptor), moment_values(ndescriptor)
        real(real64) :: direct_center(3, ndescriptor), direct_neighbors(3, ndescriptor, nneighbors)
        real(real64) :: moment_center(3), moment_neighbors(3, nneighbors)
        real(real64) :: moment_coefficients(ndescriptor), expected
        integer :: local_species(nneighbors), j, version, local_component

        do j = 1, nneighbors
            phase = 0.71_real64*real(j, real64)
            r(:, j) = [cos(phase), sin(phase), 0.13_real64*real(mod(j, 7) - 3, real64)]
            norm = sqrt(dot_product(r(:, j), r(:, j)))
            radius = 1.0_real64 + 0.11_real64*real(mod(j, 9), real64)
            r(:, j) = radius*r(:, j)/norm
            local_species(j) = 1 + mod(j, 2)
        end do
        moment_coefficients = [(0.017_real64*real(j, real64), j = 1, ndescriptor)]
        do version = 0, 10
            if (version > 1 .and. version < 10) cycle
            call initialize_config(moment_config, 2, 4.0_real64, 4, 3.5_real64, 6, &
                                   version=version, central_type_index=2)
            call evaluate_atom_with_derivatives(moment_config, r, local_species, direct_values, &
                                                direct_center, direct_neighbors)
            call evaluate_atom(moment_config, r, local_species, moment_values)
            if (maxval(abs(moment_values - direct_values)) > 2.0e-11_real64) &
                call fail("moment values differ from direct pair values")
            call contract_atom_derivatives(moment_config, r, local_species, moment_coefficients, &
                                           moment_center, moment_neighbors)
            do local_component = 1, 3
                expected = dot_product(direct_center(local_component, :), moment_coefficients)
                call assert_close("moment center contraction", moment_center(local_component), &
                                  expected, 2.0e-11_real64)
                do j = 1, nneighbors
                    expected = dot_product(direct_neighbors(local_component, :, j), moment_coefficients)
                    call assert_close("moment neighbor contraction", moment_neighbors(local_component, j), &
                                      expected, 2.0e-11_real64)
                end do
            end do
        end do
    end subroutine test_moment_path

    subroutine assert_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: actual, expected, tolerance
        if (abs(actual - expected) > tolerance) then
            failures = failures + 1
            write(*, "(A,2(1X,ES24.16))") trim(name), actual, expected
        end if
    end subroutine assert_close

    subroutine fail(name)
        character(len=*), intent(in) :: name
        failures = failures + 1
        write(*, "(A)") trim(name)
    end subroutine fail

end program test_descriptor
