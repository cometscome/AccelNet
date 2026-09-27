module accelnet_lj
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: cutoff_value, cutoff_derivative, &
        validate_cutoff_parameters, CUTOFF_HARD
    implicit none
    private

    real(real64), parameter :: EPS_DISTANCE2 = 1.0e-24_real64

    type, public :: lj_config
        real(real64) :: radial_rc = 0.0_real64
        integer :: num_species = 0
        integer :: cutoff_type = CUTOFF_HARD
        real(real64) :: cutoff_alpha = 0.0_real64
    contains
        procedure :: num_descriptors => lj_num_descriptors
    end type lj_config

    public :: initialize_lj_config
    public :: evaluate_lj_values
    public :: evaluate_lj_values_derivatives

contains

    subroutine initialize_lj_config(config, num_species, radial_rc, cutoff_type, cutoff_alpha)
        type(lj_config), intent(out) :: config
        integer, intent(in) :: num_species
        real(real64), intent(in) :: radial_rc
        integer, intent(in), optional :: cutoff_type
        real(real64), intent(in), optional :: cutoff_alpha

        if (num_species < 1) error stop "LJ num_species must be positive"
        if (radial_rc <= 0.0_real64) error stop "LJ radial_rc must be positive"
        config%num_species = num_species
        config%radial_rc = radial_rc
        if (present(cutoff_type)) config%cutoff_type = cutoff_type
        if (present(cutoff_alpha)) config%cutoff_alpha = cutoff_alpha
        call validate_cutoff_parameters(config%cutoff_type, config%cutoff_alpha)
    end subroutine initialize_lj_config

    integer function lj_num_descriptors(self) result(n)
        class(lj_config), intent(in) :: self
        n = 2*self%num_species
    end function lj_num_descriptors

    include 'legacy_lj_evaluation.inc'

end module accelnet_lj
