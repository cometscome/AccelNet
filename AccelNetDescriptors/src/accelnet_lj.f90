module accelnet_lj
    use iso_fortran_env, only: real64
    implicit none
    private

    real(real64), parameter :: EPS_DISTANCE2 = 1.0e-24_real64

    type, public :: lj_config
        real(real64) :: radial_rc = 0.0_real64
        integer :: num_species = 0
    contains
        procedure :: num_descriptors => lj_num_descriptors
    end type lj_config

    public :: initialize_lj_config
    public :: evaluate_lj_values
    public :: evaluate_lj_values_derivatives

contains

    subroutine initialize_lj_config(config, num_species, radial_rc)
        type(lj_config), intent(out) :: config
        integer, intent(in) :: num_species
        real(real64), intent(in) :: radial_rc

        if (num_species < 1) error stop "LJ num_species must be positive"
        if (radial_rc <= 0.0_real64) error stop "LJ radial_rc must be positive"
        config%num_species = num_species
        config%radial_rc = radial_rc
    end subroutine initialize_lj_config

    integer function lj_num_descriptors(self) result(n)
        class(lj_config), intent(in) :: self
        n = 2*self%num_species
    end function lj_num_descriptors

    subroutine evaluate_lj_values(config, displacements, neighbor_species, values)
        type(lj_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: values(:)
        integer :: j, species_index, output_index
        real(real64) :: r2, inverse_r2, inverse_r6, inverse_r12, cutoff2

        if (size(displacements, 1) /= 3) error stop "LJ displacements must have shape (3,n)"
        if (size(displacements, 2) /= size(neighbor_species)) error stop "LJ neighbor arrays differ"
        if (size(values) < config%num_descriptors()) error stop "LJ values array too small"
        values(1:config%num_descriptors()) = 0.0_real64
        cutoff2 = config%radial_rc*config%radial_rc

        do j = 1, size(neighbor_species)
            species_index = neighbor_species(j)
            if (species_index < 1 .or. species_index > config%num_species) &
                error stop "LJ neighbor species out of range"
            r2 = dot_product(displacements(:, j), displacements(:, j))
            if (r2 <= cutoff2 .and. r2 > EPS_DISTANCE2) then
                inverse_r2 = 1.0_real64/r2
                inverse_r6 = inverse_r2*inverse_r2*inverse_r2
                inverse_r12 = inverse_r6*inverse_r6
                output_index = 1 + 2*(species_index - 1)
                values(output_index) = values(output_index) + inverse_r6
                values(output_index + 1) = values(output_index + 1) + inverse_r12
            end if
        end do
    end subroutine evaluate_lj_values

    subroutine evaluate_lj_values_derivatives(config, displacements, neighbor_species, values, &
                                              derivative_center, derivative_neighbors)
        type(lj_config), intent(in) :: config
        real(real64), intent(in) :: displacements(:, :)
        integer, intent(in) :: neighbor_species(:)
        real(real64), intent(out) :: values(:)
        real(real64), intent(out) :: derivative_center(:, :)
        real(real64), intent(out) :: derivative_neighbors(:, :, :)
        integer :: j, species_index, output_index
        real(real64) :: r2, inverse_r2, inverse_r6, inverse_r12, cutoff2
        real(real64) :: derivative6(3), derivative12(3)

        if (size(derivative_center, 1) < 3 .or. &
            size(derivative_center, 2) < config%num_descriptors()) &
            error stop "LJ center derivative array too small"
        if (size(derivative_neighbors, 1) < 3 .or. &
            size(derivative_neighbors, 2) < config%num_descriptors() .or. &
            size(derivative_neighbors, 3) < size(neighbor_species)) &
            error stop "LJ neighbor derivative array too small"
        values(1:config%num_descriptors()) = 0.0_real64
        derivative_center(:, 1:config%num_descriptors()) = 0.0_real64
        derivative_neighbors(:, 1:config%num_descriptors(), :) = 0.0_real64
        cutoff2 = config%radial_rc*config%radial_rc

        do j = 1, size(neighbor_species)
            species_index = neighbor_species(j)
            if (species_index < 1 .or. species_index > config%num_species) &
                error stop "LJ neighbor species out of range"
            r2 = dot_product(displacements(:, j), displacements(:, j))
            if (r2 <= cutoff2 .and. r2 > EPS_DISTANCE2) then
                inverse_r2 = 1.0_real64/r2
                inverse_r6 = inverse_r2*inverse_r2*inverse_r2
                inverse_r12 = inverse_r6*inverse_r6
                output_index = 1 + 2*(species_index - 1)
                values(output_index) = values(output_index) + inverse_r6
                values(output_index + 1) = values(output_index + 1) + inverse_r12
                derivative6 = -6.0_real64*inverse_r6*inverse_r2*displacements(:, j)
                derivative12 = -12.0_real64*inverse_r12*inverse_r2*displacements(:, j)
                derivative_neighbors(:, output_index, j) = derivative6
                derivative_neighbors(:, output_index + 1, j) = derivative12
                derivative_center(:, output_index) = derivative_center(:, output_index) - derivative6
                derivative_center(:, output_index + 1) = &
                    derivative_center(:, output_index + 1) - derivative12
            end if
        end do
    end subroutine evaluate_lj_values_derivatives

end module accelnet_lj
