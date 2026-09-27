program test_neighbor_images
    use iso_fortran_env, only: real64
    use accelnet_descriptors, only: atomic_structure, neighbor_data, build_neighbor_list
    implicit none
    type(atomic_structure) :: s
    type(neighbor_data) :: nb
    integer, parameter :: n=40, bound=4
    real(real64), parameter :: rc=3.1_real64
    real(real64) :: lengths(3), displacement(3), fraction(3)
    logical :: seen(n,-bound:bound,-bound:bound,-bound:bound), expected
    integer :: example,i,j,k,x,y,z,shift(3),count_expected

    s%natoms=n; s%pbc=.true.
    allocate(s%positions(3,n),s%species(n)); s%species=1
    ! Large/mixed/small linked-cell grids, periodic self images and unwrapped atoms.
    do example=1,3
        select case(example)
        case(1)
            lengths=[12.0_real64,13.0_real64,14.0_real64]
        case(2)
            lengths=[3.5_real64,7.1_real64,11.0_real64]
        case(3)
            lengths=[2.1_real64,2.4_real64,2.8_real64]
        end select
        s%lattice=0
        do i=1,3
            s%lattice(i,i)=lengths(i)
        end do
        do i=1,n
            fraction=modulo(real(i,real64)*[0.61803_real64,0.41421_real64,0.73205_real64],1.0_real64)
            fraction(1)=fraction(1)+mod(i,3)-1
            s%positions(:,i)=lengths*fraction
        end do
        call build_neighbor_list(s,rc,nb,minimum_distance=0.1_real64)
        do i=1,n
            seen=.false.
            do k=nb%offsets(i),nb%offsets(i+1)-1
                j=nb%atom_indices(k); shift=nb%image_shifts(:,k)
                if (any(abs(shift)>bound)) error stop 'image outside test bounds'
                x=shift(1); y=shift(2); z=shift(3)
                if (seen(j,x,y,z)) error stop 'duplicate neighbor image'
                seen(j,x,y,z)=.true.
                displacement=s%positions(:,j)+lengths*shift-s%positions(:,i)
                if (maxval(abs(displacement-nb%displacements(:,k)))>1e-12_real64) &
                    error stop 'incorrect neighbor displacement'
            end do
            count_expected=0
            do z=-bound,bound
                do y=-bound,bound
                    do x=-bound,bound
                        do j=1,n
                            displacement=s%positions(:,j)+lengths*[x,y,z]-s%positions(:,i)
                            expected=sum(displacement**2)<=(rc+1e-3_real64)**2
                            if (i==j .and. x==0 .and. y==0 .and. z==0) expected=.false.
                            if (expected .neqv. seen(j,x,y,z)) error stop 'neighbor image coverage'
                            if (expected) count_expected=count_expected+1
                        end do
                    end do
                end do
            end do
            if (count_expected/=nb%count_for_atom(i)) error stop 'neighbor count'
        end do
    end do
    print *, 'PASS: orthogonal periodic image coverage, uniqueness and unwrapped displacements'
end program
