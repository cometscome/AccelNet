module probe
use omp_lib
contains
subroutine run(a)
real, contiguous :: a(:,:)
integer :: i
!$omp target enter data map(to:a)
!$omp target data map(alloc:a)
!$omp target teams distribute parallel do
 do i=1,size(a,2)
 a(:,i)=a(:,i)+1
 end do
!$omp end target teams distribute parallel do
!$omp end target data
!$omp target update from(a)
!$omp target exit data map(delete:a)
end subroutine
end module
program main
use probe
real,allocatable :: a(:,:)
allocate(a(3,100))
a=0
call run(a)
if(any(a/=1)) error stop 'result'
deallocate(a)
print *, 'array runtime probe passed'
end program
