program gnu_target_client
  use iso_c_binding
  implicit none
  interface
    function create(directory,device,mode,nspecies,cutoff,handle,message) result(status) &
        bind(C,name='accelnet_target_create_n2p2')
      import
      character(c_char), intent(in) :: directory(*)
      integer(c_int), value :: device,mode
      integer(c_int), intent(out) :: nspecies
      real(c_double), intent(out) :: cutoff
      type(c_ptr), intent(out) :: handle
      character(c_char), intent(out) :: message(*)
      integer(c_int) :: status
    end function
    function compute(handle,natoms,nrows,nedges,species,centers,offsets,indices,dr,e,f,w,message) &
        result(status) bind(C,name='accelnet_target_compute')
      import
      type(c_ptr), value :: handle
      integer(c_int), value :: natoms,nrows,nedges
      integer(c_int), intent(in) :: species(*),centers(*),offsets(*),indices(*)
      real(c_double), intent(in) :: dr(*)
      real(c_double), intent(out) :: e(*)
      real(c_double), intent(inout) :: f(*),w(*)
      character(c_char), intent(out) :: message(*)
      integer(c_int) :: status
    end function
    subroutine destroy(handle) bind(C,name='accelnet_target_destroy')
      import
      type(c_ptr), value :: handle
    end subroutine
  end interface
  type(c_ptr) :: handle
  character(len=4096) :: directory,reference
  character(len=16) :: backend
  character(c_char) :: message(512)
  character(len=16) :: symbol
  integer(c_int) :: status,n,device,species(4),centers(4),offsets(5),indices(12)
  integer :: i,j,k,unit,expected_n
  real(c_double) :: cutoff,expected_cutoff,x(3,4),dr(3,12),e(4),f(3,4),w(3,3),expected(25),actual(25)
  call get_command_argument(1,directory)
  call get_command_argument(2,reference)
  call get_command_argument(3,backend)
  device=0
  if (trim(backend)=='host') device=-1
  status=create(trim(directory)//c_null_char,device,0_c_int,n,cutoff,handle,message)
  if (status/=0) error stop 'create failed'
  x=reshape([0.d0,0.d0,0.d0,1.1d0,.2d0,.1d0,.3d0,1.3d0,-.2d0,1.5d0,1.1d0,.4d0],[3,4])
  centers=[1,2,3,4]; species=[1,2,1,2]; offsets=[1,4,7,10,13]; k=0
  do i=1,4
    do j=1,4
      if(i==j) cycle
      k=k+1; indices(k)=j; dr(:,k)=x(:,j)-x(:,i)
    end do
  end do
  f=0; w=0
  status=compute(handle,4_c_int,4_c_int,12_c_int,species,centers,offsets,indices,dr,e,f,w,message)
  if(status/=0) error stop 'compute failed'
  actual=[e,reshape(f,[12]),reshape(w,[9])]
  open(newunit=unit,file=trim(reference),status='old')
  read(unit,*) expected_n
  do i=1,expected_n
    read(unit,*) symbol
  end do
  read(unit,*) expected_cutoff
  do i=1,25
    read(unit,*) expected(i)
  end do
  close(unit)
  if(n/=expected_n.or.abs(cutoff-expected_cutoff)>1.d-12) error stop 'metadata mismatch'
  if(any(abs(actual-expected)>2.d-10+2.d-10*abs(expected))) error stop 'E/F/W mismatch'
  call destroy(handle)
  print *, 'Fortran C ABI client passed; max E/F/W error:',maxval(abs(actual-expected))
end program
