module spaceparams
implicit none
type spacepars
include 'spacedecl.gen'
#ifdef USEMPI
!
! administration of the lay-out of the distributed matrices
! In the following: p is the MPI-process number: p = 0 .. numprocs-1
! The array indices start with 1, so, for example, is(2) describes
! the situation for process rank 1.
! a is the global matrix, b is a local matrix on process p.
!
! is(:) and js(:) describe the location of the submatrix in process p:
!    b(1,1) coincides with a(is(p+1),js(p+1))
! 
! lm(:) and ln(:) describe the extend of b:
!    the dimensions of b on p are (lm(p+1),ln(p+1))
!    b concides with 
!       a(is(p+1):is(p+1)+lm(p+1)-1,js(p+1):js(p+1)+ln(p+1)-1)
!
! isleft(:), isright(:), istop(:), isbot(:) tell if 
!    matrix b respectively map from a:
!    the first column
!    the last  column
!    the first row
!    lhe last  row
! 
integer, dimension(:), pointer  :: is      => NULL()   
integer, dimension(:), pointer  :: js      => NULL()
integer, dimension(:), pointer  :: lm      => NULL()
integer, dimension(:), pointer  :: ln      => NULL()
logical, dimension(:), pointer  :: isleft  => NULL()
logical, dimension(:), pointer  :: isright => NULL()
logical, dimension(:), pointer  :: istop   => NULL()
logical, dimension(:), pointer  :: isbot   => NULL()
#endif
end type                                         
                                                 
#ifdef USEMPI

interface space_distribute
  module procedure space_distribute_matrix_real8
  module procedure space_distribute_matrix_integer
  module procedure space_distribute_block_real8
  module procedure space_distribute_block_integer
  module procedure space_distribute_block4_real8
  module procedure space_distribute_vector
  module procedure space_distribute_block_vector
end interface space_distribute

interface space_shift_borders
  module procedure space_shift_borders_matrix_real8
  module procedure space_shift_borders_block_real8
end interface space_shift_borders

interface space_collect
  module procedure space_collect_block_real8
  module procedure space_collect_block_integer
  module procedure space_collect_block4_real8
  module procedure space_collect_matrix_real8
  module procedure space_collect_matrix_integer
end interface space_collect

interface compare
  module procedure comparei2
  module procedure comparer2
  module procedure comparer3
end interface compare

#endif

interface printsum
  module procedure printsum0
  module procedure printsum1
  module procedure printsum2
  module procedure printsum3
  module procedure printsum4
  module procedure printsumi0
  module procedure printsumi1
  module procedure printsumi2
  module procedure printsumi3
  module procedure printsumi4
end interface printsum


contains                                         

subroutine doeniets(sg,sl,par)
  use params
  implicit none
  type (spacepars) :: sg,sl
  type(parameters) :: par

  type(spacepars) :: gg,ll
  type(parameters) :: ppp

  gg=sg
  ll=sl
  ppp=par

end subroutine doeniets

subroutine indextos(s,index,t)
  use mnemmodule
  use xmpi_module
  use logging_module
  implicit none
  type (spacepars), intent(in)        :: s
  integer, intent(in)                 :: index
  type(arraytype), intent(out)        :: t

  if (index .lt. 1 .or. index .gt. numvars) then
    call writelog('els','ai0a','invalid index ',index,' in indextos. Program will stop')
    call halt_program
  endif

  select case(index)
  include 'indextos.gen'
  end select

end subroutine indextos

! Generated subroutine to allocate the scalars in s
subroutine space_alloc_scalars(s)
  use mnemmodule
  implicit none
  type(spacepars),intent(inout)  :: s

  include 'space_alloc_scalars.gen'

end subroutine space_alloc_scalars

! Generated subroutine to allocate all arrays in s
subroutine space_alloc_arrays(s,par)
  use mnemmodule
  use params
  implicit none
  type(spacepars),intent(inout)  :: s
  type(parameters),intent(in)    :: par

  include 'space_alloc_arrays.gen'

end subroutine space_alloc_arrays

subroutine grid_bathy(s,par)                         
          
  use params                
  use xmpi_module
  use general_mpi_module
  use readkey_module
                                                           
  implicit none                                    
                                                   
  type(spacepars)                     :: s         
  type(parameters)                    :: par         

  character*80                        :: fnameh,fnamex,fnamey
  integer                             :: i
  integer                             :: j
  integer                             :: itheta
  real*8                              :: degrad

  !                     Input file  Keyword Default  Minimum  Maximum
  s%nx    = readkey_int('params.txt','nx',     50,      2,     10000)
  s%ny    = readkey_int('params.txt','ny',      2,      2,     10000)
  s%dx    = readkey_dbl('params.txt','dx',    0.d0,   -1d9,      1d9)
  s%dy    = readkey_dbl('params.txt','dy',    0.d0,   -1d9,      1d9)
  s%xori  = readkey_dbl('params.txt','xori',  0.d0,   -1d9,      1d9)
  s%yori  = readkey_dbl('params.txt','yori',  0.d0,   -1d9,      1d9)
  s%alfa  = readkey_dbl('params.txt','alfa',  0.d0,   -360.d0,   360.d0)
  if (s%alfa.lt.0) then 
  s%alfa = 360.d0+s%alfa
  endif

  s%alfa  = s%alfa*atan(1.0d0)/45.d0
  s%posdwn= readkey_dbl('params.txt','posdwn',1.d0,   -1.d0,     1.d0)
  s%posdwn = s%posdwn*sign(s%posdwn,1.d0)
  s%vardx = readkey_int('params.txt','vardx',   0,      0,         1)     !Jaap

  if (xmaster) then
    allocate(s%x(1:s%nx+1,1:s%ny+1))
    allocate(s%y(1:s%nx+1,1:s%ny+1))
    allocate(s%xz(1:s%nx+1))
    allocate(s%yz(1:s%ny+1))
    allocate(s%xu(1:s%nx+1))
    allocate(s%yv(1:s%ny+1))
    allocate(s%xw(1:s%nx+1,1:s%ny+1))
    allocate(s%yw(1:s%nx+1,1:s%ny+1))
    allocate(s%zb(1:s%nx+1,1:s%ny+1))
    allocate(s%zb0(1:s%nx+1,1:s%ny+1))
	allocate(s%dzbdx(1:s%nx+1,1:s%ny+1))
    allocate(s%dzbdy(1:s%nx+1,1:s%ny+1))
  endif
  !
  ! Create grid
  !
  if(s%vardx==0)then

    if (xmaster) then
      call readkey('params.txt','depfile',fnameh)
      open(31,file=fnameh)
      do j=1,s%ny+1
          read(31,*)(s%zb(i,j),i=1,s%nx+1)
      end do
      close(31)
      s%zb=-s%zb*s%posdwn
	  ! Make sure that at the lateral boundaries the bathymetry is alongshore uniform
	  s%zb(:,1) = s%zb(:,2)
	  s%zb(:,s%ny+1) = s%zb(:,s%ny)
      do j=1,s%ny+1
         do i=1,s%nx+1
           s%x(i,j)=(i-1)*s%dx
           s%y(i,j)=(j-1)*s%dy
         end do
      end do
    endif

  elseif(s%vardx==1)then

    if (xmaster) then
      call readkey('params.txt','depfile',fnameh)
      call readkey('params.txt','xfile',fnamex)
      call readkey('params.txt','yfile',fnamey)

      open(31,file=fnameh)
      open(32,file=fnamex)
      open(33,file=fnamey)
      read(31,*)((s%zb(i,j),i=1,s%nx+1),j=1,s%ny+1)
      read(32,*)((s%x(i,j),i=1,s%nx+1),j=1,s%ny+1)
      read(33,*)((s%y(i,j),i=1,s%nx+1),j=1,s%ny+1)
      close(31)
      close(32)
      close(33)
      if (abs(s%x(1,2)-s%x(1,1))>par%eps) then
         ! Apparently input grid is at an angle, therefore defined in world coordinates
         ! Find out grid orientation
         s%alfa=atan2(s%y(2,1)-s%y(1,1),s%x(2,1)-s%x(1,1))

         s%xori=s%x(1,1)
         s%yori=s%y(1,1)
         s%xw=s%x
         s%yw=s%y
         s%x= cos(s%alfa)*(s%xw-s%xori)+sin(s%alfa)*(s%yw-s%yori)
         s%y=-sin(s%alfa)*(s%xw-s%xori)+cos(s%alfa)*(s%yw-s%yori)
      endif
    s%zb=-s%zb*s%posdwn
	! Make sure that at the lateral boundaries the bathymetry is alongshore uniform
	s%zb(:,1) = s%zb(:,2)
	s%zb(:,s%ny+1) = s%zb(:,s%ny)
    endif

  endif

  if(xmaster) then
    s%xz = s%x(:,1)
    s%yz = s%y(1,:)
    s%xu(1:s%nx) = 0.5d0*(s%xz(1:s%nx)+s%xz(2:s%nx+1))
    s%xu(s%nx+1) = s%xz(s%nx+1)+0.5d0*(s%xz(s%nx+1)-s%xz(s%nx))
    s%yv(1:s%ny) = 0.5d0*(s%yz(1:s%ny)+s%yz(2:s%ny+1))
    s%yv(s%ny+1) = s%yz(s%ny+1)+0.5d0*(s%yz(s%ny+1)-s%yz(s%ny))

    s%xw=s%xori+s%x*cos(s%alfa)-s%y*sin(s%alfa)
    s%yw=s%yori+s%x*sin(s%alfa)+s%y*cos(s%alfa)

    s%zb0  = s%zb
    !
    ! Specify theta-grid
    !
    !
    ! from Nautical wave directions in degrees to Cartesian wave directions in radian !!!
    !
    s%theta0=(1.5d0*par%px-s%alfa)-par%dir0*atan(1.d0)/45.d0 ! Updated in waveparams.f90 for instat 4,5,6,7
    if (s%theta0<-par%px) s%theta0=s%theta0+2.d0*par%px
    if (s%theta0> par%px) s%theta0=s%theta0-2.d0*par%px
    degrad=par%px/180.d0
    if (par%thetanaut==1) then  
       s%thetamin=(270-par%thetamax)*degrad-s%alfa
       s%thetamax=(270-par%thetamin)*degrad-s%alfa
       if (s%thetamax>par%px) then
          s%thetamax=s%thetamax-2*par%px
          s%thetamin=s%thetamin-2*par%px
       endif
       if (s%thetamin<-par%px) then
          s%thetamax=s%thetamax+2*par%px
          s%thetamin=s%thetamin+2*par%px
       endif
    else
       s%thetamin=par%thetamin*degrad
       s%thetamax=par%thetamax*degrad
    endif
    s%dtheta=par%dtheta*degrad
    s%ntheta=(s%thetamax-s%thetamin)/s%dtheta

    allocate(s%theta(1:s%ntheta))
    allocate(s%cxsth(1:s%ntheta))
    allocate(s%sxnth(1:s%ntheta))

    do itheta=1,s%ntheta
        s%theta(itheta)=s%thetamin+s%dtheta/2+s%dtheta*(itheta-1)
    end do

    s%cxsth=dcos(s%theta)
    s%sxnth=dsin(s%theta)
  endif

  if (xmaster) then
    ! Initialize dzbdx, dzbdy
    do j=1,s%ny+1
       do i=1,s%nx
          s%dzbdx(i,j)=(s%zb(i+1,j)-s%zb(i,j))/(s%xz(i+1)-s%xz(i))
       enddo
    enddo
    ! dummy, needed to keep compiler happy
    s%dzbdx(s%nx+1,:)=s%dzbdx(s%nx,:)

    do j=1,s%ny
       do i=1,s%nx+1
          s%dzbdy(i,j)=(s%zb(i,j+1)-s%zb(i,j))/(s%yz(j+1)-s%yz(j))
       enddo
    enddo
    s%dzbdy(:,s%ny+1)=s%dzbdy(:,s%ny)
  endif
end subroutine grid_bathy                         

#ifdef USEMPI

!
! consistency check: in general, all distributed matrices and blocks
! should have the following properties:
!
!  a  is distributed matrix in this process
!  al is distributed matrix in left neighbour process
!  ar is distributed matrix in right neighbour process
!  at is distributed matrix in top neighbour process
!  ab is distributed matrix in bottom neighbour process
!
!  a(:,2) = al(:,ny+1)
!  a(2,:) = at(nx+1,:)
!  except the first and last elements of these arrays
!  
!
!  When all these checks are ok on every process, automatically care 
! has been taken for the other neighbours

!
! constencycheck for mnem
! if mnem = 'ALL' then all
subroutine space_consistency(s,mnem)
use mnemmodule
implicit none
type(spacepars)  :: s
character(len=*) :: mnem
integer          :: j,jmin,jmax
type(arraytype)  :: t

if(mnem .eq. 'ALL') then
  jmin = 1
  jmax = numvars
else
  jmin = chartoindex(mnem)
  jmax = jmin
endif
do j=jmin,jmax
  call indextos(s,j,t)
  select case(t%type)
    case('r')
      select case (t%rank)
        case(0)
          !call compare(t%r0,t%name)
        case(1)
          !call compare(t%r1,t%name)
        case(2)
          call compare(t%r2,t%name)
        case(3)
          call compare(t%r3,t%name)
        case(4)
          !call compare(t%r4,t%name)
      end select ! rank
    case('i')
      select case (t%rank)
        case(0)
          !call compare(t%i0,t%name)
        case(1)
          !call compare(t%i1,t%name)
        case(2)
          call compare(t%i2,t%name)
        case(3)
          !call compare(t%i3,t%name)
        case(4)
          !call compare(t%i4,t%name)
      end select ! rank
  end select ! type
enddo

end subroutine space_consistency

subroutine comparer2(x,s)
  use xmpi_module
  use mnemmodule
  implicit none
  real*8, dimension(:,:) :: x
  character(len=*)       :: s

  real*8, parameter     :: eps=1.0d-60
  integer               :: m
  integer               :: n
  real*8, dimension(:), allocatable :: c
  real*8, dimension(:), allocatable :: r
  real*8, dimension(2)  :: dif,difmax
  character*100         :: warning

  select case(s)
    case (mnem_tideinpz)
      return
  end select
  m=size(x,1)
  n=size(x,2)
  allocate(c(m))
  allocate(r(n))
  c = x(:,1)
  call xmpi_shift(x,':1')

  dif(1) = sum(abs(c(2:m-1)-x(2:m-1,1)))
  dif(2) = sum(abs(c(1:m  )-x(1:m  ,1)))
  x(:,1) = c

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (:,1) '//trim(s)//': ',difmax,trim(warning)
  endif

  c = x(:,n)
  call xmpi_shift(x,':n')

  dif(1) = sum(abs(c(2:m-1)-x(2:m-1,n)))
  dif(2) = sum(abs(c(1:m  )-x(1:m  ,n)))
  x(:,n) = c

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (:,n) '//trim(s)//': ',difmax,trim(warning)
  endif
  r = x(1,:)
  call xmpi_shift(x,'1:')

  dif(1) = sum(abs(r(2:n-1)-x(1,2:n-1)))
  dif(2) = sum(abs(r(1:n  )-x(1,1:n  )))
  x(1,:) = r

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (1,:) '//trim(s)//': ',difmax,trim(warning)
  endif

  r = x(m,:)
  call xmpi_shift(x,'m:')

  dif(1) = sum(abs(r(2:n-1)-x(m,2:n-1)))
  dif(2) = sum(abs(r(1:n  )-x(m,1:n  )))
  x(m,:) = r

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (m,:) '//trim(s)//': ',difmax,trim(warning)
  endif
end subroutine comparer2

subroutine comparei2(x,s)
  use xmpi_module
  implicit none
  integer, dimension(:,:) :: x
  character(len=*)       :: s
  integer, parameter    :: eps=0

  integer               :: m
  integer               :: n
  integer, dimension(:), allocatable :: c
  integer, dimension(:), allocatable :: r
  integer, dimension(2) :: dif,difmax
  character*100         :: warning

  m=size(x,1)
  n=size(x,2)
  allocate(c(m))
  allocate(r(n))
  c = x(:,1)
  call xmpi_shift(x,':1')

  dif(1) = sum(abs(c(2:m-1)-x(2:m-1,1)))
  dif(2) = sum(abs(c(1:m  )-x(1:m  ,1)))
  x(:,1) = c

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (:,1) '//trim(s)//': ',difmax,trim(warning)
  endif

  c = x(:,n)
  call xmpi_shift(x,':n')

  dif(1) = sum(abs(c(2:m-1)-x(2:m-1,n)))
  dif(2) = sum(abs(c(1:m  )-x(1:m  ,n)))
  x(:,n) = c

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (:,n) '//trim(s)//': ',difmax,trim(warning)
  endif
  r = x(1,:)
  call xmpi_shift(x,'1:')

  dif(1) = sum(abs(r(2:n-1)-x(1,2:n-1)))
  dif(2) = sum(abs(r(1:n  )-x(1,1:n  )))
  x(1,:) = r

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (1,:) '//trim(s)//': ',difmax,trim(warning)
  endif

  r = x(m,:)
  call xmpi_shift(x,'m:')

  dif(1) = sum(abs(r(2:n-1)-x(m,2:n-1)))
  dif(2) = sum(abs(r(1:n  )-x(m,1:n  )))
  x(m,:) = r

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (m,:) '//trim(s)//': ',difmax,trim(warning)
  endif
end subroutine comparei2

subroutine comparer3(x,s)
use xmpi_module
  implicit none
  real*8, dimension(:,:,:) :: x
  character(len=*)       :: s

  real*8, parameter     :: eps=1.0d-60
  integer               :: m,n,l
  real*8, dimension(:,:), allocatable :: c
  real*8, dimension(:,:), allocatable :: r
  real*8, dimension(2)  :: dif,difmax
  character*100         :: warning

  m=size(x,1)
  n=size(x,2)
  l=size(x,3)

  allocate(c(m,l))
  allocate(r(n,l))
  c = x(:,1,:)
  call xmpi_shift(x,':1')

  dif(1) = sum(abs(c(2:m-1,:)-x(2:m-1,1,:)))
  dif(2) = sum(abs(c(1:m,:)  -x(1:m  ,1,:)))
  x(:,1,:) = c

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (:,1) '//trim(s)//': ',difmax,trim(warning)
  endif

  c = x(:,n,:)
  call xmpi_shift(x,':n')

  dif(1) = sum(abs(c(2:m-1,:)-x(2:m-1,n,:)))
  dif(2) = sum(abs(c(1:m  ,:)-x(1:m  ,n,:)))
  x(:,n,:) = c

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (:,n) '//trim(s)//': ',difmax,trim(warning)
  endif

  r = x(1,:,:)
  call xmpi_shift(x,'1:')

  dif(1) = sum(abs(r(2:n-1,:)-x(1,2:n-1,:)))
  dif(2) = sum(abs(r(1:n  ,:)-x(1,1:n  ,:)))
  x(1,:,:) = r

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (1,:) '//trim(s)//': ',difmax,trim(warning)
  endif

  r = x(m,:,:)
  call xmpi_shift(x,'m:')

  dif(1) = sum(abs(r(2:n-1,:)-x(m,2:n-1,:)))
  dif(2) = sum(abs(r(1:n  ,:)-x(m,1:n  ,:)))
  x(m,:,:) = r

  call xmpi_reduce(dif,difmax,MPI_SUM)

  if(xmaster) then
    warning=' '
    if (sum(difmax) .gt. eps) then
      warning = '<===++++++++++++++++++++++'
    endif
    write (*,*) 'compare (m,:) '//trim(s)//': ',difmax,trim(warning)
  endif
end subroutine comparer3

! copies scalars from sg to sl on xmaster, and distributes
! them 
subroutine space_copy_scalars(sg,sl)
  use mnemmodule
  implicit none
  type(spacepars),intent(inout)  :: sg,sl

  type(arraytype)                :: tg,tl
  integer                        :: j

  do j = 1,numvars
    call indextos(sg,j,tg)
    if (tg%rank .eq. 0) then
      call indextos(sl,j,tl)
      select case (tg%type)
        case('i')
          tl%i0 = tg%i0
        case('r')
          tl%r0 = tg%r0
      end select
    endif
  enddo
end subroutine space_copy_scalars

! The scalars are needed for allocating the arrays,
! we distribute them here:
subroutine space_distribute_scalars(sl)
  use mnemmodule
  use xmpi_module
  type (spacepars) :: sl

  type (arraytype) :: tl
  integer          :: i

  do i=1,numvars
    call indextos(sl,i,tl)
    if (tl%rank .eq. 0) then
      if (tl%type .eq. 'i') then
        call xmpi_bcast(tl%i0)
      else
        call xmpi_bcast(tl%r0)
      endif
    endif
  enddo

end subroutine space_distribute_scalars

subroutine space_distribute_matrix_real8(sl,a,b)
  use xmpi_module
  use general_mpi_module
  implicit none
  type (spacepars), intent(inout)     :: sl
  real*8, dimension(:,:), intent(in)  :: a
  real*8, dimension(:,:), intent(out) :: b

  call matrix_distr(a,b,sl%is,sl%lm,sl%js,sl%ln,xmpi_master,xmpi_comm)

end subroutine space_distribute_matrix_real8

subroutine space_distribute_matrix_integer(sl,a,b)
  use xmpi_module
  use general_mpi_module
  implicit none
  type (spacepars), intent(inout)      :: sl
  integer, dimension(:,:), intent(in)  :: a
  integer, dimension(:,:), intent(out) :: b

  call matrix_distr(a,b,sl%is,sl%lm,sl%js,sl%ln,xmpi_master,xmpi_comm)

end subroutine space_distribute_matrix_integer

subroutine space_distribute_block_real8(sl,a,b)
  use xmpi_module
  use general_mpi_module
  implicit none
  type (spacepars), intent(inout)                :: sl
  real*8, dimension(:,:,:), intent(in)           :: a
  real*8, dimension(:,:,:), intent(out)          :: b

  integer                                        :: i

  do i=1,size(b,3)   ! assuming that b is allocated on all processes
    call matrix_distr(a(:,:,i),b(:,:,i),sl%is,sl%lm,sl%js,sl%ln,xmpi_master,xmpi_comm)
  enddo

end subroutine space_distribute_block_real8

subroutine space_distribute_block_integer(sl,a,b)
  use xmpi_module
  use general_mpi_module
  implicit none
  type (spacepars), intent(inout)                :: sl
  integer, dimension(:,:,:), intent(in)          :: a
  integer, dimension(:,:,:), intent(out)         :: b

  integer                                        :: i

  do i=1,size(b,3)   ! assuming that b is allocated on all processes
    call matrix_distr(a(:,:,i),b(:,:,i),sl%is,sl%lm,sl%js,sl%ln,xmpi_master,xmpi_comm)
  enddo

end subroutine space_distribute_block_integer



subroutine space_distribute_block4_real8(sl,a,b)
  use xmpi_module
  use general_mpi_module
  implicit none
  type (spacepars), intent(inout)                :: sl
  real*8, dimension(:,:,:,:), intent(in)         :: a
  real*8, dimension(:,:,:,:), intent(out)        :: b

  integer                                        :: i,j

  do i=1,size(b,3)
    do j=1,size(b,4)
      call matrix_distr(a(:,:,i,j),b(:,:,i,j),sl%is,sl%lm,sl%js,sl%ln,xmpi_master,xmpi_comm)
    enddo
  enddo

end subroutine space_distribute_block4_real8

subroutine space_distribute_vector(xy,sl,a,b)
  use xmpi_module
  use general_mpi_module
  implicit none
  character, intent(in)             :: xy
  type (spacepars), intent(inout)   :: sl
  real*8, dimension(:), intent(in)  :: a
  real*8, dimension(:), intent(out) :: b

  integer, dimension(:), pointer    :: ijs,lmn

  if(xy .eq.'x') then
    ijs => sl%is
    lmn => sl%lm
  else
    ijs => sl%js
    lmn => sl%ln
  endif
  call vector_distr_send(a,b,ijs,lmn,xmpi_master,xmpi_comm)

end subroutine space_distribute_vector

subroutine space_distribute_block_vector(xy,sl,a,b)
  use xmpi_module
  use general_mpi_module
  implicit none
  character, intent(in)               :: xy
  type (spacepars), intent(inout)     :: sl
  real*8, dimension(:,:), intent(in)  :: a
  real*8, dimension(:,:), intent(out) :: b

  integer                             :: i

  do i=1,sl%ntheta
    call space_distribute(xy,sl,a(:,i),b(:,i))
  enddo

end subroutine space_distribute_block_vector

subroutine space_distribute_space(sg,sl,par)
  use xmpi_module
  use logging_module
  use general_mpi_module
  use params
  use mnemmodule

  implicit none
  type(spacepars), intent(inout)  :: sg
  type(spacepars), intent(inout)  :: sl 
  type(parameters)                :: par

  integer                         :: i,j,lid,eid
  real*8, pointer, dimension(:)   :: umean1, umean2
  type (arraytype)                :: tg, tl

  !
  ! This subroutine takes care that all contents of the global
  ! space is distributed to the local space.
  !

  ! allocate scalars

  call space_alloc_scalars(sl)

  ! copy scalars to sl, only on master
  ! distributing will take place later

  if(xmaster) then
    call space_copy_scalars(sg,sl)
  endif

  ! copy scalars to all processes, nx and ny will be adapted later
  call space_distribute_scalars(sl)

  ! Also, the isleft, isright, istop and isbot logicals from
  ! the xmpi module are put in sg and sl.
  !

  ! 
  ! Distribute is,js,ln,lm,isleft,isright,istop,isbot
  !

  if(xmaster) then
    allocate(sg%is(xmpi_size))
    allocate(sg%js(xmpi_size))
    allocate(sg%lm(xmpi_size))
    allocate(sg%ln(xmpi_size))
    allocate(sg%isleft(xmpi_size))
    allocate(sg%isright(xmpi_size))
    allocate(sg%istop(xmpi_size))
    allocate(sg%isbot(xmpi_size))
  endif


  allocate(sl%is(xmpi_size))
  allocate(sl%js(xmpi_size))
  allocate(sl%lm(xmpi_size))
  allocate(sl%ln(xmpi_size))
  allocate(sl%isleft(xmpi_size))
  allocate(sl%isright(xmpi_size))
  allocate(sl%istop(xmpi_size))
  allocate(sl%isbot(xmpi_size))

  if(xmaster) then
  call det_submatrices(sg%nx+1, sg%ny+1, xmpi_m, xmpi_n, &
                             sg%is, sg%lm, sg%js, sg%ln, &
                             sg%isleft, sg%isright, sg%istop, sg%isbot)
    call writelog('sl','','Distribution of matrix on processors')
    call writelog('sl','',' proc   is   lm   js   ln')
    do i=1,xmpi_size
       call writelog('sl','(i0,i0,i0,i0,i0)',i-1,sg%is(i),sg%lm(i),sg%js(i),sg%ln(i))
    enddo
    call writelog('ls','',' proc   left right top bot')
    do i=1,xmpi_size
       call writelog('ls','',i-1,sg%isleft(i),sg%isright(i),sg%istop(i),sg%isbot(i))
    enddo
  endif

  if (xmaster) then
    sl%is       = sg%is
    sl%js       = sg%js
    sl%lm       = sg%lm
    sl%ln       = sg%ln
    sl%isleft   = sg%isleft
    sl%isright  = sg%isright
    sl%istop    = sg%istop
    sl%isbot    = sg%isbot
  endif

  call xmpi_bcast(sl%is)
  call xmpi_bcast(sl%js)
  call xmpi_bcast(sl%lm)
  call xmpi_bcast(sl%ln)
  call xmpi_bcast(sl%isleft)
  call xmpi_bcast(sl%isright)
  call xmpi_bcast(sl%istop)
  call xmpi_bcast(sl%isbot)

  !
  ! compute the values for local nx and ny
  !

  sl%nx = sl%lm(xmpi_rank+1) - 1 
  sl%ny = sl%ln(xmpi_rank+1) - 1


  !
  ! allocate all arrays in sl
  !

  call space_alloc_arrays(sl,par)

  ! for each variable is sg, find out how to distribute it
  ! and distribute
  !
  do i = 1,numvars

    if(xmaster) then
      call indextos(sg,i,tg)
    endif
    call indextos(sl,i,tl)
    
    select case (tl%btype)
      case('b')           ! have to broadcast this
        select case (tl%type)
          case('i') 
            select case(tl%rank)
              case(0) ! do nothing, scalars are already in place
                ! Have to be prudent here. all scalars are broadcasted, exept nx
                ! and ny
                !if (tl%name .ne. mnem_nx .and. tl%name .ne. mnem_ny) then
                !  if(xmaster) then
                !    tl%i0 = tg%i0
                !  endif
                !  call xmpi_bcast(tl%i0)
                !endif
              case(1)
                if(xmaster) then
                  tl%i1 = tg%i1
                endif
                call xmpi_bcast(tl%i1)
              case default
                goto 100
            end select   ! rank
          case('r')
            select case(tl%rank)
              case(0) ! do nothing here, scalars are already in place
                ! Have to be prudent here. all scalars are broadcasted, exept nx
                ! and ny
                ! here only real*8 is broadcasted, so no check necessary here
                !if(xmaster) then
                !  tl%r0 = tg%r0
                !endif
                !call xmpi_bcast(tl%r0)
              case(1)
                if(xmaster) then
                  tl%r1 = tg%r1
                endif
                call xmpi_bcast(tl%r1)
              case(2)
                if(xmaster) then
                  tl%r2 = tg%r2
                endif
                call xmpi_bcast(tl%r2)
              case default
                goto 100
            end select   ! rank
          case default
            goto 100
        end select     ! type
      case('d')        ! have to distribute this
        select case(tl%type)
          case ('i')
            select case(tl%rank)
              case(2)
                call space_distribute(sl,tg%i2,tl%i2)
				  case(3)
				    call space_distribute(sl,tg%i3,tl%i3)
              case default
                goto 100
            end select  ! rank
          case ('r')
            select case(tl%rank)
              case(1)
                select case(tl%name)
                  case(mnem_xz,mnem_xu)
                    call space_distribute_vector('x',sl,tg%r1,tl%r1)
                  case(mnem_yz, mnem_yv, mnem_bi)
                    call space_distribute_vector('y',sl,tg%r1,tl%r1)
                  case default
                    goto 100
                end select
              case(2)
                call space_distribute(sl,tg%r2,tl%r2)
              case(3)
                call space_distribute(sl,tg%r3,tl%r3)
              case(4)
                call space_distribute(sl,tg%r4,tl%r4)
              case default
            end select  ! rank
          case default
            goto 100
        end select       ! type
      case('2')    ! the umean case
                   ! dimension 2,s%ny+1
        if(xmaster) then
          allocate(umean1(sg%ny+1))
        endif
        allocate(umean2(sl%ny+1))
        do j=1,2
          if (xmaster) then
            umean1 = tg%r2(j,:)
          endif
          call space_distribute("y", sl,umean1,umean2)
          tl%r2(j,:) = umean2
        enddo
        if(xmaster) then
          deallocate(umean1)
        endif
        deallocate(umean2)
      case default
        goto 100
    end select  ! btype
  enddo  ! numvars
  return

  100 continue
  call writelog('sel','',xmpi_rank,': Error in space_distribute_space, trying to distribute:')
  call get_logfileid(lid,eid)
  call printvar(tl,lid,eid)
  call halt_program
  return

end subroutine space_distribute_space

subroutine space_shift_borders_matrix_real8(a)
  use general_mpi_module
  use xmpi_module
  implicit none
  real*8, intent(inout),dimension(:,:) :: a
    call shift_borders_matrix_real8(a,xmpi_left,xmpi_right, &
             xmpi_top,xmpi_bot,xmpi_comm)
end subroutine space_shift_borders_matrix_real8

subroutine space_shift_borders_block_real8(a)
  use general_mpi_module
  use xmpi_module
  implicit none
  real*8, intent(inout),dimension(:,:,:) :: a

  integer i

  do i=1,size(a,3)
    call shift_borders_matrix_real8(a(:,:,i),xmpi_left,xmpi_right,&
              xmpi_top,xmpi_bot,xmpi_comm)
  enddo
end subroutine space_shift_borders_block_real8

!
! wwvv a subtle point with the collect subroutines: the second
! argument: the matrix wherein the submatrices are to be collected,
! does not have to be available on the non-master processes, so
! the dimensions are not defined. 
! The third argument is always defined, on master and non-master
! processes, so its dimensions (notably the 3rd in the block subroutines
! are available
! 
! parameters of the space_collect routines:
! s: spacepars: LOCAL s 
! a: output: in this matrix the submatrices are collected
! b: input:  the local submatrix

subroutine space_collect_block_real8(s,a,b)
  use general_mpi_module
  use xmpi_module
  implicit none
  type(spacepars), intent(in)            :: s
  real*8, dimension(:,:,:), intent(out)  :: a
  real*8, dimension(:,:,:), intent(in)   :: b

  integer i
  do i = 1,size(b,3)
    call matrix_coll(a(:,:,i),b(:,:,i),s%is,s%lm,s%js,s%ln, &
                   s%isleft,s%isright,s%istop,s%isbot, &
                   xmpi_master,xmpi_comm)
  enddo

end subroutine space_collect_block_real8

subroutine space_collect_block_integer(s,a,b)
  use general_mpi_module
  use xmpi_module
  implicit none
  type(spacepars), intent(in)            :: s
  integer, dimension(:,:,:), intent(out)  :: a
  integer, dimension(:,:,:), intent(in)   :: b

  integer i

  real*8, dimension(:,:,:), allocatable :: ra,rb
  integer                             :: m,n,o

  m = size(b,1)
  n = size(b,2)
  o = size(b,3)

  allocate(rb(m,n,o))

  if (xmaster) then
    m = size(a,1)
    n = size(a,2)
	 o = size(a,3)
    allocate(ra(m,n,o))
  else
    allocate(ra(1,1,1))
  endif

  rb = b

  do i = 1,o
    call matrix_coll(ra(:,:,i),rb(:,:,i),s%is,s%lm,s%js,s%ln, &
                   s%isleft,s%isright,s%istop,s%isbot, &
                   xmpi_master,xmpi_comm)
  enddo

  if (xmaster) then
    a = ra
  endif

  deallocate(ra,rb)

end subroutine space_collect_block_integer

subroutine space_collect_block4_real8(s,a,b)
  use general_mpi_module
  use xmpi_module
  implicit none
  type(spacepars), intent(in)              :: s
  real*8, dimension(:,:,:,:), intent(out)  :: a
  real*8, dimension(:,:,:,:), intent(in)   :: b

  integer i,j
  do j = 1,size(b,4)
    do i = 1,size(b,3)
      call matrix_coll(a(:,:,i,j),b(:,:,i,j),s%is,s%lm,s%js,s%ln, &
                   s%isleft,s%isright,s%istop,s%isbot, &
                   xmpi_master,xmpi_comm)
    enddo
  enddo

end subroutine space_collect_block4_real8

subroutine space_collect_matrix_real8(s,a,b)
  use general_mpi_module
  use xmpi_module
  implicit none
  type(spacepars), intent(in)          :: s
  real*8, dimension(:,:), intent(out)  :: a
  real*8, dimension(:,:), intent(in)   :: b

  call matrix_coll(a,b,s%is,s%lm,s%js,s%ln, &
                   s%isleft,s%isright,s%istop,s%isbot, &
                   xmpi_master,xmpi_comm)

end subroutine space_collect_matrix_real8

subroutine space_collect_matrix_integer(s,a,b)
  use general_mpi_module
  use xmpi_module
  implicit none
  type(spacepars), intent(in)           :: s
  integer, dimension(:,:), intent(out)  :: a
  integer, dimension(:,:), intent(in)   :: b

  ! not used often, so we convert the integers to real*8,
  ! collect and convert back. 
  ! if this routine becomes heavily used, than a special
  ! matrix_coll has to be made. (Now we understand why
  ! C++ has templates)

  real*8, dimension(:,:), allocatable :: ra,rb
  integer                             :: m,n

  m = size(b,1)
  n = size(b,2)

  allocate(rb(m,n))

  if (xmaster) then
    m = size(a,1)
    n = size(a,2)
    allocate(ra(m,n))
  else
    allocate(ra(1,1))
  endif

  rb = b

  call matrix_coll(ra,rb,s%is,s%lm,s%js,s%ln, &
                   s%isleft,s%isright,s%istop,s%isbot, &
                   xmpi_master,xmpi_comm)

  if (xmaster) then
    a = ra
  endif

  deallocate(ra,rb)

end subroutine space_collect_matrix_integer

#endif


#ifdef USEMPI
!
!  collects data from processes in master
!  using the index number of the variable to be
!  collected
!
subroutine space_collect_index(sg,sl,index)
  use xmpi_module
  use mnemmodule
  use logging_module
  type(spacepars)                 :: sg
  type(spacepars), intent(in)     :: sl
  integer, intent(in)             :: index
  integer                         :: lid,eid

  type(arraytype)                 :: tg,tl


#ifdef USEMPI
logical, dimension(numvars)         :: avail      ! .true.: this item is collected, used to
                                                  ! prevent double space_collect
                                                  ! calls for the same item
                                                  ! 
#endif

#ifdef USEMPI
  avail = .false.
#endif


#ifdef USEMPE
  call MPE_Log_event(event_coll_start,0,'cstart')
#endif

  if(avail(index)) then
    return
  endif

  call indextos(sl,index,tl)
  if(xmaster) then
    call indextos(sg,index,tg)
  endif

  select case(tl%type)
    case('i')
      select case(tl%rank)
        case(0)             ! nothing to do
        case default     ! case 1, 2, 3 and 4 are not handled
          goto 100
      end select   ! rank
    case('r')
      select case(tl%rank)
        case(0)             ! nothing to do
        case(2)
          if (tl%name .eq. mnem_umean) then
            goto 100
          endif
          call space_collect(sl,tg%r2,tl%r2)
        case(3)
          call space_collect(sl,tg%r3,tl%r3)
        case(4)
          call space_collect(sl,tg%r4,tl%r4)
        case default
      end select   ! rank
    case default
  end select   ! type

  avail(index) = .true.

#ifdef USEMPE
  call MPE_Log_event(event_coll_start,0,'cend')
#endif

  return

  100 continue
  call writelog('lse','','Problem in space_collect_index with variable ',trim(tg%name ) )
  call writelog('lse','','Don''t know how to collect that on the masternode')
  call get_logfileid(lid,eid)
  call printvar(tl,lid,eid)
  call halt_program

end subroutine space_collect_index
#endif

! printsum* for debugging only
!
subroutine printsum0(f,str,id,val)
  implicit none
  integer, intent(in) :: f
  character(*),intent(in) :: str
  integer, intent(in) :: id
  real*8, intent(in) :: val
  write(f,*) 'printsum ',id,' ',str,':',val
end subroutine printsum0

subroutine printsum1(f,str,id,val)
  implicit none
  integer, intent(in) :: f
  character(*), intent(in) :: str
  integer, intent(in) :: id
  real*8, pointer, dimension(:) :: val
  if (associated(val) ) then
    write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
  else
    write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
  endif
end subroutine printsum1

subroutine printsum2(f,str,id,val)
  implicit none
  integer, intent(in) :: f
  character(*), intent(in) :: str
  integer, intent(in) :: id
  real*8, pointer, dimension(:,:) :: val
  if (associated(val) ) then
    write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
  else
    write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
  endif

end subroutine printsum2

subroutine printsum3(f,str,id,val)
  implicit none
  integer, intent(in) :: f
  character(*), intent(in) :: str
  integer, intent(in) :: id
  real*8, pointer, dimension(:,:,:) :: val
  if (associated(val) ) then
    write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
  else
    write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
  endif
end subroutine printsum3

subroutine printsum4(f,str,id,val)
  implicit none
  integer, intent(in) :: f
  character(*), intent(in) :: str
  integer, intent(in) :: id
  real*8, pointer, dimension(:,:,:,:) :: val
  if (associated(val) ) then
    write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
  else
    write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
  endif
end subroutine printsum4

subroutine printsumi0(f,str,id,val)
  implicit none
  integer, intent(in) :: f
  character(*),intent(in) :: str
  integer, intent(in) :: id
  integer, intent(in) :: val
  write(f,*) 'printsum ',id,' ',str,':',val
end subroutine printsumi0

subroutine printsumi1(f,str,id,val)
  implicit none
  integer, intent(in) :: f
  character(*), intent(in) :: str
  integer, intent(in) :: id
  integer, pointer, dimension(:) :: val
  if (associated(val) ) then
    write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
  else
    write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
  endif
end subroutine printsumi1

subroutine printsumi2(f,str,id,val)
  implicit none
  integer, intent(in) :: f
  character(*), intent(in) :: str
  integer, intent(in) :: id
  integer, pointer, dimension(:,:) :: val
  if (associated(val) ) then
    write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
  else
    write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
  endif

end subroutine printsumi2

subroutine printsumi3(f,str,id,val)
  implicit none
  integer, intent(in) :: f
  character(*), intent(in) :: str
  integer, intent(in) :: id
  integer, pointer, dimension(:,:,:) :: val
  if (associated(val) ) then
    write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
  else
    write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
  endif
end subroutine printsumi3

subroutine printsumi4(f,str,id,val)
  implicit none
  integer, intent(in) :: f
  character(*), intent(in) :: str
  integer, intent(in) :: id
  integer, pointer, dimension(:,:,:,:) :: val
  if (associated(val) ) then
    write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
  else
    write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
  endif
end subroutine printsumi4

subroutine printssums(s,str)
#ifdef USEMPI
  use xmpi_module
#endif
  use mnemmodule
  type (spacepars), intent(in) :: s
  character(*), intent(in) :: str
  integer :: id, f,i
  type(arraytype) :: t

#ifdef USEMPI
  id = xmpi_rank
  if (id .gt. 0 ) then
    return
  endif
#else
  id=0
#endif

  f = 100+4*numvars+id+1

  write(f, *) 'printsum: ',id,'Start of printssums ',str
  do i=1,numvars
    call indextos(s,i,t)
    select case (t%rank)
      case (0)
        if (t%type .eq. 'i') then
          call printsum(f,t%name,id,t%i0)
        else
          call printsum(f,t%name,id,t%r0)
        endif
      case (1)
        if (t%type .eq. 'i') then
          call printsum(f,t%name,id,t%i1)
        else
          call printsum(f,t%name,id,t%r1)
        endif
      case (2)
        if (t%type .eq. 'i') then
          call printsum(f,t%name,id,t%i2)
        else
          call printsum(f,t%name,id,t%r2)
        endif
      case (3)
        if (t%type .eq. 'i') then
          call printsum(f,t%name,id,t%i3)
        else
          call printsum(f,t%name,id,t%r3)
        endif
      case (4)
        if (t%type .eq. 'i') then
          call printsum(f,t%name,id,t%i4)
        else
          call printsum(f,t%name,id,t%r4)
        endif
      end select
  enddo
#ifdef USEMPI
  call printsum(f,'s%is',id,s%is) 
  call printsum(f,'s%js',id,s%js) 
  call printsum(f,'s%lm',id,s%lm) 
  call printsum(f,'s%ln',id,s%ln) 
#endif
end subroutine printssums

subroutine printssumso(s)
  type (spacepars), intent(in) :: s

  write(*,*)'Start of printssumso'

  write(*,*)'s%xw',sum(s%xw)
  write(*,*)'s%yw',sum(s%yw)
  write(*,*)'s%x',sum(s%x)
  write(*,*)'s%y',sum(s%y)
  write(*,*)'s%zs',sum(s%zs)
  write(*,*)'s%u',sum(s%u)
  write(*,*)'s%v',sum(s%v)
  write(*,*)'s%ue',sum(s%ue)
  write(*,*)'s%ve',sum(s%ve)
  write(*,*)'s%H',sum(s%H)
  write(*,*)'s%urms',sum(s%urms)
  write(*,*)'s%zb',sum(s%zb)
  write(*,*)'s%hh',sum(s%hh)
  write(*,*)'s%Fx',sum(s%Fx)
  write(*,*)'s%Fy',sum(s%Fy)
  write(*,*)'s%E',sum(s%E)
  write(*,*)'s%R',sum(s%R)
  write(*,*)'s%D',sum(s%D)
end subroutine printssumso
end module spaceparams
