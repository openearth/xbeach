module spaceparamsdef
   use mnemmodule
  implicit none
   save
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
     ! the values are determined in subroutine space_distribute_space
     !
     integer, dimension(:), pointer  :: is      => NULL()   
     integer, dimension(:), pointer  :: js      => NULL()
     integer, dimension(:), pointer  :: lm      => NULL()
     integer, dimension(:), pointer  :: ln      => NULL()
     logical, dimension(:), pointer  :: isleft  => NULL()
     logical, dimension(:), pointer  :: isright => NULL()
     logical, dimension(:), pointer  :: istop   => NULL()
     logical, dimension(:), pointer  :: isbot   => NULL()

      ! The following are determined in spaceparams, look there for a description
      integer, dimension(:), allocatable :: icgs(:)
      integer, dimension(:), allocatable :: icge(:)
      integer, dimension(:), allocatable :: jcgs(:)
      integer, dimension(:), allocatable :: jcge(:)
      integer, dimension(:), allocatable :: icls(:)
      integer, dimension(:), allocatable :: icle(:)
      integer, dimension(:), allocatable :: jcls(:)
      integer, dimension(:), allocatable :: jcle(:)

      logical, dimension(numvars) :: collected
      logical, dimension(numvars) :: precollected

#endif
  end type spacepars

  include 'ranges.inc'
  end module spaceparamsdef
