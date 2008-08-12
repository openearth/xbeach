module xmpi_module
#ifdef USEMPI
use mpi
implicit none
integer                         :: xmpi_rank    ! mpi rank of this process
integer                         :: xmpi_size    ! number of mpi processes
integer                         :: xmpi_comm    ! mpi communicator to use
integer                         :: xmpi_master  ! rank of master process
integer                         :: xmpi_m       ! 1st dimension of processor
                                                ! grid
integer                         :: xmpi_n       ! 2nd dimension of processor
                                                ! grid
integer                         :: xmpi_left    ! left neighbour
integer                         :: xmpi_right   ! right neighbour
integer                         :: xmpi_bot     ! bottom neighbour
integer                         :: xmpi_top     ! top neighbour
logical                         :: xmpi_isleft  ! submatrix is at the left side
                                                ! ie: shares first column with
                                                ! global matrix
logical                         :: xmpi_isright ! submatrix is at the right side
                                                ! ie: shares last column with
                                                ! global matrix
logical                         :: xmpi_istop   ! submatrix is at the top side
                                                ! ie: shares first row with
                                                ! global matrix
logical                         :: xmpi_isbot   ! submatrix is at the bottom side
                                                ! ie: shares first row with
                                                ! global matrix
logical                         :: xmaster      ! .true. if this process reads
                                                !  and writes files
#ifdef USEMPE
integer                         :: event_output_start 
integer                         :: event_output_end 
integer                         :: event_write_start 
integer                         :: event_write_end 
integer                         :: event_coll_start
integer                         :: event_coll_end
#endif

#else
implicit none
integer, parameter              :: xmpirank     = 0
integer, parameter              :: xmpisize     = 1
logical, parameter              :: xmaster      = .true.
logical, parameter              :: xmpi_isleft  = .true.
logical, parameter              :: xmpi_isright = .true.
logical, parameter              :: xmpi_istop   = .true.
logical, parameter              :: xmpi_isbot   = .true.
#endif
#ifdef USEMPI
interface xmpi_bcast
module procedure xmpi_bcast_int, xmpi_bcast_real8, xmpi_bcast_int8
module procedure xmpi_bcast_real4
module procedure xmpi_bcast_logical
module procedure xmpi_bcast_complex16
module procedure xmpi_bcast_array_real8
module procedure xmpi_bcast_array_integer
module procedure xmpi_bcast_matrix_integer
module procedure xmpi_bcast_matrix_real8
module procedure xmpi_bcast_char
end interface xmpi_bcast

interface xmpi_allreduce
module procedure xmpi_allreduce_real8
end interface xmpi_allreduce
#endif
contains
#ifdef USEMPI

#define USE_BARRIER
#ifdef USE_BARRIER
#define BARRIER call MPI_Barrier(xmpi_comm,ierror)
#endif
subroutine xmpi_initialize
! initialize mpi environment
  implicit none
  integer ierr
  call MPI_Init(ierr)
#ifdef USEMPE
  call MPE_Log_get_solo_eventid(event_output_start)
  call MPE_Log_get_solo_eventid(event_output_end)
  call MPE_Log_get_solo_eventid(event_write_start)
  call MPE_Log_get_solo_eventid(event_write_end)
  call MPE_Log_get_solo_eventid(event_coll_start)
  call MPE_Log_get_solo_eventid(event_coll_end)
  call MPE_Describe_event(event_output_start,'out_start','red')
  call MPE_Describe_event(event_output_end,'out_end','blue')
  call MPE_Describe_event(event_write_start,'write_start','orange')
  call MPE_Describe_event(event_write_end,'write_end','green')
  call MPE_Describe_event(event_coll_start,'coll_start','white')
  call MPE_Describe_event(event_coll_end,'coll_end','white')
#endif
  xmpi_comm = MPI_COMM_WORLD
  call MPI_Comm_rank(xmpi_comm,xmpi_rank,ierr)
  call MPI_Comm_size(xmpi_comm,xmpi_size,ierr)
  xmpi_master = 0
  xmaster  = (xmpi_rank .eq. xmpi_master)
end subroutine xmpi_initialize

subroutine xmpi_finalize
! ends mpi environment, collective subroutine
  implicit none
  integer ierr
  call MPI_Finalize(ierr)
end subroutine xmpi_finalize

subroutine xmpi_abort
! to be used if program has to end, and there is no way
! to tell other processes about this fact
  implicit none
  integer ierr
  call MPI_Abort(xmpi_comm,1,ierr)
end subroutine xmpi_abort

subroutine xmpi_determine_processor_grid
implicit none

! given the number of processes (xmpi_size), compute numbers
! xmpi_n and xmpi_m such that 
!   xmpi_n * xmpi_m = xmpi_size
!   xmpi_n and xmpi_m are as equal as possible
! Also, determine left, right upper and lower neighbours
!   and set isleft etc accordingly

  xmpi_n = 0
  xmpi_m = xmpi_size
  do while (xmpi_n .lt. xmpi_m)
    xmpi_n = xmpi_n + 1
    xmpi_m = xmpi_size/xmpi_n
    do while (xmpi_n * xmpi_m .ne. xmpi_size)
      xmpi_n = xmpi_n + 1
      xmpi_m  = xmpi_size/xmpi_n
    enddo
  enddo

! The layout of the processors is as follows:
! example 12 processors, xmpi_m=3, xmpi_n=4:

! 0 3 6  9
! 1 4 7 10 
! 2 5 8 11
  
!        top
!   left  X  right
!        bot

! Left neighbour: 
  
  xmpi_left   = xmpi_rank - xmpi_m
  xmpi_isleft = .false.
  if (xmpi_left .lt. 0) then
    xmpi_left   = MPI_PROC_NULL
    xmpi_isleft = .true.
  endif
 
! Right neighbour:

  xmpi_right   = xmpi_rank + xmpi_m
  xmpi_isright = .false.
  if (xmpi_right .ge. xmpi_size) then
    xmpi_right   = MPI_PROC_NULL
    xmpi_isright = .true.
  endif

! upper neighbour:

  if (mod (xmpi_rank,xmpi_m) .eq. 0) then
     xmpi_top   = MPI_PROC_NULL
     xmpi_istop = .true.
  else
     xmpi_top = xmpi_rank - 1
     xmpi_istop = .false.
  endif

! Lower neighbour:

  if (mod (xmpi_rank+1,xmpi_m) .eq. 0) then
    xmpi_bot   = MPI_PROC_NULL
    xmpi_isbot = .true.
  else
    xmpi_bot   = xmpi_rank + 1
    xmpi_isbot = .false.
  endif

end subroutine xmpi_determine_processor_grid

subroutine xmpi_bcast_logical(x)
  implicit none
  logical, dimension(:) :: x
  integer ierror,l
  l = size(x)
  call MPI_Bcast(x, l, MPI_LOGICAL, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_logical

subroutine xmpi_bcast_array_real8(x)
  implicit none
  real*8, dimension(:) :: x
  integer ierror,l
  l = size(x)
  call MPI_Bcast(x, l, MPI_DOUBLE_PRECISION, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_array_real8

subroutine xmpi_bcast_matrix_real8(x)
  implicit none
  real*8, dimension(:,:) :: x
  integer ierror,l
  l = size(x)
  call MPI_Bcast(x, l, MPI_DOUBLE_PRECISION, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_matrix_real8

subroutine xmpi_bcast_matrix_integer(x)
  implicit none
  integer, dimension(:,:) :: x
  integer ierror,l
  l = size(x)
  call MPI_Bcast(x, l, MPI_INTEGER, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_matrix_integer

subroutine xmpi_bcast_array_integer(x)
  implicit none
  integer, dimension(:) :: x
  integer ierror,l
  l = size(x)
  call MPI_Bcast(x, l, MPI_INTEGER, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_array_integer

subroutine xmpi_bcast_real4(x)
  implicit none
  real*4 x
  integer ierror
  call MPI_Bcast(x, 1, MPI_REAL, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_real4

subroutine xmpi_bcast_real8(x)
  implicit none
  real*8 x
  integer ierror
  call MPI_Bcast(x, 1, MPI_DOUBLE_PRECISION, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_real8

subroutine xmpi_bcast_int(x)
  implicit none
  integer x
  integer ierror
  call MPI_Bcast(x, 1, MPI_INTEGER, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_int

subroutine xmpi_bcast_int8(x)
  implicit none
  integer *8 x
  integer ierror
  call MPI_Bcast(x, 1, MPI_INTEGER8, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_int8

subroutine xmpi_bcast_complex16(x)
  implicit none
  complex*16 x
  integer ierror
  call MPI_Bcast(x, 1, MPI_DOUBLE_COMPLEX, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_complex16

subroutine xmpi_bcast_char(x)
  implicit none
  character(len=*) :: x
  integer          :: ierror,l

  if (xmpi_rank .eq. xmpi_master) then
    l = len(x)
  endif
  call xmpi_bcast(l)
  call MPI_Bcast(x, l, MPI_CHARACTER, xmpi_master, xmpi_comm, ierror)
  BARRIER
end subroutine xmpi_bcast_char

subroutine xmpi_allreduce_real8(x,op)
  implicit none
  real*8 x
  integer op

  real*8 y
  integer ierror
  y = x
  if(xmpi_comm .ne. MPI_COMM_WORLD) then
  100 continue
  goto 100
  endif
  call MPI_Allreduce(y,x,1,MPI_DOUBLE_PRECISION,op,xmpi_comm,ierror)
  BARRIER
end subroutine xmpi_allreduce_real8

subroutine xmpi_barrier
  implicit none
  integer ierror
  call MPI_Barrier(xmpi_comm,ierror)
end subroutine xmpi_barrier

#endif
subroutine halt_program
#ifdef USEMPI
  call xmpi_abort
#else
  stop
#endif
  end subroutine halt_program
end module xmpi_module
