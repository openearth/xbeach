module xmpi_module
#ifdef USEMPI
use mpi
implicit none
#ifndef HAVE_MPI_WTIME
real*8, external                :: MPI_Wtime
#endif
integer                         :: xmpi_rank    ! mpi rank of this process
integer                         :: xmpi_size    ! number of mpi processes
integer                         :: xmpi_comm    ! mpi communicator to use
integer                         :: xmpi_master  ! rank of master process
integer                         :: xmpi_m       ! 1st dimension of processor
                                                ! grid
integer                         :: xmpi_n       ! 2nd dimension of processor
                                                ! grid
integer                         :: xmpi_pcol    ! my column in processor grid (starting at 1)
integer                         :: xmpi_prow    ! my row    in processor grid (starting at 1)
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
logical                         :: xmpi_bckey=.true. ! true if parameter read in readkey has to be broadcast
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
integer, parameter              :: xmpi_pcol    = 1
integer, parameter              :: xmpi_prow    = 1
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
module procedure xmpi_allreduce_r0
end interface xmpi_allreduce

interface xmpi_reduce
module procedure xmpi_reduce_r0
module procedure xmpi_reduce_i0
module procedure xmpi_reduce_r1
module procedure xmpi_reduce_i1
end interface xmpi_reduce

interface xmpi_sendrecv
module procedure xmpi_sendrecv_r1
module procedure xmpi_sendrecv_r2
module procedure xmpi_sendrecv_i1
module procedure xmpi_sendrecv_i2
end interface xmpi_sendrecv

interface xmpi_shift
module procedure xmpi_shift_r2
module procedure xmpi_shift_i2
module procedure xmpi_shift_r3
module procedure xmpi_shift_i3
end interface xmpi_shift

#endif
contains
#ifdef USEMPI

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

subroutine xmpi_determine_processor_grid(m,n)
implicit none
integer, intent(in) :: m,n  ! the dimensions of the global domain

integer mm,nn, borderlength, min_borderlength

! determine the processor grid (xmpi_m ampi_n), such that 
! - xmpi_m * xmpi_n = xmpi_size
! - the total length of the internal borders is minimal

  min_borderlength = 1000000000
  do mm = 1,xmpi_size
    nn = xmpi_size/mm
    if (mm * nn .eq. xmpi_size) then
      borderlength = (mm - 1)*n + (nn -1)*m
      if (borderlength .lt. min_borderlength) then
        xmpi_m = mm
        xmpi_n = nn
        min_borderlength = borderlength
      endif
    endif
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

! my row and column (starting with 1,1)
  
  xmpi_pcol = xmpi_rank/xmpi_m
  xmpi_prow = xmpi_rank - xmpi_pcol*xmpi_m
  xmpi_pcol = xmpi_pcol+1
  xmpi_prow = xmpi_prow+1

end subroutine xmpi_determine_processor_grid

subroutine xmpi_bcast_logical(x)
  implicit none
  logical, dimension(:) :: x
  integer ierror,l
  l = size(x)
  call MPI_Bcast(x, l, MPI_LOGICAL, xmpi_master, xmpi_comm, ierror)
end subroutine xmpi_bcast_logical

subroutine xmpi_bcast_array_real8(x)
  implicit none
  real*8, dimension(:) :: x
  integer ierror,l
  l = size(x)
  call MPI_Bcast(x, l, MPI_DOUBLE_PRECISION, xmpi_master, xmpi_comm, ierror)
end subroutine xmpi_bcast_array_real8

subroutine xmpi_bcast_matrix_real8(x)
  implicit none
  real*8, dimension(:,:) :: x
  integer ierror,l
  l = size(x)
  call MPI_Bcast(x, l, MPI_DOUBLE_PRECISION, xmpi_master, xmpi_comm, ierror)
end subroutine xmpi_bcast_matrix_real8

subroutine xmpi_bcast_matrix_integer(x)
  implicit none
  integer, dimension(:,:) :: x
  integer ierror,l
  l = size(x)
  call MPI_Bcast(x, l, MPI_INTEGER, xmpi_master, xmpi_comm, ierror)
end subroutine xmpi_bcast_matrix_integer

subroutine xmpi_bcast_array_integer(x)
  implicit none
  integer, dimension(:) :: x
  integer ierror,l
  l = size(x)
  call MPI_Bcast(x, l, MPI_INTEGER, xmpi_master, xmpi_comm, ierror)
end subroutine xmpi_bcast_array_integer

subroutine xmpi_bcast_real4(x)
  implicit none
  real*4 x
  integer ierror
  call MPI_Bcast(x, 1, MPI_REAL, xmpi_master, xmpi_comm, ierror)
end subroutine xmpi_bcast_real4

subroutine xmpi_bcast_real8(x)
  implicit none
  real*8 x
  integer ierror
  call MPI_Bcast(x, 1, MPI_DOUBLE_PRECISION, xmpi_master, xmpi_comm, ierror)
end subroutine xmpi_bcast_real8

subroutine xmpi_bcast_int(x)
  implicit none
  integer x
  integer ierror
  call MPI_Bcast(x, 1, MPI_INTEGER, xmpi_master, xmpi_comm, ierror)
end subroutine xmpi_bcast_int

subroutine xmpi_bcast_int8(x)
  implicit none
  integer *8 x
  integer ierror
  call MPI_Bcast(x, 1, MPI_INTEGER8, xmpi_master, xmpi_comm, ierror)
end subroutine xmpi_bcast_int8

subroutine xmpi_bcast_complex16(x)
  implicit none
  complex*16 x
  integer ierror
  call MPI_Bcast(x, 1, MPI_DOUBLE_COMPLEX, xmpi_master, xmpi_comm, ierror)
end subroutine xmpi_bcast_complex16

subroutine xmpi_bcast_char(x)
!
! wwvv convert string to integer array, 
!  broadcast integer array and convert back
!  Not very efficient, but this routine is seldom called
!  and now it works for every taste of fortran90
!
  implicit none
  character(len=*)                   :: x

  integer                            :: l,i
  integer, allocatable, dimension(:) :: sx

  if (xmaster) then
    l = len_trim(x)
  endif

  call xmpi_bcast(l)
  allocate(sx(l))

  if (xmaster) then
    do i = 1,l
      sx(i) = ichar(x(i:i))
    enddo
  endif

  call xmpi_bcast(sx)

  if ( .not. xmaster) then
    x = ' '
    do i = 1,l
      x(i:i) = char(sx(i))
    enddo
  endif

  deallocate(sx)

end subroutine xmpi_bcast_char

subroutine xmpi_sendrecv_r1(sendbuf,n,dest,recvbuf,source)
  implicit none
  real*8, dimension(:), intent(in)  :: sendbuf
  real*8, dimension(:), intent(out) :: recvbuf
  integer, intent(in)               :: n           ! number of elements
  integer, intent(in)               :: dest,source

  integer                           :: ierror

  call MPI_Sendrecv(sendbuf,n,MPI_DOUBLE_PRECISION,dest,100,   &
                    recvbuf,n,MPI_DOUBLE_PRECISION,source,100, &
                    xmpi_comm,MPI_STATUS_IGNORE,ierror)

end subroutine xmpi_sendrecv_r1

subroutine xmpi_sendrecv_r2(sendbuf,n,dest,recvbuf,source)
  implicit none
  real*8, dimension(:,:), intent(in)  :: sendbuf
  real*8, dimension(:,:), intent(out) :: recvbuf
  integer, intent(in)                 :: n           ! number of elements
  integer, intent(in)                 :: dest,source

  integer                             :: ierror

  call MPI_Sendrecv(sendbuf,n,MPI_DOUBLE_PRECISION,dest,101,   &
                    recvbuf,n,MPI_DOUBLE_PRECISION,source,101, &
                    xmpi_comm,MPI_STATUS_IGNORE,ierror)

end subroutine xmpi_sendrecv_r2

subroutine xmpi_sendrecv_i1(sendbuf,n,dest,recvbuf,source)
  implicit none
  integer, dimension(:), intent(in)  :: sendbuf
  integer, dimension(:), intent(out) :: recvbuf
  integer, intent(in)                :: n           ! number of elements
  integer, intent(in)                :: dest,source

  integer                            :: ierror

  call MPI_Sendrecv(sendbuf,n,MPI_INTEGER,dest,102,   &
                    recvbuf,n,MPI_INTEGER,source,102, &
                    xmpi_comm,MPI_STATUS_IGNORE,ierror)

end subroutine xmpi_sendrecv_i1

subroutine xmpi_sendrecv_i2(sendbuf,n,dest,recvbuf,source)
  implicit none
  integer, dimension(:,:), intent(in)  :: sendbuf
  integer, dimension(:,:), intent(out) :: recvbuf
  integer, intent(in)                  :: n           ! number of elements
  integer, intent(in)                  :: dest,source

  integer                              :: ierror

  call MPI_Sendrecv(sendbuf,n,MPI_INTEGER,dest,103,   &
                    recvbuf,n,MPI_INTEGER,source,103, &
                    xmpi_comm,MPI_STATUS_IGNORE,ierror)

end subroutine xmpi_sendrecv_i2

subroutine xmpi_allreduce_r0(x,op)
  implicit none
  real*8,intent(inout)  :: x
  integer,intent(in)    :: op

  real*8  :: y
  integer :: ierror
  y = x
  call MPI_Allreduce(y,x,1,MPI_DOUBLE_PRECISION,op,xmpi_comm,ierror)
end subroutine xmpi_allreduce_r0

subroutine xmpi_reduce_r0(x,y,op)
  implicit none
  real*8, intent(in)   :: x
  real*8, intent(out)  :: y
  integer, intent(in)  :: op

  integer :: ierror
  call MPI_Reduce(x,y,1,MPI_DOUBLE_PRECISION,op,xmpi_master,xmpi_comm,ierror)
end subroutine xmpi_reduce_r0

subroutine xmpi_reduce_r1(x,y,op)
  implicit none
  real*8,dimension(:), intent(in)  :: x
  real*8,dimension(:), intent(out) :: y
  integer, intent(in)              :: op

  integer :: ierror
  call MPI_Reduce(x,y,size(x),MPI_DOUBLE_PRECISION,op,xmpi_master,xmpi_comm,ierror)
end subroutine xmpi_reduce_r1

subroutine xmpi_reduce_i0(x,y,op)
  implicit none
  integer, intent(in)   :: x
  integer, intent(out)  :: y
  integer, intent(in)   :: op

  integer :: ierror
  call MPI_Reduce(x,y,1,MPI_INTEGER,op,xmpi_master,xmpi_comm,ierror)
end subroutine xmpi_reduce_i0

subroutine xmpi_reduce_i1(x,y,op)
  implicit none
  integer,dimension(:),intent(in)  :: x
  integer,dimension(:),intent(out) :: y
  integer,intent(in)               :: op

  integer :: ierror
  call MPI_Reduce(x,y,size(x),MPI_INTEGER,op,xmpi_master,xmpi_comm,ierror)
end subroutine xmpi_reduce_i1

! 
! shift routines:  x(m,n) is the matrix in this process
! direction = 'u': shift up,    send to top   x(2,:) ,  receive from bot   x(m,:)
! direction = 'd': shift down,  send to bot   x(m-1,:), receive from top   x(1,:)
! direction = 'l': shift left,  send to left  x(:,2),   receive from right x(:,n)
! direction = 'r': shift right, send to right x(:,n-1), receive from left  x(:,1)
!
! also 'm:', '1:', ':n' and ':1' can be used, easier to remember:
! 'm:' :  x(m,:) will be replaced, except for a far bottom process
! '1:' :  x(1,:) will be replaced, except for a far top process
! ':n' :  x(:,n) will be replaced, except for a far right process
! ':1' :  x(:,1) will be replaced, except for a far left process

subroutine xmpi_shift_r2(x,direction)
  implicit none
  real*8,dimension(:,:),intent(inout) :: x
  character(len=*),intent(in)         :: direction

  integer :: m,n

  m = size(x,1)
  n = size(x,2)

  select case(direction)
    case('u','m:')
      call xmpi_sendrecv(x(2,:),n,xmpi_top,    x(m,:),xmpi_bot)
    case('d','1:')
      call xmpi_sendrecv(x(m-1,:),n,xmpi_bot,  x(1,:),xmpi_top)
    case('l',':n')
      call xmpi_sendrecv(x(:,2),m,xmpi_left,   x(:,n),xmpi_right)
    case('r',':1')
      call xmpi_sendrecv(x(:,n-1),m,xmpi_right,x(:,1),xmpi_left)
    case default
      if(xmaster) then
        write (*,*) 'Invalid direction parameter for xmpi_shift_r2: "'// &
                    direction//'"'
        call halt_program
      endif
  end select
  
end subroutine xmpi_shift_r2

subroutine xmpi_shift_i2(x,direction)
  implicit none
  integer, dimension (:,:), intent(inout) :: x
  character(len=*), intent(in)            :: direction

  integer :: m,n

  m = size(x,1)
  n = size(x,2)

  select case(direction)
    case('u','m:')
      call xmpi_sendrecv(x(2,:),n,xmpi_top,    x(m,:),xmpi_bot)
    case('d','1:')
      call xmpi_sendrecv(x(m-1,:),n,xmpi_bot,  x(1,:),xmpi_top)
    case('l',':n')
      call xmpi_sendrecv(x(:,2),m,xmpi_left,   x(:,n),xmpi_right)
    case('r',':1')
      call xmpi_sendrecv(x(:,n-1),m,xmpi_right,x(:,1),xmpi_left)
    case default
      if(xmaster) then
        write (*,*) 'Invalid direction parameter for xmpi_shift_r2: "'// &
                    direction//'"'
        call halt_program
      endif
  end select
  
end subroutine xmpi_shift_i2

subroutine xmpi_shift_r3(x,direction)
  implicit none
  real*8, dimension (:,:,:), intent(inout) :: x
  character(len=*), intent(in)             :: direction

  integer :: m,n,l

  m = size(x,1)
  n = size(x,2)
  l = size(x,3)

  select case(direction)
    case('u','m:')
      call xmpi_sendrecv(x(2,:,:),  l*n,xmpi_top,    x(m,:,:),xmpi_bot)
    case('d','1:')
      call xmpi_sendrecv(x(m-1,:,:),l*n,xmpi_bot,  x(1,:,:),xmpi_top)
    case('l',':n')
      call xmpi_sendrecv(x(:,2,:),  l*m,xmpi_left,   x(:,n,:),xmpi_right)
    case('r',':1')
      call xmpi_sendrecv(x(:,n-1,:),l*m,xmpi_right,x(:,1,:),xmpi_left)
    case default
      if(xmaster) then
        write (*,*) 'Invalid direction parameter for xmpi_shift_r2: "'// &
                    direction//'"'
        call halt_program
      endif
  end select
  
end subroutine xmpi_shift_r3

subroutine xmpi_shift_i3(x,direction)
  implicit none
  integer, dimension (:,:,:), intent(inout) :: x
  character(len=*), intent(in)              :: direction

  integer :: m,n,l

  m = size(x,1)
  n = size(x,2)
  l = size(x,3)

  select case(direction)
    case('u','m:')
      call xmpi_sendrecv(x(2,:,:),  l*n,xmpi_top,    x(m,:,:),xmpi_bot)
    case('d','1:')
      call xmpi_sendrecv(x(m-1,:,:),l*n,xmpi_bot,  x(1,:,:),xmpi_top)
    case('l',':n')
      call xmpi_sendrecv(x(:,2,:),  l*m,xmpi_left,   x(:,n,:),xmpi_right)
    case('r',':1')
      call xmpi_sendrecv(x(:,n-1,:),l*m,xmpi_right,x(:,1,:),xmpi_left)
    case default
      if(xmaster) then
        write (*,*) 'Invalid direction parameter for xmpi_shift_r2: "'// &
                    direction//'"'
        call halt_program
      endif
  end select
  
end subroutine xmpi_shift_i3

subroutine xmpi_barrier
  implicit none
  integer ierror
  call MPI_Barrier(xmpi_comm,ierror)
end subroutine xmpi_barrier

!
! get a row from a matrix in the same processor column
!
subroutine xmpi_getrow(a,n,l,prow,b)
  implicit none
  real*8, dimension(:,:), intent(in) :: a    ! the matrix
  integer, intent(in)                :: n    ! number of elements in the row
  character, intent(in)              :: l    ! '1': get first row
                                             ! 'm': get last row
  integer, intent(in)                :: prow ! row number of process
                                             ! to get the row from
  real*8, dimension(:), intent(out)  :: b    ! the row from process prow

  ! Note: a and l are only needed at the sending process

  integer                            :: row,ierror
  integer                            :: ll,dest,source,tag,r
  real*8, dimension(n)               :: rowdata
  integer, allocatable, dimension(:) :: requests

  source = (xmpi_pcol-1)*xmpi_m + prow -1 ! Sending process
  ll = 1

  if (source .eq. xmpi_rank) then  ! the sending process
    select case(l)
    case('1')
      ll = 1
    case('m')
      ll = size(a,1)
    case default
      write(*,*) 'Error in xmpi_getrow, l="'//l//'"'
      call halt_program
    end select

    rowdata = a(ll,:)
    allocate(requests(xmpi_m-1))
    dest = (xmpi_pcol-1)*xmpi_m    ! First receiving process
    r = 0
    do row = 1,xmpi_m
      if (dest .eq. xmpi_rank) then       ! do no send to myself
        b = rowdata                       ! but copy
      else
        r = r + 1
        tag = dest
        call MPI_Isend(rowdata(1), n, MPI_DOUBLE_PRECISION,  &
          dest, tag, xmpi_comm, requests(r), ierror)
      endif
      dest = dest + 1                     ! next receiving process
    enddo
    call MPI_Waitall(r,requests,MPI_STATUSES_IGNORE,ierror)
    deallocate(requests)
  else                                    ! receiving process
    tag = xmpi_rank
    call MPI_Recv(b, n, MPI_DOUBLE_PRECISION,  &
       source, tag, xmpi_comm, MPI_STATUS_IGNORE, ierror)
  endif

  ! wwvv if this is needed often, than a neat subroutine, using
  !  column communicators and MPI_Bcast would be appropriate 

end subroutine xmpi_getrow


#endif
subroutine halt_program
#ifdef USEMPI
  call xmpi_abort
#else
  stop
#endif
  end subroutine halt_program
end module xmpi_module
