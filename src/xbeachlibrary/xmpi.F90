module xmpi_module

   !#define SHIFT_TIMER

#ifdef USEMPI
   use mpi
   implicit none
   save
#ifndef HAVE_MPI_WTIME
   real*8, external                :: MPI_Wtime
#endif
   logical :: bcast_backtrace = .false.
   integer, parameter              :: MPIBOUNDARY_Y    = 1
   integer, parameter              :: MPIBOUNDARY_X    = 2
   integer, parameter              :: MPIBOUNDARY_AUTO = 3
   integer, parameter              :: MPIBOUNDARY_MAN  = 4
   ! We use two communicators:
   !  xmpi_ocomm and xmpi_comm
   !  xmpi_comm will contain all computing processes
   !  xmpi_ocomm will contain all computing processes + one process
   !     that will do most of the output, hoping that the output can be
   !     done in parallel with computing
   integer                         :: xmpi_orank   ! mpi rank of this process   in xmpi_ocomm
   integer                         :: xmpi_osize   ! number of mpi processes    in xmpi_ocomm
   integer                         :: xmpi_ocomm   ! mpi communicator to use  for the I/O
   integer                         :: xmpi_omaster ! rank of master process = rank of io-process in xmpi_ocomm
   integer                         :: xmpi_imaster ! rank of master process of computational processes in xmpi_ocomm
   !
   integer                         :: xmpi_rank    ! mpi rank of this process   in xmpi_comm
   integer                         :: xmpi_size    ! number of mpi processes    in xmpi_comm
   integer                         :: xmpi_comm    ! mpi communicator to use    in xmpi_comm
   integer                         :: xmpi_master  ! rank of master process     in xmpi_comm
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
   logical                         :: xmaster      ! .true. if process is the master process
   !                                               ! of the computing processes
   logical                         :: xomaster     ! .true if this process is the output process
   !                                               !  in xmpi_ocomm
   logical                         :: xcompute     ! .true. if this is a compute process
   !                                               !
   !
   !         1 2 3 4 5 6 7   y-axis
   !  X   1  x x x x x x x
   !  a   2  x x x x x x x
   !  x   3  x x x x x x x
   !  i   4  x x x x x x x
   !  s   5  x x x x x x x
   !
   !
   integer, parameter              :: SHIFT_X_U = 1  ! shift in x direction from high to low: from bottom to top
   integer, parameter              :: SHIFT_X_D = 2  ! shift in x direction from low to high: from top to bottom
   integer, parameter              :: SHIFT_Y_R = 3  ! shift in y direction from low to high: from left to right
   integer, parameter              :: SHIFT_Y_L = 4  ! shift in y direction from high to low: from right to left
   integer, parameter              :: SHIFT_R   = 5  ! shift in 1d array from left to right
   integer, parameter              :: SHIFT_L   = 6  ! shift in 1d array from right to left
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
   save
   integer, parameter              :: xmpi_rank     = 0
   integer, parameter              :: xmpi_orank   = 0
   integer, parameter              :: xmpi_size     = 1
   logical, parameter              :: xmaster      = .true.
   logical, parameter              :: xomaster     = .true.
   logical, parameter              :: xcompute     = .true.
   logical, parameter              :: xmpi_isleft  = .true.
   logical, parameter              :: xmpi_isright = .true.
   logical, parameter              :: xmpi_istop   = .true.
   logical, parameter              :: xmpi_isbot   = .true.
   integer, parameter              :: xmpi_pcol    = 1
   integer, parameter              :: xmpi_prow    = 1
#endif

#ifdef USEMPI
   interface xmpi_bcast
      module procedure xmpi_bcast_array_integer
      module procedure xmpi_bcast_array_integer_3
      module procedure xmpi_bcast_array_logical
      module procedure xmpi_bcast_array_logical_3
      module procedure xmpi_bcast_array_real8
      module procedure xmpi_bcast_array_real8_3
      module procedure xmpi_bcast_char
      module procedure xmpi_bcast_char_3
      module procedure xmpi_bcast_complex16
      module procedure xmpi_bcast_complex16_3
      module procedure xmpi_bcast_integer
      module procedure xmpi_bcast_integer_3
      module procedure xmpi_bcast_integer8
      module procedure xmpi_bcast_integer8_3
      module procedure xmpi_bcast_logical
      module procedure xmpi_bcast_logical_3
      module procedure xmpi_bcast_matrix_integer
      module procedure xmpi_bcast_matrix_integer_3
      module procedure xmpi_bcast_matrix_real8
      module procedure xmpi_bcast_matrix_real8_3
      module procedure xmpi_bcast_real4
      module procedure xmpi_bcast_real4_3
      module procedure xmpi_bcast_real8
      module procedure xmpi_bcast_real8_3
   end interface xmpi_bcast

   interface xmpi_allreduce
      module procedure xmpi_allreduce_r0
      module procedure xmpi_allreduce_r1
      module procedure xmpi_allreduce_i0
   end interface xmpi_allreduce

   interface xmpi_reduce
      module procedure xmpi_reduce_r0
      module procedure xmpi_reduce_i0
      module procedure xmpi_reduce_r1
      module procedure xmpi_reduce_r2
      module procedure xmpi_reduce_i1
   end interface xmpi_reduce

   interface xmpi_sendrecv
      module procedure xmpi_sendrecv_r1
      module procedure xmpi_sendrecv_r2
      module procedure xmpi_sendrecv_r3
      module procedure xmpi_sendrecv_i1
      module procedure xmpi_sendrecv_i2
   end interface xmpi_sendrecv

   interface xmpi_shift
      module procedure xmpi_shift_r2
      module procedure xmpi_shift_i2
      module procedure xmpi_shift_r3
      module procedure xmpi_shift_i3

      module procedure xmpi_shift_r2_l
      module procedure xmpi_shift_r3_l
   end interface xmpi_shift

   interface xmpi_shift_ee
      module procedure xmpi_shift_ee_r2
      module procedure xmpi_shift_ee_r3
   end interface xmpi_shift_ee

   interface xmpi_shift_uu
      module procedure xmpi_shift_uu_r2
      module procedure xmpi_shift_uu_r3
   end interface xmpi_shift_uu

   interface xmpi_shift_vv
      module procedure xmpi_shift_vv_r2
      module procedure xmpi_shift_vv_r3
   end interface xmpi_shift_vv

   interface xmpi_shift_zs
      module procedure xmpi_shift_zs_r2
      module procedure xmpi_shift_zs_r3
   end interface xmpi_shift_zs

   interface xmpi_send
      module procedure xmpi_send_r0
      module procedure xmpi_send_i0
      module procedure xmpi_send_l0
      module procedure xmpi_send_r1
      module procedure xmpi_send_r2
      module procedure xmpi_send_r3
      module procedure xmpi_send_r4
      module procedure xmpi_send_i1
      module procedure xmpi_send_l1
   end interface xmpi_send

#endif
contains
#ifdef USEMPI

   subroutine xmpi_initialize
      ! initialize mpi environment
      implicit none
      integer ierr,color,errhandler,r,comm_world
      external comm_errhandler
      ierr = 0
      ! Message buffers in openmpi are not initialized so this call can give a vallgrind error
      ! http://www.open-mpi.org/community/lists/users/2009/06/9566.php
      ! http://valgrind.org/docs/manual/manual-core.html#manual-core.suppress
      call MPI_Init(ierr)
      ! Use comm_world as a variable for compatibility with SGI MPI
      comm_world = MPI_COMM_WORLD
      call MPI_Comm_create_errhandler(comm_errhandler,errhandler,ierr)
      call MPI_Comm_set_errhandler(comm_world,errhandler,ierr)
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
      xmpi_ocomm = MPI_COMM_WORLD
      call MPI_Comm_rank(xmpi_ocomm,xmpi_orank,ierr)
      call MPI_Comm_size(xmpi_ocomm,xmpi_osize,ierr)
      if (xmpi_osize < 2) then
         print *,'Number of MPI processes must be 2 or greater, but is:',xmpi_osize
         call halt_program(.false.)
      endif
      xmpi_omaster = 0
      xomaster     = (xmpi_orank == xmpi_omaster)
      xcompute     = .not. xomaster
      !
      ! Create the compute communicator. This will contain all
      ! processes in xmpi_iocomm except process xmpi_master
      !
      color = 1
      if(xmpi_orank == xmpi_omaster) color = 0
      call MPI_Comm_split(xmpi_ocomm,color,1,xmpi_comm,ierr)
      call MPI_Comm_rank(xmpi_comm,xmpi_rank,ierr)
      call MPI_Comm_size(xmpi_comm,xmpi_size,ierr)
      call MPI_Comm_set_name(xmpi_comm,'xmpi_comm',ierr)
      call MPI_Comm_set_name(xmpi_ocomm,'xmpi_ocomm',ierr)
      xmpi_master = 0
      xmaster     = (xmpi_rank == xmpi_master)
      !
      ! on the output process, xmpi_comm and xmpi_rank
      ! are of no use, so give them values that will trigger
      ! errors when used:
      !
      if (xomaster) then
         xmpi_comm = MPI_COMM_NULL
         xmpi_rank = -123
         xmaster   = .false.
      endif
      !
      ! Let the I/O process know how many computational processes there are:
      !

      if(xomaster) then
         xmpi_size = xmpi_osize - 1
      endif

      xmpi_imaster = 1

      ! sanity check for xmpi_orank_to_rank and xmpi_rank_to_orank

      r = xmpi_orank_to_rank(xmpi_orank)
      if (r .ne. xmpi_rank) then
         print *,'Wrong conversion from xmpi_orank',xmpi_orank,'to xmpi_rank',r
         call halt_program
      endif

      r = xmpi_rank_to_orank(xmpi_rank)
      if (r .ne. xmpi_orank) then
         print *,'Wrong conversion from xmpi_rank',xmpi_rank,'to xmpi_orank',r
         call halt_program
      endif

   end subroutine xmpi_initialize

   integer function xmpi_orank_to_rank(r)
      ! given rank r in xmpi_ocomm, return rank in xmpi_comm
      integer, intent(in) :: r
      integer rr
      rr = r-1
      if (rr .lt. 0) rr = -123
      xmpi_orank_to_rank = rr
   end function xmpi_orank_to_rank

   integer function xmpi_rank_to_orank(r)
      ! given rank r in xmpi_comm, return rank in xmpi_ocomm
      integer, intent(in) :: r
      if (r .eq. -123) then
         xmpi_rank_to_orank = 0
      else
         xmpi_rank_to_orank = r+1
      endif
   end function xmpi_rank_to_orank

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
      stop 1
   end subroutine xmpi_abort

   !____________________________________________________________________________

   subroutine xmpi_determine_processor_grid(m,n,divtype,mmanual,nmanual,cyclic,error)
      implicit none
      integer, intent(in) :: m,n,mmanual,nmanual  ! the dimensions of the global domain
      integer, intent(in)     :: cyclic
      integer, intent(in)     :: divtype
      integer, intent(out)    :: error

      integer mm,nn, borderlength, min_borderlength

      ! determine the processor grid (xmpi_m ampi_n), such that
      ! - xmpi_m * xmpi_n = xmpi_size
      ! - the total length of the internal borders is minimal

      min_borderlength = 1000000000
      error = 0
      select case(divtype)
       case(MPIBOUNDARY_Y)   ! Force all subdivisions to run along y-lines
         xmpi_m=xmpi_size
         xmpi_n=1
       case(MPIBOUNDARY_X)  ! Force all subdivisions to run along x-lines
         xmpi_m=1
         xmpi_n=xmpi_size
       case (MPIBOUNDARY_AUTO)
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
       case (MPIBOUNDARY_MAN)
         if (mmanual * nmanual .ne. xmpi_size) then
            error = 2
            return
         endif
         xmpi_m = mmanual
         xmpi_n = nmanual
       case default
         error = 1
         return
      end select

       ! Check whether each MPI subdomain is large enough (at least 4 cells in x and y), preferably more
       ! Errors
       if(m+1<4*xmpi_m) then
          error = 3
          return
       endif
       if (n+1<4*xmpi_n .and. n>2) then
          error = 4
          return
       endif
       ! "Warning"
       if (m+1<8*xmpi_m .or. (n+1<8*xmpi_n .and. n>2)) then
          if (m+1<8*xmpi_m) error = 5
          if (n+1<8*xmpi_n .and. n>2) error = 6
       endif


      ! Pieter: TODO: Check whether the grid distribution did not lead to domains with to few grid cells (mpi implementation needs at least 2 at each boundary). Specify error and give a message if this is the case.

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
         select case (cyclic)
          case(0)
            xmpi_left   = MPI_PROC_NULL
            xmpi_isleft = .true.
          case(1)
            xmpi_left   = xmpi_left + xmpi_m*xmpi_n
         end select
      endif

      ! Right neighbour:

      xmpi_right   = xmpi_rank + xmpi_m
      xmpi_isright = .false.
      if (xmpi_right .ge. xmpi_size) then
         select case (cyclic)
          case(0)
            xmpi_right   = MPI_PROC_NULL
            xmpi_isright = .true.
          case(1)
            xmpi_right   = xmpi_right - xmpi_m*xmpi_n
         end select
      endif

      ! Upper neighbour:

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
   !____________________________________________________________________________

   subroutine xmpi_bcast_array_logical(x,toall)
      implicit none
      logical, dimension(:)         :: x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_array_logical

   subroutine xmpi_bcast_array_logical_3(x,src,comm)
      implicit none
      logical, dimension(:) :: x
      integer, intent(in)           :: src
      integer, intent(in)           :: comm
      integer ierror,l
      l = size(x)
      call MPI_Bcast(x, l, MPI_LOGICAL, src, comm, ierror)
   end subroutine xmpi_bcast_array_logical_3

   subroutine xmpi_bcast_logical(x,toall)
      implicit none
      logical x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_logical

   subroutine xmpi_bcast_logical_3(x,src,comm)
      implicit none
      logical x
      integer, intent(in) :: src
      integer, intent(in) :: comm
      integer ierror
      call MPI_Bcast(x, 1, MPI_LOGICAL, src, comm, ierror)
   end subroutine xmpi_bcast_logical_3

   subroutine xmpi_bcast_array_real8(x,toall)
      implicit none
      real*8, dimension(:) :: x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_array_real8

   subroutine xmpi_bcast_array_real8_3(x,src,comm)
      implicit none
      real*8, dimension(:):: x
      integer, intent(in) :: src,comm
      integer ierror,l
      l = size(x)
      call MPI_Bcast(x, l, MPI_DOUBLE_PRECISION, src, comm, ierror)
   end subroutine xmpi_bcast_array_real8_3

   subroutine xmpi_bcast_matrix_real8(x,toall)
      implicit none
      real*8, dimension(:,:) :: x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_matrix_real8

   subroutine xmpi_bcast_matrix_real8_3(x,src,comm)
      implicit none
      real*8, dimension(:,:)  :: x
      integer, intent(in)     :: src,comm
      integer ierror,l
      l = size(x)
      call MPI_Bcast(x, l, MPI_DOUBLE_PRECISION, src, comm, ierror)
   end subroutine xmpi_bcast_matrix_real8_3

   subroutine xmpi_bcast_matrix_integer(x,toall)
      implicit none
      integer, dimension(:,:)   :: x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_matrix_integer

   subroutine xmpi_bcast_matrix_integer_3(x,src,comm)
      implicit none
      integer, dimension(:,:) :: x
      integer, intent(in)       :: src,comm
      integer ierror,l
      l = size(x)
      call MPI_Bcast(x, l, MPI_INTEGER, src, comm, ierror)
   end subroutine xmpi_bcast_matrix_integer_3

   subroutine xmpi_bcast_array_integer(x,toall)
      implicit none
      integer, dimension(:)   :: x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_array_integer

   subroutine xmpi_bcast_array_integer_3(x,src,comm)
      implicit none
      integer, dimension(:) :: x
      integer, intent(in)     :: src,comm
      integer ierror,l
      l = size(x)
      call MPI_Bcast(x, l, MPI_INTEGER, src, comm, ierror)
   end subroutine xmpi_bcast_array_integer_3

   subroutine xmpi_bcast_real4(x,toall)
      implicit none
      real*4              :: x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_real4

   subroutine xmpi_bcast_real4_3(x,src,comm)
      implicit none
      real*4              :: x
      integer, intent(in) :: src,comm
      integer ierror
      call MPI_Bcast(x, 1, MPI_REAL, src, comm, ierror)
   end subroutine xmpi_bcast_real4_3

   subroutine xmpi_bcast_real8(x,toall)
      implicit none
      real*8              :: x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_real8

   subroutine xmpi_bcast_real8_3(x,src,comm)
      implicit none
      real*8              :: x
      integer, intent(in) :: src,comm
      integer ierror
      call MPI_Bcast(x, 1, MPI_DOUBLE_PRECISION, src, comm, ierror)
   end subroutine xmpi_bcast_real8_3

   subroutine xmpi_bcast_integer(x,toall)
      implicit none
      integer x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_integer

   subroutine xmpi_bcast_integer_3(x,src,comm)
      implicit none
      integer             :: x
      integer, intent(in) :: src,comm
      integer ierror
      call MPI_Bcast(x, 1, MPI_INTEGER, src, comm, ierror)
   end subroutine xmpi_bcast_integer_3

   subroutine xmpi_bcast_integer8(x,toall)
      implicit none
      integer*8                     :: x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_integer8

   subroutine xmpi_bcast_integer8_3(x,src,comm)
      implicit none
      integer*8           :: x
      integer, intent(in) :: src,comm
      integer ierror
      call MPI_Bcast(x, 1, MPI_INTEGER8, src, comm, ierror)
   end subroutine xmpi_bcast_integer8_3

   subroutine xmpi_bcast_complex16(x,toall)
      implicit none
      complex*16          :: x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_complex16

   subroutine xmpi_bcast_complex16_3(x,src,comm)
      implicit none
      complex*16          :: x
      integer, intent(in) :: src,comm
      integer ierror
      call MPI_Bcast(x, 1, MPI_DOUBLE_COMPLEX, src, comm, ierror)
   end subroutine xmpi_bcast_complex16_3

   subroutine xmpi_bcast_char(x,toall)
      implicit none
      character(len=*)         :: x
      include 'xmpi_bcast.inc'
   end subroutine xmpi_bcast_char

   subroutine xmpi_bcast_char_3(x,src,comm)
      !
      ! wwvv convert string to integer array,
      !  broadcast integer array and convert back
      !  Not very efficient, but this routine is seldom called
      !  and now it works for every taste of fortran90
      !
      implicit none
      character(len=*)                   :: x
      integer, intent(in)      :: src,comm

      integer                  :: l,i,rank,ierr
      integer, allocatable, dimension(:) :: sx

      call MPI_Comm_rank(comm,rank,ierr)
      if (rank == src) then
         l = len_trim(x)
      endif

      call xmpi_bcast(l,src,comm)
      allocate(sx(l))

      if (rank == src) then
         do i = 1,l
            sx(i) = ichar(x(i:i))
         enddo
      endif

      call xmpi_bcast(sx,src,comm)

      if ( rank /= src) then
         x = ' '
         do i = 1,l
            x(i:i) = char(sx(i))
         enddo
      endif

      deallocate(sx)

   end subroutine xmpi_bcast_char_3


   subroutine xmpi_sendrecv_r1(sendbuf,dest,recvbuf,source)
      implicit none
      real*8, dimension(:), intent(in)  :: sendbuf
      real*8, dimension(:), intent(out) :: recvbuf
      integer, intent(in)               :: dest,source

      integer                           :: ierror
      integer                           :: n

      n = size(sendbuf)

      call MPI_Sendrecv(sendbuf,n,MPI_DOUBLE_PRECISION,dest,100,   &
      recvbuf,n,MPI_DOUBLE_PRECISION,source,100, &
      xmpi_comm,MPI_STATUS_IGNORE,ierror)

   end subroutine xmpi_sendrecv_r1

   subroutine xmpi_sendrecv_r2(sendbuf,dest,recvbuf,source)
      implicit none
      real*8, dimension(:,:), intent(in)  :: sendbuf
      real*8, dimension(:,:), intent(out) :: recvbuf
      integer, intent(in)                 :: dest,source

      integer                             :: ierror
      integer                             :: n

      n = size(sendbuf)

      call MPI_Sendrecv(sendbuf,n,MPI_DOUBLE_PRECISION,dest,101,   &
      recvbuf,n,MPI_DOUBLE_PRECISION,source,101, &
      xmpi_comm,MPI_STATUS_IGNORE,ierror)

   end subroutine xmpi_sendrecv_r2

   subroutine xmpi_sendrecv_r3(sendbuf,dest,recvbuf,source)
      implicit none
      real*8, dimension(:,:,:), intent(in)  :: sendbuf
      real*8, dimension(:,:,:), intent(out) :: recvbuf
      integer, intent(in)                   :: dest,source

      integer                               :: ierror
      integer                               :: n

      n = size(sendbuf)

      call MPI_Sendrecv(sendbuf,n,MPI_DOUBLE_PRECISION,dest,101,   &
      recvbuf,n,MPI_DOUBLE_PRECISION,source,101, &
      xmpi_comm,MPI_STATUS_IGNORE,ierror)

   end subroutine xmpi_sendrecv_r3

   subroutine xmpi_sendrecv_i1(sendbuf,dest,recvbuf,source)
      implicit none
      integer, dimension(:), intent(in)  :: sendbuf
      integer, dimension(:), intent(out) :: recvbuf
      integer, intent(in)                :: dest,source

      integer                            :: ierror
      integer                            :: n

      n = size(sendbuf)

      call MPI_Sendrecv(sendbuf,n,MPI_INTEGER,dest,102,   &
      recvbuf,n,MPI_INTEGER,source,102, &
      xmpi_comm,MPI_STATUS_IGNORE,ierror)

   end subroutine xmpi_sendrecv_i1

   subroutine xmpi_sendrecv_i2(sendbuf,dest,recvbuf,source)
      implicit none
      integer, dimension(:,:), intent(in)  :: sendbuf
      integer, dimension(:,:), intent(out) :: recvbuf
      integer, intent(in)                  :: dest,source

      integer                              :: ierror
      integer                              :: n

      n = size(sendbuf)

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

   subroutine xmpi_allreduce_r1(x,op)
      implicit none
      real*8,dimension(:), intent(inout)  :: x
      real*8,dimension(:), allocatable    :: y
      integer,intent(in)    :: op

      integer :: ierror
      allocate(y(size(x)))
      y = x
      call MPI_Allreduce(y,x,size(x),MPI_DOUBLE_PRECISION,op,xmpi_comm,ierror)
      deallocate(y)
   end subroutine xmpi_allreduce_r1

   subroutine xmpi_allreduce_i0(x,op)
      implicit none
      integer,intent(inout)  :: x
      integer,intent(in)    :: op

      integer :: y
      integer :: ierror
      y = x
      call MPI_Allreduce(y,x,1,MPI_INTEGER,op,xmpi_comm,ierror)
   end subroutine xmpi_allreduce_i0

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

   subroutine xmpi_reduce_r2(x,y,op)
      implicit none
      real*8,dimension(:,:), intent(in)  :: x
      real*8,dimension(:,:), intent(out) :: y
      integer, intent(in)                :: op

      integer :: ierror
      call MPI_Reduce(x,y,size(x),MPI_DOUBLE_PRECISION,op,xmpi_master,xmpi_comm,ierror)
   end subroutine xmpi_reduce_r2

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
         call xmpi_sendrecv(x(2,:),xmpi_top,    x(m,:),xmpi_bot)
       case('d','1:')
         call xmpi_sendrecv(x(m-1,:),xmpi_bot,  x(1,:),xmpi_top)
       case('l',':n')
         call xmpi_sendrecv(x(:,2),xmpi_left,   x(:,n),xmpi_right)
       case('r',':1')
         call xmpi_sendrecv(x(:,n-1),xmpi_right,x(:,1),xmpi_left)
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
         call xmpi_sendrecv(x(2,:),   xmpi_top,  x(m,:),xmpi_bot)
       case('d','1:')
         call xmpi_sendrecv(x(m-1,:), xmpi_bot,  x(1,:),xmpi_top)
       case('l',':n')
         call xmpi_sendrecv(x(:,2),   xmpi_left, x(:,n),xmpi_right)
       case('r',':1')
         call xmpi_sendrecv(x(:,n-1), xmpi_right,x(:,1),xmpi_left)
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
         call xmpi_sendrecv(x(2,:,:),  xmpi_top,    x(m,:,:),xmpi_bot)
       case('d','1:')
         call xmpi_sendrecv(x(m-1,:,:),xmpi_bot,  x(1,:,:),xmpi_top)
       case('l',':n')
         call xmpi_sendrecv(x(:,2,:),  xmpi_left,   x(:,n,:),xmpi_right)
       case('r',':1')
         call xmpi_sendrecv(x(:,n-1,:),xmpi_right,x(:,1,:),xmpi_left)
       case default
         if(xmaster) then
            write (*,*) 'Invalid direction parameter for xmpi_shift_r3: "'// &
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
         call xmpi_sendrecv(x(2,:,:),  xmpi_top,  x(m,:,:),xmpi_bot)
       case('d','1:')
         call xmpi_sendrecv(x(m-1,:,:),xmpi_bot,  x(1,:,:),xmpi_top)
       case('l',':n')
         call xmpi_sendrecv(x(:,2,:),  xmpi_left, x(:,n,:),xmpi_right)
       case('r',':1')
         call xmpi_sendrecv(x(:,n-1,:),xmpi_right,x(:,1,:),xmpi_left)
       case default
         if(xmaster) then
            write (*,*) 'Invalid direction parameter for xmpi_shift_i3: "'// &
            direction//'"'
            call halt_program
         endif
      end select

   end subroutine xmpi_shift_i3

   !  translation from nx+1, ny+1 to m and n
   !
   !  m=nx+1  nx=m-1
   !  n=ny+1  ny=n-1
   !
   !  variabele   l->r      r->l        b->t       t->b
   !              r         l           u          d      shift
   !          SHIFT_Y_R  SHIFT_Y_L  SHIFT_X_U  SHIFT_X_D


   !    ee,rr    m-3:m-2,:  3:4,:      :,n-3:n-2  :,3:4
   !             1:2,:      m-1:m,:    :,1:2      :,n-1:n
   !
   !    uu       m-3,:      2:3,:      :,n-3:n-2  :,3:4
   !             1,:        m-2,m-1,:  :,1:2      :,n-1:n
   !
   !    vv       m-3:m-2,:  3:4,:      :,n-3      :,2:3
   !             1:2,:      m-1:m,:    :,1        :,n-2:n-1
   !
   !  zs,ccg,zb  m-2,:      3,:        :,n-2      :,3
   !             2,:        m-1,:      :,2        :,n-1
   !
   !  Dus:
   !        i     p    j     q
   !
   !        1 <-> m-3  1 <-> n-3
   !        2 <-> m-2  2 <-> n-2
   !        3 <-> m-1  3 <-> n-1
   !        4 <-> m    4 <-> n
   !
   !        p=m-4+i    q=n-4+j

   !   About the _l subroutines:
   !
   !   xmpi_shift(x,SHIFT_Y_R,1,2): columns 1:2   are filled in (aka updated)
   !   xmpi_shift(x,SHIFT_Y_L,3,4): columns n-1,n are filled in (aka updated)
   !   xmpi_shift(x,SHIFT_X_D,1,2): rows    1:2   are filled in (aka updated)
   !   xmpi_shift(x,SHIFT_X_U,3,4): rows    m-1,m are filled in (aka updated)
   !
   subroutine xmpi_shift_r2_l(x,direction,i1,i2)
      implicit none
      real*8,dimension(:,:),intent(inout) :: x
      integer, intent(in)                 :: direction
      integer, intent(in)                 :: i1,i2

      integer            :: m,n,s1,s2,r1,r2
      integer, parameter :: nbord = 2, nover = 2*nbord

      ! nover is the number of overlapping rows/columns

      ! s1,s2 will contain the indices of the first and last row/column to send
      ! r1,r2 will contain the indices of the first and last row/column to send
      ! mn will contain the 1st or 2nd dimension as appropriate:
      !   shift in x direction: mn = m (= nx+1)
      !   shift in y direction: mn = n (= ny+1)

      m = size(x,1)
      n = size(x,2)

      ! sanity check
      select case(direction)
       case(SHIFT_Y_R,SHIFT_Y_L,SHIFT_X_U,SHIFT_X_D)
         continue
       case default
         if (xmaster) then
            write(*,*) 'Invalid value for direction in xmpi_shift_r2_l ',direction
            call halt_program
         endif
      endselect

      select case(direction)
       case(SHIFT_Y_R)
         s1 = n - nover + i1
         s2 = n - nover + i2
         r1 = i1
         r2 = i2
         call xmpi_sendrecv(x(:,s1:s2),xmpi_right,x(:,r1:r2),xmpi_left)
       case(SHIFT_X_D)
         s1 = m - nover + i1
         s2 = m - nover + i2
         r1 = i1
         r2 = i2
         call xmpi_sendrecv(x(s1:s2,:),xmpi_bot,  x(r1:r2,:),xmpi_top)
       case(SHIFT_Y_L)
         s1 = i1
         s2 = i2
         r1 = n - nover + i1
         r2 = n - nover + i2
         call xmpi_sendrecv(x(:,s1:s2),xmpi_left, x(:,r1:r2),xmpi_right)
       case(SHIFT_X_U)
         s1 = i1
         s2 = i2
         r1 = m - nover + i1
         r2 = m - nover + i2
         call xmpi_sendrecv(x(s1:s2,:),xmpi_top,  x(r1:r2,:),xmpi_bot)
      endselect

   end subroutine xmpi_shift_r2_l

   subroutine xmpi_shift_r3_l(x,direction,i1,i2)
      !
      implicit none
      real*8,dimension(:,:,:),intent(inout) :: x
      integer, intent(in)                   :: direction
      integer, intent(in)                   :: i1,i2

      integer            :: m,n,s1,s2,r1,r2
      integer, parameter :: nbord = 2, nover = 2*nbord


      ! nover is the number of overlapping rows/columns

      ! s1,s2 will contain the indices of the first and last row/column to send
      ! r1,r2 will contain the indices of the first and last row/column to send
      ! mn will contain the 1st or 2nd dimension as appropriate:
      !   shift in x direction: mn = m (= nx+1)
      !   shift in y direction: mn = n (= ny+1)

      m = size(x,1)
      n = size(x,2)

      ! sanity check
      select case(direction)
       case(SHIFT_Y_R,SHIFT_Y_L,SHIFT_X_U,SHIFT_X_D)
         continue
       case default
         if (xmaster) then
            write(*,*) 'Invalid value for direction in xmpi_shift_r2_l ',direction
            call halt_program
         endif
      endselect

      select case(direction)
       case(SHIFT_Y_R)
         s1 = n - nover + i1
         s2 = n - nover + i2
         r1 = i1
         r2 = i2
         call xmpi_sendrecv(x(:,s1:s2,:),xmpi_right,x(:,r1:r2,:),xmpi_left)
       case(SHIFT_X_D)
         s1 = m - nover + i1
         s2 = m - nover + i2
         r1 = i1
         r2 = i2
         call xmpi_sendrecv(x(s1:s2,:,:),xmpi_bot,  x(r1:r2,:,:),xmpi_top)
       case(SHIFT_Y_L)
         s1 = i1
         s2 = i2
         r1 = n - nover + i1
         r2 = n - nover + i2
         call xmpi_sendrecv(x(:,s1:s2,:),xmpi_left, x(:,r1:r2,:),xmpi_right)
       case(SHIFT_X_U)
         s1 = i1
         s2 = i2
         r1 = m - nover + i1
         r2 = m - nover + i2
         call xmpi_sendrecv(x(s1:s2,:,:),xmpi_top,  x(r1:r2,:,:),xmpi_bot)
      endselect
   end subroutine xmpi_shift_r3_l

   subroutine xmpi_shift_ee_r2(x)
      implicit none
      real*8,dimension(:,:),intent(inout) :: x
#ifdef SHIFT_TIMER
      real*8 ttt
      call xmpi_barrier
      ttt=MPI_Wtime()
#endif
      call xmpi_shift(x,SHIFT_Y_R,1,2)
      call xmpi_shift(x,SHIFT_Y_L,3,4)
      call xmpi_shift(x,SHIFT_X_U,3,4)
      call xmpi_shift(x,SHIFT_X_D,1,2)
#ifdef SHIFT_TIMER
      call xmpi_barrier
      if(xmaster)print *,'shift_ee_r2:',MPI_Wtime()-ttt,size(x,1),size(x,2)
#endif
   end subroutine xmpi_shift_ee_r2

   subroutine xmpi_shift_ee_r3(x)
      implicit none
      real*8,dimension(:,:,:),intent(inout) :: x
#ifdef SHIFT_TIMER
      real*8 ttt
      call xmpi_barrier
      ttt=MPI_Wtime()
#endif
      call xmpi_shift(x,SHIFT_Y_R,1,2)
      call xmpi_shift(x,SHIFT_Y_L,3,4)
      call xmpi_shift(x,SHIFT_X_U,3,4)
      call xmpi_shift(x,SHIFT_X_D,1,2)
#ifdef SHIFT_TIMER
      call xmpi_barrier
      if(xmaster)print *,'shift_ee_r3:',MPI_Wtime()-ttt
#endif
   end subroutine xmpi_shift_ee_r3

   subroutine xmpi_shift_uu_r2(x)
      implicit none
      real*8,dimension(:,:),intent(inout) :: x
#ifdef SHIFT_TIMER
      real*8 ttt
      call xmpi_barrier
      ttt=MPI_Wtime()
#endif
      call xmpi_shift(x,SHIFT_Y_R,1,2)
      call xmpi_shift(x,SHIFT_Y_L,3,4)
      call xmpi_shift(x,SHIFT_X_U,2,3)
      call xmpi_shift(x,SHIFT_X_D,1,1)
#ifdef SHIFT_TIMER
      call xmpi_barrier
      if(xmaster)print *,'shift_uu_r2:',MPI_Wtime()-ttt
#endif
   end subroutine xmpi_shift_uu_r2

   subroutine xmpi_shift_uu_r3(x)
      implicit none
      real*8,dimension(:,:,:),intent(inout) :: x
#ifdef SHIFT_TIMER
      real*8 ttt
      call xmpi_barrier
      ttt=MPI_Wtime()
#endif
      call xmpi_shift(x,SHIFT_Y_R,1,2)
      call xmpi_shift(x,SHIFT_Y_L,3,4)
      call xmpi_shift(x,SHIFT_X_U,2,3)
      call xmpi_shift(x,SHIFT_X_D,1,1)
#ifdef SHIFT_TIMER
      call xmpi_barrier
      if(xmaster)print *,'shift_uu_r3:',MPI_Wtime()-ttt
#endif
   end subroutine xmpi_shift_uu_r3

   subroutine xmpi_shift_vv_r2(x)
      implicit none
      real*8,dimension(:,:),intent(inout) :: x
#ifdef SHIFT_TIMER
      real*8 ttt
      call xmpi_barrier
      ttt=MPI_Wtime()
#endif
      call xmpi_shift(x,SHIFT_Y_R,1,1)
      call xmpi_shift(x,SHIFT_Y_L,2,3)
      call xmpi_shift(x,SHIFT_X_U,3,4)
      call xmpi_shift(x,SHIFT_X_D,1,2)
#ifdef SHIFT_TIMER
      call xmpi_barrier
      if(xmaster)print *,'shift_vv_r2:',MPI_Wtime()-ttt
#endif
   end subroutine xmpi_shift_vv_r2

   subroutine xmpi_shift_vv_r3(x)
      implicit none
      real*8,dimension(:,:,:),intent(inout) :: x
#ifdef SHIFT_TIMER
      real*8 ttt
      call xmpi_barrier
      ttt=MPI_Wtime()
#endif
      call xmpi_shift(x,SHIFT_Y_R,1,1)
      call xmpi_shift(x,SHIFT_Y_L,2,3)
      call xmpi_shift(x,SHIFT_X_U,3,4)
      call xmpi_shift(x,SHIFT_X_D,1,2)
#ifdef SHIFT_TIMER
      call xmpi_barrier
      if(xmaster)print *,'shift_vv_r3:',MPI_Wtime()-ttt
#endif
   end subroutine xmpi_shift_vv_r3

   subroutine xmpi_shift_zs_r2(x)
      implicit none
      real*8,dimension(:,:),intent(inout) :: x
#ifdef SHIFT_TIMER
      real*8 ttt
      call xmpi_barrier
      ttt=MPI_Wtime()
#endif
      call xmpi_shift(x,SHIFT_Y_R,2,2)
      call xmpi_shift(x,SHIFT_Y_L,3,3)
      call xmpi_shift(x,SHIFT_X_U,3,3)
      call xmpi_shift(x,SHIFT_X_D,2,2)
#ifdef SHIFT_TIMER
      call xmpi_barrier
      if(xmaster)print *,'shift_zs_r2:',MPI_Wtime()-ttt
#endif
   end subroutine xmpi_shift_zs_r2

   subroutine xmpi_shift_zs_r3(x)
      implicit none
      real*8,dimension(:,:,:),intent(inout) :: x
#ifdef SHIFT_TIMER
      real*8 ttt
      call xmpi_barrier
      ttt=MPI_Wtime()
#endif
      call xmpi_shift(x,SHIFT_Y_R,2,2)
      call xmpi_shift(x,SHIFT_Y_L,3,3)
      call xmpi_shift(x,SHIFT_X_U,3,3)
      call xmpi_shift(x,SHIFT_X_D,2,2)
#ifdef SHIFT_TIMER
      call xmpi_barrier
      if(xmaster)print *,'shift_zs_r3:',MPI_Wtime()-ttt
#endif
   end subroutine xmpi_shift_zs_r3
   !________________________________________________________________________________

   subroutine xmpi_send_r0(from,to,x)
      ! use this to send to one process, in communicator xmpi_ocomm
      ! receiver and sender call this same subroutine with the same
      ! from and to
      integer, intent(in)    :: from,to
      real*8, intent(inout)  :: x

      integer ier
      ! MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)

      if (from .eq. to) return

      if (xmpi_orank .eq. from) then
         call MPI_Send(x, 1, MPI_DOUBLE_PRECISION, to, 1011, xmpi_ocomm, ier)
      elseif (xmpi_orank .eq. to) then
         call MPI_Recv(x, 1, MPI_DOUBLE_PRECISION, from, 1011, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
      endif
   end subroutine xmpi_send_r0
   !________________________________________________________________________________

   subroutine xmpi_send_i0(from,to,x)
      ! use this to send to one process, in communicator xmpi_ocomm
      ! receiver and sender call this same subroutine with the same
      ! from and to
      integer, intent(in)     :: from,to
      integer, intent(inout)  :: x

      integer ier
      ! MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)

      if (from .eq. to) return

      if (xmpi_orank .eq. from) then
         call MPI_Send(x, 1, MPI_INTEGER, to, 1012, xmpi_ocomm, ier)
      elseif (xmpi_orank .eq. to) then
         call MPI_Recv(x, 1, MPI_INTEGER, from, 1012, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
      endif
   end subroutine xmpi_send_i0
   !________________________________________________________________________________

   subroutine xmpi_send_l0(from,to,x)
      ! use this to send to one process, in communicator xmpi_ocomm
      ! receiver and sender call this same subroutine with the same
      ! from and to
      integer, intent(in)     :: from,to
      logical, intent(inout)  :: x

      integer ier
      ! MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)

      if (from .eq. to) return

      if (xmpi_orank .eq. from) then
         call MPI_Send(x, 1, MPI_LOGICAL, to, 1013, xmpi_ocomm, ier)
      elseif (xmpi_orank .eq. to) then
         call MPI_Recv(x, 1, MPI_LOGICAL, from, 1013, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
      endif
   end subroutine xmpi_send_l0
   !________________________________________________________________________________

   subroutine xmpi_send_r1(from,to,x)
      ! use this to send to one process, in communicator xmpi_ocomm
      ! receiver and sender call this same subroutine with the same
      ! from and to
      integer, intent(in)    :: from,to
      real*8, intent(inout)  :: x(:)

      integer ier
      ! MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)

      if (from .eq. to) return

      if (xmpi_orank .eq. from) then
         call MPI_Send(x, size(x), MPI_DOUBLE_PRECISION, to, 1014, xmpi_ocomm, ier)
      elseif (xmpi_orank .eq. to) then
         call MPI_Recv(x, size(x), MPI_DOUBLE_PRECISION, from, 1014, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
      endif
   end subroutine xmpi_send_r1
   !________________________________________________________________________________

   subroutine xmpi_send_r2(from,to,x)
      ! use this to send to one process, in communicator xmpi_ocomm
      ! receiver and sender call this same subroutine with the same
      ! from and to
      integer, intent(in)    :: from,to
      real*8, intent(inout)  :: x(:,:)

      integer ier
      ! MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)

      if (from .eq. to) return

      if (xmpi_orank .eq. from) then
         call MPI_Send(x, size(x), MPI_DOUBLE_PRECISION, to, 1018, xmpi_ocomm, ier)
      elseif (xmpi_orank .eq. to) then
         call MPI_Recv(x, size(x), MPI_DOUBLE_PRECISION, from, 1018, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
      endif
   end subroutine xmpi_send_r2
   !________________________________________________________________________________

   subroutine xmpi_send_r3(from,to,x)
      ! use this to send to one process, in communicator xmpi_ocomm
      ! receiver and sender call this same subroutine with the same
      ! from and to
      integer, intent(in)    :: from,to
      real*8, intent(inout)  :: x(:,:,:)

      integer ier
      ! MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)

      if (from .eq. to) return

      if (xmpi_orank .eq. from) then
         call MPI_Send(x, size(x), MPI_DOUBLE_PRECISION, to, 1019, xmpi_ocomm, ier)
      elseif (xmpi_orank .eq. to) then
         call MPI_Recv(x, size(x), MPI_DOUBLE_PRECISION, from, 1019, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
      endif
   end subroutine xmpi_send_r3
   !________________________________________________________________________________

   subroutine xmpi_send_r4(from,to,x)
      ! use this to send to one process, in communicator xmpi_ocomm
      ! receiver and sender call this same subroutine with the same
      ! from and to
      integer, intent(in)    :: from,to
      real*8, intent(inout)  :: x(:,:,:,:)

      integer ier
      ! MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)

      if (from .eq. to) return

      if (xmpi_orank .eq. from) then
         call MPI_Send(x, size(x), MPI_DOUBLE_PRECISION, to, 1020, xmpi_ocomm, ier)
      elseif (xmpi_orank .eq. to) then
         call MPI_Recv(x, size(x), MPI_DOUBLE_PRECISION, from, 1020, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
      endif
   end subroutine xmpi_send_r4
   !________________________________________________________________________________

   subroutine xmpi_send_i1(from,to,x)
      ! use this to send to one process, in communicator xmpi_ocomm
      ! receiver and sender call this same subroutine with the same
      ! from and to
      integer, intent(in)     :: from,to
      integer, intent(inout)  :: x(:)

      integer ier
      ! MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)

      if (from .eq. to) return

      if (xmpi_orank .eq. from) then
         call MPI_Send(x, size(x), MPI_INTEGER, to, 1015, xmpi_ocomm, ier)
      elseif (xmpi_orank .eq. to) then
         call MPI_Recv(x, size(x), MPI_INTEGER, from, 1015, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
      endif
   end subroutine xmpi_send_i1
   !________________________________________________________________________________

   subroutine xmpi_send_l1(from,to,x)
      ! use this to send to one process, in communicator xmpi_ocomm
      ! receiver and sender call this same subroutine with the same
      ! from and to
      integer, intent(in)     :: from,to
      logical, intent(inout)  :: x(:)

      integer ier
      ! MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)

      if (from .eq. to) return

      if (xmpi_orank .eq. from) then
         call MPI_Send(x, size(x), MPI_LOGICAL, to, 1016, xmpi_ocomm, ier)
      elseif (xmpi_orank .eq. to) then
         call MPI_Recv(x, size(x), MPI_LOGICAL, from, 1016, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
      endif
   end subroutine xmpi_send_l1
   !________________________________________________________________________________

   ! following sends one integer from 'from' to 'to'.
   ! The receive by 'to' is done using a sleep/MPI_Test cycle
   ! to prevent cpu usage when waiting for the message.
   ! 'from' and 'to' in xmpi_ocomm

   subroutine xmpi_send_sleep(from,to)
      use sleeper
      integer, intent(in) :: from,to

      integer             :: ier, request
      integer             :: buf = 0
      logical             :: flag

      if (xmpi_orank .eq. from) then
         call MPI_Send(buf, 1, MPI_INTEGER, to, 1017, xmpi_ocomm, ier)
      elseif (xmpi_orank .eq. to) then
         call MPI_Irecv(buf, 1, MPI_INTEGER, from, 1017, xmpi_ocomm, request, ier)
         do
            call MPI_test(request, flag, MPI_STATUS_IGNORE, ier)
            if (flag) exit
            call myusleep(2000) ! sleep 2 milliseconds
         enddo
      endif
   end subroutine xmpi_send_sleep
   !________________________________________________________________________________

   subroutine testje
      implicit none
      real*8, dimension(10,10) :: x
      real*8, dimension(10,10,10) :: y
      call xmpi_shift_ee(x)
      call xmpi_shift_uu(x)
      call xmpi_shift_vv(x)
      call xmpi_shift_zs(x)
      call xmpi_shift_ee(y)
      call xmpi_shift_uu(y)
      call xmpi_shift_vv(y)
      call xmpi_shift_zs(y)
   end subroutine testje

   subroutine xmpi_barrier(toall)
      implicit none
      logical, intent(in), optional :: toall
      integer ierror,comm
      comm = xmpi_comm
      if (present(toall)) then
         if(toall) then
            comm = xmpi_ocomm
         endif
      endif
      call MPI_Barrier(comm,ierror)
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
   subroutine halt_program(normal)
      logical,intent(in),optional :: normal
      logical                     :: lnormal
      integer                     :: ierr

      if(present(normal)) then
         lnormal = normal
      else
         lnormal = .true.
      endif

      write(0,*) 'halt_program called by process', xmpi_orank
      if (lnormal) then
         call backtrace
      endif

#ifdef USEMPI
      if (lnormal) then
         call xmpi_abort
      else
         call MPI_Abort(xmpi_ocomm,1,ierr)
      endif
#else
      stop 1
#endif
   end subroutine halt_program

end module xmpi_module


#ifdef USEMPI
subroutine comm_errhandler(comm,error_code)
   use mpi
   use xmpi_module
   implicit none
   integer comm,error_code

   character(MPI_MAX_ERROR_STRING) e
   character(MPI_MAX_OBJECT_NAME) o
   integer r,ier,ra

   call MPI_Comm_get_name(comm,o,r,ier)
   call MPI_Comm_rank(MPI_COMM_WORLD,ra,ier)
   write(0,*) 'MPI process #',ra,'in ',o(1:r),' generated an error:',error_code
   call MPI_Error_string(error_code,e,r,ier)
   write(0,*) e(1:r)
   call halt_program
end subroutine comm_errhandler
#endif

#define MBACKTRACE
#ifdef __GFORTRAN__
#if __GNUC__ >=4
#if __GNUC_MINOR__ >=8
#undef MBACKTRACE
#endif
#endif
#endif

#ifdef __INTEL_COMPILER
#undef MBACKTRACE
subroutine backtrace
   use ifcore
   call tracebackqq('traceback:',-1)
end subroutine backtrace
#endif

#ifdef MBACKTRACE
subroutine backtrace
   implicit none
   integer, pointer :: x => null()
   write(0,*) 'no backtrace for this compiler'
   write(0,*) 'forcing a segmentation fault ...'
   x = 0
   print *,x
end subroutine backtrace
#endif
