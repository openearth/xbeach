program testgenmodule
use mpi
use xmpi_module
use general_mpi_module
implicit none
real*8, external :: fill,fill1
integer, parameter :: ma=101
integer, parameter :: na=501
real *8, allocatable, dimension(:,:) :: a,b
integer, allocatable, dimension(:,:) :: ia,ib
real *8, allocatable, dimension(:) :: vs,vr
integer, allocatable, dimension(:) :: is,js,lm,ln
integer i,j,p,ii,jj,errors,total_errors,ierror
integer mlocal,nlocal
integer, allocatable, dimension(:) :: leftn,rightn,topn,botn
logical, allocatable, dimension(:) :: isleft,isright,istop,isbot
real*8 t0,t1
integer maal,times

times=1000

call xmpi_initialize

allocate (leftn(xmpi_size))
allocate (rightn(xmpi_size))
allocate (topn(xmpi_size))
allocate (botn(xmpi_size))

allocate (isleft(xmpi_size))
allocate (isright(xmpi_size))
allocate (istop(xmpi_size))
allocate (isbot(xmpi_size))

if (xmaster) then
  write(*,*)'xmpi_size:',xmpi_size
  ! call flush
endif
call xmpi_determine_processor_grid
call MPI_Gather((/xmpi_left/),1,MPI_INTEGER,&
                leftn,1,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)
call MPI_Gather((/xmpi_right/),1,MPI_INTEGER,&
                rightn,1,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)
call MPI_Gather((/xmpi_top/),1,MPI_INTEGER,&
                topn,1,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)
call MPI_Gather((/xmpi_bot/),1,MPI_INTEGER,&
                botn,1,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)

if(xmaster) then
  write(*,*)ma,na,xmpi_size,xmpi_m, xmpi_n
  write(*,*)'proc left right top bot'
  do i=1,xmpi_size
    write(*,*)i-1,leftn(i),rightn(i),topn(i),botn(i)
  enddo
  ! call flush()
endif

if (xmpi_size .ne. xmpi_m*xmpi_n) then
  write(*,*)'xmpisize:',xmpi_size,' m:',xmpi_m,' n',xmpi_n
  call xmpi_finalize
  stop
endif

allocate(is(xmpi_size),lm(xmpi_size),js(xmpi_size),ln(xmpi_size))

if (xmpi_rank .eq. xmpi_master) then

  call det_submatrices(ma,na,xmpi_m,xmpi_n,is,lm,js,ln, &
                isleft,isright,istop,isbot)
  write(*,*)'verdeling'
  i=0
  j=1
  write(*,*)'    p    i    j   is   lm   js   ln left right  top  bot'
  do p=0,xmpi_size-1
    i=i+1
    if (i .gt. xmpi_m) then
      j=j+1
      i=1
    endif

    write(*,"(' ',7i5,4l5)")p,i,j,is(p+1),lm(p+1),js(p+1),ln(p+1), &
          isleft(p+1),isright(p+1),istop(p+1),isbot(p+1)
    ! call flush()
  enddo

  allocate(a(ma,na),ia(ma,na))
  do j=1,na
    do i=1,ma
      a(i,j)  =  fill(i,j)
      ia(i,j) =  fill(i,j)
    enddo
  enddo
else
  allocate(a(1,1),ia(1,1))
endif

call MPI_Scatter(ln(1), 1, MPI_INTEGER, nlocal, 1, MPI_INTEGER, &
                xmpi_master, xmpi_comm, ierror)
call MPI_Scatter(lm(1), 1, MPI_INTEGER, mlocal, 1, MPI_INTEGER, &
                xmpi_master, xmpi_comm, ierror)
call MPI_Bcast(ln, xmpi_size, MPI_INTEGER, xmpi_master, &
               xmpi_comm, ierror)
call MPI_Bcast(lm, xmpi_size, MPI_INTEGER, xmpi_master, &
               xmpi_comm, ierror)
call MPI_Bcast(isleft, xmpi_size, MPI_LOGICAL, xmpi_master, &
               xmpi_comm, ierror)
call MPI_Bcast(isright, xmpi_size, MPI_LOGICAL, xmpi_master, &
               xmpi_comm, ierror)
call MPI_Bcast(istop, xmpi_size, MPI_LOGICAL, xmpi_master, &
               xmpi_comm, ierror)
call MPI_Bcast(isbot, xmpi_size, MPI_LOGICAL, xmpi_master, &
               xmpi_comm, ierror)
write(*,*)xmpi_rank,'mlocal,nlocal:',mlocal,nlocal
! call flush()
allocate(b(mlocal,nlocal),ib(mlocal,nlocal))

call MPI_Barrier(xmpi_comm,ierror)
t0 = MPI_Wtime()
do maal = 1,times
  call matrix_distr(a,b,is,lm,js,ln,xmpi_master,xmpi_comm)
enddo
call MPI_Barrier(xmpi_comm,ierror)
t1 = MPI_Wtime()
if (xmpi_master .eq. xmpi_rank) then
   write(*,*)'Time for real*8 scattering:',(t1-t0)/times
endif

call MPI_Barrier(xmpi_comm,ierror)
t0 = MPI_Wtime()
do maal = 1,times
  call matrix_distr(ia,ib,is,lm,js,ln,xmpi_master,xmpi_comm)
enddo
call MPI_Barrier(xmpi_comm,ierror)
t1 = MPI_Wtime()
if (xmpi_master .eq. xmpi_rank) then
   write(*,*)'Time for integer scattering:',(t1-t0)/times
endif

call MPI_Bcast(is,xmpi_size,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)
call MPI_Bcast(lm,xmpi_size,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)
call MPI_Bcast(js,xmpi_size,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)
call MPI_Bcast(ln,xmpi_size,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)
errors=0
do j=1,nlocal
  do i=1,mlocal
    ii = is(xmpi_rank+1)+i-1
    jj = js(xmpi_rank+1)+j-1
    if (b(i,j) .ne. fill(ii,jj)) then
      errors = errors+1
      write(*,*)xmpi_rank,'scat real*8 error:', i,j,ii,jj,b(i,j),fill(ii,jj)
      ! call flush()
    endif
  enddo
enddo

call MPI_Reduce(errors,total_errors,1,MPI_INTEGER,MPI_SUM,&
                xmpi_master,xmpi_comm,ierror)
if (xmaster) then
  write(*,*)'Number of real*8 scatter errors is ',total_errors
  ! call flush()
endif

errors=0
do j=1,nlocal
  do i=1,mlocal
    ii = is(xmpi_rank+1)+i-1
    jj = js(xmpi_rank+1)+j-1
    if (ib(i,j) .ne. fill(ii,jj)) then
      errors = errors+1
      write(*,*)xmpi_rank,'scat integer error:', i,j,ii,jj,ib(i,j),fill(ii,jj)
      ! call flush()
    endif
  enddo
enddo

call MPI_Reduce(errors,total_errors,1,MPI_INTEGER,MPI_SUM,&
                xmpi_master,xmpi_comm,ierror)
if (xmaster) then
  write(*,*)'Number of integer scatter errors is ',total_errors
  ! call flush()
endif

if (xmpi_master .eq. xmpi_rank) then
  a=-1
endif


  if ( .not. xmpi_isleft) then
    do i=1,mlocal
      b(i,1) = -1
    enddo
  endif
  if (.not. xmpi_isright) then
    do i=1,mlocal
      b(i,nlocal) = -1
    enddo
  endif
  if (.not. xmpi_istop) then
    do j=1,nlocal
      b(1,j) = -1
    enddo
  endif
  if (.not. xmpi_isbot) then
    do j=1,nlocal
      b(mlocal,j) = -1
    enddo
  endif

call MPI_Barrier(xmpi_comm, ierror)
t0=MPI_Wtime()

do maal=1,times
!call matrix_coll(a,b,is,lm,js,ln,xmpi_isleft,xmpi_isright,&
!    xmpi_istop,xmpi_isbot,xmpi_master,xmpi_comm)
call matrix_coll(a,b,is,lm,js,ln,isleft,isright,&
    istop,isbot,xmpi_master,xmpi_comm)
enddo
t1=MPI_Wtime()
if(xmaster) then

  write(*,*)'Time for collecting:',(t1-t0)/times
  ! call flush
  errors=0
  do j=1,na
    do i=1,ma
      if (fill(i,j).ne. a(i,j)) then
        errors=errors+1
        write(*,*)'coll error:',i,j,a(i,j),fill(i,j)
      endif
    enddo
  enddo
  write(*,*)'Number of collect errors is',errors
  ! call flush()
endif

! testing the test_borders subroutine

! fill again the global array a:
if(xmaster) then
  do j=1,na
    do i=1,ma
      a(i,j) =  fill(i,j)
    enddo
  enddo
endif

! distribute the matrix
call matrix_distr(a,b,is,lm,js,ln,xmpi_master,xmpi_comm)

! distribute is,lm,js,ln:
call MPI_Bcast(is,xmpi_size,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)
call MPI_Bcast(lm,xmpi_size,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)
call MPI_Bcast(js,xmpi_size,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)
call MPI_Bcast(ln,xmpi_size,MPI_INTEGER,xmpi_master,xmpi_comm,ierror)

! fill the colomuns and rows that are to be received from the
! neighbours (the borders, or ghost cells) with -1. After the
! shift, they should be filled as usual again
! Take care of processes ath the borders of the domain.
! The processes on the left will not receive their first column,
! for example.
!
! Also be prudent about the elements that are the corners of the
! whole domain.
if (xmpi_top .ne. MPI_PROC_NULL) then
  b(1,2:nlocal-1) = -1
endif
if (xmpi_left .ne. MPI_PROC_NULL) then
  b(2:mlocal-1,1) = -1
endif
if (xmpi_bot .ne. MPI_PROC_NULL) then
  b(mlocal,2:nlocal-1) = -1
endif
if (xmpi_right .ne. MPI_PROC_NULL) then
  b(2:mlocal-1,nlocal) = -1
endif

if (xmpi_top .ne. MPI_PROC_NULL .and. xmpi_left .ne. MPI_PROC_NULL) then
  b(1,1) = -1
endif
if (xmpi_left .ne. MPI_PROC_NULL .and. xmpi_bot .ne. MPI_PROC_NULL) then
  b(mlocal,1) = -1
endif
if (xmpi_bot .ne. MPI_PROC_NULL .and. xmpi_right .ne. MPI_PROC_NULL) then
  b(mlocal,nlocal) = -1
endif
if (xmpi_right .ne. MPI_PROC_NULL .and. xmpi_top .ne. MPI_PROC_NULL) then
  b(1,nlocal) = -1
endif

do maal=1,times
  call shift_borders(b,xmpi_left,xmpi_right,xmpi_top,xmpi_bot,xmpi_comm)
enddo

errors=0
do j=1,nlocal
  do i=1,mlocal
    ii = is(xmpi_rank+1)+i-1
    jj = js(xmpi_rank+1)+j-1
    if (b(i,j) .ne. fill(ii,jj)) then
      errors = errors+1
      write(*,*)xmpi_rank,'shift error:', i,j,ii,jj,b(i,j),fill(ii,jj)
      ! call flush()
    endif
  enddo
enddo

call MPI_Reduce(errors,total_errors,1,MPI_INTEGER,MPI_SUM,&
                xmpi_master,xmpi_comm,ierror)
if (xmaster) then
  write(*,*)'Number of shift errors is ',total_errors
  ! call flush()
endif

! test vector_distr_send

! allocate vector-to-send

if (xmpi_master .eq. xmpi_rank) then
  allocate (vs(ma),vr(mlocal))
  do i=1,ma
    vs(i) = fill1(i)
  enddo
else
  allocate(vs(1),vr(mlocal))
endif

vr = -1

do maal=1,times
call vector_distr_send(vs,vr,is,lm,xmpi_master,xmpi_comm)
enddo

errors=0
do i=1, mlocal
  if (vr(i) .ne. fill1(is(xmpi_rank+1)+i-1)) then
    errors = errors+1
    write(*,*)xmpi_rank,': vector distribute error',i,fill(is(xmpi_rank+1)),vr(i)
  endif
enddo

call MPI_Reduce(errors,total_errors,1,MPI_INTEGER,MPI_SUM,&
                xmpi_master,xmpi_comm,ierror)
if (xmaster) then
  write(*,*)'Number of vector errors is ',total_errors
  ! call flush()
endif

deallocate(is,lm,js,ln,a,b)
deallocate (leftn)
deallocate (rightn)
deallocate (topn)
deallocate (botn)
deallocate (vs,vr)
goto 100
100 continue
if (xmpi_master .eq. xmpi_rank) then
  write(*,*)xmpi_rank,': calling xmpi_finalize:'
endif
call xmpi_finalize

end

real*8 function fill(i,j)
  integer, intent(in) :: i,j

  fill = 10000*i+j

end function fill
real*8 function fill1(i)
  integer, intent(in) :: i

  fill1 = i*10

end function fill1

