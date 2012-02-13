module testgenmod
  implicit none
  integer, parameter ::  times=100
  integer, parameter :: ma=101
  integer, parameter :: na=501
  integer, allocatable, dimension(:) :: is,js,lm,ln
  integer, allocatable, dimension(:) :: leftn,rightn,topn,botn
  logical, allocatable, dimension(:) :: isleft,isright,istop,isbot
  integer mlocal,nlocal

contains

real*8 function fill(i,j)
  implicit none
  integer, intent(in) :: i,j

  fill = 10000*i+j

end function fill
real*8 function fill1(i)
  implicit none
  integer, intent(in) :: i

  fill1 = i*10

end function fill1

subroutine scattertest(a,ia,b,ib)
  use mpi
  use xmpi_module
  use general_mpi_module
  implicit none
  real*8, dimension(:,:) :: a,b
  integer, dimension(:,:) :: ia,ib
  integer i,maal,ii,jj,j,ierror
  integer errors, total_errors
  real*8 t0,t1
  if(xmaster) then
    do j=1,size(a,2)
      do i=1,size(a,1)
        a(i,j)  =  fill(i,j)
        ia(i,j) =  fill(i,j)
      enddo
    enddo
  endif
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
end subroutine scattertest

subroutine colltest(a,b)
  use xmpi_module
  use general_mpi_module
  implicit none
  real*8, dimension(:,:) :: a,b
  real*8 t0,t1
  integer i,j,maal,errors,ierror

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
end subroutine colltest

subroutine shifttest(a,b)
  use xmpi_module
  use general_mpi_module
  implicit none
  real*8, dimension(:,:) :: a,b
  integer i,j,ii,jj,maal,errors,ierror,total_errors
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

  ! fill the columns and rows that are to be received from the
  ! neighbours (the borders, or ghost cells) with -1. After the
  ! shift, they should be filled as usual again
  ! Take care of processes at the borders of the domain.
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
end subroutine shifttest

subroutine vscattertest(vs,vr)
  use xmpi_module
  use general_mpi_module
  implicit none
  real*8, dimension(:) :: vs,vr
  integer i,maal,errors,total_errors,ierror

  if(xmaster) then
    do i=1,size(vs)
      vs(i) = fill1(i)
    enddo
  endif
  vr = -1

  do maal=1,times
    call vector_distr_send(vs,vr,is,lm,xmpi_master,xmpi_comm)
  enddo

  errors=0
  do i=1, mlocal
    if (vr(i) .ne. fill1(is(xmpi_rank+1)+i-1)) then
      errors = errors+1
      write(*,*)xmpi_rank,': vector distribute error',i,fill1(is(xmpi_rank+1)),vr(i)
    endif
  enddo

  call MPI_Reduce(errors,total_errors,1,MPI_INTEGER,MPI_SUM,&
                  xmpi_master,xmpi_comm,ierror)
  if (xmaster) then
    write(*,*)'Number of vector errors is ',total_errors
    ! call flush()
  endif
end subroutine vscattertest

subroutine getrowtest(b,x)
  use mpi
  use xmpi_module
  implicit none
  real*8, dimension(:,:) :: b
  real*8, dimension(:)   :: x
  character, dimension(2) :: rows
  integer i,j,k,errors,total_errors,ierror
  if(xmaster) write(*,*) 'testing xmpi_getrow'
  rows(1) = '1'
  rows(2) = 'm'
  errors=0
  do i = 1,xmpi_m         !loop over source processes
    do j = 1,2            !loop over top, bot
      b = xmpi_prow
      x=-1
      call xmpi_getrow(b,nlocal,rows(j),i,x)
      do k=1,nlocal
        if (x(k) .ne. i) then
          errors = errors+1
          write(*,*) xmpi_rank,':error in getrow i,j,k,x(k)',i,j,k,x(k)
        endif
      enddo
    enddo
  enddo

  call MPI_Reduce(errors,total_errors,1,MPI_INTEGER,MPI_SUM,&
                  xmpi_master,xmpi_comm,ierror)
  if (xmaster) then
    write(*,*)'Number of getrow errors is ',total_errors
    ! call flush()
  endif
end subroutine getrowtest

subroutine xshiftest(b,ib,b3,ib3)
  use mpi
  use xmpi_module
  implicit none
  real*8, dimension(:,:) :: b
  real*8, dimension(:,:,:) :: b3
  integer, dimension(:,:) :: ib
  integer, dimension(:,:,:) :: ib3

  character*2, dimension(4) :: w = (/'m:','1:',':n',':1'/)
  character *2              :: ww
  integer :: i, mb,nb,errors,j,total_errors,ierror,iw,l,lb

  mb = size(b,1)
  nb = size(b,2)
  errors=0
  do iw=1,4
    ww = w(iw)
    if (xmaster) then
!      write(*,*) 'testing b',mb,nb,ww
    endif
    b = -1
    select case(ww)
    case ('m:')
      do j=1,nb
        b(2,j) = xmpi_rank+j
      enddo
    case ('1:')
      do j=1,nb
        b(mb-1,j) = xmpi_rank+j
      enddo
    case (':n')
      do i=1,mb
        b(i,2)  = xmpi_rank+i
      enddo
    case (':1')
      do i=1,mb
        b(i,nb-1)  = xmpi_rank+i
      enddo
    end select
    call xmpi_shift(b,ww)
    select case(ww)
    case('m:')
      if (xmpi_isbot) then
        do j=1,nb
          if (b(mb,j) .ne. -1) then
            write (*,*) xmpi_rank,mb,j,' ',ww,': Error in xshiftest, expected -1 got ',b(mb,j)
            errors = errors + 1
          endif
        enddo
      else
        do j=1,nb
          if (b(mb,j) .ne. xmpi_bot+j) then
            write (*,*) xmpi_rank,mb,j,' ',ww,': Error in xshiftest, expected',xmpi_bot+j,&
                                  ' got ',b(mb,j)
            errors = errors + 1
          endif
        enddo
      endif
    case('1:')
      if (xmpi_istop) then
        do j=1,nb
          if (b(1,j) .ne. -1) then
            write (*,*) xmpi_rank,1,j,' ',ww,': Error in xshiftest, expected -1 got ',b(1,j)
            errors = errors + 1
          endif
        enddo
      else
        do j=1,nb
          if (b(1,j) .ne. xmpi_top+j) then
            write (*,*) xmpi_rank,1,j,' ',ww,': Error in xshiftest, expected',xmpi_top+j,&
                                  ' got ',b(1,j)
            errors = errors + 1
          endif
        enddo
      endif
    case(':n')
      if (xmpi_isright) then
        do i=1,mb
          if (b(i,nb) .ne. -1) then
            write (*,*) xmpi_rank,i,nb,' ',ww,': Error in xshiftest, expected -1 got ',b(i,nb)
            errors = errors + 1
          endif
        enddo
      else
        do i=1,mb
          if (b(i,nb) .ne. xmpi_right+i) then
            write (*,*) xmpi_rank,i,nb,' ',ww,': Error in xshiftest, expected',xmpi_right+i,&
                                  ' got ',b(i,nb)
            errors = errors + 1
          endif
        enddo
      endif
    case(':1')
      if (xmpi_isleft) then
        do i=1,mb
          if (b(i,1) .ne. -1) then
            write (*,*) xmpi_rank,i,1,' ',ww,': Error in xshiftest, expected -1 got ',b(i,1)
            errors = errors + 1
          endif
        enddo
      else
        do i=1,mb
          if (b(i,1) .ne. xmpi_left+i) then
            write (*,*) xmpi_rank,i,1,' ',ww,': Error in xshiftest, expected',xmpi_left+i,&
                                  ' got ',b(i,1)
            errors = errors + 1
          endif
        enddo
      endif
    end select
  enddo
  mb = size(ib,1)
  nb = size(ib,2)
  errors=0
  do iw=1,4
    ww = w(iw)
    if (xmaster) then
!      write(*,*) 'testing ib',mb,nb,ww
    endif
    ib = -1
    select case(ww)
    case ('m:')
      do j=1,nb
        ib(2,j) = xmpi_rank+j
      enddo
    case ('1:')
      do j=1,nb
        ib(mb-1,j) = xmpi_rank+j
      enddo
    case (':n')
      do i=1,mb
        ib(i,2)  = xmpi_rank+i
      enddo
    case (':1')
      do i=1,mb
        ib(i,nb-1)  = xmpi_rank+i
      enddo
    end select
    call xmpi_shift(ib,ww)
    select case(ww)
    case('m:')
      if (xmpi_isbot) then
        do j=1,nb
          if (ib(mb,j) .ne. -1) then
            write (*,*) xmpi_rank,mb,j,' ',ww,': Error in xshiftest, expected -1 got ',ib(mb,j)
            errors = errors + 1
          endif
        enddo
      else
        do j=1,nb
          if (ib(mb,j) .ne. xmpi_bot+j) then
            write (*,*) xmpi_rank,mb,j,' ',ww,': Error in xshiftest, expected',xmpi_bot+j,&
                                  ' got ',ib(mb,j)
            errors = errors + 1
          endif
        enddo
      endif
    case('1:')
      if (xmpi_istop) then
        do j=1,nb
          if (ib(1,j) .ne. -1) then
            write (*,*) xmpi_rank,1,j,' ',ww,': Error in xshiftest, expected -1 got ',ib(1,j)
            errors = errors + 1
          endif
        enddo
      else
        do j=1,nb
          if (ib(1,j) .ne. xmpi_top+j) then
            write (*,*) xmpi_rank,1,j,' ',ww,': Error in xshiftest, expected',xmpi_top+j,&
                                  ' got ',ib(1,j)
            errors = errors + 1
          endif
        enddo
      endif
    case(':n')
      if (xmpi_isright) then
        do i=1,mb
          if (ib(i,nb) .ne. -1) then
            write (*,*) xmpi_rank,i,nb,' ',ww,': Error in xshiftest, expected -1 got ',ib(i,nb)
            errors = errors + 1
          endif
        enddo
      else
        do i=1,mb
          if (ib(i,nb) .ne. xmpi_right+i) then
            write (*,*) xmpi_rank,i,nb,' ',ww,': Error in xshiftest, expected',xmpi_right+i,&
                                  ' got ',ib(i,nb)
            errors = errors + 1
          endif
        enddo
      endif
    case(':1')
      if (xmpi_isleft) then
        do i=1,mb
          if (ib(i,1) .ne. -1) then
            write (*,*) xmpi_rank,i,1,' ',ww,': Error in xshiftest, expected -1 got ',ib(i,1)
            errors = errors + 1
          endif
        enddo
      else
        do i=1,mb
          if (ib(i,1) .ne. xmpi_left+i) then
            write (*,*) xmpi_rank,i,1,' ',ww,': Error in xshiftest, expected',xmpi_left+i,&
                                  ' got ',ib(i,1)
            errors = errors + 1
          endif
        enddo
      endif
    end select
  enddo

  mb = size(b3,1)
  nb = size(b3,2)
  lb = size(b3,3)
  errors=0
  do iw=1,4
    ww = w(iw)
    if (xmaster) then
!      write(*,*) 'testing b3',mb,nb,lb,ww
    endif
    b3 = -1
    select case(ww)
    case ('m:')
      do j=1,nb
        do l=1,lb
          b3(2,j,l) = xmpi_rank+j+l
        enddo
      enddo
    case ('1:')
      do j=1,nb
        do l=1,lb
          b3(mb-1,j,l) = xmpi_rank+j+l
        enddo
      enddo
    case (':n')
      do i=1,mb
        do l=1,lb
          b3(i,2,l)  = xmpi_rank+i+l
        enddo
      enddo
    case (':1')
      do i=1,mb
        do l=1,lb
          b3(i,nb-1,l)  = xmpi_rank+i+l
        enddo
      enddo
    end select
    call xmpi_shift(b3,ww)
    select case(ww)
    case('m:')
      if (xmpi_isbot) then
        do j=1,nb
          do l=1,lb
            if (b3(mb,j,l) .ne. -1) then
              write (*,*) xmpi_rank,mb,j,l,' ',ww,': Error in xshiftest, expected -1 got ',b3(mb,j,l)
              errors = errors + 1
            endif
          enddo
        enddo
      else
        do j=1,nb
          do l=1,lb
            if (b3(mb,j,l) .ne. xmpi_bot+j+l) then
              write (*,*) xmpi_rank,mb,j,l,' ',ww,': Error in xshiftest, expected',xmpi_bot+j+l,&
                                    ' got ',b3(mb,j,l)
              errors = errors + 1
            endif
          enddo
        enddo
      endif
    case('1:')
      if (xmpi_istop) then
        do j=1,nb
          do l=1,lb
            if (b3(1,j,l) .ne. -1) then
              write (*,*) xmpi_rank,1,j,l,' ',ww,': Error in xshiftest, expected -1 got ',b3(1,j,l)
              errors = errors + 1
            endif
          enddo
        enddo
      else
        do j=1,nb
          do l=1,lb
            if (b3(1,j,l) .ne. xmpi_top+j+l) then
              write (*,*) xmpi_rank,1,j,l,' ',ww,': Error in xshiftest, expected',xmpi_top+j+l,&
                                    ' got ',b3(1,j,l)
              errors = errors + 1
            endif
          enddo
        enddo
      endif
    case(':n')
      if (xmpi_isright) then
        do i=1,mb
          do l=1,lb
            if (b3(i,nb,l) .ne. -1) then
              write (*,*) xmpi_rank,i,nb,l,' ',ww,': Error in xshiftest, expected -1 got ',b3(i,nb,l)
              errors = errors + 1
            endif
          enddo
        enddo
      else
        do i=1,mb
          do l=1,lb
            if (b3(i,nb,l) .ne. xmpi_right+i+l) then
              write (*,*) xmpi_rank,i,nb,l,' ',ww,': Error in xshiftest, expected',xmpi_right+i+l,&
                                    ' got ',b3(i,nb,l)
              errors = errors + 1
            endif
          enddo
        enddo
      endif
    case(':1')
      if (xmpi_isleft) then
        do i=1,mb
          do l=1,lb
            if (b3(i,1,l) .ne. -1) then
              write (*,*) xmpi_rank,i,1,l,' ',ww,': Error in xshiftest, expected -1 got ',b3(i,1,l)
              errors = errors + 1
            endif
          enddo
        enddo
      else
        do i=1,mb
          do l=1,lb
            if (b3(i,1,l) .ne. xmpi_left+i+l) then
              write (*,*) xmpi_rank,i,1,l,' ',ww,': Error in xshiftest, expected',xmpi_left+i+l,&
                                    ' got ',b3(i,1,l)
              errors = errors + 1
            endif
          enddo
        enddo
      endif
    end select
  enddo
  mb = size(ib3,1)
  nb = size(ib3,2)
  lb = size(ib3,3)
  errors=0
  do iw=1,4
    ww = w(iw)
    if (xmaster) then
!      write(*,*) 'testing ib3',mb,nb,lb,ww
    endif
    ib3 = -1
    select case(ww)
    case ('m:')
      do j=1,nb
        do l=1,lb
          ib3(2,j,l) = xmpi_rank+j+l
        enddo
      enddo
    case ('1:')
      do j=1,nb
        do l=1,lb
          ib3(mb-1,j,l) = xmpi_rank+j+l
        enddo
      enddo
    case (':n')
      do i=1,mb
        do l=1,lb
          ib3(i,2,l)  = xmpi_rank+i+l
        enddo
      enddo
    case (':1')
      do i=1,mb
        do l=1,lb
          ib3(i,nb-1,l)  = xmpi_rank+i+l
        enddo
      enddo
    end select
    call xmpi_shift(ib3,ww)
    select case(ww)
    case('m:')
      if (xmpi_isbot) then
        do j=1,nb
          do l=1,lb
            if (ib3(mb,j,l) .ne. -1) then
              write (*,*) xmpi_rank,mb,j,l,' ',ww,': Error in xshiftest, expected -1 got ',ib3(mb,j,l)
              errors = errors + 1
            endif
          enddo
        enddo
      else
        do j=1,nb
          do l=1,lb
            if (ib3(mb,j,l) .ne. xmpi_bot+j+l) then
              write (*,*) xmpi_rank,mb,j,l,' ',ww,': Error in xshiftest, expected',xmpi_bot+j+l,&
                                    ' got ',ib3(mb,j,l)
              errors = errors + 1
            endif
          enddo
        enddo
      endif
    case('1:')
      if (xmpi_istop) then
        do j=1,nb
          do l=1,lb
            if (ib3(1,j,l) .ne. -1) then
              write (*,*) xmpi_rank,1,j,l,' ',ww,': Error in xshiftest, expected -1 got ',ib3(1,j,l)
              errors = errors + 1
            endif
          enddo
        enddo
      else
        do j=1,nb
          do l=1,lb
            if (ib3(1,j,l) .ne. xmpi_top+j+l) then
              write (*,*) xmpi_rank,1,j,l,' ',ww,': Error in xshiftest, expected',xmpi_top+j+l,&
                                    ' got ',ib3(1,j,l)
              errors = errors + 1
            endif
          enddo
        enddo
      endif
    case(':n')
      if (xmpi_isright) then
        do i=1,mb
          do l=1,lb
            if (ib3(i,nb,l) .ne. -1) then
              write (*,*) xmpi_rank,i,nb,l,' ',ww,': Error in xshiftest, expected -1 got ',ib3(i,nb,l)
              errors = errors + 1
            endif
          enddo
        enddo
      else
        do i=1,mb
          do l=1,lb
            if (ib3(i,nb,l) .ne. xmpi_right+i+l) then
              write (*,*) xmpi_rank,i,nb,l,' ',ww,': Error in xshiftest, expected',xmpi_right+i+l,&
                                    ' got ',ib3(i,nb,l)
              errors = errors + 1
            endif
          enddo
        enddo
      endif
    case(':1')
      if (xmpi_isleft) then
        do i=1,mb
          do l=1,lb
            if (ib3(i,1,l) .ne. -1) then
              write (*,*) xmpi_rank,i,1,l,' ',ww,': Error in xshiftest, expected -1 got ',ib3(i,1,l)
              errors = errors + 1
            endif
          enddo
        enddo
      else
        do i=1,mb
          do l=1,lb
            if (ib3(i,1,l) .ne. xmpi_left+i+l) then
              write (*,*) xmpi_rank,i,1,l,' ',ww,': Error in xshiftest, expected',xmpi_left+i+l,&
                                    ' got ',ib3(i,1,l)
              errors = errors + 1
            endif
          enddo
        enddo
      endif
    end select
  enddo

  call MPI_Reduce(errors,total_errors,1,MPI_INTEGER,MPI_SUM,&
                  xmpi_master,xmpi_comm,ierror)
  if (xmaster) then
    write(*,*)'Number of xshiftest errors is ',total_errors
    ! call flush()
  endif

end subroutine xshiftest

end module testgenmod

program testgenmodule
use mpi
use xmpi_module
use general_mpi_module
use testgenmod
implicit none
!integer, parameter :: ma=51
!integer, parameter :: na=61
real *8, allocatable, dimension(:,:) :: a,b
real *8, allocatable, dimension(:,:,:) :: b3,bb3
integer, allocatable, dimension(:,:,:) :: ib3,ibb3
integer, allocatable, dimension(:,:) :: ia,ib
real *8, allocatable, dimension(:,:) :: aa,bb
integer, allocatable, dimension(:,:) :: iaa,ibb
real *8, allocatable, dimension(:) :: vs,vr,x
real *8, allocatable, dimension(:) :: vvs,vvr,xx
integer i,j,p,ierror
real*8 MPI_Wtime

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
  allocate(aa(ma+1,na+1),iaa(ma+1,na+1))
else
  allocate(a(1,1),ia(1,1))
  allocate(aa(1,1),iaa(1,1))
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
allocate(bb(mlocal+1,nlocal+1),ibb(mlocal+1,nlocal+1))

if (xmaster) then
  write(*,*)'scattertest with contiguous matrices:'
endif
if (xmaster) then
  aa(1:ma,1:na) = a
  iaa(1:ma,1:na) = ia
endif
call scattertest(a,ia,b,ib)
if (xmaster) then
  write(*,*)'scattertest with non-contiguous matrices:'
endif
call scattertest(aa(1:ma,1:na),iaa(1:ma,1:na),&
                 bb(1:mlocal,1:nlocal),ibb(1:mlocal,1:nlocal))

if (xmaster) then
  write(*,*)'colltest with contiguous matrices:'
endif
call colltest(a,b)
if (xmaster) then
  write(*,*)'colltest with non-contiguous matrices:'
endif
call colltest(aa(1:ma,1:na),bb(1:mlocal,1:nlocal))

if (xmaster) then
  write(*,*)'shifttest with contiguous matrices:'
endif
call shifttest(a,b)

if (xmaster) then
  write(*,*)'test with non-contiguous matrices:'
endif
call shifttest(aa(1:ma,1:na),bb(1:mlocal,1:nlocal))
! test vector_distr_send

! allocate vector-to-send
if (xmpi_master .eq. xmpi_rank) then
  allocate (vs(ma),vr(mlocal))
  allocate (vvs(2*ma),vvr(2*mlocal))
else
  allocate(vs(1),vr(mlocal))
  allocate(vvs(2),vvr(2*mlocal))
endif

if (xmaster) then
  write(*,*)'vector scattertest with contiguous vectors:'
endif
call vscattertest(vs,vr)

if (xmaster) then
  write(*,*)'vector scattertest with non-contiguous vectors:'
endif
call vscattertest(vvs(1:2*ma:2),vvr(1:2*mlocal:2))

allocate(x(nlocal))
allocate(xx(2*nlocal))

if (xmaster) then
  write(*,*)'getrowtest with contiguous matrices/vectors:'
endif
call getrowtest(b,x)

if (xmaster) then
  write(*,*)'getrowtest with non-contiguous matrices/vectors:'
endif
call getrowtest(bb(1:mlocal,1:nlocal),xx(1:2*nlocal:2))

allocate (b3(mlocal,nlocal,10),ib3(mlocal,nlocal,10))
if (xmaster) then
  write(*,*)'xshiftest with contiguous matrices/vectors:'
endif

call xshiftest(b,ib,b3,ib3)

allocate (bb3(mlocal*2,nlocal+1,20),ibb3(mlocal*2,nlocal+1,20))
if (xmaster) then
  write(*,*)'xshiftest with non-contiguous matrices/vectors:'
endif

call xshiftest(bb(1:mlocal,1:nlocal),ibb(1:mlocal,1:nlocal),&
               bb3(1:2*mlocal:2,1:nlocal,1:20:2),ibb3(1:2*mlocal:2,1:nlocal,1:20:2))

if (xmpi_master .eq. xmpi_rank) then
  write(*,*)xmpi_rank,': calling xmpi_finalize:'
endif
call xmpi_finalize

end program testgenmodule
