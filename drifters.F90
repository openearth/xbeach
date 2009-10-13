module drifter_module
contains
subroutine drifter(s,par)
    use params
    use spaceparams
    use readkey_module
    use xmpi_module
    use mnemmodule

    IMPLICIT NONE

    type(spacepars), target     :: s
    type(parameters)            :: par

    integer                     :: i
    integer                     :: j
    integer                     :: reclen,wordsize
    integer, save               :: it_drifter=0
    integer, save               :: ndrifter

    integer, dimension(:)  ,allocatable,save  :: iudrift,judrift,ivdrift,jvdrift
    real*8 , dimension(:)  ,allocatable,save  :: xdrift,ydrift,xwdrift,ywdrift,releasetime,retrievaltime
    logical, save                             :: first_drifter=.true.
    character*14                              :: fname
    character*78                              :: drifterfile
    real*8                                    :: facux,facuy,facvx,facvy,udrift,vdrift

    include 's.ind'
    include 's.inp'

if (first_drifter) then
   first_drifter=.false.
   ! check how many drifters are present
   if (xmaster) then
      ndrifter   = readkey_int     ('params.txt','ndrifter',    0,         0,        50)
   endif
#ifdef USEMPI
       call xmpi_bcast(ndrifter)
#endif
   if (ndrifter>0) then
      allocate(iudrift(ndrifter))
      allocate(judrift(ndrifter))
      allocate(ivdrift(ndrifter))
      allocate(jvdrift(ndrifter))
      allocate(xdrift(ndrifter))
      allocate(ydrift(ndrifter))
      allocate(xwdrift(ndrifter))
      allocate(ywdrift(ndrifter))
      allocate(releasetime(ndrifter))
      allocate(retrievaltime(ndrifter))
      if (xmaster) then
         call readkey('params.txt','drifterfile',drifterfile)
         open(10,file=drifterfile)
         do i=1,ndrifter
            read(10,*)xwdrift(i),ywdrift(i),releasetime(i),retrievaltime(i)
         enddo
         close(10)
      endif
#ifdef USEMPI
      call xmpi_bcast(xwdrift)
      call xmpi_bcast(ywdrift)
      call xmpi_bcast(releasetime)
      call xmpi_bcast(retrievaltime)
#endif
      if (xmaster) then
         do i=1,ndrifter
            xdrift(i)= cos(alfa)*(xwdrift(i)-xori)+sin(alfa)*(ywdrift(i)-yori)
            ydrift(i)=-sin(alfa)*(xwdrift(i)-xori)+cos(alfa)*(ywdrift(i)-yori)
         enddo
         inquire(iolength=wordsize) 1.d0
         reclen=wordsize*3
         do i=1,ndrifter
            write(fname(7:10),'(i4)')i+1000
            fname(1:7)='drifter'
            fname(11:14)='.dat'
            open(700+i,file=fname,form='unformatted',access='direct',recl=reclen)
         enddo
      endif
   endif
else
   if (ndrifter>0) then

      do i=1,ndrifter
         if (par%t>releasetime(i).and.par%t<retrievaltime(i)) then
!        
!        |     |     |
!        
!        +  -  +  -  +  -
!             o 
!  yv(jv)|     |     |
!        
!  y(ju) +  -  +  -  +  -
!         xu(iu)
!      x(iv)
            ! Only update position if within domain
            if (   xdrift(i)>xu(1) .and. xdrift(i)<=x(nx+1,1) .and.  &
       &           ydrift(i)>yv(1) .and. ydrift(i)<=y(1,ny+1) ) then 
               call hunt(xu,nx+1,xdrift(i),iudrift(i))
               call hunt(y(1,:),ny+1,ydrift(i),judrift(i))
               call hunt(x(:,1),nx+1,xdrift(i),ivdrift(i))
               call hunt(yv,ny+1,ydrift(i),jvdrift(i))
               facux=(xu(iudrift(i)+1)-xdrift(i))/(xu(iudrift(i)+1)-xu(iudrift(i)))
               facuy=(y(1,judrift(i)+1)-ydrift(i))/(y(1,judrift(i)+1)-y(1,judrift(i)))
               facvx=(x(ivdrift(i)+1,1)-xdrift(i))/(x(ivdrift(i)+1,1)-x(ivdrift(i),1))
               facvy=(yv(jvdrift(i)+1)-ydrift(i))/(yv(jvdrift(i)+1)-yv(jvdrift(i)))
               udrift=     facuy *(facux*uu(iudrift(i),judrift(i)  )+(1.d0-facux)*uu(iudrift(i)+1,judrift(i)  )) &
                  & +(1.d0-facuy)*(facux*uu(iudrift(i),judrift(i)+1)+(1.d0-facux)*uu(iudrift(i)+1,judrift(i)+1))
               vdrift=     facvx *(facvy*vv(ivdrift(i)  ,jvdrift(i))+(1.d0-facvy)*vv(ivdrift(i)  ,jvdrift(i)+1)) &
                  & +(1.d0-facvx)*(facvy*vv(ivdrift(i)+1,jvdrift(i))+(1.d0-facvy)*vv(ivdrift(i)+1,jvdrift(i)+1))
               xdrift(i)=xdrift(i)+udrift*par%dt
               ydrift(i)=ydrift(i)+vdrift*par%dt
#ifdef USEMPI
            else     ! In case of MPI set drifter coordinates to huge if outside domain
                     ! This allows the right coordinates to be communicated by allreduce
               xdrift(i)=huge(0.0d0)
               ydrift(i)=huge(0.0d0)  
#endif
            endif
#ifdef USEMPI

               call xmpi_allreduce(xdrift,MPI_MIN)
               call xmpi_allreduce(ydrift,MPI_MIN)
#endif            
         endif
      enddo
   endif
   if (xmaster) then
      if (abs(mod(par%t,par%tintp))<1.d-6) then
         it_drifter=it_drifter+1
         do i=1,ndrifter
            if (par%t>releasetime(i).and.par%t<retrievaltime(i)) then
               xwdrift(i)=xori + cos(alfa)*xdrift(i) - sin(alfa)*ydrift(i)
               ywdrift(i)=yori + sin(alfa)*xdrift(i) + cos(alfa)*ydrift(i)
            else
               xwdrift(i)=-999.d0
               ywdrift(i)=-999.d0
            endif
            write(700+i,rec=it_drifter) xwdrift(i),ywdrift(i),par%t
         enddo
      endif 
   endif     
endif


end subroutine drifter










subroutine hunt(xx        ,n         ,x         ,jlo       )
!
    implicit none
!
! Global variables
!
    integer                                  :: jlo
    integer                    , intent(in)  :: n
    real*8                     , intent(in)  :: x
    real*8       , dimension(n), intent(in)  :: xx
!
! Local variables
!
    integer :: inc
    integer :: jhi
    integer :: jm
    logical :: ascnd
!
!! executable statements -------------------------------------------------------
!
    ascnd = xx(n)>=xx(1)
    if (jlo<=0 .or. jlo>n) then
       jlo = 0
       jhi = n + 1
       goto 3
    endif
    inc = 1
    if (x>=xx(jlo) .eqv. ascnd) then
    1  continue
       jhi = jlo + inc
       if (jhi>n) then
          jhi = n + 1
       elseif (x>=xx(jhi) .eqv. ascnd) then
          jlo = jhi
          inc = inc + inc
          goto 1
       else
       endif
    else
       jhi = jlo
    2  continue
       jlo = jhi - inc
       if (jlo<1) then
          jlo = 0
       elseif (x<xx(jlo) .eqv. ascnd) then
          jhi = jlo
          inc = inc + inc
          goto 2
       else
       endif
    endif
    3 continue
    if (jhi - jlo==1) then
       return
    endif
    jm = (jhi + jlo)/2
    if (x>xx(jm) .eqv. ascnd) then
       jlo = jm
    else
       jhi = jm
    endif
    goto 3
end subroutine hunt
end module drifter_module
