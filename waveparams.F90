module waveparams
implicit none
type waveparameters

integer                                 :: K, Npy, Nr
integer, dimension(:), pointer          :: index_vector

real*8                                  :: mainang,dang
real*8                                  :: hm0gew, df
real*8                                  :: Ly, dt, rt
real*8,dimension(:),pointer             :: S0, dthetafin, fgen, theta0
real*8,dimension(:),pointer             :: Sf, Dd, f, theta, window
real*8,dimension(:,:),pointer           :: S_array

double complex,dimension(:,:),pointer   :: CompFn

end type

contains

! --------------------------------------------------------------
! ------------- Sorting for calling functions ------------------
! --------------------------------------------------------------
subroutine makebcf(par,s,wp)

use params
use spaceparams
use readkey_module
use xmpi_module

IMPLICIT NONE

! Input / output variables
type(parameters), INTENT(INOUT)             :: par
type(spacepars), intent(IN)                 :: s

type(waveparameters)                        :: wp
real*8                                      :: h0t0
real*8,save                                 :: bcendtime
character*80                                :: fname,Ebcfname,qbcfname
character*8                                 :: testc
logical                                     :: makefile
integer,save                                :: reuse  ! = 0 used to be in code


makefile=.false.
! First time tests
if (abs(par%t-par%dt)<1.d-6) then
    bcendtime=0
    par%listline=0
    if(xmaster) then
      call readkey('params.txt','bcfile',fname)
	  call checkbcfilelength(par,fname)
      open(74,file=fname,form='formatted')
      if (par%instat/=41) read(74,*)testc
    endif
    if(xmaster) then
      open(53,file='ebcflist.bcf',form='formatted',status='replace')  ! Make new files, don't add to existing files
      open(54,file='qbcflist.bcf',form='formatted',status='replace')
      close(53)
      close(54)
    endif
    if (testc=='FILELIST' .or. par%instat==41) then     ! New BCF-files are made when needed
        reuse=0
    else
        reuse=1                     ! Same file is reused
        if(xmaster) then
          close(74)    ! only close if this is not the list of files
        endif
    end if
end if

if (par%t>=(par%tstop-par%dt)) then
    return                          ! Hard code to prevent recalculation of bc for last timestep
end if                              ! par%listline is not increased, therfore, first line of current bcf is used (i.e. all zeros)


! Lookup rt and dt for bcf-files and create names for E and q bcf-files
if (reuse==0) then
    if(xmaster) then
      if (par%instat/=41) then
         read(74,*)wp%rt,wp%dt,fname
		 wp%rt = wp%rt / max(par%morfac,1.d0)
      endif
    endif
#ifdef USEMPI
    !Dano call xmpi_bcast(wp%rt)
    !Dano call xmpi_bcast(wp%dt)
    !Dano call xmpi_bcast(fname)
#endif
    Ebcfname='E_'//fname
    qbcfname='q_'//fname
else 
    wp%rt = readkey_dbl ('params.txt','rt'      , 3600.d0, 1200.d0,7200.d0,bcast=.false.)
    wp%rt = wp%rt / max(par%morfac,1.d0)
    wp%dt = readkey_dbl ('params.txt','dtbc', 0.1d0,0.01d0,1.0d0,bcast=.false.)
    Ebcfname='E_reuse.bcf'
    qbcfname='q_reuse.bcf'
end if

! Check to (re)make bcf files only if t=0 or (t>0 and reuse=0)
if (abs(par%t-par%dt)<1.d-6) then
    makefile=.true.
else
    if (reuse==0) then
        makefile=.true.
    end if
end if

! This code is used to generate bcf-files, only if makefile==.true.
if (makefile) then
    h0t0=sum(s%hh(1,:))/(s%ny+1)
    if (par%instat==4.or.par%instat==41) then
        call build_jonswap(par,s,wp,fname)
        call build_etdir(par,s,wp,h0t0,Ebcfname)
        call build_boundw(par,s,wp,h0t0,qbcfname)
    elseif (par%instat==5) then
        call swanreader(par,s,wp,fname)
        call build_etdir(par,s,wp,h0t0,Ebcfname)
        call build_boundw(par,s,wp,h0t0,qbcfname)
    elseif (par%instat==6) then
        call vardensreader(par,s,wp,fname)
        call build_etdir(par,s,wp,h0t0,Ebcfname)
        call build_boundw(par,s,wp,h0t0,qbcfname)
    endif
end if

! Keep track of which line in ebcflist.bcf and qbcflist.bcf should be read in boudaryconditions.f90
! Calculate when boudaryconditions.f90 should start making new bcf-files
par%listline=par%listline+1
bcendtime=bcendtime+wp%rt

if(xmaster) then
  ! Keep index/list of bcf file names and time information (ebcflist.bcf and qbcflist.bcf).
  open(53,file='ebcflist.bcf',form='formatted',position='append')
  open(54,file='qbcflist.bcf',form='formatted',position='append')
  write(53,'(f12.3,a,f12.3,a,f9.3,a,f9.5,a,f11.5,a)') & 
   & bcendtime,' ',wp%rt,' ',wp%dt,' ',par%Trep,' ',wp%mainang,' '//trim(Ebcfname)
  write(54,'(f12.3,a,f12.3,a,f9.3,a,f9.5,a,f11.5,a)') &
   & bcendtime,' ',wp%rt,' ',wp%dt,' ',par%Trep,' ',wp%mainang,' '//trim(qbcfname)
  close(53)
  close(54)
endif

end subroutine makebcf


! --------------------------------------------------------------
! ----------------------Read spectrum files --------------------
! --------------------------------------------------------------
subroutine build_jonswap(par,s,wp,fname)

use readkey_module
use params
use spaceparams
use xmpi_module

IMPLICIT NONE

! Input / output variables
type(parameters), INTENT(INout)         :: par
type(spacepars), intent(IN)             :: s
type(waveparameters), INTENT(INOUT)     :: wp
character(len=*)                        :: fname

! Internal variables
integer                                 :: i=0,ii,nang,nfreq,ier
integer                                 :: firstp, lastp
real*8,dimension(:),allocatable         :: temp, x, y
real*8                                  :: t1, dfj, fnyq, fp
real*8                                  :: gam, scoeff
character(len=80)                       :: dummystring

! Start program
if(xmaster) then
  if (par%instat /= 41) then
     write(*,*)'waveparams: Reading from ',fname,' ...'
  else
     call readkey('params.txt','bcfile',fname)
     write(*,*)'waveparams: reading from table',fname,' ...'
     read(74,*,iostat=ier)wp%hm0gew,fp,wp%mainang,gam,scoeff,wp%rt,wp%dt
	 wp%rt = wp%rt/max(par%morfac,1.d0)
	 fp=1.d0/fp
	 fnyq = 3.d0*fp
	 dfj= fp/20.d0
  endif
endif
! Read JONSWAP file 
!                             Input file  Keyword      Default       Minimum     Maximum   
if (par%instat/=41) then               
   wp%hm0gew           = readkey_dbl (fname,'Hm0'     ,   0.0d0,       0.00d0,       5.0d0,bcast=.false.)
   fp                  = readkey_dbl (fname,'fp'      ,   0.08d0,      0.0625d0,     0.4d0,bcast=.false.)
   fnyq                = readkey_dbl (fname,'fnyq'    ,   0.3d0,       0.2d0,        1.0d0,bcast=.false.)
   dfj                 = readkey_dbl (fname,'dfj'     ,   fnyq/200,    fnyq/1000,  fnyq/20,bcast=.false.)
   gam                 = readkey_dbl (fname,'gammajsp',   3.3d0,       1.0d0,        5.0d0,bcast=.false.)
   scoeff              = readkey_dbl (fname,'s'       ,   10.0d0,      1.0d0,     1000.0d0,bcast=.false.)
   wp%mainang          = readkey_dbl (fname,'mainang' ,   270.0d0,     0.0d0,      360.0d0,bcast=.false.)
   if(xmaster) then
     call readkey(fname,'checkparams',dummystring)
   endif
endif


wp%Npy=s%ny+1
write(*,*)'Hm0 = ',wp%hm0gew,'Tp = ',1.d0/fp,'dir = ',wp%mainang,'duration = ',wp%rt
!par%Trep=0.8345d0*(1/fp)                      ! approximation from Coastal Engineering: Processes, Theory and Design Practice
                                            ! Dominic Reeve, Andrew Chadwick 2004
allocate(temp(ceiling((fnyq-dfj)/dfj)))
temp=(/(i,i=1,size(temp))/)

allocate(wp%f(size(temp)))
wp%f=temp*dfj
deallocate (temp)

allocate(x(size(wp%f)))
x=wp%f/fp

allocate(y(size(wp%f)))
call jonswapgk(x,gam,y)


y=(wp%hm0gew/(4.d0*sqrt(sum(y)*dfj)))**2*y
deallocate (x)
        
t1=-(par%px)/2.d0

allocate(temp(101))
allocate(wp%theta(101))
temp=(/(i,i=0,100)/)
wp%theta=temp*((par%px)/100.d0)+t1
deallocate (temp)

wp%dang=wp%theta(2)-wp%theta(1)


allocate (wp%Dd(size(wp%theta)))
wp%mainang=(1.5d0*par%px-s%alfa)-wp%mainang*atan(1.d0)/45.0d0
if (wp%mainang>2.d0*par%px) then
    wp%mainang=wp%mainang-2.d0*par%px
elseif (wp%mainang<-2.d0*par%px) then
    wp%mainang=wp%mainang+2.d0*par%px
endif
wp%Dd = cos((wp%theta-wp%mainang)/2.d0)**(2.d0*scoeff)
wp%Dd = wp%Dd / (sum(wp%Dd)*wp%dang)


nang=size(wp%theta)
nfreq=size(y)
allocate(wp%S_array(nfreq,nang))

do i=1,nang                             ! Fill S_array
    do ii=1,nfreq
        wp%S_array(ii,i)=y(ii)*wp%Dd(i)
    end do  
end do
deallocate (y)
                                        ! integrate S_array over angles
allocate(wp%Sf(size(wp%f)))
wp%Sf = sum(wp%S_array, DIM = 2)*wp%dang 
                                        ! dimension [length f]

call tpDcalc(wp%Sf,wp%f,par%Trep)
par%Trep=1.d0/par%Trep

allocate(temp(size(wp%Sf)))
temp=0.d0
call frange(par, wp%Sf,firstp,lastp,temp)
deallocate (temp)

!!!!! ja/ap wp%K=max(100,ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1))

wp%K=ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1)  !!! this has changed

return

end subroutine build_jonswap
!
! --------------------------------------------------------------------
!
subroutine swanreader(par,s,wp,fname)

use params
use spaceparams
use readkey_module
use math_tools
use xmpi_module

IMPLICIT NONE

! Input / output variables
type(parameters), INTENT(INout)         :: par
type(spacepars), intent(IN)             :: s
type(waveparameters), INTENT(INOUT)     :: wp
character(len=*), INTENT(IN)            :: fname

! Internal variables
character(6)                            :: rtext
real*8                                  :: factor,exc,m0,dthetaS_XB 
integer                                 :: nfreq, ndir, switch, i, flipped
integer                                 :: firstp,lastp,nt,Ashift
real*8, dimension(:),allocatable        :: temp, findline
real*8, dimension(:,:),allocatable      :: tempA

dthetaS_XB = readkey_dbl ('params.txt','dthetaS_XB', 0.0d0, -360.0d0, 360.0d0,bcast=.false.)

flipped=0
wp%Npy=s%ny+1

switch = 0
if(xmaster) then
  write(*,*)'Reading from SWAN file ',fname,' ...'
  open(44,file=fname,form='formatted',status='old')


! Read file until RFREQ or AFREQ is found

  do while (switch==0)
      read(44,'(a)')rtext
      if (rtext == 'RFREQ ') then
          switch = 1
      elseif (rtext == 'AFREQ ') then
          switch = 2
      end if
  end do

! Read nfreq and f

  read(44,*)nfreq
endif
#ifdef USEMPI
!Dano call xmpi_bcast(nfreq)
!Dano call xmpi_bcast(switch)
#endif
allocate(wp%f(nfreq))
if(xmaster) then
  do i=1,nfreq
      read(44,*)wp%f(i)
  end do
endif
#ifdef USEMPI
!Dano call xmpi_bcast(wp%f)
#endif

! Convert to absolute frequencies

if (switch == 1) then
    wp%f = wp%f
else 
    wp%f = wp%f
end if

! Read CDIR or NDIR

if(xmaster) then
  read(44,'(a)')rtext
  if (rtext == 'NDIR  ') then
      switch = 1
  elseif (rtext == 'CDIR  ') then
      switch = 2
  else
      write(*,*) 'SWAN directional bins keyword not found'
#ifdef USEMPI
      call xmpi_abort
#else
      stop
#endif
  end if

  ! Read ndir, theta

  read(44,*)ndir
endif
#ifdef USEMPI
!Dano call xmpi_bcast(ndir)
#endif
allocate(wp%theta(ndir))

if(xmaster) then
do i=1,ndir
    read(44,*)wp%theta(i)
end do
endif

#ifdef USEMPI
!Dano call xmpi_bcast(wp%theta)
#endif

! Convert angles to XBeach angles and radians

if (switch == 1) then
    wp%theta = 270-s%alfa-wp%theta
else
    wp%theta = wp%theta-dthetaS_XB                  ! dthetaS_XB is the angle in the degrees to rotate the x-axis in SWAN to the
                                                    ! x-axis in XBeach (in Cartesian terms) (Have Fun :-))
end if

! Ensure angles are increasing instead of decreasing
if ((wp%theta(2)-wp%theta(1))<0) then
    call flipv(wp%theta,size(wp%theta))
    flipped=1
end if

nt = 0
Ashift = 0
! Make sure that all angles are in -180 to 180 degrees
if(minval(wp%theta)<-180)then
  allocate (temp(ndir))
  Ashift=-1
  temp=0
  do i=1,ndir
      if (wp%theta(i)<-180) then
          wp%theta(i)=wp%theta(i)+360.0d0
          nt = nt+1
      endif
  enddo
  !temp(1:ndir-nt) = wp%theta(nt+1:ndir)
  !temp(ndir-nt+1:ndir) = wp%theta(1:nt)
  temp(1:nt) = wp%theta(ndir-nt+1:ndir)
  temp(nt+1:ndir) = wp%theta(1:ndir-nt)
  wp%theta=temp
  deallocate(temp)
elseif(maxval(wp%theta)>180.0d0)then
  allocate(temp(ndir))
  Ashift=1.0d0
  temp=0
  do i=1,ndir
      if (wp%theta(i)>180.0d0) then
          wp%theta(i)=wp%theta(i)-360.0d0
          nt = nt+1
      endif
  enddo
  !temp(1:ndir-nt) = wp%theta(nt+1:ndir)
  !temp(ndir-nt+1:ndir) = wp%theta(1:nt)
  temp(1:nt) = wp%theta(ndir-nt+1:ndir)
  temp(nt+1:ndir) = wp%theta(1:ndir-nt)
  wp%theta=temp
  deallocate(temp)
endif

wp%theta=wp%theta*par%px/180
wp%dang=wp%theta(2)-wp%theta(1)

! Skip Quant, next line, read VaDens or EnDens
if(xmaster) then
  read(44,'(a)')rtext
  read(44,'(a)')rtext
  read(44,'(a)')rtext
  if (rtext == 'VaDens') then
      switch = 1
  elseif (rtext == 'EnDens') then
      switch = 2
  else
      write(*,*) 'SWAN VaDens/EnDens keyword not found'
#ifdef USEMPI
      call xmpi_abort
#else
      stop
#endif
  end if
endif
#ifdef USEMPI
!Dano call xmpi_bcast(switch)
#endif

if(xmaster) then
  read(44,'(a)')rtext
  read(44,*)exc
endif

#ifdef USEMPI
!Dano call xmpi_bcast(rtext)
!Dano call xmpi_bcast(exc)
#endif

if(xmaster) then
  i=0
  ! Find FACTOR keyword
  do while (i==0)
      read(44,'(a)')rtext
      if (rtext == 'FACTOR') then
          i=1
      elseif (rtext == 'ZERO  ') then
          write(*,*) 'Zero energy density input for this point'
#ifdef USEMPI
          call xmpi_abort
#else
          stop
#endif
      elseif (rtext == 'NODATA') then
          write(*,*) 'SWAN file has no data for this point'
#ifdef USEMPI
          call xmpi_abort
#else
          stop
#endif
      end if
  end do
  read(44,*)factor
endif

#ifdef USEMPI
!Dano call xmpi_bcast(factor)
#endif
! Read S_array
allocate(wp%S_array(nfreq,ndir))

if(xmaster) then
  do i=1,nfreq
      read(44,*)wp%S_array(i,:)
  end do
endif

#ifdef USEMPI
!Dano call xmpi_bcast(wp%S_array)
#endif

where (wp%S_array == exc) wp%S_array =0

! If angles were decreasing, flip S_array as also dir is flipped
if (flipped == 1) then
    flipped=2
    call flipa(wp%S_array,nfreq,ndir,flipped)
end if


! Make sure that all wave variance is between -180 to 180 degrees range
if(Ashift==-1)then
    allocate(tempA(nfreq,ndir))
    tempA=0
!   tempA(:,ndir-nt+1:ndir) = wp%S_array(:,1:nt)
!   tempA(:,1:ndir-nt) = wp%S_array(:,nt+1:ndir)
    tempA(:,1:nt) = wp%S_array(:,ndir-nt+1:ndir)
    tempA(:,nt+1:ndir) = wp%S_array(:,1:ndir-nt)
    wp%S_array=tempA
    deallocate(tempA)
elseif (Ashift==1) then
    allocate(tempA(nfreq,ndir))
    tempA=0
!   tempA(:,ndir-nt+1:ndir) = wp%S_array(:,1:nt)
!   tempA(:,1:ndir-nt) = wp%S_array(:,nt+1:ndir)
    tempA(:,1:nt) = wp%S_array(:,ndir-nt+1:ndir)
    tempA(:,nt+1:ndir) = wp%S_array(:,1:ndir-nt)
    wp%S_array=tempA
    deallocate(tempA)
endif

wp%S_array=wp%S_array*factor

if(xmaster) then
  close(44)                               ! Finished reading file
endif


! Convert to m2/Hz/rad

wp%S_array=wp%S_array*180/par%px

! Convert from energy density to variance density

if (switch == 2) then
    wp%S_array=wp%S_array/(par%rho*par%g)
end if

allocate(wp%Sf(nfreq))  
wp%Sf = sum(wp%S_array, DIM = 2)*wp%dang

! Find main wave direction
allocate (temp(ndir))
temp=sum(wp%S_array, DIM = 1)
i=maxval(maxloc(temp))
wp%mainang=wp%theta(i)
deallocate(temp)

allocate (temp(nfreq+1))
temp(1)=0
temp(2:nfreq)=0.5*wp%f(1:nfreq-1)+0.5*wp%f(2:nfreq)
temp(nfreq+1)=wp%f(nfreq)
! Calculate zero-order moment
m0=0
!m1=0
do i=1,nfreq
    m0=m0+wp%Sf(i)*(temp(i+1)-temp(i))
!    m1=m1+wp%f(i)*wp%Sf(i)*(temp(i+1)-temp(i))
end do
deallocate (temp)

wp%hm0gew=4*sqrt(m0)
!par%Trep=1/(m1/m0)

call tpDcalc(wp%Sf,wp%f,par%Trep)
par%Trep=1.d0/par%Trep

allocate(findline(size(wp%Sf)))

firstp=0
lastp=0
call frange(par, wp%Sf,firstp,lastp,findline)
deallocate(findline)


!!!!! ja/ap wp%K=max(100,ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1))

wp%K=ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1)  !!! this has changed


allocate(wp%Dd(ndir))

wp%Dd=sum(wp%S_array, DIM = 1)

end subroutine swanreader


!
! --------------------------------------------------------------------
!
subroutine vardensreader(par,s,wp,fname)

use params
use spaceparams
use xmpi_module

IMPLICIT NONE

! Input / output variables
type(parameters), INTENT(INout)         :: par
type(spacepars), intent(IN)             :: s
type(waveparameters), INTENT(INOUT)     :: wp
character(len=*), INTENT(IN)            :: fname

! Internal variables
real*8, external                        :: readkey_dbl
real*8                                  :: m0 
integer                                 :: nfreq, ndir,i
integer                                 :: firstp,lastp
real*8, dimension(:),allocatable        :: temp, findline

wp%Npy=s%ny+1

if(xmaster) then
  write(*,*)'Reading from VarDens file ',fname,' ...'
  open(44,file=fname,form='formatted',status='old')

  read(44,*)nfreq
endif
#ifdef USEMPI
!Dano call xmpi_bcast(nfreq)
#endif
allocate(wp%f(nfreq))

if(xmaster) then
  do i=1,nfreq
      read(44,*)wp%f(i)
  end do

  read(44,*)ndir
endif
#ifdef USEMPI
!Dano call xmpi_bcast(wp%f)
!Dano call xmpi_bcast(ndir)
#endif
allocate(wp%theta(ndir))

if(xmaster) then
  do i=1,ndir
      read(44,*)wp%theta(i)
  end do
endif
#ifdef USEMPI
!Dano call xmpi_bcast(wp%theta)
#endif

wp%theta=wp%theta*par%px/180
wp%dang=wp%theta(2)-wp%theta(1)

! Read S_array
allocate(wp%S_array(nfreq,ndir))

if(xmaster) then
  do i=1,nfreq
      read(44,*)wp%S_array(i,:)
  end do

  close(44)                               ! Finished reading file
endif

#ifdef USEMPI
!Dano call xmpi_bcast(wp%S_array)
#endif

! Convert to m2/Hz/rad

wp%S_array=wp%S_array*180/par%px

allocate(wp%Sf(nfreq))  
wp%Sf = sum(wp%S_array, DIM = 2)*wp%dang

allocate (temp(nfreq+1))
temp(1)=0
temp(2:nfreq)=0.5*wp%f(1:nfreq-1)+0.5*wp%f(2:nfreq)
temp(nfreq+1)=wp%f(nfreq)
! Calculate zero-order moment
m0=0.0d0
!m1=0.0d0
do i=1,nfreq
    m0=m0+wp%Sf(i)*(temp(i+1)-temp(i))
!    m1=m1+wp%f(i)*wp%Sf(i)*(temp(i+1)-temp(i))
end do
deallocate (temp)

! Find main wave direction
allocate (temp(ndir))
temp=sum(wp%S_array, DIM = 1)
i=maxval(maxloc(temp))
wp%mainang=wp%theta(i)
deallocate(temp)

wp%hm0gew=4*sqrt(m0)
!par%Trep=1/(m1/m0)
call tpDcalc(wp%Sf,wp%f,par%Trep)
par%Trep=1.d0/par%Trep
!
allocate(findline(size(wp%Sf)))
firstp=0
lastp=0
call frange(par, wp%Sf,firstp,lastp,findline)
deallocate(findline)

!!!!! ja/ap wp%K=max(100,ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1))

wp%K=ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1)  !!! this has changed



allocate(wp%Dd(ndir))

wp%Dd=sum(wp%S_array, DIM = 1)

end subroutine vardensreader


! --------------------------------------------------------------
! ---------------------- Build E_tdir file ---------------------
! --------------------------------------------------------------

subroutine build_etdir(par,s,wp,h,Ebcfname)

use params
use math_tools
use spaceparams
use interp
use xmpi_module

IMPLICIT NONE

! Input / output variables

type(parameters), INTENT(IN)            :: par
type(waveparameters), INTENT(INOUT)     :: wp
type(spacepars), INTENT(IN)             :: s
real*8, INTENT(IN)                      :: h
character(len=*), INTENT(IN)            :: Ebcfname
! Internal variables

! help integers
integer                                 :: Ns ! Number of theta bins

! counters
integer                                 :: i, ii, iii, stepf, stepang, index2
integer                                 :: firstp, lastp, M

! Nothings
integer                                 :: F2
real*8                                  :: pp

! help single variables with meaning
real*8                                  :: TT, kmax
real*8                                  :: hm0now, s1, s2, modf, modang

! help vectors
integer, dimension(wp%K)                :: Nbin
real*8,dimension(size(wp%Sf))           :: findline
real*8,dimension(size(wp%Dd))           :: Dmean, P
real*8,dimension(wp%K)                  :: P0, k,phase, Sf0, A, Sf0org,S0org
real*8,dimension(wp%K*2)                                :: randummy
real*8,dimension(:),allocatable         :: temp, temp2, t, Nbox
real*8,dimension(1:401)                 :: ktemp, ftemp

! help arrays
real*8,dimension(:,:), allocatable      :: D
real*8,dimension(:,:,:), allocatable    :: zeta, Ampzeta, E_tdir

! Complex help variables
! double complex                          :: ctemp
! wwvv double complex,dimension(:),allocatable  :: Gn, Comptemp
complex(fftkind),dimension(:),allocatable   :: Gn, Comptemp

! Check to see is Nyquest frequencydoes not interfere with wp%dt
pp=maxval(wp%f)*2.d0
if (wp%dt>(1.d0/pp)) then
    wp%dt=1.d0/pp
    if (xmaster) then
       write(*,'(a,f6.4,a)')'Changing dt in wave boundary conditions to satisfy Nyquist condition. New dt = ',wp%dt,' s.'
    endif
endif

! Start program
if(xmaster) then
  write(*,*)'Calculating wave energy at boundary'
endif

findline=0.0d0
call frange(par, wp%Sf,firstp,lastp,findline)
M=sum(findline)

allocate (temp(size(wp%Sf)))
temp=1.0d0                                  ! turn into ONES(length(findline),1)

allocate (D(size(findline),size(wp%Dd)))
D=0.0d0

do i=1,size(wp%Dd)                      ! Build D
    do ii=1,size(findline)
        D(ii,i)=temp(ii)*wp%Dd(i)
    end do
end do

deallocate (temp)

do i=1,size(wp%Dd)
    D(:,i)=D(:,i)*findline              ! Information in D is now destroyed, but not further necessary
end do

Dmean=sum(D, DIM=1)/M

allocate (temp(wp%K))                   ! choose K wave components in range around peak
temp=(/(i,i=0,wp%K-1)/)
                                        ! fgen = frequency array of generated waves
allocate(wp%fgen(wp%K))
wp%fgen=temp*((wp%f(lastp)-wp%f(firstp))/(wp%K-1))+wp%f(firstp)  
deallocate(temp)

wp%df=(wp%fgen(wp%K)-wp%fgen(1))/(dble(wp%K)-1.d0)
! Avoid leakage
!!!!wp%fgen=floor(wp%fgen/wp%df)*wp%df  !taken out because this gives double values of fgen (Pieter, Jaap and Ap 28/5)                  
TT=1/wp%df
    
! Make table of dispersion relation
kmax=((2*(par%px)*wp%f(lastp))**2)/par%g    
allocate(temp(401))
temp=(/(i,i=0,400)/)                    
ktemp=temp*(kmax/200)
ftemp=sqrt((par%g)*ktemp*tanh(ktemp*h))/(2*par%px)      
deallocate (temp)

do i=1,size(wp%fgen)
    call LINEAR_INTERP(ftemp, ktemp,401, wp%fgen(i),pp,F2)
    k(i)=pp
end do

pp=1/(sum(Dmean)*wp%dang)               ! Normalization facort for wave variance
do i=1,size(wp%theta)
    P(i)=sum(Dmean(1:i))*wp%dang*pp     ! Cummulative normalized wave variance in directional space used as probability density function
end do  

if (par%random==1) CALL RANDOM_SEED                        ! Call random seed
!call random_number(P0)
call random_number(randummy)
if (xmaster) then
!  open(555,file='tempout1')  ! wwvv todo: better to use the same unit number for
!                             ! short I/O sequences like this. 555 cab interfere
!                             ! with varoutput unit numbers
!  write(555,*)randummy
!  close(555)
endif
P0=randummy(1:wp%K)


P0=0.95*P0+0.05/2                       ! Pick K wave directions

allocate(wp%theta0(wp%K))
do i=1,size(P0)
    call LINEAR_INTERP(P(1:size(P)-1),wp%theta(1:size(P)-1),size(P)-1,P0(i),pp,F2)
    wp%theta0(i)=pp
end do

F2=nint(TT/wp%dt)
!if (mod(real(F2),2.)/=0) then
if (mod(F2,2)/=0) then
    F2=F2+1
end if

allocate(t(F2))
do i=1,F2
    t(i)=wp%dt*i                        ! time axis
end do

wp%Nr=nint(TT/wp%dt)    
!if (mod(real(wp%Nr),2.d0)/=0) then        ! Ensure even number
if (mod(wp%Nr,2)/=0) then        ! Ensure even number
    wp%Nr=wp%Nr+1
end if

!call random_number(phase)               ! Random phase
phase=randummy(wp%K+1:2*wp%K)

phase=2*phase*par%px

! Interp across main diagonal of S_array
allocate(wp%S0(wp%K))
do i=1,size(wp%fgen)
    ! Find in which 'box' interpolation point is
    ! In what frequency step is the point?
    ! -1 to move to row number instead of step number
    allocate(temp(size(wp%f)))
    allocate(temp2(size(wp%f)))
    temp2=(/(ii,ii=1,size(wp%f))/)
    temp=1
    where (wp%f < wp%fgen(i) ) temp=0
    temp=temp*temp2
    if (sum(temp)==0) then
        stepf=size(wp%f)-1
    else
        stepf=max(nint(minval(temp, MASK = temp .gt. 0) -1),1)
    end if
    modf=(wp%fgen(i)-wp%f(stepf))/(wp%f(stepf+1)-wp%f(stepf))
    deallocate(temp,temp2)
    ! In what angle step is the point?
    ! -1 to move to row number instead of step number
    allocate(temp(size(wp%theta)))
    allocate(temp2(size(wp%theta)))
    temp2=(/(ii,ii=1,size(wp%theta))/)
    temp=1
    where (wp%theta < wp%theta0(i) ) temp=0
    temp=temp*temp2
    if (wp%theta0(i)==wp%theta(1)) then
        stepang=1
    else
        stepang=nint(minval(temp, MASK = temp .gt. 0) -1)
    end if
    modang=(wp%theta0(i)-wp%theta(stepang))/(wp%theta(stepang+1)-wp%theta(stepang))
    deallocate(temp,temp2)
    ! Linear interpolation for two lines in theta direction
    s1=(1.d0-modang)*wp%S_array(stepf,stepang)+modang*wp%S_array(stepf,stepang+1)
    s2=(1.d0-modang)*wp%S_array(stepf+1,stepang)+modang*wp%S_array(stepf+1,stepang+1)
    
    ! Linear interpolation between s1 and s2
    wp%S0(i)=max(tiny(0.d0),(1.d0-modf)*s1+modf*s2)                     ! Robert, in case no energy is drawn
!       wp%S0(i)=(1.d0-modf)*s1+modf*s2
end do

do i=1,size(wp%fgen)
    call Linear_interp(wp%f,wp%Sf,size(wp%f),wp%fgen(i),pp,F2)
    Sf0(i)=pp
end do
    
! correction for Hm0gem 6/8/2001
hm0now = 4*sqrt(sum(Sf0)*wp%df)         ! not S0!!!


Sf0org = Sf0
S0org=wp%S0
!!!wp%S0 = (wp%hm0gew/hm0now)**2*wp%S0
!!!!Sf0 = (wp%hm0gew/hm0now)**2*Sf0

allocate(wp%dthetafin(wp%K))
wp%dthetafin = Sf0/wp%S0                ! is chosen such that the spreading function is one.                                                  


                         
! amplitude of components
A = sqrt(2*wp%S0*wp%df*wp%dthetafin)    

Sf0=Sf0org
wp%S0=S0org
! Fourier coefficients
allocate(wp%CompFn(wp%Npy,wp%Nr))                   
wp%CompFn=0.d0
! index_vector of Fourier components associated with fgen
allocate(wp%index_vector(wp%K))
wp%index_vector = floor(wp%f(firstp)/wp%df)+1+nint((wp%fgen-wp%f(firstp))/wp%df)

do i=1,wp%K
    wp%CompFn(1,wp%index_vector(i)) = A(i)/2*exp(par%compi*phase(i))
enddo

allocate(Comptemp(size(wp%CompFn(1,wp%Nr/2+2:wp%Nr))))
Comptemp = conjg(wp%CompFn(1,2:wp%Nr/2))
call flipiv(Comptemp,size(Comptemp))    ! Flip left-right routine
wp%CompFn(1,wp%Nr/2+2:wp%Nr)=Comptemp

do index2=2,wp%Npy                      ! construct Fourier components for every y position
    wp%CompFn(index2,wp%index_vector)=wp%CompFn(1,wp%index_vector) & 
                            *exp(-par%compi*k*sin(wp%theta0)*s%yz(index2))
    Comptemp = conjg(wp%CompFn(index2,2:wp%Nr/2))
    call flipiv(Comptemp,size(Comptemp))
    wp%CompFn(index2,wp%Nr/2+2:wp%Nr)=Comptemp
end do
deallocate(Comptemp)
Ns=s%ntheta
allocate(temp(Ns+1))
temp=(/(i,i=0,Ns)/)
! ensure all theta=themamax is included
temp(Ns+1)=temp(Ns+1)+epsilon(1.0)          
allocate(Nbox(Ns+1))
Nbox=(par%thetamin*par%px/180)+temp*(par%dtheta*par%px/180)
deallocate (temp)

! Histc function on theta0 with Nbox edges
do i=1,size(wp%theta0)
    ! Relate each theta0 to a bin number
    Nbin(i)=ceiling((wp%theta0(i)-Nbox(1))/(par%dtheta*par%px/180.d0))
    if (mod((wp%theta0(i)-Nbox(1)),(par%dtheta*par%px/180.d0))==0) then
            Nbin(i)=Nbin(i)+1           ! To ensure binlow<=x<binhigh
    end if
end do
i=(maxval(Nbin))
!! Change theta0 direction to middle of wave energy bin : Robert
if (par%nspr==1) then
        do i=1,wp%K
                ! If component is outside wave distribution bins move to outer bins (so include energy)
                                if (Nbin(i)<=0) then 
                                        Nbin(i)=1
                                elseif (Nbin(i)>Ns) then
                                        Nbin(i)=Ns
                                endif
                                wp%theta0(i)=s%theta(Nbin(i))
        enddo
endif


deallocate(Nbox)

allocate(wp%window(size(t)))
allocate(temp(size(t)))
temp=t
where(t>wp%rt)temp=0        ! ensure window makes 0 at t=0 and at t=rt
wp%window=1  
! fc = 2*96
wp%window=wp%window*(tanh(192.d0*temp/maxval(temp))**2)*(tanh(192.d0*(1.d0-temp/maxval(temp)))**2)
deallocate(temp)

allocate(zeta(wp%Npy,wp%Nr,Ns))
allocate(Ampzeta(wp%Npy,wp%Nr,Ns))
zeta=0
Ampzeta=0

! Nr = time, wp%Npy = y-axis points, Ns= theta bins
do ii=1,Ns
    if(xmaster) then
      write(*,'(A,I0,A,I0)')'Calculating wave energy for theta bin ',ii,' of ',Ns
    endif
    do index2=1,wp%Npy
            
        allocate(Gn(wp%Nr))
        Gn=0
        allocate(temp(size(Nbin)))
        temp=0
        where (Nbin==ii)
            temp=1.
        end where
        
        F2=nint(sum(temp))
                
        if (F2/=0) then                 ! Check for 0-length vectors

            allocate(temp2(F2))
            temp2=0
            
            do i=1,F2
                iii=maxval(maxloc(temp))
                temp(iii)=0
                temp2(i)=wp%index_vector(iii)
            end do
            
            
            Gn(int(temp2))=wp%CompFn(index2,int(temp2))
            deallocate(temp2)
            
            allocate(Comptemp(size(Gn(wp%Nr/2+2:wp%Nr))))
            Comptemp = conjg(Gn(2:wp%Nr/2))
                                
            call flipiv(Comptemp,size(Comptemp))

            Gn(wp%Nr/2+2:wp%Nr)=Comptemp
            deallocate(Comptemp)

            allocate(Comptemp(size(Gn)))
            Comptemp=Gn
            F2=0
            Comptemp=fft(Comptemp,inv=.true.,stat=F2)
            
            ! Scale factor
            Comptemp=Comptemp/sqrt(real(size(Comptemp)))    

            zeta(index2,:,ii)=dble(Comptemp*wp%Nr)*wp%window
                        
            Comptemp=zeta(index2,:,ii)

            call hilbert(Comptemp,size(Comptemp))
        
            Ampzeta(index2,:,ii)=abs(Comptemp)
            
            if(xmaster) then
              if (F2/=0) then
                  write(*,'(A,I0,A,I0,A,I0)')'Y-point ',index2,' of ',wp%Npy,' done. Error code: ',F2
              else
                  write(*,'(A,I0,A,I0,A)')'Y-point ',index2,' of ',wp%Npy,' done.'
              end if          
            endif
                    
            deallocate(Comptemp)
        
        else
            if(xmaster) then
              write(*,'(A,I0,A)')'Theta bin ',ii,' empty at this point. Continuing to next point'
            endif

        end if

        deallocate(temp)
        deallocate(Gn)
    end do
end do

if(xmaster) then
  write(*,fmt='(a)')'writing wave energy to ',trim(Ebcfname),' ...'
endif
allocate(E_tdir(wp%Npy,wp%Nr,Ns))
E_tdir=0.0d0
E_tdir= 0.5d0*(par%rho)*(par%g)*Ampzeta**2
E_tdir=E_tdir/s%dtheta

if(xmaster) then
  open(12,file=Ebcfname,form='unformatted')
!  open(12,file=Ebcfname,form='binary')
  do i=1,wp%Nr
      write(12)E_tdir(:,i,:)
  end do
  do i=1,4                ! Ensure the file is always full to the end
      write(12)0.d0*E_tdir(:,1,:)
  end do
  close(12)
  write(*,*)'file done'
endif

deallocate (D,t,zeta,Ampzeta,E_tdir)

return

end subroutine build_etdir

! --------------------------------------------------------------
! ----------------------- Bound long wave ----------------------
! --------------------------------------------------------------
subroutine build_boundw(par,s,wp,h,qbcfname)

use params
use spaceparams
use math_tools
use xmpi_module

IMPLICIT NONE


! Input / output variables

type(parameters), INTENT(IN)            :: par
type(spacepars), INTENT(IN)             :: s
type(waveparameters), INTENT(INOUT)     :: wp
real*8, INTENT(IN)                      :: h
character(len=*), INTENT(IN)            :: qbcfname

! internal variables
integer                                 :: K, m, index1, Npy, Nr, i=0, jj
integer,dimension(:),allocatable        :: index2

real*8                                  :: g
real*8                                  :: df, deltaf
real*8,dimension(:), allocatable        :: w1, k1
real*8,dimension(:), allocatable        :: term1, term2, chk1, chk2
real*8,dimension(:,:),allocatable       :: Eforc, D, deltheta, KKx, KKy, theta3
real*8,dimension(:,:),allocatable       :: dphi3, k3, cg3, Abnd, qx

double complex,dimension(:),allocatable     :: Comptemp, Comptemp2
! wwvv double complex,dimension(:,:),allocatable    :: Gn, Ftemp, Ftemp2
complex(fftkind),dimension(:,:),allocatable :: Gn, Ftemp, Ftemp2

g=par%g
K=wp%K
df=wp%df
index1=wp%index_vector(1)
Npy=wp%Npy
Nr=wp%Nr

if(xmaster) then
  write(*,*) 'Calculating flux at boundary'
endif

allocate(Eforc(K-1,K),D(K-1,K),deltheta(K-1,K),KKx(K-1,K),KKy(K-1,K))
allocate(dphi3(K-1,K),k3(K-1,K),cg3(K-1,K))

! Set to zeros
Eforc = 0                               ! Herbers et al. (1994) eq. 1
                                        ! rows = diff freq, 
                                        ! columns is interaction

D = 0                                   ! Herbers eq. A5.
                                                                                
deltheta = 0                            ! difference angle between two primary wave components

KKx = 0                                 ! x component of wave number of difference wave
KKy = 0                                 ! y component

dphi3 = 0                               ! phase shift of diff wave

k3 = 0                                  ! wavenumber of diff wave

cg3 = 0                                 ! alternative def of group speed

allocate(w1(size(wp%fgen)),k1(size(wp%fgen)))
w1=0
k1=0

do m=1,K-1
    
    deltaf=m*df                         ! difference frequency

    w1=2*par%px*wp%fgen                 ! wave numbers of primary waves

    call bc_disper(k1,w1,size(w1),h,g)  ! wave numbers of primary waves

    deltheta(m,1:K-m) = abs(wp%theta0(m+1:K)-wp%theta0(1:K-m))+par%px

    KKy(m,1:K-m)=k1(m+1:K)*sin(wp%theta0(m+1:K))-k1(1:K-m)*sin(wp%theta0(1:K-m))

    KKx(m,1:K-m)=k1(m+1:K)*cos(wp%theta0(m+1:K))-k1(1:K-m)*cos(wp%theta0(1:K-m))
    
    k3(m,1:K-m) =sqrt(k1(1:K-m)**2+k1(m+1:K)**2+2*k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m)))

    cg3(m,1:K-m)= 2.d0*par%px*deltaf/k3(m,1:K-m)

    ! build transferfunction D (Herbers eq. A5)
    allocate(term1(K-m),term2(K-m),chk1(K-m),chk2(K-m))
    
    term1 = (-w1(1:K-m))*w1(m+1:K)
    term2 = (-w1(1:K-m))+w1(m+1:K)
    chk1  = cosh(k1(1:K-m)*h)
    chk2  = cosh(k1(m+1:K)*h)

!previous version with bug:    D(m,1:K-m) = -g*k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m))/2/term1+g*term2*(chk1*chk2)/ &
!                ((g*k3(m,1:K-m)*tanh(k3(m,1:K-m)*h)-(term2)**2)*term1*cosh(k3(m,1:K-m)*h))* &
!                (term2*((term1)**2/g/g - k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m))) - &
!            0.5*((-w1(1:K-m))*k1(m+1:K)**2/(chk1**2)+w1(m+1:K)*k1(1:K-m)**2/(chk2**2)))

    D(m,1:K-m) = -g*k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m))/2.d0/term1+g*term2*(chk1*chk2)/ &
                ((g*k3(m,1:K-m)*tanh(k3(m,1:K-m)*h)-(term2)**2)*term1*cosh(k3(m,1:K-m)*h))* &
                (term2*((term1)**2/g/g - k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m))) - &
            0.50d0*((-w1(1:K-m))*k1(m+1:K)**2/(chk2**2)+w1(m+1:K)*k1(1:K-m)**2/(chk1**2)))


    deallocate(term1,term2,chk1,chk2)

    ! correction for surface elevation input and output instead of bottom pressure
    D(m,1:K-m) = D(m,1:K-m)*cosh(k3(m,1:K-m)*h)/(cosh(k1(1:K-m)*h)*cosh(k1(m+1:K)*h))  

    ! exclude interactions where f<deltaf (== lower limit Herbers eq. 1)    
    
    where(wp%fgen<=m*df) D(m,:)=0.d0

    if (m*df<=par%fcutoff) D(m,:)=0.d0

    ! S0 is energy density of the surface elevation, as measured
    ! by bottom pressure sensors in FRF array.  Herbers eq. 1.

    Eforc(m,1:K-m) = 2*D(m,1:K-m)**2*wp%S0(1:K-m)*wp%S0(m+1:K)* & 
                     wp%dthetafin(1:K-m)*wp%dthetafin(m+1:K)*df

    ! phase of longwave, so that the longwave is 180 degrees out of
    ! phase with the pair of primary waves.

    allocate(Comptemp(K-m),Comptemp2(K-m))
    
    Comptemp=conjg(wp%CompFn(1,index1+m:index1+K-1))
    Comptemp2=conjg(wp%CompFn(1,index1:index1+K-m-1))
    
    ! angle defined as imaginary part of teh natural log of a complex number
    ! angle1 = imag(log(Comptemp))
    ! angle2 = imag(log(Comptemp2))
    ! Works as long as the complex number is not zero
    dphi3(m,1:K-m) = par%px+imag(log(Comptemp))-imag(log(Comptemp2))
    deallocate (Comptemp,Comptemp2)

end do

allocate(theta3(K-1,K))
theta3 = atan(KKy/(KKx+0.0001d0))        ! angle of longwave

allocate(Gn(Npy,Nr))
allocate(Abnd(K-1,K))

Gn=0                                    ! Fourier coefficient of lo freq.
Abnd = sqrt(2*Eforc*df)                 ! Amplitude of low frequency

allocate(index2(K-1))
index2=(/(i,i=1,K-1)/)

allocate(Ftemp(K-1,K))                  ! everything in fluxes
                                        ! this is still of function of diff freq AND interaction pair
Ftemp = Abnd/2*exp(-1*par%compi*dphi3)*cg3*cos(theta3)
    
Gn(1,index2+1)=(sum(Ftemp,DIM=2))       ! sum over components

allocate (Comptemp(Nr/2-1))             ! deallocated later
Comptemp=conjg(Gn(1,2:Nr/2))
call flipiv(Comptemp,Nr/2-1)
Gn(1,Nr/2+2:Nr)=Comptemp

allocate(qx(Npy,Nr))
qx=0.0d0
allocate(Comptemp2(Nr)) 
if(xmaster) then
  write(*,'(A,I0)')'Flux 1 of ',Npy
endif
Comptemp2=fft(Gn(1,:),inv=.true.)
Comptemp2=Comptemp2/sqrt(real(Nr))
qx(1,:)=real(Comptemp2*Nr)*wp%window    ! flux as a function of time

! determine flux at every y position
allocate(Ftemp2(K-1,K))

do jj=2,Npy
    ! phase shift
    Ftemp2 = Ftemp*exp(-1*par%compi*KKy*s%yz(jj))
    Gn(jj,index2+1) = (sum(Ftemp2,DIM=2))
    Comptemp = conjg(Gn(jj,2:Nr/2))
    call flipiv(Comptemp,Nr/2-1)
    Gn(jj,Nr/2+2:Nr) = Comptemp

    if(xmaster) then
      write(*,'(A,I0,A,I0)')'Flux ',jj,' of ',Npy
    endif
    Comptemp2=fft(Gn(jj,:),inv=.true.)
    Comptemp2=Comptemp2/sqrt(real(Nr))
    qx(jj,:)=real(Comptemp2*Nr)*wp%window

end do

deallocate(Comptemp)
deallocate(Comptemp2)
deallocate(Ftemp)
deallocate(Ftemp2)
deallocate(index2)

if(xmaster) then
  write(*,fmt='(a)')'writing long wave mass flux to ',trim(qbcfname),' ...'
  open(21,file=qbcfname,form='unformatted')
!  open(21,file=qbcfname,form='binary')
  do i=1,wp%Nr
      write(21)qx(:,i)
  end do
  do i=1,4                ! Ensure the file is always full to the end
      write(21)0.d0*qx(:,1)
  end do
  close(21)
  write(*,*)'file done'
endif

! Clean up memory
deallocate(wp%index_vector,wp%S0,wp%dthetafin,wp%fgen,wp%theta0,wp%window)
deallocate(wp%Sf,wp%Dd,wp%f,wp%theta,wp%S_array,wp%CompFn)

end subroutine build_boundw




! -----------------------------------------------------------
! --------- JONSWAP  unscaled JONSWAP spectrum --------------
! ----------------(used by build_jonswap)--------------------
subroutine jonswapgk(x,gam,y)

IMPLICIT NONE
! Required input: - x           : nondimensional frequency, divided by the peak frequency
!                 - gam         : peak enhancement factor, optional parameter (DEFAULT 3.3)
!                 - y is output : nondimensional relative spectral density, equal to one at the peak

real*8, INTENT(IN)                  :: gam
real*8,dimension(:), INTENT(IN)     :: x
real*8,dimension(:), INTENT(INOUT)  :: y

! Internal variables
real*8,dimension(size(x))           :: xa, sigma, fac1, fac2, fac3, temp

xa=abs(x)

where (xa==0) 
    xa=1e-20
end where

sigma=xa

where (sigma<1.) 
    sigma=0.07
end where

where (sigma>=1.) 
    sigma=0.09
end where

temp=0*xa+1

fac1=xa**(-5)
fac2=exp(-1.25*(xa**(-4)))
fac3=(gam*temp)**(exp(-((xa-1)**2)/(2*(sigma**2))))

y=fac1*fac2*fac3
y=y/maxval(y)

return

end subroutine jonswapgk

! -----------------------------------------------------------
! ---- Small subroutine to determine f-range round peak -----
! ----(used by build_jonswap, swanreader, vardensreader)-----
subroutine frange(par,Sf,firstp,lastp,findlineout)

use params

implicit none

type(parameters)                        :: par


real*8, dimension(:), intent(in)        :: Sf
integer, intent(out)                    :: firstp, lastp

real*8, dimension(:), intent(out)       :: findlineout
real*8, dimension(:),allocatable        :: temp, findline
integer                                 :: i = 0 

allocate(findline(size(Sf)))
findline=0*Sf                           ! find frequency range around peak

where (Sf>par%sprdthr*maxval(Sf)) 
    findline=1
end where


firstp=maxval(maxloc(findline))         ! Picks the first "1" in temp

allocate (temp(size(findline)))
temp=(/(i,i=1,size(findline))/)
lastp=maxval(maxloc(temp*findline))     ! Picks the last "1" in temp

findlineout=findline
deallocate(temp, findline)

end subroutine frange


! -----------------------------------------------------------
! ----------- Small subroutine to determine tpD -------------
! ----(used by build_jonswap, swanreader, vardensreader)-----
subroutine tpDcalc(Sf,f,Trep)

implicit none

real*8, dimension(:), intent(in)        :: Sf, f
real*8, intent(out)                     :: Trep

real*8, dimension(:),allocatable        :: temp

allocate(temp(size(Sf)))
temp=0.d0
where (Sf>0.8*maxval(Sf)) 
    temp=1.d0
end where

Trep=sum(temp*Sf*f)/sum(temp*Sf)

end subroutine tpDcalc


! --------------------------------------------------------------
! --------------------- Dispersion relation --------------------
! ----------------- (used only by build_boundw) ----------------
subroutine bc_disper(k1,w1,m,h,g)
!          k  = wave number             (2 * pi / wave length)
!          w  = wave angular frequency  (2 * pi / wave period)
!          m  = size k and w vectors
!          h  = water depth
!          g  = gravitational acceleration constant, optional (DEFAULT 9.81)
!
!          absolute error in k*h < 5.0e-16 for all k*h
!
!
!          original Matlab code by: G. Klopman, Delft Hydraulics, 6 Dec 1994

integer, intent(in)                     :: m
real*8,dimension(m),intent(in)          :: w1
real*8,dimension(m),intent(out)         :: k1
real*8, intent(in)                      :: h, g

! internal variables

real*8,dimension(m)                     :: w2,q,thq,thq2,a,b,c,arg,sign
integer                                 :: j
real*8                                  :: hu

w2 = w1**2*(h/g)
q = w2/(1.0d0-exp(-(w2**(5.0d0/4.0d0))))*(2.0d0/5.0d0)

do j=1,4
    thq  = tanh(q)
    thq2 = 1.0d0-thq**2
    a    = (1.0d0-q*thq)*thq2
    b    = thq + q*thq2
    c    = q*thq-w2
    where (abs(a*c)<(b**2*1.0e-8))
        arg = -c/b
        elsewhere
                arg  = (b**2)-4.0d0*a*c
                arg  = (-b + sqrt(arg))/(2.0d0*a)
        endwhere
    q    = q+arg
end do

where (w1>0.0d0) sign=1.0d0

where (w1==0.0d0) sign=0.0d0

where (w1<0.0d0) sign=-1.0d0

k1 = sign*q/h

where (k1==huge(hu)) k1=0.0d0

where (k1==-1.0d0*huge(hu)) k1=0.0d0

return

end subroutine bc_disper



subroutine checkbcfilelength(par,filename)

use params
use xmpi_module

IMPLICIT NONE

type(parameters), INTENT(IN)             :: par
character*80      :: filename,dummy
character*8       :: testc
character*1       :: ch
integer           :: i,ier,nlines,filetype
real*8            :: t,dt,total,d1,d2,d3,d4,d5


if (xmaster) then
	open(741,file=filename)
	i=0
	do while (ier==0)
	   read(741,'(a)',iostat=ier)ch
	   if (ier==0)i=i+1
	enddo
	nlines=i
	rewind(741)    
		
	if (par%instat==4 .or. par%instat==5 .or. par%instat==6) then 
		read(741,*)testc
		if (testc=='FILELIST') then
			filetype = 1
			nlines=nlines-1
		else
			filetype = 0
		endif
	elseif (par%instat==40 .or. par%instat==41) then
        filetype = 2
	endif
    
	total=0.d0
	select case (filetype)
	    case(0)
			total=2.d0*par%tstop*max(par%morfac,1.d0)
		case(1)
			do i=1,nlines
				read(741,*)t,dt,dummy
				total=total+t
			enddo
		case(2)
			do i=1,nlines
				read(741,*)d1,d2,d3,d4,d5,t,dt
				total=total+t
			enddo
	end select
	total=total/max(par%morfac,1.d0)
	
	close(741)
	if (total<par%tstop) then
		write(*,*)'Error !!!! Wave boundary condition time series too short. Stopping calculation !!!!'
		call halt_program
	endif

endif

end subroutine checkbcfilelength

end module waveparams
