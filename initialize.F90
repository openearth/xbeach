module initialize 
contains
subroutine wave_init (s,par)
use params
use spaceparams
use readkey_module
use xmpi_module
use interp
use wave_timestep_module

IMPLICIT NONE

type(spacepars)                     :: s
type(parameters)                    :: par

integer                             :: i
integer                             :: j,ig
integer                             :: itheta,indt

real*8,dimension(:),allocatable     :: yzs0,szs0

  allocate(s%thetamean(1:s%nx+1,1:s%ny+1))
  allocate(s%Fx(1:s%nx+1,1:s%ny+1))
  allocate(s%Fy(1:s%nx+1,1:s%ny+1))
  allocate(s%Sxx(1:s%nx+1,1:s%ny+1))
  allocate(s%Sxy(1:s%nx+1,1:s%ny+1))
  allocate(s%Syy(1:s%nx+1,1:s%ny+1))
  allocate(s%n(1:s%nx+1,1:s%ny+1))
  allocate(s%H(1:s%nx+1,1:s%ny+1))
  allocate(s%hh(1:s%nx+1,1:s%ny+1))
  allocate(s%cgx(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%cgy(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%cx(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%cy(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%ctheta(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%thet(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%costhet(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%sinthet(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%sigt(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%ee(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%rr(1:s%nx+1,1:s%ny+1,1:s%ntheta))
  allocate(s%sigm(1:s%nx+1,1:s%ny+1))
  allocate(s%c(1:s%nx+1,1:s%ny+1))
  allocate(s%cg(1:s%nx+1,1:s%ny+1))
  allocate(s%k(1:s%nx+1,1:s%ny+1))
  allocate(s%ui(1:s%nx+1,1:s%ny+1))
  allocate(s%E(1:s%nx+1,1:s%ny+1)) 
  allocate(s%R(1:s%nx+1,1:s%ny+1)) 
  allocate(s%urms(1:s%nx+1,1:s%ny+1)) 
  allocate(s%D(1:s%nx+1,1:s%ny+1)) 
  allocate(s%Qb(1:s%nx+1,1:s%ny+1)) 
  allocate(s%ust(1:s%nx+1,1:s%ny+1)) 
  allocate(s%tm(1:s%nx+1,1:s%ny+1)) 
  allocate(s%uwf(1:s%nx+1,1:s%ny+1)) 
  allocate(s%vwf(1:s%nx+1,1:s%ny+1)) 
  allocate(s%ustr(1:s%nx+1,1:s%ny+1)) 
  allocate(s%usd(1:s%nx+1,1:s%ny+1))
  allocate(s%zs0(1:s%nx+1,1:s%ny+1)) 
  allocate(szs0(1:2)) 
  allocate(yzs0(1:2)) 
  allocate(s%bi(1:s%ny+1))
  allocate(s%DR(1:s%nx+1,1:s%ny+1)) 
  allocate(s%umean(2,1:s%ny+1))
  allocate(s%vmean(2,1:s%ny+1))
  allocate(s%umwci       (1:s%nx+1,1:s%ny+1))
  allocate(s%vmwci       (1:s%nx+1,1:s%ny+1))
  allocate(s%zswci       (1:s%nx+1,1:s%ny+1))
  allocate(s%BR(1:s%nx+1,1:s%ny+1))

  s%umean=0.
  s%vmean=0.

  if (par%tideloc>0) then

    par%zs01=s%tideinpz(1,1)

    if(par%tideloc.eq.1) par%zs02=par%zs01

    if(par%tideloc.eq.2 .and. par%paulrevere.eq.0) then
          par%zs03=s%tideinpz(1,2)
          par%zs02=par%zs01
          par%zs04=par%zs03
    endif 

    if(par%tideloc.eq.2 .and. par%paulrevere.eq.1) then
          par%zs02=s%tideinpz(1,2)
          par%zs03=0.d0
          par%zs04=0.d0
    endif

    if(par%tideloc.eq.4) then
          par%zs02=s%tideinpz(1,2)
          par%zs03=s%tideinpz(1,3)
          par%zs04=s%tideinpz(1,4)
    endif

    if(par%tideloc.eq.1) s%zs0 = s%x*0.0d0 + par%zs01

    if(par%tideloc.eq.2 .and. par%paulrevere.eq.1) then
          yzs0(1)=s%y(1,1)
          yzs0(2)=s%y(1,s%ny+1)
          szs0(1)=par%zs01
          szs0(2)=par%zs02
          do i = 1,s%ny+1
                  call LINEAR_INTERP(yzs0, szs0, 2, s%y(1,i), s%zs0(1,i), indt)
          enddo
          do j = 1,s%ny+1 
                  do i = 1,s%nx+1
                          s%zs0(i,j) = s%zs0(1,j)
                  enddo
          enddo 
    endif

    if(par%tideloc.eq.2 .and. par%paulrevere.eq.0) then
          yzs0(1)=s%x(1,1)
          yzs0(2)=s%x(s%nx+1,1)
          szs0(1)=par%zs01
          szs0(2)=par%zs04
          s%zs0(1,:)=par%zs01
          s%zs0(s%nx+1,:)=par%zs03
          do j = 1,s%ny+1 
                  do i = 1,s%nx+1
                          if (s%zb(i,j).gt.s%zs0(1,j)+par%eps) goto 302
                          s%zs0(i,j) = s%zs0(1,j)
                  enddo
                  
  302             if (i.lt.s%nx+1) then
                          do ig=i+1,s%nx+1
                                  s%zs0(ig,j) = s%zs0(s%nx+1,j) 
                          enddo
                  else
                          do i = 1,s%nx+1
                                  call LINEAR_INTERP(yzs0, szs0, 2, s%x(i,j), s%zs0(i,j), indt)
                          enddo
                  endif
          enddo
    endif

    if(par%tideloc.eq.4) then
          yzs0(1)=s%y(1,1)
          yzs0(2)=s%y(1,s%ny+1)
          szs0(1)=par%zs01
          szs0(2)=par%zs02
          do i = 1,s%ny+1
                  call LINEAR_INTERP(yzs0, szs0, 2, s%y(1,i), s%zs0(1,i), indt)
          enddo
          yzs0(1)=s%y(s%nx+1,1)
          yzs0(2)=s%y(s%nx+1,s%ny+1)
          szs0(1)=par%zs04
          szs0(2)=par%zs03
          do i = 1,s%ny+1
                  call LINEAR_INTERP(yzs0, szs0, 2, s%y(s%nx+1,i), s%zs0(s%nx+1,i), indt)
          enddo

          do j = 1,s%ny+1 
                  do i = 1,s%nx+1
                          if (s%zb(i,j).gt.s%zs0(1,j)+par%eps) goto 303
                          s%zs0(i,j) = s%zs0(1,j)
                  enddo
                  
  303             if (i.lt.s%nx+1) then
                          do ig=i+1,s%nx+1
                                  s%zs0(ig,j) = s%zs0(s%nx+1,j) 
                          enddo
                  else
                          do i = 1,s%nx+1
                                  call LINEAR_INTERP(yzs0, szs0, 2, s%x(i,j), s%zs0(i,j), indt)
                          enddo
                  endif
          enddo
    endif

  else

     s%zs0 = par%zs01

  endif

  s%zs0 = max(s%zs0,s%zb)
  !Dano s%hh = max(s%zs0-s%zb,par%eps);
  s%hh = s%zs0-s%zb
  !
  ! Initial condition
  !
  do itheta=1,s%ntheta
      do j=1,s%ny+1
          do i=1,s%nx+1
              s%ee(i,j,itheta)=0.d0
          end do
      end do
  end do

  s%thetamean = 0.d0
  s%Fx        = 0.d0
  s%Fy        = 0.d0
  s%Sxx       = 0.d0
  s%Sxy       = 0.d0
  s%Syy       = 0.d0
  s%n         = 0.d0
  s%H         = 0.d0
  s%cgx       = 0.d0
  s%cgy       = 0.d0
  s%cx        = 0.d0
  s%cy        = 0.d0
  s%ctheta    = 0.d0
  s%thet      = 0.d0
  s%sigt      = 0.d0
  s%rr        = 0.d0
  s%sigm      = 0.d0
  s%c         = 0.d0
  s%cg        = 0.d0
  s%k         = 0.d0
  s%ui        = 0.d0
  s%E         = 0.d0
  s%R         = 0.d0
  s%urms      = 0.d0
  s%D         = 0.d0
  s%Qb        = 0.d0
  s%ust       = 0.d0
  s%tm        = 0.d0
  s%uwf       = 0.d0
  s%vwf       = 0.d0
  s%ustr      = 0.d0
  s%usd       = 0.d0
  s%bi        = 0.d0
  s%DR        = 0.d0
  s%BR        = par%Beta


  ! added for bound long wave bc Ad 27 march 2006
  do itheta=1,s%ntheta
      s%thet(:,:,itheta) = s%theta(itheta)
      s%costhet(:,:,itheta)= cos(s%theta(itheta))
      s%sinthet(:,:,itheta)= sin(s%theta(itheta))
  end do

  ! introduce intrinsic frequencies for wave action
  if (par%instat==4 .or. par%instat==5 .or. par%instat==6 .or. par%instat==7 .or. par%instat==8) par%Trep=10.d0 !Robert
  ! incorrect values are computed below for instat = 4/5/6/7
  ! in this case right values are computed in wave params.f90
  do itheta=1,s%ntheta
      s%sigt(:,:,itheta) = 2*par%px/par%Trep
  end do
  s%sigm = sum(s%sigt,3)/s%ntheta
  call dispersion(par,s)


end subroutine wave_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine flow_init (s,par)
use params
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars)                     :: s
type(parameters)                    :: par

integer*4                           :: i,j
logical                             :: exists
integer*4                           :: iUnit

  allocate(s%zs(1:s%nx+1,1:s%ny+1))
  allocate(s%dzsdt(1:s%nx+1,1:s%ny+1))
  allocate(s%dzsdx(1:s%nx+1,1:s%ny+1))
  allocate(s%dzsdy(1:s%nx+1,1:s%ny+1))
  allocate(s%dzbdt(1:s%nx+1,1:s%ny+1))
  allocate(s%uu(1:s%nx+1,1:s%ny+1))
  allocate(s%vv(1:s%nx+1,1:s%ny+1))
  allocate(s%qx(1:s%nx+1,1:s%ny+1))
  allocate(s%qy(1:s%nx+1,1:s%ny+1))
  allocate(s%sedero(1:s%nx+1,1:s%ny+1))
  allocate(s%ueu(1:s%nx+1,1:s%ny+1)) 
  allocate(s%vev(1:s%nx+1,1:s%ny+1)) 
  allocate(s%vmagu(1:s%nx+1,1:s%ny+1)) 
  allocate(s%vmagv(1:s%nx+1,1:s%ny+1)) 
  allocate(s%vmageu(1:s%nx+1,1:s%ny+1)) 
  allocate(s%vmagev(1:s%nx+1,1:s%ny+1)) 
  allocate(s%u(1:s%nx+1,1:s%ny+1)) 
  allocate(s%v(1:s%nx+1,1:s%ny+1)) 
  allocate(s%ue(1:s%nx+1,1:s%ny+1)) 
  allocate(s%ve(1:s%nx+1,1:s%ny+1)) 
  allocate(s%hold(1:s%nx+1,1:s%ny+1)) 
  allocate(s%wetu(1:s%nx+1,1:s%ny+1)) 
  allocate(s%wetv(1:s%nx+1,1:s%ny+1)) 
  allocate(s%wetz(1:s%nx+1,1:s%ny+1))
  allocate(s%hu(1:s%nx+1,1:s%ny+1)) 
  allocate(s%hv(1:s%nx+1,1:s%ny+1))
  allocate(s%hum(1:s%nx+1,1:s%ny+1)) 
  allocate(s%hvm(1:s%nx+1,1:s%ny+1))
  allocate(s%vu(1:s%nx+1,1:s%ny+1))
  allocate(s%uv(1:s%nx+1,1:s%ny+1))
  allocate(s%maxzs(1:s%nx+1,1:s%ny+1))
  allocate(s%minzs(1:s%nx+1,1:s%ny+1))
  
!  if (par%nonh==1) then   
! Robert: required to stop MPI crash
    allocate(s%ws(1:s%nx+1,1:s%ny+1))
    allocate(s%wb(1:s%nx+1,1:s%ny+1))
    allocate(s%pres(1:s%nx+1,1:s%ny+1))
    s%ws   = 0.0d0
    s%wb   = 0.0d0
    s%pres = 0.0d0
!  endif

  ! cjaap: replaced par%hmin by par%eps
  s%hh=max(s%zs0-s%zb,par%eps)
  s%zs=0.d0
  s%zs=max(s%zb,s%zs0)
  
  !For certain tests I need to prescribe the initial water elevation (Pieter)
  !Disabled in the case of mpi

#ifndef USEMPI  
  inquire(file=trim('zsfile.ini'),EXIST=exists)
  if (exists) then
    !
    iUnit = 10000
    !Check if unit is free
    inquire(unit=iUnit,OPENED=exists)
    if (.not. exists) then
      open(iUnit,file='zsfile.ini')
      do j=1,s%ny+1
        !
        read(iUnit,*)(s%zs(i,j),i=1,s%nx+1)
        !
      enddo
      close(iUnit)
      s%hh = max(s%zs-s%zb,par%eps)
    endif  
    !
  endif
#endif

  
  !Initialize hu correctly to prevent spurious initial flow (Pieter)
  do j=1,s%ny+1
    do i=1,s%nx
      s%hu(i,j) = max(s%zs(i,j),s%zs(i+1,j))-max(s%zb(i,j),s%zb(i+1,j))
    enddo
  enddo
  s%hu(s%nx+1,:)=s%hu(s%nx,:)

  do j=1,s%ny
    do i=1,s%nx+1
      s%hv(i,j) = max(s%zs(i,j),s%zs(i,j+1))-max(s%zb(i,j),s%zb(i,j+1))
    enddo
  enddo
  s%hv(:,s%ny+1)=s%hv(:,s%ny+1)

  s%hum(1:s%nx,:) = 0.5d0*(s%hh(1:s%nx,:)+s%hh(2:s%nx+1,:))
  s%hum(s%nx+1,:)=s%hh(s%nx+1,:)

  s%dzsdt=0.d0
  s%dzsdx=0.d0
  s%dzsdy=0.d0
  s%dzbdt=0.d0
  s%uu=0.d0
  s%u=0.d0
  s%vv=0.d0
  s%v=0.d0
  s%vu=0.d0
  s%uv=0.d0
  s%qx=0.d0
  s%qy=0.d0
  s%sedero=0.d0
  s%vmagu=0.d0
  s%vmagv=0.d0
  s%vmageu=0.d0
  s%vmagev=0.d0
  s%maxzs=-999.d0
  s%minzs=999.d0



end subroutine flow_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine sed_init (s,par)
use params
use spaceparams
use readkey_module
use xmpi_module

IMPLICIT NONE

type(spacepars)                     :: s
type(parameters)                    :: par

integer                             :: i,j,m,jg,start,rgd
character*80                        :: fnameg,line
character*4                         :: tempc
real*8                              :: tempr



   rgd = par%ngd



   allocate(s%ccg(1:s%nx+1,1:s%ny+1,par%ngd))
	allocate(s%dcbdy(1:s%nx+1,1:s%ny+1))
	allocate(s%dcbdx(1:s%nx+1,1:s%ny+1))
	allocate(s%dcsdy(1:s%nx+1,1:s%ny+1))
	allocate(s%dcsdx(1:s%nx+1,1:s%ny+1))
	allocate(s%Tsg(1:s%nx+1,1:s%ny+1,par%ngd)) 
	allocate(s%Susg(1:s%nx+1,1:s%ny+1,par%ngd)) 
	allocate(s%Svsg(1:s%nx+1,1:s%ny+1,par%ngd)) 
	allocate(s%Subg(1:s%nx+1,1:s%ny+1,par%ngd)) 
	allocate(s%Svbg(1:s%nx+1,1:s%ny+1,par%ngd)) 
	allocate(s%vmag(1:s%nx+1,1:s%ny+1)) 
	allocate(s%ceqsg(1:s%nx+1,1:s%ny+1,par%ngd))
	allocate(s%ceqbg(1:s%nx+1,1:s%ny+1,par%ngd))
	allocate(s%D50(1:par%ngd))
	allocate(s%D90(1:par%ngd))
	allocate(s%sedcal(1:par%ngd))
	allocate(s%ucrcal(1:par%ngd))
	allocate(s%nd(1:s%nx+1,1:s%ny+1))
	allocate(s%dzbed(1:s%nx+1,1:s%ny+1,1:max(par%nd,2))) 
	allocate(s%pbbed(1:s%nx+1,1:s%ny+1,1:max(par%nd,2),1:par%ngd)) 
	allocate(s%z0bed(1:s%nx+1,1:s%ny+1))
	allocate(s%ureps(1:s%nx+1,1:s%ny+1))
	allocate(s%urepb(1:s%nx+1,1:s%ny+1))
	allocate(s%vreps(1:s%nx+1,1:s%ny+1))
	allocate(s%vrepb(1:s%nx+1,1:s%ny+1))
	allocate(s%dzbdx(1:s%nx+1,1:s%ny+1))
	allocate(s%dzbdy(1:s%nx+1,1:s%ny+1))
	allocate(s%ero(1:s%nx+1,1:s%ny+1,1:par%ngd))
	allocate(s%depo_ex(1:s%nx+1,1:s%ny+1,1:par%ngd))
	allocate(s%depo_im(1:s%nx+1,1:s%ny+1,1:par%ngd))
    allocate(s%kb(1:s%nx+1,1:s%ny+1))
    allocate(s%Tbore(1:s%nx+1,1:s%ny+1))
    allocate(s%ua(1:s%nx+1,1:s%ny+1))  
    allocate(s%dzav(1:s%nx+1,1:s%ny+1))  
    allocate(s%Sk(1:s%nx+1,1:s%ny+1))
    allocate(s%As(1:s%nx+1,1:s%ny+1))
    allocate(s%kturb(1:s%nx+1,1:s%ny+1))
    allocate(s%rolthick(1:s%nx+1,1:s%ny+1))


  

 !
 ! Set grain size(s)
 !
 call readkey('params.txt','D50',line)
 if (line=='') then
    s%D50=0.0002d0   ! Default
 else
    read(line,*) s%D50(1:rgd)
 endif
 if (par%form<4 .and. maxval(s%D50)>0.002) write(*,*) 'D50 > 2mm, out of validity range'
 call readkey('params.txt','D90',line)
 if (line=='') then
    s%D90=0.0003d0   ! Default
 else
    read(line,*) s%D90(1:rgd)
 endif
 call readkey('params.txt','sedcal',line)
 if (line=='') then
    s%sedcal=1.d0    ! Default
 else
    read(line,*) s%sedcal(1:rgd)
 endif
 call readkey('params.txt','ucrcal',line)
 if (line=='') then
    s%ucrcal=1.d0    ! Default
 else
    read(line,*) s%ucrcal(1:rgd)
 endif

if (rgd==1) then
! No multi sediment, but we do need some data to keep the script running
   s%pbbed=1.d0
	if (par%nd==1) then
	   par%nd_var=1
		s%dzbed = 99999.d0
	else
	   par%nd_var=1
		s%dzbed = par%dzg1
	endif
else
   !
   ! Fill s%pbed en s%dzbed
   ! 
   do jg=1,rgd
	  write(tempc,'(i4)')jg
      start=4-floor(log10(real(jg)))
	  write(fnameg,'(a,a,a)')'gdist',tempc(start:4),'.inp'
      open(31,file=fnameg)
	  do m=1,par%nd
        do j=1,s%ny+1
	        read(31,*)(s%pbbed(i,j,m,jg),i=1,s%nx+1)
		  enddo
	  enddo
	  close(31)
   enddo
   ! Rework pbbed so that sum fractions = 1
   do m=1,par%nd
	   do j=2,s%ny
		   do i=2,s%nx
		       
			   tempr=sum(s%pbbed(i,j,m,:))
			   if (abs(1.d0-tempr)>0.d0) then
				   write(*,*)' Warning: Resetting sum of sediment fractions in point (',&
					                     i,',',j,') layer ,',m,&
												' to equal unity.'
				   s%pbbed(i,j,m,:)=s%pbbed(i,j,m,:)/tempr
				endif
			enddo
		enddo
	enddo
	! boundary neumann
	s%pbbed(1,:,:,:)=s%pbbed(2,:,:,:)
	s%pbbed(s%nx+1,:,:,:)=s%pbbed(s%nx,:,:,:)
	s%pbbed(:,1,:,:)=s%pbbed(:,2,:,:)
	s%pbbed(:,s%ny+1,:,:)=s%pbbed(:,s%ny,:,:)
	! sediment thickness				   
	s%dzbed(:,:,1:par%nd_var-1)       = par%dzg1
	s%dzbed(:,:,par%nd_var)           = par%dzg2
	s%dzbed(:,:,par%nd_var+1:par%nd)  = par%dzg3
endif

! bottom of sediment model
s%z0bed = s%zb - sum(s%dzbed,DIM=3)

s%nd = max(par%nd,2)

s%ureps      = 0.d0
s%vreps      = 0.d0
s%ccg        = 0.d0
s%ceqbg      = 0.d0
s%ceqsg      = 0.d0
s%Susg       = 0.d0
s%Svsg       = 0.d0
s%Subg       = 0.d0
s%Svbg       = 0.d0
s%dcsdx      = 0.d0
s%dcsdy      = 0.d0
s%dcbdx      = 0.d0
s%dcbdy      = 0.d0
s%ero        = 0.d0
s%depo_im    = 0.d0
s%depo_ex    = 0.d0
s%kb         = 0.d0
s%Tbore      = 0.d0
s%ua         = 0.d0
s%dzav       = 0.d0
s%Sk         = 0.d0
s%As         = 0.d0
s%kturb      = 0.d0
s%rolthick   = 0.d0


end subroutine sed_init

end module initialize
