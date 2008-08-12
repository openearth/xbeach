module initialize 
contains
subroutine wave_init (s,par)
use params
use spaceparams
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

if (xmaster) then
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
  allocate(s%BR(1:s%nx+1,1:s%ny+1))
  allocate(s%kb(1:s%nx+1,1:s%ny+1))
  allocate(s%Tbore(1:s%nx+1,1:s%ny+1))
  allocate(s%uon(1:s%nx+1,1:s%ny+1))
  allocate(s%uoff(1:s%nx+1,1:s%ny+1))
  allocate(s%ua(1:s%nx+1,1:s%ny+1))  
  allocate(s%dzav(1:s%nx+1,1:s%ny+1))  

  s%umean=0.

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
  s%ust       = 0.d0
  s%tm        = 0.d0
  s%uwf       = 0.d0
  s%vwf       = 0.d0
  s%ustr      = 0.d0
  s%usd       = 0.d0
  s%bi        = 0.d0
  s%DR        = 0.d0
  s%BR        = par%Beta
  s%kb        = 0.d0
  s%Tbore     = 0.d0
  s%uon       = 0.d0
  s%uoff      = 0.d0
  s%ua        = 0.d0
  s%dzav      = 0.d0


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

endif

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

if(xmaster) then
  allocate(s%zs(1:s%nx+1,1:s%ny+1))
  allocate(s%dzsdt(1:s%nx+1,1:s%ny+1))
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


  ! cjaap: replaced par%hmin by par%eps
  s%hh=max(s%zs0-s%zb,par%eps)
  s%hu = s%hh !Jaap
  s%hv = s%hh !Jaap
  s%hum(1:s%nx,:) = 0.5d0*(s%hh(1:s%nx,:)+s%hh(2:s%nx+1,:))
  s%hum(s%nx+1,:)=s%hh(s%nx+1,:)
  s%zs=0.d0
  s%zs=max(s%zb,s%zs0)
  s%dzsdt=0.d0
  s%dzbdt=0.d0
  s%uu=0.d0
  s%u=0.d0
  s%vv=0.d0
  s%vu=0.d0
  s%qx=0.d0
  s%qy=0.d0
  s%sedero=0.d0
  s%vmagu=0.d0
  s%vmagv=0.d0
  s%vmageu=0.d0
  s%vmagev=0.d0
  s%maxzs=-999.d0
  s%minzs=999.d0

endif


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

integer                             :: i,j,m
character*80                        :: fnameg

allocate(s%ccg(1:s%nx+1,1:s%ny+1,par%ngd))
allocate(s%dcdy(1:s%nx+1,1:s%ny+1))
allocate(s%dcdx(1:s%nx+1,1:s%ny+1))
allocate(s%Tsg(1:s%nx+1,1:s%ny+1,par%ngd)) 
allocate(s%Sug(1:s%nx+1,1:s%ny+1,par%ngd)) 
allocate(s%Svg(1:s%nx+1,1:s%ny+1,par%ngd)) 
allocate(s%vmag(1:s%nx+1,1:s%ny+1)) 
allocate(s%ceqg(1:s%nx+1,1:s%ny+1,par%ngd))
allocate(s%graindistr(1:s%nx+1,1:s%ny+1,1:par%nd,1:par%ngd))
allocate(s%D50(1:par%ngd))
allocate(s%D90(1:par%ngd))
allocate(s%sedcal(1:par%ngd))
!
! Specify Grain Distribution
!
if(xmaster) then
  if (par%ngd>1) then
     call readkey('params.txt','gdist',fnameg)
     open(34,file=fnameg)
     do j = 2,s%ny
        do m = 1,par%nd
               read(34,*)(s%graindistr(i,j,m,1),i=1,s%nx+1)
            enddo
     enddo
     close(34)
     s%graindistr(:,1,:,:)      = s%graindistr(:,2,:,:)
     s%graindistr(:,s%ny+1,:,:) = s%graindistr(:,s%ny,:,:)
     s%graindistr(:,:,:,2) = abs(1.d0-s%graindistr(:,:,:,1))
     s%graindistr(:,:,:,3) = 0.0d0
  else
     s%graindistr = 1.0 
  endif
endif  ! xmaster

s%D50(1)     = par%D501
s%D90(1)     = par%D901
s%sedcal(1)  = par%sedcal1 

if (par%ngd>1) then
   s%D50(2)     = par%D502
   s%D90(2)     = par%D902
   s%sedcal(2)  = par%sedcal2 
   s%D50(3)     = par%D503
   s%D90(3)     = par%D903
   s%sedcal(3)  = par%sedcal3 
endif

s%ccg        = 0.d0
s%Sug        = 0.d0
s%Svg        = 0.d0
s%dcdx       = 0.d0
s%dcdy       = 0.d0

end subroutine sed_init

end module initialize
