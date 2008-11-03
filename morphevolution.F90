module morphevolution
contains
subroutine transus(s,par)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Copyright (C) 2007 UNESCO-IHE, WL|Delft Hydraulics and Delft University !
! Dano Roelvink, Ap van Dongeren, Ad Reniers, Jamie Lescinski,            !
! Jaap van Thiel de Vries, Robert McCall                                  !       
!                                                                         !
! d.roelvink@unesco-ihe.org                                               !
! UNESCO-IHE Institute for Water Education                                !
! P.O. Box 3015                                                           !
! 2601 DA Delft                                                           !
! The Netherlands                                                         !
!                                                                         !
! This library is free software; you can redistribute it and/or           !
! modify it under the terms of the GNU Lesser General Public              !
! License as published by the Free Software Foundation; either            !
! version 2.1 of the License, or (at your option) any later version.      !
!                                                                         !
! This library is distributed in the hope that it will be useful,         !
! but WITHOUT ANY WARRANTY; without even the implied warranty of          !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        !
! Lesser General Public License for more details.                         !
!                                                                         !
! You should have received a copy of the GNU Lesser General Public        !
! License along with this library; if not, write to the Free Software     !
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307     !
! USA                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use params
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars),target                   :: s
type(parameters)                         :: par

integer                                  :: i
integer                                  :: j,jg

real*8,dimension(:,:),allocatable,save   :: vmag2,uau,uav
real*8 , dimension(s%nx+1,s%ny+1)        :: dzbx,dzby,dzremain,source
real*8,dimension(:,:),allocatable,save   :: cc,cu,cv,Su,Sv,Dc,termh

include 's.ind'
include 's.inp'

if (.not. allocated(vmag2)) then
   allocate(vmag2 (nx+1,ny+1))
   allocate(uau (nx+1,ny+1))
   allocate(uav (nx+1,ny+1))
   allocate(cu (nx+1,ny+1))
   allocate(cv (nx+1,ny+1))
   allocate(cc (nx+1,ny+1))
   allocate(Su (nx+1,ny+1))
   allocate(Sv (nx+1,ny+1))
   allocate(Dc (nx+1,ny+1))
   allocate(termh (nx+1,ny+1))
endif
! use eulerian velocities
vmag2    = ue**2+ve**2
dcdx     = 0.0d0
dcdy     = 0.0d0
! calculate equilibrium concentration
if (par%form==1) then           ! soulsby van Rijn
        call sb_vr(s,par)
elseif (par%form==2) then       ! Van Thiel de Vries & Reniers 2008
    call sednew(s,par)
end if

! compute diffusion coefficient

if (par%nuhfac==1) then
   termh = hh/max(H,.01d0)
   ! termh = max(hh(i,j)/2.d0,0.01d0)
   termh = min(termh,10.d0);
   ! Dc = 5*(DR/par%rho)**(1.d0/3.d0)*hh/(exp(termh)-1.d0)
   Dc = par%nuh+par%nuhfac*hh*(DR/par%rho)**(1.d0/3.d0)
   ! Dc = par%dico
else
   Dc = par%dico
end if

do jg = 1,par%ngd
   cc = ccg(:,:,jg)
!   cc = ceqg(:,:,jg) ! Can be used to test total transport mode
   !
   ! X-direction
   do j=1,ny+1
      do i=1,nx
         if(ueu(i,j)>0.d0) then
   !test          cu(i,j)=cc(i,j)
           cu(i,j)=par%thetanum*cc(i,j)+(1.d0-par%thetanum)*cc(min(i+1,nx),j)
         elseif(ueu(i,j)<0.d0) then
            cu(i,j)=par%thetanum*cc(i+1,j)+(1.d0-par%thetanum)*cc(max(i,2),j)
         else
            cu(i,j)=0.5d0*(cc(i,j)+cc(i+1,j))
         endif
         dcdx(i,j)=(cc(i+1,j)-cc(i,j))/(xz(i+1)-xz(i)) !Jaap
      enddo 
   enddo
   ! wwvv dcdx(nx:1,:) is still untouched, correct this ofr the parallel case
#ifdef USEMPI
   call xmpi_shift(dcdx,'m:')
#endif
   cu(nx+1,:) = cc(nx+1,:) !Robert
   ! wwvv fix this in parallel case
#ifdef USEMPI
   call xmpi_shift(cu,'m:')
#endif
   !
   ! Bed slope terms
   !
   dzbx=0.d0
   do j=1,ny+1
       do i=1,nx
           dzbx(i,j)=(zb(i+1,j)-zb(i,j))/(xz(i+1)-xz(i))
       enddo
   enddo
   ! wwvv in parallel version, there will be a discrepancy between the values
   ! of dzbx(nx+1,:).
   !wwvv so fix that
#ifdef USEMPI
   call xmpi_shift(dzbx,'m:')
#endif
   ! Jaap: get ua in u points and split out in u and v direction
   uau(1:nx,:) = 0.5*(ua(1:nx,:)*cos(theta0)+ua(2:nx+1,:)*cos(theta0))
   uav(1:nx,:) = 0.5*(ua(1:nx,:)*sin(theta0)+ua(2:nx+1,:)*sin(theta0))
   ! Jaap: compute vmagu including ua
   vmagu = sqrt((uu+uau)**2+(vu+uav)**2)
   !
   Su=(cu*(ueu+uau)*hu-Dc*hu*dcdx-par%facsl*cu*vmagu*hu*dzbx)*wetu   !

   !
   ! Y-direction
   !
   do j=1,ny
      do i=1,nx+1
         if(vev(i,j)>0) then
        !   cv(i,j)=cc(i,j)
            cv(i,j)=par%thetanum*cc(i,j)+(1.d0-par%thetanum)*cc(i,min(j+1,ny))
          else if(vev(i,j)<0) then
        !   cv(i,j)=cc(i,j+1)
            cv(i,j)=par%thetanum*cc(i,j+1)+(1.d0-par%thetanum)*cc(i,max(j,2))
        else
         cv(i,j)=0.5d0*(cv(i,j)+cv(i,j+1))
         end if
         dcdy(i,j)=(cc(i,j+1)-cc(i,j))/(yz(j+1)-yz(j)) !Jaap
      end do 
    end do
    ! wwvv dcdy(:,ny+1) is not filled in, so in parallel case:
#ifdef USEMPI
    call xmpi_shift(dcdy,':n')
#endif
    cv(:,ny+1) = cc(:,ny+1) !Robert
    !
    ! Bed slope terms
    !
    dzby=0.d0
    do j=1,ny
        do i=1,nx+1
            dzby(i,j)=(zb(i,j+1)-zb(i,j))/(yz(j+1)-yz(j))
        enddo
    enddo
   ! wwvv in parallel version, there will be a discrepancy between the values
   ! of dzby(:,ny+1).
    ! wwvv so fix that
#ifdef USEMPI
    call xmpi_shift(dzby,':n')
#endif
	! Jaap: get ua in v points and split out in u and v direction
	uau(:,1:ny) = 0.5*(ua(:,1:ny)*sin(theta0)+ua(:,2:ny+1)*sin(theta0))
    uav(:,1:ny) = 0.5*(ua(:,1:ny)*cos(theta0)+ua(:,2:ny+1)*cos(theta0))
    ! Jaap: compute vmagu including ua
    vmagv = sqrt((uv+uau)**2+(vv+uav)**2)
    !
    Sv=(cv*(vev+uav)*hv-Dc*hv*dcdy-par%facsl*cv*vmagv*hv*dzby)*wetv
    ! Jaap: compute remaining sediment thickness above hard layer
    dzremain = max(0.d0,dzlayer+sedero)
    source = (1-par%por)*dzremain/par%dt/max(par%morfac,1.d0)
    ! Jaap: we now already how much sand is picked up from the bed
    dzbdt=0.0d0
    ! dzbdt = -1.d0/(1-par%por)*min(source,hold*(ceqg(:,:,jg)*graindistr(:,:,1,jg)-cc)/Tsg(:,:,jg))
    ! Jaap: First compute cc before updating sediment transports...
    do j=2,ny+1
      do i=2,nx+1
         cc(i,j) = hold(i,j)*cc(i,j)-par%dt*((Su(i,j)-Su(i-1,j))/(xu(i)-xu(i-1))+&
                                             (Sv(i,j)-Sv(i,j-1))/(yv(j)-yv(j-1))-&
	  									     !hold(i,j)*(ceqg(i,j,jg)*graindistr(i,j,1,jg)-cc(i,j))/Tsg(i,j,jg))
											 !Jaap: set source to zero in case of hard layer near surface...
                                             min(source(i,j),hold(i,j)*(ceqg(i,j,jg)*graindistr(i,j,1,jg)-cc(i,j))/Tsg(i,j,jg)))

         cc(i,j)=max(cc(i,j),0.0d0)
      enddo
    enddo
    do j=1,ny+1
       do i=1,nx+1
          if(hh(i,j)>=par%hmin) then 
             cc(i,j)=cc(i,j)/hh(i,j)  
          else
             cc(i,j)=0.d0
          end if
       end do
    end do
    cc(1,:)=cc(2,:)
    cc(:,1)=cc(:,2)
    cc(nx+1,:)=cc(nx+1-1,:)
    cc(:,ny+1)=cc(:,ny+1-1)
    ! wwvv fix the first and last rows and columns of cc in parallel case
#ifdef USEMPI
    call xmpi_shift(cc,'1:')
    call xmpi_shift(cc,'m:')
    call xmpi_shift(cc,':1')
    call xmpi_shift(cc,':n')
#endif
    !Jaap
    cc=cc*wetz
    !
    ! wwvv border columns and rows of ccg Svg and Sug have to be communicated
    ccg(:,:,jg) = cc
    Svg(:,:,jg) = Sv
    Sug(:,:,jg) = Su
end do

vmag=sqrt(max(vmag2,par%umin))

end subroutine transus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bed_update(s,par)
use params
use spaceparams

IMPLICIT NONE

type(spacepars),target              :: s
type(parameters)                    :: par

integer                             :: i,ii
integer                             :: j,jg,jd
real*8                              :: dzb,dzmax,dz
real*8 , dimension(s%nx+1,s%ny+1)   :: dzbx,dzby,Su,Sv,dzremain
real*8 , dimension(s%nx+1,s%ny+1)   :: dzbtot,fact0,fact1,fact2,fact3,fact4,fact5
real*8 , dimension(s%nx+1,s%ny+1,par%ngd) :: dzbdtg
real*8 , dimension(s%nx+1,s%ny+1,par%nd)  :: graindistrm
logical                                   :: aval

include 's.ind'
include 's.inp'

dzbdtg = 0.d0
dzbtot = 0.d0

! Update bed level using continuity eq.
dzb=0.d0 !!!Ap
dzbdt=0.0d0

if (par%t>=par%morstart .and. par%morfac > .999d0) then

do jg = 1,par%ngd
   ! Update bed level using continuity eq.
   Su = Sug(:,:,jg)
   Sv = Svg(:,:,jg)

   ! Ap and Jaap: Bottom update is morfac dependent....
   do j=2,ny !Jaap nx instead of ny+1
      do i=2,nx  !Jaap nx instead of nx+1
         dzbdt(i,j)=1/(1-par%por)*        &
             (-(Su(i,j)-Su(i-1,j))/(xu(i)-xu(i-1)) &
              -(Sv(i,j)-Sv(i,j-1))/(yv(j)-yv(j-1)))
         zb(i,j)=zb(i,j)+dzbdt(i,j)*par%morfac*par%dt
         sedero(i,j)=sedero(i,j)+dzbdt(i,j)*par%morfac*par%dt           ! total bed level changes for all sediment clases over all time steps so far
         dzbdtg(i,j,jg) = dzbdtg(i,j,jg) + dzbdt(i,j)*par%morfac*par%dt ! total bed level change for each sediment class at this time step
         dzbtot(i,j) = dzbtot(i,j) +  dzbdt(i,j)*par%morfac*par%dt      ! total bed level change for all sediment classes at this time step
       enddo
   enddo
enddo

fact0 = 0.d0
fact1 = 0.d0
fact2 = 0.d0
fact3 = 0.d0
fact4 = 0.d0
fact5 = 0.d0

! update sediment fractions in sediment layers
! cjaap
if (par%ngd>1) then

do j=2,ny 
   do i=2,nx  
      if (dzbtot(i,j)<0.d0) then
         fact0(i,j) = 0.d0      
         fact1(i,j) = 1.d0
         fact2(i,j) = -dzbtot(i,j)/par%dzg
         fact3(i,j) = 0.d0
         fact4(i,j) = (par%dzg+dzbtot(i,j))/par%dzg
         fact5(i,j) = -dzbtot(i,j)/par%dzg
      else
         fact0(i,j) = 1.d0  ! only in case of vegetated layer 
         fact1(i,j) = (par%dzg-dzbtot(i,j))/par%dzg
         fact2(i,j) = 0.d0
         fact3(i,j) = dzbtot(i,j)/par%dzg
         fact4(i,j) = (par%dzg-dzbtot(i,j))/par%dzg
         fact5(i,j) = 0.d0
      endif
   enddo
enddo

do jd = 1,par%nd-1
   do jg = 1,par%ngd
      if (jd == 1) then
         if (jg == 1) then
            graindistr(:,:,jd,jg) = dzbdtg(:,:,jg)/par%dzg + fact1*graindistr(:,:,jd,jg) + &
            fact2*graindistr(:,:,jd+1,jg) + max(dzbdtg(:,:,jg+1)/par%dzg,0.0d0)
         elseif (jg==2) then
            !graindistr(:,:,jd,jg) = dzbdtg(:,:,jg)/par%dzg + fact1*graindistr(:,:,jd,jg) + fact2*graindistr(:,:,jd+1,jg)
            graindistr(:,:,jd,jg) = min(dzbdtg(:,:,jg)/par%dzg,0.0d0) + fact1*graindistr(:,:,jd,jg) + &
            fact2*graindistr(:,:,jd+1,jg)
         else
            graindistr(:,:,jd,jg) = dzbdtg(:,:,jg)/par%dzg + fact1*graindistr(:,:,jd,jg) + fact2*graindistr(:,:,jd+1,jg) 
         endif
      else
         graindistr(:,:,jd,jg) = fact3*graindistr(:,:,jd-1,jg) + fact4*graindistr(:,:,jd,jg) + fact5*graindistr(:,:,jd+1,jg)
      endif
   enddo
enddo

! ensure total fractions equal 1
graindistrm = 0.d0
do jd = 1,par%nd
    do jg = 1,par%ngd
       graindistrm(:,:,jd) = graindistrm(:,:,jd) + graindistr(:,:,jd,jg)
    enddo
    do jg = 1,par%ngd
       graindistr(:,:,jd,jg) = graindistr(:,:,jd,jg)/max(graindistrm(:,:,jd),.0010d0)
    enddo
enddo

endif
!
! bed boundary conditions
! 
! Fix bed at back boundary, but allow sediment to pass though
!zb(nx+1,:) = zb(nx+1,:)+dzbdt(nx,:)*par%dt  !Ap
!zb(nx+1,:) = zb(nx+1,:)+dzbdt(nx,:)*par%dt  !Ap
!sedero(nx+1,:) = sedero(nx,:) !Ap
zb(:,1) = zb(:,2)
sedero(:,1) = sedero(:,2)
zb(:,ny+1) = zb(:,ny)
sedero(:,ny+1) = sedero(:,ny)
!
! Avalanching
!
do ii=1,nint(par%morfac)
   ! Jaap; update dzremain before avalanching
   dzremain = max(0.d0,dzlayer+sedero)
   
   aval=.false.
   dzbx=0.d0
   dzby=0.d0
   do j=1,ny+1
      do i=1,nx
         dzbx(i,j)=(zb(i+1,j)-zb(i,j))/(xz(i+1)-xz(i))
      enddo
   enddo
   !
   Fimpact = 0.0d0
   do i=2,nx-1
      do j=1,ny+1
         if(max(hh(i,j),hh(i+1,j))>par%hswitch+par%eps) then
            dzmax=par%wetslp;
         else
            dzmax=par%dryslp;
         end if
         if(abs(dzbx(i,j))>dzmax ) then
		    aval=.true.     
            dzb=sign(1.0d0,dzbx(i,j))*(abs(dzbx(i,j))-dzmax)*(xz(i+1)-xz(i));
			! Jaap: compute dzmax using Overton and Fisher
			if (par%impact==1) then
			   Fimpact(i,j) = min(hu(i-1,j),max(0.d0,zb(i+1,j)-zb(i,j)))*par%rho*(uu(i-1,j)**2 + 0.125d0*(urms(i,j)+urms(i-1,j))**2)
			   dz = min(par%dzmax,par%CE/((xz(i+1)-xz(i))**2*par%rhos*(1-par%por)*par%g)*Fimpact(i,j))
			else
			   dz = par%dzmax
		    endif	
            ! Jaap: limit dzb with dz and with dzremain
			if (dzb>=0) then
			   dzb=min(dzb,dzremain(i+1,j),dz*par%dt*(xz(i+1)-xz(i)))
			else
               dzb=max(dzb,-1.d0*dzremain(i,j)*(xu(i)-xu(i-1))/(xu(i+1)-xu(i)),-1.d0*dz*par%dt*(xz(i+1)-xz(i)))
			endif
            
			zb(i,j)=zb(i,j)+dzb*(xu(i+1)-xu(i))/(xu(i)-xu(i-1)) !Jaap make sure there is continuity of sediment in non uniform grids;
            sedero(i,j) = sedero(i,j)+dzb*(xu(i+1)-xu(i))/(xu(i)-xu(i-1))
            dzbtot(i,j) = dzbtot(i,j)+dzb*(xu(i+1)-xu(i))/(xu(i)-xu(i-1))  
            zs(i,j)=zs(i,j)+dzb*(xu(i+1)-xu(i))/(xu(i)-xu(i-1))
			dzav(i,j)= dzav(i,j)+dzb*(xu(i+1)-xu(i))/(xu(i)-xu(i-1))
			
			zb(i+1,j)=zb(i+1,j)-dzb
            sedero(i+1,j)=sedero(i+1,j)-dzb
            dzbtot(i+1,j) = dzbtot(i+1,j)-dzb
            zs(i+1,j)=zs(i+1,j)-dzb
			dzav(i+1,j)= dzav(i+1,j)-dzb  !Jaap compute total bed level change due to avalanching
		 end if
      end do
   end do

   !JJ: update y slopes after avalanching in X-direction seems more appropriate
   do j=1,ny
      do i=1,nx+1
         dzby(i,j)=(zb(i,j+1)-zb(i,j))/(yz(j+1)-yz(j))
      enddo
   enddo

   do j=2,ny-1;
        do i=1,nx+1
            if(max(hh(i,j),hh(i,j+1))>par%hswitch+par%eps) then
                dzmax=par%wetslp
            else
                dzmax=par%dryslp
            end if
            if(abs(dzby(i,j))>dzmax ) then 
               aval=.true. 
               dzb=sign(1.0d0,dzby(i,j))*(abs(dzby(i,j))-dzmax)*(yz(j+1)-yz(j));
               dzb=sign(1.0d0,dzb)*min(abs(dzb),par%dzmax*par%dt*(yz(j+1)-yz(j))); !0.005d0
               ! Jaap: limit dzb with dz and with dzremain
			   if (dzb>=0) then
			      dzb=min(dzb,dzremain(i,j+1),dz*par%dt*(yz(j+1)-yz(j)))
			   else
                  dzb=max(dzb,-1.d0*dzremain(i,j)*(yv(j)-yv(j-1))/(yv(j+1)-yv(j)),-1.d0*dz*par%dt*(yz(j+1)-yz(j)))
			   endif
		    
			   zb(i,j)=zb(i,j)+dzb*(yv(j+1)-yv(j))/(yv(j)-yv(j-1)) !Jaap make sure there is continuity of sediment in non uniform grids;
               sedero(i,j) = sedero(i,j)+dzb*(yv(j+1)-yv(j))/(yv(j)-yv(j-1))
               dzbtot(i,j) = dzbtot(i,j)+dzb*(yv(j+1)-yv(j))/(yv(j)-yv(j-1)) 
               zs(i,j)=zs(i,j)+dzb*(yv(j+1)-yv(j))/(yv(j)-yv(j-1))!
			   dzav(i,j)= dzav(i,j)+dzb*(yv(j+1)-yv(j))/(yv(j)-yv(j-71))
			
			   zb(i,j+1)=zb(i,j+1)-dzb
               sedero(i,j+1)=sedero(i,j+1)-dzb
               dzbtot(i,j+1) = dzbtot(i,j+1)-dzb
               zs(i,j+1)=zs(i,j+1)-dzb!
		       dzav(i,j+1)= dzav(i,j+1)-dzb  !Jaap compute total bed level change due to avalanching
            end if
        end do
   end do
   if (.not.aval) exit
end do
endif

end subroutine bed_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sb_vr(s,par)
use params
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars),target              :: s
type(parameters)                    :: par

integer                             :: i
integer                             :: j,jg
real*8                              :: dster,twothird
real*8                              :: m1,m2,m3,m4,m5,m6
real*8                              :: z0,Ass,delta
real*8                              :: Te,kvis,Sster,c1,c2,wster


real*8 , dimension(:,:),allocatable,save   :: vmag2,Cd,Asb,dhdx,dhdy,Ts,Ur,Bm,B1
real*8 , dimension(:,:),allocatable,save   :: urms2,Ucr,term1,term2
real*8 , dimension(:,:),allocatable,save   :: uandv,b,fslope,hloc,ceq

include 's.ind'
include 's.inp'

if (.not. allocated(vmag2)) then
   allocate (vmag2 (nx+1,ny+1))
   allocate (Cd    (nx+1,ny+1))
   allocate (Asb   (nx+1,ny+1))
   allocate (dhdx  (nx+1,ny+1))   ! not used wwvv
   allocate (dhdy  (nx+1,ny+1))   ! not used wwvv
   allocate (urms2 (nx+1,ny+1))
   allocate (Ucr   (nx+1,ny+1))
   allocate (term1 (nx+1,ny+1))
   allocate (term2 (nx+1,ny+1))
   allocate (uandv (nx+1,ny+1))  ! not used wwvv
   allocate (b     (nx+1,ny+1))  ! not used wwvv
   allocate (fslope(nx+1,ny+1))  ! not used wwvv
   allocate (hloc  (nx+1,ny+1))
   allocate (kb    (nx+1,ny+1))
   allocate (Ts    (nx+1,ny+1))
   allocate (ceq   (nx+1,ny+1))
   allocate (Ur    (nx+1,ny+1))
   allocate (Bm    (nx+1,ny+1))
   allocate (B1    (nx+1,ny+1))
endif
! Soulsby van Rijn sediment transport formula
! Ad Reniers april 2006
!
! z is defined positive upward 
! x is defined positive toward the shore
!  Formal parameters:
!  ------------------
!
!   Var. I/O  Type Dimensions
!   -------------------------
!

hloc   = max(hh,par%hmin) !Jaap par%hmin instead of par%eps
twothird=2.d0/3.d0
delta = (par%rhos-par%rho)/par%rho
! use eulerian velocities
! cjaap: add turbulence near bottom
do j=1,ny+1 
    do i=1,nx+1
       kb(i,j) = (DR(i,j)/par%rho)**twothird/ &
         (exp( min( hloc(i,j)/max(H(i,j),0.1d0) ,100.d0)) -1.d0)
       vmag2(i,j) = ue(i,j)**2+ve(i,j)**2
    enddo
enddo

vmag2  = ue**2+ve**2  ! wwvv todo just to be sure ?
urms2  = urms**2+0.50d0*kb

do jg = 1,par%ngd

   ! cjaap: compute fall velocity with simple expression from Ahrens (2000)
   Te    = 20.d0
   kvis  = 4.d0/(20.d0+Te)*1d-5 ! Van rijn, 1993 
   Sster = D50(jg)/(4*kvis)*sqrt((par%rhos/par%rho-1)*par%g*D50(jg))
   c1    = 1.06d0*tanh(0.064d0*Sster*exp(-7.5d0/Sster**2))
   c2    = 0.22d0*tanh(2.34d0*Sster**(-1.18d0)*exp(-0.0064d0*Sster**2))
   wster = c1+c2*Sster
   par%w = wster*sqrt((par%rhos/par%rho-1.d0)*par%g*D50(jg))

   Ts       = par%tsfac*hh/par%w
   !Ts       = max(Ts,0.2d0)
   Tsg(:,:,jg) = max(Ts,0.2d0) 

   dster=25296*D50(jg)
   
   ! calculate treshold velocity Ucr
   if(D50(jg)<=0.0005d0) then
     Ucr=0.19d0*D50(jg)**0.1d0*log10(4*hloc/D90(jg))
   else if(D50(jg)<0.002d0) then
     Ucr=8.5d0*D50(jg)**0.6d0*log10(4*hloc/D90(jg))
   else
#ifdef USEMPI
     write(*,'(a,i4)') 'In process',xmpi_rank
#endif
     write(*,*) '  Remark from sb_vr: D50(jg) > 2mm, out of validity range'
   end if
   ! drag coefficient
   z0 = par%z0
   Cd=(0.40d0/(log(max(hh,par%hmin)/z0)-1.0d0))**2 !Jaap
   !Cd = par%g/par%C**2;   ! consistent with flow modelling 

   ! Diane Foster and Robert: limit Shields to par%smax -> vmag2 for transp. limited
  ! vmag2=min(vmag2,par%smax*par%C**2*D50(jg)*delta)
   vmag2=min(vmag2,par%smax*par%g/par%cf*D50(jg)*delta)       ! In terms of cf

   ! transport parameters
   Asb=0.005d0*hloc*(D50(jg)/hloc/(delta*par%g*D50(jg)))**1.2d0     ! bed load coefficent
   Ass=0.012d0*D50(jg)*dster**(-0.6d0)/(delta*par%g*D50(jg))**1.2d0 ! suspended load coeffient
   term1=sqrt(vmag2+0.018d0/Cd*urms2)     ! nearbed-velocity
      
   term2 = 0
   do j=1,ny+1
      do i=1,nx
         if(term1(i,j)>Ucr(i,j) .and. hh(i,j)>par%eps) then
            term2(i,j)=(term1(i,j)-Ucr(i,j))**2.4d0
!           term2(i,j)=(term1(i,j)-Ucr(i,j))**(1.7-0.7*tanh((ue(i,j)/sqrt(par%g*hh(i,j))-0.5)*10))
         end if
      end do
   end do
   ! wwvv in parallel version, there will be a discrepancy between the values
   ! of term2. term2(nx+1,:) is zero, while the corresponding row in the process
   ! below term2(2,:) has some value, different from zero.
   ! so we fix this:
#ifdef USEMPI
   call xmpi_shift(term2,'m:')
#endif
   ceq =(Asb+Ass)*term2  
   ceq = min(ceq,0.2d0)       ! equilibrium concentration
   ceq = ceq/hloc
   ceqg(:,:,jg) = ceq*sedcal(jg)
enddo  ! end og grain size classes

m1 = 0;       ! a = 0
m2 = 0.7939;  ! b = 0.79 +/- 0.023
m3 = -0.6065; ! c = -0.61 +/- 0.041
m4 = 0.3539;  ! d = -0.35 +/- 0.032 
m5 = 0.6373;  ! e = 0.64 +/- 0.025
m6 = 0.5995;  ! f = 0.60 +/- 0.043

do j=1,ny+1     
   do i=1,nx+1
      if (k(i,j)*h(i,j)<par%px/2.d0 .and. H(i,j)>0.01d0) then
         Ur(i,j) = 3.d0/8.d0*sqrt(2.d0)*H(i,j)*k(i,j)/(k(i,j)*hloc(i,j))**3.d0                       !Ursell number
         Bm(i,j) = m1+(m2-m1)/(1.d0+exp((m3-log10(Ur(i,j)))/m4))                                     !Boltzmann sigmoid (eq 6) 
         B1(i,j) = -90.d0+90.d0*tanh(m5/Ur(i,j)**m6)    
         B1(i,j) = B1(i,j)*par%px/180.d0
         Sk(i,j) = Bm(i,j)*cos(B1(i,j))                                                              !Skewness (eq 8)
         As(i,j) = Bm(i,j)*sin(B1(i,j))                                                              !Skewness (eq 9)
		 ua(i,j) = par%facua*Sk(i,j)*urms(i,j)
	  endif
   enddo
enddo

end subroutine sb_vr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sednew(s,par)
use params
use spaceparams
use readkey_module
use xmpi_module

IMPLICIT NONE

type(spacepars),target                  :: s
type(parameters)                        :: par

character*80                            :: fnamet
integer                                 :: i,nw,ii
integer                                 :: j,jg
real*8                                  :: dster,onethird,twothird,Ass,dcf,ML
real*8                                  :: Te,kvis,Sster,cc1,cc2,wster,z0
real*8                                  :: ih0,it0,ih1,it1,p,q,f0,f1,f2,f3,uad,duddtmax,dudtmax,siguref,t0fac,duddtmean,dudtmean
real*8                               ,save     :: dh,dt,nh,nt
real*8 , dimension(:,:),allocatable  ,save     :: vmg,Asb,Ts
real*8 , dimension(:,:),allocatable  ,save     :: uorb,Ucr,Ucrc,Ucrw,term1,B2,Cd
real*8 , dimension(:,:),allocatable  ,save     :: hloc,ceq,h0,t0,detadxmax,detadxmean
real*8 , dimension(:,:,:),allocatable,save     :: RF 

include 's.ind'
include 's.inp'

if (.not. allocated(vmg)) then
   allocate (vmg   (nx+1,ny+1))
   allocate (h0    (nx+1,ny+1))
   allocate (t0    (nx+1,ny+1))
   allocate (detadxmax    (nx+1,ny+1))
   allocate (detadxmean   (nx+1,ny+1))
   allocate (term1 (nx+1,ny+1))
   allocate (B2    (nx+1,ny+1))
   allocate (Cd    (nx+1,ny+1))
   allocate (Asb   (nx+1,ny+1))
   allocate (Ucr   (nx+1,ny+1))
   allocate (Ucrc  (nx+1,ny+1))
   allocate (Ucrw  (nx+1,ny+1))
   allocate (uorb  (nx+1,ny+1))
   allocate (hloc  (nx+1,ny+1))
   allocate (kb    (nx+1,ny+1))
   allocate (Ts    (nx+1,ny+1))
   allocate (ceq   (nx+1,ny+1))
   allocate (RF    (18,33,40))
endif

hloc   = max(hh,par%hmin) !Jaap par%hmin instead of par%eps
onethird=1.d0/3.d0
twothird=2.d0/3.d0

!
! compute wave shape short waves (Rienecker and Fenton + Ruessink and van Rijn to estimate weighting of sine's and cosines)
!
! read table at t=0;
if (par%t==par%dt) then
   RF = RF*0.d0
   if (xmaster) then
      call readkey('params.txt','swtable',fnamet)
      open(31,file=fnamet);
      do i=1,18
         do j=1,33
            read(31,*)(RF(i,j,ii),ii=1,40)
         enddo
      enddo
   endif
   #ifdef USEMPI
   do i=1,18
      call xmpi_bcast(RF(i,:,:))
   enddo
   #endif
   dh = 0.03d0
   dt = 1.25d0
   nh = floor(0.99d0/dh);
   nt = floor(50.d0/dt);
endif
close(31)

hloc   = max(hh,par%hmin)

! read us and duddtmax from table....
h0 = min(nh*dh,max(dh,H/hloc))
t0 = min(nt*dt,max(dt,par%Trep*sqrt(par%g/hloc)))

do j=1,ny+1
   do i=1,nx+1
      ! interpolate table values....
      ih0=floor(h0(i,j)/dh);
      it0=floor(T0(i,j)/dt);
      ih1=min(ih0+1,nh);
      it1=min(it0+1,nt);
      p=(h0(i,j)-ih0*dh)/dh;
      q=(T0(i,j)-it0*dt)/dt;

      f0=(1-p)*(1-q);
      f1=p*(1-q);
      f2=q*(1-p);
      f3=p*q;
      
      Sk(i,j) = f0*RF(13,ih0,it0)+f1*RF(13,ih1,it0)+ f2*RF(13,ih0,it1)+f3*RF(13,ih1,it1)
	  As(i,j) = f0*RF(14,ih0,it0)+f1*RF(14,ih1,it0)+ f2*RF(14,ih0,it1)+f3*RF(14,ih1,it1)
      duddtmax = f0*RF(15,ih0,it0)+f1*RF(15,ih1,it0)+ f2*RF(15,ih0,it1)+f3*RF(15,ih1,it1)
	  siguref = f0*RF(16,ih0,it0)+f1*RF(16,ih1,it0)+ f2*RF(16,ih0,it1)+f3*RF(16,ih1,it1)
	  
      ! correct slope in case 1.25>T0>50
	  if (t0(i,j)==50.d0) then
	     t0fac = 50.d0/max((par%Trep*sqrt(par%g/hloc(i,j))),50.d0) 
	  elseif (t0(i,j)==1.25)then
	     t0fac = 1.25d0/min((par%Trep*sqrt(par%g/hloc(i,j))),1.25d0) 
	  else
	     t0fac = 1.d0
	  endif
	  ! translate dimensionless duddtmax to real world dudtmax
	  !         /scale with variance and go from [-] to [m/s^2]     /tableb./dimensionless dudtmax
      dudtmax = urms(i,j)/max(par%eps,siguref)*sqrt(par%g/hloc(i,j))*t0fac*duddtmax
	  detadxmax(i,j) = dudtmax*sinh(k(i,j)*hloc(i,j))/max(c(i,j),sqrt(H(i,j)*par%g))/sigm(i,j)
      
	  uad = f0*RF(17,ih0,it0)+f1*RF(17,ih1,it0)+ f2*RF(17,ih0,it1)+f3*RF(17,ih1,it1)
	  ! ua(i,j) = par%facua*uad*urms(i,j)/max(par%eps,siguref)
	  ! Jaap: Dano's approach shoudl be the same....
	  ua(i,j) = par%facua*Sk(i,j)*urms(i,j)

      ! Jaap: use average slope over bore front in roller energy balance...
	  duddtmean = f0*RF(18,ih0,it0)+f1*RF(18,ih1,it0)+ f2*RF(18,ih0,it1)+f3*RF(18,ih1,it1)
	  dudtmean = urms(i,j)/max(par%eps,siguref)*sqrt(par%g/hloc(i,j))*t0fac*duddtmean
	  detadxmean(i,j) = dudtmean*sinh(k(i,j)*hloc(i,j))/max(c(i,j),sqrt(H(i,j)*par%g))/sigm(i,j)
   enddo
enddo
Tbore = max(par%Trep/25.d0,min(par%Trep/4.d0,H/(max(c,sqrt(H*par%g))*max(detadxmax,par%eps))))
Tbore = par%Tbfac*Tbore
if (par%rfb==1) then
   BR = par%BRfac*sin(atan(detadxmean))
else
   BR = par%beta
endif
!
! compute near bed turbulence
!
do j=1,ny+1	
	do i=1,nx+1
	   ! compute mixing length
	   ! ML = 2*R(i,j)*par%Trep/(par%rho*c(i,j)*max(H(i,j),par%eps))
	   ML = dsqrt(2*R(i,j)*par%Trep/(par%rho*c(i,j)))
	   ! ML = 0.9d0*H(i,j)
	   ML = min(ML,hloc(i,j)+0.5d0*H(i,j));
	   ! exponential decay turbulence over depth
	   dcf = min(1.d0,1.d0/(exp(hloc(i,j)/max(ML,0.1d0)) -1.d0))
	   ! dcf = min(1.d0,1.d0/(cosh(hloc(i,j)/max(ML,0.1d0)) -1.d0))
	   if (par%turb == 2) then
          kb(i,j) = (DR(i,j)/par%rho)**twothird*dcf*par%Trep/Tbore(i,j)
	   elseif (par%turb == 1) then
	      kb(i,j) = (DR(i,j)/par%rho)**twothird*dcf
	   elseif (par%turb == 0) then
	      kb(i,j) = 0.0d0
	   endif
    enddo
enddo

vmg  = dsqrt(ue**2+ve**2)
uorb = dsqrt(urms**2.d0+1.45d0*kb)

do jg = 1,par%ngd

   ! Jaap: compute fall velocity with simple expression from Ahrens (2000)
   Te    = 20.d0
   kvis  = 4.d0/(20.d0+Te)*1d-5	! Van rijn, 1993 
   Sster = s%D50(jg)/(4*kvis)*sqrt((par%rhos/par%rho-1)*par%g*s%D50(jg))
   cc1   = 1.06d0*tanh(0.064d0*Sster*exp(-7.5d0/Sster**2))
   cc2   = 0.22d0*tanh(2.34d0*Sster**-1.18d0*exp(-0.0064d0*Sster**2))
   wster = cc1+cc2*Sster
   par%w = wster*sqrt((par%rhos/par%rho-1.d0)*par%g*s%D50(jg))
   
   !where (ceqg(:,:,jg)>ccg(:,:,jg))  Ts = hh/max(dsqrt(kb),par%w)
   !where (ceqg(:,:,jg)<=ccg(:,:,jg)) Ts = hh/par%w
   !Tsg(:,:,jg) = max(Ts,0.2d0)
   Ts       = par%tsfac*hh/par%w
   Tsg(:,:,jg) = max(Ts,par%Tsmin) 
   
   dster=25296*s%D50(jg)
   !
   ! calculate treshold velocity Ucr
   !
   if(s%D50(jg)<=0.0005) then
     Ucrc=0.19d0*s%D50(jg)**0.1d0*log10(4.d0*hloc/s%D90(jg)) 
	 Ucrw=0.24d0*(1.65d0*par%g)**0.66d0*s%D50(jg)**0.33d0*par%Trep**0.33d0 
   else if(s%D50(jg)<0.002) then
     Ucrc=8.5d0*s%D50(jg)**0.6d0*log10(4.d0*hloc/s%D90(jg))                          !Shields
	 Ucrw=0.95d0*(1.65d0*par%g)**0.57d0*s%D50(jg)**0.43*par%Trep**0.14               !Komar and Miller (1975)
   else
     write(*,*) ' s%D50(jg) > 2mm, out of validity range'
   end if
   B2 = vmg/max(vmg+uorb,par%eps)
   Ucr = B2*Ucrc + (1-B2)*Ucrw                                                       ! Van Rijn 2008 (Bed load transport paper)

   ! transport parameters
   Asb=0.015d0*hloc*(s%D50(jg)/hloc)**1.2d0/(1.65d0*par%g*s%D50(jg))**0.75d0         !bed load coefficent
   Ass=0.012d0*s%D50(jg)*dster**(-0.6d0)/(1.65d0*par%g*s%D50(jg))**1.2d0             !suspended load coeffient
   
   ! Jaap: par%sws to set short wave stirring to zero
   ! Jaap: Van Rijn use Peak orbital flow velocity --> 0.64 correpsonds to 0.4 coefficient regular waves Van Rijn (2007)  
   term1= dsqrt(vmg**2+0.64d0*par%sws*uorb**2)                                       

   ! Try Soulsby van rijn approach...
   ! drag coefficient
   ! z0 = par%z0
   ! Cd=(0.40/(log(max(hh,par%hmin)/z0)-1.0))**2 !Jaap
   ! term1=(vmg**2+0.018/Cd*uorb**2)**0.5   
      
   ceq = 0*term1                                                                     !initialize ceq
   do j=1,ny+1
	  do i=1,nx
         if(term1(i,j)>Ucr(i,j) .and. hh(i,j)>par%eps) then
            ceq(i,j)=Asb(i,j)*(term1(i,j)-Ucr(i,j))**1.5 + Ass*(term1(i,j)-Ucr(i,j))**2.4
		 end if
      end do
   end do
  
   ceq = min(ceq,0.05)		      ! equilibrium concentration
   ceq = ceq/hloc
   ceqg(:,:,jg) = ceq*sedcal(jg)*wetz
enddo                             ! end of grain size classes
                            ! end of grain size classes

end subroutine sednew

end module morphevolution
