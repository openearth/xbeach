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
use interp
! use vsmumod

IMPLICIT NONE

type(spacepars),target                   :: s
type(parameters)                         :: par

integer                                  :: i,isig
integer                                  :: j,jg

real*8,dimension(:),allocatable,save     :: chain,cumchain
real*8,dimension(:,:),allocatable,save   :: vmag2,uau,uav,um,vm
real*8,dimension(:,:),allocatable,save   :: ccvt,dcdz,dsigt,aref
real*8,dimension(:,:),allocatable,save   :: cc,cu,cv,Sus,Svs,Dc
real*8,dimension(:,:),allocatable,save   :: ccb,cub,cvb,Sub,Svb,pbbedu,pbbedv
real*8,dimension(:,:),allocatable,save   :: suq3d,svq3d,eswmax,eswbed,sigs,deltas
real*8,dimension(:,:,:),allocatable,save :: dsig,ccv,sdif,cuq3d,cvq3d

include 's.ind'
include 's.inp'

if (.not. allocated(vmag2)) then
   allocate(vmag2 (nx+1,ny+1))
   allocate(uau (nx+1,ny+1))
   allocate(uav (nx+1,ny+1))
   allocate(cu  (nx+1,ny+1))
   allocate(cv  (nx+1,ny+1))
   allocate(cc  (nx+1,ny+1))
   allocate(Sus (nx+1,ny+1))
   allocate(Svs (nx+1,ny+1))
   allocate(cub (nx+1,ny+1))
   allocate(cvb (nx+1,ny+1))
   allocate(ccb (nx+1,ny+1))
   allocate(Sub (nx+1,ny+1))
   allocate(Svb (nx+1,ny+1))
   allocate(Dc  (nx+1,ny+1))
   allocate(pbbedu (nx+1,ny+1))
   allocate(pbbedv (nx+1,ny+1))
   allocate(ccvt (nx+1,ny+1))
   allocate(dcdz (nx+1,ny+1))
   allocate(dsigt (nx+1,ny+1))
   allocate(dsig(s%nx+1,s%ny+1,par%kmax))
   allocate(ccv(s%nx+1,s%ny+1,par%kmax))
   allocate(sdif(s%nx+1,s%ny+1,par%kmax))
   allocate(um (nx+1,ny+1))
   allocate(vm (nx+1,ny+1))
   allocate(deltas(nx+1,ny+1))
   allocate(sigs(nx+1,ny+1))
   allocate(eswmax(nx+1,ny+1))
   allocate(eswbed(nx+1,ny+1))
   allocate(suq3d(nx+1,ny+1))
   allocate(svq3d(nx+1,ny+1))
   allocate(cuq3d(nx+1,ny+1,par%kmax))
   allocate(cvq3d(nx+1,ny+1,par%kmax))
   allocate(aref(nx+1,ny+1))
   allocate(chain(par%kmax))
   allocate(cumchain(par%kmax))
  
   uau       = 0.d0
   ccg       = 0.d0
   uav       = 0.d0
   Dc        = 0.d0
   um        = 0.d0
   vm        = 0.d0
   chain     = 0.0d0
   cumchain  = 0.0d0
   ! generate sigma grid shape...
   do isig=2,par%kmax
	  chain(isig) = par%sigfac**(isig-1)
	  cumchain(isig) = cumchain(isig-1)+chain(isig)
   enddo
endif
! use eulerian velocities
vmag2     = ue**2+ve**2
cu        = 0.0d0
cv        = 0.0d0
cub       = 0.0d0
cvb       = 0.0d0 	
dcsdx     = 0.0d0
dcsdy     = 0.0d0
dcbdx     = 0.0d0
dcbdy     = 0.0d0

! compute long wave turbulence due to breaking
if (par%lwt==1) then
   call longwaveturb(s,par)
endif

! calculate equilibrium concentration
if (par%form==1) then           ! Soulsby van Rijn
   call sb_vr(s,par)
elseif (par%form==2) then       ! Van Thiel de Vries & Reniers 2008
   call sednew(s,par)
end if

! compute long wave turbulence due to breaking
if (par%lwt==1) then
   call longwaveturb(s,par)
endif

! compute diffusion coefficient

!if (par%nuhfac==1) then
   Dc = par%nuh+par%nuhfac*hh*(DR/par%rho)**(1.d0/3.d0)
!end if

dzbdt=0.0d0

do jg = 1,par%ngd
   cc = ccg(:,:,jg)
   if (D50(jg)>0.002d0) then
      ! RJ: set ceqsg to zero for gravel.
      cc = 0.d0 ! Can be used to test total transport mode
   endif 
   !
   ! X-direction
   do j=1,ny+1
      do i=1,nx
         if(ueu(i,j)>0.d0) then
            ! test cu(i,j)=cc(i,j)
            cu(i,j)=par%thetanum*cc(i,j)+(1.d0-par%thetanum)*cc(min(i+1,nx),j)
			cub(i,j)=par%thetanum*ceqbg(i,j,jg)+(1.d0-par%thetanum)*ceqbg(min(i+1,nx),j,jg)
         elseif(ueu(i,j)<0.d0) then
            cu(i,j)=par%thetanum*cc(i+1,j)+(1.d0-par%thetanum)*cc(max(i,2),j)
			cub(i,j)=par%thetanum*ceqbg(i+1,j,jg)+(1.d0-par%thetanum)*ceqbg(max(i,2),j,jg)
         else
            cu(i,j)=0.5d0*(cc(i,j)+cc(i+1,j))
			cub(i,j)=0.5d0*(ceqbg(i,j,jg)+ceqbg(i+1,j,jg))
         endif
         dcsdx(i,j)=(cc(i+1,j)-cc(i,j))/(xz(i+1)-xz(i))
		 dcbdx(i,j)=(ceqbg(i+1,j,jg)-ceqbg(i,j,jg))/(xz(i+1)-xz(i))
      enddo 
   enddo
   ! wwvv dcdx(nx:1,:) is still untouched, correct this ofr the parallel case
#ifdef USEMPI
   call xmpi_shift(dcsdx,'m:')
   call xmpi_shift(dcbdx,'m:')
#endif
   cu(nx+1,:) = cc(nx+1,:) !Robert
   ! wwvv fix this in parallel case
#ifdef USEMPI
   call xmpi_shift(cu,'m:')
#endif
   ! wwvv in parallel version, there will be a discrepancy between the values
   ! of dzbdx(nx+1,:).
   !wwvv so fix that
#ifdef USEMPI
   call xmpi_shift(dzbdx,'m:')
#endif
   ! Get ua in u points and split out in u and v direction
   uau(1:nx,:) = 0.5*(ua(1:nx,:)*cos(thetamean(1:nx,:))+ua(2:nx+1,:)*cos(thetamean(1:nx,:)))
   uav(1:nx,:) = 0.5*(ua(1:nx,:)*sin(thetamean(1:nx,:))+ua(2:nx+1,:)*sin(thetamean(1:nx,:)))
   ! Compute vmagu including ua
   vmagu = sqrt((uu+uau)**2+(vu+uav)**2)
   !
   ! Compute sedimnent transport in u-direction
   !
   ureps = ueu
   urepb = ueu   ! RJ maybe reduce this velocity?
   !
   Sus = 0.d0
   Sub = 0.d0
   ! suspended load
   Sus=(cu*(ureps+uau)*hu-Dc*hu*dcsdx)*wetu   !No bed slope term in suspended transport?
   ! bed load
   Sub=(cub*(urepb+uau)*hu-par%facsl*cub*vmagu*hu*dzbdx)*wetu 
   ! 
   do j=1,ny+1
      do i=1,nx
         if(Sub(i,j)>0.d0) then
			pbbedu(i,j) = pbbed(i,j,1,jg)
         elseif(Sub(i,j)<0.d0) then
			pbbedu(i,j)= pbbed(i+1,j,1,jg)
         else
			pbbedu(i,j)=0.5d0*(pbbed(i,j,1,jg)+pbbed(i+1,j,1,jg))
         endif
      enddo 
   enddo
   !
   Sub = pbbedu*Sub
   !
   ! Y-direction (2DH & Q3D)
   !
   do j=1,ny
      do i=1,nx+1
         if(vev(i,j)>0) then
            ! cv(i,j)=cc(i,j)
            cv(i,j)=par%thetanum*cc(i,j)+(1.d0-par%thetanum)*cc(i,min(j+1,ny))
			cvb(i,j)=par%thetanum*ceqbg(i,j,jg)+(1.d0-par%thetanum)*ceqbg(i,min(j+1,ny),jg)
         else if(vev(i,j)<0) then
            ! cv(i,j)=cc(i,j+1)
            cv(i,j)=par%thetanum*cc(i,j+1)+(1.d0-par%thetanum)*cc(i,max(j,2))
			cvb(i,j)=par%thetanum*ceqbg(i,j+1,jg)+(1.d0-par%thetanum)*ceqbg(i,max(j,2),jg)
         else
            cv(i,j)=0.5d0*(cc(i,j)+cc(i,j+1)) !Jaap: cc instead of cv
			cvb(i,j)=0.5d0*(ceqbg(i,j,jg)+ceqbg(i,j+1,jg))
         end if
         dcsdy(i,j)=(cc(i,j+1)-cc(i,j))/(yz(j+1)-yz(j)) !Jaap
		 dcbdy(i,j)=(ceqbg(i,j+1,jg)-ceqbg(i,j,jg))/(yz(j+1)-yz(j)) !Jaap
      end do 
    end do
    ! wwvv dcdy(:,ny+1) is not filled in, so in parallel case:
#ifdef USEMPI
    call xmpi_shift(dcsdy,':n')
	call xmpi_shift(dcbdy,':n')
#endif
    cv(:,ny+1) = cc(:,ny+1) !Robert
   ! wwvv in parallel version, there will be a discrepancy between the values
   ! of dzbdy(:,ny+1).
   ! wwvv so fix that
#ifdef USEMPI
    call xmpi_shift(dzbdy,':n')
#endif
	! Jaap: get ua in v points and split out in u and v direction
    uau(:,1:ny) = 0.5*(ua(:,1:ny)*cos(thetamean(:,1:ny))+ua(:,2:ny+1)*cos(thetamean(:,1:ny)))
    uav(:,1:ny) = 0.5*(ua(:,1:ny)*sin(thetamean(:,1:ny))+ua(:,2:ny+1)*sin(thetamean(:,1:ny)))
    ! Jaap: compute vmagv including ua
    vmagv = sqrt((uv+uau)**2+(vv+uav)**2)
    !
    ! Compute sedimnent transport in v-direction
    !
    vreps = vev
	vrepb = vev   ! RJ maybe reduce this velocity?
    !
    Svs = 0.d0
    Svb = 0.d0
    ! Suspended load
    Svs=(cv*(vreps+uav)*hv-Dc*hv*dcsdy)*wetv
	! Bed load
	Svb=(cvb*(vrepb+uav)*hv-par%facsl*cvb*vmagv*hv*dzbdy)*wetv
	!
	do j=1,ny
      do i=1,nx+1
         if(Svb(i,j)>0) then
			pbbedv(i,j)=pbbed(i,j,1,jg)
         else if(Svb(i,j)<0) then
			pbbedv(i,j)=pbbed(i,j+1,1,jg)
         else
			pbbedv(i,j)=0.5d0*(pbbed(i,j,1,jg)+pbbed(i,j+1,1,jg))
         end if
      end do 
    end do 
	!
	Svb = 0.d0
	Svb = pbbedv*Svb
    !
    !BRJ: compute sources (sink must be computed after updating actual sed.conc.)
    !    
    do j=2,ny
      do i=2,nx
         ero(i,j,jg) = hh(i,j)*ceqsg(i,j,jg)*pbbed(i,j,1,jg)/Tsg(i,j,jg)    ! Changed to hh from hold by RJ (13072009)
         ! depo_ex(i,j,jg) = max(hold(i,j),0.01d0)*cc(i,j)/Tsg(i,j,jg)                    
	     !BRJ: the volume in the water column is updated and not the volume concentration.
		 !BRJ: implicit concentration update
         cc(i,j) = (par%dt*Tsg(i,j,jg))/(par%dt+Tsg(i,j,jg))* &
		                                      (hold(i,j)*cc(i,j)/par%dt -((Sus(i,j)-Sus(i-1,j))/(xu(i)-xu(i-1))+&
                                                                          (Svs(i,j)-Svs(i,j-1))/(yv(j)-yv(j-1))-&
		 									                              ero(i,j,jg)))	

         cc(i,j)=max(cc(i,j),0.0d0) ! Jaap: negative cc's are possible...
		 depo_ex(i,j,jg) = cc(i,j)/Tsg(i,j,jg)                    																		  																	  													
      enddo
    enddo

	! cc = min(1.d0,cc/max(hh,0.01d0))
	cc = cc/hh
	! do lateral bounadries...
if(xmpi_istop)then
    cc(1,:)=cc(2,:)
	ero(1,:,jg)=ero(2,:,jg)
    depo_ex(1,:,jg)=depo_ex(2,:,jg)
endif
if(xmpi_isleft)then
	cc(:,1)=cc(:,2)
    ero(:,1,jg)=ero(:,2,jg)
    depo_ex(:,1,jg)=depo_ex(:,2,jg)
endif
if(xmpi_istop)then
	cc(nx+1,:)=cc(nx+1-1,:)
    ero(nx+1,:,jg)=ero(nx,:,jg)
    depo_ex(nx+1,:,jg)=depo_ex(nx,:,jg)
endif
if(xmpi_isright)then
	cc(:,ny+1)=cc(:,ny+1-1)
	ero(:,ny+1,jg)=ero(:,ny,jg)
    depo_ex(:,ny+1,jg)=depo_ex(:,ny,jg)
endif
	
    ! wwvv fix the first and last rows and columns of cc in parallel case
#ifdef USEMPI
    call xmpi_shift(cc,'1:')
    call xmpi_shift(cc,'m:')
    call xmpi_shift(cc,':1')   ! Maybe not necessary as zb calculated from 1:ny+1, but just in case...
    call xmpi_shift(cc,':n')   ! Dito
#endif
    ! Jaap
    ! cc=cc*wetz
    !
    ! wwvv border columns and rows of ccg Svg and Sug have to be communicated
    ccg(:,:,jg) = cc
    Svsg(:,:,jg) = Svs
    Susg(:,:,jg) = Sus
	Svbg(:,:,jg) = Svb
    Subg(:,:,jg) = Sub

end do ! number of sediment fractions

vmag=sqrt(max(vmag2,par%umin))

end subroutine transus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bed_update(s,par)
use params
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars),target              :: s
type(parameters)                    :: par

integer                             :: i,j,ii,iii,jjj,jdz,ndz,ind
integer                             :: jg,nt_e,nt_d,nt_sub,nt_sub_max
real*8                              :: dzb,dzmax,dzt,dzleft
real*8 , dimension(:,:),allocatable,save :: dzbtot
real*8 , dimension(par%ngd)         :: edg,edg1,edg2,dzg
real*8 , dimension(:),pointer       :: dz
real*8 , dimension(:,:),pointer     :: pb
logical                             :: aval

include 's.ind'
include 's.inp'

if (.not. allocated(dzbtot)) allocate(dzbtot(s%nx+1,s%ny+1))

dzbtot = 0.d0
dzbdt  = 0.d0
dzb    = 0.d0

if (par%t>=par%morstart .and. par%morfac > .999d0) then
!
par%frac_dz = 0.7d0 ! relative thickness to split time step for bed updating
par%split = 1.01d0
par%merge = 0.01d0

nt_sub=1
nt_sub_max = nt_sub	
!
! bed_predict
!
do j=2,ny
   do i=2,nx
      dz=>dzbed(i,j,:)
	  pb=>pbbed(i,j,:,:)
	  ! bed level changes per fraction in this morphological time step in meters sand including pores
	  if (par%sourcesink==0) then
	  dzg=par%morfac*par%dt/(1.d0-par%por)*( & ! Dano, dz from sus transport gradients      
											 (Susg(i,j,:)-Susg(i-1,j,:))/(xu(i)-xu(i-1)) +&           
		                                     (Svsg(i,j,:)-Svsg(i,j-1,:))/(yv(j)-yv(j-1)) +&        	                                                                                                
											 ! dz from bed load transport gradients
                                             (Subg(i,j,:)-Subg(i-1,j,:))/(xu(i)-xu(i-1)) +&           
		                                     (Svbg(i,j,:)-Svbg(i,j-1,:))/(yv(j)-yv(j-1))    )        
      elseif (par%sourcesink==1) then
	  dzg=par%morfac*par%dt/(1.d0-par%por)*( & ! BRJ, dz from source-sink terms
	                                         ero(i,j,:)-depo_ex(i,j,:)                   +&        	                                                                                                
											 ! dz from bed load transport gradients
                                             (Subg(i,j,:)-Subg(i-1,j,:))/(xu(i)-xu(i-1)) +&           
		                                     (Svbg(i,j,:)-Svbg(i,j-1,:))/(yv(j)-yv(j-1))    )      
	  endif
      
	  ind = minval(minloc(dz))
      nt_e=max(1,ceiling(maxval(dzg(:)/(par%frac_dz*dz(ind)*max(pb(ind,:),0.001d0)) ) ) )

      nt_d=max(1,ceiling(abs(sum(dzg))/par%frac_dz/minval(dz)))
	  nt_sub=max(nt_e,nt_d)
	  nt_sub_max = max(nt_sub_max,nt_sub)

      !erosion/deposition rate of sand mass (m/s)
	  !positive in case of erosion
		 
	  edg = dzg*(1.d0-par%por)/par%dt           
   	 
	  call update_fractions(par,s,i,j,nt_sub,dz,pb,edg,0.d0)
 
   enddo ! nx+1
enddo ! ny+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Avalanching
!
do ii=1,nint(par%morfac)
   
   aval=.false.
   dzbdx=0.d0
   dzbdy=0.d0
   do j=1,ny+1
      do i=1,nx
         dzbdx(i,j)=(zb(i+1,j)-zb(i,j))/(xz(i+1)-xz(i))
      enddo
   enddo
   !
   do i=2,nx-1
      do j=1,ny+1
         if(max(hh(i,j)+par%delta*H(i,j),hh(i+1,j)+par%delta*H(i+1,j))>par%hswitch+par%eps) then
            dzmax=par%wetslp;
         else
            dzmax=par%dryslp;
         end if
         if(abs(dzbdx(i,j))>dzmax ) then
		    aval=.true.     
            dzb=sign(1.0d0,dzbdx(i,j))*(abs(dzbdx(i,j))-dzmax)*(xz(i+1)-xz(i));
			!
			if (dzb >= 0.d0) then
               iii = i+1                                      ! upwind erosion
			   dzb=min(dzb,par%dzmax*par%dt/(xz(i+1)-xz(i)))
			else
			   iii = i                                        ! upwind deposition
			   dzb=max(dzb,-par%dzmax*par%dt/(xz(i+1)-xz(i)))
			endif
            dz => dzbed(iii,j,:) 
	 		pb => pbbed(iii,j,:,:)
			dzleft = abs(dzb)
			ndz = max(1,1+ceiling(max(dzleft-dz(1),0.d0)/dz(2))) ! check number of affected layers by avalanching
			do jdz=1,ndz
				dzt = min(dz(jdz),dzleft)
			   dzleft = dzleft-dzt;
               do jg=1,par%ngd ! erosion deposition per fraction for the eroding point (upwind or downwind); edg is positive in case of erosion   
			      edg2(jg) = sedcal(jg)*sign(1.d0,dzb)*dzt*pb(jdz,jg)*(1.d0-par%por)/par%dt       
			      edg1(jg) = -sedcal(jg)*edg2(jg)*(xu(i+1)-xu(i))/(xu(i)-xu(i-1))
			   enddo
			   nt_d=max(1,ceiling(abs(dzb)/par%frac_dz/minval(dzbed(iii-sign(1,ceiling(dzb)),j,:)))) ! nt_d is not necessarily equal to nt_e (=ndz)
			   nt_sub = max(1,nt_d)

			   ! jaap: update dzb
			   dzb = sum(edg2)/(1.d0-par%por)*par%dt
               
			   dz=>dzbed(i+1,j,:)
	           pb=>pbbed(i+1,j,:,:)
               call update_fractions(par,s,i+1,j,nt_sub,dz,pb,edg2,-dzb) ! upwind point

			   dz=>dzbed(i,j,:)
	           pb=>pbbed(i,j,:,:)
	   	       call update_fractions(par,s,i,j,nt_sub,dz,pb,edg1,dzb*(xu(i+1)-xu(i))/(xu(i)-xu(i-1)))    ! downwind point
			enddo

		 end if
      end do
   end do
   !JJ: update y slopes after avalanching in X-direction seems more appropriate
   do j=1,ny
      do i=1,nx+1
         dzbdy(i,j)=(zb(i,j+1)-zb(i,j))/(yz(j+1)-yz(j))
      enddo
   enddo

   do j=2,ny-1
        do i=1,nx+1
            if(max(hh(i,j),hh(i,j+1))>par%hswitch+par%eps) then
                dzmax=par%wetslp
            else
                dzmax=par%dryslp
            end if
            if(abs(dzbdy(i,j))>dzmax ) then 
		    aval=.true. 
            dzb=sign(1.0d0,dzbdy(i,j))*(abs(dzbdy(i,j))-dzmax)*(yz(j+1)-yz(j));
			!
			if (dzb >= 0.d0) then
               jjj = j+1
			   dzb=min(dzb,par%dzmax*par%dt/(yz(j+1)-yz(j)))
			else
			   jjj = j;
			   dzb=max(dzb,-par%dzmax*par%dt/(yz(j+1)-yz(j)))
			endif
            
			dz => dzbed(i,jjj,:) 
			pb => pbbed(i,jjj,:,:)

			dzleft = abs(dzb)
			ndz = max(1,1+ceiling(max(dzleft-dz(1),0.d0)/dz(2))) ! check number of affected layers by avalanching

			do jdz=1,ndz
			   dzt = min(dz(jdz),dzleft)
			   dzleft = dzleft-dzt;
               do jg=1,par%ngd ! erosion deposition per fraction for the eroding point (upwind or downwind)   
			      edg2(jg) = sedcal(jg)*sign(1.d0,dzb)*dzt*pb(jdz,jg)*(1.d0-par%por)/par%dt       ! upwind sand ends up in down wind point --> use pbbed(i+1)
			      edg1(jg) = -sedcal(jg)*edg2(jg)*(yv(j+1)-yv(j))/(yv(j)-yv(j-1))
			   enddo
			   
			   nt_d=max(1,ceiling(abs(dzb)/par%frac_dz/minval(dzbed(i,jjj-sign(1,ceiling(dzb)),:)))) ! nt_d is not necessarily equal to nt_e (=ndz)
			   nt_sub = max(1,nt_d)
               
               ! jaap: update dzb
			   dzb = sum(edg2)/(1.d0-par%por)*par%dt

			   dz=>dzbed(i,j+1,:)
	           pb=>pbbed(i,j+1,:,:)

               call update_fractions(par,s,i,j+1,nt_sub,dz,pb,edg2,-dzb) ! upwind point

			   dz=>dzbed(i,j,:)
	           pb=>pbbed(i,j,:,:)

	   	       call update_fractions(par,s,i,j,nt_sub,dz,pb,edg1,dzb*(yv(j+1)-yv(j))/(yv(j)-yv(j-1)))    ! downwind point
			enddo
		 end if
      end do
   end do
   if (.not.aval) exit
end do
!
! bed boundary conditions
! 
if(xmpi_isleft) then
   zb(:,1) = zb(:,2)
   sedero(:,1) = sedero(:,2)
   pbbed(:,1,:,:)=pbbed(:,2,:,:)
   z0bed(:,1)=z0bed(:,2)
   dzbed(:,1,:)=dzbed(:,2,:)
endif

if(xmpi_isright) then
   zb(:,ny+1) = zb(:,ny)
   sedero(:,ny+1) = sedero(:,ny)
   pbbed(:,ny+1,:,:)=pbbed(:,ny,:,:)
   z0bed(:,ny+1)=z0bed(:,ny)
   dzbed(:,ny+1,:)=dzbed(:,ny,:)
endif

! Robert: in parallel version bed update must take place on internal boundaries:
#ifdef USEMPI
    call xmpi_shift(zb,'1:')
    call xmpi_shift(zb,'m:')
    call xmpi_shift(zb,':1')
    call xmpi_shift(zb,':n')
#endif

! Update representative sed.diameter at the bed for flow friction and output
if (par%ngd>1) then
   do j=2,s%ny
      do i=2,s%nx
         s%D50top(i,j) =  sum(s%pbbed(i,j,1,:)*s%D50)
         s%D90top(i,j) =  sum(s%pbbed(i,j,1,:)*s%D90)
      enddo
   enddo
endif

endif ! if par%t>par%morstart

end subroutine bed_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_fractions(par,s,i,j,nt_sub,dz,pb,edg,dzb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Copyright (C) 2009 Technische Universiteit Delft
!        Bram van Prooijen
!        b.c.vanprooijen@tudelft.nl
!	   +31(0)15 2784070   
!        Faculty of Civil Engineering and Geosciences
!        department of Hydraulic Engineering
!	   PO Box 5048
!        2600 GA Delft
!        The Netherlands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use params
use spaceparams
use xmpi_module


IMPLICIT NONE

type(spacepars),target                          :: s
type(parameters)                                :: par

integer                                         :: i,j,jg,jd,t_sub,nt_sub
real*8                                          :: ED,zbold,dzb
real*8 , dimension(par%ngd),intent(in)          :: edg
real*8 , dimension(s%nd(i,j)),intent(inout)     :: dz
real*8 , dimension(s%nd(i,j),par%ngd),intent(inout) :: pb

real*8 , dimension(:),allocatable,save                :: Ap,b
real*8 , dimension(:,:),allocatable,save              :: Sm,A

if (.not. allocated(Ap)) then
   allocate(Ap   (s%nd(1,1)))
   allocate(b    (s%nd(1,1)))
   allocate(Sm   (s%nd(1,1),par%ngd))
   allocate(A    (s%nd(1,1),3))
endif

!!!initialize Sm

ED = sum(edg)

do t_sub=1,nt_sub  !loop over subtimesteps
  do jg=1,par%ngd

   A=0.
   Ap=0.
   b=0.

   Sm(:,jg)=pb(:,jg)*dz*(1.d0-par%por)
				
   !!!build matrix A							
   select case(par%nd_var)
   case(1)
   !in this case: A=0
   case(2)
     A (1,1:3)= (/0.d0           ,  min(ED,0.d0) ,  max(ED,0.d0) /)   
     A (2,1:3)= (/-min(ED,0.d0)  , -max(ED,0.d0) ,  0.d0         /)
   case(3:1000)
     A (1,1:3)= (/0.d0           ,  min(ED,0.d0) ,  max(ED,0.d0) /)   

     A (2:par%nd_var-1,1)=-min(ED,0.d0)
     A (2:par%nd_var-1,2)=-abs(ED)
     A (2:par%nd_var-1,3)=max(ED,0.d0)

     A (par%nd_var,1:3)= (/-min(ED,0.d0) , -max(ED,0.d0) ,   0.d0    /)
   end select

   !!!determine RHS
   Ap(1) = sum(A(1,2:3)*pb(1:2,jg))
   do jd = 2,par%nd_var
      Ap(jd) = sum(A (jd,:)*pb(jd-1:jd+1,jg))
   enddo
   Ap(par%nd_var+1) = sum(A(par%nd_var+1,1:2)*pb(par%nd_var:par%nd_var+1,jg))			

   b(1) = -edg(jg)

   !!!update S
   Sm(1:par%nd_var+1,jg) = Sm(1:par%nd_var+1,jg) + par%dt/nt_sub*(Ap(1:par%nd_var+1)+b(1:par%nd_var+1))			

  enddo !fractions

  do jd=1,par%nd_var+1
   pb(jd,:) = Sm(jd,:)/sum(Sm(jd,:))
   dz(jd) = sum(Sm(jd,:))/(1.d0-par%por)
  enddo

  !!! modify grid
  !merge two upper layers in case of erosion
  if (dz(par%nd_var) .lt. par%merge*dz(par%nd_var+1)) then
    forall (jg=1:par%ngd)
      pb(par%nd_var,jg) = (dz(par%nd_var)*pb(par%nd_var,jg) + dz(par%nd_var+1)* &
	                      pb(par%nd_var+1,jg))/(dz(par%nd_var)+dz(par%nd_var+1)) 
	  pb(par%nd_var+1:s%nd(i,j)-1,jg) = pb(par%nd_var+2:s%nd(i,j),jg)
	  pb(s%nd(i,j),jg) = pb(s%nd(i,j),jg)	
    endforall
    s%z0bed(i,j) = s%z0bed(i,j)-dz(par%nd_var+1)
    dz(par%nd_var) = dz(par%nd_var+1)+dz(par%nd_var)   
  endif
  !split upper layer in case of sedimentation
  if (dz(par%nd_var)>par%split*dz(par%nd_var+1)) then
    pb(par%nd_var+1:s%nd(i,j),:) = pb(par%nd_var:s%nd(i,j)-1,:)
    s%z0bed(i,j) = s%z0bed(i,j)+dz(par%nd_var+1)
    dz(par%nd_var) = dz(par%nd_var)-dz(par%nd_var+1)  
  endif
enddo ! nt_sub

pb = max(0.d0,min(pb,1.d0))

zbold = s%zb(i,j)
s%zb(i,j) = s%z0bed(i,j)+sum(dz)
s%dzbdt(i,j) = s%dzbdt(i,j)+(s%zb(i,j)-zbold)
s%sedero(i,j) = s%sedero(i,j)+s%dzbdt(i,j)
if (abs(dzb)>0.d0) then 
  s%zs(i,j)=s%zs(i,j)+(s%zb(i,j)-zbold)
  s%dzav(i,j)=s%dzav(i,j)+dzb
endif


end subroutine update_fractions

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
real*8                              :: onethird,twothird
real*8                              :: z0,Ass,dcf,dcfin
real*8                              :: Te,kvis,Sster,c1,c2,wster
real*8,save                         :: m1,m2,m3,m4,m5,m6
real*8,save                         :: alpha,beta,delta

real*8 , dimension(:)  ,allocatable,save   :: w,dster
real*8 , dimension(:,:),allocatable,save   :: vmg,Cd,Asb,dhdx,dhdy,Ts,Ur,Bm,B1
real*8 , dimension(:,:),allocatable,save   :: urms2,Ucr,term1,term2
real*8 , dimension(:,:),allocatable,save   :: uandv,b,fslope,hloc,ceqs,ceqb

include 's.ind'
include 's.inp'

if (.not. allocated(vmg)) then
   allocate (vmg   (nx+1,ny+1))
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
   allocate (ceqs   (nx+1,ny+1))
   allocate (ceqb   (nx+1,ny+1))
   allocate (Ur    (nx+1,ny+1))
   allocate (Bm    (nx+1,ny+1))
   allocate (B1    (nx+1,ny+1))
   allocate (w     (par%ngd))
   allocate (dster (par%ngd))
   vmg = 0.d0
   m1 = 0;       ! a = 0
   m2 = 0.7939;  ! b = 0.79 +/- 0.023
   m3 = -0.6065; ! c = -0.61 +/- 0.041
   m4 = 0.3539;  ! d = -0.35 +/- 0.032 
   m5 = 0.6373;  ! e = 0.64 +/- 0.025
   m6 = 0.5995;  ! f = 0.60 +/- 0.043
   alpha = -log10(exp(1.d0))/m4
   beta  = exp(m3/m4)
   onethird=1.d0/3.d0
   ! Robert: do only once, not necessary every time step
   do jg=1,par%ngd
      ! cjaap: compute fall velocity with simple expression from Ahrens (2000)
      Te    = 20.d0
      kvis  = 4.d0/(20.d0+Te)*1d-5 ! Van rijn, 1993 
      Sster = D50(jg)/(4*kvis)*sqrt((par%rhos/par%rho-1)*par%g*D50(jg))
      c1    = 1.06d0*tanh(0.064d0*Sster*exp(-7.5d0/Sster**2))
      c2    = 0.22d0*tanh(2.34d0*Sster**(-1.18d0)*exp(-0.0064d0*Sster**2))
      wster = c1+c2*Sster
      w(jg) = wster*sqrt((par%rhos/par%rho-1.d0)*par%g*D50(jg))
	  ! RJ: for modeling gravel
	  delta = (par%rhos-par%rho)/par%rho
      dster(jg)=(delta*par%g/1.d-12)**onethird*s%D50(jg) 
   enddo
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
hloc = hh
twothird=2.d0/3.d0
! use eulerian velocities
! cjaap: add turbulence near bottom
do j=1,ny+1 
    do i=1,nx+1
	   ! exponential decay turbulence over depth
	   dcfin = exp(min(100.d0,hloc(i,j)/max(H(i,j),0.1d0)))
	   dcf = min(1.d0,1.d0/(dcfin-1.d0))
       kb(i,j) = par%nuhfac*(DR(i,j)/par%rho)**twothird*dcf
    enddo
enddo

if (par%lws==1) then
   vmg  = dsqrt(ue**2+ve**2)
elseif (par%lws==0) then
   ! vmg lags on actual mean flow; but long wave contribution to mean flow is included... 
   vmg = (1.d0-1.d0/par%cats/par%Trep*par%dt)*vmg + (1.d0/par%cats/par%Trep*par%dt)*dsqrt(ue**2+ve**2) 
endif

urms2  = urms**2+0.50d0*(kb+kturb)


do jg = 1,par%ngd

   Ts       = par%tsfac*hloc/w(jg)
   Tsg(:,:,jg) = max(Ts,par%Tsmin) 
   
   ! calculate treshold velocity Ucr
   if(D50(jg)<=0.0005d0) then
     Ucr=0.19d0*D50(jg)**0.1d0*log10(4.d0*hloc/D90(jg))
   else if(D50(jg)<0.05d0) then   !Dano see what happens with coarse material
     Ucr=8.5d0*D50(jg)**0.6d0*log10(4*hloc/D90(jg))
   else
#ifdef USEMPI
     write(*,'(a,i4)') 'In process',xmpi_rank
#endif
   end if
   ! drag coefficient
   z0 = par%z0
   Cd=(0.40d0/(log(max(hloc,10.d0*z0)/z0)-1.0d0))**2 !Jaap: max(hloc,10.d0*z0) 

   ! transport parameters
   Asb=0.005d0*hloc*(D50(jg)/hloc/(delta*par%g*D50(jg)))**1.2d0     ! bed load coefficent
   Ass=0.012d0*D50(jg)*dster(jg)**(-0.6d0)/(delta*par%g*D50(jg))**1.2d0 ! suspended load coeffient

  ! Diane Foster and Robert: limit Shields to par%smax -> vmag2 for transp. limited
  ! vmag2=min(vmag2,par%smax*par%C**2*D50(jg)*delta)
  ! vmag2=min(vmag2,par%smax*par%g/par%cf*D50(jg)*delta)            ! In terms of cf
  ! term1=sqrt(vmag2+0.018d0/Cd*urms2)     ! nearbed-velocity
  ! the two above lines are comment out and replaced by a limit on total velocity u2+urms2, robert 1/9 and ap 28/11

   term1=(vmg**2+0.018/Cd*par%sws*urms2) ! Make 0.018/Cd is always smaller than the flow friction coefficient 

   term1=min(term1,par%smax*par%g/par%cf*s%D50(jg)*delta)
   term1=sqrt(term1)      
   
   term2 = 0.d0
   do j=1,ny+1
      do i=1,nx
	     if(term1(i,j)>Ucr(i,j) .and. hloc(i,j)>par%eps) then
		   term2(i,j)=(term1(i,j)-Ucr(i,j))**2.4d0
           ! term2(i,j)=(term1(i,j)-Ucr(i,j))**(1.7-0.7*tanh((ue(i,j)/sqrt(par%g*hh(i,j))-0.5)*10))
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
   ceqb = Asb*term2                    		      
   ceqb = min(ceqb/hloc,0.05d0)             ! maximum equilibrium bed concentration
   ceqbg(:,:,jg) = ceqb*sedcal(jg)*wetz
   ceqs = Ass*term2                    
   ceqs = min(ceqs/hloc,0.05d0)             ! maximum equilibrium suspended concentration		      
   ceqsg(:,:,jg) = ceqs*sedcal(jg)*wetz

enddo  ! end og grain size classes

! Robert + Pieter : should be faster
if (abs(par%facua)>tiny(0.d0)) then     ! Robert: Very slow loop, so only do if necessary 
   Ur = 3.d0/8.d0*sqrt(2.d0)*H*k/(k*hloc)**3                  !Ursell number
   Ur = max(Ur,0.000000000001d0)
   Bm = m1 + (m2-m1)/(1.d0+beta*Ur**alpha)                    !Boltzmann sigmoid (eq 6)         
   B1 = (-90.d0+90.d0*tanh(m5/Ur**m6))*par%px/180.d0
   Sk = Bm*cos(B1)                                            !Skewness (eq 8)
   As = Bm*sin(B1)                                            !Asymmetry(eq 9)
   ua = par%facua*(Sk-As)*urms
endif

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
integer                                 :: i,ii
integer , save                          :: nh,nt    
integer                                 :: ih0,it0,ih1,it1
integer                                 :: j,jg
real*8                                  :: onethird,twothird,Ass,dcf,dcfin,ML
real*8                                  :: Te,kvis,Sster,cc1,cc2,wster
real*8                                  :: p,q,f0,f1,f2,f3,uad,duddtmax,dudtmax,siguref,t0fac,duddtmean,dudtmean
real*8 , save                           :: dh,dt,delta
real*8 , save                           :: m1,m2,m3,m4,m5,m6
real*8 , save                           :: alpha,beta,Ur,Bm,B1
real*8 , dimension(:),allocatable    ,save     :: w,dster  
real*8 , dimension(:,:),allocatable  ,save     :: vmg,Asb,Ts
real*8 , dimension(:,:),allocatable  ,save     :: urmsturb,Ucr,Ucrc,Ucrw,term1,B2,Cd
real*8 , dimension(:,:),allocatable  ,save     :: hloc,ceqs,ceqb,h0,t0,detadxmax,detadxmean
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
   allocate (urmsturb  (nx+1,ny+1))
   allocate (hloc  (nx+1,ny+1))
   allocate (kb    (nx+1,ny+1))
   allocate (Ts    (nx+1,ny+1))
   allocate (ceqs  (nx+1,ny+1))
   allocate (ceqb  (nx+1,ny+1))
   allocate (w     (par%ngd))
   allocate (dster (par%ngd))
   allocate (RF    (18,33,40))
   vmg = 0.d0
   onethird=1.d0/3.d0
   if (par%waveform==1)then
      m1 = 0;       ! a = 0
      m2 = 0.7939;  ! b = 0.79 +/- 0.023
      m3 = -0.6065; ! c = -0.61 +/- 0.041
      m4 = 0.3539;  ! d = -0.35 +/- 0.032 
      m5 = 0.6373;  ! e = 0.64 +/- 0.025
      m6 = 0.5995;  ! f = 0.60 +/- 0.043
      alpha = -log10(exp(1.d0))/m4
      beta  = exp(m3/m4)
   endif
   ! Robert: do only once, not necessary every time step
   do jg=1,par%ngd
      ! cjaap: compute fall velocity with simple expression from Ahrens (2000)
      Te    = 20.d0
      kvis  = 4.d0/(20.d0+Te)*1d-5	! Van rijn, 1993 
      Sster = s%D50(jg)/(4*kvis)*sqrt((par%rhos/par%rho-1)*par%g*s%D50(jg))
      cc1   = 1.06d0*tanh(0.064d0*Sster*exp(-7.5d0/Sster**2))
	  cc2    = 0.22d0*tanh(2.34d0*Sster**(-1.18d0)*exp(-0.0064d0*Sster**2))
      wster = cc1+cc2*Sster
      w(jg) = wster*sqrt((par%rhos/par%rho-1.d0)*par%g*s%D50(jg))
	  ! RJ: for modeling gravel
	  delta = (par%rhos-par%rho)/par%rho
      dster(jg)=(delta*par%g/1.d-12)**onethird*s%D50(jg)
   enddo
endif

! hloc   = max(hh,0.01d0) ! Jaap 
hloc = hh
onethird=1.d0/3.d0
twothird=2.d0/3.d0
!
! compute wave shape short waves (Rienecker and Fenton + Ruessink and van Rijn to estimate weighting of sine's and cosines)
!
! read table at t=0;
if (abs(par%t-par%dt)<1.d-6) then
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

! read us and duddtmax from table....
h0 = min(nh*dh,max(dh,min(H,hloc)/hloc))  ! Jaap: try this
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
      
      if (par%waveform==1) then
         Ur = 3.d0/8.d0*sqrt(2.d0)*H(i,j)*k(i,j)/(k(i,j)*hloc(i,j))**3   !Ursell number
         Ur = max(Ur,0.000000000001d0)
         Bm = m1 + (m2-m1)/(1.d0+beta*Ur**alpha)                         !Boltzmann sigmoid (eq 6)         
         B1 = (-90.d0+90.d0*tanh(m5/Ur**m6))*par%px/180.d0
         Sk(i,j) = Bm*cos(B1)                                            !Skewness (eq 8)
         As(i,j) = Bm*sin(B1)                                            !Asymmetry(eq 9)
	  elseif (par%waveform==2) then
         Sk(i,j) = f0*RF(13,ih0,it0)+f1*RF(13,ih1,it0)+ f2*RF(13,ih0,it1)+f3*RF(13,ih1,it1)
	     As(i,j) = f0*RF(14,ih0,it0)+f1*RF(14,ih1,it0)+ f2*RF(14,ih0,it1)+f3*RF(14,ih1,it1)
      endif

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
	  
	  ua(i,j) = par%sws*par%facua*(Sk(i,j)-As(i,j))*urms(i,j)

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
! due to short waves
do j=1,ny+1	
	do i=1,nx+1
	   ! compute mixing length
	   ! ML = 2*R(i,j)*par%Trep/(par%rho*c(i,j)*max(H(i,j),par%eps))
	   ML = dsqrt(2*R(i,j)*par%Trep/(par%rho*c(i,j)))
	   ! ML = 0.9d0*H(i,j)
	   ML = min(ML,hloc(i,j));
	   ! exponential decay turbulence over depth
	   dcfin = exp(min(100.d0,hloc(i,j)/max(ML,0.01d0)))
	   dcf = min(1.d0,1.d0/(dcfin-1.d0))
	   if (par%turb == 2) then
          kb(i,j) = par%nuhfac*(DR(i,j)/par%rho)**twothird*dcf*par%Trep/Tbore(i,j)
	   elseif (par%turb == 1) then
	      kb(i,j) = par%nuhfac*(DR(i,j)/par%rho)**twothird*dcf
	   elseif (par%turb == 0) then
	      kb(i,j) = 0.0d0
	   endif
    enddo
enddo

! switch to include long wave stirring
if (par%lws==1) then
   vmg  = dsqrt(ue**2+ve**2)
elseif (par%lws==0) then
   ! vmg lags on actual mean flow; but long wave contribution to mean flow is included... 
   vmg = (1.d0-1.d0/par%cats/par%Trep*par%dt)*vmg + (1.d0/par%cats/par%Trep*par%dt)*dsqrt(ue**2+ve**2) 
endif

urmsturb = dsqrt(urms**2.d0+1.45d0*(kb+kturb))

do jg = 1,par%ngd
   
   Ts       = par%tsfac*hloc/w(jg)
   Tsg(:,:,jg) = max(Ts,par%Tsmin) 
   !
   ! calculate treshold velocity Ucr
   !
   if(s%D50(jg)<=0.0005) then
     Ucrc=0.19d0*D50(jg)**0.1d0*log10(4.d0*hloc/D90(jg))                           !Shields
	 Ucrw=0.24d0*(delta*par%g)**0.66d0*s%D50(jg)**0.33d0*par%Trep**0.33d0          !Komar and Miller (1975)
   else if(s%D50(jg)<0.002) then
     Ucrc=8.5d0*D50(jg)**0.6d0*log10(4.d0*hloc/D90(jg))                            !Shields
	 Ucrw=0.95d0*(delta*par%g)**0.57d0*D50(jg)**0.43*par%Trep**0.14                !Komar and Miller (1975)
   else if(s%D50(jg)>0.002) then
     Ucrc=1.3d0*sqrt(delta*par%g*D50(jg))*(hloc/D50(jg))**(0.5d0*onethird)         !Maynord (1978) --> also Neill (1968) where 1.3d0 = 1.4d0
	 Ucrw=0.95d0*(delta*par%g)**0.57d0*D50(jg)**0.43*par%Trep**0.14                !Komar and Miller (1975)
   end if
   B2 = vmg/max(vmg+urmsturb,par%eps)
   Ucr = B2*Ucrc + (1-B2)*Ucrw                                                     !Van Rijn 2007 (Bed load transport paper)

   ! transport parameters
   Asb=0.015d0*hloc*(s%D50(jg)/hloc)**1.2d0/(delta*par%g*s%D50(jg))**0.75d0        !bed load coefficent
   Ass=0.012d0*s%D50(jg)*dster(jg)**(-0.6d0)/(delta*par%g*s%D50(jg))**1.2d0        !suspended load coeffient
   
   ! Jaap: Gravel test:

   ! Jaap: par%sws to set short wave stirring to zero
   ! Jaap: Van Rijn use Peak orbital flow velocity --> 0.64 corresponds to 0.4 coefficient regular waves Van Rijn (2007)  
   term1= dsqrt(vmg**2+0.64d0*par%sws*urmsturb**2)                                       

   ! Try Soulsby van rijn approach...
   ! drag coefficient
   ! z0 = par%z0
   ! Cd=(0.40/(log(hloc/z0)-1.0))**2 !Jaap
   ! term1=(vmg**2+0.018/Cd*uorb**2)**0.5   
      
   ceqb = 0.d0*term1                                                                     !initialize ceqb
   ceqs = 0.d0*term1                                                                     !initialize ceqs
   do j=1,ny+1
	  do i=1,nx
         if(term1(i,j)>Ucr(i,j) .and. hloc(i,j)>par%eps) then
            ceqb(i,j)=Asb(i,j)*(term1(i,j)-Ucr(i,j))**1.5 
			ceqs(i,j)=Ass*(term1(i,j)-Ucr(i,j))**2.4
		 end if
      end do
   end do
  
   ceqb = min(ceqb/hloc,0.05)		      ! maximum equilibrium bed concentration
   ceqbg(:,:,jg) = ceqb*sedcal(jg)*wetz
   ceqs = min(ceqs/hloc,0.05)		      ! maximum equilibrium suspended concentration
   ceqsg(:,:,jg) = ceqs*sedcal(jg)*wetz
enddo                                 ! end of grain size classes
                            ! end of grain size classes

end subroutine sednew

subroutine longwaveturb(s,par)
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
integer                                  :: j
real                                     :: ksource
real*8,dimension(:,:),allocatable,save   :: kturbu,kturbv,Sturbu,Sturbv,dzsdt_cr

include 's.ind'
include 's.inp'

if (.not. allocated(kturbu)) then
   allocate(kturbu (nx+1,ny+1))
   allocate(kturbv (nx+1,ny+1))
   allocate(Sturbu (nx+1,ny+1))
   allocate(Sturbv (nx+1,ny+1))
   allocate(dzsdt_cr (nx+1,ny+1))
endif
! use lagrangian velocities
kturbu       = 0.0d0  !Jaap
kturbv       = 0.0d0  !Jaap		
dzsdt_cr=par%beta*c
! Update roller thickness
rolthick=rolthick+par%dt*(abs(dzsdt)-dzsdt_cr)
rolthick=max(rolthick,0.d0)
!
!  X-direction
do j=1,ny+1
   do i=1,nx
      if(uu(i,j)>0.d0) then
         kturbu(i,j)=par%thetanum*kturb(i,j)+(1.d0-par%thetanum)*kturb(min(i+1,nx),j)
      elseif(uu(i,j)<0.d0) then
         kturbu(i,j)=par%thetanum*kturb(i+1,j)+(1.d0-par%thetanum)*kturb(max(i,2),j)
      else
         kturbu(i,j)=0.5d0*(kturb(i,j)+kturb(i+1,j))
      endif
   enddo 
enddo
kturbu(nx+1,:) = kturb(nx+1,:) !Robert
! wwvv fix this in parallel case
#ifdef USEMPI
   call xmpi_shift(kturbu,'m:')
#endif
!
Sturbu=kturbu*uu*hu*wetu   !
!
! Y-direction
!
do j=1,ny
   do i=1,nx+1
      if(vv(i,j)>0) then
         kturbv(i,j)=par%thetanum*kturb(i,j)+(1.d0-par%thetanum)*kturb(i,min(j+1,ny))
      else if(vev(i,j)<0) then
         kturbv(i,j)=par%thetanum*kturb(i,j+1)+(1.d0-par%thetanum)*kturb(i,max(j,2))
      else
         kturbv(i,j)=0.5d0*(kturbv(i,j)+kturbv(i,j+1))
      end if
   end do 
end do
kturbv(:,ny+1) = kturb(:,ny+1) !Robert
!
Sturbv=kturbv*vv*hv*wetv
do j=2,ny+1
   do i=2,nx+1
      ksource=par%g*rolthick(i,j)*par%beta*c(i,j)      !only important in shallow water, where c=sqrt(gh)  
      kturb(i,j) = hold(i,j)*kturb(i,j)-par%dt*((Sturbu(i,j)-Sturbu(i-1,j))/(xu(i)-xu(i-1))+&
                                                (Sturbv(i,j)-Sturbv(i,j-1))/(yv(j)-yv(j-1))-&
								                (ksource-par%betad*kturb(i,j)**1.5d0))
      kturb(i,j)=max(kturb(i,j),0.0d0)
   enddo
enddo
! Jaap
kturb = kturb/max(hh,0.01d0)
kturb(1,:)=kturb(2,:)
kturb(:,1)=kturb(:,2)
kturb(nx+1,:)=kturb(nx+1-1,:)
kturb(:,ny+1)=kturb(:,ny+1-1)
! wwvv fix the first and last rows and columns of kturb in parallel case
#ifdef USEMPI
    call xmpi_shift(kturb,'1:')
    call xmpi_shift(kturb,'m:')
    call xmpi_shift(kturb,':1')   ! Maybe not necessary as zb calculated from 1:ny+1, but just in case...
    call xmpi_shift(kturb,':n')   ! Dito
#endif
kturb=kturb*wetz
!

end subroutine longwaveturb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module morphevolution
