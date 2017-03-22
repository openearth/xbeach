module morphevolution
   implicit none
   save

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
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
      use paramsconst
      ! use vsmumod

      implicit none

      type(spacepars),target                   :: s
      type(parameters)                         :: par

      integer                                  :: i,isig
      integer                                  :: j,jg
      real*8                                   :: exp_ero

      real*8,dimension(:),allocatable,save     :: chain,cumchain
      real*8,dimension(:,:),allocatable,save   :: vmag2,uau,uav,um,vm,ueu_sed,uev_sed,veu_sed,vev_sed
      real*8,dimension(:,:),allocatable,save   :: ccvt,dcdz,dsigt,aref
      real*8,dimension(:,:),allocatable,save   :: cc,ccb,cu,cv,Sus,Svs
      real*8,dimension(:,:),allocatable,save   :: cub,cvb,Sub,Svb,pbbedu,pbbedv
      real*8,dimension(:,:),allocatable,save   :: suq3d,svq3d,eswmax,eswbed,sigs,deltas
      real*8,dimension(:,:,:),allocatable,save :: dsig,ccv,sdif,cuq3d,cvq3d,fac

      real*8,dimension(:,:),allocatable,save   :: sinthm,costhm

      real*8                                   :: delta,delta_x,shields,ftheta,psi_x,Sbtot ! Lodewijk: for direction of sediment transport (bed slope effect)
      real*8                                   :: Ssmtot, dzbds,  Sbmtot ! Lodewijk: for magnitude of sediment transport (bed slope effect)

      !include 's.ind'
      !include 's.inp'

      if (.not. allocated(vmag2)) then
         allocate(vmag2 (s%nx+1,s%ny+1))
         allocate(uau (s%nx+1,s%ny+1))
         allocate(uav (s%nx+1,s%ny+1))
         allocate(ueu_sed (s%nx+1,s%ny+1))
         allocate(uev_sed (s%nx+1,s%ny+1))
         allocate(veu_sed (s%nx+1,s%ny+1))
         allocate(vev_sed (s%nx+1,s%ny+1))
         allocate(cu  (s%nx+1,s%ny+1))
         allocate(cv  (s%nx+1,s%ny+1))
         allocate(cc  (s%nx+1,s%ny+1))
         allocate(ccb (s%nx+1,s%ny+1))
         allocate(fac (s%nx+1,s%ny+1,par%ngd))
         allocate(Sus (s%nx+1,s%ny+1))
         allocate(Svs (s%nx+1,s%ny+1))
         allocate(cub (s%nx+1,s%ny+1))
         allocate(cvb (s%nx+1,s%ny+1))
         allocate(Sub (s%nx+1,s%ny+1))
         allocate(Svb (s%nx+1,s%ny+1))
         allocate(pbbedu (s%nx+1,s%ny+1))
         allocate(pbbedv (s%nx+1,s%ny+1))
         allocate(ccvt (s%nx+1,s%ny+1))
         allocate(dcdz (s%nx+1,s%ny+1))
         allocate(dsigt (s%nx+1,s%ny+1))
         allocate(dsig(s%nx+1,s%ny+1,par%kmax))
         allocate(ccv(s%nx+1,s%ny+1,par%kmax))
         allocate(sdif(s%nx+1,s%ny+1,par%kmax))
         allocate(um (s%nx+1,s%ny+1))
         allocate(vm (s%nx+1,s%ny+1))
         allocate(deltas(s%nx+1,s%ny+1))
         allocate(sigs(s%nx+1,s%ny+1))
         allocate(eswmax(s%nx+1,s%ny+1))
         allocate(eswbed(s%nx+1,s%ny+1))
         allocate(suq3d(s%nx+1,s%ny+1))
         allocate(svq3d(s%nx+1,s%ny+1))
         allocate(cuq3d(s%nx+1,s%ny+1,par%kmax))
         allocate(cvq3d(s%nx+1,s%ny+1,par%kmax))
         allocate(aref(s%nx+1,s%ny+1))
         allocate(chain(par%kmax))
         allocate(cumchain(par%kmax))
         allocate (sinthm(s%nx+1,s%ny+1))
         allocate (costhm(s%nx+1,s%ny+1))
         delta_x   = 0.d0 ! Lodewijk
         shields   = 0.d0 ! Lodewijk
         ftheta    = 0.d0 ! Lodewijk
         psi_x     = 0.d0 ! Lodewijk
         Sbtot     = 0.d0 ! Lodewijk
         delta     = 0.d0 ! Lodewijk
         uau       = 0.d0
         uav       = 0.d0
         um        = 0.d0
         vm        = 0.d0
         chain     = 0.0d0
         cumchain  = 0.0d0
         fac       = 1.d0
         exp_ero   = 0.d0
         ! generate sigma grid shape...
         do isig=2,par%kmax
            chain(isig) = par%sigfac**(isig-1)
            cumchain(isig) = cumchain(isig-1)+chain(isig)
         enddo
      endif

      ! use eulerian velocities
      vmag2     = s%ue**2+s%ve**2
      cu        = 0.0d0
      cv        = 0.0d0
      cub       = 0.0d0
      cvb       = 0.0d0
      s%dcsdx     = 0.0d0
      s%dcsdy     = 0.0d0

      sinthm = sin(s%thetamean-s%alfaz)
      costhm = cos(s%thetamean-s%alfaz)

      ! short wave runup
      if (par%swrunup==1 .and. par%struct==1) then
         call hybrid(s,par)
      endif

      ! compute turbulence due to wave breaking
      if (par%lwt==1 .or. par%turb == TURB_BORE_AVERAGED .or. par%turb == TURB_WAVE_AVERAGED) then
         call waveturb(s,par)
      endif

      if (par%swave==1) then
         ! include wave skewness and assymetry in sediment advection velocity
         if (par%waveform==WAVEFORM_RUESSINK_VANRIJN)then
            call RvR(s,par)
         elseif (par%waveform==WAVEFORM_VANTHIEL) then
            call vT(s,par)
         endif

      endif

      ! calculate equilibrium concentration/sediment source
      if ((par%form==FORM_SOULSBY_VANRIJN) .or. (par%form==FORM_VANTHIEL_VANRIJN))then           ! Soulsby van Rijn
         ! Soulsby van Rijn and Van Thiel de Vries & Reniers 2008 formulations
         call sedtransform(s,par)
      end if

      ! compute reduction factor for sediment sources due to presence of hard layers
      do jg = 1,par%ngd
         do j=1,s%ny+1
            do i=1,s%nx+1
               exp_ero = par%morfac*par%dt/(1.d0-par%por)*s%hh(i,j)*(s%ceqsg(i,j,jg)*s%pbbed(i,j,1,jg)/s%Tsg(i,j,jg) &
               + s%ceqbg(i,j,jg)*s%pbbed(i,j,1,jg)/par%dt)
               ! limit erosion to available sediment on top of hard layer wwvv changed tiny into epsilon
               fac(i,j,jg) = min(1.d0,s%structdepth(i,j)*s%pbbed(i,j,1,jg)/max(epsilon(0.d0),exp_ero) )
               !if (fac(i,j,jg)*exp_ero > dzbed(i,j,1)*pbbed(i,j,1,jg)) then
               !   limit erosion to available sand in top layer
               !   fac(i,j,jg) = min(fac(i,j,jg),dzbed(i,j,1)*pbbed(i,j,1,jg)/max(tiny(0.d0),exp_ero) )
               !   write(*,*)'WARNING: expected erosion from top layer is larger than available sand in top layer'
               !endif
            enddo
         enddo
      enddo

      ! compute diffusion coefficient
      s%Dc = par%facDc*(par%nuh+par%nuhfac*s%hh*(s%DR/par%rho)**(1.d0/3.d0))

      do jg = 1,par%ngd
         cc = s%ccg(:,:,jg)
         if (s%D50(jg)>0.002d0) then
            ! RJ: set ceqsg to zero for gravel.
            ! Dano: try without this fix cc = 0.d0 ! Can be used to test total transport mode
         endif
         !
         ! X-direction
         !
         ! Get ua in u points and split out in u and v direction
         uau(1:s%nx,:) = 0.5*(s%ua(1:s%nx,:)*costhm(1:s%nx,:)+s%ua(2:s%nx+1,:)*costhm(1:s%nx,:))
         uav(1:s%nx,:) = 0.5*(s%ua(1:s%nx,:)*sinthm(1:s%nx,:)+s%ua(2:s%nx+1,:)*sinthm(1:s%nx,:))
         if (par%nz>1) then
            ueu_sed(1:s%nx,:) = 0.5*(s%ue_sed(1:s%nx,:)+s%ue_sed(2:s%nx+1,:))
         else
            ueu_sed=s%ueu
         endif
         veu_sed(1:s%nx,:) = 0.5*(s%ve_sed(1:s%nx,:)+s%ve_sed(2:s%nx+1,:))
         if (xmpi_isbot) then
            veu_sed(s%nx+1,:) = veu_sed(s%nx,:)
         endif

         ! Compute vmagu including ua
!         s%vmagu = sqrt((s%uu+uau)**2+(s%vu+uav)**2)
         s%vmagu = sqrt((ueu_sed+uau)**2+(veu_sed+uav)**2)
         ! sediment advection velocity for suspended load and bed load respectively
         ! REMARK: when vreps does not equal vv; no mass conservation
!         s%ureps = s%ueu+uau
!         s%urepb = s%ueu+uau  ! RJ maybe reduce this velocity?
         s%ureps = ueu_sed+uau
         s%urepb = ueu_sed+uau  
         !
         do j=1,s%ny+1
            do i=1,s%nx
               if(s%ureps(i,j)>0.d0) then
                  ! test cu(i,j)=cc(i,j)
                  cu(i,j)=par%thetanum*cc(i,j)+(1.d0-par%thetanum)*cc(min(i+1,s%nx),j)
                  cub(i,j)=par%thetanum*s%pbbed(i,j,1,jg)*s%ceqbg(i,j,jg)+(1.d0-par%thetanum)&
                  *s%pbbed(min(i+1,s%nx),j,1,jg)*s%ceqbg(min(i+1,s%nx),j,jg)
                  !cub(i,j)=par%thetanum*ccb(i,j)+(1.d0-par%thetanum)*ccb(min(i+1,nx),j)
               elseif(s%ureps(i,j)<0.d0) then
                  cu(i,j)=par%thetanum*cc(i+1,j)+(1.d0-par%thetanum)*cc(max(i,2),j)
                  cub(i,j)=par%thetanum*s%pbbed(i+1,j,1,jg)*s%ceqbg(i+1,j,jg)+(1.d0-par%thetanum)&
                  *s%pbbed(max(i,2),j,1,jg)*s%ceqbg(max(i,2),j,jg)
                  !cub(i,j)=par%thetanum*ccb(i+1,j)+(1.d0-par%thetanum)*ccb(max(i,2),j)
               else
                  cu(i,j)=0.5d0*(cc(i,j)+cc(i+1,j))
                  cub(i,j)=0.5d0*(s%pbbed(i,j,1,jg)*s%ceqbg(i,j,jg)+s%pbbed(i+1,j,1,jg)*s%ceqbg(i+1,j,jg))
                  !cub(i,j)=0.5d0*(ccb(i,j)+ccb(i+1,j))
               endif
               s%dcsdx(i,j)=(cc(i+1,j)-cc(i,j))/s%dsu(i,j)
            enddo
         enddo
         ! wwvv dcdx(nx:1,:) is still untouched, correct this ofr the parallel case
         cu(s%nx+1,:) = cc(s%nx+1,:) !Robert
         ! wwvv fix this in parallel case
         ! wwvv in parallel version, there will be a discrepancy between the values
         ! of dzbdx(nx+1,:).
         !wwvv so fix that
         !
         Sus = 0.d0
         Sub = 0.d0
         !
         ! suspended load, Lodewijk: no bed slope effect (yet)
         Sus=par%sus*(cu*s%ureps*s%hu-s%Dc*s%hu*s%dcsdx)*s%wetu
         ! bed load, Lodewijk: no bed slope effect (yet)
         Sub=par%bed*(cub*s%urepb*s%hu)*s%wetu
         !
         ! Originally bed slope effect of XBeach : par%bdslpeffmag = 1
         ! Original one, but only on bed load : par%bdslpeffmag = 2
         if (par%bdslpeffmag == BDSLPEFFMAG_ROELV_TOTAL) then
            if (par%bermslope>0) then
               where (s%H/s%hu>1.0d0.or.(s%hu<1.d0.and. par%instat==INSTAT_STAT_TABLE))
                  Sus = Sus-par%sus*(10.d0*par%facsl*cu*s%vmagu*s%hu*(s%dzbdx-par%bermslope))*s%wetu
               elsewhere
                  Sus = Sus-par%sus*(par%facsl*cu*s%vmagu*s%hu*s%dzbdx)*s%wetu
               end where
            else
               Sus = Sus-par%sus*(par%facsl*cu*s%vmagu*s%hu*s%dzbdx)*s%wetu
            endif
         endif

         if (par%bdslpeffmag == BDSLPEFFMAG_ROELV_TOTAL .or. par%bdslpeffmag == BDSLPEFFMAG_ROELV_BED) then
            if (par%bermslope>0) then
               where (s%H/s%hu>1.0d0 .or. (s%hu<1.d0.and. par%instat==INSTAT_STAT_TABLE))
                  Sub = Sub-par%bed*(10.0d0*par%facsl*cub*s%vmagu*s%hu*(s%dzbdx-par%bermslope))*s%wetu
               elsewhere
                  Sub = Sub-par%bed*(par%facsl*cub*s%vmagu*s%hu*s%dzbdx)*s%wetu
               end where
            else
               Sub = Sub-par%bed*(par%facsl*cub*s%vmagu*s%hu*s%dzbdx)*s%wetu 
            endif
         endif
         !
         !
         ! Y-direction
         !
         ! Jaap: get ua in v points and split out in u and v direction
         if (s%ny>0) then
            uau(:,1:s%ny) = 0.5*(s%ua(:,1:s%ny)*costhm(:,1:s%ny)+s%ua(:,2:s%ny+1)*costhm(:,1:s%ny))
            uav(:,1:s%ny) = 0.5*(s%ua(:,1:s%ny)*sinthm(:,1:s%ny)+s%ua(:,2:s%ny+1)*sinthm(:,1:s%ny))
            uau(:,s%ny+1) = uau(:,s%ny) ! Jaap
            uav(:,s%ny+1) = uav(:,s%ny) ! Jaap
            if (par%nz>1) then
               vev_sed(:,1:s%ny) = 0.5*(s%ve_sed(:,1:s%ny)+s%ve_sed(:,2:s%ny+1))
            else
               vev_sed=s%vev
            endif
            uev_sed(:,1:s%ny) = 0.5*(s%ue_sed(:,1:s%ny)+s%ue_sed(:,2:s%ny+1))
            if (xmpi_isright) then
               uev_sed(:,1:s%ny+1) = uev_sed(:,1:s%ny)
            endif
        else
            uau=s%ua*costhm
            uav=s%ua*sinthm
            uev_sed=s%ue_sed
            vev_sed=s%ve_sed
         endif
         ! Jaap: compute vmagv including ua
!         s%vmagv = sqrt((s%uv+uau)**2+(s%vv+uav)**2)
         s%vmagv = sqrt((uev_sed+uau)**2+(vev_sed+uav)**2)
         ! sediment advection velocity for suspended load and bed load respectively
         ! REMARK: when vreps does not equal vv; no mass conservation
!         s%vreps = s%vev+uav
!         s%vrepb = s%vev+uav   ! RJ maybe reduce this velocity? Should be s%vv instead of s%vev?
         s%vreps = vev_sed+uav
         s%vrepb = vev_sed+uav   ! RJ maybe reduce this velocity? Should be s%vv instead of s%vev?
         !
         if (s%ny>0) then
            do j=1,s%ny
               do i=1,s%nx+1
                  if(s%vreps(i,j)>0) then
                     cv(i,j)=par%thetanum*cc(i,j)+(1.d0-par%thetanum)*cc(i,min(j+1,s%ny))
                     cvb(i,j)=par%thetanum*s%pbbed(i,j,1,jg)*s%ceqbg(i,j,jg)+(1.d0-par%thetanum)&
                     *s%pbbed(i,min(j+1,s%ny),1,jg)*s%ceqbg(i,min(j+1,s%ny),jg)
                     !cvb(i,j)=par%thetanum*ccb(i,j)+(1.d0-par%thetanum)*ccb(i,min(j+1,ny))
                  elseif(s%vreps(i,j)<0) then
                     cv(i,j)=par%thetanum*cc(i,j+1)+(1.d0-par%thetanum)*cc(i,max(j,2))
                     cvb(i,j)=par%thetanum*s%pbbed(i,j+1,1,jg)*s%ceqbg(i,j+1,jg)+(1.d0-par%thetanum)&
                     *s%pbbed(i,max(j,2),1,jg)*s%ceqbg(i,max(j,2),jg)
                     !cvb(i,j)=par%thetanum*ccb(i,j+1)+(1.d0-par%thetanum)*ccb(i,max(j,2))
                  else
                     cv(i,j)=0.5d0*(cc(i,j)+cc(i,j+1)) !Jaap: cc instead of cv
                     cvb(i,j)=0.5d0*(s%pbbed(i,j,1,jg)*s%ceqbg(i,j,jg)+s%pbbed(i,j+1,1,jg)*s%ceqbg(i,j+1,jg))
                     !cvb(i,j)=0.5d0*(ccb(i,j)+ccb(i,j+1))
                  end if
                  s%dcsdy(i,j)=(cc(i,j+1)-cc(i,j))/s%dnv(i,j) !Jaap

               end do
            end do
            ! wwvv dcdy(:,ny+1) is not filled in, so in parallel case:
            cv(:,s%ny+1) = cc(:,s%ny+1) !Robert
            ! wwvv in parallel version, there will be a discrepancy between the values
            ! of dzbdy(:,ny+1).
            ! wwvv so fix that
         else
            cv = cc
            cvb = s%ceqbg(:,:,jg)
         endif ! s%ny>0
         !
         ! Compute sedimnent transport in v-direction
         !
         Svs = 0.d0
         Svb = 0.d0
         ! Suspended load
         Svs=par%sus*(cv*s%vreps*s%hv-s%Dc*s%hv*s%dcsdy)*s%wetv
         ! Bed load
         Svb=par%bed*(cvb*s%vrepb*s%hv)*s%wetv
         !
         ! Originally bed slope magnitude effect of XBeach : par%bdslpeffmag = 1
         ! Original one, but only on bed load : par%bdslpeffmag = 2
         if (par%bdslpeffmag == BDSLPEFFMAG_ROELV_TOTAL) then
            Svs = Svs-par%sus*(par%facsl*cv*s%vmagv*s%hv*s%dzbdy)*s%wetv
         endif
         if (par%bdslpeffmag == BDSLPEFFMAG_ROELV_TOTAL .or. par%bdslpeffmag == BDSLPEFFMAG_ROELV_BED  ) then
            Svb = Svb-par%bed*(par%facsl*cvb*s%vmagv*s%hv*s%dzbdy)*s%wetv
         endif
         !
         !
         ! Bed slope magnitude effect (as Souslby intended) and change direction transport (see Van Rijn 1993 (section 7.2.6))
         !
         !
         if (par%bdslpeffmag == BDSLPEFFMAG_SOULS_TOTAL .or. par%bdslpeffmag == BDSLPEFFMAG_SOULS_BED) then
            do j=1,s%ny+1
               do i=1,s%nx+1
                  if ((dabs(Sub(i,j)) > 0.000001d0) .or. (dabs(Svb(i,j)) > 0.000001d0)) then
                     Sbmtot = dsqrt(  Sub(i,j)**2.d0  +  Svb(i,j)**2.d0   )
                     dzbds = s%dzbdx(i,j)*Sub(i,j)/Sbmtot + s%dzbdy(i,j)*Svb(i,j)/Sbmtot
                     ! dzbdn = s%dzbdx*Svb/Sbtot + s%dzbdy*Sub/Sbtot
                     Sub(i,j) = Sub(i,j)*(1.d0 - par%facsl*dzbds)
                     Svb(i,j) = Svb(i,j)*(1.d0 - par%facsl*dzbds)
                     !
                  endif
                  if (((dabs(Sus(i,j)) > 0.000001d0) .or. (dabs(Svs(i,j)) > 0.000001d0)) &
                  .and. par%bdslpeffmag == BDSLPEFFMAG_SOULS_TOTAL) then
                     Ssmtot = dsqrt(  Sus(i,j)**2.d0  +  Svs(i,j)**2.d0  )
                     dzbds = s%dzbdx(i,j)*Sus(i,j)/Ssmtot + s%dzbdy(i,j)*Svs(i,j)/Ssmtot
                     ! dzbdn = s%dzbdx*Svb/Sbtot + s%dzbdy*Sub/Sbtot
                     Sus(i,j) = Sus(i,j)*(1.d0 - par%facsl*dzbds)
                     Svs(i,j) = Svs(i,j)*(1.d0 - par%facsl*dzbds)
                  endif
               enddo
            enddo
         endif
         !
         ! Lodewijk: modify the direction of the bed load transport based on the bed slope, see Van Rijn 1993 (section 7.2.6)
         if (par%bdslpeffdir == BDSLPEFFDIR_TALMON) then
            do j=1,s%ny+1
               do i=1,s%nx+1
                  if (((dabs(s%urepb(i,j)) > 0.0001d0) .or. (dabs(s%vrepb(i,j)) > 0.0001d0)) &
                  .and. ((dabs(s%taubx(i,j)) > 0.0001d0) .or. (dabs(s%tauby(i,j)) > 0.0001d0))) then
                     if (s%urepb(i,j) < 0.d0) then
                        delta_x = datan(s%vrepb(i,j)/s%urepb(i,j))+par%px ! Angle between fluid velocity vector and the s%x-axis
                     else
                        delta_x = datan(s%vrepb(i,j)/s%urepb(i,j))        ! Angle between fluid velocity vector and the s%x-axis
                     endif
                     delta = (par%rhos-par%rho)/par%rho
                     shields = sqrt(s%taubx(i,j)**2 + s%tauby(i,j)**2)/(delta*par%rho*par%g*s%D50(jg))
                     ! shields = (urepb(i,j)**2.d0+vrepb(i,j)**2.d0)*s%cf(i,j)/(par%g*D50(jg)*delta)
                     ftheta = 1.d0/(9.d0*(s%D50(jg)/s%hh(i,j))**0.3d0*sqrt(shields)) ! Talmon
                     psi_x = datan(  (dsin(delta_x)-ftheta*s%dzbdy(i,j))  /  (dcos(delta_x)-ftheta*s%dzbdx(i,j))  )
                     psi_x = par%bdslpeffdirfac*(psi_x - delta_x)+delta_x
                     Sbtot = dsqrt(  Sub(i,j)**2.d0  +  Svb(i,j)**2.d0  )  ! Magnitude of sediment transport without direction modifcation
                     ! Decompose the sediment transport again, know with the knowledge of the direction of the sediment transport vector
                     Sub(i,j) = Sbtot * dcos(psi_x)
                     Svb(i,j) = Sbtot * dsin(psi_x)
                  else
                     Sub(i,j) = 0.d0
                     Svb(i,j) = 0.d0
                  endif
               enddo
            enddo
         endif
         !
         !
         !
         do j=1,s%ny+1
            do i=1,s%nx
               if(Sub(i,j)>0.d0) then
                  pbbedu(i,j) = s%pbbed(i,j,1,jg)
               elseif(Sub(i,j)<0.d0) then
                  pbbedu(i,j)= s%pbbed(i+1,j,1,jg)
               else
                  pbbedu(i,j)=0.5d0*(s%pbbed(i,j,1,jg)+s%pbbed(i+1,j,1,jg))
               endif
            enddo
         enddo
         !
         Sub = pbbedu*Sub
         !
         do j=1,s%ny
            do i=1,s%nx+1
               if(Svb(i,j)>0) then
                  pbbedv(i,j)=s%pbbed(i,j,1,jg)
               else if(Svb(i,j)<0) then
                  pbbedv(i,j)=s%pbbed(i,j+1,1,jg)
               else
                  pbbedv(i,j)=0.5d0*(s%pbbed(i,j,1,jg)+s%pbbed(i,j+1,1,jg))
               end if
            end do
         end do
         !
         Svb = pbbedv*Svb
         !
         ! BRJ: implicit concentration update (compute sources first, sink must be computed after updating actual sed.conc.)
         !
         if (s%ny>0) then
            do j=2,s%ny
               do i=2,s%nx
                  ! Changed to hh from hold by RJ (13072009) !**2/max(hh(i,j),par%hmin)
                  s%ero(i,j,jg) = fac(i,j,jg)*s%hh(i,j)*s%ceqsg(i,j,jg)*s%pbbed(i,j,1,jg)/s%Tsg(i,j,jg)
                  ! depo_ex(i,j,jg) = max(hold(i,j),0.01d0)*cc(i,j)/Tsg(i,j,jg)
                  ! BRJ: the volume in the water column is updated and not the volume concentration.
                  cc(i,j) = (par%dt*s%Tsg(i,j,jg))/(par%dt+s%Tsg(i,j,jg))* &
                  (s%hold(i,j)*cc(i,j)/par%dt -((Sus(i,j)*s%dnu(i,j)-Sus(i-1,j)*s%dnu(i-1,j)+&
                  Svs(i,j)*s%dsv(i,j)-Svs(i,j-1)*s%dsv(i,j-1))*s%dsdnzi(i,j)-&
                  s%ero(i,j,jg)))

                  cc(i,j)=max(cc(i,j),0.0d0) ! Jaap: negative cc's are possible...
                  cc(i,j)=min(cc(i,j),par%cmax*s%hh(i,j))
                  s%depo_ex(i,j,jg) = cc(i,j)/s%Tsg(i,j,jg)
               enddo
            enddo

         else
            j=1
            do i=2,s%nx
               ! Changed to hh from hold by RJ (13072009) !**2/max(hh(i,j),par%hmin)
               s%ero(i,j,jg) = fac(i,j,jg)*s%hh(i,j)*s%ceqsg(i,j,jg)*s%pbbed(i,j,1,jg)/s%Tsg(i,j,jg)
               ! depo_ex(i,j,jg) = max(hold(i,j),0.01d0)*cc(i,j)/Tsg(i,j,jg)
               ! BRJ: the volume in the water column is updated and not the volume concentration.
               cc(i,j) = (par%dt*s%Tsg(i,j,jg))/(par%dt+s%Tsg(i,j,jg))* &
               (s%hold(i,j)*cc(i,j)/par%dt -((Sus(i,j)*s%dnu(i,j)-Sus(i-1,j)*s%dnu(i-1,j))*s%dsdnzi(i,j)-&
               s%ero(i,j,jg)))

               cc(i,j)=max(cc(i,j),0.0d0) ! Jaap: negative cc's are possible...
               cc(i,j)=min(cc(i,j),par%cmax*s%hh(i,j))

               s%depo_ex(i,j,jg) = cc(i,j)/s%Tsg(i,j,jg)
            enddo
         endif


         cc = cc/s%hh

         ! do lateral boundaries...
         if(xmpi_istop)then
            cc(1,:)=cc(2,:)
            s%ero(1,:,jg)=s%ero(2,:,jg)
            s%depo_ex(1,:,jg)=s%depo_ex(2,:,jg)
         endif
         if(xmpi_isleft .and. s%ny>0)then
            cc(:,1)=cc(:,2)
            s%ero(:,1,jg)=s%ero(:,2,jg)
            s%depo_ex(:,1,jg)=s%depo_ex(:,2,jg)
         endif
         if(xmpi_istop)then
            cc(s%nx+1,:)=cc(s%nx+1-1,:)
            s%ero(s%nx+1,:,jg)=s%ero(s%nx,:,jg)
            s%depo_ex(s%nx+1,:,jg)=s%depo_ex(s%nx,:,jg)
         endif
         if(xmpi_isright .and. s%ny>0)then
            cc(:,s%ny+1)=cc(:,s%ny+1-1)
            s%ero(:,s%ny+1,jg)=s%ero(:,s%ny,jg)
            s%depo_ex(:,s%ny+1,jg)=s%depo_ex(:,s%ny,jg)
         endif

         ! wwvv fix the first and last rows and columns of cc in parallel case
#ifdef USEMPI
         call xmpi_shift_ee(cc)
#endif
         !
         ! wwvv border columns and rows of ccg Svg and Sug have to be communicated
         s%ccg(:,:,jg) = cc
         s%Svsg(:,:,jg) = Svs
         s%Susg(:,:,jg) = Sus
         s%Svbg(:,:,jg) = Svb
         s%Subg(:,:,jg) = Sub

      end do ! number of sediment fractions

      s%vmag=sqrt(max(vmag2,par%umin))

   end subroutine transus

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine bed_update(s,par)
      use params
      use interp
      use spaceparams
      use xmpi_module

      implicit none

      type(spacepars),target              :: s
      type(parameters)                    :: par

      integer                                     :: i,j,j1,jg,ii,ie,id,je,jd,jdz,ndz, hinterland
      integer , dimension(:,:,:),allocatable,save :: indSus,indSub,indSvs,indSvb
      real*8                                      :: dzb,dzmax,dzt,dzleft,sdz,dzavt,fac,Savailable,dAfac
      real*8 , dimension(:,:),allocatable,save    :: Sout,hav
      real*8 , dimension(par%ngd)                 :: edg,edg1,edg2,dzg
      real*8 , dimension(:),pointer               :: dz
      real*8 , dimension(:,:),pointer             :: pb
      logical                                     :: aval

      !include 's.ind'
      !include 's.inp'

      if (.not. allocated(Sout)) then
         allocate(Sout(s%nx+1,s%ny+1))
         allocate(hav(s%nx+1,s%ny+1))
         allocate(indSus(s%nx+1,s%ny+1,par%ngd))
         allocate(indSub(s%nx+1,s%ny+1,par%ngd))
         allocate(indSvs(s%nx+1,s%ny+1,par%ngd))
         allocate(indSvb(s%nx+1,s%ny+1,par%ngd))
      endif

      ! Super fast 1D
      if (s%ny==0) then
         j1 = 1
      else
         j1 = 2
      endif
      s%dzbnow  = 0.d0
      dzb    = 0.d0
      if (par%t>=par%morstart .and. par%t < par%morstop .and. par%morfac > .999d0) then
         !
         ! bed_predict
         !
         ! reduce sediment transports when hard layer comes to surface
         ! this step is mainly necessary at the transition from hard layers to sand
         if (par%struct == 1) then

            do jg = 1,par%ngd
               indSus = 0
               indSub = 0
               indSvs = 0
               indSvb = 0
               Sout   = 0.d0
               do j=j1,s%ny+1
                  do i=2,s%nx+1
                     ! fluxes at i,j
                     if (s%Subg(i,j,jg) > 0.d0) then      ! bed load s%u-direction
                        indSub(i,j,jg) = 1
                        Sout(i,j) = Sout(i,j) + s%Subg(i,j,jg)*s%dnu(i,j)
                     endif
                     ! fluxes at i-1,j
                     if (s%Subg(i-1,j,jg) < 0.d0 ) then   ! bed load s%u-direction
                        Sout(i,j) = Sout(i,j) - s%Subg(i-1,j,jg)*s%dnu(i-1,j)
                     endif
                     if (par%sourcesink==0) then
                        ! fluxes at i,j
                        if (s%Susg(i,j,jg) > 0.d0 ) then     ! suspended load s%u-direction
                           indSus(i,j,jg) = 1
                           Sout(i,j) = Sout(i,j) + s%Susg(i,j,jg)*s%dnu(i,j)
                        endif
                        ! fluxes at i-1,j
                        if (s%Susg(i-1,j,jg) < 0.d0 ) then   ! suspended load s%u-direction
                           Sout(i,j) = Sout(i,j) - s%Susg(i-1,j,jg)*s%dnu(i-1,j)
                        endif
                     endif
                  enddo
               enddo
               if (s%ny>0) then
                  do j=j1,s%ny+1
                     do i=2,s%nx+1
                        if (s%Svbg(i,j,jg) > 0.d0 ) then     ! bed load s%v-direction
                           indSvb(i,j,jg) = 1
                           Sout(i,j) = Sout(i,j) + s%Svbg(i,j,jg)*s%dsv(i,j)
                        endif
                        ! fluxes at i,j-1
                        if (s%Svbg(i,j-1,jg) < 0.d0 ) then   ! bed load s%v-direction
                           Sout(i,j) = Sout(i,j) - s%Svbg(i,j-1,jg)*s%dsv(i,j-1)
                        endif
                        if (par%sourcesink==0) then
                           if (s%Svsg(i,j,jg) > 0.d0 ) then     ! suspended load s%v-direction
                              indSvs(i,j,jg) = 1
                              Sout(i,j) = Sout(i,j) + s%Svsg(i,j,jg)*s%dsv(i,j)
                           endif
                           ! fluxes at i,j-1
                           if (s%Svsg(i,j-1,jg) < 0.d0 ) then   ! suspended load s%v-direction
                              Sout(i,j) = Sout(i,j) - s%Svsg(i,j-1,jg)*s%dsv(i,j-1)
                           endif
                        endif ! sourcesink = 0
                     enddo !s%nx+1
                  enddo !s%ny+1
               endif !s%ny>0
               !
               do j=j1,s%ny+1
                  do i=2,s%nx+1
                     Savailable = s%structdepth(i,j)*s%pbbed(i,j,1,jg)/par%morfac/par%dt*(1.d0-par%por)/s%dsdnzi(i,j)
                     ! reduction factor for cell outgoing sediment transports wwvv changed tiny into epsilon
                     fac  = min(1.d0,Savailable/max(Sout(i,j),epsilon(0.d0)) )
                     ! fix sediment transports for the presence of a hard layer; remind indSus etc are 1 in cases of cell outgoing transports
                     ! updated S         oell outgoing transports                  cell incoming transports
                     if (fac < 1.d0)then
                        s%Subg(i,j,jg)   = fac*indSub(i,j,jg)*s%Subg(i,j,jg)         + (1-indSub(i,j,jg))*s%Subg(i,j,jg)
                        s%Subg(i-1,j,jg) = fac*(1-indSub(i-1,j,jg))*s%Subg(i-1,j,jg) + indSub(i-1,j,jg)*s%Subg(i-1,j,jg)
                        if (s%ny>0) then
                           s%Svbg(i,j,jg)   = fac*indSvb(i,j,jg)*s%Svbg(i,j,jg)         + (1-indSvb(i,j,jg))*s%Svbg(i,j,jg)
                           s%Svbg(i,j-1,jg) = fac*(1-indSvb(i,j-1,jg))*s%Svbg(i,j-1,jg) + indSvb(i,j-1,jg)*s%Svbg(i,j-1,jg)
                        endif
                        if (par%sourcesink==0) then
                           s%Susg(i,j,jg)   = fac*indSus(i,j,jg)*s%Susg(i,j,jg)         + (1-indSus(i,j,jg))*s%Susg(i,j,jg)
                           s%Susg(i-1,j,jg) = fac*(1-indSus(i-1,j,jg))*s%Susg(i-1,j,jg) + indSus(i-1,j,jg)*s%Susg(i-1,j,jg)
                           if (s%ny>0) then
                              s%Svsg(i,j,jg)   = fac*indSvs(i,j,jg)*s%Svsg(i,j,jg)         + (1-indSvs(i,j,jg))*s%Svsg(i,j,jg)
                              s%Svsg(i,j-1,jg) = fac*(1-indSvs(i,j-1,jg))*s%Svsg(i,j-1,jg) + indSvs(i,j-1,jg)*s%Svsg(i,j-1,jg)
                           endif !s%ny = 0
                        endif ! sourcesink = 0
                     endif !fac<1.d0
                  enddo ! s%nx+1
               enddo !s%ny + 1
            enddo !par%ngd
         endif !struct == 1

         if (s%ny>0) then
            do j=2,s%ny
               do i=2,s%nx

                  ! bed level changes per fraction in this morphological time step in meters sand including pores
                  ! positive in case of erosion
                  if (par%sourcesink==0) then
                     dzg=par%morfac*par%dt/(1.d0-par%por)*( & ! Dano, dz from sus transport gradients
                     ( s%Susg(i,j,:)*s%dnu(i,j)-s%Susg(i-1,j,:)*s%dnu(i-1,j) +&
                     s%Svsg(i,j,:)*s%dsv(i,j)-s%Svsg(i,j-1,:)*s%dsv(i,j-1) +&
                     ! dz from bed load transport gradients
                     s%Subg(i,j,:)*s%dnu(i,j)-s%Subg(i-1,j,:)*s%dnu(i-1,j)+&
                     s%Svbg(i,j,:)*s%dsv(i,j)-s%Svbg(i,j-1,:)*s%dsv(i,j-1) )*s%dsdnzi(i,j)    )
                  elseif (par%sourcesink==1) then
                     dzg=par%morfac*par%dt/(1.d0-par%por)*( &
                     s%ero(i,j,:)-s%depo_ex(i,j,:)   +&
                     ( s%Subg(i,j,:)*s%dnu(i,j)-s%Subg(i-1,j,:)*s%dnu(i-1,j)+&
                     s%Svbg(i,j,:)*s%dsv(i,j)-s%Svbg(i,j-1,:)*s%dsv(i,j-1) )*s%dsdnzi(i,j)    )
                  endif

                  if (par%ngd==1) then ! Simple bed update in case one fraction

                     s%zb(i,j) = s%zb(i,j)-sum(dzg)
                     s%dzbnow(i,j) = s%dzbnow(i,j)-sum(dzg) ! naamgeveing?
                     s%dzbdt(i,j) = s%dzbnow(i,j)/par%dt
                     s%sedero(i,j) = s%sedero(i,j)-sum(dzg)
                     s%structdepth(i,j) = max(0.d0,s%structdepth(i,j)-sum(dzg))

                  else
                     ! erosion/deposition rate of sand mass (m/s)
                     ! positive in case of erosion
                     edg = dzg*(1.d0-par%por)/par%dt


                     dz=>s%dzbed(i,j,:)
                     pb=>s%pbbed(i,j,:,:)

                     call update_fractions(par,s,i,j,dz,pb,edg,sum(dzg))

                  endif

               enddo ! s%nx+1
            enddo ! s%ny+1
         else
            j=1
            do i=2,s%nx
               ! bed level changes per fraction in this morphological time step in meters sand including pores
               ! positive in case of erosion
               if (par%sourcesink==0) then
                  dzg=par%morfac*par%dt/(1.d0-par%por)*( & ! Dano, dz from sus transport gradients
                  ( s%Susg(i,j,:)*s%dnu(i,j)-s%Susg(i-1,j,:)*s%dnu(i-1,j) +&
                    s%Subg(i,j,:)*s%dnu(i,j)-s%Subg(i-1,j,:)*s%dnu(i-1,j) )*s%dsdnzi(i,j) +&
                   (s%Svsg(i,j,:)+s%Svbg(i,j,:))*par%lsgrad)
               elseif (par%sourcesink==1) then
                  dzg=par%morfac*par%dt/(1.d0-par%por)*( &
                  s%ero(i,j,:)-s%depo_ex(i,j,:)       +&
                  ( s%Subg(i,j,:)*s%dnu(i,j)-s%Subg(i-1,j,:)*s%dnu(i-1,j) )*s%dsdnzi(i,j) +&
                   (s%Svsg(i,j,:)+s%Svbg(i,j,:))*par%lsgrad)
               endif

               if (par%ngd==1) then ! Simple bed update in case one fraction

                  s%zb(i,j) = s%zb(i,j)-sum(dzg)
                  s%dzbnow(i,j) = s%dzbnow(i,j)-sum(dzg)
                  s%dzbdt(i,j) = s%dzbnow(i,j)/par%dt
                  s%sedero(i,j) = s%sedero(i,j)-sum(dzg)
                  s%structdepth(i,j) = max(0.d0,s%structdepth(i,j)-sum(dzg))

               else ! multiple fractions...
                  ! erosion/deposition rate of sand mass (m/s)
                  ! positive in case of erosion
                  edg = dzg*(1.d0-par%por)/par%dt

                  dz=>s%dzbed(i,j,:)
                  pb=>s%pbbed(i,j,:,:)

                  call update_fractions(par,s,i,j,dz,pb,edg,sum(dzg))

               endif

            enddo ! s%nx+1
         endif !s%ny>0
#ifdef USEMPI
         call xmpi_shift_ee(s%zb)
#endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Avalanching
         !


         if (par%avalanching==1) then
            do ii=1,nint(par%morfac)

               aval=.false.
               s%dzbdx=0.d0
               s%dzbdy=0.d0
               do j=1,s%ny+1
                  do i=1,s%nx
                     s%dzbdx(i,j)=(s%zb(i+1,j)-s%zb(i,j))/s%dsu(i,j)
                  enddo
               enddo
               !
               hav = s%hh
               ! Fix Hav for short wave runup:
               if (par%swrunup == 1) then
                  do j = 1,s%ny+1
                     hinterland = 0
                     do i = 1, s%nx+1
                        if (hinterland == 0 .and. s%runup(j)+s%zs(nint(s%iwl(j)),j) < s%zb(i,j) + par%eps) then
                           hinterland = 1;
                        endif
                        if (s%wetz(i,j) == 1) then
                           hav(i,j) = max(par%eps,(s%hh(i,j) + s%runup(j)))
                        elseif (hinterland == 0) then
                           hav(i,j) =  max(par%eps,s%runup(j)+s%zs(nint(s%iwl(j)),j)-s%zb(i,j) )
                        else
                           hav(i,j) = par%eps;
                        endif
                     enddo
                  enddo
               endif
               !
               do i=2,s%nx !-1 Jaap -1 gives issues for bed updating at mpi boundaries
                  do j=1,s%ny+1
                     !if (max( max(hh(i,j),par%delta*H(i,j)), max(hh(i+1,j),par%delta*H(i+1,j)) )>par%hswitch+par%eps) then
                     if(max(hav(i,j),hav(i+1,j))>par%hswitch+par%eps) then ! Jaap instead of s%hh
                        dzmax=par%wetslp;
                        ! tricks: seaward of istruct (transition from sand to structure) wetslope is set to 0.03;
                        if (i>nint(s%istruct(j))) then
                           !dzmax = 0.03d0
                           dzmax = max(0.03d0,abs(s%dzbdx(i,j))*0.99d0)
                        endif
                     else
                        dzmax=par%dryslp;
                     end if

                     if(abs(s%dzbdx(i,j))>dzmax .and. s%structdepth(i+nint(max(0.d0,sign(1.d0,s%dzbdx(i,j)))),j)>par%eps) then
                        aval=.true.
                        dzb=sign(1.0d0,s%dzbdx(i,j))*(abs(s%dzbdx(i,j))-dzmax)*s%dsu(i,j)

                        if (dzb >= 0.d0) then
                           ie = i+1                                        ! index erosion point
                           id = i                                          ! index deposition point
                           dAfac = s%dsdnzi(i,j)/s%dsdnzi(i+1,j)               ! take into account varying grid sizes
                           dzb=min(dzb,par%dzmax*par%dt/s%dsu(i,j))          ! make sure dzb is not in conflict with par%dzmax
                           dzb=min(dzb,s%structdepth(i+1,j))                 ! make sure dzb is not larger than sediment layer thickness
                        else
                           ie = i                                          ! index erosion point
                           id = i+1                                        ! index deposition point
                           dAfac = s%dsdnzi(i+1,j)/s%dsdnzi(i,j)               ! take into account varying grid sizes
                           dzb=max(dzb,-par%dzmax*par%dt/s%dsu(i,j))
                           dzb=max(dzb,-s%structdepth(i,j))
                        endif


                        if (par%ngd == 1) then ! Simple bed update in case one fraction

                           dzleft = abs(dzb)

                           s%zb(id,j) = s%zb(id,j)+dzleft*dAfac
                           s%zb(ie,j) = s%zb(ie,j)-dzleft
                           s%dzbnow(id,j) = s%dzbnow(id,j)+dzleft*dAfac ! naamgeveing?
                           s%dzbnow(ie,j) = s%dzbnow(ie,j)-dzleft
                           s%sedero(id,j) = s%sedero(id,j)+dzleft*dAfac
                           s%sedero(ie,j) = s%sedero(ie,j)-dzleft
                           s%structdepth(id,j) = max(0.d0,s%structdepth(id,j)+dzleft*dAfac)
                           s%structdepth(ie,j) = max(0.d0,s%structdepth(ie,j)-dzleft)

                           s%zs(id,j)  = s%zs(id,j)+dzleft*dAfac
                           s%zs(ie,j)  = s%zs(ie,j)-dzleft
                           s%dzav(id,j)= s%dzav(id,j)+dzleft*dAfac
                           s%dzav(ie,j)= s%dzav(ie,j)-dzleft

                        else ! multiple fractions...

                           ! now fix fractions....
                           dz => s%dzbed(ie,j,:)
                           pb => s%pbbed(ie,j,:,:)

                           ! figure out how many depth layers (ndz) are eroded in point iii
                           sdz = 0
                           ndz = 0
                           do while (sdz<abs(dzb))
                              ndz = ndz+1
                              sdz = sdz+dz(ndz)
                           enddo

                           ! now update bed and fractions by stepping through each layer seperately
                           dzleft = abs(dzb)
                           dzavt  = 0.d0

                           do jdz=1,ndz

                              dzt = min(dz(jdz),dzleft)
                              dzleft = dzleft-dzt;

                              ! erosion deposition per fraction (upwind or downwind); edg is positive in case of erosion
                              do jg=1,par%ngd
                                 edg2(jg) =  s%sedcal(jg)*dzt*pb(jdz,jg)*(1.d0-par%por)/par%dt       ! erosion    (dzt always > 0 )
                                 edg1(jg) = -s%sedcal(jg)*edg2(jg)*dAfac                             ! deposition (dzt always < 0 )
                              enddo

                              dzavt = dzavt + sum(edg2)*par%dt/(1.d0-par%por)

                              call update_fractions(par,s,ie,j,s%dzbed(ie,j,:),s%pbbed(ie,j,:,:),edg2,dzavt)           ! update bed in eroding point

                              call update_fractions(par,s,id,j,s%dzbed(id,j,:),s%pbbed(id,j,:,:),edg1,-dzavt*dAfac)    ! update bed in deposition point

                           enddo

                           ! update water levels and dzav
                           s%zs(ie,j)  = s%zs(ie,j)-dzavt
                           s%dzav(ie,j)= s%dzav(ie,j)-dzavt

                           s%zs(id,j)  = s%zs(id,j)+dzavt*dAfac
                           s%dzav(id,j)= s%dzav(id,j)+dzavt*dAfac

                        end if ! yes/no multiple fractions
                     end if
                  end do
               end do
               !JJ: update y slopes after avalanching in X-direction seems more appropriate
               do j=1,s%ny
                  do i=1,s%nx+1
                     s%dzbdy(i,j)=(s%zb(i,j+1)-s%zb(i,j))/s%dnv(i,j)
                  enddo
               enddo

               do j=2,s%ny !-1 Jaap -1 gives issues for bed updating at mpi boundaries
                  do i=1,s%nx+1
                     if(max(s%hh(i,j),s%hh(i,j+1))>par%hswitch+par%eps) then
                        dzmax=par%wetslp
                     else
                        dzmax=par%dryslp
                     end if
                     if(abs(s%dzbdy(i,j))>dzmax .and. s%structdepth(i,j+nint(max(0.d0,sign(1.d0,s%dzbdy(i,j)))))>par%eps) then ! Jaap
                        aval=.true.
                        dzb=sign(1.0d0,s%dzbdy(i,j))*(abs(s%dzbdy(i,j))-dzmax)*s%dnv(i,j)
                        !
                        if (dzb >= 0.d0) then
                           je = j+1                                        ! index erosion point
                           jd = j                                          ! index deposition point
                           dAfac = s%dsdnzi(i,j)/s%dsdnzi(i,j+1)               ! take into account varying grid sizes
                           dzb=min(dzb,par%dzmax*par%dt/s%dnv(i,j))
                           dzb=min(dzb,s%structdepth(i,j+1))
                        else
                           je = j                                          ! index erosion point
                           jd = j+1                                        ! index deposition point
                           dAfac = s%dsdnzi(i,j+1)/s%dsdnzi(i,j)               ! take into account varying grid sizes
                           dzb=max(dzb,-par%dzmax*par%dt/s%dnv(i,j))
                           dzb=max(dzb,-s%structdepth(i,j))
                        endif

                        if (par%ngd == 1) then ! Simple bed update in case one fraction

                           dzleft = abs(dzb)

                           s%zb(i,jd) = s%zb(i,jd)+dzleft*dAfac
                           s%zb(i,je) = s%zb(i,je)-dzleft
                           s%dzbnow(i,jd) = s%dzbnow(i,jd)+dzleft*dAfac ! naamgeveing?
                           s%dzbnow(i,je) = s%dzbnow(i,je)-dzleft
                           s%sedero(i,jd) = s%sedero(i,jd)+dzleft*dAfac
                           s%sedero(i,je) = s%sedero(i,je)-dzleft
                           s%structdepth(i,jd) = max(0.d0,s%structdepth(i,jd)+dzleft*dAfac)
                           s%structdepth(i,je) = max(0.d0,s%structdepth(i,je)-dzleft)

                           s%zs(i,jd)  = s%zs(i,jd)+dzleft*dAfac
                           s%zs(i,je)  = s%zs(i,je)-dzleft
                           s%dzav(i,jd)= s%dzav(i,jd)+dzleft*dAfac
                           s%dzav(i,je)= s%dzav(i,je)-dzleft

                        else ! multiple fractions...

                           dz => s%dzbed(i,je,:)
                           pb => s%pbbed(i,je,:,:)

                           ! figure out how many depth layers (ndz) are affected
                           sdz = 0
                           ndz = 0
                           do while (sdz<abs(dzb))
                              ndz = ndz+1
                              sdz = sdz+dz(ndz)
                           enddo

                           ! now update bed and fractions by stepping through each layer seperately
                           dzleft = abs(dzb)
                           dzavt  = 0.d0

                           do jdz=1,ndz
                              dzt = min(dz(jdz),dzleft)
                              dzleft = dzleft-dzt;

                              ! erosion deposition per fraction (upwind or downwind); edg is positive in case of erosion
                              do jg=1,par%ngd
                                 edg2(jg) = s%sedcal(jg)*dzt*pb(jdz,jg)*(1.d0-par%por)/par%dt        ! erosion    (dzt always > 0 )
                                 edg1(jg) = -s%sedcal(jg)*edg2(jg)*dAfac                             ! deposition (dzt always < 0 )
                              enddo

                              dzavt = dzavt + sum(edg2)*par%dt/(1.d0-par%por)

                              call update_fractions(par,s,i,je,s%dzbed(i,je,:),s%pbbed(i,je,:,:),edg2,dzavt)           ! upwind point

                              call update_fractions(par,s,i,jd,s%dzbed(i,jd,:),s%pbbed(i,jd,:,:),edg1,-dzavt*dAfac)    ! downwind point

                           enddo

                           ! update water levels and dzav
                           s%zs(i,je)  = s%zs(i,je)-dzavt
                           s%dzav(i,je)= s%dzav(i,je)-dzavt

                           s%zs(i,jd)  = s%zs(i,jd)+dzavt*dAfac
                           s%dzav(i,jd)= s%dzav(i,jd)+dzavt*dAfac

                        endif !yes/no multiple fractions
                     end if
                  end do
               end do
               if (.not.aval) exit
            end do
         else
            s%dzbdx=0.d0
            s%dzbdy=0.d0
            do j=1,s%ny+1
               do i=1,s%nx
                  s%dzbdx(i,j)=(s%zb(i+1,j)-s%zb(i,j))/s%dsu(i,j)
               enddo
            enddo
            do j=1,s%ny
               do i=1,s%nx+1
                  s%dzbdy(i,j)=(s%zb(i,j+1)-s%zb(i,j))/s%dnv(i,j)
               enddo
            enddo
         end if
         !
         ! bed boundary conditions
         !
         if(xmpi_isleft .and. s%ny>0) then
            s%zb(:,1) = s%zb(:,2)
            s%dzbdt(:,1) = s%dzbdt(:,2)
            s%dzbnow(:,1) = s%dzbnow(:,2)
            s%sedero(:,1) = s%sedero(:,2)
            s%structdepth(:,1) = s%structdepth(:,2)
            s%pbbed(:,1,:,:)=s%pbbed(:,2,:,:)
            s%z0bed(:,1)=s%z0bed(:,2)
            s%dzbed(:,1,:)=s%dzbed(:,2,:)
         endif

         if(xmpi_isright .and. s%ny>0) then
            s%zb(:,s%ny+1) = s%zb(:,s%ny)
            s%dzbdt(:,s%ny+1) = s%dzbdt(:,s%ny)
            s%dzbnow(:,s%ny+1) = s%dzbnow(:,s%ny)
            s%sedero(:,s%ny+1) = s%sedero(:,s%ny)
            s%structdepth(:,s%ny+1) = s%structdepth(:,s%ny)
            s%pbbed(:,s%ny+1,:,:)=s%pbbed(:,s%ny,:,:)
            s%z0bed(:,s%ny+1)=s%z0bed(:,s%ny)
            s%dzbed(:,s%ny+1,:)=s%dzbed(:,s%ny,:)
         endif

         ! Robert: in parallel version bed update must take place on internal boundaries:
#ifdef USEMPI
         call xmpi_shift_ee(s%zb)
#endif

         ! Update representative sed.diameter at the bed for flow friction and output
         if (par%ngd>1) then
            do j=j1,max(s%ny,1)
               do i=2,s%nx
                  s%D50top =  sum(s%pbbed(i,j,1,:)*s%D50)
                  s%D90top =  sum(s%pbbed(i,j,1,:)*s%D90)
               enddo
            enddo
         endif

      endif ! if par%t>par%morstart

   end subroutine bed_update

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine update_fractions(par,s,i,j,dz,pb,edg,dzb)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !     Copyright (C) 2009 Technische Universiteit Delft
      !        Bram van Prooijen
      !        b.c.vanprooijen@tudelft.nl
      !      +31(0)15 2784070
      !        Faculty of Civil Engineering and Geosciences
      !        department of Hydraulic Engineering
      !      PO Box 5048
      !        2600 GA Delft
      !        The Netherlands
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use params
      use spaceparams
      use xmpi_module


      implicit none

      type(spacepars),target                          :: s
      type(parameters)                                :: par

      integer                                         :: i,j,jg,jd
      real*8                                          :: ED,zbold,dzbt,fac,dzb_loc
      real*8, intent(in)                              :: dzb
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
      !TODO, dzb_loc is not initialized can be Nan, leading to infinite loop
      dzb_loc = dzb

      !!!initialize Sm

      ED = sum(edg)

      ! do t_sub=1,nt_sub  !loop over subtimesteps
      do while (abs(dzb_loc) .gt. 0.d0)
         dzbt     = min(dzb_loc,dz(par%nd_var))                 ! make sure erosion (dzg is positive) is limited to thickness of variable layer
         dzbt     = max(dzbt,-par%frac_dz*dz(par%nd_var+1))     ! make sure deposition (dzg is negative) is limited to thickness of first layer below variable layer

         fac = dzbt/dzb_loc                                     ! factor over mass change in cells to limit erosion and deposition

         dzb_loc = dzb_loc-dzbt                                 ! update dzb

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

            ! Ap makes sure that layer nd_var varies in thickness in time
            ! Ap = 0 with single fraction (in case(1))
            Ap(1) = sum(A(1,2:3)*pb(1:2,jg))
            do jd = 2,par%nd_var
               Ap(jd) = sum(A (jd,:)*pb(jd-1:jd+1,jg))
            enddo
            Ap(par%nd_var+1) = sum(A(par%nd_var+1,1:2)*pb(par%nd_var:par%nd_var+1,jg))

            ! b represents the actual erosion and deposition in the top layer.
            ! However, the thickness of the top layer remains constant and instead layer nd_var will breath in thickness
            b(1) = -edg(jg)

            !!!update Sm

            ! Sm is the sediment mass per fraction per layer
            ! Sm(1:par%nd_var+1,jg) = Sm(1:par%nd_var+1,jg) + par%dt/nt_sub*(Ap(1:par%nd_var+1)+b(1:par%nd_var+1))
            Sm(1:par%nd_var+1,jg) = Sm(1:par%nd_var+1,jg) + par%dt*fac*(Ap(1:par%nd_var+1)+b(1:par%nd_var+1))

         enddo !fractions

         ! From Sm we can compute the new fraction ratios per layer and the layer thickness...
         do jd=1,par%nd_var+1
            if (sum(Sm(jd,:))>0) then
               pb(jd,:) = Sm(jd,:)/sum(Sm(jd,:))
            else
               pb(jd,:) = pb(jd,:)
            endif
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
      s%dzbnow(i,j) = s%dzbnow(i,j)+(s%zb(i,j)-zbold)
      s%sedero(i,j) = s%sedero(i,j)+(s%zb(i,j)-zbold)
      s%structdepth(i,j) = max(0.d0,s%structdepth(i,j)+(s%zb(i,j)-zbold))

   end subroutine update_fractions

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine sedtransform(s,par)
      use params
      use spaceparams
      use readkey_module
      use xmpi_module

      implicit none

      type(spacepars),target                  :: s
      type(parameters)                        :: par

      integer                                 :: i,j,jg
      real*8                                  :: z0,dcf,dcfin,ML
      real*8                                  :: Te,Sster,cc1,cc2,wster,Ass
      real*8                                  :: kl,alpha,alpha1,alpha2,beta,psi
      real*8 , save                           :: delta,kvis,onethird,twothird,phi
      real*8 , dimension(:),allocatable    ,save     :: dster,ws0
      real*8 , dimension(:,:),allocatable  ,save     :: vmg,Cd,Asb,dhdx,dhdy,Ts,hfac
      real*8 , dimension(:,:),allocatable  ,save     :: urms2,Ucr,Ucrc,Ucrw,term1,B2,srfTotal,srfRhee,vero,Ucrb,Ucrs
      real*8 , dimension(:,:),allocatable  ,save     :: uandv,b,fslope,hloc,ceqs,ceqb,fallvelredfac
      real*8 , dimension(:,:,:),allocatable,save     :: w

      !include 's.ind'
      !include 's.inp'

      if (.not. allocated(vmg)) then
         allocate (vmg   (s%nx+1,s%ny+1))
         allocate (term1 (s%nx+1,s%ny+1))
         allocate (B2    (s%nx+1,s%ny+1))
         allocate (Cd    (s%nx+1,s%ny+1))
         allocate (Asb   (s%nx+1,s%ny+1))
         allocate (Ucr   (s%nx+1,s%ny+1))
         allocate (Ucrc  (s%nx+1,s%ny+1))
         allocate (Ucrw  (s%nx+1,s%ny+1))
         allocate (urms2 (s%nx+1,s%ny+1))
         allocate (hloc  (s%nx+1,s%ny+1))
         allocate (Ts    (s%nx+1,s%ny+1))
         allocate (ceqs  (s%nx+1,s%ny+1))
         allocate (ceqb  (s%nx+1,s%ny+1))
         allocate (srfTotal(s%nx+1,s%ny+1))       ! Lodewijk
         allocate (Ucrb  (s%nx+1,s%ny+1))         ! Lodewijk
         allocate (Ucrs  (s%nx+1,s%ny+1))         ! Lodewijk
         allocate (srfRhee(s%nx+1,s%ny+1))        ! Lodewijk
         allocate (vero  (s%nx+1,s%ny+1))         ! Lodewijk
         allocate (fallvelredfac(s%nx+1,s%ny+1))  ! Lodewijk
         allocate (w     (s%nx+1,s%ny+1,par%ngd)) ! Lodewijk
         allocate (dster (par%ngd))
         allocate (ws0   (par%ngd))
         allocate (dhdx  (s%nx+1,s%ny+1))
         allocate (dhdy  (s%nx+1,s%ny+1))
         allocate (uandv (s%nx+1,s%ny+1))
         allocate (b     (s%nx+1,s%ny+1))
         allocate (fslope(s%nx+1,s%ny+1))
         allocate (hfac  (s%nx+1,s%ny+1))
         vmg = 0.d0
         onethird=1.d0/3.d0
         twothird=2.d0/3.d0
         phi = par%reposeangle/180*par%px ! Angle of internal friction
         ! Robert: do only once, not necessary every time step
         do jg=1,par%ngd
            ! cjaap: compute fall velocity with simple expression from Ahrens (2000)
            Te    = 20.d0
            kvis  = 4.d0/(20.d0+Te)*1d-5 ! Van rijn, 1993
            Sster = s%D50(jg)/(4*kvis)*sqrt((par%rhos/par%rho-1)*par%g*s%D50(jg))
            cc1   = 1.06d0*tanh(0.064d0*Sster*exp(-7.5d0/Sster**2))
            cc2    = 0.22d0*tanh(2.34d0*Sster**(-1.18d0)*exp(-0.0064d0*Sster**2))
            wster = cc1+cc2*Sster
            ws0(jg) = wster*sqrt((par%rhos/par%rho-1.d0)*par%g*s%D50(jg))
            ! RJ: for modeling gravel
            delta = (par%rhos-par%rho)/par%rho
            dster(jg)=(delta*par%g/1.d-12)**onethird*s%D50(jg)
            if (par%fallvelred==0) then
               w(:,:,jg) = ws0(jg)
            endif
         enddo
      endif
      

      ! Lodewijk: calculate fall velocity reduction coefficient based on concentration of previous time step
      ! Rowe (1987) made an estimation of the exponent alpha by fitting a logarithmic function on a dataset of Richardson and Zaki (1954).
      if (par%fallvelred==1) then
         do jg=1,par%ngd
            alpha = 2.35d0*(2.d0+0.175d0*(ws0(jg)*s%D50(jg)/kvis)**0.75d0)/(1.d0+0.175d0*(ws0(jg)*s%D50(jg)/kvis)**0.75d0)
            fallvelredfac = (1.d0-s%ccg(:,:,jg))**(alpha)
            w(:,:,jg) = ws0(jg)*fallvelredfac
         enddo
      endif
      !
      ! hloc   = max(hh,0.01d0) ! Jaap
      hloc = max(s%hh,0.01)
      ! Compute mean fall velocity 
      par%ws=0.d0
      do jg=1,par%ngd
          par%ws=par%ws+w(1,1,jg)
      enddo
      par%ws=par%ws/par%ngd
          
      !
      ! compute near bed turbulence
      !
      ! due to short waves

      if (par%swave==1) then
         ! wave breaking induced turbulence due to short waves
         do j=1,s%ny+1
            do i=1,s%nx+1
               ! compute mixing length
               ! ML = 2*R(i,j)*par%Trep/(par%rho*c(i,j)*max(H(i,j),par%eps))
               ML = dsqrt(2*s%R(i,j)*par%Trep/(par%rho*max(s%c(i,j),1d-10)))
               ! ML = 0.9d0*H(i,j)
               ML = min(ML,hloc(i,j));
               ! exponential decay turbulence over depth
               dcfin = exp(min(100.d0,hloc(i,j)/max(ML,0.01d0)))
               dcf = min(1.d0,1.d0/(dcfin-1.d0))
               ! Jaap: new approach: compute kb based on waveturb result
               s%kb(i,j) = s%kturb(i,j)*dcf
               !
               if (par%turb == TURB_BORE_AVERAGED) then
                  s%kb(i,j) = s%kb(i,j)*par%Trep/s%Tbore(i,j)
               endif

            enddo
            ! Jaap: rundown jet creating additional turbulence
            if (par%swrunup==1)then
               s%kb(nint(s%istruct(j)),j) = s%kb(nint(s%istruct(j)),j) + par%jetfac* &
               (s%E(nint(s%istruct(j)),j)*s%strucslope(j)*sqrt(par%g/s%hh(nint(s%istruct(j)),j)))**twothird
            endif
         enddo
      elseif (par%nonh==1) then
         do j=1,s%ny+1
            do i=1,s%nx+1
               !ML=max(s%rolthick(i,j),0.07d0*max(s%hh(i,j),par%hmin))
               !s%kb(i,j) = s%kturb(i,j)*ML/max(s%hh(i,j),par%hmin) ! Simpler expression 
                s%kb(i,j) = s%kturb(i,j)  ! even simpler expression :-)
            enddo
         enddo
      endif !par%swave == 1

      ! switch to include long wave stirring
      if (par%lws==1) then
         vmg  = dsqrt(s%ue**2+s%ve**2)
      elseif (par%lws==0) then
         ! vmg lags on actual mean flow; but long wave contribution to mean flow is included...
         vmg = (1.d0-1.d0/par%cats/par%Trep*par%dt)*vmg + (1.d0/par%cats/par%Trep*par%dt)*dsqrt(s%ue**2+s%ve**2)
      endif

      urms2 = s%urms**2.d0+1.45d0*s%kb

      do jg = 1,par%ngd

         Ts       = par%tsfac*hloc/w(:,:,jg)
         s%Tsg(:,:,jg) = max(Ts,par%Tsmin)
         !
         ! calculate treshold velocity Ucr
         !
         if (par%form==FORM_SOULSBY_VANRIJN) then       ! Soulsby van Rijn
            if(s%D50(jg)<=0.0005d0) then
               Ucr=0.19d0*s%D50(jg)**0.1d0*log10(4.d0*hloc/s%D90(jg))
            else   !Dano see what happens with coarse material
               Ucr=8.5d0*s%D50(jg)**0.6d0*log10(4.d0*hloc/s%D90(jg))
            end if
         elseif (par%form==FORM_VANTHIEL_VANRIJN) then  ! Van Thiel de Vries & Reniers 2008
            if(s%D50(jg)<=0.0005) then
               Ucrc=0.19d0*s%D50(jg)**0.1d0*log10(4.d0*hloc/s%D90(jg))                           !Shields
               Ucrw=0.24d0*(delta*par%g)**0.66d0*s%D50(jg)**0.33d0*par%Trep**0.33d0          !Komar and Miller (1975)
            else if(s%D50(jg)<=0.002) then
               Ucrc=8.5d0*s%D50(jg)**0.6d0*log10(4.d0*hloc/s%D90(jg))                            !Shields
               Ucrw=0.95d0*(delta*par%g)**0.57d0*s%D50(jg)**0.43*par%Trep**0.14                !Komar and Miller (1975)
            else if(s%D50(jg)>0.002) then
               Ucrc=1.3d0*sqrt(delta*par%g*s%D50(jg))*(hloc/s%D50(jg))**(0.5d0*onethird)         !Maynord (1978) --> also Neill (1968) where 1.3d0 = 1.4d0
               Ucrw=0.95d0*(delta*par%g)**0.57d0*s%D50(jg)**0.43*par%Trep**0.14                !Komar and Miller (1975)
            end if
            B2 = vmg/max(vmg+dsqrt(urms2),par%eps)
            Ucr = B2*Ucrc + (1-B2)*Ucrw                                                     !Van Rijn 2007 (Bed load transport paper)
         end if

         ! Lodewijk: implementation Rhee (2010), reduction sediment transport due dilatancy
         srfRhee(:,:)  = 0.d0
         srfTotal(:,:) = 1.d0
         if (par%dilatancy == 1) then
            vero(:,:) = max(0.d0,-1.d0*s%dzbdt)      ! Erosion velocity, for now asume it equal to -s%dzbdt of the previous time step
            kl = par%g/(160.d0*kvis)*(s%D15(jg)**2.d0)*((par%por**3.d0)/(1.d0-par%por)**2.d0) ! Permeability, Adel 1987, which is based on the Ergun equation
            ! bed porosity, estimated here as maximum porosity, to be added to user input
            ! A=3/4 for single particles and A=1/(1-n0) for a continuum
            ! Reduction factor on the critical Shields parameter by dilatancy (Van Rhee, 2010)
            srfRhee(:,:) = vero(:,:)/kl*(par%pormax-par%por)/(1.d0-par%pormax)*par%rheeA/delta
         endif
         ! Lodewijk: implementation bed slope reduction on critical flow velocity
         if (par%bdslpeffini == BDSLPEFFINI_NONE) then
            srfTotal(:,:) = 1.d0 + srfRhee(:,:)
         elseif (par%bdslpeffini == BDSLPEFFINI_TOTAL .or. par%bdslpeffini == BDSLPEFFINI_BED) then
            do j=1,s%ny+1
               do i=1,s%nx+1
                  ! Prevent NaN values if too small values
                  if  ((dabs(s%ue(i,j)) > 0.000001d0 .or. dabs(s%ve(i,j)) > 0.000001d0) .and.  &
                  (dabs(s%dzbdx(i,j))>0.000001d0 .or. dabs(s%dzbdy(i,j))>0.000001d0)) then
                     ! Angle between the x-axis and the flow velocity
                     ! REMARK: also include waves in the velocity?
                     if (s%ue(i,j) < 0.d0) then
                        alpha1 = datan(s%ve(i,j)/s%ue(i,j)) + par%px
                     else
                        alpha1 = datan(s%ve(i,j)/s%ue(i,j))
                     endif
                     ! Angle between the x-axis and the bed slope vector directed in down-slope direction, derived in thesis Lodewijk
                     if (s%dzbdy(i,j) >= 0.d0) then
                        alpha2 = -datan(s%dzbdx(i,j)/s%dzbdy(i,j))+1.5d0*par%px
                     else
                        alpha2 = -datan(s%dzbdx(i,j)/s%dzbdy(i,j))+0.5d0*par%px
                     endif
                     psi = alpha1-(alpha2-par%px) ! Angle between the flow direction and the up-slope directed vector
                     if (dabs(s%dzbdx(i,j))<0.000001d0) then ! A smaller slope could result in a NaN for beta
                        ! Beta purely based on dzbdy
                        beta = datan(dabs(s%dzbdy(i,j)))
                     else
                        beta = datan(dabs(s%dzbdx(i,j)/dsin(datan(s%dzbdx(i,j)/s%dzbdy(i,j)))))     ! Maximum absolute bed slope angle, derived in thesis Lodewijk
                     endif
                     beta = min(beta,phi) ! Take min to exclude NaN's
                     if (par%dilatancy == 1) then
                        srfTotal(i,j) = (dcos(psi)*dsin(beta)+dsqrt( &
                        ((srfRhee(i,j))**2+2*srfRhee(i,j)*dcos(beta)+dcos(beta)**2 &
                        )*dtan(phi)**2-dsin(psi)**2*dsin(beta)**2 &
                        ))/dtan(phi) ! Soulsby (1997), modified by Lodewijk (see Thesis)
                     else
                        srfTotal(i,j) = (dcos(psi)*dsin(beta)+ &
                        dsqrt(dcos(beta)**2*dtan(phi)**2-dsin(psi)**2*dsin(beta)**2))/dtan(phi) ! Soulsby (1997)
                     endif
                  endif
               end do
            end do
         endif
         ! Calculate the new critical velocity based on the modification factors on the Shields parameter, bed slope only on bed load
         Ucrb(:,:) = Ucr(:,:)*sqrt(srfTotal) ! Lodewijk
         if (par%bdslpeffini == BDSLPEFFINI_TOTAL) then
            Ucrs = Ucrb
         else
            Ucrs = Ucr*(1.d0+sqrt(srfRhee)) ! Lodewijk, no effect on suspended load by bed slope
         endif


         if (par%form==FORM_SOULSBY_VANRIJN) then       ! Soulsby van Rijn
            !
            ! drag coefficient
            z0 = par%z0
            Cd=(0.40d0/(log(max(hloc,10.d0*z0)/z0)-1.0d0))**2
            !
            ! transport parameters
            Asb=0.005d0*hloc*(s%D50(jg)/hloc/(delta*par%g*s%D50(jg)))**1.2d0         ! bed load coefficent
            Ass=0.012d0*s%D50(jg)*dster(jg)**(-0.6d0)/(delta*par%g*s%D50(jg))**1.2d0 ! suspended load coeffient
            !
            term1=(vmg**2+0.018/Cd*par%sws*urms2) ! Make 0.018/Cd is always smaller than the flow friction coefficient
            !
            ! reduce sediment suspensions for (inundation) overwash conditions with critical flow velocitties
            ! vmag2=min(vmag2,par%smax*par%C**2*D50(jg)*delta)
            ! vmag2=min(vmag2,par%smax*par%g/par%cf*D50(jg)*delta)            ! In terms of cf
            ! term1=sqrt(vmag2+0.018d0/Cd*urms2)     ! nearbed-velocity
            ! the two above lines are comment out and replaced by a limit on total velocity u2+urms2, robert 1/9 and ap 28/11
            !
            term1=min(term1,par%smax*par%g/s%cf*s%D50(jg)*delta)
            term1=sqrt(term1)
            !
            ceqb = 0.d0*term1                                                                     !initialize ceqb
            ceqs = 0.d0*term1                                                                     !initialize ceqs
            do j=1,s%ny+1
               do i=1,s%nx
                  ! Lodewijk: sepperate bed load from suspended load since Ucr not anymore the same
                  if(term1(i,j)>Ucrb(i,j) .and. hloc(i,j)>par%eps) then
                     ceqb(i,j)=Asb(i,j)*(term1(i,j)-Ucrb(i,j))**2.4d0
                  end if
                  if(term1(i,j)>Ucrs(i,j) .and. hloc(i,j)>par%eps) then
                     ceqs(i,j)=Ass*(term1(i,j)-Ucrs(i,j))**2.4d0
                  end if
               end do
            end do
            !
         elseif (par%form==FORM_VANTHIEL_VANRIJN) then  ! Van Thiel de Vries & Reniers 2008
            !
            ! transport parameters
            Asb=0.015d0*hloc*(s%D50(jg)/hloc)**1.2d0/(delta*par%g*s%D50(jg))**0.75d0        !bed load coefficent
            Ass=0.012d0*s%D50(jg)*dster(jg)**(-0.6d0)/(delta*par%g*s%D50(jg))**1.2d0        !suspended load coeffient
            !
            ! Jaap: par%sws to set short wave stirring to zero
            ! Jaap: Van Rijn use Peak orbital flow velocity --> 0.64 corresponds to 0.4 coefficient regular waves Van Rijn (2007)
            term1=vmg**2+0.64d0*par%sws*urms2
            ! reduce sediment suspensions for (inundation) overwash conditions with critical flow velocitties
            term1=min(term1,par%smax*par%g/max(s%cf,1d-10)*s%D50(jg)*delta)
            term1=sqrt(term1)
            !
            ceqb = 0.d0*term1                                                                     !initialize ceqb
            ceqs = 0.d0*term1                                                                     !initialize ceqs
            do j=1,s%ny+1
               do i=1,s%nx
                  ! Lodewijk: sepperate bed load from suspended load since Ucr not anymore the same
                  if(term1(i,j)>Ucrb(i,j) .and. hloc(i,j)>par%eps) then
                     ceqb(i,j)=Asb(i,j)*(term1(i,j)-Ucrb(i,j))**1.5
                  end if
                  if(term1(i,j)>Ucrs(i,j) .and. hloc(i,j)>par%eps) then
                     ceqs(i,j)=Ass*(term1(i,j)-Ucrs(i,j))**2.4
                  end if
               end do
            end do
            !
         end if
         !
         ceqb = min(ceqb/hloc,par%cmax/2) ! maximum equilibrium bed concentration
         s%ceqbg(:,:,jg) = (1-par%bulk)*ceqb*s%sedcal(jg)*s%wetz
         ceqs = min(ceqs/hloc,par%cmax/2) ! maximum equilibrium suspended concentration
         s%ceqsg(:,:,jg) = (ceqs+par%bulk*ceqb)*s%sedcal(jg)*s%wetz
         ! Jaap: old brute method to prevent strong coastline erosion
         ! where (hloc<=par%hmin) ceqsg(:,:,jg) = 0.d0
      enddo                                 ! end of grain size classes
      ! end of grain size classes

   end subroutine sedtransform

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine waveturb(s,par)
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

      implicit none

      type(spacepars),target                   :: s
      type(parameters)                         :: par

      integer                                  :: i
      integer                                  :: j
      real*8                                   :: ML, disturb
      real*8, save                             :: twothird
      real*8,dimension(:,:),allocatable,save   :: ksource, kturbu,kturbv,Sturbu,Sturbv,dzsdt_cr

      !include 's.ind'
      !include 's.inp'

      if (.not. allocated(kturbu)) then
         allocate(ksource (s%nx+1,s%ny+1))
         allocate(kturbu (s%nx+1,s%ny+1))
         allocate(kturbv (s%nx+1,s%ny+1))
         allocate(Sturbu (s%nx+1,s%ny+1))
         allocate(Sturbv (s%nx+1,s%ny+1))
         allocate(dzsdt_cr (s%nx+1,s%ny+1))
         twothird=2.d0/3.d0
      endif
      ! use lagrangian velocities
      kturbu       = 0.0d0  !Jaap
      kturbv       = 0.0d0  !Jaap
!     dzsdt_cr=par%beta*s%c
      dzsdt_cr=par%beta*sqrt(par%g*s%hh)
      ! Update roller thickness
      !s%rolthick=s%rolthick+par%dt*(abs(s%dzsdt)-dzsdt_cr) Dano: not abs!
      s%rolthick=s%rolthick+par%dt*(s%dzsdt-dzsdt_cr)
      s%rolthick=max(s%rolthick,0.d0)
      s%rolthick=min(s%rolthick,s%hh)
      where (s%wetz==0)
          s%rolthick=0.d0
      endwhere

      ! Jaap compute sources and sinks
      ! long wave source
      ksource = 0
      if (par%lwt == 1) then
         ksource = par%g*s%rolthick*par%beta*sqrt(par%g*s%hh)      !only important in shallow water, where s%c=sqrt(gh)
      endif
      if (par%turbadv == TURBADV_NONE) then
         s%kturb = (s%DR/par%rho)**twothird           ! See Battjes, 1975 / 1985
      else
         ksource = ksource + s%DR/par%rho


         ! Jaap do long wave turb approach for both short waves and long waves
         !
         ! Turbulence in uu-points
         !
         do j=1,s%ny+1
            do i=1,s%nx
               if(s%uu(i,j)>0.d0) then
                  kturbu(i,j)=par%thetanum*s%kturb(i,j)+(1.d0-par%thetanum)*s%kturb(min(i+1,s%nx),j)
               elseif(s%uu(i,j)<0.d0) then
                  kturbu(i,j)=par%thetanum*s%kturb(i+1,j)+(1.d0-par%thetanum)*s%kturb(max(i,2),j)
               else
                  kturbu(i,j)=0.5d0*(s%kturb(i,j)+s%kturb(i+1,j))
               endif
            enddo
         enddo
         if (xmpi_isbot) kturbu(s%nx+1,:) = s%kturb(s%nx+1,:)
         !
         ! Turbulence in vv-points
         !
         do j=1,s%ny
            do i=1,s%nx+1
               if(s%vv(i,j)>0) then
                  kturbv(i,j)=par%thetanum*s%kturb(i,j)+(1.d0-par%thetanum)*s%kturb(i,min(j+1,s%ny))
               else if(s%vv(i,j)<0) then
                  kturbv(i,j)=par%thetanum*s%kturb(i,j+1)+(1.d0-par%thetanum)*s%kturb(i,max(j,2))
               else
                  kturbv(i,j)=0.5d0*(kturbv(i,j)+kturbv(i,j+1))
               endif
            enddo
         enddo
         kturbv(:,s%ny+1) = s%kturb(:,s%ny+1) !Robert
         !
         ! Turbulence advection in X and Y direction
         !
         if (par%turbadv == TURBADV_LAGRANGIAN) then
            Sturbu=kturbu*s%uu*s%hu*s%wetu
            Sturbv=kturbv*s%vv*s%hv*s%wetv
         elseif (par%turbadv == TURBADV_EULERIAN) then
            Sturbu=kturbu*s%ueu*s%hu*s%wetu
            Sturbv=kturbv*s%vev*s%hv*s%wetv
         endif
         !
         ! Update turbulence
         !
         if (s%ny>0) then
            do j=2,s%ny+1
               do i=2,s%nx+1
                  if (par%betad>0) then
                     disturb=par%betad*s%kturb(i,j)**1.5d0
                  else
                     ML=max(s%rolthick(i,j),0.07d0*max(s%hh(i,j),par%hmin))
                     disturb=0.08d0*max(s%hh(i,j),par%hmin)/ML*s%kturb(i,j)**1.5d0
                  endif
                  s%kturb(i,j) = (s%hold(i,j)*s%kturb(i,j)-par%dt*(       &
                  (Sturbu(i,j)*s%dnu(i,j)-Sturbu(i-1,j)*s%dnu(i-1,j)+&
                  Sturbv(i,j)*s%dsv(i,j)-Sturbv(i,j-1)*s%dsv(i,j-1))*s%dsdnzi(i,j)-&
                  (ksource(i,j)-disturb)))/max(s%hh(i,j),par%hmin)
                  s%kturb(i,j)=max(s%kturb(i,j),0.0d0)

               enddo
            enddo
         else
            j=1
            do i=2,s%nx+1

               s%kturb(i,j) = (s%hold(i,j)*s%kturb(i,j)-par%dt*(       &
               (Sturbu(i,j)*s%dnu(i,j)-Sturbu(i-1,j)*s%dnu(i-1,j))*s%dsdnzi(i,j)-&
               (ksource(i,j)-par%betad*s%kturb(i,j)**1.5d0)))/max(s%hh(i,j),par%hmin)
               s%kturb(i,j)=max(s%kturb(i,j),0.0d0)

            enddo
         endif

         s%kturb = s%kturb/max(s%hh,0.01d0)

         ! Jaap only required for advection mode?
         if (xmpi_istop) s%kturb(1,:)=s%kturb(2,:)
         if (xmpi_isbot) s%kturb(s%nx+1,:)=s%kturb(s%nx+1-1,:)
         if (s%ny>0) then
            if (xmpi_isleft)  s%kturb(:,1)=s%kturb(:,2)
            if (xmpi_isright) s%kturb(:,s%ny+1)=s%kturb(:,s%ny+1-1)
         endif

#ifdef USEMPI
         call xmpi_shift_ee(s%kturb)
#endif

      endif ! turbadv == 'none'

   end subroutine waveturb

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine RvR(s,par)

      use params
      use spaceparams
      use xmpi_module

      implicit none

      type(spacepars),target                   :: s
      type(parameters)                         :: par

      real*8 , save                            :: m1,m2,m3,m4,m5,m6,alpha,beta

      real*8 , dimension(:,:),allocatable,save   :: Urs,Bm,B1

      !include 's.ind'
      !include 's.inp'

      ! only in first timestep..
      if (.not. allocated(Urs)) then

         allocate (Urs    (s%nx+1,s%ny+1))
         allocate (Bm    (s%nx+1,s%ny+1))
         allocate (B1    (s%nx+1,s%ny+1))

         m1 = 0;       ! a = 0
         m2 = 0.7939;  ! b = 0.79 +/- 0.023
         m3 = -0.6065; ! c = -0.61 +/- 0.041
         m4 = 0.3539;  ! d = -0.35 +/- 0.032
         m5 = 0.6373;  ! e = 0.64 +/- 0.025
         m6 = 0.5995;  ! f = 0.60 +/- 0.043
         alpha = -log10(exp(1.d0))/m4
         beta  = exp(m3/m4)

      endif

      Urs = 3.d0/8.d0*sqrt(2.d0)*s%H*s%k/(s%k*s%hh)**3                    !Ursell number
      Urs = max(Urs,0.000000000001d0)
      Bm = m1 + (m2-m1)/(1.d0+beta*Urs**alpha)                    !Boltzmann sigmoid (eq 6)
      B1 = (-90.d0+90.d0*tanh(m5/Urs**m6))*par%px/180.d0
      s%Sk = Bm*cos(B1)                                            !Skewness (eq 8)
      s%As = Bm*sin(B1)                                            !Asymmetry(eq 9)
      s%ua = par%sws*(par%facSk*s%Sk-par%facAs*s%As)*s%urms

      ! multiply Sk and As with wetz to get zeros at h = 0 for output
      s%Sk = s%Sk*s%wetz
      s%As = s%As*s%wetz

   end subroutine RvR

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine vT(s,par)

      use params
      use spaceparams
      use readkey_module
      use xmpi_module

      implicit none

      type(spacepars),target                   :: s
      type(parameters)                         :: par

      integer                                  :: i,j
      integer , save                           :: nh,nt
      integer                                  :: ih0,it0,ih1,it1
      real*8                                   :: p,q,f0,f1,f2,f3
      real*8                                   :: t0fac,siguref,duddtmax,dudtmax,duddtmean,dudtmean,detadxmean
      real*8 , save                            :: dh,dt

      real*8 , dimension(:,:),allocatable  ,save     :: h0,t0,detadxmax
      ! Robert: RF table now included in source code, rather than read from file
      ! Rienecker Fenton table with amongst others amplitudes non-linear components obtained with stream function theory
      include 'RF.inc'
      ! Robert: 'RF.inc' contains definition of RF as real*8(5,33,40) with "parameter" attribute
      !   so RF values may not be modified! To save memory, only rows 13,14,15,16 and 18 of the
      !   original matrix are stored. So new row 1 corresponds with old row 13, etc.

      !include 's.ind'
      !include 's.inp'


      ! only in first timestep..
      if (.not. allocated(h0)) then
         allocate (h0    (s%nx+1,s%ny+1))
         allocate (t0    (s%nx+1,s%ny+1))
         allocate (detadxmax    (s%nx+1,s%ny+1))
         dh = 0.03d0
         dt = 1.25d0
         nh = floor(0.99d0/dh);
         nt = floor(50.d0/dt);
      endif

      ! non-linearity of short waves is listed in table as function of dimensionless wave height h0 and dimensionless wave period t0

      ! compute dimensionless wave height and wave period in each grid point..
      h0 = min(nh*dh,max(dh,min(s%H,s%hh)/s%hh))
      t0 = min(nt*dt,max(dt,par%Trep*sqrt(par%g/s%hh)))

      ! estimate Sk, As and ua by interpolating table values
      do j=1,s%ny+1
         do i=1,s%nx+1
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

            ! Skewness and assymetry
            s%Sk(i,j) = f0*RF(1,ih0,it0)+f1*RF(1,ih1,it0)+ f2*RF(1,ih0,it1)+f3*RF(1,ih1,it1)
            s%As(i,j) = f0*RF(2,ih0,it0)+f1*RF(2,ih1,it0)+ f2*RF(2,ih0,it1)+f3*RF(2,ih1,it1)

            ! Sediment advection velocity from Skewness and Assymetry
            ! ua(i,j) = par%sws*par%facua*(Sk(i,j)-As(i,j))*urms(i,j)
            s%ua(i,j) = par%sws*(par%facSk*s%Sk(i,j)-par%facAs*s%As(i,j))*s%urms(i,j)

            ! Estimate bore period Tbore and mean slope bore front to feeded back in roller energy balance

            ! correct slope in case 1.25>T0>50
            if (t0(i,j)==50.d0) then
               t0fac = 50.d0/max((par%Trep*sqrt(par%g/s%hh(i,j))),50.d0)
            elseif (t0(i,j)==1.25)then
               t0fac = 1.25d0/min((par%Trep*sqrt(par%g/s%hh(i,j))),1.25d0)
            else
               t0fac = 1.d0
            endif

            ! detadxmax for Tbore...
            ! dimnesionless maximum acceleration under bore front
            duddtmax = f0*RF(3,ih0,it0)+f1*RF(3,ih1,it0)+ f2*RF(3,ih0,it1)+f3*RF(3,ih1,it1)
            siguref = f0*RF(4,ih0,it0)+f1*RF(4,ih1,it0)+ f2*RF(4,ih0,it1)+f3*RF(4,ih1,it1)
            ! translate dimensionless duddtmax to real world dudtmax
            !         /scale with variance and go from [-] to [m/s^2]     /tableb./dimensionless dudtmax
            dudtmax = s%urms(i,j)/max(par%eps,siguref)*sqrt(par%g/s%hh(i,j))*t0fac*duddtmax
            detadxmax(i,j) = dudtmax*sinh(s%k(i,j)*s%hh(i,j))/max(max(s%c(i,j),sqrt(s%H(i,j)*par%g)),1d-10)/s%sigm(i,j)

            ! detadxmean for roller energy balance dissipation...
            if (par%rfb==1) then
               duddtmean = f0*RF(5,ih0,it0)+f1*RF(5,ih1,it0)+ f2*RF(5,ih0,it1)+f3*RF(5,ih1,it1)
               dudtmean = s%urms(i,j)/max(par%eps,siguref)*sqrt(par%g/s%hh(i,j))*t0fac*duddtmean
               detadxmean = dudtmean*sinh(s%k(i,j)*s%hh(i,j))/max(s%c(i,j),sqrt(s%H(i,j)*par%g))/s%sigm(i,j)
               s%BR(i,j) = par%BRfac*sin(atan(detadxmean))
            endif

         enddo
      enddo

      s%Tbore = max(par%Trep/25.d0,min(par%Trep/4.d0,s%H/(max(max(s%c,sqrt(s%H*par%g)),1d-10)*max(detadxmax,par%eps))))
      s%Tbore = par%Tbfac*s%Tbore

      ! multiply Sk and As with wetz to get zeros at h = 0 for output
      s%Sk = s%Sk*s%wetz
      s%As = s%As*s%wetz

   end subroutine vT

   subroutine setbathy_update(s, par)

      use params
      use spaceparams
      use interp

      implicit none

      type(spacepars)                     :: s
      type(parameters)                    :: par

      integer                             :: i,j,dummy
      real*8,dimension(s%nx+1,s%ny+1)     :: zbnew

      ! interpolate from file
      do j=1,s%ny+1
         do i=1,s%nx+1
            call LINEAR_INTERP(s%tsetbathy,s%setbathy(i,j,:),par%nsetbathy, &
            par%t,zbnew(i,j),dummy)
         enddo
      enddo
      ! update water level
      s%zs = s%zs+zbnew-s%zb
      ! update bed level
      s%zb = zbnew
      ! update wet and dry cells
      where (s%zs<s%zb+par%eps)
         s%wetz=0
         s%zs = s%zb+par%eps
         s%hh = par%eps
      endwhere

   end subroutine setbathy_update

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine hybrid(s,par)

      use params
      use interp
      use spaceparams
      use xmpi_module

      implicit none

      type(spacepars),target                   :: s
      type(parameters)                         :: par

      integer                                  :: i,j,j1,indx,first, nIter, maxIter
      integer , dimension(:), allocatable,save :: slopeind
      real*8 , dimension(:), allocatable,save  :: hav1d
      real*8                                   :: irrb,runup_old

      !include 's.ind'
      !include 's.inp'

      if (.not. allocated(hav1d)) then
         allocate(hav1d (s%nx+1))
         allocate(slopeind (s%nx+1))
      endif

      do j=1,s%ny+1
         indx = s%nx+1
         first = 0
         do i=1,s%nx
            if (s%wetz(i,j)-s%wetz(max(i-1,1),j)==-1 .and. first==0) then ! transition from wet to dry
               ! only consider first dry point
               first = 1
               ! find wave height for runup at L1 meter from water line
               call linear_interp(s%xz(:,j),s%H(:,j),s%nx+1,s%xz(i-1,j)-s%L1(i-1,j),s%Hrunup(j),indx)
               ! Find toe of runup slope if present (dzbdx > 0.15).
               ! If not present Hrunup will converge to H at the water line (where H = 0 per definition)
               do j1=indx,i-1
                  ! TODO: Strange condition. In case of no toe, or no steep slope,
                  !   there will still be extra turbulence at L1 meter from the water line...
                  ! cross shore location structure toe
                  if (s%dzbdx(j1,j)<0.15d0 .or. s%structdepth(j1,j)>0.1d0) then
                     indx = j1
                  endif
               enddo
               ! update Hrunup and runup x-location
               s%Hrunup(j) = s%H(indx,j)       ! short wave height at revetment toe
               s%xHrunup(j) = s%xz(indx,j)     ! cross-shore location revetment toe
               s%istruct(j) = indx*1.d0      ! cross-shore index revetment toe
               s%iwl(j) = (i-1)*1.d0         ! cross-shore location waterline (inlcuding lw-s%runup)

               ! now iteratively compute runup
               hav1d = s%hh(:,j)
               runup_old = huge(0.d0)
               s%runup(j) = 0;
               nIter = 0;
               maxIter = 50;
               do while (abs(s%runup(j)-runup_old)>0.01d0 .and. nIter < maxIter)
                  nIter = nIter +1;
                  runup_old = s%runup(j)
                  slopeind = 0
                  where (hav1d>par%eps .and. s%dzbdx(:,j)>0.15d0)
                     slopeind = 1
                  endwhere
                  s%strucslope(j) = sum(s%dzbdx(indx:s%nx,j)*s%dsu(indx:s%nx,j)*slopeind(indx:s%nx))/ &
                  max(par%eps,sum(s%dsu(indx:s%nx,j)*slopeind(indx:s%nx)))
                  if (s%strucslope(j) > 0.d0) then
                     irrb = s%strucslope(j)/sqrt(2*par%px*max(s%Hrunup(j),par%eps)/par%g/par%Trep**2)
                     s%runup(j) = par%facrun*min(irrb,2.3d0)*s%Hrunup(j)*cos(2*par%px/par%Trep*par%t)
                  else
                     s%runup(j) = 0.d0;
                  endif
                  ! This triggers s%runup to be calculated for the complete hinterland (also behind a dune).
                  ! Only values calculated for the first dune front are used in the avalanching algorithm (morphevolution)
                  hav1d =  s%wetz(:,j)*max(par%eps,(s%hh(:,j) + s%runup(j))) + &
                  (1.d0-s%wetz(:,j))*max(par%eps,s%runup(j)+s%zs(i-1,j)-s%zb(:,j) )
               enddo
            endif
         enddo
      enddo



   end subroutine hybrid

end module morphevolution
