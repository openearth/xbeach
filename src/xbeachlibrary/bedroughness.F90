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
module bedroughness_module

   implicit none
   save

   real*8,dimension(:,:),allocatable,private  :: kru,krv,kru50,krv50,kru90,krv90
   real*8,dimension(:,:),allocatable,private  :: urms_upd,u2_upd
   real*8,dimension(:,:),allocatable,private  :: facbl,blphi,infilb,Ubed,Ventilation
   real*8,dimension(:,:),allocatable,private  :: ueuf,uevf,vevf,veuf
   real*8,dimension(:,:),allocatable,private  :: ueuold,uevold,vevold,veuold
   real*8,dimension(:,:),allocatable,private  :: dudtsmooth,dvdtsmooth
   real*8,dimension(:,:),allocatable,private  :: shieldsu,shieldsv
   real*8,dimension(2),              private  :: dtold
   real*8,                           private  :: delta,rhogdelta


   public bedroughness_init
   public bedroughness_update

contains

   subroutine bedroughness_init(s,par)
      use params
      use spaceparams
      use paramsconst

      implicit none

      type(parameters),intent(in)                 :: par
      type(spacepars),target                      :: s
      integer                                     :: i,j

      !include 's.ind'
      !include 's.inp'

      s%taubx_add = 0.d0
      s%tauby_add = 0.d0
      
      select case (par%bedfriction)

       case(BEDFRICTION_CHEZY)
         s%cfu = par%g/s%bedfriccoef**2
         s%cfv = par%g/s%bedfriccoef**2
       case(BEDFRICTION_CF)
         s%cfu = s%bedfriccoef
         s%cfv = s%bedfriccoef
       case (BEDFRICTION_MANNING)
         where(s%wetu==1)
            s%cfu = par%g*s%bedfriccoef**2/s%hu**(1.d0/3)
            s%cfu = min(s%cfu,0.1d0)
         elsewhere
            s%cfu = 0.1d0
         endwhere
         where(s%wetv==1)
            s%cfv = par%g*s%bedfriccoef**2/s%hv**(1.d0/3)
            s%cfv= min(s%cfv,0.1d0)
         elsewhere
            s%cfv = 0.1d0
         endwhere
       case (BEDFRICTION_WHITE_COLEBROOK_GRAINSIZE)
         allocate (kru(s%nx+1,s%ny+1))
         allocate (krv(s%nx+1,s%ny+1))
         ! compute initial D90 at cell interface
         if (par%ngd==1) then
            kru = 3.d0*s%D90top
            krv = 3.d0*s%D90top
         else
            ! Rougness height in u-points
            do j=1,s%ny+1
               do i=1,s%nx
                  kru(i,j) = 1.5d0*(s%D90top(i,j)+s%D90top(i+1,j))
               enddo
            enddo
            kru(s%nx+1,:) = kru(s%nx,:)
            ! Roughness height in v-points
            if (par%ny==0) then
               ! v-point is central on top of z-point
               krv = 3*s%D90top
            else
               do j=1,s%ny
                  do i=1,s%nx+1
                     krv(i,j) = 1.5d0*(s%D90top(i,j)+s%D90top(i,j+1))
                  enddo
               enddo
               krv(:,s%ny+1) = krv(:,s%ny)
            endif
         endif
         where(s%wetu==1)
            s%cfu = par%g/(18*log10(12*max(s%hu,kru)/kru))**2
            s%cfu = min(s%cfu,0.1d0)
         elsewhere
            s%cfu = 0.1d0
         endwhere
         where(s%wetv==1)
            s%cfv = par%g/(18*log10(12*max(s%hv,krv)/krv))**2
            s%cfv = min(s%cfv,0.1d0)
         elsewhere
            s%cfv = 0.1d0
         endwhere
       case (BEDFRICTION_WHITE_COLEBROOK)
         allocate (kru(s%nx+1,s%ny+1))
         allocate (krv(s%nx+1,s%ny+1))
         kru = s%bedfriccoef
         krv = s%bedfriccoef
         where(s%wetu==1)
            s%cfu = par%g/(18*log10(12*max(s%hu,kru)/kru))**2
            s%cfu = min(s%cfu,0.1d0)
         elsewhere
            s%cfu = 0.1d0
         endwhere
         where(s%wetv==1)
            s%cfv = par%g/(18*log10(12*max(s%hv,krv)/krv))**2
            s%cfv = min(s%cfv,0.1d0)
         elsewhere
            s%cfv = 0.1d0
         endwhere
      end select



   end subroutine bedroughness_init


   subroutine bedroughness_update(s,par)

      use params
      use spaceparams
      use paramsconst

      implicit none

      type(parameters),intent(in)                 :: par
      type(spacepars),target                      :: s
      integer                                     :: i,j

      select case (par%bedfriction)

       case(BEDFRICTION_CHEZY)
         ! do nothing, this is constant bed roughness, already set in initialise
       case(BEDFRICTION_CF)
         ! do nothing, this is constant bed roughness, already set in initialise
       case (BEDFRICTION_MANNING)
         ! C = H**(1/6)/n
         ! cf = g/C**2 = g/(hu**(1/6)/n)**2 = g*n**2/hu**(1/3)
         where(s%wetu==1)
            s%cfu = par%g*s%bedfriccoef**2/s%hu**(1.d0/3)
            s%cfu = min(s%cfu,0.1d0)
         elsewhere
            s%cfu = 0.1d0
         endwhere
         where(s%wetv==1)
            s%cfv = par%g*s%bedfriccoef**2/s%hv**(1.d0/3)
            s%cfv = min(s%cfv,0.1d0)
         elsewhere
            s%cfv = 0.1d0
         endwhere
       case (BEDFRICTION_WHITE_COLEBROOK_GRAINSIZE)
         if (par%ngd>1) then
            ! Rougness height in u-points
            do j=1,s%ny+1
               do i=1,s%nx
                  if (s%wetu(i,j)==1) then
                     kru(i,j) = 1.5d0*(s%D90top(i,j)+s%D90top(i+1,j))
                  endif
               enddo
            enddo
            kru(s%nx+1,:) = kru(s%nx,:)
            ! Roughness height in v-points
            if (par%ny==0) then
               ! v-point is central on top of z-point
               where (s%wetv==1)
                  krv = 3*s%D90top
               endwhere
            else
               do j=1,s%ny
                  do i=1,s%nx+1
                     if (s%wetv(i,j)==1) then
                        krv(i,j) = 1.5d0*(s%D90top(i,j)+s%D90top(i,j+1))
                     endif
                  enddo
               enddo
               krv(:,s%ny+1) = krv(:,s%ny)
            endif
         endif
         where(s%wetu==1)
            s%cfu = par%g/(18*log10(12*max(s%hu,kru)/kru))**2
            s%cfu = min(s%cfu,0.1d0)
         elsewhere
            s%cfu = 0.1d0
         endwhere
         where(s%wetv==1)
            s%cfv = par%g/(18*log10(12*max(s%hv,krv)/krv))**2
            s%cfv = min(s%cfv,0.1d0)
         elsewhere
            s%cfv = 0.1d0
         endwhere
       case (BEDFRICTION_WHITE_COLEBROOK)
         where(s%wetu==1)
            s%cfu = par%g/(18*log10(12*max(s%hu,kru)/kru))**2
            s%cfu = min(s%cfu,0.1d0)
         elsewhere
            s%cfu = 0.1d0
         endwhere
         where(s%wetv==1)
            s%cfv = par%g/(18*log10(12*max(s%hv,krv)/krv))**2
            s%cfv = min(s%cfv,0.1d0)
         elsewhere
            s%cfv = 0.1d0
         endwhere
      end select

      ! Additional bed friction terms due to acceleration
      s%taubx_add = 0.d0
      s%tauby_add = 0.d0
   
      ! Increased bed shear due to acceleration terms (phase-shift type or Morison type)
      if (par%friction_acceleration==CF_ACC_NIELSEN) then
         call acceleration_boundary_layer_effect_Nielsen(s,par)
      elseif (par%friction_acceleration==CF_ACC_MCCALL) then
         call acceleration_boundary_layer_effect(s,par)
      endif
      
      ! turbulence and acceleration add to shear stress, but do
      ! not necessarily also need enhanced cf from infiltration
      if (par%friction_turbulence==1) then
         call turbulence_boundary_layer_effect(s,par)   
      endif
      
      ! Infiltration enhances normal bed shear stress through 
      ! increase of cfu and cfv
      if (par%friction_infiltration==1) then
         call infiltration_boundary_layer_effect(s,par)
      endif
      

      s%cfu = min(s%cfu,1.d0)    ! check why needed
      s%cfv = min(s%cfv,1.d0)


   end subroutine bedroughness_update

   subroutine infiltration_boundary_layer_effect(s,par)

      use params
      use spaceparams
   
      IMPLICIT NONE

      type(parameters),intent(in)                 :: par
      type(spacepars)                             :: s
      integer                                     :: i,j
      real*8,parameter                            :: epsVentilation = 1.d0
      real*8,parameter                            :: maxEnhancement1 = 3.d0 ! Pers. Corr. Conley 2014
      real*8,parameter                            :: maxEnhancement2 = 0.1d0 ! Pers. Corr. Conley 2014
   
   
      if (.not.allocated(facbl)) then
         ! arrays available in bed friction module only
         allocate (facbl(s%nx+1,s%ny+1))
         allocate (blphi(s%nx+1,s%ny+1))
         allocate (infilb(s%nx+1,s%ny+1))
         allocate (Ubed(s%nx+1,s%ny+1))
         allocate (Ventilation(s%nx+1,s%ny+1))
      endif
   
   
      ! Include infiltration effect on boundary layer thinning
      ! Taken from eq 13 Butt, Russell, Turner 2001, based on
      ! Conley and Inman 1994
      ! In C&I derivation f = 2*cf and w is positive up, and b=0.9-2.0 so modified to:
      ! Tau_w/Tau_0 = C/(2*cf) * -infil/abs(u)
      ! where C is a constant of ~2
      !
      ! Can re-write to taubx_total = taubx + taubx_add
      ! 
      ! Tau ~ cf
      ! cf_w/cf_0 = Tau_w/Tau_0
      !
      ! Therefore: taubx_total = cf_0*rho*u*abs(u) + (1-cf_w/cf_0)*cf*rho*u*abs(u)
      ! (here we neglect the urms part to the drag in taubx_add for surf-beat, as intrawave velocities and intrawave infiltration 
      !  not resolved in surf-beat/stationary type approach)

      ! U-points
      do j=1,s%ny+1
         do i=1,s%nx
            if(s%wetu(i,j)==1) then
               infilb(i,j) = 0.5d0*(s%infil(i,j)+s%infil(i+1,j))
            else
               infilb(i,j) = 0.d0
            endif
         enddo
      enddo
      infilb(s%nx+1,:) = infilb(s%nx,:)
      where(s%wetu==1)
         !!!Ubed = sqrt(cfu)*abs(vmageu) older
         !Ubed = abs(s%uu) ! original XBeach-G
         Ubed = abs(s%vmageu)
         Ventilation  = -infilb/max(Ubed,1.d-6)
         Ventilation  = max(min(Ventilation,epsVentilation),-epsVentilation)
         blphi = 0.9d0/2/s%cfu*Ventilation
         ! need to account for devision by zero
         where(infilb>=0.d0)
            blphi = min(blphi,-1d-4)
            facbl = blphi/(exp(blphi)-1.d0)
            facbl = min(facbl,maxEnhancement1) 
         elsewhere
            blphi = max(blphi,1d-4)
            facbl = blphi/(exp(blphi)-1.d0)
            facbl = max(facbl,maxEnhancement2)
         endwhere
         !s%cfu = s%cfu*facbl
         !s%cfu = min(s%cfu,1.d0)
         ! backward compatibility maximum value (almost never met)
         facbl = min(facbl,1.d0/s%cfu) ! maximizes cfu at 1.0 (=huge)
         s%taubx_add = s%taubx_add + (facbl-1.d0)*s%cfu*par%rho*s%ueu*s%vmageu
      endwhere
      
      
      ! V-points
      if(par%ny>0) then
         do j=1,s%ny
            do i=1,s%nx+1
               if(s%wetv(i,j)==1) then
                  infilb(i,j) = 0.5d0*(s%infil(i,j)+s%infil(i,j+1))
               else
                  infilb(i,j) = 0.d0
               endif
            enddo
         enddo
         infilb(:,s%ny+1) = s%infil(:,s%ny)
      else
         ! v-points centered on cell centre
         infilb = s%infil
      endif
      where(s%wetv==1)
         !!!Ubed = sqrt(cfv)*abs(vmagev)  ! older
         ! Ubed = abs(vv) ! orginal XBeach-G
         Ubed = abs(s%vmagev)
         Ventilation  = -infilb/max(Ubed,1.d-6)
         Ventilation  = max(min(Ventilation,epsVentilation),-epsVentilation)
         blphi = 0.9d0/2/s%cfv*Ventilation
         ! need to account for devision by zero
         where(infilb>=0.d0)
            blphi = min(blphi,-1d-4)
            facbl = blphi/(exp(blphi)-1.d0)
            facbl = min(facbl,maxEnhancement1)
         elsewhere
            blphi = max(blphi,1d-4)
            facbl = blphi/(exp(blphi)-1.d0)
            facbl = max(facbl,maxEnhancement2)
         endwhere
         !s%cfv = s%cfv*facbl
         !s%cfv = min(s%cfv,1.d0)
         facbl = min(facbl,1.d0/s%cfv) 
         s%tauby_add = s%tauby_add + (facbl-1.d0)*s%cfv*par%rho*s%vev*s%vmagev
      endwhere

   end subroutine infiltration_boundary_layer_effect

   subroutine acceleration_boundary_layer_effect(s,par)

      use params
      use paramsconst
      use spaceparams
   
      IMPLICIT NONE

      type(parameters),intent(in)                 :: par
      type(spacepars)                             :: s
   
      real*8                                      :: Tsmooth,factime
      real*8,dimension(:,:),allocatable,save      :: Fi
 
   
      if(.not.allocated(dudtsmooth)) then
         ! arrays available in bed friction module only
         allocate (dudtsmooth(s%nx+1,s%ny+1))
         allocate (dvdtsmooth(s%nx+1,s%ny+1))
         allocate (ueuold    (s%nx+1,s%ny+1))
         allocate (uevold    (s%nx+1,s%ny+1))
         allocate (vevold    (s%nx+1,s%ny+1))
         allocate (veuold    (s%nx+1,s%ny+1))
         allocate (ueuf      (s%nx+1,s%ny+1))
         allocate (uevf      (s%nx+1,s%ny+1))
         allocate (vevf      (s%nx+1,s%ny+1))
         allocate (veuf      (s%nx+1,s%ny+1))
         allocate (Fi        (s%nx+1,s%ny+1))
         ueuold = 0.d0
         uevold = 0.d0
         vevold = 0.d0
         veuold = 0.d0
         ueuf   = 0.d0
         uevf   = 0.d0
         vevf   = 0.d0
         veuf   = 0.d0
         dtold = par%dt
      endif
   
      ! Problem with derivation of DUdt in U and V points
      ! ueu and veu are not known until end of flow timestep, so
      ! both are based on previous time step. The temporal gradient
      ! can only be found by going back one more time step.
      ! We need to keep track of par%dt so the correct time step for
      ! the gradient is used:
      !
      !    i-2            i-1           i (now)          i+1 (future)
      !         dtold(1)       dtold(2)          par%dt   
      !   ueuold          ueuf         unknown
      !
      ! Note that this assumes that the acceleration term does not change significantly
      ! over 2 timesteps and that the location mask wetu can still be used,
      ! even though it is not the wetu mask of 2 timesteps ago.
      Tsmooth = par%Trep/20
      factime = min(dtold(1)/Tsmooth,1.d0)

      where(s%wetu==1)
         ueuf = (1-factime)*ueuf + factime*s%ueu ! note, because taub_add has linear relation with total acceleration (dvmag/dt)
                                                 ! the x-component taubx_add can be computed from only the u-component of the 
                                                 ! acceleration
         dudtsmooth = (ueuf-ueuold)/dtold(1)
         !dudtsmooth = (ueuf-ueuold)/dtold(1)+ududx*wetu_ududx  ! Following Baldock attempt to correct du/dt to dp/dx
         ! cut off rediculous values
         where (dudtsmooth>100*par%g)
            dudtsmooth = 100*par%g
         elsewhere (dudtsmooth<-100*par%g)
            dudtsmooth = -100*par%g
         endwhere
         s%taubx_add = s%taubx_add + par%ci*par%rho*min(s%D50top,s%hu)*dudtsmooth
      elsewhere
         ueuf = 0.d0
         veuf = 0.d0
      endwhere
      
      ! now set the old variable to the i-1 variables
      ueuold = ueuf
      veuold = veuf
       
      ! in some cases we don't need to compute v-acceleration
      ! - 1D non-h has only u components
      ! - 1D surfbeat / stationary has no v-components if boundary set to wall
      if (par%ny==0 .and. ((par%left==LR_WALL .and. par%right==LR_WALL) .or. par%swave==0)) then
         ! do nothing (don't change tauby_add)
      else
         where(s%wetv==1)
            vevf = (1-factime)*vevf + factime*s%vev
            dvdtsmooth = (vevf-vevold)/dtold(1)
            where (dvdtsmooth>100*par%g)
               dvdtsmooth = 100*par%g
            elsewhere (dvdtsmooth<-100*par%g)
               dvdtsmooth = -100*par%g
            endwhere
            s%tauby_add = s%tauby_add + par%ci*par%rho*min(s%D50top,s%hv)*dvdtsmooth 
         elsewhere
            vevf = 0.d0
            uevf = 0.d0
         endwhere
         
         ! now set the old variable to the i-1 variables
         uevold = uevf
         vevold = vevf
      endif
      
      ! update timekeeping
      dtold(1) = dtold(2)
      dtold(2) = par%dt 
   
   end subroutine acceleration_boundary_layer_effect

   subroutine acceleration_boundary_layer_effect_Nielsen(s,par)

      use params
      use spaceparams
   
      IMPLICIT NONE

      type(parameters),intent(in)                 :: par
      type(spacepars)                             :: s
      real*8                                      :: omegap,Tsmooth,factime,iomegap
      real*8,save                                 :: phirad
   
   
      if(.not.allocated(dudtsmooth)) then
         allocate (dudtsmooth(s%nx+1,s%ny+1))
         allocate (dvdtsmooth(s%nx+1,s%ny+1))
         allocate (ueuold    (s%nx+1,s%ny+1))
         !allocate (uevold    (s%nx+1,s%ny+1))
         allocate (vevold    (s%nx+1,s%ny+1))
         !allocate (veuold    (s%nx+1,s%ny+1))
         allocate (ueuf      (s%nx+1,s%ny+1))
         !allocate (uevf      (s%nx+1,s%ny+1))
         allocate (vevf      (s%nx+1,s%ny+1))
         !allocate (veuf      (s%nx+1,s%ny+1))
         ueuold = 0.d0
         !uevold = 0.d0
         vevold = 0.d0
         !veuold = 0.d0
         ueuf   = 0.d0
         !uevf   = 0.d0
         vevf   = 0.d0
         !veuf   = 0.d0
         dtold = par%dt
         phirad  = par%phit/180*par%px
      endif
   
      omegap = 2*par%px/par%Trep
      iomegap = 1/omegap
      ! Problem with derivation of DUdt in U and V points
      ! ueu and veu are not known until end of flow timestep, so
      ! both are based on previous time step. The temporal gradient
      ! can only be found by going back one more time step.
      ! We need to keep track of par%dt so the correct time step for
      ! the gradient is used:
      !
      !    i-2            i-1           i (now)          i+1 (future)
      !         dtold(1)       dtold(2)          par%dt   
      !   ueuold          ueuf         unknown
      !
      ! Note that this assumes that the acceleration term does not change significantly
      ! over 2 timesteps and that the location mask wetu can still be used,
      ! even though it is not the wetu mask of 2 timesteps ago.
      Tsmooth = par%Trep/20
      factime = min(dtold(1)/Tsmooth,1.d0)
      where(s%wetu==1)
         !ueuf = (1-factime)*ueuf + factime*s%ueu
         !veuf = (1-factime)*veuf + factime*s%veu
         !dudtsmooth = sqrt(((ueuf-ueuold)/dtold(1))**2 + &
         !                  ((veuf-veuold)/dtold(1))**2)
         !
         ! Robert: stick to 1D approach, similar to acceleration_boundary_layer_effect subroutine
         ueuf = (1-factime)*ueuf + factime*s%ueu
         dudtsmooth = (ueuf-ueuold)/dtold(1)
      elsewhere
         ueuf = 0.d0
         !veuf = 0.d0
         dudtsmooth = 0.d0
      endwhere
      
      ! in some cases we don't need to compute v-acceleration
      ! - 1D cases without short wave forcing (e.g., non-h) has only u components
      ! - 1D surfbeat / stationary with short wave forcing has no v-components if boundary set to wall
      if (par%ny==0 .and. ((par%left==LR_WALL .and. par%right==LR_WALL) .or. par%swave==0)) then
         dvdtsmooth = 0.d0
      else
         where(s%wetv==1)
            vevf = (1-factime)*vevf + factime*s%vev
            dvdtsmooth = (vevf-vevold)/dtold(1)
            !uevf = (1-factime)*uevf + factime*uev
            !dvdtsmooth = sqrt(((vevf-vevold)/dtold(1))**2 + &
            !                  ((uevf-uevold)/dtold(1))**2)
         elsewhere
            vevf = 0.d0
            !uevf = 0.d0
            dvdtsmooth = 0.d0
         endwhere
      endif
      ! now set the old variable to the i-1 variables
      ueuold = ueuf
      !uevold = uevf
      vevold = vevf
      !veuold = veuf
      dtold(1) = dtold(2)
      dtold(2) = par%dt
      ! New approach:
      ! we want to use Nielsen to compute bed shear stress (taubx_total):
      ! taubx_total = cf*rho*(cos(phi)*ue + 1/w*sin(phi)*dudt)^2*sign(u)
      !
      ! in flow timestep the "regular" taubx is computed
      ! (note, here we again ignore urms component, Nielsen not meant for wave-averaged approach):
      ! taubx = cf*rho*ueu*umageu
      !
      ! therefore:
      ! taubx_add = taubx_total-taubx = cf*rho*((cos(phi)*ue + 1/w*sin(phi)*dudt)^2*sign(u)-ueu*umageu)
      !
      ! there is inconsistency in definitions of ue and umageu (v-component), which is not resolved in Nielsen paper, but if we
      ! rewrite to maintain consistency for stationary case (phi = 0 and dudt = 0) then:
      !
      ! taubx_total = cf*rho* (cos^2(phi)*ueu*vmageu + 2*cos(phi)*ueu*sin(phi)*1/w*dudt + (sin(phi)*1/w*dudt)^2*sign(ueu))
      ! and
      ! taubx_add = cf*rho * ((cos^2(phi)-1)*ueu*vmageu + 2*cos(phi)*ueu*sin(phi)*1/w*dudt + (sin(phi)*1/w*dudt)^2*sign(ueu))
      ! (ignoring urms component in taubx)
      where(s%wetu==1)
         s%taubx_add = s%taubx_add + s%cfu*par%rho*( ((cos(phirad))**2-1)*s%ueu*s%vmageu &
                                                    + 2*cos(phirad)*s%ueu*sin(phirad)*iomegap*dudtsmooth &
                                                    + (sin(phirad)*iomegap*dudtsmooth)**2*sign(1.d0,s%ueu) &
                                                   )
      endwhere
      if (par%ny==0 .and. ((par%left==LR_WALL .and. par%right==LR_WALL) .or. par%swave==0)) then
         ! do nothing, leave tauby_add as is
      else
         where(s%wetv==1)
            s%tauby_add = s%tauby_add + s%cfv*par%rho*( ((cos(phirad))**2-1)*s%vev*s%vmagev &
                                                       + 2*cos(phirad)*s%vev*sin(phirad)*iomegap*dvdtsmooth &
                                                       + (sin(phirad)*iomegap*dvdtsmooth)**2*sign(1.d0,s%vev) &
                                                      )
         endwhere
      endif
      

   end subroutine acceleration_boundary_layer_effect_Nielsen

   subroutine turbulence_boundary_layer_effect(s,par)

      use params
      use spaceparams
   
      IMPLICIT NONE

      type(parameters),intent(in)                 :: par
      type(spacepars)                             :: s
      real*8,dimension(:,:),allocatable,save      :: kbl
      integer                                     :: i,j
   
      ! Following Reniers et al (2004)
      ! ubed_turb = sqrt(ubed^2+gamma*kb) :: gamma ~ 1
      ! 
      ! New way: compute taubx_add and tauby_add
      ! taubx_total = cf*rho*(ubed^2+gamma*kb)
      ! taubx = cf*rho*(ubed^2)
      ! taubx_add = cf*rho*(gamma*kb)
  
      if (.not.allocated(kbl)) then
         allocate (kbl(s%nx+1,s%ny+1))
      endif
    
       ! U-points
      do j=1,s%ny+1
         do i=1,s%nx
            if(s%wetu(i,j)==1) then
               kbl(i,j) = 0.5d0*(s%kb(i,j)+s%kb(i+1,j))
            else
               kbl(i,j) = 0.d0
            endif
         enddo
      enddo
      kbl(s%nx+1,:) = kbl(s%nx,:)
      where(s%wetu==1 .and. kbl>0.d0)
         ! New way:
         s%taubx_add = s%taubx_add + s%cfu*par%rho*par%gamma_turb*kbl*sign(1.d0,s%ue)
      
         ! Old way:
         !cfu = cfu*(1.d0+(kbl/max(uu**2,1.d-6)))
         !cfu = cfu*(1.d0+(2.d0**(1.25d0)*sqrt(kbl))/max(abs(uu),1.d-6)+(sqrt(2.d0)*kbl)/max(uu**2,1.d-6))
      endwhere
      ! V-points
      if(par%ny>0) then
         do j=1,s%ny
            do i=1,s%nx+1
               if(s%wetv(i,j)==1) then
                  kbl(i,j) = 0.5d0*(s%kb(i,j)+s%kb(i,j+1))
               else
                  kbl(i,j) = 0.d0
               endif
            enddo
         enddo
         kbl(:,s%ny+1) = kbl(:,s%ny)
      else
         ! v-points centered on cell centre
         kbl = s%kb
      endif
      where(s%wetv==1 .and. kbl>0.d0)
          ! New way:
         s%tauby_add = s%tauby_add + s%cfv*par%rho*par%gamma_turb*kbl*sign(1.d0,s%ve)
      
         ! Old way:
         !cfv = cfv*(1.d0+(kbl/max(vv**2,1.d-6)))
         !cfv = cfv*(1.d0+(2.d0**(1.25d0)*sqrt(kbl))/max(abs(vv),1.d-6)+(sqrt(2.d0)*kbl)/max(vv**2,1.d-6))
      endwhere
  
   end subroutine turbulence_boundary_layer_effect


end module bedroughness_module
