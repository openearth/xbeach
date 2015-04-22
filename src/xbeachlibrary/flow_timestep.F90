module flow_timestep_module
   implicit none
   save
contains
   subroutine flow(s,par)
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
      use boundaryconditions
      use flow_secondorder_module
      use nonh_module
      use bedroughness_module

      IMPLICIT NONE

      type(spacepars),target                  :: s
      type(parameters)                        :: par

      integer                                 :: i
      integer                                 :: j,j1,jp1
      real*8,dimension(:,:),allocatable,save  :: vv_old                   !Velocity at previous timestep
      real*8,dimension(:,:),allocatable,save  :: uu_old                   !Velocity at previous timestep
      real*8,dimension(:,:),allocatable,save  :: zs_old
      real*8,dimension(:,:),allocatable,save  :: vsu,usu,vsv,usv,veu,uev
      real*8,dimension(:,:),allocatable,save  :: dudx,dvdy
      real*8,dimension(:,:),allocatable,save  :: us,vs
      real*8                                  :: nuh1,nuh2
      real*8                                  :: dudx1,dudx2,dudy1,dudy2
      real*8                                  :: dvdy1,dvdy2,dvdx1,dvdx2  !Jaap
      real*8                                  :: dalfa                    !difference in grid angles
      real*8                                  :: uin,vin                  !s%u resp s%v-velocity corrected for grid angle change
      real*8                                  :: qin                      !specific discharge entering cell
      real*8                                  :: dzsdnavg                 !alongshore water level slope
      real*8,save                             :: fc
      real*8,dimension(:,:),allocatable,save  :: sinthm,costhm

      integer                                 :: imax,jmax,jmin


      if (.not. allocated(vsu) ) then
         allocate (   vsu(s%nx+1,s%ny+1))
         allocate (   usu(s%nx+1,s%ny+1))
         allocate (   vsv(s%nx+1,s%ny+1))
         allocate (   usv(s%nx+1,s%ny+1))
         allocate (   veu(s%nx+1,s%ny+1))
         allocate (   uev(s%nx+1,s%ny+1))
         allocate (  dudx(s%nx+1,s%ny+1))
         allocate (  dvdy(s%nx+1,s%ny+1))
         allocate (    us(s%nx+1,s%ny+1))
         allocate (    vs(s%nx+1,s%ny+1))
         allocate (sinthm(s%nx+1,s%ny+1))
         allocate (costhm(s%nx+1,s%ny+1))

         if (par%secorder == 1) then
            allocate(vv_old(s%nx+1,s%ny+1)); vv_old = s%vv
            allocate(uu_old(s%nx+1,s%ny+1)); uu_old = s%uu
            allocate(zs_old(s%nx+1,s%ny+1)); zs_old = s%zs
         endif

         s%vu      =0.d0
         vsu     =0.d0
         usu     =0.d0
         vsv     =0.d0
         usv     =0.d0
         s%uv      =0.d0
         veu     =0.d0
         uev     =0.d0
         s%ueu     =0.d0
         s%vev     =0.d0
         dudx    =0.d0
         s%ududx   =0.d0
         dvdy    =0.d0
         s%vdvdy   =0.d0
         s%udvdx   =0.d0
         s%vdudy   =0.d0
         s%viscu   =0.d0
         s%viscv   =0.d0
         us      =0.d0
         vs      =0.d0
         s%hum     =0.d0
         s%hvm     =0.d0
         s%u       =0.d0
         s%v       =0.d0
         s%ue      =0.d0
         s%ve      =0.d0
         fc      =2.d0*par%wearth*sin(par%lat)

         call bedroughness_init(s,par) ! note, this is not yet designed for initialisation
         ! on sglobal, so don't call from initialize.F90
      endif

      ! Super fast 1D
      if (s%ny==0) then
         j1 = 1
      else
         j1 = 2
      endif

      ! Add vertical discharges
      call discharge_boundary_v(s,par)

      !
      ! zs=zs*wetz
      ! Water level slopes
      do j=1,s%ny+1
         do i=2,s%nx
            s%dzsdx(i,j)=(s%zs(i+1,j)+s%ph(i+1,j)-s%zs(i,j)-s%ph(i,j))/s%dsu(i,j)
         end do
      end do
      !    do j=2,ny
      do j=1,s%ny ! Dano need to get correct slope on boundary s%y=0
         do i=1,s%nx+1
            s%dzsdy(i,j)=(s%zs(i,j+1)+s%ph(i,j+1)-s%zs(i,j)-s%ph(i,j))/s%dnv(i,j)
         end do
      end do

      ! wwvv in the next lines
      !  hu(i,j) is more or less a function of hh(i+1,j)
      !  In the parallel case, this needs some action because
      !  for processes not at the bottom, the last row of
      !  hu (hu(nx+1,:)) has to be taken from the neighbour below
      !  hu is used later on in this subroutine, so we have to insert
      !  an mpi call.
      !  The same for hum
      ! Of course, no action is necessary if hu(nx+1:) is never used...
      do j=1,s%ny+1
         do i=1,s%nx+1 !Ap
            ! Water depth in u-points do momentum equation: mean
            !!
            !! ARBJ: mean water depth or weighted water depth? How to deal with this in curvi-linear?
            !!
            s%hum(i,j)=max(.5d0*(s%hh(i,j)+s%hh(min(s%nx,i)+1,j)),par%eps)
         end do
      end do
      ! wwvv here the mpi code to communicate a row of hu
      ! we send to the neighbour above and receive from the neighbour
      ! below:
      ! Wetting and drying criterion (only do momentum balance)
      do j=1,s%ny+1
         do i=1,s%nx+1
            if(s%hu(i,j)>par%eps .and. s%hum(i,j)>par%eps) then  ! Jaap and Pieter: If you want to compute correct advection term
               s%wetu(i,j)=1                                  ! then both s%hu and s%hum should be larger than par%eps. It is not
            else                                             ! necessarily true that if s%hu>par%eps also s%hum>par%eps.
               s%wetu(i,j)=0
            end if
         end do
      end do
      ! wwvv about the same for hv, only in the left-right direction
      ! hv(i,j) is more or less a function of hh(i,j+1)
      ! so in the parallel case, hv(:,ny+1) has to be collected
      ! from the right neighbour
      ! the same for hvm
      do j=1,s%ny+1
         do i=1,s%nx+1
            ! Water depth in v-points do momentum equation: mean
            s%hvm(i,j)=max(.5d0*(s%hh(i,j)+s%hh(i,min(s%ny,j)+1)),par%eps)
         end do
      end do
      ! Wetting and drying criterion (only do momentum balance)
      do j=1,s%ny+1
         do i=1,s%nx+1
            if(s%hv(i,j)>par%eps .and. s%hvm(i,j)>par%eps) then
               s%wetv(i,j)=1
            else
               s%wetv(i,j)=0
            end if
         end do
      end do

      ! Jaap Wetting and drying criterion eta points
      do j=1,s%ny+1
         do i=1,s%nx+1
            !A eta point is wet if any of the surrounding velocity points is wet...
            !            wetz(i,j) = min(1,wetu(max(i,2)-1,j)+wetu(i,j)+wetv(i,j)+wetv(i,max(j,2)-1))
            if(s%hh(i,j)>par%eps) then
               s%wetz(i,j)=1
            else
               s%wetz(i,j)=0
            end if
         end do
      end do

      !
      ! Compute velocity gradients for viscosity terms.
      ! Robert: Check whether should be same gradients as advection terms?
      do j=j1,max(s%ny,1)
         do i=2,s%nx+1
            dudx(i,j) = (s%uu(i,j)-s%uu(i-1,j))/s%dsz(i,j)
         enddo
      enddo
      ! wwvv: added: xmpi_istop
      if (xmpi_istop) then
         dudx(1,:) = 0.d0 ! Robert: by defintion of Neumann boundary
      endif
      if (s%ny>2) then
         do j=2,s%ny+1
            do i=1,s%nx+1
               dvdy(i,j) = (s%vv(i,j)-s%vv(i,j-1))/s%dnz(i,j)
            enddo
         enddo
         ! wwvv: added: xmpi_isleft
         if (xmpi_isleft) then
            dvdy(:,1) = 0.d0 ! Robert: by defintion of Neumann boundary
         endif
      else
         dvdy = 0.d0  ! Robert: by definition of 1D model
      endif

      ! wwvv: added shift_ee
      ! R+D: ToDo, check needed
#ifdef USEMPI
      call xmpi_shift_ee(dvdy)
      call xmpi_shift_ee(dudx)
#endif

      ! Update bed roughness coefficient
      call bedroughness_update(s,par)
      !cf = cfu : Robert: cf is not used anymore

      !
      ! X-direction
      !
      do j=j1,max(s%ny,1)
         do i=2,s%nx
            s%ududx(i,j)        = 0.d0
            qin               = .5d0*(s%qx(i,j)+s%qx(i-1,j))
            if (qin>0) then
               dalfa          = s%alfau(i,j)-s%alfau(i-1,j)
               uin            = s%uu(i-1,j)*cos(dalfa) + s%vu(i-1,j)*sin(dalfa)
               if ((s%uu(i,j)-s%uu(i-1,j))>par%eps_sd) then
                  ! Conservation of energy head
                  s%ududx(i,j)     = s%ududx(i,j) + 0.5d0*(s%uu(i-1,j)+s%uu(i,j))*(s%uu(i,j)-uin)*s%dnz(i,j)*s%dsdnui(i,j)
               else
                  ! Conservation of momentum
                  s%ududx(i,j)     = s%ududx(i,j) +        qin/s%hum(i,j)      *(s%uu(i,j)-uin)*s%dnz(i,j)*s%dsdnui(i,j)
               endif
            endif
            qin               = -.5d0*(s%qx(i,j)+s%qx(i+1,j))
            if (qin>0) then
               dalfa          = s%alfau(i,j)-s%alfau(i+1,j)
               uin            = s%uu(i+1,j)*cos(dalfa) + s%vu(i+1,j)*sin(dalfa)
               if ((s%uu(i+1,j)-s%uu(i,j))>par%eps_sd) then
                  ! Conservation of energy head
                  s%ududx(i,j)     = s%ududx(i,j) - 0.5d0*(s%uu(i+1,j)+s%uu(i,j))*(s%uu(i,j)-uin)*s%dnz(i+1,j)*s%dsdnui(i,j)
               else
                  ! Conservation of momentum
                  s%ududx(i,j)     = s%ududx(i,j) +        qin/s%hum(i,j)      *(s%uu(i,j)-uin)*s%dnz(i+1,j)*s%dsdnui(i,j)
               endif
            endif
         end do
      end do
      do j=2,s%ny
         do i=1,s%nx
            s%vdudy(i,j)        = 0.d0
            qin               = .5d0*(s%qy(i,j-1)+s%qy(i+1,j-1))
            if (qin>0) then
               dalfa          = s%alfau(i,j)-s%alfau(i,j-1)
               uin            = s%uu(i,j-1)*cos(dalfa) + s%vu(i,j-1)*sin(dalfa)
               if ((s%vv(i,j)-s%vv(i,j-1))>par%eps_sd) then
                  ! Conservation of energy head
                  s%vdudy(i,j)     = s%vdudy(i,j) + 0.5d0*(s%vv(i,j-1)+s%vv(i+1,j-1))*(s%uu(i,j)-uin)*s%dsc(i,j-1)*s%dsdnui(i,j)
               else
                  ! Conservation of momentum
                  s%vdudy(i,j)     = s%vdudy(i,j) +        qin/s%hum(i,j)          *(s%uu(i,j)-uin)*s%dsc(i,j-1)*s%dsdnui(i,j)
               endif
            endif
            qin               = -.5d0*(s%qy(i,j)+s%qy(i+1,j))
            if (qin>0) then
               dalfa          = s%alfau(i,j)-s%alfau(i,j+1)
               uin            = s%uu(i,j+1)*cos(dalfa) + s%vu(i,j+1)*sin(dalfa)
               if ((s%vv(i,j+1)-s%vv(i,j))>par%eps_sd) then
                  ! Conservation of energy head
                  s%vdudy(i,j)     = s%vdudy(i,j) - 0.5d0*(s%vv(i,j)+s%vv(i+1,j))*(s%uu(i,j)-uin)*s%dsc(i,j)*s%dsdnui(i,j)
               else
                  ! Conservation of momentum
                  s%vdudy(i,j)     = s%vdudy(i,j) +        qin/s%hum(i,j)       *(s%uu(i,j)-uin)*s%dsc(i,j)*s%dsdnui(i,j)
               endif
            endif
         end do
      end do
      !

      ! Jaap: Slightly changes approach; 1) background viscosity is user defined or obtained from Smagorinsky, 2) nuh = max(nuh,roller induced viscosity)
      if (par%smag == 1) then
         !Use smagorinsky subgrid model
         call visc_smagorinsky(s,par)
      else
         s%nuh = par%nuh
      endif
      ! Add viscosity for wave breaking effects
      if (par%swave == 1) then
         do j=j1,max(s%ny,1)
            do i=2,s%nx
               s%nuh(i,j) = max(s%nuh(i,j),par%nuhfac*s%hh(i,j)*(s%DR(i,j)/par%rho)**(1.0d0/3.0d0)) ! Ad: change to max
            end do
         end do
      elseif (par%swave==0 .and. par%nonh==1) then
         select case (par%nhbreaker)
          case (1)
            where (s%breaking/=0)
               s%nuh = par%breakviscfac*s%nuh
            endwhere
          case (2)
            where (s%breaking==1)
               s%nuh = s%nuh + (par%nuh*par%breakvisclen*s%hh)**2*sqrt(dudx**2+dvdy**2)
            endwhere
          case (3)
            ! Ad en Arnold: Battjes 1975 formulation to smoothen front of lf wave bore in the swash
            ! compute (long) wave turbulence due to breaking
            !nuh = max(par%nuh, par%avis*hloc*sqrt(kturb))
          case default
            ! do nothing to increase viscosity
         end select
      endif

      do j=j1,max(s%ny,1)
         do i=2,s%nx
            !write(*,*)i,j,2
            dudx1 = s%nuh(i+1,j)*s%hh(i+1,j)*(s%uu(i+1,j)-s%uu(i,j))/s%dsz(i+1,j)
            dudx2 = s%nuh(i,j)  *s%hh(i  ,j)*(s%uu(i,j)-s%uu(i-1,j))/s%dsz(i,j)
            s%viscu(i,j) = (1.0d0/s%hum(i,j))*( 2*(dudx1-dudx2)/(s%dsz(i,j)+s%dsz(i+1,j)) )
         end do
      end do

      if (par%smag == 1) then
         !
         !   For non constant eddy viscosity the stress terms read:
         !
         !   d           d           d      d         d
         !   -- [ 2* mu -- (U) ]  + --[ mu --(U) +mu --(V) ]
         !   dx         dx          dy     dy        dx
         !
         !                d         d
         !   Only when   -- [mu] = --[mu] = 0 we have (for incompressible flow)
         !               dx        dy
         !
         !         2            2
         !        d            d
         !    mu --2 [U] + mu --2 [U]
         !       dx           dy
         !
         s%viscu = 2.0d0*s%viscu
      endif

      do j=2,s%ny
         do i=2,s%nx
            !Nuh is defined at eta points, interpolate from four surrounding points
            nuh1  = .25d0*(s%nuh(i,j)+s%nuh(i+1,j)+s%nuh(i+1,j+1)+s%nuh(i,j+1))
            nuh2  = .25d0*(s%nuh(i,j)+s%nuh(i+1,j)+s%nuh(i+1,j-1)+s%nuh(i,j-1))

            dudy1 = nuh1 *.5d0*(s%hvm(i,j  )+s%hvm(i+1,j  ))*(s%uu(i,j+1)-s%uu(i,j))/s%dnc(i,j)
            dudy2 = nuh2 *.5d0*(s%hvm(i,j-1)+s%hvm(i+1,j-1))*(s%uu(i,j)-s%uu(i,j-1))/s%dnc(i,j-1)
            s%viscu(i,j) = s%viscu(i,j) + (1.0d0/s%hum(i,j))* &
            ( 2.0d0*(dudy1-dudy2)/(s%dnc(i,j)+s%dnc(i,j-1)) )*s%wetu(i,j+1)*s%wetu(i,j-1)
         end do
      end do

      if (par%smag == 1) then
         do j=2,s%ny
            do i=2,s%nx
               !Nuh is defined at eta points, interpolate from four surrounding points
               nuh1  = .25d0*(s%nuh(i,j)+s%nuh(i+1,j)+s%nuh(i+1,j+1)+s%nuh(i,j+1))
               nuh2  = .25d0*(s%nuh(i,j)+s%nuh(i+1,j)+s%nuh(i+1,j-1)+s%nuh(i,j-1))

               dvdx1 = nuh1*.5d0*(s%hvm(i,j  )+s%hvm(i+1,j  ))*(s%vv(i+1,j  )-s%vv(i,j  ))/s%dsc(i,j)
               dvdx2 = nuh2*.5d0*(s%hvm(i,j-1)+s%hvm(i+1,j-1))*(s%vv(i+1,j-1)-s%vv(i,j-1))/s%dsc(i,j-1)
               s%viscu(i,j) = s%viscu(i,j) + (1.d0/s%hum(i,j))*(dvdx1-dvdx2)/s%dnz(i,j)*real(s%wetv(i+1,j) &
               * s%wetv(i,j)*s%wetv(i+1,j-1)*s%wetv(i,j-1),8)
            enddo
         enddo
      endif !smag ==1 and s%ny>0
      !
      ! Bed friction term
      where (s%wetu==1)
         s%taubx=s%cfu*par%rho*s%ueu*sqrt((1.16d0*s%urms)**2+s%vmageu**2) !Ruessink et al, 2001
      elsewhere
         s%taubx = 0.d0
      endwhere
      !
      ! Explicit Euler step momentum u-direction
      !
      do j=j1,max(s%ny,1)
         ! do i=2,nx-1   ! wwvv uu(nx,:) is never computed in this subroutine, is that ok?
         if (xmpi_isbot) then
            imax = s%nx-1
         else
            imax=s%nx
         endif
         do i=2,imax ! wwvv with this modification, parallel and serial version
            ! give the same results. If this modification is not ok, then
            ! we have a problem
            if(s%wetu(i,j)==1) then
               s%uu(i,j)=s%uu(i,j)-par%dt*(s%ududx(i,j)+s%vdudy(i,j)-s%viscu(i,j) & !Ap,Robert,Jaap
               + par%g*s%dzsdx(i,j) &
               + s%taubx(i,j)/(par%rho*s%hu(i,j)) &  ! Dano: changed s%hum to s%hu NOT s%cf volume approach
               + s%Fvegu(i,j) &
               - par%lwave*s%Fx(i,j)/(par%rho*max(s%hum(i,j),par%hmin)) &
               - fc*s%vu(i,j) &
               - par%rhoa*par%Cd*s%windsu(i,j)*sqrt(s%windsu(i,j)**2+s%windnv(i,j)**2)/(par%rho*s%hum(i,j)))    ! Kees: wind correction
            else
               s%uu(i,j)=0.0d0
            end if
         end do
      end do
      ! Lateral boundary conditions for uu
      if (s%ny>0) then
         if (xmpi_isleft) then !Dano/Robert only on outer boundary
            s%uu(1:s%nx+1,1)=s%uu(1:s%nx+1,2) ! RJ: can also be done after continuity but more appropriate here
         endif
         ! Lateral boundary at y=ny*dy
         if (xmpi_isright) then !Dano/Robert only at outer boundary
            s%uu(1:s%nx+1,s%ny+1)=s%uu(1:s%nx+1,s%ny) ! RJ: can also be done after continuity but more appropriate here
         endif
      endif
#ifdef USEMPI
      call xmpi_shift_ee(s%uu)
#endif
      !
      ! Y-direction
      !
      ! Robert: Complicated
      !         In overall model (not mpi-subdomains) vv(:,ny+1) is never calculated and has dummy values
      !         In overall model vv(:,1) and vv(:,ny) are used as boundary conditions
      !                          vv(:,1) and vv(:,ny) = 0 for wall boundary conditions
      !                          vv(:,1) and vv(:,ny) = copied from internal in Neumann conditions
      !                          vv(:,1) and vv(:,ny) = calculated without advection terms in no_advec conditions
      !                          vv(:,1) and vv(:,ny) = calculated with advection terms in free conditions
      !         In overall model (or in non-MPI) vv is calculated over j=2,ny-1
      !
      !         In mpi subdomain vv(:,1) and vv(:,ny+1) are communicated and no boundary conditions are given.
      !         In mpi subdomain vv(:,ny+1) is NOT a dummy and can be used to calculate advection terms
      !         In mpi subdomain vv is calculated over j=2,ny

      ! Advection term vdvdy
      ! Robert: in MPI every model subdomain has vv(:,ny+1) from other domain
      !         at overall model boundaries vv(:,ny+1) assumed to be vv(:,ny) and vv(:,0) = vv(:,1)
      !         so calculate vdvdy where we can
      if (xmpi_isright .and. s%ny>0) then  ! no such condition needed for _isleft, because s%vdvdy(:,1) not needed in mpi subdomain
         jmax = s%ny-1
      elseif (xmpi_isright .and. s%ny==0) then
         jmax = 1
      else
         jmax = s%ny
      endif
      s%vdvdy=0.d0
      ! calculate true vdvdy up to ny in central domains and up to ny-1 on isright
      do j=2,jmax
         do i=2,s%nx
            qin               = .5d0*(s%qy(i,j)+s%qy(i,j-1))
            if (qin>0) then
               dalfa          = s%alfav(i,j)-s%alfav(i,j-1)
               vin            = s%vv(i,j-1)*cos(dalfa) - s%uv(i,j-1)*sin(dalfa)
               if ((s%vv(i,j)-s%vv(i,j-1))>par%eps_sd) then
                  ! Conservation of energy head
                  s%vdvdy(i,j)     = s%vdvdy(i,j) + 0.5d0*(s%vv(i,j-1)+s%vv(i,j))*(s%vv(i,j)-vin)*s%dsz(i,j)*s%dsdnvi(i,j)
               else
                  ! Conservation of momentum
                  s%vdvdy(i,j)     = s%vdvdy(i,j) +        qin/s%hvm(i,j)      *(s%vv(i,j)-vin)*s%dsz(i,j)*s%dsdnvi(i,j)
               endif
            endif
            qin               = -.5d0*(s%qy(i,j)+s%qy(i,j+1))
            if (qin>0) then
               dalfa          = s%alfav(i,j)-s%alfav(i,j+1)
               vin            = s%vv(i,j+1)*cos(dalfa) - s%uv(i,j+1)*sin(dalfa)
               if ((s%vv(i,j+1)-s%vv(i,j))>par%eps_sd) then
                  ! Conservation of energy head
                  s%vdvdy(i,j)     = s%vdvdy(i,j) - 0.5d0*(s%vv(i,j+1)+s%vv(i,j))*(s%vv(i,j)-vin)*s%dsz(i,j+1)*s%dsdnvi(i,j)
               else
                  ! Conservation of momentum
                  s%vdvdy(i,j)     = s%vdvdy(i,j) +        qin/s%hvm(i,j)      *(s%vv(i,j)-vin)*s%dsz(i,j+1)*s%dsdnvi(i,j)
               endif
            endif
         enddo
      enddo
      if (s%ny>0) then
         ! Global boundary conditions for vdvdy(:,1) and vdvdy(:,ny), global vdvdy(:,ny+1) not needed anywhere
         if (xmpi_isleft) then
            ! (vv(:,1)-vv(:,0))/dy == 0 so only second part of the vdvdy equation:
            do i=2,s%nx
               qin            = -.5d0*(s%qy(i,1)+s%qy(i,2))
               if (qin>0) then
                  dalfa       = s%alfav(i,1)-s%alfav(i,2)
                  vin         = s%vv(i,2)*cos(dalfa) - s%uv(i,2)*sin(dalfa)
                  s%vdvdy(i,1)  = s%vdvdy(i,1) + qin*(s%vv(i,1)-vin)*s%dsz(i,2)/s%hvm(i,1)*s%dsdnvi(i,1)
               endif
            enddo
         endif
         if (xmpi_isright) then
            ! (vv(:,ny+1)-vv(:,ny))/dy == 0 so only first part of the vdvdy equation:
            do i=2,s%nx
               qin            = .5d0*(s%qy(i,s%ny)+s%qy(i,s%ny-1))
               if (qin>0) then
                  dalfa       = s%alfav(i,s%ny)-s%alfav(i,s%ny-1)
                  vin         = s%vv(i,s%ny-1)*cos(dalfa) - s%uv(i,s%ny-1)*sin(dalfa)
                  s%vdvdy(i,s%ny) = s%vdvdy(i,s%ny) + qin*(s%vv(i,s%ny)-vin)*s%dsz(i,s%ny)/s%hvm(i,s%ny)*s%dsdnvi(i,s%ny)
               endif
            enddo
         endif
      endif

      s%udvdx=0.d0
      if (s%ny>0) then
         ! Robert: udvdx not usually needed at j = 1
         do j=1,s%ny !1,s%ny instead of 2,s%ny
            do i=2,s%nx
               qin            = .5d0*(s%qx(i-1,j)+s%qx(i-1,j+1))
               if (qin>0) then
                  dalfa       = s%alfav(i,j)-s%alfav(i-1,j)
                  vin         = s%vv(i-1,j)*cos(dalfa) - s%uv(i-1,j)*sin(dalfa)
                  if ((s%uu(i,j)-s%uu(i-1,j))>par%eps_sd) then
                     ! Conservation of energy head
                     s%udvdx(i,j)     = s%udvdx(i,j) + 0.5d0*(s%uu(i-1,j)+s%uu(i-1,j+1))*(s%vv(i,j)-vin)*s%dnc(i-1,j)*s%dsdnvi(i,j)
                  else
                     ! Conservation of momentum
                     s%udvdx(i,j)     = s%udvdx(i,j) +        qin/s%hvm(i,j)          *(s%vv(i,j)-vin)*s%dnc(i-1,j)*s%dsdnvi(i,j)
                  endif
               endif
               qin            = -.5d0*(s%qx(i,j)+s%qx(i,j+1))
               if (qin>0) then
                  dalfa       = s%alfav(i,j)-s%alfav(i+1,j)
                  vin         = s%vv(i+1,j)*cos(dalfa) - s%uv(i+1,j)*sin(dalfa)
                  if ((s%uu(i+1,j)-s%uu(i,j))>par%eps_sd) then
                     ! Conservation of energy head
                     s%udvdx(i,j)  = s%udvdx(i,j) - 0.5d0*(s%uu(i,j)+s%uu(i,j+1))*(s%vv(i,j)-vin)*s%dnc(i,j)*s%dsdnvi(i,j)
                  else
                     ! Conservation of momentum
                     s%udvdx(i,j)  = s%udvdx(i,j) +        qin/s%hvm(i,j)      *(s%vv(i,j)-vin)*s%dnc(i,j)*s%dsdnvi(i,j)
                  endif
               endif
            end do
         end do
      else
         do i=2,s%nx
            qin               = s%qx(i-1,1)
            if (qin>0) then
               dalfa          = s%alfav(i,1)-s%alfav(i-1,1)
               vin            = s%vv(i-1,1)*cos(dalfa) - s%uv(i-1,1)*sin(dalfa)
               s%udvdx(i,1)     = s%udvdx(i,1) + qin*(s%vv(i,1)-vin)*s%dnc(i-1,1)/s%hvm(i,1)*s%dsdnvi(i,1)
            endif
            qin               = -s%qx(i,1)
            if (qin>0) then
               dalfa          = s%alfav(i,1)-s%alfav(i+1,1)
               vin            = s%vv(i+1,1)*cos(dalfa) - s%uv(i+1,1)*sin(dalfa)
               s%udvdx(i,1)     = s%udvdx(i,1) + qin*(s%vv(i,1)-vin)*s%dnc(i,1)/s%hvm(i,1)*s%dsdnvi(i,1)
            endif
         enddo
      endif
      !
      s%viscv =0.d0
      do j=2,s%ny
         do i=2,s%nx
            dvdy1 = s%nuh(i,j+1)*s%hh(i,j+1)*(s%vv(i,j+1)-s%vv(i,j))/s%dnz(i,j+1)
            dvdy2 = s%nuh(i,j)  *s%hh(i,j  )*(s%vv(i,j)-s%vv(i,j-1))/s%dnz(i,j)
            s%viscv(i,j) = (1.0d0/s%hvm(i,j))* 2*(dvdy1-dvdy2)/(s%dnz(i,j)+s%dnz(i,j+1))*s%wetv(i,j+1)*s%wetv(i,j-1)
         end do
      end do
      ! Robert: global boundary at (:,1) edge
      if (s%ny>0) then
         if (xmpi_isleft) then
            s%viscv(:,1) = s%viscv(:,2)
         endif
         if (xmpi_isright) then
            s%viscv(:,s%ny) = s%viscv(:,s%ny-1)
         endif
      endif
      !
      ! Viscosity
      !
      if (par%smag == 1) then
         s%viscv = 2.0d0*s%viscv
      endif

      s%nuh = par%nuhv*s%nuh !Robert en Ap: increase s%nuh interaction in d2v/dx2
      do j=1,max(s%ny,1)
         jp1 = min(j+1,s%ny+1)
         do i=2,s%nx
            !Nuh is defined at eta points, interpolate from four surrounding points
            nuh1  = .25d0*(s%nuh(i,j)+s%nuh(i+1,j)+s%nuh(i+1,jp1)+s%nuh(i,jp1))
            nuh2  = .25d0*(s%nuh(i,j)+s%nuh(i-1,j)+s%nuh(i-1,jp1)+s%nuh(i,jp1))

            dvdx1 = nuh1*.5d0*(s%hum(i  ,j)+s%hum(i  ,jp1))*(s%vv(i+1,j)-s%vv(i,j))/s%dsc(i,j)
            dvdx2 = nuh2*.5d0*(s%hum(i-1,j)+s%hum(i-1,jp1))*(s%vv(i,j)-s%vv(i-1,j))/s%dsc(i-1,j)
            s%viscv(i,j) = s%viscv(i,j) + (1.0d0/s%hvm(i,j))*( 2*(dvdx1-dvdx2)/(s%dsc(i-1,j)+s%dsc(i,j)) ) &
            *s%wetv(i+1,j)*s%wetv(i-1,j)
         end do
      end do
      !
      if (par%smag == 1) then
         do j=1,max(s%ny,1)
            jp1 = min(j+1,s%ny+1)
            do i=2,s%nx
               !Nuh is defined at eta points, interpolate from four surrounding points
               nuh1  = .25d0*(s%nuh(i,j)+s%nuh(i+1,j)+s%nuh(i+1,jp1)+s%nuh(i,jp1))
               nuh2  = .25d0*(s%nuh(i,j)+s%nuh(i-1,j)+s%nuh(i-1,jp1)+s%nuh(i,jp1))

               dudy1 = nuh1 *.5d0*(s%hum(i  ,j)+s%hum(i  ,jp1))*(s%uu(i,jp1  )-s%uu(i,j  ))/s%dnc(i,j)
               dudy2 = nuh2 *.5d0*(s%hum(i-1,j)+s%hum(i-1,jp1))*(s%uu(i-1,jp1)-s%uu(i-1,j))/s%dnc(i-1,j)

               s%viscv(i,j) = s%viscv(i,j) + (1.d0/s%hvm(i,j))*(dudy1-dudy2)/s%dsz(i,j)  &
               * real(s%wetu(i,jp1)*s%wetu(i,j)*s%wetu(i-1,jp1)*s%wetv(i-1,j),8)
            enddo
         enddo
      endif
      !
      ! Bed friction term
      !
      where (s%wetv==1)
         s%tauby=s%cfv*par%rho*s%vev*sqrt((1.16d0*s%urms)**2+s%vmagev**2) !Ruessink et al, 2001
      elsewhere
         s%tauby = 0.d0
      endwhere
      !
      ! Explicit Euler step momentum v-direction
      !
      if (s%ny==0) then
         jmin = 1
      else
         jmin = 2
         if (s%ny==2) then
            jmax = 2 ! Robert: very special case of s%ny=2 and xmpi_isright would otherwise lead to no calculation of s%vv with Neumann boundaries
         endif
      endif
      !
      do j=jmin,jmax
         do i=2,s%nx !jaap instead of s%nx+1
            if(s%wetv(i,j)==1) then
               ! Robert: ensure taubx always has the same sign as uu (always decelerates)
               ! Dano: I don't agree
               s%vv(i,j)=s%vv(i,j)-par%dt*(s%udvdx(i,j)+s%vdvdy(i,j)-s%viscv(i,j)& !Ap,Robert,Jaap
               + par%g*s%dzsdy(i,j)&
               + s%tauby(i,j)/(par%rho*s%hv(i,j)) &  ! Dano: s%hv instead of s%hvm, NOT s%cf volume approach
               + s%Fvegv(i,j) &
               - par%lwave*s%Fy(i,j)/(par%rho*max(s%hvm(i,j),par%hmin)) &
               + fc*s%uv(i,j) &
               - par%rhoa*par%Cd*s%windnv(i,j)*sqrt(s%windsu(i,j)**2+s%windnv(i,j)**2)/(par%rho*s%hvm(i,j)))    ! Kees: wind correction
            else
               s%vv(i,j)=0.0d0
            end if
         end do
      end do
      ! Communicate vv at internal boundaries
#ifdef USEMPI
      call xmpi_shift_ee(s%vv)
#endif
      ! Robert: Boundary conditions along the global boundaries
      ! Function flow_lat_bc located in boundaryconditions.F90
      ! function call takes care of 1D vs 2D models and boundary condition types
      if (s%ny>0) then
         if (xmpi_isleft) then
            s%vv(:,1)=flow_lat_bc(s,par,par%right,1,2,s%udvdx(:,1),s%vdvdy(:,1),s%viscv(:,1))
         endif
         if (xmpi_isright) then
            s%vv(:,s%ny)=flow_lat_bc(s,par,par%left,s%ny,s%ny-1,s%udvdx(:,s%ny),s%vdvdy(:,s%ny),s%viscv(:,s%ny))
         endif
      endif

      if (par%nonh==1) then
         !Do explicit predictor step with pressure
         call nonh_explicit(s,par)

#ifdef USEMPI
         call xmpi_shift_ee(s%uu)
         call xmpi_shift_ee(s%vv)
#endif
      end if

      if (par%secorder==1) then
         !Call second order correction to the advection
         call flow_secondorder_advUV(s,par,uu_old,vv_old)
#ifdef USEMPI
         call xmpi_shift_ee(s%uu)
         call xmpi_shift_ee(s%vv)
#endif
      end if

      ! Pieter and Jaap: update hu en hv for continuity
      do j=1,s%ny+1
         do i=1,s%nx+1 !Ap
            ! Water depth in u-points do continuity equation: upwind
            if (s%uu(i,j)>par%umin) then
               if (par%oldhu == 1) then
                  s%hu(i,j)=s%hh(i,j)
               else
                  s%hu(i,j)=s%zs(i,j)-max(s%zb(i,j),s%zb(min(s%nx,i)+1,j))
               endif
            elseif (s%uu(i,j)<-par%umin) then
               if (par%oldhu == 1) then
                  s%hu(i,j)=s%hh(min(s%nx,i)+1,j)
               else
                  s%hu(i,j)=s%zs(min(s%nx,i)+1,j)-max(s%zb(i,j),s%zb(min(s%nx,i)+1,j))
               endif
            else
               s%hu(i,j)=max(max(s%zs(i,j),s%zs(min(s%nx,i)+1,j))-max(s%zb(i,j),s%zb(min(s%nx,i)+1,j)),par%eps)
            end if
         end do
      end do

      s%hu = max(s%hu,0.d0)

      do j=1,s%ny+1
         do i=1,s%nx+1
            ! Water depth in v-points do continuity equation: upwind
            if (s%vv(i,j)>par%umin) then
               if (par%oldhu == 1) then
                  s%hv(i,j)=s%hh(i,j)
               else
                  s%hv(i,j)=s%zs(i,j)-max(s%zb(i,j),s%zb(i,min(s%ny,j)+1))
               endif
            elseif (s%vv(i,j)<-par%umin) then
               if (par%oldhu == 1) then
                  s%hv(i,j)=s%hh(i,min(s%ny,j)+1)
               else
                  s%hv(i,j)=s%zs(i,min(s%ny,j)+1)-max(s%zb(i,j),s%zb(i,min(s%ny,j)+1))
               endif
            else
               s%hv(i,j)=max(max(s%zs(i,j),s%zs(i,min(s%ny,j)+1))-max(s%zb(i,j),s%zb(i,min(s%ny,j)+1)),par%eps)
            end if
         end do
      end do
      s%hv = max(s%hv,0.d0)

      if (par%nonh==1) then
         ! do non-hydrostatic pressure compensation to solve short waves
         call nonh_cor(s,par)
         ! note: MPI shift in subroutine nonh_cor
      end if

      ! Flux in u-point
      s%qx=s%uu*s%hu
      ! Flux in v-points
      ! first column of qy is used later, and it is defined in the loop above
      ! no communication  necessary at this point
      s%qy=s%vv*s%hv
      !
      ! Add horizontal discharges
      !
      call discharge_boundary_h(s,par)
      !
      ! Update water level using continuity eq.
      !
      if (xmpi_isright) then
         jmax=s%ny
      else
         jmax=s%ny+1
      endif
      if (xmpi_isbot) then
         imax=s%nx
      else
         imax=s%nx+1
      endif
      if (s%ny>0) then
         do j=2,jmax
            do i=2,imax
               s%dzsdt(i,j) = (-1.d0)*( s%qx(i,j)*s%dnu(i,j)-s%qx(i-1,j)*s%dnu(i-1,j)  &
               + s%qy(i,j)*s%dsv(i,j)-s%qy(i,j-1)*s%dsv(i,j-1) )*s%dsdnzi(i,j) &
               - s%infil(i,j)
            end do
         end do
         s%zs(2:s%nx,2:s%ny) = s%zs(2:s%nx,2:s%ny)+s%dzsdt(2:s%nx,2:s%ny)*par%dt !Jaap s%nx instead of s%nx+1
      else
         j=1
         do i=2,imax
            s%dzsdt(i,j) = (-1.d0)*( s%qx(i,j)*s%dnu(i,j)-s%qx(i-1,j)*s%dnu(i-1,j) )*s%dsdnzi(i,j) &
            - s%infil(i,j)
         end do
         s%zs(2:s%nx,1) = s%zs(2:s%nx,1)+s%dzsdt(2:s%nx,1)*par%dt !Jaap s%nx instead of s%nx+1
      endif !s%ny>0
      ! call discharge_boundary(s,par)
      !

      if (par%secorder == 1) then
         !Second order correction
         call flow_secondorder_con(s,par,zs_old)
      endif

      !
      ! Lateral boundary conditions
      !
      if (s%ny>0) then
         ! RJ: Neumann water levels in case of right = 1 or right = 0
         do i=1,s%nx+1
            ! Jaap multiply with wetz(i,ny+1)*wetz(i,1) here to prevent prsssure grdaient over land
            dzsdnavg=s%wetz(i,s%ny+1)*s%wetz(i,1)*(s%zs0(i,s%ny+1)-s%zs0(i,1))/(s%ndist(i,s%ny+1)-s%ndist(i,1))
            ! Lateral boundary at y=0
            if (xmpi_isleft) then
               s%zs(i,1)=max(s%zs(i,2) - dzsdnavg*s%dnv(i,1),s%zb(i,1))
            endif
            ! Lateral boundary at y=ny+1
            if (xmpi_isright) then
               s%zs(i,s%ny+1)=max(s%zs(i,s%ny) + dzsdnavg*s%dnv(i,s%ny),s%zb(i,s%ny+1))
            endif
         enddo
      endif !s%ny>0


      ! wwvv zs, uu, vv have to be communicated now, because they are used later on
#ifdef USEMPI
      call xmpi_shift_ee(s%zs)
#endif

      if (par%secorder == 1) then
         vv_old = s%vv
         uu_old = s%uu
         zs_old = s%zs
      endif

      ! offshore boundary
      !
      ! U and V in cell centre; do output and sediment stirring
      !
      s%u(2:s%nx,:)=0.5d0*(s%uu(1:s%nx-1,:)+s%uu(2:s%nx,:))
      if(xmpi_istop) then
         s%u(1,:)=s%uu(1,:)
      endif
      if(xmpi_isbot) then
         s%u(s%nx+1,:)=s%u(s%nx,:)
      endif

      if (s%ny>0) then
         s%v(:,2:s%ny)=0.5d0*(s%vv(:,1:s%ny-1)+s%vv(:,2:s%ny))
         if(xmpi_isleft) then
            s%v(:,1)=s%vv(:,1)
         endif
         if(xmpi_isright) then
            s%v(:,s%ny+1)=s%v(:,s%ny)        ! bas: need this for calculation of s%ee in wci routine
         endif
         !Ap
         s%v(s%nx+1,:)=s%v(s%nx,:)
      else ! Dano
         s%v=s%vv
      endif !s%ny>0

      ! Robert + Jaap: compute derivatives of u and v
      !

      sinthm = sin(s%thetamean-s%alfaz)
      costhm = cos(s%thetamean-s%alfaz)

      ! V-velocities at u-points
      if (s%ny>0) then
         s%vu(1:s%nx,2:s%ny)= 0.25d0*(s%vv(1:s%nx,1:s%ny-1)+s%vv(1:s%nx,2:s%ny)+ &
         s%vv(2:s%nx+1,1:s%ny-1)+s%vv(2:s%nx+1,2:s%ny))
         ! how about boundaries?
         if(xmpi_isleft) then
            s%vu(:,1) = s%vu(:,2)
         endif
         if(xmpi_isright) then
            s%vu(:,s%ny+1) = s%vu(:,s%ny)
         endif
      else
         s%vu(1:s%nx,1)= 0.5d0*(s%vv(1:s%nx,1)+s%vv(2:s%nx+1,1))
      endif !s%ny>0
      ! wwvv fill in vu(:1) and vu(:ny+1) for non-left and non-right processes
      !  and vu(nx+1,:)
      s%vu=s%vu*s%wetu
      ! V-stokes velocities at U point
      if (s%ny>0) then
         vsu(1:s%nx,2:s%ny)=0.5d0*(s%ust(1:s%nx,2:s%ny)*sinthm(1:s%nx,2:s%ny)+ &
         s%ust(2:s%nx+1,2:s%ny)*sinthm(2:s%nx+1,2:s%ny))
         if(xmpi_isleft) then
            vsu(:,1)=vsu(:,2)
         endif
         if(xmpi_isright) then
            vsu(:,s%ny+1) = vsu(:,s%ny)
         endif
      else
         vsu(1:s%nx,1)=0.5d0*(s%ust(1:s%nx,1)*sinthm(1:s%nx,1)+ &
         s%ust(2:s%nx+1,1)*sinthm(2:s%nx+1,1))
      endif !s%ny>0
      ! wwvv same for vsu
      vsu = vsu*s%wetu
      ! U-stokes velocities at U point
      if (s%ny>0) then
         usu(1:s%nx,2:s%ny)=0.5d0*(s%ust(1:s%nx,2:s%ny)*costhm(1:s%nx,2:s%ny)+ &
         s%ust(2:s%nx+1,2:s%ny)*costhm(2:s%nx+1,2:s%ny))
         if(xmpi_isleft) then
            usu(:,1)=usu(:,2)
         endif
         if(xmpi_isright) then
            usu(:,s%ny+1)=usu(:,s%ny)
         endif
      else
         usu(1:s%nx,1)=0.5d0*(s%ust(1:s%nx,1)*costhm(1:s%nx,1)+ &
         s%ust(2:s%nx+1,1)*costhm(2:s%nx+1,1))
      endif !s%ny>0
      ! wwvv same for usu
      usu=usu*s%wetu

      ! V-euler velocities at u-point
      veu = s%vu - vsu
      ! U-euler velocties at u-point
      s%ueu = s%uu - usu
      ! Velocity magnitude at u-points
      if (par%sedtrans == 0) then
         s%vmagu=sqrt(s%uu**2+s%vu**2)
      endif
      ! Eulerian velocity magnitude at u-points
      s%vmageu=sqrt(s%ueu**2+veu**2)

      ! U-velocities at v-points
      if (s%ny>0) then
         s%uv(2:s%nx,1:s%ny)= .25d0*(s%uu(1:s%nx-1,1:s%ny)+s%uu(2:s%nx,1:s%ny)+ &
         s%uu(1:s%nx-1,2:s%ny+1)+s%uu(2:s%nx,2:s%ny+1))
         ! boundaries?
         ! wwvv and what about uv(:,1) ?
         if(xmpi_isright) then
            s%uv(:,s%ny+1) = s%uv(:,s%ny)
         endif
      else
         s%uv(2:s%nx,1)= .5d0*(s%uu(1:s%nx-1,1)+s%uu(2:s%nx,1))
      endif !s%ny>0
      ! wwvv fix uv(:,ny+1) for non-right processes
      ! uv(1,:) and uv(nx+1,:) need to be filled in for
      ! non-bot or top processes
      s%uv=s%uv*s%wetv
      ! V-stokes velocities at V point
      if (s%ny>0) then
         vsv(2:s%nx,1:s%ny)=0.5d0*(s%ust(2:s%nx,1:s%ny)*sinthm(2:s%nx,1:s%ny)+&
         s%ust(2:s%nx,2:s%ny+1)*sinthm(2:s%nx,2:s%ny+1))
         if(xmpi_isleft) then
            vsv(:,1) = vsv(:,2)
         endif
         if(xmpi_isright) then
            vsv(:,s%ny+1) = vsv(:,s%ny)
         endif
      else
         vsv(2:s%nx,1)= s%ust(2:s%nx,1)*sinthm(2:s%nx,1)
      endif !s%ny>0
      ! wwvv fix vsv(:,1) and vsv(:,ny+1) and vsv(1,:) and vsv(nx+1,:)

      vsv=vsv*s%wetv
      ! U-stokes velocities at V point
      if (s%ny>0) then
         usv(2:s%nx,1:s%ny)=0.5d0*(s%ust(2:s%nx,1:s%ny)*costhm(2:s%nx,1:s%ny)+&
         s%ust(2:s%nx,2:s%ny+1)*costhm(2:s%nx,2:s%ny+1))
         if(xmpi_isleft) then
            usv(:,1) = usv(:,2)
         endif
         if(xmpi_isright) then
            usv(:,s%ny+1) = usv(:,s%ny)
         endif
      else
         usv(2:s%nx,1)=s%ust(2:s%nx,1)*costhm(2:s%nx,1)
      endif !s%ny>0
      ! wwvv fix usv(:,1) and usv(:,ny+1) and usv(1,:) and usv(nx+1,:)
      usv=usv*s%wetv

      ! V-euler velocities at V-point
      s%vev = s%vv - vsv
      ! U-euler velocties at V-point
      uev = s%uv - usv
      ! Velocity magnitude at v-points
      if (par%sedtrans==0) then
         s%vmagv=sqrt(s%uv**2+s%vv**2)
      endif
      ! Eulerian velocity magnitude at v-points
      s%vmagev=sqrt(uev**2+s%vev**2)


      ! Ue and Ve in cell centre; do output and sediment stirring
      s%ue(2:s%nx,:)=0.5d0*(s%ueu(1:s%nx-1,:)+s%ueu(2:s%nx,:))
      s%ue(1,:)=s%ueu(1,:)
      ! wwvv ue(nx+1,:) ?

      if (s%ny>0) then
         s%ve(:,2:s%ny)=0.5d0*(s%vev(:,1:s%ny-1)+s%vev(:,2:s%ny)) !Jaap s%ny+1
         s%ve(:,1)=s%vev(:,1)
      else
         s%ve(:,1) = s%vev(:,1)
      endif !s%ny>0
      ! wwvv vev(nx+1,:) ?
      !
      s%hold =s%hh    ! wwvv ?  s%hold is never else used
      !
      s%hh=max(s%zs-s%zb,par%eps)

      s%maxzs=max(s%zs,s%maxzs)
      s%minzs=min(s%zs,s%minzs)

   end subroutine flow


   subroutine visc_smagorinsky(s,par)
      use params
      use spaceparams
      use xmpi_module

      IMPLICIT NONE

      ! DATE               AUTHOR               CHANGES
      !
      ! December 2010       Pieter Bart Smit     New Subroutine
      ! March    2010       Pieter Bart Smit     Changed formulation to standard smag. model

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------
      !
      !  Calculates the turbulent viscocity coefficient nuh according to the smagorinsky
      !  subgrid model.
      !
      !--------------------------        METHOD          ----------------------------
      !
      ! The turbulent viscocity is given as:
      !
      ! nuh = C^2*dx*dy*Tau
      !
      ! Tau =2^(1/2) *  [ (du/dx)^2 + (dv/dy)^2 + 1/2 * (du/dy + dv/dx)^2 ] ^ (1/2)
      !
      ! Where
      !
      ! dx,dy : grid size
      ! C     : Constant ~0.15 (set by par%nuh)
      ! Tau   : Measure for the magnitude of the turbulent stresses
      !
      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars),target                   ,intent(inout) :: s
      type(parameters)                         ,intent(in)    :: par

      !--------------------------     LOCAL VARIABLES    ----------------------------

      real*8                                                  :: dudx !U Velocity gradient in x-dir
      real*8                                                  :: dudy !U Velocity gradient in y-dir
      real*8                                                  :: dvdx !V Velocity gradient in x-dir
      real*8                                                  :: dvdy !V Velocity gradient in y-dir
      real*8                                                  :: Tau  !Measure for magnitude viscous stresses
      real*8                                                  :: l    !Local gridcell area
      integer                                                 :: i    !Index variable
      integer                                                 :: j    !Index variable


      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------

      !MPI WARNING -> Check loop indices
      if (s%ny>2) then
         do j=2,s%ny
            do i=2,s%nx
               dudx = (s%uu(i,j)-s%uu(i-1,j))/s%dsz(i,j)
               dudy = .5d0*(s%uu(i,j+1) - s%uu(i,j-1) + s%uu(i-1,j+1) - s%uu(i-1,j-1))/(s%dnv(i,j)+s%dnv(i,j-1))
               dvdx = .5d0*(s%vv(i+1,j) - s%vv(i-1,j) + s%vv(i+1,j-1) - s%vv(i-1,j-1))/(s%dsu(i,j)+s%dsu(i-1,j))
               dvdy = (s%vv(i,j)-s%vv(i,j-1))/s%dnz(i,j)
               Tau  = sqrt(2.0d0 * dudx**2+2.0d0 * dvdy**2 + (dvdx+dudy)**2)
               l    = 1.d0/s%dsdnzi(i,j)
               s%nuh(i,j) = par%nuh**2 * l * Tau * real(s%wetu(i,j)*s%wetu(i-1,j)*s%wetv(i,j)*s%wetv(i,j-1),kind=8)
            enddo
         enddo

      else

         j = max(s%ny,1)

         do i=2,s%nx
            dudx = (s%uu(i,j)-s%uu(i-1,j))/s%dsz(i,j)
            dvdx = (s%vv(i+1,j) - s%vv(i-1,j) )/(s%dsu(i,j)+s%dsu(i-1,j))
            Tau  = sqrt(2.0d0 * dudx**2 + dvdx**2)

            if (par%dy > -1.d0) then
               l = s%dsz(i,j)*par%dy
            else
               l = s%dsz(i,j)**2
            endif

            s%nuh(i,j) = par%nuh**2 * l * Tau * real(s%wetu(i,j)*s%wetu(i-1,j),kind=8)
         enddo

      endif !s%ny>2

      if (s%ny>0) then
         if (xmpi_isleft)    s%nuh(:,1)      = s%nuh(:,2)
         if (xmpi_isright)   s%nuh(:,s%ny+1)   = s%nuh(:,s%ny)
      endif

      if (xmpi_istop)       s%nuh(1,:)      = s%nuh(2,:)
      if (xmpi_isbot)       s%nuh(s%nx+1,:)   = s%nuh(s%nx,:)

   end subroutine visc_smagorinsky

end module flow_timestep_module
