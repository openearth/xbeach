module flow_timestep_module
contains
  subroutine flow_timestep(s,par)
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

#ifndef USEMPI
    use flow_secondorder_module
    use nonh_module
#endif


    IMPLICIT NONE

    type(spacepars),target                  :: s
    type(parameters)                        :: par

    integer                                 :: i
    integer                                 :: j,j1,jp1
    real*8,dimension(:,:),allocatable,save  :: vv_old                   !Velocity at previous timestep
    real*8,dimension(:,:),allocatable,save  :: uu_old                   !Velocity at previous timestep
    real*8,dimension(:,:),allocatable,save  :: zs_old
    real*8,dimension(:,:),allocatable,save  :: vsu,usu,vsv,usv,veu,uev
    real*8,dimension(:,:),allocatable,save  :: ududx,vdvdy,udvdx,vdudy
    real*8,dimension(:,:),allocatable,save  :: viscu,viscv
    real*8,dimension(:,:),allocatable,save  :: us,vs
    real*8                                  :: nuh1,nuh2
    real*8                                  :: dudx1,dudx2,dudy1,dudy2
    real*8                                  :: dvdy1,dvdy2,dvdx1,dvdx2  !Jaap
	real*8                                  :: dalfa                    !difference in grid angles
	real*8                                  :: uin,vin                  !u resp v-velocity corrected for grid angle change
	real*8                                  :: qin                      !specific discharge entering cell
	real*8                                  :: dzsdnavg                 !alongshore water level slope
	real*8,save                             :: fc

    integer                                 :: imax,jmax,jmin

    include 's.ind'
    include 's.inp'

    if (.not. allocated(vsu) ) then
       allocate (   vsu(nx+1,ny+1))
       allocate (   usu(nx+1,ny+1))
       allocate (   vsv(nx+1,ny+1))
       allocate (   usv(nx+1,ny+1))
       allocate (   veu(nx+1,ny+1))
       allocate (   uev(nx+1,ny+1))
       allocate ( ududx(nx+1,ny+1))
       allocate ( vdvdy(nx+1,ny+1))
       allocate ( udvdx(nx+1,ny+1))
       allocate ( vdudy(nx+1,ny+1))
       allocate ( viscu(nx+1,ny+1))
       allocate ( viscv(nx+1,ny+1))
       allocate (    us(nx+1,ny+1))
       allocate (    vs(nx+1,ny+1))

       if (par%secorder == 1) then
          allocate(vv_old(nx+1,ny+1)); vv_old = vv
          allocate(uu_old(nx+1,ny+1)); uu_old = uu
          allocate(zs_old(nx+1,ny+1)); zs_old = zs
       endif

       vu      =0.d0
       vsu     =0.d0
       usu     =0.d0
       vsv     =0.d0
       usv     =0.d0
       uv      =0.d0
       veu     =0.d0
       uev     =0.d0
       ueu     =0.d0
       vev     =0.d0
       ududx   =0.d0
       vdvdy   =0.d0
       udvdx   =0.d0
       vdudy   =0.d0
       viscu   =0.d0
       viscv   =0.d0
       us      =0.d0
       vs      =0.d0
       hum     =0.d0
       hvm     =0.d0
       u       =0.d0
       v       =0.d0
       ue      =0.d0
       ve      =0.d0
       fc      =2.d0*par%wearth*sin(par%lat)
    endif
    ! update bedfriction coefficient cf
    if (trim(par%bedfriction)=='white-colebrook') then
       cf = par%g/(18.d0*log10(4.*hh/min(hh,D90top)))**2 ! cf = g/C^2 where C = 18*log(4*hh/D90)
    endif

    ! Super fast 1D
    if (ny==0) then
      j1 = 1
    else
      j1 = 2
    endif
    
    ! Add vertical discharges
    call discharge_boundary_v(s,par)
    
    !
    ! zs=zs*wetz
    ! Water level slopes
    do j=1,ny+1
       do i=2,nx
          dzsdx(i,j)=(zs(i+1,j)-zs(i,j))/dsu(i,j)
       end do
    end do
    !    do j=2,ny
    do j=1,ny ! Dano need to get correct slope on boundary y=0
       do i=1,nx+1
          dzsdy(i,j)=(zs(i,j+1)-zs(i,j))/dnv(i,j)
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
    do j=1,ny+1
       do i=1,nx+1 !Ap
          ! Water depth in u-points do momentum equation: mean
		  !!
		  !! ARBJ: mean water depth or weighted water depth? How to deal with this in curvi-linear?
		  !!
          hum(i,j)=max(.5d0*(hh(i,j)+hh(min(nx,i)+1,j)),par%eps) 
       end do
    end do
    ! wwvv here the mpi code to communicate a row of hu
    ! we send to the neighbour above and receive from the neighbour
    ! below:
#ifdef USEMPI
    call xmpi_shift(hum,'m:')
#endif
    ! Wetting and drying criterion (only do momentum balance)
    do j=1,ny+1
       do i=1,nx+1
          if(hu(i,j)>par%eps .and. hum(i,j)>par%eps) then  ! Jaap and Pieter: If you want to compute correct advection term
             wetu(i,j)=1                                  ! then both hu and hum should be larger than par%eps. It is not
          else                                             ! necessarily true that if hu>par%eps also hum>par%eps.
             wetu(i,j)=0
          end if
       end do
    end do
    ! wwvv about the same for hv, only in the left-right direction
    ! hv(i,j) is more or less a function of hh(i,j+1)
    ! so in the parallel case, hv(:,ny+1) has to be collected
    ! from the right neighbour
    ! the same for hvm
    do j=1,ny+1
       do i=1,nx+1
          ! Water depth in v-points do momentum equation: mean
          hvm(i,j)=max(.5d0*(hh(i,j)+hh(i,min(ny,j)+1)),par%eps)
       end do
    end do
    ! send to the left, read from the right
#ifdef USEMPI
    call xmpi_shift(hvm,':n')
#endif
    ! Wetting and drying criterion (only do momentum balance)
    do j=1,ny+1
       do i=1,nx+1
          if(hv(i,j)>par%eps .and. hvm(i,j)>par%eps) then
             wetv(i,j)=1
          else
             wetv(i,j)=0
          end if
       end do
    end do

    ! Jaap Wetting and drying criterion eta points
    do j=1,ny+1
       do i=1,nx+1
          !A eta point is wet if any of the surrounding velocity points is wet...
          !            wetz(i,j) = min(1,wetu(max(i,2)-1,j)+wetu(i,j)+wetv(i,j)+wetv(i,max(j,2)-1))
          if(hh(i,j)>par%eps) then
             wetz(i,j)=1
          else
             wetz(i,j)=0
          end if
       end do
    end do

    !
    ! X-direction
    !
    do j=j1,max(ny,1)
       do i=2,nx
          ududx(i,j)=0.d0
          qin=.5d0*(qx(i,j)+qx(i-1,j))
          if (qin>0) then
             dalfa=alfau(i,j)-alfau(i-1,j)
             uin=uu(i-1,j)*cos(dalfa)+vu(i-1,j)*sin(dalfa)
             ududx(i,j)=ududx(i,j)+qin*(uu(i,j  )-uin)*dnz(i,j)/ hum(i,j)*dsdnui(i,j)
          endif
          qin=-.5d0*(qx(i,j)+qx(i+1,j))
          if (qin>0) then
             dalfa=alfau(i,j)-alfau(i+1,j)
             uin=uu(i+1,j)*cos(dalfa)+vu(i+1,j)*sin(dalfa)
             ududx(i,j)=ududx(i,j)+qin*(uu(i,j  )-uin)*dnz(i+1,j)/ hum(i,j)*dsdnui(i,j)
          endif
       end do
    end do
    ! wwvv fix border rows and columns of ududx
#ifdef USEMPI
    call xmpi_shift(ududx,'1:')
    call xmpi_shift(ududx,'m:')
    call xmpi_shift(ududx,':1')
    call xmpi_shift(ududx,':n')
#endif
    do j=2,ny
       do i=1,nx
          vdudy(i,j)=0.d0
          qin=.5d0*(qy(i,j-1)+qy(i+1,j-1))
          if (qin>0) then
             dalfa=alfau(i,j)-alfau(i,j-1)
             uin=uu(i,j-1)*cos(dalfa)+vu(i,j-1)*sin(dalfa)
             vdudy(i,j)=vdudy(i,j)+qin*(uu(i,j  )-uin)*dsc(i,j-1)/ hum(i,j)*dsdnui(i,j)
          endif
          qin=-.5d0*(qy(i,j)+qy(i+1,j))
          if (qin>0) then
             dalfa=alfau(i,j)-alfau(i,j+1)
             uin=uu(i,j+1)*cos(dalfa)+vu(i,j+1)*sin(dalfa)
             vdudy(i,j)=vdudy(i,j)+qin*(uu(i,j  )-uin)*dsc(i,j)/ hum(i,j)*dsdnui(i,j)
          endif
       end do
    end do
    ! wwvv fix border rows and columns of vdudy
#ifdef USEMPI
    ! call xmpi_shift(vdudy,'1:') ! not necessary wwvv
    call xmpi_shift(vdudy,'m:')
    call xmpi_shift(vdudy,':1')
    call xmpi_shift(vdudy,':n')
#endif
    !

    ! Jaap: Slightly changes approach; 1) background viscosity is user defined or obtained from Smagorinsky, 2) nuh = max(nuh,roller induced viscosity)
    if (par%smag == 1) then
       !Use smagorinsky subgrid model
       call visc_smagorinsky(s,par)
    else
       nuh = par%nuh
    endif

    do j=j1,max(ny,1)
       do i=2,nx
          nuh(i,j) = max(nuh(i,j),par%nuhfac*hh(i,j)*(DR(i,j)/par%rho)**(1.0d0/3.0d0)) ! Ad: change to max
       end do
    end do

    do j=j1,max(ny,1)
       do i=2,nx
       !write(*,*)i,j,2
          dudx1 = nuh(i+1,j)*hh(i+1,j)*(uu(i+1,j)-uu(i,j))/dsz(i+1,j)         
          dudx2 = nuh(i,j)  *hh(i  ,j)*(uu(i,j)-uu(i-1,j))/dsz(i,j)           
          viscu(i,j) = (1.0d0/hum(i,j))*( 2*(dudx1-dudx2)/(dsz(i,j)+dsz(i+1,j)) )
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
       viscu = 2.0d0*viscu
    endif

    do j=2,ny
       do i=2,nx
          !Nuh is defined at eta points, interpolate from four surrounding points
          nuh1  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,j+1)+nuh(i,j+1))
          nuh2  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,j-1)+nuh(i,j-1))

          dudy1 = nuh1 *.5d0*(hvm(i,j  )+hvm(i+1,j  ))*(uu(i,j+1)-uu(i,j))/dnc(i,j)
          dudy2 = nuh2 *.5d0*(hvm(i,j-1)+hvm(i+1,j-1))*(uu(i,j)-uu(i,j-1))/dnc(i,j-1)
          viscu(i,j) = viscu(i,j) + (1.0d0/hum(i,j))*( 2.0d0*(dudy1-dudy2)/(dnc(i,j)+dnc(i,j-1)) )*wetu(i,j+1)*wetu(i,j-1)
       end do
    end do

    if (par%smag == 1) then
       do j=2,ny
          do i=2,nx
             !Nuh is defined at eta points, interpolate from four surrounding points
             nuh1  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,j+1)+nuh(i,j+1))
             nuh2  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,j-1)+nuh(i,j-1))

             dvdx1 = nuh1*.5d0*(hvm(i,j  )+hvm(i+1,j  ))*(vv(i+1,j  )-vv(i,j  ))/dsc(i,j)
             dvdx2 = nuh2*.5d0*(hvm(i,j-1)+hvm(i+1,j-1))*(vv(i+1,j-1)-vv(i,j-1))/dsc(i,j-1)
             viscu(i,j) = viscu(i,j) + (1.d0/hum(i,j))*(dvdx1-dvdx2)/dnz(i,j)*real(wetv(i+1,j) &
                  * wetv(i,j)*wetv(i+1,j-1)*wetv(i,j-1),8)
          enddo
       enddo
    endif !smag ==1 and ny>0
    !
    ! Bed friction term
	where (wetu==1)
	   taubx=cf*par%rho*ueu*sqrt((1.16d0*urms)**2+vmageu**2) !Ruessink et al, 2001
	elsewhere
	   taubx = 0.d0
	endwhere
    !
    ! Explicit Euler step momentum u-direction
    !
    do j=j1,max(ny,1)
       ! do i=2,nx-1   ! wwvv uu(nx,:) is never computed in this subroutine, is that ok?
       if (xmpi_isbot) then
          imax = nx-1
       else
          imax=nx
       endif
       do i=2,imax ! wwvv with this modification, parallel and serial version
          ! give the same results. If this modification is not ok, then
          ! we have a problem
          if(wetu(i,j)==1) then
             uu(i,j)=uu(i,j)-par%dt*(ududx(i,j)+vdudy(i,j)-viscu(i,j) & !Ap,Robert,Jaap
                  + par%g*dzsdx(i,j) &
                  + taubx(i,j)/(par%rho*hu(i,j)) &  ! Dano: changed hum to hu NOT cf volume approach
                  - par%lwave*Fx(i,j)/(par%rho*max(hum(i,j),par%hmin)) &
                  - fc*vu(i,j) &
                  - par%rhoa*par%Cd*windsu(i,j)**2/(par%rho*hum(i,j)))
          else
             uu(i,j)=0.0d0
          end if
       end do
    end do
   	! Lateral boundary conditions for uu
   	if (ny>0) then
	   if (xmpi_isleft) then !Dano/Robert only on outer boundary
          uu(1:nx+1,1)=uu(1:nx+1,2) ! RJ: can also be done after continuity but more appropriate here
       endif
	   ! Lateral boundary at y=ny*dy 
	   if (xmpi_isright) then !Dano/Robert only at outer boundary
          uu(1:nx+1,ny+1)=uu(1:nx+1,ny) ! RJ: can also be done after continuity but more appropriate here
       endif
    endif
#ifdef USEMPI
    ! wwvv qx is used later on, also the first row, fix the first row
    ! of uu first
    call xmpi_shift(uu,'1:')
    call xmpi_shift(uu,'m:')
    call xmpi_shift(uu,':1')
    call xmpi_shift(uu,':n')
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
    if (xmpi_isright .and. ny>0) then  ! no such condition needed for _isleft, because vdvdy(:,1) not needed in mpi subdomain
       jmax = ny-1
    elseif (xmpi_isright .and. ny==0) then
	   jmax = 1
	else
       jmax = ny
    endif
	vdvdy=0.d0
	! calculate true vdvdy up to ny in central domains and up to ny-1 on isright
	do j=2,jmax
       do i=2,nx
          qin=.5d0*(qy(i,j)+qy(i,j-1))
          if (qin>0) then
             dalfa=alfav(i,j)-alfav(i,j-1)
             vin=vv(i,j-1)*cos(dalfa)-uv(i,j-1)*sin(dalfa)
             vdvdy(i,j)=vdvdy(i,j)+qin*(vv(i,j  )-vin)*dsz(i,j)/ hvm(i,j)*dsdnvi(i,j)
          endif
          qin=-.5d0*(qy(i,j)+qy(i,j+1))
          if (qin>0) then
             dalfa=alfav(i,j)-alfav(i,j+1)
             vin=vv(i,j+1)*cos(dalfa)-uv(i,j+1)*sin(dalfa)
             vdvdy(i,j)=vdvdy(i,j)+qin*(vv(i,j  )-vin)*dsz(i,j+1)/ hvm(i,j)*dsdnvi(i,j)
          endif
       enddo
    enddo
	if (ny>0) then 
       ! Global boundary conditions for vdvdy(:,1) and vdvdy(:,ny), global vdvdy(:,ny+1) not needed anywhere
       if (xmpi_isleft) then
          ! (vv(:,1)-vv(:,0))/dy == 0 so only second part of the vdvdy equation:
          do i=2,nx
		     qin=-.5d0*(qy(i,1)+qy(i,2))
             if (qin>0) then
                dalfa=alfav(i,1)-alfav(i,2)
                vin=vv(i,2)*cos(dalfa)-uv(i,2)*sin(dalfa)
                vdvdy(i,1)=vdvdy(i,1)+qin*(vv(i,1  )-vin)*dsz(i,2)/ hvm(i,1)*dsdnvi(i,1)
             endif
		  enddo
       endif
	   if (xmpi_isright) then
	      ! (vv(:,ny+1)-vv(:,ny))/dy == 0 so only first part of the vdvdy equation:
		  do i=2,nx
		     qin=.5d0*(qy(i,ny)+qy(i,ny-1))
		     if (qin>0) then
                dalfa=alfav(i,ny)-alfav(i,ny-1)
                vin=vv(i,ny-1)*cos(dalfa)-uv(i,ny-1)*sin(dalfa)
                vdvdy(i,ny)=vdvdy(i,ny)+qin*(vv(i,ny  )-vin)*dsz(i,ny)/ hvm(i,ny)*dsdnvi(i,ny)
             endif
          enddo
       endif
    endif

    udvdx=0.d0
    if (ny>0) then
	   ! Robert: udvdx not usually needed at j = 1
	   do j=1,ny !1,ny instead of 2,ny
          do i=2,nx
             qin=.5d0*(qx(i-1,j)+qx(i-1,j+1))
             if (qin>0) then
                dalfa=alfav(i,j)-alfav(i-1,j)
                vin=vv(i-1,j)*cos(dalfa)-uv(i-1,j)*sin(dalfa)
                udvdx(i,j)=udvdx(i,j)+qin*(vv(i,j  )-vin)*dnc(i-1,j)/ hvm(i,j)*dsdnvi(i,j)
             endif
             qin=-.5d0*(qx(i,j)+qx(i,j+1))
             if (qin>0) then
                dalfa=alfav(i,j)-alfav(i+1,j)
                vin=vv(i+1,j)*cos(dalfa)-uv(i+1,j)*sin(dalfa)
                udvdx(i,j)=udvdx(i,j)+qin*(vv(i,j  )-vin)*dnc(i,j)/ hvm(i,j)*dsdnvi(i,j)
             endif
          end do
       end do
    else
	   do i=2,nx
	      qin=qx(i-1,1)
          if (qin>0) then
             dalfa=alfav(i,1)-alfav(i-1,1)
             vin=vv(i-1,1)*cos(dalfa)-uv(i-1,1)*sin(dalfa)
             udvdx(i,1)=udvdx(i,1)+qin*(vv(i,1  )-vin)*dnc(i-1,1)/ hvm(i,1)*dsdnvi(i,1)
          endif
          qin=-qx(i,1)
          if (qin>0) then
             dalfa=alfav(i,1)-alfav(i+1,1)
             vin=vv(i+1,1)*cos(dalfa)-uv(i+1,1)*sin(dalfa)
             udvdx(i,1)=udvdx(i,1)+qin*(vv(i,1  )-vin)*dnc(i,1)/ hvm(i,1)*dsdnvi(i,1)
          endif
       enddo
    endif
	!
	viscv =0.d0	
    do j=2,ny
       do i=2,nx
          dvdy1 = nuh(i,j+1)*hh(i,j+1)*(vv(i,j+1)-vv(i,j))/dnz(i,j+1)
          dvdy2 = nuh(i,j)  *hh(i,j  )*(vv(i,j)-vv(i,j-1))/dnz(i,j)
          viscv(i,j) = (1.0d0/hvm(i,j))* 2*(dvdy1-dvdy2)/(dnz(i,j)+dnz(i,j+1))*wetv(i,j+1)*wetv(i,j-1)
       end do
    end do
	! Robert: global boundary at (:,1) edge
	if (ny>0) then
       if (xmpi_isleft) then
          viscv(i,1) = viscv(i,2)
       endif
       if (xmpi_isright) then
          viscv(i,ny) = viscv(i,ny-1)
       endif
    endif
	!
	! Viscosity
	!
    if (par%smag == 1) then
       viscv = 2.0d0*viscv
    endif

    nuh = par%nuhv*nuh !Robert en Ap: increase nuh interaction in d2v/dx2
    do j=1,max(ny,1)
       jp1 = min(j+1,ny+1)
       do i=2,nx
          !Nuh is defined at eta points, interpolate from four surrounding points
          nuh1  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,jp1)+nuh(i,jp1))
          nuh2  = .25d0*(nuh(i,j)+nuh(i-1,j)+nuh(i-1,jp1)+nuh(i,jp1))

          dvdx1 = nuh1*.5d0*(hum(i  ,j)+hum(i  ,jp1))*(vv(i+1,j)-vv(i,j))/dsc(i,j)
          dvdx2 = nuh2*.5d0*(hum(i-1,j)+hum(i-1,jp1))*(vv(i,j)-vv(i-1,j))/dsc(i-1,j)
          viscv(i,j) = viscv(i,j) + (1.0d0/hvm(i,j))*( 2*(dvdx1-dvdx2)/(dsc(i-1,j)+dsc(i,j)) )*wetv(i+1,j)*wetv(i-1,j)
       end do
    end do
    !
    if (par%smag == 1) then
       do j=1,max(ny,1)
          jp1 = min(j+1,ny+1)
          do i=2,nx
             !Nuh is defined at eta points, interpolate from four surrounding points
             nuh1  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,jp1)+nuh(i,jp1))
             nuh2  = .25d0*(nuh(i,j)+nuh(i-1,j)+nuh(i-1,jp1)+nuh(i,jp1))

             dudy1 = nuh1 *.5d0*(hum(i  ,j)+hum(i  ,jp1))*(uu(i,jp1  )-uu(i,j  ))/dnc(i,j)
             dudy2 = nuh2 *.5d0*(hum(i-1,j)+hum(i-1,jp1))*(uu(i-1,jp1)-uu(i-1,j))/dnc(i-1,j)
             
             viscv(i,j) = viscv(i,j) + (1.d0/hvm(i,j))*(dudy1-dudy2)/dsz(i,j)  &
                  * real(wetu(i,jp1)*wetu(i,j)*wetu(i-1,jp1)*wetv(i-1,j),8)
          enddo
       enddo
    endif
	!
    ! Bed friction term
    !
	where (wetv==1)
	   tauby=cf*par%rho*vev*sqrt((1.16d0*urms)**2+vmagev**2) !Ruessink et al, 2001
	elsewhere
	   tauby = 0.d0
	endwhere
    !
    ! Explicit Euler step momentum v-direction
    !
	if (ny==0) then
	   jmin = 1
	else
	   jmin = 2
	   if (ny==2) then
	      jmax = 2 ! Robert: very special case of ny=2 and xmpi_isright would otherwise lead to no calculation of vv with Neumann boundaries
	   endif
	endif
    !
    do j=jmin,jmax
       do i=2,nx !jaap instead of nx+1
          if(wetv(i,j)==1) then
             ! Robert: ensure taubx always has the same sign as uu (always decelerates)
             ! Dano: I don't agree 
             vv(i,j)=vv(i,j)-par%dt*(udvdx(i,j)+vdvdy(i,j)-viscv(i,j)& !Ap,Robert,Jaap
                  + par%g*dzsdy(i,j)&
                  + tauby(i,j)/(par%rho*hv(i,j)) &  ! Dano: hv instead of hvm, NOT cf volume approach
                  - par%lwave*Fy(i,j)/(par%rho*max(hvm(i,j),par%hmin)) &
                  + fc*uv(i,j) &
                  - par%rhoa*par%Cd*windnv(i,j)**2/(par%rho*hvm(i,j)))
          else
             vv(i,j)=0.0d0
          end if
       end do
    end do
    ! Communicate vv at internal boundaries
#ifdef USEMPI
    call xmpi_shift(vv,'1:')
    call xmpi_shift(vv,'m:')
    call xmpi_shift(vv,':1')
    call xmpi_shift(vv,':n')
#endif
    ! Robert: Boundary conditions along the global boundaries
	! Function flow_lat_bc located in boundaryconditions.F90
	! function call takes care of 1D vs 2D models and boundary condition types
	if (ny>0) then
       if (xmpi_isleft) then
	      vv(:,1)=flow_lat_bc(s,par,par%right,1,2,udvdx(:,1),vdvdy(:,1),viscv(:,1))
	   endif
       if (xmpi_isright) then
	      vv(:,ny)=flow_lat_bc(s,par,par%left,ny,ny-1,udvdx(:,ny),vdvdy(:,ny),viscv(:,ny))
	   endif
	endif

#ifndef USEMPI
    if (par%nonh==1) then
       !Do explicit predictor step with pressure
       call nonh_explicit(s,par,nuh)
    end if

    if (par%secorder==1) then
       !Call second order correction to the advection
       call flow_secondorder_advUV(s,par,uu_old,vv_old)
    end if
#endif

! Robert: include again when 2nd order returned
#ifdef USEMPI
    call xmpi_shift(uu,':1')
    call xmpi_shift(uu,':n')
    call xmpi_shift(uu,'1:')
    call xmpi_shift(uu,'m:')
    call xmpi_shift(vv,':1')
    call xmpi_shift(vv,':n')
    call xmpi_shift(vv,'1:')
    call xmpi_shift(vv,'m:')
#endif

    ! Pieter and Jaap: update hu en hv for continuity
    do j=1,ny+1
       do i=1,nx+1 !Ap
          ! Water depth in u-points do continuity equation: upwind
          if (uu(i,j)>par%umin) then
             if (par%oldhu == 1) then 
			    hu(i,j)=hh(i,j)
			 else
                hu(i,j)=zs(i,j)-max(zb(i,j),zb(min(nx,i)+1,j))
			 endif
          elseif (uu(i,j)<-par%umin) then
		     if (par%oldhu == 1) then 
                hu(i,j)=hh(min(nx,i)+1,j)
			 else
                hu(i,j)=zs(min(nx,i)+1,j)-max(zb(i,j),zb(min(nx,i)+1,j))
		     endif
          else
             hu(i,j)=max(max(zs(i,j),zs(min(nx,i)+1,j))-max(zb(i,j),zb(min(nx,i)+1,j)),par%eps)
          end if
       end do
    end do

    hu = max(hu,0.d0)

    do j=1,ny+1
       do i=1,nx+1
          ! Water depth in v-points do continuity equation: upwind
          if (vv(i,j)>par%umin) then
		     if (par%oldhu == 1) then 
                hv(i,j)=hh(i,j)
			 else
                hv(i,j)=zs(i,j)-max(zb(i,j),zb(i,min(ny,j)+1))
			 endif
          elseif (vv(i,j)<-par%umin) then
		     if (par%oldhu == 1) then 
                hv(i,j)=hh(i,min(ny,j)+1)
			 else
                hv(i,j)=zs(i,min(ny,j)+1)-max(zb(i,j),zb(i,min(ny,j)+1))
			 endif
          else
             hv(i,j)=max(max(zs(i,j),zs(i,min(ny,j)+1))-max(zb(i,j),zb(i,min(ny,j)+1)),par%eps)
          end if
       end do
    end do
    hv = max(hv,0.d0)

#ifdef USEMPI
    call xmpi_shift(hu ,'m:')
    call xmpi_shift(hv ,':n')
#endif


#ifndef USEMPI
    if (par%nonh==1) then
       ! do non-hydrostatic pressure compensation to solve short waves
       call nonh_cor(s,par)
    end if
#endif

    ! Flux in u-point
    qx=uu*hu
    ! Flux in v-points
    ! first column of qy is used later, and it is defined in the loop above
    ! no communication  necessary at this point
    qy=vv*hv
    !
    ! Add horizontal discharges
	!
    call discharge_boundary_h(s,par)
    !
    ! Update water level using continuity eq.
    !
    if (xmpi_isright) then
       jmax=ny
    else
       jmax=ny+1
    endif
    if (xmpi_isbot) then
       imax=nx
    else
       imax=nx+1
    endif
    if (ny>0) then
      do j=2,jmax
        do i=2,imax
          dzsdt(i,j) = (-1.d0)*( qx(i,j)*dnu(i,j)-qx(i-1,j)*dnu(i-1,j)  &
                               + qy(i,j)*dsv(i,j)-qy(i,j-1)*dsv(i,j-1) )*dsdnzi(i,j) &
                        - gww(i,j)
        end do
      end do
      zs(2:nx,2:ny) = zs(2:nx,2:ny)+dzsdt(2:nx,2:ny)*par%dt !Jaap nx instead of nx+1
    else
       j=1
       do i=2,imax
          dzsdt(i,j) = (-1.d0)*( qx(i,j)*dnu(i,j)-qx(i-1,j)*dnu(i-1,j) )*dsdnzi(i,j) &
                        - gww(i,j)
       end do
       zs(2:nx,1) = zs(2:nx,1)+dzsdt(2:nx,1)*par%dt !Jaap nx instead of nx+1
    endif !ny>0
    ! call discharge_boundary(s,par)
    !

#ifndef USEMPI
    if (par%secorder == 1) then
       !Second order correction
       call flow_secondorder_con(s,par,zs_old)
    endif
#endif
    !
    ! Lateral boundary conditions
    !
    if (ny>0) then
      ! RJ: Neumann water levels in case of right = 1 or right = 0
      do i=1,nx+1
         dzsdnavg=(zs0(i,ny+1)-zs0(i,1))/ndist(i,ny+1)
         ! Lateral boundary at y=0
         if (xmpi_isleft) then
           zs(i,1)=max(zs(i,2) - dzsdnavg*dnv(i,1),zb(i,1))
         endif
		 ! Lateral boundary at y=ny+1
         if (xmpi_isright) then
            zs(i,ny+1)=max(zs(i,ny) + dzsdnavg*dnv(i,ny),zb(i,ny+1))
         endif
       enddo
    endif !ny>0


    ! wwvv zs, uu, vv have to be communicated now, because they are used later on
#ifdef USEMPI
    call xmpi_shift(zs,':1')
    call xmpi_shift(zs,':n')
    call xmpi_shift(zs,'1:')
    call xmpi_shift(zs,'m:')
    call xmpi_shift(dzsdt,':1')  ! wwvv dzsdt maybe not necessary because first and last columns are not used
    call xmpi_shift(dzsdt,':n')
    call xmpi_shift(dzsdt,'1:')
    call xmpi_shift(dzsdt,'m:')
    call xmpi_shift(qx,':1')  ! wwvv qx maybe not necessary because first and last columns are not used
    call xmpi_shift(qx,':n')
    call xmpi_shift(qx,'1:')
    call xmpi_shift(qx,'m:')
    call xmpi_shift(qy,':1')
    call xmpi_shift(qy,':n')
    call xmpi_shift(qy,'1:')! wwvv qy maybe not necessary because first and last rows are not used
    call xmpi_shift(qy,'m:')
#endif

    if (par%secorder == 1) then
       vv_old = vv
       uu_old = uu
       zs_old = zs
    endif

    ! offshore boundary
    !
    ! U and V in cell centre; do output and sediment stirring
    !
    u(2:nx,:)=0.5d0*(uu(1:nx-1,:)+uu(2:nx,:))
    if(xmpi_istop) then
       u(1,:)=uu(1,:)
    endif
    if(xmpi_isbot) then
       u(nx+1,:)=u(nx,:)
    endif
    
#ifdef USEMPI
    call xmpi_shift(u,'1:')
    call xmpi_shift(u,'m:')
    call xmpi_shift(u,':1')
    call xmpi_shift(u,':n')
#endif

    if (ny>0) then
        v(:,2:ny)=0.5d0*(vv(:,1:ny-1)+vv(:,2:ny))
        if(xmpi_isleft) then
           v(:,1)=vv(:,1)
        endif
        if(xmpi_isright) then
           v(:,ny+1)=v(:,ny)        ! bas: need this for calculation of ee in wci routine
        endif
        !Ap
        v(nx+1,:)=v(nx,:)
    endif !ny>0

#ifdef USEMPI
    call xmpi_shift(v,':1')
    call xmpi_shift(v,':n')
    call xmpi_shift(v,'1:')
    call xmpi_shift(v,'m:')
#endif
    ! Robert + Jaap: compute derivatives of u and v
    !
    ! V-velocities at u-points
    if (ny>0) then
      vu(1:nx,2:ny)= 0.25d0*(vv(1:nx,1:ny-1)+vv(1:nx,2:ny)+ &
                             vv(2:nx+1,1:ny-1)+vv(2:nx+1,2:ny))
      ! how about boundaries?
      if(xmpi_isleft) then
        vu(:,1) = vu(:,2)
      endif
      if(xmpi_isright) then
        vu(:,ny+1) = vu(:,ny)
      endif
    else 
      vu(1:nx,1)= 0.5d0*(vv(1:nx,1)+vv(2:nx+1,1))
    endif !ny>0
    ! wwvv fill in vu(:1) and vu(:ny+1) for non-left and non-right processes
    !  and vu(nx+1,:)
#ifdef USEMPI
    call xmpi_shift(vu,':1')
    call xmpi_shift(vu,':n')
    call xmpi_shift(vu,'m:')
    ! Jaap: whu not vu,'1:' --> seesm to be necessary to compute vmagu that is used to compute Su
#endif
    vu=vu*wetu
    ! V-stokes velocities at U point
    if (ny>0) then
      vsu(1:nx,2:ny)=0.5d0*(ust(1:nx,2:ny)*sin(thetamean(1:nx,2:ny))+ &
                            ust(2:nx+1,2:ny)*sin(thetamean(2:nx+1,2:ny)))
      if(xmpi_isleft) then
        vsu(:,1)=vsu(:,2)
      endif
      if(xmpi_isright) then
        vsu(:,ny+1) = vsu(:,ny)
      endif
    else
      vsu(1:nx,1)=0.5d0*(ust(1:nx,1)*sin(thetamean(1:nx,1))+ &
                         ust(2:nx+1,1)*sin(thetamean(2:nx+1,1)))
    endif !ny>0
    ! wwvv same for vsu
#ifdef USEMPI
    call xmpi_shift(vsu,':1')
    call xmpi_shift(vsu,':n')
    call xmpi_shift(vsu,'m:')
#endif
    vsu = vsu*wetu
    ! U-stokes velocities at U point
    if (ny>0) then
      usu(1:nx,2:ny)=0.5d0*(ust(1:nx,2:ny)*cos(thetamean(1:nx,2:ny))+ &
                            ust(2:nx+1,2:ny)*cos(thetamean(2:nx+1,2:ny)))
      if(xmpi_isleft) then
        usu(:,1)=usu(:,2)
      endif
      if(xmpi_isright) then
        usu(:,ny+1)=usu(:,ny)
      endif
    else
      usu(1:nx,1)=0.5d0*(ust(1:nx,1)*cos(thetamean(1:nx,1))+ &
                         ust(2:nx+1,1)*cos(thetamean(2:nx+1,1)))
    endif !ny>0
    ! wwvv same for usu   
#ifdef USEMPI
    call xmpi_shift(usu,':1')
    call xmpi_shift(usu,':n')
    call xmpi_shift(usu,'m:')
#endif
    usu=usu*wetu

    ! V-euler velocities at u-point
    veu = vu - vsu
    ! U-euler velocties at u-point
    ueu = uu - usu
    ! Velocity magnitude at u-points
    ! vmagu=sqrt(uu**2+vu**2)
    ! Eulerian velocity magnitude at u-points
    vmageu=sqrt(ueu**2+veu**2)

    ! U-velocities at v-points
    if (ny>0) then
      uv(2:nx,1:ny)= .25d0*(uu(1:nx-1,1:ny)+uu(2:nx,1:ny)+ &
                            uu(1:nx-1,2:ny+1)+uu(2:nx,2:ny+1))
      ! boundaries?
      ! wwvv and what about uv(:,1) ?
      if(xmpi_isright) then
        uv(:,ny+1) = uv(:,ny)
      endif
    else
      uv(2:nx,1)= .5d0*(uu(1:nx-1,1)+uu(2:nx,1))
    endif !ny>0
    ! wwvv fix uv(:,ny+1) for non-right processes
    ! uv(1,:) and uv(nx+1,:) need to be filled in for
    ! non-bot or top processes
#ifdef USEMPI
    call xmpi_shift(uv,':n')
    call xmpi_shift(uv,':1')
    call xmpi_shift(uv,'1:')
    call xmpi_shift(uv,'m:')
#endif
    uv=uv*wetv
    ! V-stokes velocities at V point
    if (ny>0) then
      vsv(2:nx,1:ny)=0.5d0*(ust(2:nx,1:ny)*sin(thetamean(2:nx,1:ny))+&
                            ust(2:nx,2:ny+1)*sin(thetamean(2:nx,2:ny+1)))
      if(xmpi_isleft) then
        vsv(:,1) = vsv(:,2)
      endif
      if(xmpi_isright) then
        vsv(:,ny+1) = vsv(:,ny)
      endif
    else
      vsv(2:nx,1)= ust(2:nx,1)*sin(thetamean(2:nx,1))
    endif !ny>0
    ! wwvv fix vsv(:,1) and vsv(:,ny+1) and vsv(1,:) and vsv(nx+1,:)
#ifdef USEMPI
    call xmpi_shift(vsv,':n')
    call xmpi_shift(vsv,':1')
    call xmpi_shift(vsv,'1:')
    call xmpi_shift(vsv,'m:')
#endif

    vsv=vsv*wetv
    ! U-stokes velocities at V point
    if (ny>0) then
      usv(2:nx,1:ny)=0.5d0*(ust(2:nx,1:ny)*cos(thetamean(2:nx,1:ny))+&
                            ust(2:nx,2:ny+1)*cos(thetamean(2:nx,2:ny+1)))
      if(xmpi_isleft) then
        usv(:,1) = usv(:,2)
      endif
      if(xmpi_isleft) then
        usv(:,ny+1) = usv(:,ny)
      endif
    else
      usv(2:nx,1)=ust(2:nx,1)*cos(thetamean(2:nx,1))
    endif !ny>0
    ! wwvv fix usv(:,1) and usv(:,ny+1) and usv(1,:) and usv(nx+1,:)
#ifdef USEMPI
    call xmpi_shift(usv,':n')
    call xmpi_shift(usv,':1')
    call xmpi_shift(usv,'1:')
    call xmpi_shift(usv,'m:')
#endif
    usv=usv*wetv

    ! V-euler velocities at V-point
    vev = vv - vsv
    ! U-euler velocties at V-point
    uev = uv - usv
    ! Velocity magnitude at v-points
    ! vmagv=sqrt(uv**2+vv**2)
    ! Eulerian velocity magnitude at v-points
    vmagev=sqrt(uev**2+vev**2)


    ! Ue and Ve in cell centre; do output and sediment stirring
    ue(2:nx,:)=0.5d0*(ueu(1:nx-1,:)+ueu(2:nx,:))
    ue(1,:)=ueu(1,:)
    ! wwvv ue(nx+1,:) ?
#ifdef USEMPI
    call xmpi_shift(ue,':1')
    call xmpi_shift(ue,':n')
    call xmpi_shift(ue,'1:')
    call xmpi_shift(ue,'m:')
#endif

    if (ny>0) then
      ve(:,2:ny)=0.5d0*(vev(:,1:ny-1)+vev(:,2:ny)) !Jaap ny+1
      ve(:,1)=vev(:,1)
    else
      ve(:,1) = vev(:,1)
    endif !ny>0
    ! wwvv vev(nx+1,:) ?
#ifdef USEMPI
    call xmpi_shift(ve,':1')
    call xmpi_shift(ve,':n')
    call xmpi_shift(ve,'1:')
    call xmpi_shift(ve,'m:')
#endif
    !
    hold =hh    ! wwvv ?  hold is never else used
    !
    hh=max(zs-zb,par%eps)

    maxzs=max(zs,maxzs)
    minzs=min(zs,minzs)

  end subroutine flow_timestep


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

  include 's.ind'
  include 's.inp'

!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------

  !MPI WARNING -> Check loop indices
  if (ny>0) then
    do j=2,ny
      do i=2,nx
        dudx = (uu(i,j)-uu(i-1,j))/dsz(i,j)
        dudy = .5d0*(uu(i,j+1) - uu(i,j-1) + uu(i-1,j+1) - uu(i-1,j-1))/(dnv(i,j)+dnv(i,j-1))
        dvdx = .5d0*(vv(i+1,j) - vv(i-1,j) + vv(i+1,j-1) - vv(i-1,j-1))/(dsu(i,j)+dsu(i-1,j))
        dvdy = (vv(i,j)-vv(i,j-1))/dnz(i,j)
        Tau  = sqrt(2.0d0 * dudx**2+2.0d0 * dvdy**2 + (dvdx+dudy)**2)
        l    = 1.d0/dsdnzi(i,j)
        nuh(i,j) = par%nuh**2 * l * Tau * real(wetu(i,j)*wetu(i-1,j)*wetv(i,j)*wetv(i,j-1),kind=8)
      enddo
    enddo
    
    if (xmpi_isleft)  nuh(:,1)      = nuh(:,2)      ! Bas+Jaap: changed from boundaries=0.d0 to neumann
    if (xmpi_isright) nuh(:,ny+1)   = nuh(:,ny)
  
  else
  
    j = 1
    do i=2,nx
        dudx = (uu(i,j)-uu(i-1,j))/dsz(i,j)
        dvdx = (vv(i+1,j) - vv(i-1,j) )/(dsu(i,j)+dsu(i-1,j))
        Tau  = sqrt(2.0d0 * dudx**2 + dvdx**2)
        l    = dsz(i,j)**2
        nuh(i,j) = par%nuh**2 * l * Tau * real(wetu(i,j)*wetu(i-1,j),kind=8)
    enddo
  
  endif !ny>0

if (xmpi_istop)   nuh(1,:)      = nuh(2,:)      ! Bas+Jaap: changed from boundaries=0.d0 to neumann
if (xmpi_isbot)   nuh(nx+1,:)   = nuh(nx,:)

#ifdef USEMPI
  call xmpi_shift(nuh,'1:')
  call xmpi_shift(nuh,'m:')
  call xmpi_shift(nuh,':1')
  call xmpi_shift(nuh,':n')
#endif  
end subroutine visc_smagorinsky

end module flow_timestep_module
