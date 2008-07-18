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
#ifdef USEMPI
use mpi
#endif

IMPLICIT NONE

type(spacepars),target                  :: s
type(parameters)                        :: par

integer                                 :: i
integer                                 :: j
real*8,dimension(:,:),allocatable,save  :: vsu,usu,vsv,usv,veu,uev
real*8,dimension(:,:),allocatable,save  :: dzsdx,dzsdy
real*8,dimension(:,:),allocatable,save  :: ududx,vdvdy,udvdx,vdudy
real*8,dimension(:,:),allocatable,save  :: viscu,viscv
real*8,dimension(:,:),allocatable,save  :: us,vs  
real*8,dimension(:,:),allocatable,save  :: nuh
real*8                                  :: dudx1,dudx2,dudy1,dudy2
real*8                                  :: dvdy1,dvdy2,dvdx1,dvdx2  !Jaap

include 's.ind'
include 's.inp'

if (.not. allocated(vsu) ) then
   allocate (   vsu(s%nx+1,s%ny+1))
   allocate (   usu(s%nx+1,s%ny+1))
   allocate (   vsv(s%nx+1,s%ny+1))
   allocate (   usv(s%nx+1,s%ny+1))
   allocate (   veu(s%nx+1,s%ny+1))
   allocate (   uev(s%nx+1,s%ny+1))
   allocate ( dzsdx(s%nx+1,s%ny+1))
   allocate ( dzsdy(s%nx+1,s%ny+1))
   allocate ( ududx(s%nx+1,s%ny+1))
   allocate ( vdvdy(s%nx+1,s%ny+1))
   allocate ( udvdx(s%nx+1,s%ny+1))
   allocate ( vdudy(s%nx+1,s%ny+1))
   allocate ( viscu(s%nx+1,s%ny+1))
   allocate ( viscv(s%nx+1,s%ny+1))
   allocate (    us(s%nx+1,s%ny+1))
   allocate (    vs(s%nx+1,s%ny+1))
   allocate (   nuh(s%nx+1,s%ny+1))

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
    dzsdx   =0.d0
    dzsdy   =0.d0
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
    nuh     =0.d0
endif


    
    ! Water depth
    ! hh=zs-zb
    ! hh=max(hh,par%eps)
    ! Jaap Wetting and drying criterion eta points
    do j=1,ny+1
        do i=1,nx+1
            if(hh(i,j)>par%eps) then
                wetz(i,j)=1
            else
                wetz(i,j)=0
            end if
        end do
    end do
    ! zs=zs*wetz
    ! Water level slopes
    do j=1,ny+1
        do i=2,nx
                dzsdx(i,j)=(zs(i+1,j)-zs(i,j))/(xz(i+1)-xz(i))  
        end do
    end do
    do j=2,ny
        do i=1,nx+1
                dzsdy(i,j)=(zs(i,j+1)-zs(i,j))/(yz(j+1)-yz(j))    
        end do
    end do 
    do j=1,ny+1
        do i=1,nx+1 !Ap
            ! Water depth in u-points do momentum equation: mean
            hum(i,j)=max(.5d0*(hh(i,j)+hh(min(nx,i)+1,j)),par%eps)
            ! Water depth in u-points do continuity equation: upwind
            if (uu(i,j)>par%umin) then
                hu(i,j)=hh(i,j)
            elseif (uu(i,j)<-par%umin) then
                !Ap
                hu(i,j)=hh(min(nx,i)+1,j)  
            else
                !Ap
                hu(i,j)=max(max(zs(i,j),zs(min(nx,i)+1,j))-max(zb(i,j),zb(min(nx,i)+1,j)),par%eps)
                !hu(i,j) = (hh(i,j)+hh(min(nx,i)+1,j))/2
            end if           
        end do 
    end do
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
    do j=1,ny+1
        do i=1,nx+1
            ! Water depth in v-points do momentum equation: mean
            hvm(i,j)=max(.5d0*(hh(i,j)+hh(i,min(ny,j)+1)),par%eps)
            ! Water depth in v-points do continuity equation: upwind
            if (vv(i,j)>par%umin) then
                hv(i,j)=hh(i,j)
            elseif (vv(i,j)<-par%umin) then
                hv(i,j)=hh(i,min(ny,j)+1)
            else
                hv(i,j)=max(max(zs(i,j),zs(i,min(ny,j)+1))-max(zb(i,j),zb(i,min(ny,j)+1)),par%eps)
                !hv(i,j) = (hh(i,j)+hh(i,min(ny,j)+1))/2
            end if           
        end do 
    end do
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
    !
    ! X-direction
    !
    do j=2,ny ! Jaap 2,ny instead of 1,ny+1
        do i=2,nx
            ! Advection terms (momentum conserving method)
            if (uu(i,j)>=0.d0) then
               ududx(i,j)=.5d0*(hu(i,j)*uu(i,j)+hu(i-1,j)*uu(i-1,j))/&
                          hum(i,j)*(uu(i,j)-uu(i-1,j))/(xu(i)-xu(i-1))
            else
               ududx(i,j)=.5d0*(hu(i,j)*uu(i,j)+hu(i+1,j)*uu(i+1,j))/&
                          hum(i,j)*(uu(i+1,j)-uu(i,j))/(xu(i+1)-xu(i))
            end if           
        end do 
    end do
    !do j=2,ny
    !    do i=1,nx
    !        ! Advection terms (momentum conserving method) --> not exactly as in Stelling and Duinmeijer
    !        vdudy(i,j)=vu(i,j)*(uu(i,j+1)-uu(i,j-1))/(yz(j+1)-yz(j-1))
    !    end do 
    !end do
    ! Robert & Jaap: Let's do Stelling & Duinmeijer for all advection terms
    do j=2,ny
        do i=1,nx
            if (vu(i,j)>=0.d0) then
               vdudy(i,j)=.5d0*(hu(i,j)*vu(i,j)+hu(i,j-1)*vu(i,j-1))/&
                          hum(i,j)*(uu(i,j)-uu(i,j-1))/(yz(j)-yz(j-1))
            else
               vdudy(i,j)=.5d0*(hu(i,j)*vu(i,j)+hu(i,j+1)*vu(i,j+1))/&
                          hum(i,j)*(uu(i,j+1)-uu(i,j))/(yz(j+1)-yz(j))
            endif
        end do 
    end do
    !
    do j=2,ny
        do i=2,nx
            nuh(i,j) = par%nuh + par%nuhfac*hh(i,j)*(DR(i,j)/par%rho)**(1.0d0/3.0d0) !Robert en Jaap; increase eddy viscosity by wave induced breaking as in Reniers 2004
            dudx1 = (uu(i+1,j)-uu(i,j))/(xu(i+1)-xu(i))
            dudx2 = (uu(i,j)-uu(i-1,j))/(xu(i)-xu(i-1)) 
            viscu(i,j) = nuh(i,j)*( 2*(dudx1-dudx2)/(xu(i+1)-xu(i-1)) )*wetu(i+1,j)*wetu(i-1,j)  !Set viscu = 0.0 near water line            
        end do 
    end do
    do j=2,ny
        do i=2,nx
            dudy1 = (uu(i,j+1)-uu(i,j))/(yz(j+1)-yz(j))
            dudy2 = (uu(i,j)-uu(i,j-1))/(yz(j)-yz(j-1))
            viscu(i,j) = viscu(i,j) + nuh(i,j)*( 2*(dudy1-dudy2)/(yz(j+1)-yz(j-1)) )*wetu(i,j+1)*wetu(i,j-1)
        end do 
    end do
    !
    ! Explicit Euler step momentum u-direction
    !
    do j=2,ny
        do i=2,nx-1
            if(wetu(i,j)==1) then
                uu(i,j)=uu(i,j)-par%dt*(ududx(i,j)+vdudy(i,j)-viscu(i,j) & !Ap,Robert,Jaap 
                    + par%g*dzsdx(i,j) &
!                   + par%g/par%C**2.d0/hu(i,j)*vmageu(i,j)*ueu(i,j) & 
                                + par%cf/hu(i,j)*ueu(i,j)*sqrt((1.16d0*s%urms(i,j))**2+vmageu(i,j)**2) &    
                    - Fx(i,j)/par%rho/hu(i,j) &
                                        - par%fc*vu(i,j))
            else
                uu(i,j)=0.0d0
            end if
        end do 
    end do
    ! Flux in u-point
    qx=uu*hu
    !
    ! Y-direction
    !
    do j=1,ny
        do i=2,nx+1 !
            ! Advection terms (momentum conserving method)
            if (vv(i,j)>=0.d0) then
                                if (j>1) then
                  vdvdy(i,j)=.5d0*(hv(i,j)*vv(i,j)+hv(i,j-1)*vv(i,j-1))/hvm(i,j)*(vv(i,j)-vv(i,j-1))&
                  /(yv(j)-yv(j-1))
                                else
                                                vdvdy(i,j) = 0.0d0
                                endif
            else
               ! Dano no need to set vdvdy to 0 for v>0
               if (j<ny) then
                  vdvdy(i,j)=.5d0*(hv(i,j)*vv(i,j)+hv(i,j+1)*vv(i,j+1))/hvm(i,j)*(vv(i,j+1)-vv(i,j))&
                  /(yv(j+1)-yv(j))
               else
                  vdvdy(i,j)=0.d0
               endif
            end if           
        end do 
    end do
    !do j=1,ny !Jaap 1,ny instead of 1,ny+1 (v's at ny+1 are dummy's)
    !    do i=2,nx
    !        ! Advection terms (momentum conserving method) --> not exactly as in Stelling and Duinmeijer
    !        udvdx(i,j)=uv(i,j)*(vv(i+1,j)-vv(i-1,j))/(xz(i+1)-xz(i-1))
    !    end do 
    !end do
    ! Robert & Jaap: Let's do Stelling & Duinmeijer for all advection terms
    do j=1,ny
        do i=2,nx
            if (uv(i,j)>=0.d0) then
               udvdx(i,j)=.5d0*(hv(i,j)*uv(i,j)+hv(i-1,j)*uv(i-1,j))/&
                          hvm(i,j)*(vv(i,j)-vv(i-1,j))/(xz(i)-xz(i-1))
            else
               udvdx(i,j)=.5d0*(hv(i,j)*uv(i,j)+hv(i+1,j)*uv(i+1,j))/&
                          hvm(i,j)*(vv(i+1,j)-vv(i,j))/(xz(i+1)-xz(i))
            endif
        end do 
    end do
    do j=2,ny
        do i=2,nx
          dvdy1 = (vv(i,j+1)-vv(i,j))/(yv(j+1)-yv(j))
          dvdy2 = (vv(i,j)-vv(i,j-1))/(yv(j)-yv(j-1))
          viscv(i,j) = nuh(i,j)*( 2*(dvdy1-dvdy2)/(yv(j+1)-yv(j-1)) )*wetv(i,j+1)*wetv(i,j-1)
        end do 
    end do
    do j=2,ny
        do i=2,nx
            nuh(i,j) = par%nuhv*nuh(i,j)         !Robert en Ap: increase nuh interaction in d2v/dx2
            dvdx1 = (vv(i+1,j)-vv(i,j))/(xz(i+1)-xz(i))
            dvdx2 = (vv(i,j)-vv(i-1,j))/(xz(i)-xz(i-1))
            viscv(i,j) = viscv(i,j) + nuh(i,j)*( 2*(dvdx1-dvdx2)/(xz(i+1)-xz(i-1)) )*wetv(i+1,j)*wetv(i-1,j)
        end do 
    end do  

    udvdx(nx+1,:)=0.0d0       !Jaap udvdx(nx+1,:) is not defined but is used to compute vv(nx+1,:)
    viscv(nx+1,:)=viscv(nx,:) !Jaap viscv(nx+1,:) is not defined but is used to compute vv(nx+1,:)

    viscv(:,1)=viscv(:,2)   
    viscv(:,ny)=viscv(:,ny-1)
    !
    ! Explicit Euler step momentum v-direction
    !
    do j=1,ny
        do i=2,nx !jaap instead of nx+1        
            if(wetv(i,j)==1) then
                vv(i,j)=vv(i,j)-par%dt*(udvdx(i,j)+vdvdy(i,j)-viscv(i,j)& !Ap,Robert,Jaap 
                    + par%g*dzsdy(i,j)&
!                                       + par%g/par%C**2/hv(i,j)*vmagev(i,j)*vev(i,j)&
                                        + par%cf/hv(i,j)*vev(i,j)*sqrt((1.16d0*s%urms(i,j))**2+vmagev(i,j)**2) &   !Ruessink et al 2001
                                        - Fy(i,j)/par%rho/hv(i,j) &
                                        + par%fc*uv(i,j))
            else
                vv(i,j)=0.0d0
            end if
        end do 
    end do
    ! Flux in v-points
    qy=vv*hv
    !
        ! do non-hydrostatic pressure compensation to solve short waves
        !if (par%nonh==1) then
        !       call nhcorrection(s,par)
    !    qx=uu*hu
        !        qy=vv*hv
    !end if
        !
    ! Update water level using continuity eq.
    !
    do j=2,ny
        do i=2,nx
           dzsdt(i,j) = (-1.d0)*((qx(i,j)-qx(i-1,j))/(xu(i)-xu(i-1)) &
                               +(qy(i,j)-qy(i,j-1))/(yv(j)-yv(j-1)))
         end do
    end do
    !
    zs(2:nx,2:ny) = zs(2:nx,2:ny)+dzsdt(2:nx,2:ny)*par%dt !Jaap nx instead of nx+1
    !    
    ! Output
    !
    !
    ! Lateral boundary at y=0;
    !
    zs(1:nx+1,1)=zs(1:nx+1,2) - (zs0(:,ny+1)-zs0(:,1))/(yz(ny+1)-yz(1))*(yz(2)-yz(1))
    uu(1:nx+1,1)=uu(1:nx+1,2) !Jaap nx instead of nx+1
    if (par%right==1) then
       vv(2:nx+1,1) = 0.d0
    endif
    !
    ! Lateral boundary at y=ny*dy;
    !
    zs(1:nx+1,ny+1)=zs(1:nx+1,ny) + (zs0(:,ny+1)-zs0(:,1))/(yz(ny+1)-yz(1))*(yz(ny+1)-yz(ny))
    uu(1:nx+1,ny+1)=uu(1:nx+1,ny) !Jaap nx instead of nx+1
    if (par%left==1) then
       vv(2:nx+1,ny) = 0.d0
    endif
        ! offshore boundary
    !
    ! U and V in cell centre; do output and sediment stirring
    !
    u(2:nx,:)=0.5d0*(uu(1:nx-1,:)+uu(2:nx,:))
    u(1,:)=uu(1,:)
    !Ap
    u(nx+1,:)=u(nx,:)
    v(:,2:ny)=0.5d0*(vv(:,1:ny-1)+vv(:,2:ny)) !Jaap: ny+1
    v(:,1)=vv(:,1)
    
    !Ap
    v(nx+1,:)=v(nx,:)
        ! Robert + Jaap: compute derivaties of u and v
        !
    ! V-velocities at u-points
    vu(1:nx,2:ny)= &                         
        0.25d0*(vv(1:nx,1:ny-1)+vv(1:nx,2:ny)+ &
        vv(2:nx+1,1:ny-1)+vv(2:nx+1,2:ny))
    ! how about boundaries?
    vu(:,1) = vu(:,2)
    vu(:,ny+1) = vu(:,ny)
    vu=vu*wetu
    ! V-stokes velocities at U point
    vsu(1:nx,2:ny)=0.5d0*(ust(1:nx,2:ny)*sin(thetamean(1:nx,2:ny))+ &
        ust(2:nx+1,2:ny)*sin(thetamean(2:nx+1,2:ny)))
    vsu(:,1)=vsu(:,2)
    vsu(:,ny+1) = vsu(:,ny)
    vsu = vsu*wetu
    ! U-stokes velocities at U point
    usu(1:nx,2:ny)=0.5d0*(ust(1:nx,2:ny)*cos(thetamean(1:nx,2:ny))+ &
        ust(2:nx+1,2:ny)*cos(thetamean(2:nx+1,2:ny)))
    usu(:,1)=usu(:,2)
    usu(:,ny+1)=usu(:,ny)
    usu=usu*wetu
    
    ! V-euler velocities at u-point
    veu = vu - vsu
    ! U-euler velocties at u-point
    ueu = uu - usu
    ! Velocity magnitude at u-points
    vmagu=sqrt(uu**2+vu**2)
    ! Eulerian velocity magnitude at u-points
    vmageu=sqrt(ueu**2+veu**2)

    ! U-velocities at v-points
    uv(2:nx,1:ny)= &
        .25d0*(uu(1:nx-1,1:ny)+uu(2:nx,1:ny)+ &
        uu(1:nx-1,2:ny+1)+uu(2:nx,2:ny+1))
    ! boundaries?
    uv(:,ny+1) = uv(:,ny)
    uv=uv*wetv
     ! V-stokes velocities at V point
    vsv(2:nx,1:ny)=0.5d0*(ust(2:nx,1:ny)*sin(thetamean(2:nx,1:ny))+&
        ust(2:nx,2:ny+1)*sin(thetamean(2:nx,2:ny+1)))
    vsv(:,1) = vsv(:,2)
    vsv(:,ny+1) = vsv(:,ny)
    vsv=vsv*wetv
    ! U-stokes velocities at V point
    usv(2:nx,1:ny)=0.5d0*(ust(2:nx,1:ny)*cos(thetamean(2:nx,1:ny))+&
        ust(2:nx,2:ny+1)*cos(thetamean(2:nx,2:ny+1)))
    usv(:,1) = usv(:,2)
    usv(:,ny+1) = usv(:,ny)
    usv=usv*wetv

    ! V-euler velocities at V-point
    vev = vv - vsv
    ! U-euler velocties at V-point
    uev = uv - usv
    ! Velocity magnitude at v-points
    vmagv=sqrt(uv**2+vv**2)
     ! Eulerian velocity magnitude at v-points
    vmagev=sqrt(uev**2+vev**2)


    ! Ue and Ve in cell centre; do output and sediment stirring
    ue(2:nx,:)=0.5d0*(ueu(1:nx-1,:)+ueu(2:nx,:))
    ue(1,:)=ueu(1,:)
    ve(:,2:ny)=0.5d0*(vev(:,1:ny-1)+vev(:,2:ny)) !Jaap ny+1
    ve(:,1)=vev(:,1)
    !
    hold =hh
    !
    hh=max(zs-zb,par%eps)

        maxzs=max(zs,maxzs)
        minzs=min(zs,minzs)

end subroutine flow_timestep

end module flow_timestep_module
