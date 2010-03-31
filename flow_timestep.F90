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
integer                                 :: j
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

integer                                 :: imax,jmax

include 's.ind'
include 's.inp'

if (.not. allocated(vsu) ) then
   allocate (   vsu(s%nx+1,s%ny+1))
   allocate (   usu(s%nx+1,s%ny+1))
   allocate (   vsv(s%nx+1,s%ny+1))
   allocate (   usv(s%nx+1,s%ny+1))
   allocate (   veu(s%nx+1,s%ny+1))
   allocate (   uev(s%nx+1,s%ny+1))
   allocate ( ududx(s%nx+1,s%ny+1))
   allocate ( vdvdy(s%nx+1,s%ny+1))
   allocate ( udvdx(s%nx+1,s%ny+1))
   allocate ( vdudy(s%nx+1,s%ny+1))
   allocate ( viscu(s%nx+1,s%ny+1))
   allocate ( viscv(s%nx+1,s%ny+1))
   allocate (    us(s%nx+1,s%ny+1))
   allocate (    vs(s%nx+1,s%ny+1))

   if (par%secorder == 1) then
      allocate(vv_old(s%nx+1,s%ny+1)); vv_old = s%vv
      allocate(uu_old(s%nx+1,s%ny+1)); uu_old = s%uu
      allocate(zs_old(s%nx+1,s%ny+1)); zs_old = s%zs
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
endif
    ! update bedfriction coefficient cf
	if (par%cf < 0.d0) then
	   s%cf = par%g/(18.d0*log10(4.*s%hh/min(s%hh,s%D90top)))**2 ! cf = g/C^2 where C = 18*log(4*hh/D90) 
    endif

    ! zs=zs*wetz
    ! Water level slopes
    do j=1,ny+1
        do i=2,nx
                dzsdx(i,j)=(zs(i+1,j)-zs(i,j))/(xz(i+1)-xz(i))  
        end do
    end do
!    do j=2,ny
    do j=1,ny ! Dano need to get correct slope on boundary y=0
        do i=1,nx+1
                dzsdy(i,j)=(zs(i,j+1)-zs(i,j))/(yz(j+1)-yz(j))    
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
!            s%wetz(i,j) = min(1,s%wetu(max(i,2)-1,j)+s%wetu(i,j)+s%wetv(i,j)+s%wetv(i,max(j,2)-1))
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
    do j=2,ny ! Jaap 2,ny instead of 1,ny+1
        do i=2,nx         
			ududx(i,j)=max(0.d0,.5d0*(qx(i,j)+qx(i-1,j)) ) / hum(i,j)*(uu(i,j  )-uu(i-1,j))/(xu(i  )-xu(i-1)) + &
			           min(0.d0,.5d0*(qx(i,j)+qx(i+1,j)) ) / hum(i,j)*(uu(i+1,j)-uu(i,j  ))/(xu(i+1)-xu(i  ))
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
            vdudy(i,j)=max(0.d0,.5d0*(qy(i,j-1)+qy(i+1,j-1)) ) / hum(i,j)*(uu(i,j  )-uu(i,j-1))/(yz(j  )-yz(j-1)) + &
                       min(0.d0,.5d0*(qy(i,j  )+qy(i+1,j  )) ) / hum(i,j)*(uu(i,j+1)-uu(i,j  ))/(yz(j+1)-yz(j  ))
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
#ifndef USEMPI    
    if (par%smag == 1) then
      !Use smagorinsky subgrid model
      call visc_smagorinsky(s,par)    
    else
#endif    
      s%nuh = par%nuh
#ifndef USEMPI      
    endif  
#endif

    do j=2,ny
       do i=2,nx
          nuh(i,j) = max(s%nuh(i,j),par%nuhfac*hh(i,j)*(DR(i,j)/par%rho)**(1.0d0/3.0d0)) ! Ad: change to max
       end do 
    end do

    do j=2,ny
        do i=2,nx
            dudx1 = nuh(i+1,j)*s%hh(i+1,j)*(uu(i+1,j)-uu(i,j))/(xu(i+1)-xu(i))
            dudx2 = nuh(i,j)  *s%hh(i  ,j)*(uu(i,j)-uu(i-1,j))/(xu(i)-xu(i-1)) 
            viscu(i,j) = (1.0d0/s%hum(i,j))*( 2*(dudx1-dudx2)/(xu(i+1)-xu(i-1)) )*wetu(i+1,j)*wetu(i-1,j)  !Set viscu = 0.0 near water line   wwvv: viscu is overwritten in next loopnest before it is used
        end do 
    end do
    
#ifndef USEMPI    
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
#endif        
    
    do j=2,ny
        do i=2,nx
            !Nuh is defined at eta points, interpolate from four surrounding points
            nuh1  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,j+1)+nuh(i,j+1))
            nuh2  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,j-1)+nuh(i,j-1))
            
            dudy1 = nuh1 *.5d0*(s%hvm(i,j  )+s%hvm(i+1,j  ))*(uu(i,j+1)-uu(i,j))/(yz(j+1)-yz(j))
            dudy2 = nuh2 *.5d0*(s%hvm(i,j-1)+s%hvm(i+1,j-1))*(uu(i,j)-uu(i,j-1))/(yz(j)-yz(j-1))
            viscu(i,j) = viscu(i,j) + (1.0d0/s%hum(i,j))*( 2.0d0*(dudy1-dudy2)/(yz(j+1)-yz(j-1)) )*wetu(i,j+1)*wetu(i,j-1)
        end do 
    end do
    
#ifndef USEMPI    
  if (par%smag == 1) then
    do j=2,ny
        do i=2,nx
            !Nuh is defined at eta points, interpolate from four surrounding points
            nuh1  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,j+1)+nuh(i,j+1))
            nuh2  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,j-1)+nuh(i,j-1))
            
            dvdx1 = nuh1*.5d0*(s%hvm(i,j  )+s%hvm(i+1,j  ))*(s%vv(i+1,j  )-s%vv(i,j  ))/(s%xz(i+1)-s%xz(i))
            dvdx2 = nuh1*.5d0*(s%hvm(i,j-1)+s%hvm(i+1,j-1))*(s%vv(i+1,j-1)-s%vv(i,j-1))/(s%xz(i+1)-s%xz(i))
            viscu(i,j) = viscu(i,j) + (1.d0/s%hum(i,j))*(dvdx1-dvdx2)/(s%yv(j)-s%yv(j-1))*real(wetv(i+1,j) &
                 * wetv(i,j)*wetv(i+1,j-1)*wetv(i,j-1),8)
        enddo
    enddo    
  endif  
#endif      
    
    
    !
    ! Explicit Euler step momentum u-direction
    !
    do j=2,ny
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
			    taubx(i,j)=s%cf(i,j)*par%rho*ueu(i,j)*sqrt((1.16d0*s%urms(i,j))**2+vmageu(i,j)**2)
                uu(i,j)=uu(i,j)-par%dt*(ududx(i,j)+vdudy(i,j)-viscu(i,j) & !Ap,Robert,Jaap 
                    + par%g*dzsdx(i,j) &
!                   + par%g/par%C**2.d0/hu(i,j)*vmageu(i,j)*ueu(i,j) & 
                    + (taubx(i,j) - par%lwave*Fx(i,j))/(par%rho*hu(i,j)) &
                    - par%fc*vu(i,j) &
					- par%rhoa*par%Cd*cos(s%winddirnow)*s%windvnow**2)
            else
                uu(i,j)=0.0d0
            end if
        end do 
    end do
! wwvv since the loops range from 2:nx-1 , 2:ny we have to do something in the
! parallel case
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
    do j=1,ny
        do i=2,nx+1 !
            ! Advection terms (momentum conserving method)
            if (j>1 .and. j<ny) then
               vdvdy(i,j)=max(0.d0,.5d0*(qy(i,j)+qy(i,j-1)) ) / hvm(i,j)*(vv(i,j)-vv(i,j-1))/(yv(j)-yv(j-1)) + &
                          min(0.d0,.5d0*(qy(i,j)+qy(i,j+1)) ) / hvm(i,j)*(vv(i,j+1)-vv(i,j))/(yv(j+1)-yv(j))	          
            elseif (j==1) then
               vdvdy(i,j)=min(0.d0,.5d0*(qy(i,j)+qy(i,j+1)) ) / hvm(i,j)*(vv(i,j+1)-vv(i,j))/(yv(j+1)-yv(j))	                      
		      	elseif (j==ny) then
               vdvdy(i,j)=max(0.d0,.5d0*(qy(i,j)+qy(i,j-1)) ) / hvm(i,j)*(vv(i,j)-vv(i,j-1))/(yv(j)-yv(j-1))
            endif
         end do
    end do
    ! Robert & Jaap: Let's do Stelling & Duinmeijer for all advection terms
    do j=1,ny
        do i=2,nx
           udvdx(i,j)=max(0.d0,.5d0*(qx(i-1,j)+qx(i-1,j+1)) ) / hvm(i,j)*(vv(i  ,j)-vv(i-1,j))/(xz(i)-xz(i-1)) +&
                      min(0.d0,.5d0*(qx(i  ,j)+qx(i  ,j+1)) ) / hvm(i,j)*(vv(i+1,j)-vv(i  ,j))/(xz(i+1)-xz(i))
        end do 
    end do
    do j=2,ny
        do i=2,nx
          dvdy1 = nuh(i,j+1)*s%hh(i,j+1)*(vv(i,j+1)-vv(i,j))/(yv(j+1)-yv(j))
          dvdy2 = nuh(i,j)  *s%hh(i,j  )*(vv(i,j)-vv(i,j-1))/(yv(j)-yv(j-1))
          viscv(i,j) = (1.0d0/s%hvm(i,j))*( 2*(dvdy1-dvdy2)/(yv(j+1)-yv(j-1)) )*wetv(i,j+1)*wetv(i,j-1)
        end do 
    end do
    
#ifndef USEMPI    
    if (par%smag == 1) then
      viscv = 2.0d0*viscv
    endif    
#endif      
    
    nuh = par%nuhv*nuh !Robert en Ap: increase nuh interaction in d2v/dx2
    do j=2,ny
        do i=2,nx   
            !Nuh is defined at eta points, interpolate from four surrounding points
            nuh1  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,j+1)+nuh(i,j+1))
            nuh2  = .25d0*(nuh(i,j)+nuh(i-1,j)+nuh(i-1,j+1)+nuh(i,j+1))            
            
            dvdx1 = nuh1*.5d0*(s%hum(i  ,j)+s%hum(i  ,j+1))*(vv(i+1,j)-vv(i,j))/(xz(i+1)-xz(i))
            dvdx2 = nuh2*.5d0*(s%hum(i-1,j)+s%hum(i-1,j+1))*(vv(i,j)-vv(i-1,j))/(xz(i)-xz(i-1))
            viscv(i,j) = viscv(i,j) + (1.0d0/s%hvm(i,j))*( 2*(dvdx1-dvdx2)/(xz(i+1)-xz(i-1)) )*wetv(i+1,j)*wetv(i-1,j)
        end do 
    end do
    
#ifndef USEMPI    
    if (par%smag == 1) then
      do j=2,ny
          do i=2,nx
              !Nuh is defined at eta points, interpolate from four surrounding points
              nuh1  = .25d0*(nuh(i,j)+nuh(i+1,j)+nuh(i+1,j+1)+nuh(i,j+1))
              nuh2  = .25d0*(nuh(i,j)+nuh(i-1,j)+nuh(i-1,j+1)+nuh(i,j+1)) 
              
              dudy1 = nuh1 *.5d0*(s%hum(i  ,j)+s%hum(i  ,j+1))*(s%uu(i,j+1  )-s%uu(i,j  ))/(s%yz(j+1)-s%yz(j))
              dudy2 = nuh2 *.5d0*(s%hum(i-1,j)+s%hum(i-1,j+1))*(s%uu(i-1,j+1)-s%uu(i-1,j))/(s%yz(j+1)-s%yz(j))
              viscv(i,j) = viscv(i,j) + (1.d0/s%hvm(i,j))*(dudy1-dudy2)/(s%xu(i)-s%xu(i-1)) &
                   * real(wetu(i,j+1)*wetu(i,j)*wetu(i-1,j+1)*wetv(i-1,j),8)
          enddo
      enddo    
    endif  
#endif       

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
                tauby(i,j)=s%cf(i,j)*par%rho*vev(i,j)*sqrt((1.16d0*s%urms(i,j))**2+vmagev(i,j)**2) !Ruessink et al, 2001
                vv(i,j)=vv(i,j)-par%dt*(udvdx(i,j)+vdvdy(i,j)-viscv(i,j)& !Ap,Robert,Jaap 
                       + par%g*dzsdy(i,j)&
                      ! + par%g/par%C**2/hv(i,j)*vmagev(i,j)*vev(i,j)&
                       + (tauby(i,j)- par%lwave*Fy(i,j))/(par%rho*hv(i,j)) &
                       + par%fc*uv(i,j) &
                       - par%rhoa*par%Cd*sin(s%winddirnow)*s%windvnow**2)

            else
                vv(i,j)=0.0d0
            end if
        end do 
    end do
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

    
    ! Pieter and Jaap: update hu en hv for continuity
    do j=1,ny+1
      do i=1,nx+1 !Ap
        ! Water depth in u-points do continuity equation: upwind
        if (uu(i,j)>par%umin) then
       !   hu(i,j)=hh(i,j)
		      hu(i,j)=zs(i,j)-max(zb(i,j),zb(min(nx,i)+1,j))
        elseif (uu(i,j)<-par%umin) then
      !    hu(i,j)=hh(min(nx,i)+1,j)  
		       hu(i,j)=zs(min(nx,i)+1,j)-max(zb(i,j),zb(min(nx,i)+1,j))
        else
          hu(i,j)=max(max(zs(i,j),zs(min(nx,i)+1,j))-max(zb(i,j),zb(min(nx,i)+1,j)),par%eps)          
        end if
		  end do
	  end do
	
    do j=1,ny+1
      do i=1,nx+1			 
			! Water depth in v-points do continuity equation: upwind
        if (vv(i,j)>par%umin) then
          !hv(i,j)=hh(i,j)
		      hv(i,j)=zs(i,j)-max(zb(i,j),zb(i,min(ny,j)+1))
        elseif (vv(i,j)<-par%umin) then
          !hv(i,j)=hh(i,min(ny,j)+1)
		      hv(i,j)=zs(i,min(ny,j)+1)-max(zb(i,j),zb(i,min(ny,j)+1))
        else
          hv(i,j)=max(max(zs(i,j),zs(i,min(ny,j)+1))-max(zb(i,j),zb(i,min(ny,j)+1)),par%eps)
        end if           
      end do 
    end do
#ifdef USEMPI
    call xmpi_shift(hu ,'m:')
	call xmpi_shift(hv ,':n')
#endif	

    ! 
    ! Dano in case of closed boundaries we need to set vv to 0 before computing qv,
    ! to avoid mass errors
	! Lateral boundary at y=0
    if (xmpi_isleft) then !Dano/Robert only on outer boundary
       if (par%right==1) then
         vv(2:nx+1,1) = 0.d0
	   else
	     vv(2:nx+1,1) = vv(2:nx+1,2) ! RJ: 
       endif
       uu(1:nx+1,1)=uu(1:nx+1,2) ! RJ: can also be done after continuity but more appropriate here
    endif
	! Lateral boundary at y=ny*dy 
	if (xmpi_isright) then !Dano/Robert only at outer boundary
       if (par%left==1) then
         vv(2:nx+1,ny) = 0.d0       
       else
	     vv(2:nx+1,ny) = vv(2:nx+1,ny-1) ! RJ
       endif
	   uu(1:nx+1,ny+1)=uu(1:nx+1,ny)
	endif
	
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

    call discharge_boundary(s,par)
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
    do j=2,jmax
        do i=2,imax
           dzsdt(i,j) = (-1.d0)*((qx(i,j)-qx(i-1,j))/(xu(i)-xu(i-1))  &
                               + (qy(i,j)-qy(i,j-1))/(yv(j)-yv(j-1))) &
							   - s%gww(i,j)
         end do
    end do
    ! call discharge_boundary(s,par)
    !
    zs(2:nx,2:ny) = zs(2:nx,2:ny)+dzsdt(2:nx,2:ny)*par%dt !Jaap nx instead of nx+1
    
#ifndef USEMPI
	if (par%secorder == 1) then
      !Second order correction
      call flow_secondorder_con(s,par,zs_old)
    endif	
#endif
    !    
    ! Output
    !
	! RJ: Neumann water levels in case of right = 1 or right = 0
    ! Lateral boundary at y=0

if(xmpi_isleft) then
    zs(1:nx+1,1)=max(zs(1:nx+1,2) - (zs0(:,ny+1)-zs0(:,1))/(yz(ny+1)-yz(1))*(yz(2)-yz(1)),zb(1:nx+1,1))
endif
if(xmpi_isright) then
    ! Lateral boundary at y=ny*dy
    zs(1:nx+1,ny+1)=max(zs(1:nx+1,ny) + (zs0(:,ny+1)-zs0(:,1))/(yz(ny+1)-yz(1))*(yz(ny+1)-yz(ny)),zb(1:nx+1,ny))
endif
    
   
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
      vv_old = s%vv
      uu_old = s%uu
      zs_old = s%zs
    endif

        ! offshore boundary
    !
    ! U and V in cell centre; do output and sediment stirring
    !
    u(2:nx,:)=0.5d0*(uu(1:nx-1,:)+uu(2:nx,:))
if(xmpi_istop) then
      u(1,:)=uu(1,:)
endif
    !Ap
if(xmpi_isbot) then
      u(nx+1,:)=u(nx,:)
endif
#ifdef USEMPI
    call xmpi_shift(u,'1:')
    call xmpi_shift(u,'m:')
    call xmpi_shift(u,':1')
    call xmpi_shift(u,':n')
#endif

    v(:,2:ny)=0.5d0*(vv(:,1:ny-1)+vv(:,2:ny)) !Jaap: ny+1
    v(:,1)=vv(:,1)
    
    !Ap
    v(nx+1,:)=v(nx,:)   ! wwvv symmetry: why not v(:,ny+1) = v(:,ny)

#ifdef USEMPI
    call xmpi_shift(v,':1')
    call xmpi_shift(v,':n')
    call xmpi_shift(v,'1:') 
    call xmpi_shift(v,'m:') 
#endif
    ! Robert + Jaap: compute derivatives of u and v
    !
    ! V-velocities at u-points
    vu(1:nx,2:ny)= &                         
        0.25d0*(vv(1:nx,1:ny-1)+vv(1:nx,2:ny)+ &
        vv(2:nx+1,1:ny-1)+vv(2:nx+1,2:ny))
    ! how about boundaries?
if(xmpi_isleft) then
    vu(:,1) = vu(:,2)
endif
if(xmpi_isright) then
    vu(:,ny+1) = vu(:,ny)
endif
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
    vsu(1:nx,2:ny)=0.5d0*(ust(1:nx,2:ny)*sin(thetamean(1:nx,2:ny))+ &
        ust(2:nx+1,2:ny)*sin(thetamean(2:nx+1,2:ny)))
if(xmpi_isleft) then
    vsu(:,1)=vsu(:,2)
endif
if(xmpi_isright) then
    vsu(:,ny+1) = vsu(:,ny)
    ! wwvv same for vsu
endif
#ifdef USEMPI
    call xmpi_shift(vsu,':1')
    call xmpi_shift(vsu,':n')
    call xmpi_shift(vsu,'m:')
#endif
    vsu = vsu*wetu
    ! U-stokes velocities at U point
    usu(1:nx,2:ny)=0.5d0*(ust(1:nx,2:ny)*cos(thetamean(1:nx,2:ny))+ &
        ust(2:nx+1,2:ny)*cos(thetamean(2:nx+1,2:ny)))
if(xmpi_isleft) then
    usu(:,1)=usu(:,2)
endif
if(xmpi_isright) then
    usu(:,ny+1)=usu(:,ny)
endif
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
    uv(2:nx,1:ny)= &
        .25d0*(uu(1:nx-1,1:ny)+uu(2:nx,1:ny)+ &
        uu(1:nx-1,2:ny+1)+uu(2:nx,2:ny+1))
    ! boundaries?
    ! wwvv and what about uv(:,1) ?
if(xmpi_isright) then
    uv(:,ny+1) = uv(:,ny)
endif
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
    vsv(2:nx,1:ny)=0.5d0*(ust(2:nx,1:ny)*sin(thetamean(2:nx,1:ny))+&
        ust(2:nx,2:ny+1)*sin(thetamean(2:nx,2:ny+1)))
if(xmpi_isleft) then
    vsv(:,1) = vsv(:,2)
endif
if(xmpi_isright) then
    vsv(:,ny+1) = vsv(:,ny)
endif
    ! wwvv fix vsv(:,1) and vsv(:,ny+1) and vsv(1,:) and vsv(nx+1,:)
#ifdef USEMPI
    call xmpi_shift(vsv,':n')
    call xmpi_shift(vsv,':1')
    call xmpi_shift(vsv,'1:')
    call xmpi_shift(vsv,'m:')
#endif

    vsv=vsv*wetv
    ! U-stokes velocities at V point
    usv(2:nx,1:ny)=0.5d0*(ust(2:nx,1:ny)*cos(thetamean(2:nx,1:ny))+&
        ust(2:nx,2:ny+1)*cos(thetamean(2:nx,2:ny+1)))
if(xmpi_isleft) then
    usv(:,1) = usv(:,2)
endif
if(xmpi_isleft) then
    usv(:,ny+1) = usv(:,ny)
endif
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

    ve(:,2:ny)=0.5d0*(vev(:,1:ny-1)+vev(:,2:ny)) !Jaap ny+1
    ve(:,1)=vev(:,1)
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

#ifndef USEMPI

subroutine visc_smagorinsky(s,par)
  use params
  use spaceparams

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
  real*8                                                  :: dx   !Local gridsize in x-dir
  real*8                                                  :: dy   !Local gridsize in y-dir
  real*8                                                  :: l   !Local gridsize in y-dir  
  integer                                                 :: i    !Index variable
  integer                                                 :: j    !Index variable

!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------      
      
  !MPI WARNING -> Check loop indices
  do j=2,s%ny
    do i=2,s%nx     
      dx   = (s%xu(i)-s%xu(i-1))
      dy   = (s%yv(j)-s%yv(j-1))
      dudx = (s%uu(i,j)-s%uu(i-1,j))/dx
      dudy = .5d0*(s%uu(i,j+1) - s%uu(i,j-1) + s%uu(i-1,j+1) - s%uu(i-1,j-1))/(s%yz(j+1)-s%yz(j-1))
      dvdx = .5d0*(s%vv(i+1,j) - s%vv(i-1,j) + s%vv(i+1,j-1) - s%vv(i-1,j-1))/(s%xz(i+1)-s%xz(i-1))
      dvdy = (s%vv(i,j)-s%vv(i,j-1))/dy
      Tau  = sqrt(2.0d0 * dudx**2+2.0d0 * dvdy**2 + (dvdx+dudy)**2)      
      l    = dx * dy
      s%nuh(i,j) = par%nuh**2 * l * Tau * real(s%wetu(i,j)*s%wetu(i-1,j)*s%wetv(i,j)*s%wetv(i,j-1),kind=8)
    enddo
  enddo
  
  s%nuh(1,:)      = 0.0d0
  s%nuh(:,s%ny+1) = 0.0d0
  s%nuh(:,1)      = 0.0d0
  s%nuh(s%nx+1,:) = 0.0d0
end subroutine visc_smagorinsky

#endif
end module flow_timestep_module
