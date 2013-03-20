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
module groundwaterflow
  use typesandkinds
contains



  subroutine gw_init(s,par)
    use params
    use xmpi_module
    use spaceparams
    use readkey_module
    use logging_module

    IMPLICIT NONE

    type(parameters)                            :: par
    type(spacepars)                             :: s

    character(slen)                              :: fname
    real*8                                      :: aquiferbot,temp
    integer                                     :: i,j,ier


    allocate (s%gwhead(s%nx+1,s%ny+1))
    allocate (s%gwlevel(s%nx+1,s%ny+1))
    allocate (s%gwheight(s%nx+1,s%ny+1))
    allocate (s%gwu(s%nx+1,s%ny+1))
    allocate (s%gwv(s%nx+1,s%ny+1))
    allocate (s%gww(s%nx+1,s%ny+1))
    allocate (s%gwbottom(s%nx+1,s%ny+1))
    allocate (s%dinfil(s%nx+1,s%ny+1))
    allocate (s%gw0back(2,s%ny+1))

    s%gww=0.d0
    s%gw0back=0.d0
    s%gwlevel=0.d0
    s%gwheight=0.d0
    s%gwhead=0.d0
    s%gwu=0.d0
    s%gwv=0.d0
    s%gww=0.d0
    s%dinfil=0.d0
    s%gwbottom=0.d0

    if (par%gwflow==1) then
       fname = readkey_name('params.txt','aquiferbotfile',bcast=.false.)
       if (fname=='') then     ! Not a filename
          temp = minval(s%zb)
          aquiferbot = readkey_dbl('params.txt','aquiferbot',temp-3.d0,-100.d0,100.d0,bcast=.false.)
          s%gwbottom=aquiferbot
       else
          open(31,file=fname)
          do j=1,s%ny+1
             read(31,*,iostat=ier)(s%gwbottom(i,j),i=1,s%nx+1)
             if (ier .ne. 0) then
                call report_file_read_error(fname)
             endif 
          end do
          close(31)
       endif

       fname = readkey_name('params.txt','gw0file',bcast=.false.)
       if (fname=='') then     ! Not a filename
          temp = readkey_dbl('params.txt','gw0',0.d0,-5.d0,5.d0,bcast=.false.)
          s%gwhead=temp
       else
          open(31,file=fname)
          do j=1,s%ny+1
             read(31,*,iostat=ier)(s%gwhead(i,j),i=1,s%nx+1)
             if (ier .ne. 0) then
                call report_file_read_error(fname)
             endif
          end do
          close(31)
       endif

       s%gw0back=s%gwhead(s%nx:s%nx+1,:)
       s%gwlevel=min(s%zb,s%gwhead)
       s%gwlevel=max(s%gwlevel,s%gwbottom+par%eps)
       s%gwheight=s%gwlevel-s%gwbottom
       s%gwu=0.d0
       s%gwv=0.d0
       s%gww=0.d0
       s%dinfil=max(par%dwetlayer/3.d0,0.02)   ! Centroid of area influenced instantly by free surface level lies at dwetlayer/3
    endif
  end subroutine gw_init



  subroutine gw_bc(s,par)

    use params
    use xmpi_module
    use spaceparams

    IMPLICIT NONE

    type(parameters)                            :: par
    type(spacepars)                             :: s


    if(xmpi_istop) then
       s%gwhead(1,:)=s%zs(1,:)
    elseif (xmpi_isbot) then
       if (par%tideloc==4 .or. (par%tideloc==2 .and. trim(par%paulrevere)=='land')) then
          s%gwhead(s%nx+1,:)=s%zs(s%nx+1,:)
       else
          s%gwhead(s%nx+1,:)=s%gw0back(2,:)
       endif
    endif

    if (s%ny>0) then 
       s%gwhead(:,1)=s%gwhead(:,2)
       s%gwhead(:,s%ny+1)=s%gwhead(:,s%ny)
    endif

#ifdef USEMPI
    call xmpi_shift(s%gwhead,':1')
    call xmpi_shift(s%gwhead,':n')
    call xmpi_shift(s%gwhead,'1:')
    call xmpi_shift(s%gwhead,'m:')
#endif

    s%gwlevel(1,:)=min(s%gwhead(1,:),s%zb(1,:))
    s%gwlevel(s%nx+1,:)=min(s%gwhead(s%nx+1,:),s%zb(s%nx+1,:))
    if (s%ny>0) then
       s%gwlevel(:,1)=min(s%gwhead(:,1),s%zb(:,1))
       s%gwlevel(:,s%ny+1)=min(s%gwhead(:,s%ny+1),s%zb(:,s%ny+1))
    endif

    s%gwbottom=min(s%gwbottom,s%zb-par%eps)

  end subroutine gw_bc



  subroutine gwflow(s,par)

    use params
    use xmpi_module
    use spaceparams

    IMPLICIT NONE

    type(parameters)                            :: par
    type(spacepars)                             :: s

    integer                                     :: i,j
    real*8,dimension(:,:),allocatable,save      :: dheaddx,dheaddy,dleveldt,gwqx,gwqy,gwhu,gwhv
    real*8,dimension(:,:),allocatable,save      :: infilhorgw,infilhorsw
    real*8,dimension(:,:),allocatable,save      :: c1,c2,r,fsh

    if (.not. allocated(dheaddx)) then
       allocate(dheaddx(s%nx+1,s%ny+1))
       allocate(dheaddy(s%nx+1,s%ny+1))
       allocate(dleveldt(s%nx+1,s%ny+1))
       allocate(gwqx(s%nx+1,s%ny+1))
       allocate(gwqy(s%nx+1,s%ny+1))
       allocate(gwhu(s%nx+1,s%ny+1))
       allocate(gwhv(s%nx+1,s%ny+1))
       allocate(c1(s%nx+1,s%ny+1))
       allocate(c2(s%nx+1,s%ny+1))
       allocate(r(s%nx+1,s%ny+1))
       allocate(fsh(s%nx+1,s%ny+1))
       allocate(infilhorgw(1:s%nx+1,1:s%ny+1))
       allocate(infilhorsw(1:s%nx+1,1:s%ny+1))
    endif

    ! Free surface head and ratio free surface to groundwater head to be used
    fsh=(s%zs-s%zb)
    r=(s%zb-s%gwhead)/(par%dwetlayer)
    where (r<0.d0)
       r=0.d0
    elsewhere (r>1.d0)
       r=1.d0
    endwhere


    ! Momentum balance
    ! Determine pressure gradients

    ! Update groundwater head
    ! note smoothing over 2 timesteps
    where (s%wetz==1 .and. s%gwlevel>s%zb-par%dwetlayer)
       s%gwhead=(s%gwhead+s%gwlevel+(s%zs-s%gwlevel)*(1.d0-r))/2
    elsewhere
       s%gwhead=(s%gwhead+s%gwlevel)/2
    endwhere

    dheaddx=0.d0
    dheaddy=0.d0

    do j=1,s%ny+1
       do i=1,s%nx
          dheaddx(i,j)=(s%gwhead(i+1,j)-s%gwhead(i,j))/s%dsu(i,j)  
       end do
    end do
    
    do j=1,s%ny
       do i=1,s%nx+1
          dheaddy(i,j)=(s%gwhead(i,j+1)-s%gwhead(i,j))/s%dnv(i,j)   
       end do
    end do
 
    s%gwheight=s%gwlevel-s%gwbottom

    ! Determine intermediate aquifer depths
    gwhu(1:s%nx,:)=0.5d0*(s%gwheight(1:s%nx,:)+s%gwheight(2:s%nx+1,:))
    if (s%ny>0) then
       gwhv(:,1:s%ny)=0.5d0*(s%gwheight(:,1:s%ny)+s%gwheight(:,2:s%ny+1))
    else
       gwhv = 0.d0
    endif

    ! Determine fluxes
    s%gwu=-par%kx*dheaddx
    s%gwv=-par%ky*dheaddy

    ! Limit for stability in case of very high kx values
    !s%gwu(1:s%nx,1:s%ny)=min(s%gwu(1:s%nx,1:s%ny),0.5d0*(s%xz(2:s%nx+1,1:s%ny)-s%xz(1:s%nx,1:s%ny))/par%dt)
    !s%gwv(1:s%nx,1:s%ny)=min(s%gwv(1:s%nx,1:s%ny),0.5d0*(s%yz(1:s%nx,2:s%ny+1)-s%yz(1:s%nx,1:s%ny))/par%dt)

    gwqx=s%gwu*gwhu
    gwqy=s%gwv*gwhv

    ! Stop cells from drying up
    if (s%ny>0) then
       where(s%gwlevel(2:s%nx,2:s%ny)<=s%gwbottom(2:s%nx,2:s%ny)+par%eps)
          gwqx(2:s%nx,2:s%ny)=min(gwqx(2:s%nx,2:s%ny),0.d0)
          gwqx(1:s%nx-1,2:s%ny)=max(gwqx(1:s%nx-1,2:s%ny),0.d0)
          gwqy(2:s%nx,2:s%ny)=min(gwqy(2:s%nx,2:s%ny),0.d0)
          gwqy(2:s%nx,1:s%ny-1)=max(gwqy(2:s%nx,1:s%ny-1),0.d0)
       end where
    else
       where(s%gwlevel(2:s%nx,1)<=s%gwbottom(2:s%nx,1)+par%eps)
          gwqx(2:s%nx,1)=min(gwqx(2:s%nx,1),0.d0)
          gwqx(1:s%nx-1,1)=max(gwqx(1:s%nx-1,1),0.d0)
       end where
    endif

    ! Based on old groundwater levels, interaction with free water calculated
    ! This could be done by a double do-loop with if statements, but logical indexing prob. faster
    c1=0.d0
    c2=0.d0
    !c3=0.d0

    where (s%gwlevel>=s%zb)   ! Water permeates out
       c1=1.d0
    elsewhere (s%wetz==1)    ! Water can permeate in
       c2=1.d0
    endwhere

    ! Assume that infiltration layers with no water on top drain out of the way of subsequent infiltration actions
    ! But maintain minimum layer thickness to prevent numerical instability
    where(s%wetz==0)
       s%dinfil=par%dwetlayer/3.d0        ! Centroid of area influenced instantly by free surface level lies at dwetlayer/3
    elsewhere
       s%dinfil=min(s%dinfil,s%zb-s%gwlevel)
       s%dinfil=max(s%dinfil,par%dwetlayer/3.d0)
    endwhere

    s%gww=0.d0 ! w defined positive from sea to groundwater in volumes of surface water (no pores).
    s%gww=par%por*(&
         -( c1* (s%gwlevel-s%zb) / par%dt                                      )&
         +( c2* ((1.d0-r)*(s%zb-s%gwlevel) / par%dt + r*(par%kz*(1.d0 + fsh/s%dinfil))/par%por ) ) ) ! Jaap: add effect gravity for computing gww
    !             +(case2*((1.d0-r)*max((s%zb-s%gwhead),0.d0)+r*(fsh*gw%kper)))&
    ! +(case3*fsh*gw%kper)&

    ! ensure that water extracted from surface layer is not more than available
    where (s%gww*par%dt>s%hh)
       s%gww=s%hh/par%dt
    endwhere

    dleveldt=0.d0
    ! Mass balance
    if (s%ny>0) then
       do j=2,s%ny
          do i=2,s%nx
             dleveldt(i,j)=-1.d0*(gwqx(i,j)*s%dnu(i,j)-gwqx(i-1,j)*s%dnu(i-1,j) + &
                  gwqy(i,j)*s%dsv(i,j)-gwqy(i,j-1)*s%dsv(i,j-1))*s%dsdnzi(i,j)/par%por &
                  +1.d0*s%gww(i,j)/par%por
          enddo
       enddo
    else
       do i=2,s%nx
          dleveldt(i,1)=-1.d0*(gwqx(i,1)-gwqx(i-1,1))/s%dsz(i,1)/par%por &
               +1.d0*s%gww(i,1)/par%por
       enddo
    endif

    s%gwlevel=s%gwlevel+dleveldt*par%dt

    ! Update quasi vertical model infiltration layer thickness
    s%dinfil=s%dinfil+s%gww*par%dt/par%por
    
    if (par%gwhorinfil==1) then
       ! post-update horizontal infiltration computation
       dleveldt = 0.d0
       infilhorgw = 0.d0
       infilhorsw = 0.d0
       call gw_horizontal_infil_exfil(s,par,infilhorgw,infilhorsw,par%kx,par%ky)
       ! update water levels
       s%gwlevel=s%gwlevel+infilhorgw*par%dt
       ! update for mass balance in flow
       s%gww = s%gww+infilhorsw
    endif

  end subroutine gwflow


  subroutine gw_horizontal_infil_exfil(s,par,infilhorgw,infilhorsw,Kx,Ky)
     ! compute horizontal part of infiltration/exfiltration
     use params
     use xmpi_module
     use spaceparams
     
     type(parameters),intent(in)                 :: par
     type(spacepars)                             :: s
     real*8,intent(in)                           :: Kx,Ky
     real*8,dimension(s%nx+1,s%ny+1),intent(out) :: infilhorgw,infilhorsw
     ! internal
     real*8,dimension(s%nx+1,s%ny+1)             :: dheaddx,dheaddy
     integer                                     :: i,j
     real*8                                      :: dz,dx,dy,vel,hsurf,hgw,dhead
     real*8                                      :: scaleareab,scaleareas
     
     !
     ! compute horizontal flow velocity between bed and surface water VERTICAL interfaces
     ! 
     ! loop over all u-interfaces
     do i=1,s%nx
        do j=1,s%ny+1
          dz = s%zb(i+1,j)-s%zb(i,j)
          if (dz>0.d0 .and. s%wetz(i,j)==1) then ! jump up
             hsurf = s%zs(i,j)
             hgw = s%gwhead(i+1,j)
             dx = s%dsz(i+1,j)
             dhead = hsurf-hgw         ! positive gradient if hsurf>hgw
                                       ! negative velocitity if infiltration
             vel = dhead/dx*Kx         ! convert negative down velocity into positive if infiltration,
                                       ! negative if exfiltration 
             ! horizontal infiltration occurs along small side section, so scale the horizontal
             ! infiltration area by the vertical area of the cell used in the mass balance equation
             ! and infiltration velocities
             ! 
             ! Scale area for bed cell and surface water cell
             scaleareab = dz/s%dsz(i+1,j)
             scaleareas = dz/s%dsz(i,j)
             ! Limit  flow to that available in space in ground and 
             ! in surface water volume
             if (vel>0.d0) then  ! infiltration
                vel = min(vel, &
                          s%hh(i,j)*s%dsz(i,j)/par%dt/dz, &
                          max(s%zb(i+1,j)-s%gwlevel(i+1,j),0.d0)*par%por*s%dsz(i+1,j)/par%dt/dz)
             else
                vel = max(vel, &
                          -(s%gwlevel(i+1,j)-s%gwbottom(i+1,j))*s%dsz(i+1,j)*par%por/par%dt/dz)
             endif                                           
             infilhorgw(i+1,j) = infilhorgw(i+1,j)+vel*scaleareab
             infilhorsw(i,j) = infilhorsw(i,j)+vel*scaleareas
          elseif (dz<0.d0 .and. s%wetz(i+1,j)==1) then  ! jump down
             dz = -dz
             hsurf = s%zs(i+1,j)
             hgw = s%gwhead(i,j)
             dx = s%dsz(i,j)
             dhead = hsurf-hgw         ! positive gradient if hsurf>hgw
                                       ! negative velocitity if infiltration
             vel = dhead/dx*Kx         ! convert negative down velocity into positive if infiltration,
                                       ! negative if exfiltration 
             ! horizontal infiltration occurs along small side section, so scale the horizontal
             ! infiltration area by the vertical area of the cell used in the mass balance equation
             ! and infiltration velocities
             ! 
             ! Scale area for bed cell and surface water cell
             scaleareab = dz/s%dsz(i,j)
             scaleareas = dz/s%dsz(i+1,j)
             ! Limit  flow to that available in space in ground and 
             ! in surface water volume
             if (vel>0.d0) then  ! infiltration
                vel = min(vel, &
                          s%hh(i+1,j)*s%dsz(i+1,j)/par%dt/dz, &
                          max(s%zb(i,j)-s%gwlevel(i,j),0.d0)*par%por*s%dsz(i,j)/par%dt/dz)
             else
                vel = max(vel, &
                          -(s%gwlevel(i,j)-s%gwbottom(i,j))*s%dsz(i,j)*par%por/par%dt/dz)
             endif                                           
             infilhorgw(i,j) = infilhorgw(i,j)+vel*scaleareab
             infilhorsw(i+1,j) = infilhorsw(i+1,j)+vel*scaleareas
          else                ! equal bed height
             ! nothing
          endif
       enddo
    enddo
    !
    !
    ! loop over all v-interfaces
    do i=1,s%nx+1
        do j=1,s%ny
          dz = s%zb(i,j+1)-s%zb(i,j)
          if (dz>0.d0 .and. s%wetz(i,j)==1) then ! jump up
             hsurf = s%zs(i,j)
             hgw = s%gwhead(i,j+1)
             dy = s%dnz(i,j+1)
             dhead = hsurf-hgw         ! positive gradient if hsurf>hgw
                                       ! negative velocitity if infiltration
             vel = dhead/dy*Ky         ! convert negative down velocity into positive if infiltration,
                                       ! negative if exfiltration 
             ! horizontal infiltration occurs along small side section, so scale the horizontal
             ! infiltration area by the vertical area of the cell used in the mass balance equation
             ! and infiltration velocities
             ! 
             ! Scale area for bed cell and surface water cell
             scaleareab = dz/s%dnz(i,j+1)
             scaleareas = dz/s%dnz(i,j)
             ! Limit  flow to that available in space in ground and 
             ! in surface water volume
             if (vel>0.d0) then  ! infiltration
                vel = min(vel, &
                          s%hh(i,j)*s%dnz(i,j)/par%dt/dz, &
                          max(s%zb(i,j+1)-s%gwlevel(i,j+1),0.d0)*par%por*s%dnz(i,j+1)/par%dt/dz)
             else
                vel = max(vel, &
                          -(s%gwlevel(i,j+1)-s%gwbottom(i,j+1))*s%dnz(i,j+1)*par%por/par%dt/dz)
             endif                                           
             infilhorgw(i,j+1) = infilhorgw(i,j+1)+vel*scaleareab
             infilhorsw(i,j) = infilhorsw(i,j)+vel*scaleareas
          elseif (dz<0.d0 .and. s%wetz(i,j+1)==1) then  ! jump down
             dz = -dz
             hsurf = s%zs(i,j+1)
             hgw = s%gwhead(i,j)
             dy = s%dsz(i,j)
             dhead = hsurf-hgw         ! positive gradient if hsurf>hgw
                                       ! negative velocitity if infiltration
             vel = dhead/dy*Ky         ! convert negative down velocity into positive if infiltration,
                                       ! negative if exfiltration 
             ! horizontal infiltration occurs along small side section, so scale the horizontal
             ! infiltration area by the vertical area of the cell used in the mass balance equation
             ! and infiltration velocities
             ! 
             ! Scale area for bed cell and surface water cell
             scaleareab = dz/s%dnz(i,j)
             scaleareas = dz/s%dnz(i,j+1)
             ! Limit  flow to that available in space in ground and 
             ! in surface water volume
             if (vel>0.d0) then  ! infiltration
                vel = min(vel, &
                          s%hh(i,j+1)*s%dnz(i,j+1)/par%dt/dz, &
                          max(s%zb(i,j)-s%gwlevel(i,j),0.d0)*par%por*s%dnz(i,j)/par%dt/dz)
             else
                vel = max(vel, &
                          -(s%gwlevel(i,j)-s%gwbottom(i,j))*s%dnz(i,j)*par%por/par%dt/dz)
             endif                                           
             infilhorgw(i,j) = infilhorgw(i,j)+vel*scaleareab
             infilhorsw(i,j+1) = infilhorsw(i,j+1)+vel*scaleareas
          else                ! equal bed height
             ! nothing
          endif
       enddo
    enddo

  end subroutine gw_horizontal_infil_exfil

end module groundwaterflow
