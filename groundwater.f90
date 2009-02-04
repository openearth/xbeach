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

contains



subroutine gwinit(par,s)
use params
use xmpi_module
use spaceparams
use readkey_module

IMPLICIT NONE

type(parameters)                            :: par
type(spacepars)								:: s

character*80                                :: fname
real*8                                      :: aquiferbot,temp
integer										:: i,j


if (xmaster) then
	allocate (s%gwhead(s%nx+1,s%ny+1))
	allocate (s%gwlevel(s%nx+1,s%ny+1))
	allocate (s%gwheight(s%nx+1,s%ny+1))
	allocate (s%gwu(s%nx+1,s%ny+1))
	allocate (s%gwv(s%nx+1,s%ny+1))
	allocate (s%gww(s%nx+1,s%ny+1))
    allocate (s%gwbottom(s%nx+1,s%ny+1))
	allocate (s%dinfil(s%nx+1,s%ny+1))
	allocate (s%gw0back(s%ny+1))


    call readkey('params.txt','aquiferbotfile',fname)
	if (fname=='') then     ! Not a filename
	   temp = minval(s%zb)
	   aquiferbot = readkey_dbl('params.txt','aquiferbot',temp-3.d0,-100.d0,100.d0,bcast=.false.)
       s%gwbottom=aquiferbot
	else
	  open(31,file=fname)
      do j=1,s%ny+1
          read(31,*)(s%gwbottom(i,j),i=1,s%nx+1)
      end do
      close(31)
    endif

    call readkey('params.txt','gw0file',fname)
	if (fname=='') then     ! Not a filename
	   temp = readkey_dbl('params.txt','gw0',0.d0,-5.d0,5.d0,bcast=.false.)
       s%gwhead=temp
	else
	  open(31,file=fname)
      do j=1,s%ny+1
          read(31,*)(s%gwhead(i,j),i=1,s%nx+1)
      end do
      close(31)
    endif

	s%gw0back=s%gwhead(s%nx+1,:)
	s%gwlevel=min(s%zb,s%gwhead)
    s%gwlevel=max(s%gwlevel,s%gwbottom+par%eps)
	s%gwheight=s%gwlevel-s%gwbottom
	s%gwu=0.d0
	s%gwv=0.d0
	s%gww=0.d0
	s%dinfil=max(par%dwetlayer/3.d0,0.02)   ! Centroid of area influenced instantly by free surface level lies at dwetlayer/3
endif
end subroutine



subroutine gwbc(par,s)

use params
use xmpi_module
use spaceparams

IMPLICIT NONE

type(parameters)                            :: par
type(spacepars)								:: s


if(xmpi_istop) then
   s%gwhead(1,:)=s%zs(1,:)
elseif (xmpi_isbot) then
   if (par%tideloc==4 .or. (par%tideloc==2 .and. par%paulrevere==0)) then
       s%gwhead(s%nx+1,:)=s%zs(s%nx+1,:)
   else
      s%gwhead(s%nx+1,:)=s%gw0back
   endif
endif

s%gwhead(:,1)=s%gwhead(:,2)
s%gwhead(:,s%ny+1)=s%gwhead(:,s%ny)

#ifdef USEMPI
   call xmpi_shift(s%gwhead,':1')
   call xmpi_shift(s%gwhead,':n')
   call xmpi_shift(s%gwhead,'1:')
   call xmpi_shift(s%gwhead,'m:')
#endif

s%gwlevel(1,:)=min(s%gwhead(1,:),s%zb(1,:))
s%gwlevel(s%nx+1,:)=min(s%gwhead(s%nx+1,:),s%zb(s%nx+1,:))
s%gwlevel(:,1)=min(s%gwhead(:,1),s%zb(:,1))
s%gwlevel(:,s%ny+1)=min(s%gwhead(:,s%ny+1),s%zb(:,s%ny+1))

s%gwbottom=min(s%gwbottom,s%zb-par%eps)

end subroutine



subroutine gwflow(par,s)

use params
use xmpi_module
use spaceparams

IMPLICIT NONE

type(parameters)                            :: par
type(spacepars)								:: s

integer										:: i,j
real*8,dimension(:,:),allocatable           :: dheaddx,dheaddy,dleveldt,gwqx,gwqy,gwhu,gwhv
real*8,dimension(:,:),allocatable           :: c1,c2,r,fsh

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

where (s%wetz==1 .and. s%gwlevel>s%zb-par%dwetlayer)
    s%gwhead=s%gwlevel+(s%zs-s%gwlevel)*(1.d0-r)
elsewhere
    s%gwhead=s%gwlevel
endwhere

dheaddx=0.d0
dheaddy=0.d0

do j=1,s%ny+1
    do i=1,s%nx
        dheaddx(i,j)=(s%gwhead(i+1,j)-s%gwhead(i,j))/(s%xz(i+1)-s%xz(i))  
    end do
end do

do j=1,s%ny
    do i=1,s%nx+1
        dheaddy(i,j)=(s%gwhead(i,j+1)-s%gwhead(i,j))/(s%yz(j+1)-s%yz(j))    
    end do
end do

s%gwheight=s%gwlevel-s%gwbottom

! Determine intermediate aquifer depths
gwhu(1:s%nx,:)=0.5d0*(s%gwheight(1:s%nx,:)+s%gwheight(2:s%nx+1,:))
gwhv(:,1:s%ny)=0.5d0*(s%gwheight(:,1:s%ny)+s%gwheight(:,2:s%ny+1))

! Determine fluxes
s%gwu=-par%kx*dheaddx
s%gwv=-par%ky*dheaddy

! Limit for stability in case of very high kx values
!s%gwu(1:s%nx,1:s%ny)=min(s%gwu(1:s%nx,1:s%ny),0.5d0*(s%x(2:s%nx+1,1:s%ny)-s%x(1:s%nx,1:s%ny))/par%dt)
!s%gwv(1:s%nx,1:s%ny)=min(s%gwv(1:s%nx,1:s%ny),0.5d0*(s%y(1:s%nx,2:s%ny+1)-s%y(1:s%nx,1:s%ny))/par%dt)
	   
gwqx=s%gwu*gwhu
gwqy=s%gwv*gwhv

! Stop cells from drying up
where(s%gwlevel(2:s%nx,2:s%ny)<=s%gwbottom(2:s%nx,2:s%ny)+par%eps)
   gwqx(2:s%nx,2:s%ny)=min(gwqx(2:s%nx,2:s%ny),0.d0)
   gwqx(1:s%nx-1,2:s%ny)=max(gwqx(1:s%nx-1,2:s%ny),0.d0)
   gwqy(2:s%nx,2:s%ny)=min(gwqy(2:s%nx,2:s%ny),0.d0)
   gwqy(2:s%nx,1:s%ny-1)=max(gwqy(2:s%nx,1:s%ny-1),0.d0)
end where

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


s%gww=0.d0						! w defined positive from sea to groundwater in volumes of surface water (no pores).
s%gww=par%por*(&
              -( c1* (s%gwlevel-s%zb) / par%dt                                      )&
              +( c2* ((1.d0-r)*(s%zb-s%gwlevel) / par%dt + r*(par%kz*fsh/s%dinfil)) )&
!             +(case2*((1.d0-r)*max((s%zb-s%gwhead),0.d0)+r*(fsh*gw%kper)))&
!			  +(case3*fsh*gw%kper)&
			  )


! ensure that water extracted from surface layer is not more than available
where (s%gww*par%dt>s%hh)
	s%gww=s%hh/par%dt
endwhere

dleveldt=0.d0
! Mass balance
do j=2,s%ny
    do i=2,s%nx
	    dleveldt(i,j)=-1.d0*(gwqx(i,j)-gwqx(i-1,j))/(s%xu(i)-s%xu(i-1)) &
		              -1.d0*(gwqy(i,j)-gwqy(i,j-1))/(s%yv(j)-s%yv(j-1)) &
				      +1.d0*s%gww(i,j)/par%por
    enddo
enddo

s%gwlevel=s%gwlevel+dleveldt*par%dt

! Update quasi vertical model infiltration layer thickness
s%dinfil=s%dinfil+s%gww*par%dt/par%por






end subroutine


!subroutine gwoutput(par,s,it,gw)
!
!use params
!use xmpi_module
!use spaceparams
!
!IMPLICIT NONE
!
!type(parameters)                            :: par
!type(gwpars)                                :: gw
!type(spacepars)								:: s
!
!integer										:: i,j,reclen,it
!
!
!
!reclen=2*(s%nx+1)*(s%ny+1)
!
!if (it==1) then
!    open(1102,file='gwhead.dat',form='unformatted',access='direct',recl=reclen)
!	open(1103,file='gwu.dat',form='unformatted',access='direct',recl=reclen)
!	open(1104,file='gwv.dat',form='unformatted',access='direct',recl=reclen)
!	open(1105,file='gww.dat',form='unformatted',access='direct',recl=reclen)
!	open(1106,file='gwlevel.dat',form='unformatted',access='direct',recl=reclen)
!endif
!
!write(1102,rec=it)s%gwhead
!write(1103,rec=it)s%gwu
!write(1104,rec=it)s%gwv
!write(1105,rec=it)s%gww
!write(1106,rec=it)s%gwlevel
!
!if(par%t>=par%tstop) then
!    close(1102)
!    close(1103)
!	close(1104)
!	close(1105)
!	close(1106)
!endif
!
!end subroutine


end module




