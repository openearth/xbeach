module roelvink_module

    interface janssen_battjes
        module procedure janssen_battjes_1D
        module procedure janssen_battjes_2D
    end interface janssen_battjes
    
contains
subroutine roelvink1(E,hh,Trep,alpha,gamma,n,rho,g,delta,D,ntot,break)
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
IMPLICIT NONE

integer                         :: ntot
character(24)                   :: break
real*8,dimension(ntot)          :: E,hh,D
real*8                          :: fac,rho,g,delta,alpha,gamma,n,Trep


real*8,dimension(:),allocatable,save :: H,Qb,hroelvink

if (.not. allocated(H)) then
   allocate(H       (ntot))
   allocate(Qb      (ntot))
   allocate(hroelvink(ntot))
endif

! Dissipation acc. to Roelvink (1993)

fac=8.0d0/rho/g
H=sqrt(fac*E)
hroelvink=hh+delta*H  !add breakerdelay here, to correct for waterdepth, and cancel subroutine breakerdelay
					  !danojaaprobertap 10/2, also in baldock, roelvink

Qb=min(1-exp(-(H/gamma/hroelvink)**n),1.0d0) !gebruik hier Hb = (0.88d0/k)*tanh(gamma*kh/0.88d0)

! D=Qb*2.*alpha/Trep*E

! cjaap : two options:
if (trim(break)=='roelvink1') then
   D=Qb*2*alpha/Trep*E
elseif (trim(break)=='roelvink2') then
   ! D=Qb*2.*alpha/Trep*E*H/hroelvink; !breaker delay depth for dissipation in bore and fraction of breaking waves
   D=Qb*2*alpha/Trep*E*H/hh            !breaker delay depth for fraction breaking waves and actual water depth for disipation in bore
end if


end subroutine roelvink1

subroutine roelvink(par,s,km)
use params
use spaceparams

IMPLICIT NONE

type(spacepars)                 :: s
type(parameters)                :: par

real*8                          :: fac
integer                         :: i,j
real*8,dimension(s%nx+1,s%ny+1) :: km
real*8,dimension(:,:),allocatable,save :: H,hroelvink,arg,kmr

if (.not. allocated(H)) then
   allocate(H       (s%nx+1,s%ny+1))
   allocate(arg     (s%nx+1,s%ny+1))
   allocate(kmr     (s%nx+1,s%ny+1))
   allocate(hroelvink(s%nx+1,s%ny+1))
endif

! Dissipation acc. to Roelvink (1993)

fac=8.0d0/par%rho/par%g
H=sqrt(fac*s%E)

hroelvink=s%hh+par%delta*s%H

kmr=km
where(km<0.01d0)
   kmr=0.01d0
elsewhere(km>100.d0)
   kmr=100.d0
endwhere

if (trim(par%break)/='roelvink_daly') then
   if (par%wci==1) then
      arg = -(H/(par%gamma*tanh(min(max(km,0.01d0),100.d0)*hroelvink)/min(max(km,0.01d0),100.d0)))**par%n
      s%Qb = min(1.d0-exp(max(arg,-100.d0)),1.d0)
   else
      s%Qb=min(1-exp(-(H/par%gamma/hroelvink)**par%n),1.0d0)
   endif
else
   do j=1,s%ny+1
      do i=1,s%nx+1
         if (H(i,j)>par%gamma *hroelvink(i,j)) s%Qb(i,j)=1.d0
         if (H(i,j)<par%gamma2*hroelvink(i,j)) s%Qb(i,j)=0.d0
      enddo
   enddo
   s%Qb=max(s%Qb,0.d0)
   
endif

! cjaap : two options:
if (trim(par%break)=='roelvink1') then
   if (par%wci==1) then
      s%D=s%Qb*2.d0*par%alpha*s%sigm*s%E/2.d0/par%px;
!     s%D=Qb*par%alpha*km*s%E*H/par%px*sqrt(par%g*km/tanh(max(km,0.001)*hroelvink));  ! according to CK2002
   else
      s%D=s%Qb*2*par%alpha/par%Trep*s%E
   endif
elseif (trim(par%break)=='roelvink2' .or. trim(par%break)=='roelvink_daly') then
   ! s%D=s%Qb*2.*par%alpha/par%Trep*s%E*H/hroelvink; !breaker delay depth for dissipation in bore and fraction of breaking waves
   s%D=s%Qb*2*par%alpha/par%Trep*s%E*H/s%hh        !breaker delay depth for fraction breaking waves and actual water depth for disipation in bore
end if

end subroutine roelvink

subroutine baldock1(E,hh,k,Trep,alpha,gamma,rho,g,delta,D,ntot)
  IMPLICIT NONE

  integer                         :: ntot
  real*8,dimension(ntot)          :: E,hh,D,k,arg
  real*8                          :: fac,rho,g,delta,alpha,gamma,Trep


  real*8,dimension(:),allocatable,save :: H,Hb,Qb,kh,tkh,hbaldock
  
  if (.not. allocated(H)) then
     allocate(H       (ntot))
     allocate(Hb      (ntot))
     allocate(Qb      (ntot))
     allocate(kh      (ntot))
     allocate(tkh     (ntot))
     allocate(hbaldock(ntot))
  endif

  ! Dissipation acc. to Baldock et al. 1998

  fac=8.0d0/rho/g
  H=sqrt(fac*E)

  hbaldock=hh+delta*H

  kh=k*hbaldock
  ! tkh=tanh(kh)   ! tkh not used

  ! Wave dissipation acc. to Baldock et al. 1998
  Hb = (0.88d0/k)*tanh(gamma*kh/0.88d0)
  arg = -(Hb/max(H,0.00001d0))**2
  Qb = exp(max(arg,-100.d0))
  D = 0.25d0*alpha*Qb*rho*(1.0d0/Trep)*g*(Hb**2+H**2)

end subroutine baldock1

subroutine baldock(par,s,km)
  use params
  use spaceparams

  IMPLICIT NONE

  type(spacepars)                 :: s
  type(parameters)                :: par

  real*8                          :: fac
  integer                         :: alpha
  real*8,dimension(s%nx+1,s%ny+1) :: km

  real*8,dimension(:,:),allocatable,save :: H,Hb,kh,hbaldock,arg
  if (.not. allocated(H)) then
     allocate(H       (s%nx+1,s%ny+1))
     allocate(Hb      (s%nx+1,s%ny+1))
     allocate(kh      (s%nx+1,s%ny+1))
     allocate(arg     (s%nx+1,s%ny+1))
     allocate(hbaldock(s%nx+1,s%ny+1))
  endif

  ! Dissipation acc. to Baldock et al. 1998

  fac=8.d0/par%rho/par%g
  alpha=1
  hbaldock=s%hh+par%delta*s%H

  H=sqrt(fac*s%E)

  kh=s%k*hbaldock
  !tkh=tanh(kh)     ! tkh not used

  ! Wave dissipation acc. to Baldock et al. 1998
  !arg = -(Hb/max(H,0.00001d0))**2
  if (par%wci==1) then
     Hb = (0.88d0/km)*tanh(par%gamma*kh/0.88d0)
     !s%Qb = exp(max(arg,-100.d0))                                          ! bas: seems not right since arg is computed before Hb
     s%Qb = exp(-(Hb/max(H,0.00001d0))**2)
     s%D = 0.25d0*alpha*s%Qb*par%rho*s%sigm*par%g*(Hb**2+H**2)/2.d0/par%px
  else
     Hb = (0.88d0/s%k)*tanh(par%gamma*kh/0.88d0)
     s%Qb = exp(-(Hb/max(H,0.00001d0))**2)
     s%D = 0.25d0*alpha*s%Qb*par%rho*(1.d0/par%Trep)*par%g*(Hb**2+H**2)
  endif

end subroutine baldock

subroutine janssen_battjes_1D(par,s,km,i)
    
    use params
    use spaceparams

    implicit none

    type(spacepars)                                 :: s
    type(parameters)                                :: par
    
    integer                                         :: i
    real*8                                          :: B
    
    real*8, dimension(s%ny+1)                       :: km,kh,f,k,H,Hb,R

    ! Dissipation according to Janssen and Battjes 2007
    
    B   = par%alpha
    
    if (par%wci==1) then
        f = s%sigm(i,:) / 2.d0 / par%px
        k = km
    else
        f = 1.d0 / par%Trep
        k = s%k(i,:)
    endif
    
    kh  = s%k(i,:) * (s%hh(i,:) + par%delta*s%H(i,:))
    
    H   = sqrt(8.d0/par%rho/par%g*s%E(i,:))
    Hb  = tanh(par%gamma*kh/0.88d0)*(0.88d0/k)
    R   = Hb/max(H,0.00001d0)
    
    s%Qb(i,:)   = 1 + 4/(3*sqrt(par%px)) * (R**3 + 3/2*R) * exp(-R**2) - derf(R)
    s%D (i,:)   = 3*sqrt(par%px)/16      * B * f * par%rho * par%g * H**3/s%hh(i,:) * s%Qb(i,:)
    
end subroutine janssen_battjes_1D

subroutine janssen_battjes_2D(par,s,km)

    use params
    use spaceparams

    implicit none

    type(spacepars)                                 :: s
    type(parameters)                                :: par
    
    integer                                         :: i
    
    real*8, dimension(s%nx+1,s%ny+1)                :: km
    
    do i = 1,s%nx+1
        call janssen_battjes_1D(par,s,km(i,:),i)
    enddo
    
end subroutine janssen_battjes_2D

end module roelvink_module
