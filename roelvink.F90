module roelvink_module
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

integer                         :: ntot,break
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
if (break==1) then
   D=Qb*2*alpha/Trep*E
elseif (break==3) then
   ! D=Qb*2.*alpha/Trep*E*H/hroelvink; !breaker delay depth for dissipation in bore and fraction of breaking waves
   D=Qb*2*alpha/Trep*E*H/hh        !breaker delay depth for fraction breaking waves and actual water depth for disipation in bore
end if


end subroutine roelvink1

subroutine roelvink(par,s)
use params
use spaceparams

IMPLICIT NONE

type(spacepars)                 :: s
type(parameters)                :: par

real*8                          :: fac



real*8,dimension(:,:),allocatable,save :: H,Qb,hroelvink

if (.not. allocated(H)) then
   allocate(H       (s%nx+1,s%ny+1))
   allocate(Qb      (s%nx+1,s%ny+1))
   allocate(hroelvink(s%nx+1,s%ny+1))
endif

! Dissipation acc. to Roelvink (1993)

fac=8.0d0/par%rho/par%g
hroelvink=s%hh+par%delta*s%H

H=sqrt(fac*s%E)
Qb=min(1-exp(-(H/par%gamma/hroelvink)**par%n),1.0d0)

! cjaap : two options:
if (par%break==1) then
   s%D=Qb*2*par%alpha/par%Trep*s%E
elseif (par%break==3) then
   ! s%D=Qb*2.*par%alpha/par%Trep*s%E*H/hroelvink; !breaker delay depth for dissipation in bore and fraction of breaking waves
   s%D=Qb*2*par%alpha/par%Trep*s%E*H/s%hh        !breaker delay depth for fraction breaking waves and actual water depth for disipation in bore
end if

end subroutine roelvink

subroutine baldock1(E,hh,k,Trep,alpha,gamma,rho,g,delta,D,ntot)
  IMPLICIT NONE

  integer                         :: ntot
  real*8,dimension(ntot)          :: E,hh,D,k
  real*8                          :: fac,rho,g,delta,alpha,gamma,Trep


  real*8,dimension(:),allocatable,save :: H,Qb,Hb,kh,tkh,hbaldock
  
  if (.not. allocated(H)) then
     allocate(H       (ntot))
     allocate(Qb      (ntot))
     allocate(Hb      (ntot))
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
  Qb = exp(-(Hb/max(H,0.00001d0))**2)
  D = 0.25d0*alpha*Qb*rho*(1.0d0/Trep)*g*(Hb**2+H**2)

end subroutine baldock1

subroutine baldock(par,s)
  use params
  use spaceparams

  IMPLICIT NONE

  type(spacepars)                 :: s
  type(parameters)                :: par

  real*8                          :: fac
  integer                         :: alpha


  real*8,dimension(:,:),allocatable,save :: H,Qb,Hb,kh,tkh,hbaldock
  if (.not. allocated(H)) then
     allocate(H       (s%nx+1,s%ny+1))
     allocate(Qb      (s%nx+1,s%ny+1))
     allocate(Hb      (s%nx+1,s%ny+1))
     allocate(kh      (s%nx+1,s%ny+1))
     allocate(tkh     (s%nx+1,s%ny+1))
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
  Hb = (0.88d0/s%k)*tanh(par%gamma*kh/0.88d0)
  Qb = exp(-(Hb/max(H,0.00001d0))**2)
  s%D = 0.25d0*alpha*Qb*par%rho*(1.d0/par%Trep)*par%g*(Hb**2+H**2)

end subroutine baldock

end module roelvink_module
