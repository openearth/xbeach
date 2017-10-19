module roelvink_module
   implicit none
   save
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

   interface roelvink
      module procedure roelvink_1D
      module procedure roelvink_2D
   end interface roelvink

   interface baldock
      module procedure baldock_1D
      module procedure baldock_2D
   end interface baldock

   interface janssen_battjes
      module procedure janssen_battjes_1D
      module procedure janssen_battjes_2D
   end interface janssen_battjes

contains

   subroutine roelvink_1D(par,s,i)

      use params
      use spaceparams
      use paramsconst

      implicit none

      type(spacepars)                                 :: s
      type(parameters)                                :: par
      integer,intent(in)                              :: i

      integer                                         :: j

      real*8, dimension(s%ny+1)                       ::kmr,hr,H,arg

      ! Dissipation according to Roelvink (1993)

      H   = s%H(i,:)
      hr  = s%hhw(i,:)

      kmr = min(max(s%k(i,:), 0.01d0), 100.d0)

      if (par%break /= BREAK_ROELVINK_DALY) then
         if (par%wci==1) then
            arg = -( H / (par%gamma*tanh(kmr*hr)/kmr))**par%n
         else
            arg = -( H / (par%gamma*hr              ))**par%n
         endif

         s%Qb(i,:) = min(1.d0 - exp(max(arg,-100.d0)), 1.d0)
      else
         do j=1,s%ny+1
            if (H(j) > par%gamma *hr(j)) s%Qb(i,j) = 1.d0
            if (H(j) < par%gamma2*hr(j)) s%Qb(i,j) = 0.d0
         enddo

         s%Qb(i,:) = max(s%Qb(i,:), 0.d0)
      endif

      s%D(i,:) = s%Qb(i,:) * 2.d0 * par%alpha * s%E(i,:)

      if (par%break == BREAK_ROELVINK1) then
         if (par%wci==1) then
            s%D(i,:) = s%D(i,:) * s%sigm(i,:)/2.d0/par%px;
         else
            s%D(i,:) = s%D(i,:) / par%Trep
         endif
      elseif (par%break == BREAK_ROELVINK2 .or. par%break == BREAK_ROELVINK_DALY) then
         ! Jaap: also wci switch for roelvink2
         if (par%wci==1) then
            s%D(i,:) = s%D(i,:) * s%sigm(i,:)/2.d0/par%px * H/s%hh(i,:);
         else
            s%D(i,:) = s%D(i,:) / par%Trep * H/s%hh(i,:)
         endif
      end if

   end subroutine roelvink_1D

   subroutine roelvink_2D(par,s)

      use params
      use spaceparams

      implicit none

      type(spacepars)                                 :: s
      type(parameters)                                :: par

      integer                                         :: i

      do i = 1,s%nx+1
         call roelvink_1D(par,s,i)
      enddo

   end subroutine roelvink_2D

   subroutine baldock_1D(par,s,i)

      use params
      use spaceparams

      implicit none

      type(spacepars)                 :: s
      type(parameters)                :: par

      integer,intent(in)              :: i

      real*8, dimension(s%ny+1)       :: kh,f,k,H,Hb,R,gamma

      ! Dissipation according to Baldock et al. (1998)

      k = s%k(i,:)
      if (par%wci==1) then
         f = s%sigm(i,:) / 2.d0 / par%px
      else
         f = 1.d0 / par%Trep
      endif

      kh  = s%k(i,:) * s%hhw(i,:)

      if (par%wci == 1) then
         gamma = 0.76d0*kh + 0.29d0 !Jaap: spatial varying gamma according to Ruessink et al., 1998
      else
         gamma = par%gamma
      endif

      H   = s%H(i,:)
      Hb  = tanh(gamma*kh/0.88d0)*(0.88d0/max(k,1d-10))
      R   = Hb/max(H,0.00001d0)

      s%Qb(i,:)   = exp(-R**2)
      s%D (i,:)   = 0.25d0 * par%alpha * f * par%rho * par%g * (Hb**2+H**2) * s%Qb(i,:)

   end subroutine baldock_1D

   subroutine baldock_2D(par,s)

      use params
      use spaceparams

      implicit none

      type(spacepars)                                 :: s
      type(parameters)                                :: par

      integer                                         :: i

      do i = 1,s%nx+1
         call baldock_1D(par,s,i)
      enddo

   end subroutine baldock_2D

   subroutine janssen_battjes_1D(par,s,i)

      use params
      use spaceparams
      use math_tools, only: xerf

      implicit none

      type(spacepars)                                 :: s
      type(parameters)                                :: par

      integer,intent(in)                              :: i
      real*8                                          :: B

      real*8, dimension(s%ny+1)                       :: kh,f,k,H,Hb,R

      ! Dissipation according to Janssen and Battjes (2007)

      B   = par%alpha

      k = s%k(i,:)
      if (par%wci==1) then
         f = s%sigm(i,:) / 2.d0 / par%px
      else
         f = 1.d0 / par%Trep
      endif

      kh  = s%k(i,:) * s%hhw(i,:)

      H   = s%H(i,:)
      Hb  = tanh(par%gamma*kh/0.88d0)*(0.88d0/k)
      R   = Hb/max(H,0.00001d0)

      s%Qb(i,:)   = 1 + 4/(3*sqrt(par%px)) * (R**3 + 3/2*R) * exp(-R**2) - xerf(R)
      s%D (i,:)   = 3*sqrt(par%px)/16      * B * f * par%rho * par%g * H**3/s%hh(i,:) * s%Qb(i,:)

   end subroutine janssen_battjes_1D

   subroutine janssen_battjes_2D(par,s)

      use params
      use spaceparams

      implicit none

      type(spacepars)                                 :: s
      type(parameters)                                :: par

      integer                                         :: i

      do i = 1,s%nx+1
         call janssen_battjes_1D(par,s,i)
      enddo

   end subroutine janssen_battjes_2D

end module roelvink_module
