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

real*8,dimension(:,:),allocatable,save,private  :: kru,krv,kru50,krv50,kru90,krv90
real*8,dimension(:,:),allocatable,save,private  :: urms_upd,u2_upd
real*8,dimension(:,:),allocatable,save,private  :: facbl,blphi,infilb,Ubed,Ventilation
real*8,dimension(:,:),allocatable,save,private  :: ueuf,uevf,vevf,veuf
real*8,dimension(:,:),allocatable,save,private  :: ueuold,uevold,vevold,veuold
real*8,dimension(:,:),allocatable,save,private  :: dudtsmooth,dvdtsmooth
real*8,dimension(:,:),allocatable,save,private  :: shieldsu,shieldsv
real*8,dimension(2),save,private                :: dtold
real*8,save,private                             :: delta,rhogdelta


public bedroughness_init
public bedroughness_update

contains

subroutine bedroughness_init(s,par)
   use params
   use spaceparams
   use paramsconst
   
   IMPLICIT NONE

   type(parameters),intent(in)                 :: par
   type(spacepars),target                      :: s
   integer                                     :: i,j

   !include 's.ind'
   !include 's.inp'
      
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
   
   IMPLICIT NONE

   type(parameters),intent(in)                 :: par
   type(spacepars),target                      :: s
   integer                                     :: i,j
   real*8                                      :: ubed,dudtbed,hda,klocal
   
   !include 's.ind'
   !include 's.inp'
   
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
   
   s%cfu = min(s%cfu,1.d0)
   s%cfv = min(s%cfv,1.d0)
   
   
end subroutine bedroughness_update 
   


end module bedroughness_module
