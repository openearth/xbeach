module wetcells_module
contains
  subroutine compute_wetcells(s,par)
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

    IMPLICIT NONE

    type(spacepars),target                  :: s
    type(parameters),intent(in)             :: par

    integer                                 :: i
    integer                                 :: j
   
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
  end subroutine compute_wetcells
end module wetcells_module