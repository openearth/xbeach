module readtide_module
contains
subroutine readtide (s,par)
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
use readkey_module

IMPLICIT NONE

type(spacepars),target                  :: s
type(parameters)                        :: par

character*80                        :: fnamezs0
integer                             :: i
integer                             :: j
real*8,dimension(:,:),allocatable   :: tidedummy

include 's.ind'
include 's.inp'

allocate(tidedummy(par%tidelen,par%tideloc+1))
allocate(s%tideinpt(par%tidelen))
allocate(s%tideinpz(par%tidelen,par%tideloc))

if (xreader) then
  write(*,*) 'tidelen=', par%tidelen , par%tideloc
endif

call readkey('params.txt','zs0file',fnamezs0)

if (xreader) then
  open(31,file=fnamezs0)
    do i=1,par%tidelen
      if (xreader) then
        read(31,*)(tidedummy(i,j),j=1,par%tideloc+1)
      endif
    end do
  close(31)
endif

#ifdef USEMPI
call xmpi_bcast(tidedummy)
#endif
do i=1,par%tidelen
    s%tideinpt(i) = tidedummy(i,1)
end do

do j=1,par%tideloc
    do i=1,par%tidelen
        s%tideinpz(i,j) = tidedummy(i,j+1)
    end do
end do

deallocate(tidedummy)

end subroutine readtide

end module readtide_module
