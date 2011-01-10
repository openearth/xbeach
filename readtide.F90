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
use logging_module

IMPLICIT NONE

type(spacepars),target                  :: s
type(parameters)                        :: par

character*80                        :: fnamezs0
integer                             :: i
integer                             :: io,ntide
real*8,dimension(par%tideloc+1)     :: temp

include 's.ind'
include 's.inp'

! this must only work on master
if (xmaster) then
  if (par%tideloc .eq. 0) then
     s%tidelen = 2
     allocate(s%tideinpt(s%tidelen))
     allocate(s%tideinpz(s%tidelen,par%tideloc))
     return
  endif

  io = 0
  ntide = 0

  call readkey('params.txt','zs0file',fnamezs0)
  call writelog('ls','','readtide: reading tide time series from ',fnamezs0,' ...')
  open(31,file=fnamezs0)
  do while (io==0)
     ntide=ntide+1
     read(31,*,IOSTAT=io) temp
  enddo
  rewind(31)

  s%tidelen=ntide-1

  allocate(s%tideinpt(s%tidelen))
  allocate(s%tideinpz(s%tidelen,par%tideloc))
  s%tideinpz = 0.0d0
  do i=1,s%tidelen
     read(31,*)s%tideinpt(i),s%tideinpz(i,:)
  end do
  close(31)

  if (par%morfacopt==1) s%tideinpt = s%tideinpt / max(par%morfac,1.d0)
  if (s%tideinpt(s%tidelen)<par%tstop) then
     call writelog('els','','Tide condition time series too short. Stopping calculation')
     call halt_program
  endif

endif

end subroutine readtide

end module readtide_module
