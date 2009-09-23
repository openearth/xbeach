module readwind_module
contains
subroutine readwind (s,par)
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

character*80                        :: fname
integer                             :: i
integer                             :: nwind,io
real*8,dimension(3)						:: temp


include 's.ind'
include 's.inp'

io    = 0
nwind = 0

if (xmaster) then

  par%rhoa    = readkey_dbl ('params.txt','rhoa',   1.25d0,     1.0d0,   2.0d0,bcast=.false.)
  par%Cd      = readkey_dbl ('params.txt','Cd',    0.002d0,  0.0001d0,  0.01d0,bcast=.false.)

  call readkey('params.txt','windfile',fname)

  if (fname=='') then   ! Stationary wind
      par%windv   = readkey_dbl ('params.txt','windv',   0.0d0,     0.0d0, 200.0d0,bcast=.false.)
      par%windth  = readkey_dbl ('params.txt','windth', 270.0d0,  -360.0d0, 360.0d0,bcast=.false.)

      par%windlen=2
      allocate(s%windinpt(par%windlen))
      allocate(s%windvel (par%windlen))
      allocate(s%winddir (par%windlen))

		s%windinpt(1)=0
      s%windinpt(2)=par%tstop
		s%windvel=par%windv
		s%winddir=(270.d0-par%windth-s%alfa)*par%px/180.d0
	else                 ! Non-stationary wind
     write(*,*)'readwind: reading wind time series from ',fname,' ...'
     open(31,file=fname)
     do while (io==0)
	    nwind=nwind+1
       read(31,*,IOSTAT=io) temp
     enddo
     rewind(31)
     par%windlen=nwind-1
     allocate(s%windinpt(par%windlen))
     allocate(s%windvel (par%windlen))
     allocate(s%winddir (par%windlen))
     do i=1,par%windlen
       read(31,*,IOSTAT=io) s%windinpt(i),s%windvel(i),s%winddir(i)
     enddo
     close(31)
	  s%windinpt = s%windinpt / max(par%morfac,1.d0)
	  if (s%windinpt(par%windlen)<par%tstop) then
	      write(*,*)'Error !!!! Wind condition time series too short. Stopping calculation !!!'
         call halt_program
	  endif
	endif
endif


end subroutine readwind

end module readwind_module
