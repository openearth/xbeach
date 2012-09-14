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
    use logging_module

    IMPLICIT NONE

    type(spacepars),target                  :: s
    type(parameters), intent(in)            :: par

    integer                             :: i
    integer                             :: nwind,io
    real*8,dimension(3)                 :: temp


    include 's.ind'
    include 's.inp'

    io    = 0
    nwind = 0

    if (xmaster) then
       allocate(s%windsu(s%nx+1,s%ny+1))
       allocate(s%windnv(s%nx+1,s%ny+1))
       s%windsu=0.d0
       s%windnv=0.d0
       if (par%windfile==' ') then   ! Stationary wind

          s%windlen=1
          allocate(s%windinpt(1))
          allocate(s%windvelts (1))
          allocate(s%winddirts (1))
          allocate(s%windxts(1))
          allocate(s%windyts(1))

          s%windinpt=par%tstop
          s%windvelts=par%windv
          s%winddirts=(270.d0-par%windth-s%alfa)*par%px/180.d0
          ! Alfa is East to X still ?? 
          s%windxts = dcos(s%winddirts)*s%windvelts
          s%windyts = dsin(s%winddirts)*s%windvelts
          ! If stationary wind then we will calculate now what the wind velocity is everywhere now
          s%windsu = s%windxts(1)*dcos(alfau) + s%windyts(1)*dsin(s%alfau)
          s%windnv = s%windyts(1)*dcos(s%alfav-0.5d0*par%px) - s%windxts(1)*dsin(s%alfav-0.5d0*par%px)
       else                 ! Non-stationary wind
          call writelog('ls','','readwind: reading wind time series from ',trim(par%windfile),' ...')
          open(31,file=par%windfile)
          do while (io==0)
             nwind=nwind+1
             read(31,*,IOSTAT=io) temp
          enddo
          rewind(31)
          s%windlen=nwind-1
          allocate(s%windinpt(s%windlen))
          allocate(s%windvelts(s%windlen))
          allocate(s%winddirts(s%windlen))
          allocate(s%windxts(s%windlen))
          allocate(s%windyts(s%windlen))
          do i=1,s%windlen
             read(31,*,IOSTAT=io) s%windinpt(i),s%windvelts(i),s%winddirts(i)
             if (io .ne. 0) then
                call report_file_read_error(par%windfile)
             endif
          enddo
          ! to cartesian radians
          s%winddirts=(270.d0-s%winddirts-s%alfa)*par%px/180.d0
          s%windxts = dcos(s%winddirts)*s%windvelts
          s%windyts = dsin(s%winddirts)*s%windvelts
          close(31)
          if (par%morfacopt==1) s%windinpt = s%windinpt / max(par%morfac,1.d0)
          if (s%windinpt(s%windlen)<par%tstop) then
             call writelog('els','','Wind condition time series too short. Stopping calculation')
             call halt_program
          endif
       endif
    endif


  end subroutine readwind

end module readwind_module
