module wave_timestep_module
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

contains
  subroutine wave(s,par)
    
    use params
    use spaceparams
    use wave_stationary_module
    use wave_instationary_module
    
    implicit none
    
    type(spacepars)                                 :: s
    type(parameters)                                :: par
  
    if (trim(par%instat) == 'stat' .or. trim(par%instat) == 'stat_table') then

#ifdef USEMPI
        call wave_instationary(s,par)
        if ((abs(mod(par%t,par%wavint))<0.000001d0) .or. par%newstatbc==1) then
            par%newstatbc   = 0
        endif
#else
        if ((abs(mod(par%t,par%wavint))<0.000001d0) .or. par%newstatbc==1) then
            call wave_stationary(s,par)
            par%newstatbc   = 0
        endif
#endif

    else
        par%newstatbc       = 0
        call wave_instationary(s,par)
    endif
    
  end subroutine

end module wave_timestep_module
