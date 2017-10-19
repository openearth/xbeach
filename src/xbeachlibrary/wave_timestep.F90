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
   implicit none
   save

contains
   subroutine wave(s,par)

      use params
      use spaceparams
      use wave_stationary_module
      use wave_directions_module
      use wave_instationary_module
      use paramsconst
      use wave_functions_module, only: wave_dispersion,update_means_wave_flow

      implicit none

      type(spacepars)                                 :: s
      type(parameters)                                :: par
      
      ! set basic water depth for all wave calculations: note that depending on wci setting, other water depths may be used in 
      ! dispersion routine etc. 
      if (par%delta>0.d0) then
         s%hhw = max(s%hh+par%delta*s%H,par%eps)
      else
         s%hhw = max(s%hh,par%eps) ! hh can be less than eps after morphevolution?
      endif
      
      select case (par%wavemodel)
      case(WAVEMODEL_STATIONARY)
         ! call update_means_wave_flow not required: stationary always uses instantaneous water depth s%hhw
         if ((abs(mod(par%t,par%wavint))<0.001d0*par%dt) .or. s%newstatbc==1) then
            call wave_dispersion(s,par,0)  ! use instantaneous water depth (and velocity)
            call wave_stationary(s,par)
            s%newstatbc   = 0
         endif
      case(WAVEMODEL_SURFBEAT)
         if (par%single_dir==1) then
            ! always need to update depths for wave directions model, as well as in case of wci (dealt with in subroutine)
            call update_means_wave_flow(s,par)  
            !
            ! if necessary, update wave directions
            if ((abs(mod(par%t,par%wavint))<0.001d0*par%dt) .or. s%newstatbc==1 .or. par%t==par%dt) then
               call wave_dispersion(s,par,1)  ! use s%hhws water depth (and velocity)
               call wave_directions(s,par)
               s%newstatbc   = 0
            endif
            !
            s%newstatbc       = 0 ! not sure if this is needed every timestep, but no overhead to keep in ...
            !
            ! need to call dispersion again, because dispersion is different in timestep and direction computation
            if (par%wci==1) then
               call wave_dispersion(s,par,2) ! use s%hhwcins water depth (and velocity)
            else
               call wave_dispersion(s,par,0) ! use s%hhw water depth
            endif
            call wave_instationary(s,par)
         else
            s%newstatbc       = 0
            if (par%wci==1) then
               ! in this case we only need to update flow depth in case of wci
               call update_means_wave_flow(s,par)
               call wave_dispersion(s,par,2) ! use s%hhwcins water depth (and velocity)
            else
               call wave_dispersion(s,par,0) ! use s%hhw water depth
            endif
            call wave_instationary(s,par)
         endif
      end select
         

   end subroutine wave

end module wave_timestep_module
