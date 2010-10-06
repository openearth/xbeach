module output_module
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef USENETCDF
use ncoutput_module
#endif
use params
use spaceparams
use timestep_module
use logging_module
use fortoutput_module
! IFDEF used in case netcdf support is not compiled, f.i. Windows (non-Cygwin)


contains
subroutine output_init(sglobal, slocal, par, tpar)
  type(spacepars), target, intent(in)  :: sglobal
  type(spacepars), target, intent(in)  :: slocal
  type(parameters), intent(in)         :: par
  type(timepars), intent(in)           :: tpar

  
! initialize the correct output module (clean this up?, move to another module?)
if (par%outputformat=='fortran') then
   ! only fortran
   call writelog('ls', '', 'Fortran outputformat')
   call var_output_init(sglobal,slocal,par,tpar)
elseif (par%outputformat=='netcdf') then
   ! only netcdf, stop if it's not build
   call writelog('ls', '', 'NetCDF outputformat')
#ifdef USENETCDF
   call ncoutput_init(sglobal,slocal,par,tpar)
#else
   call writelog('lse', '', 'This xbeach executable has no netcdf support. Rebuild with netcdf or outputformat=fortran')
   call halt_program
#endif
elseif (par%outputformat=='debug') then
   call writelog('ls', '', 'Debug outputformat, writing both netcdf and fortran output')
#ifdef USENETCDF
   call ncoutput_init(sglobal,slocal,par,tpar)
#endif
   call var_output_init(sglobal,slocal,par,tpar)
endif

end subroutine output_init
subroutine output(s,sglobal,par,tpar)
   
   implicit none
   
   type(spacepars)                     :: s,sglobal
   type(parameters)                    :: par
   type(timepars)                      :: tpar
        
   ! update output times
   call outputtimes_update(par, tpar)
   ! update log
   call log_progress(par)
   ! Output
   if (par%outputformat=='fortran') then
      call var_output(sglobal,s,par,tpar)
   elseif (par%outputformat=='netcdf') then
#ifdef USENETCDF
      call ncoutput(sglobal,s,par, tpar)
#endif
   elseif (par%outputformat=='debug') then
#ifdef USENETCDF
      call ncoutput(sglobal,s,par, tpar)
#endif
      call var_output(sglobal,s,par,tpar)
   endif
   
end subroutine output

subroutine log_progress(par)

   type(parameters)                    :: par
   logical,save                        :: firsttime = .true.
   integer,save                        :: day
   real*8,save                         :: tprev,percprev
   real*8                              :: tnow,percnow,tpredicted
   integer,dimension(8)                :: datetime
   
   
    if (firsttime) then 
       day=0
       call date_and_time(VALUES=datetime)
       tprev = day*24.d0*3600.d0+datetime(5)*3600.d0+60.d0*datetime(6)+1.d0*datetime(7)+0.001d0*datetime(8)
       percprev = 0.d0
       firsttime = .false.
    else    
       if (par%timings .ne. 0) then
          
          call date_and_time(VALUES=datetime)
          ! Current time in seconds
          tnow=day*24.d0*3600.d0+datetime(5)*3600.d0+60.d0*datetime(6)+1.d0*datetime(7)+0.001d0*datetime(8)
          ! Is it time to update the log again?
          if (tnow>=tprev+5.d0) then
             percnow = 100.d0*par%t/par%tstop
             ! Predict time based on percentage change rate in the last 5 seconds. 
             tpredicted = 100.d0*(1.d0-par%t/par%tstop)/(max(percnow-percprev,0.01d0)/(tnow-tprev))
             ! Percentage complete:
             call writelog('ls','(a,f5.1,a)','Simulation ',percnow,' percent complete')
             if (tpredicted>=3600) then 
                call writelog('ls','(a,I0,a,I2,a)','Time remaining ',&
                     floor(tpredicted/3600.0d0),' hours and ',&
                     nint((tpredicted-3600.0d0*floor(tpredicted/3600.0d0))/60.0d0),&
                     ' minutes')
             elseif (tpredicted>=600) then
                call writelog('ls','(a,I2,a)','Time remaining ',&
                     floor(tpredicted/60.0d0),' minutes')
             elseif (tpredicted>=60) then
                call writelog('ls','(a,I2,a,I2,a)','Time remaining ',&
                     floor(tpredicted/60.0d0),' minutes and ',&
                     nint((tpredicted-60.0d0*floor(tpredicted/60.0d0))),' seconds')
             else
                call writelog('ls','(a,I2,a)','Time remaining ',nint(tpredicted),' seconds')
             endif
             tprev=tnow
             percprev=percnow
          elseif (tnow<tprev-60.d0) then  ! It's probably the next day 
             day=day+1
          endif
       endif ! timings on
    endif ! firsttime logical
end subroutine log_progress
   
end module output_module

