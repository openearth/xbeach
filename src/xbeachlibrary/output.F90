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
  use paramsconst
  ! IFDEF used in case netcdf support is not compiled, f.i. Windows (non-Cygwin)


contains
  subroutine output_init(sglobal, slocal, par, tpar)
    type(spacepars), target, intent(in)  :: sglobal
    type(spacepars), target, intent(in)  :: slocal
    type(parameters), intent(in)         :: par
    type(timepars), intent(in)           :: tpar


    ! initialize the correct output module (clean this up?, move to another module?)
    if (par%outputformat==OUTPUTFORMAT_FORTRAN) then
       ! only fortran
       call writelog('ls', '', 'Fortran outputformat')
       call var_output_init(sglobal,slocal,par,tpar)
    elseif (par%outputformat==OUTPUTFORMAT_NETCDF) then
       ! only netcdf, stop if it's not build
       call writelog('ls', '', 'NetCDF outputformat')
#ifdef USENETCDF
       call ncoutput_init(sglobal,slocal,par,tpar)
#else
       call writelog('lse', '', 'This xbeach executable has no netcdf support. Rebuild with netcdf or outputformat=fortran')
       call halt_program
#endif
    elseif (par%outputformat==OUTPUTFORMAT_DEBUG) then
       call writelog('ls', '', 'Debug outputformat, writing both netcdf and fortran output')
#ifdef USENETCDF
       call ncoutput_init(sglobal,slocal,par,tpar)
#endif
       call var_output_init(sglobal,slocal,par,tpar)
    endif

  end subroutine output_init

  subroutine output(sglobal,s,par,tpar,update)

    use means_module

    implicit none

    type(spacepars)                     :: s,sglobal
    type(parameters)                    :: par
    type(timepars)                      :: tpar

    logical, optional                   :: update
    logical                             :: lupdate

    if (present(update)) then
       lupdate = update
    else
       lupdate = .true.
    endif

    ! update output times
    if (lupdate) call outputtimes_update(par, tpar)

    ! update log
    call log_progress(par)

    ! update meanvars in current averaging period with current timestep
    if (par%nmeanvar/=0) then
       if (par%t>tpar%tpm(1) .and. par%t<=tpar%tpm(size(tpar%tpm))) then
          call makeaverage(s,par)
       endif
    endif

    ! output
    if (par%outputformat==OUTPUTFORMAT_FORTRAN) then
       call var_output(sglobal,s,par,tpar)
    elseif (par%outputformat==OUTPUTFORMAT_NETCDF) then
#ifdef USENETCDF
       call ncoutput(sglobal,s,par, tpar)
#endif
    elseif (par%outputformat==OUTPUTFORMAT_DEBUG) then
#ifdef USENETCDF
       call ncoutput(sglobal,s,par, tpar)
#endif
       call var_output(sglobal,s,par,tpar)
    endif

    ! clear averages after output of means
    if (tpar%outputm .and. tpar%itm>1) then
       call clearaverage(par)
    endif

  end subroutine output

  subroutine output_error(s, sglobal, par, tpar)

    use logging_module
    use xmpi_module, only: halt_program

    implicit none

    type(spacepars)                     :: s,sglobal
    type(parameters)                    :: par
    type(timepars)                      :: tpar

    call output(s, sglobal, par, tpar, update=.false.)
    call writelog('lse','','An extra output timestep is created to inquire the last timestep')
    call writelog('lse','','    before an error occured')
    call halt_program

  end subroutine output_error

  subroutine log_progress(par)

    type(parameters)                    :: par
    logical,save                        :: firsttime = .true.
    integer,save                        :: day, ndt
    real*8,save                         :: tprev,percprev,sumdt,dtavg,t0
    real*8                              :: tnow,percnow,tpredicted,tpredicted2,tpredmean
    integer,dimension(8)                :: datetime


    if (firsttime) then 
       day=0
       call date_and_time(VALUES=datetime)
       tprev = day*24.d0*3600.d0+datetime(5)*3600.d0+60.d0*datetime(6)+1.d0*datetime(7)+0.001d0*datetime(8)
       percprev = 0.d0
       t0 = tprev
       firsttime = .false.
       ndt = 0
       sumdt = 0.d0
    else    
       if (par%timings .ne. 0) then
          ndt=ndt+1
          sumdt=sumdt+par%dt
          call date_and_time(VALUES=datetime)
          ! Current time in seconds
          tnow=day*24.d0*3600.d0+datetime(5)*3600.d0+60.d0*datetime(6)+1.d0*datetime(7)+0.001d0*datetime(8)
          ! Is it time to update the log again?
          if (tnow>=tprev+5.d0) then
             ! Percentage complete
             percnow = 100.d0*par%t/par%tstop
             ! Average time step in the last 5 seconds
             dtavg = sumdt/ndt
             call writelog('ls','(a,f5.1,a,f7.3,a)','Simulation ',percnow,' percent complete. Average dt ',sumdt/ndt,' seconds')
             ndt=0
             sumdt=0.d0
             ! Predict time based on percentage change rate in the last 5 seconds. 
             tpredicted = 100.d0*(1.d0-par%t/par%tstop)/(max(percnow-percprev,0.01d0)/(tnow-tprev))
             tpredicted2 = 100.d0*(1.d0-par%t/par%tstop)/(max(percnow,0.01d0)/(tnow-t0))
             tpredmean = (tpredicted+tpredicted2)/2.d0
             ! Percentage complete:
             if (tpredmean>=3600) then 
                call writelog('ls','(a,I3,a,I3,a)','Time remaining',&
                     floor(tpredmean/3600.0d0),' hours and ',&
                     nint((tpredmean-3600.0d0*floor(tpredmean/3600.0d0))/60.0d0),&
                     ' minutes')
             elseif (tpredmean>=600) then
                call writelog('ls','(a,I3,a)','Time remaining ',&
                     floor(tpredmean/60.0d0),' minutes')
             elseif (tpredmean>=60) then
                call writelog('ls','(a,I3,a,I3,a)','Time remaining ',&
                     floor(tpredmean/60.0d0),' minutes and ',&
                     nint((tpredmean-60.0d0*floor(tpredmean/60.0d0))),' seconds')
             else
                call writelog('ls','(a,I3,a)','Time remaining ',nint(tpredmean),' seconds')
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

