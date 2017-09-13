module output_module

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
   use ncoutput_module
   use params
   use spaceparams
   use timestep_module
   use logging_module
   use paramsconst
   ! IFDEF used in case netcdf support is not compiled, f.i. Windows (non-Cygwin)
   implicit none
   save

contains
   !_________________________________________________________________________________

   subroutine output_init(sglobal, slocal, par, tpar)
      implicit none
      type(spacepars), target, intent(inout)  :: sglobal
      type(spacepars), target, intent(in)  :: slocal
      type(parameters), intent(in)         :: par
      type(timepars), intent(in)           :: tpar


      ! get xz and yz in sglobal
      ! and initialize sglobal%collected and sglobal%precollected

#ifdef USEMPI
      sglobal%collected = .false.

      call space_collect_mnem(sglobal,slocal,par,mnem_xz)

      call space_collect_mnem(sglobal,slocal,par,mnem_yz)

      if(xomaster) then
         sglobal%precollected = sglobal%collected
      endif
#endif

      ! initialize the correct output module (clean this up?, move to another module?)
      call points_output_init(sglobal,par)

      select case(par%outputformat)

       case(OUTPUTFORMAT_FORTRAN)
         ! only fortran
         call writelog('ls', '', 'Fortran outputformat')
         !call var_output_init(sglobal,slocal,par,tpar)
         call fortoutput_init(sglobal,par,tpar)
       case(OUTPUTFORMAT_NETCDF)
         ! only netcdf, stop if it's not build
         call writelog('ls', '', 'NetCDF outputformat')
#ifdef USENETCDF
         call ncoutput_init(sglobal,slocal,par,tpar)
#else
         call writelog('lse', '', 'This xbeach executable has no netcdf support. ', &
         'Rebuild with netcdf or run with outputformat=fortran')
         call halt_program
#endif
       case(OUTPUTFORMAT_DEBUG)
         call writelog('ls', '', 'Debug outputformat, writing both netcdf and fortran output')
#ifdef USENETCDF
         call writelog('ls', '', 'NetCDF outputformat')
         call ncoutput_init(sglobal,slocal,par,tpar)
#endif
         !call var_output_init(sglobal,slocal,par,tpar)
         call writelog('ls', '', 'Fortran outputformat')
         call fortoutput_init(sglobal,par,tpar)
      endselect

   end subroutine output_init
   !_________________________________________________________________________________

   subroutine output(sglobal,s,par,tpar,update)

      use means_module
      use postprocessmod

      implicit none

      type(spacepars)                     :: s,sglobal
      type(parameters)                    :: par
      type(timepars)                      :: tpar

      logical, optional                   :: update
      logical                             :: lupdate

      logical                               :: end_program
#ifdef USEMPI
      logical                               :: toall = .true.
#endif


      if (present(update)) then
         lupdate = update
      else
         lupdate = .true.
      endif

      ! update output times
      if(xcompute) then
         if (lupdate) call outputtimes_update(par, tpar)
      endif
      ! update log
      call log_progress(par)


      ! update meanvars in current averaging period with current timestep
      if(par%nmeanvar .gt. 0) then
         if(xcompute) then
            if (par%t>tpar%tpm(1) .and. par%t<=tpar%tpm(size(tpar%tpm))) then
               call makeaverage(s,par)
            endif
         endif
      endif

      end_program = .false.
#ifdef USEMPI
      if(xcompute) then
         if(tpar%output) then
            call xmpi_send_sleep(xmpi_imaster,xmpi_omaster) ! wake up omaster
            call xmpi_bcast(end_program,toall) ! matching the xmpi_bcast
            !                                  ! in the do loop a few lines below
            call tell_xomaster_what_time_it_is ! matching the call a few lines below
         else
            return
         endif
      endif
#endif

      do
         ! xomaster will not leave this loop, except
         ! at the very end of the program
#ifdef USEMPI
         if(xomaster) then
            call xmpi_send_sleep(xmpi_imaster,xmpi_omaster)
            call xmpi_bcast(end_program,toall) ! matching the xmpi_bcast
            !                                  ! above or the xmpi_bcast
            !                                  ! in finalize
            if(end_program) then
               call xmpi_barrier(toall)
               call xmpi_finalize
               stop
            endif
            call tell_xomaster_what_time_it_is ! matching the call a few lines above
         endif
#endif

         ! output

         call ncoutput(sglobal,s,par, tpar)

         ! clear averages after output of means
         if (tpar%outputm .and. tpar%itm>1) then
            call clearaverage(par)
         endif
         if (xcompute) exit
      enddo

#ifdef USEMPI
   contains
      subroutine tell_xomaster_what_time_it_is

         call xmpi_send(xmpi_imaster,xmpi_omaster,tpar%tnext)

         call xmpi_send(xmpi_imaster,xmpi_omaster,tpar%itg)
         call xmpi_send(xmpi_imaster,xmpi_omaster,tpar%itp)
         call xmpi_send(xmpi_imaster,xmpi_omaster,tpar%itm)
         call xmpi_send(xmpi_imaster,xmpi_omaster,tpar%itw)

         call xmpi_send(xmpi_imaster,xmpi_omaster,tpar%outputg)
         call xmpi_send(xmpi_imaster,xmpi_omaster,tpar%outputp)
         call xmpi_send(xmpi_imaster,xmpi_omaster,tpar%outputm)
         call xmpi_send(xmpi_imaster,xmpi_omaster,tpar%outputw)
         call xmpi_send(xmpi_imaster,xmpi_omaster,tpar%output)

         call xmpi_send(xmpi_imaster,xmpi_omaster,par%t)

      end subroutine tell_xomaster_what_time_it_is
#endif
   end subroutine output
   !_________________________________________________________________________________

   subroutine output_error(s, sglobal, par, tpar)

      use logging_module
      use xmpi_module, only: halt_program

      implicit none

      type(spacepars)                     :: s,sglobal
      type(parameters)                    :: par
      type(timepars)                      :: tpar

      !call output(s, sglobal, par, tpar, update=.false.)
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

