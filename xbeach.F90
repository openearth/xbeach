program xbeach

use params
use spaceparams
use xmpi_module
use initialize
use boundaryconditions
use drifter_module
use flow_timestep_module
use morphevolution
use outputmod
use readtide_module
use readwind_module
use wave_stationary_module
use wave_timestep_module
use timestep_module
use readkey_module
use groundwaterflow
use logging_module
! IFDEF used in case netcdf support is not compiled, f.i. Windows (non-Cygwin)
#ifdef USENETCDF
use ncoutput_module
#endif


IMPLICIT NONE

type(parameters)         :: par
type(timepars)           :: tpar
type(spacepars), pointer :: s
type(spacepars), target  :: sglobal
type(spacepars), target  :: slocal
character(len=80)        :: dummystring
character(len=8)         :: date
character(len=10)        :: time
character(len=5)         :: zone

logical                  :: newstatbc

integer                  :: it,error
real*8                   :: tbegin,tend
#ifdef USEMPI
real*8                   :: t0,t01,t1
#endif

character(len=155)       :: cwd ! for printing the working dir

! ----------------------------
! Initialize program
! ----------------------------

! autotools 
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

! subversion information
include 'version.def'
include 'version.dat'

! Setup of MPI 
#ifdef USEMPI
s=>slocal
call xmpi_initialize
t0 = MPI_Wtime()
#endif

! Start up log files
call start_logfiles(error)
if (error==1) then
   write(*,*) 'Error: not able to open log file. Please contact XBeach team. Stopping simulation'
   stop
endif
! 
call cpu_time(tbegin)
call DATE_AND_TIME(DATE=date, TIME=time, ZONE=zone)
! only run this on linux
call getcwd(cwd)

if (xmaster) then
  call writelog('ls','','**********************************************************')
  call writelog('ls','','                   Welcome to XBeach                      ')
  call writelog('ls','','                                                          ')
  call writelog('ls','','            revision ',trim(Build_Revision)                )
  call writelog('ls','','            date ',trim(Build_Date)                        )
  call writelog('ls','',' URL: ',trim(Build_URL)                                    )
  call writelog('ls','','**********************************************************')
  call writelog('ls','','                                                          ')
  call writelog('ls','','Simulation started: YYYYMMDD    hh:mm:ss     time zone (UTC)')
  call writelog('ls','','                    '//date //'  '//time(1:2)//':'//time(3:4)//':'//time(5:6)//'     '//zone)
  call writelog('ls','','                                                          ')
  call writelog('ls','',' running in: ',cwd )
  call writelog('ls','','General Input Module')
#ifdef USEMPI
  if(xmaster) then
    call writelog('ls','','MPI version, running on ',xmpi_size,'processes')
  endif
#endif
endif


! ----------------------------
! Initialize simulation
! ----------------------------

! TODO: move these to a params structure?
it=0
newstatbc=.true. ! This really shouldn't be here.

! General input per module
!
! This routine does need all processes, so not just xmaster ! Robert
call all_input(par)
! Do check of params.txt to spot errors 
! TODO: This shouldn't be in the main
if (xmaster) call readkey('params.txt','checkparams',dummystring) 
! TODO: We're not stepping into the timeloop just yet....
call writelog('ls','','Stepping into the time loop ....')   ! writelog is xmaster aware

#ifdef USEMPI
call distribute_par(par)
#endif

if (xmaster) then
  call writelog('l','' ,'------------------------------------')
  call writelog('ls','','Building Grid and Bathymetry and....')
  call writelog('ls','','Distributing wave energy across the directional space ....')
  call writelog('l','', '------------------------------------')
endif
! Grid and bathymetry

!
! grid_bathy will allocate x,y,xz,yz,xu,yv,xw,yw,zb,zb0 only
! on master process 
!
call space_alloc_scalars(sglobal)
s => sglobal
call grid_bathy(s,par)  ! s%nx and s%ny are available now
#ifdef USEMPI
call xmpi_determine_processor_grid(s%nx,s%ny,par%mpiboundary,error)
if(xmaster) then
  if (error==1) then
     call writelog('els','','Unknown mpi division ',par%mpiboundary)
	 call halt_program
  else
     call writelog('ls','','processor grid: ',xmpi_m,' X ',xmpi_n)
  endif
endif
#endif

! initialize timesteps
call timestep_init(par, tpar)

if (xmaster) then
  ! Jump into subroutine readtide
  call readtide (s,par)  !Ap 15/10 ! runs oonly on master wwvv
  call readwind (s,par)  !Robert 8/7/2009 only on master

  call writelog('ls','','Initializing .....')
! Initialisations
  call wave_init (s,par)  ! Always do this       wave_init only works on master process
endif
#ifdef USEMPI
! some of par has been changed, so:
call distribute_par(par)
#endif
if (xmaster) then
   call flow_init (s,par)  ! Always do this      works only on master process
   call gwinit(par,s)      ! works only on master process
   call sed_init (s,par)   ! works only on master process
endif

! initialize the correct output module (clean this up?, move to another module?)
if (par%outputformat=='fortran') then
   ! only fortran
   call writelog('ls', '', 'Fortran outputformat')
   call output_init(sglobal,slocal,par,tpar)
elseif (par%outputformat=='netcdf') then
   ! only netcdf, stop if it's not build
   call writelog('ls', '', 'NetCDF outputformat')
#ifdef USENETCDF
   call ncoutput_init(sglobal,slocal,par,tpar)
#else
   call writelog('lse', '', 'This xbeach executable has no netcdf support. Rebuild with netcdf or outputformat=fortran')
   stop 1
#endif
elseif (par%outputformat=='debug') then
   call writelog('ls', '', 'Debug outputformat, writing both netcdf and fortran output')
#ifdef USENETCDF
   call ncoutput_init(sglobal,slocal,par,tpar)
#endif
   call output_init(sglobal,slocal,par,tpar)
endif


#ifdef USEMPI
! some par has changed, so:
call distribute_par(par)
#endif

#ifdef USEMPI
s => slocal
!
!  determine how to divide the submatrices on the processor grid
!  distribute all values in sglobal to slocal
!  nx and ny will be adjusted in slocal
!  arrays is,js,lm,ln (describing the distribution) will
!  be filled in slocal
!  Note: slocal is available on all nodes, including master
!
call space_distribute_space(sglobal,slocal,par)
!call space_consistency(slocal,'ALL')
#endif

! update times at which we need output
call outputtimes_update(par, tpar)
! Store first timestep (always)
if (par%outputformat=='fortran') then
   call var_output(it,sglobal,s,par,tpar)
elseif (par%outputformat=='netcdf') then
#ifdef USENETCDF
   call nc_output(sglobal,s,par, tpar)
#endif
elseif (par%outputformat=='debug') then
#ifdef USENETCDF
   call nc_output(sglobal,s,par, tpar)
#endif
   call var_output(it,sglobal,s,par,tpar)
endif

! ----------------------------
! This is the main time loop
! ----------------------------
do while (par%t<par%tstop)
   ! Calculate timestep
   call timestep(s,par,tpar, it)
   ! Wave boundary conditions
   if (par%swave==1) call wave_bc (sglobal,slocal,par,newstatbc)
   ! Flow boundary conditions
   if (par%gwflow==1) call gwbc(par,s)
   if (par%flow+par%nonh>0) call flow_bc (s,par)
   !Dano moved here, after (long) wave bc generation
   if (it==0) then
#ifdef USEMPI
      t01 = MPI_Wtime()
#endif
   endif
   ! Wave timestep
   if (par%swave==1) then
      if (trim(par%instat) == 'stat' .or. trim(par%instat) == 'stat_table') then
         if ((abs(mod(par%t,par%wavint))<0.000001d0).or.newstatbc) then
            call wave_stationary(s,par)
            newstatbc=.false.
         endif
      else
	     newstatbc=.false.
         call wave_timestep(s,par)
      endif
   endif
   ! Flow timestep
   if (par%gwflow==1) call gwflow(par,s)
   if (par%flow+par%nonh>0) call flow_timestep (s,par)
   if (par%ndrifter>0) call drifter (s,par)
   ! Suspended transport
   if(par%sedtrans==1) call transus(s,par)
   ! Bed level update
   if (par%morphology==1) call bed_update(s,par)
   ! Calculate new output times, so we know when to stop
   call outputtimes_update(par, tpar)
   ! Output
   if (par%outputformat=='fortran') then
      call var_output(it,sglobal,s,par,tpar)
   elseif (par%outputformat=='netcdf') then
#ifdef USENETCDF
      call nc_output(sglobal,s,par, tpar)
#endif
   elseif (par%outputformat=='debug') then
#ifdef USENETCDF
      call nc_output(sglobal,s,par, tpar)
#endif
      call var_output(it,sglobal,s,par,tpar)
   endif
enddo
! ------------------------
! End of main time loop
! ------------------------
!
! Finish files
if(xmaster) then
   call cpu_time(tend)
   call writelog('ls','','Total calculation time: ',tend-tbegin,' seconds')
#ifdef USEMPI
   t1 = MPI_Wtime()
   call writelog('ls','','Timing: procs: ',xmpi_size,' seconds: total:',t1-t0,'loop: ',t1-t01)
#endif
   call writelog('ls','','End of program xbeach')
endif
call close_logfiles
#ifdef USEMPI
call xmpi_finalize
#endif

end program


