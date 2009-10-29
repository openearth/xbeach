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


IMPLICIT NONE

type(parameters)         :: par
type(spacepars), pointer :: s
type(spacepars), target  :: sglobal
type(spacepars), target  :: slocal
character(len=80)        :: dummystring
character(len=10)        :: date,time,zone
logical                  :: newstatbc

integer                  :: it
real*8                   :: tbegin,tend
#ifdef USEMPI
real*8                   :: t0,t01,t1
#endif

! subversion information
include 'version.def'
include 'version.dat'

!
!build_revision = '$Revision$'
!build_date = '$Date$'
!build_url = '$HeadURL$'

! autotools 
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif



#ifdef USEMPI
s=>slocal
call xmpi_initialize
t0 = MPI_Wtime()
#endif

call cpu_time(tbegin)
call DATE_AND_TIME(DATE=date, TIME=time, ZONE=zone)

if (xmaster) then
  write(*,*)'**********************************************************'
  write(*,*)'                   Welcome to XBeach                      '
  write(*,*)'                                                          '
  write(*,*)'            revision ',trim(Build_Revision)
  write(*,*)'            date $Date$'
  write(*,*)' URL: $HeadURL$ '
  write(*,*)'**********************************************************'
  write(*,*)'                                                          '
  write(*,*)'Simulation started: YYYYMMDD    hh:mm:ss     time zone (UTC)'
  write(*,*)'                    ',date(1:10),'  ',time(1:2),':',time(3:4),':',time(5:6),'     ',zone(1:5)
  write(*,*)'                                                          '
  write(*,*)'General Input Module'
#ifdef USEMPI
  if(xmaster) then
    write(*,*) 'MPI version, running on ',xmpi_size,'processes'
  endif
#endif
endif

! General input per module
!
! the basic input routines, used by the following three subroutines
! are MPI-aware, no need to do something special here
!

par%t=0.d0
it=0
newstatbc=.true.

call wave_input(par)
call flow_input(par)
call sed_input(par)

#ifdef USEMPI
call distribute_par(par)
#endif

if (xmaster) then
  write(*,*) 'Building Grid and Bathymetry and....'
  write(*,*) 'Distributing wave energy across the directional space ....'     
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
call xmpi_determine_processor_grid(s%nx,s%ny)
if(xmaster) then
  write(*,*) 'processor grid: ',xmpi_m,' X ',xmpi_n
endif
#endif

if (xmaster) then
  ! Jump into subroutine readtide
  call readtide (s,par)  !Ap 15/10 ! runs oonly on master wwvv
  call readwind (s,par)  !Robert 8/7/2009 only on master

  write(*,*) 'Initializing .....'
! Initialisations
  call wave_init (s,par)  ! wave_init only works on master process
endif
#ifdef USEMPI
! some of par has been changed, so:
call distribute_par(par)
#endif
if (xmaster) then
   call flow_init (s,par)  ! works only on master process
   call gwinit(par,s)      ! works only on master process
   call sed_init (s,par)   ! works only on master process
endif
call init_output(sglobal,slocal,par,it)
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
call printit(sglobal,slocal,par,it,'after space_distribute_space')

! If wanted, produce output at t=0
if (it==1) call var_output(it,sglobal,s,par)

if (xmaster) then
  ! Do check of params.txt to spot errors
  call readkey('params.txt','checkparams',dummystring) 
  write(*,*) 'Stepping into the time loop ....'  
endif

!#ifdef USEMPI
!t01 = MPI_Wtime()
!#endif
do while (par%t<par%tstop)
    ! Calculate timestep
    call timestep(s,par,it)
    ! Wave boundary conditions
    call wave_bc (sglobal,slocal,par,newstatbc)
    call printit(sglobal,slocal,par,it,'after wave_bc')
    ! Flow boundary conditions
	if (par%gwflow==1) call gwbc(par,s)
	call flow_bc (s,par)
    !Dano moved here, after (long) wave bc generation
	if (it==1) then
#ifdef USEMPI
       t01 = MPI_Wtime()
#endif
    endif
    call printit(sglobal,slocal,par,it,'after flow_bc')
#ifdef USEMPI
    !call space_consistency(slocal,'ALL')
#endif
    ! Wave timestep
    if (par%instat==0.or.par%instat==40) then
       if ((abs(mod(par%t,par%wavint))<0.000001d0).or.newstatbc) then
          call wave_stationary(s,par)
		  newstatbc=.false.
          call printit(sglobal,slocal,par,it,'after wave_stationary')
       endif
    else
      call wave_timestep(s,par)
          call printit(sglobal,slocal,par,it,'after wave_timestep')
    endif
    ! Flow timestep
	if (par%gwflow==1) call gwflow(par,s)
    call flow_timestep (s,par)
    call drifter (s,par)
    call printit(sglobal,slocal,par,it,'after flow_timestep')
    ! Suspended transport
    call transus(s,par)
    call printit(sglobal,slocal,par,it,'after transus')
    ! Bed level update
    call bed_update(s,par)
    call printit(sglobal,slocal,par,it,'after bed_update')
    ! Output
    call var_output(it,sglobal,s,par)
#ifdef USEMPI
!   varoutput changed some in parameters, so:
    call distribute_par(par)
#endif

enddo


if(xmaster) then
call cpu_time(tend)
write(*,*)'Total calculation time: ',tend-tbegin,' seconds'
#ifdef USEMPI
    t1 = MPI_Wtime()
    write(*,*)'Timing: procs: ',xmpi_size,' seconds: total:',t1-t0,&
            'loop: ',t1-t01
#endif
  write(*,*)'End of program xbeach'
endif
#ifdef USEMPI
call xmpi_finalize
#endif

end program

subroutine printit(sglobal,slocal,par,it,s)
  use spaceparams
  use params
  use xmpi_module
  IMPLICIT none
  type(spacepars)          :: sglobal,slocal
  type(parameters)         :: par
  integer                  :: it
  character(len=*)         :: s
  integer,save             :: iter=0 
  return
  iter = iter+1
#ifdef USEMPI
  write(*,*) par%t,xmpi_rank,trim(s)
  call space_collect(slocal,sglobal%H,slocal%H)
  call space_collect(slocal,sglobal%zs,slocal%zs)
  call space_collect(slocal,sglobal%zs0,slocal%zs0)
  call space_collect(slocal,sglobal%u,slocal%u)
  call space_collect(slocal,sglobal%uu,slocal%uu)
  call space_collect(slocal,sglobal%ui,slocal%ui)
  call space_collect(slocal,sglobal%hh,slocal%hh)
  call space_collect(slocal,sglobal%vu,slocal%vu)
  call space_collect(slocal,sglobal%v,slocal%v)
  !call space_consistency(slocal,'ALL')
#else
  slocal%nx = slocal%nx  ! to prevent compiler warning about
                         ! unused slocal
#endif
  if(xmaster) call printsum(6,'H',1000*it+iter,sglobal%H)
  if(xmaster) call printsum(6,'zs',1000*it+iter,sglobal%zs)
  if(xmaster) call printsum(6,'zs0',1000*it+iter,sglobal%zs0)
  if(xmaster) call printsum(6,'u',1000*it+iter,sglobal%u)
  if(xmaster) call printsum(6,'uu',1000*it+iter,sglobal%uu)
  if(xmaster) call printsum(6,'ui',1000*it+iter,sglobal%ui)
  if(xmaster) call printsum(6,'hh',1000*it+iter,sglobal%hh)
  if(xmaster) call printsum(6,'vu',1000*it+iter,sglobal%vu)
  if(xmaster) call printsum(6,'v',1000*it+iter,sglobal%v)

  if(xmaster) print *,'par%t:',par%t
  if(xmaster) print *,'par%zs01:',par%zs01
#ifdef USEMPI
  if(xmaster) print *,'s%tideinpt:',slocal%tideinpt
  if(xmaster) print *,'s%tideinpz:',slocal%tideinpz(:,1)
#else
  if(xmaster) print *,'s%tideinpt:',sglobal%tideinpt
  if(xmaster) print *,'s%tideinpz:',sglobal%tideinpz(:,1)
#endif
end subroutine printit
