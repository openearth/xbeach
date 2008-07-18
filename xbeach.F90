program xbeach

use params
use spaceparams
use xmpi_module
use initialize
use boundaryconditions
use flow_timestep_module
use morphevolution
use outputmod
use readtide_module
use wave_stationary_module
use wave_timestep_module
use timestep_module


IMPLICIT NONE

type(parameters)         :: par
type(spacepars), pointer :: s
type(spacepars), target  :: sglobal
type(spacepars), target  :: slocal

integer                  :: it
real*8                   :: tbegin,tend
#ifdef USEMPI
real*8                   :: t0,t01,t1
#endif

#ifdef USEMPI
call xmpi_initialize
t0 = MPI_Wtime()
if(xreader) then
  write(*,*) 'MPI version, running on ',xmpi_size,'processes'
endif
call xmpi_determine_processor_grid
if(xreader) then
  write(*,*) 'processor grid: ',xmpi_m,' X ',xmpi_n
endif
#endif

call cpu_time(tbegin)

if (xprinter) then
  write(*,*) 'Welcome to Xbeach'
  write(*,*) 'General Input Module'
endif

! General input per module
!
! the basic input routines, used by the following three subroutines
! are MPI-aware, no need to do something special here
!
call wave_input(par)
call flow_input(par)
call sed_input(par)

#ifdef USEMPI
call distribute_par(par)
#endif

if (xprinter) then
  write(*,*) 'Building Grid and Bathymetry and....'
  write(*,*) 'Distributing wave energy across the directional space ....'     
endif
! Grid and bathymetry

!
! grid_bathy will allocate x,y,xz,yz,xu,yv,xw,yw,zb,zb0 only
! on master process 
!
s => sglobal
call grid_bathy(s,par)

! Jump into subroutine readtide
if(par%tideloc>=1)then 
   call readtide (s,par)  !Ap 15/10
end if

par%t=0.d0

call init_output(s,par)

if (xprinter) then
  write(*,*) 'Initializing .....'
endif
! Initialisations
call wave_init (s,par)  ! wave_init only works on master process
#ifdef USEMPI
! some of par has been changed, so:
call distribute_par(par)
#endif
call flow_init (s,par)  ! works only on master process
call sed_init (s,par)       ! works on all processes

#ifdef USEMPI

s => slocal

call space_distribute_space(sglobal,slocal,par)

#endif

par%t=0.0d0
par%tnext=par%tint
it=0

if (xprinter) then
  write(*,*) 'Stepping into the time loop ....'   
endif

#ifdef USEMPI
t01 = MPI_Wtime()
#endif
do while (par%t<par%tstop)
    ! Calculate timestep
        call timestep(s,par,it)
        ! Wave boundary conditions
    if(xreader) then
      !print *,'calling wave_bc'
    endif
    call wave_bc (sglobal,slocal,par)
#ifdef USEMPI
       !DANO Communicate ee,rr
       call space_shift_borders(s,s%ee)  
       call space_shift_borders(s,s%rr)  
       !WILLEM Communicate ui
       !call space_shift_borders(s,s%ui)
#endif
    ! Flow boundary conditions
    call flow_bc (s,par)
#ifdef USEMPI
    !DANO Communicate uu,vv,zs; hh uitzoeken
    call space_shift_borders(s,s%uu)
    call space_shift_borders(s,s%vv)
    call space_shift_borders(s,s%zs)
    call space_shift_borders(s,s%hh)
#endif
    ! Wave timestep
    if (par%instat==0) then
       if (mod(par%t,real(par%wavint))==0) then
          if(xreader) then
            !print *,'calling wave_stationary'
          endif
          call wave_stationary(s,par)
       endif
    else
      if(xreader) then
        !print *,'calling wave_timestep'
      endif
      call wave_timestep(s,par)
    endif
    ! Flow timestep
    if(xreader) then
      !print *,'calling flow_timestep'
    endif
    call flow_timestep (s,par)
    ! Suspended transport
    call transus(s,par)
#ifdef USEMPI
    !DANO communicate cc
    !wwvv cc does not exists anymore, communicate ccg instead
    call space_shift_borders(s,s%ccg)
#endif
    ! Bed level update
    call bed_update(s,par)
#ifdef USEMPI
    !DANO communicate zb
    call space_shift_borders(s,s%zb)
#endif
    ! Output
    call var_output(it,s,par)

enddo


if(xreader) then
call cpu_time(tend)
write(*,*)'Total calculation time: ',tend-tbegin,' seconds'
#ifdef USEMPI
    t1 = MPI_Wtime()
    print *,'Timing: procs: ',xmpi_size,' seconds: total:',t1-t0,&
            'loop: ',t1-t01
#endif
  print *,'End of program xbeach'
endif
#ifdef USEMPI
call xmpi_finalize
#endif

end program
