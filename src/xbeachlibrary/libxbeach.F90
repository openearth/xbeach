module libxbeach_module
  use iso_c_binding
  use mnemiso_module
  use mnemmodule
  use getkey_module
  use params
  use spaceparams
  use xmpi_module
  use initialize
  use boundaryconditions
  use drifter_module
  use flow_timestep_module
  use morphevolution
  use readtide_module
  use readwind_module
  use wave_timestep_module
  use timestep_module
  use readkey_module
  use beachwizard_module
  use groundwaterflow
  use logging_module
  use means_module
  use output_module
  use ship_module
  implicit none

  type(parameters), save                              :: par
  type(timepars), save                                     :: tpar
  type(spacepars), pointer                            :: s
  type(spacepars), target, save                       :: sglobal
  type(ship), dimension(:), pointer                   :: sh

  integer                                             :: n,it,error
  real*8                                              :: tbegin

#ifdef USEMPI
  type(spacepars), target, save                       :: slocal
  real*8                                              :: t0,t01
#endif


  !startinit

  !-----------------------------------------------------------------------------!
  ! Initialize program                                                          !
  !-----------------------------------------------------------------------------!


contains
  integer(c_int) function init()

    error   = 0

    ! setup of MPI 
#ifdef USEMPI
    s=>slocal
    call xmpi_initialize
    t0 = MPI_Wtime()
#endif

    ! create log files
    call start_logfiles(error)

    ! set starting time and date
    call cpu_time(tbegin)

    ! show statup message
    call writelog_startup()

    !-----------------------------------------------------------------------------!
    ! Initialize simulation                                                       !
    !-----------------------------------------------------------------------------!

    ! initialize time counter
    it      = 0

    ! read input from params.txt
    call all_input(par)

    ! allocate space scalars
    call space_alloc_scalars(sglobal)
    s => sglobal

    ! read grid and bathymetry
    call grid_bathy(s,par)

    ! distribute grid over processors
#ifdef USEMPI
    call xmpi_determine_processor_grid(s%nx,s%ny,par%mpiboundary,error)
    call writelog_mpi(par%mpiboundary,error)
#endif

    ! initialize timestep
    call timestep_init(par, tpar)

    if (xmaster) then

       call writelog('ls','','Initializing .....')

       ! initialize physics
       call readtide           (s,par)
       call readwind           (s,par)

       call flow_init          (s,par)
       call discharge_init     (s,par)
       call drifter_init       (s,par)
       call wave_init          (s,par)
       call gw_init            (s,par)
       ! TODO, fix ordening of arguments....
       call bwinit             (s,par)		! works only on master process 

       call sed_init           (s,par)
       if (par%ships==1)        call ship_init (s,par,sh)

    endif

#ifdef USEMPI
    call distribute_par(par)
    s => slocal
    call space_distribute_space (sglobal,s,par     )
#endif

    ! initialize output
    call means_init             (sglobal,s,par     )
    call output_init            (sglobal,s,par,tpar)


    ! store first timestep
    call output(sglobal,s,par,tpar)
    init = 0
  end function init

  integer(c_int) function outputext()
    ! store first timestep
    call output(sglobal,s,par,tpar)
    outputext = 0
  end function outputext
  !-----------------------------------------------------------------------------!
  ! Start simulation                                                            !
  !-----------------------------------------------------------------------------!

  integer(c_int) function executestep()


    ! This is now called in output, but timestep depends on it...
    ! call outputtimes_update(par, tpar)
    !n = 0
    executestep = -1
    !do while (par%t<par%tstop)
    ! determine timestep
    call timestep(s,par,tpar,it,ierr=error)
    if (error==1) call output_error(s,sglobal,par,tpar)

    ! boundary conditions
    call wave_bc        (sglobal,s,par)
    if (par%gwflow==1)       call gw_bc          (s,par)
    if (par%flow+par%nonh>0) call flow_bc        (s,par)

#ifdef USEMPI
    if (it==0) t01 = MPI_Wtime()
#endif

    ! compute timestep
    if (par%ships==1)        call shipwave       (s,par,sh)
    if (par%swave==1)        call wave           (s,par)
    if (par%gwflow==1)       call gwflow         (s,par)
    if (par%flow+par%nonh>0) call flow           (s,par)
    if (par%ndrifter>0)      call drifter        (s,par)
    if (par%sedtrans==1)     call transus        (s,par)
    ! Beach wizard
    if (par%bchwiz>0)        call assim          (s,par)
    ! Bed level update
    if ((par%morphology==1) .and. (.not. par%bchwiz == 1)) call bed_update     (s,par)
    if (par%bchwiz>0)        call assim_update   (s, par)
    ! n = n + 1
    executestep = 0
    ! enddo
  end function executestep



  integer(c_int) function finalize()


    !-----------------------------------------------------------------------------!
    ! Finalize simulation                                                         !
    !-----------------------------------------------------------------------------!

#ifdef USEMPI
    call writelog_finalize(tbegin,n,par%t,par%nx,par%ny,t0,t01)
    call xmpi_finalize
#else
    call writelog_finalize(tbegin,n,par%t,par%nx,par%ny)
#endif
    finalize = 0
  end function finalize
end module libxbeach_module
