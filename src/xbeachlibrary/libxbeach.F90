module libxbeach_module
   use loopcounters_module
   use iso_c_binding
   use mnemiso_module
   use mnemmodule
   use getkey_module
   use params
   use spaceparams
   use xmpi_module
   use initialize_module
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
   use nonh_module
   use vegetation_module

   implicit none
   save

   type(parameters)                     :: par
   type(timepars)                       :: tpar
   type(spacepars), pointer             :: s
   type(spacepars), target              :: sglobal
   type(ship), dimension(:), pointer    :: sh

   integer                              :: n,it,error
   real*8                               :: tbegin

#ifdef USEMPI
   type(spacepars), target              :: slocal
   real*8                               :: t0,t01
   logical                              :: toall = .true.
   logical                              :: end_program
   integer                              :: nxbak, nybak
#endif


   !startinit

   !-----------------------------------------------------------------------------!
   ! Initialize program                                                          !
   !-----------------------------------------------------------------------------!


contains
   integer(c_int) function init()
#ifdef USEMPI
      integer, dimension(12)               :: info
      character(256)                       :: line
      integer                              :: rank,i
#endif

      error   = 0
      n = 0

      ! setup of MPI
#ifdef USEMPI
      s=>slocal
      call xmpi_initialize
      call xmpi_barrier(toall)
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
      params_inio = .false.
      call all_input(par)

      ! allocate space scalars
      call space_alloc_scalars(sglobal)
      s => sglobal

      ! read grid and bathymetry
      call grid_bathy(s,par)

      ! distribute grid over processors
#ifdef USEMPI
      call xmpi_determine_processor_grid(s%nx,s%ny,par%mpiboundary,par%mmpi,par%nmpi,par%cyclic,error)
#if 0
      ! print information about the neighbours of the processes
      info                      = 0
      info(1)                   = xmpi_orank
      info(2)                   = xmpi_rank
      info(3)                   = xmpi_prow
      info(4)                   = xmpi_pcol
      info(5)                   = xmpi_left
      info(6)                   = xmpi_right
      info(7)                   = xmpi_top
      info(8)                   = xmpi_bot
      if(xmpi_isleft)  info(9)  = 1
      if(xmpi_isright) info(10) = 1
      if(xmpi_istop)   info(11) = 1
      if(xmpi_isbot)   info(12) = 1

      do i=5,8
         if(info(i) .eq. MPI_PROC_NULL) then
            info(i) = -99
         endif
      enddo

      call writelog("ls"," "," ranks and neigbours (-99 means: no neighbour):")
      call writelog("ls"," ",' ')
      call writelog("ls"," ","  orank rank pcol prow left right  top  bot isleft isright istop isbot")

      do rank = 0,xmpi_osize-1
         if (rank .ne. xmpi_omaster) then
            call xmpi_send(rank,xmpi_imaster,info)
            if (xmaster) then
               write(line,'(i5,i5,i5,i5,i5,i6,i5,i5,i7,i8,i6,i6)') info
            endif
            call writelog("ls"," ",trim(line))
         endif
      enddo
      call writelog("ls"," ",' ')
#endif
      call writelog_mpi(par%mpiboundary,error)
#endif

      ! initialize timestep
      call timestep_init(par, tpar)

      if (xmaster) then

         call writelog('ls','','Initializing .....')
      endif
      call setbathy_init      (s,par)
      ! initialize physics
      call readtide           (s,par)
      call readwind           (s,par)
      call flow_init          (s,par)
      call discharge_init     (s,par)
      call drifter_init       (s,par)
      call wave_init          (s,par)
      call gw_init            (s,par)
      ! TODO, fix ordening of arguments....
      call bwinit             (s)          ! works only on master process

      call sed_init           (s,par)

      call ship_init          (s,par,sh)   ! always need to call initialise in order
      ! to reserve memory on MPI subprocesses.
      ! Note: if par%ships==0 then don't allocate
      ! and read stuff for sh structures
      call vegie_init         (s,par)

#ifdef USEMPI
      call distribute_par(par)
      s => slocal
      !
      ! here an hack to ensure that sglobal is populated, also on
      ! the not-(o)master processes, just to get valid addresses.
      !
      if (.not. xmaster .and. .not. xomaster) then
         nxbak = sglobal%nx
         nybak = sglobal%ny
         sglobal%nx=0
         sglobal%ny=0
         call space_alloc_arrays(sglobal,par)
         sglobal%nx = nxbak
         sglobal%ny = nybak
      endif
      call space_distribute_space (sglobal,s,par     )
#endif

      call ranges_init(s)

      ! nonh_init does not always need to be called
      if (par%nonh==1) call nonh_init(s,par)

      ! initialize output
      call means_init             (sglobal,s,par     )

      call output_init            (sglobal,s,par,tpar)


      ! store first timestep
      ! from this point on, xomaster will hang in subroutine output
      ! until a broadcast .true. is received
      call output(sglobal,s,par,tpar)
      init = 0
   end function init

   integer(c_int) function outputext()
      ! store first timestep
      call output(sglobal,s,par,tpar,.false.)
      outputext = 0
   end function outputext
   !-----------------------------------------------------------------------------!
   ! Start simulation                                                            !
   !-----------------------------------------------------------------------------!

   !_____________________________________________________________________________

  integer(c_int) function executestep(dt)

    real*8, optional :: dt

#ifdef USEMPI
      if (execute_counter .eq. 1) then
         ! exclude first pass from time measurement
         call xmpi_barrier
         t01 = MPI_Wtime()
      endif
#endif
      execute_counter = execute_counter + 1

      call outputtimes_update(par, tpar)
      executestep = -1

      ! determine timestep
      if(xcompute) then

         ! determine time step
         call timestep(s,par,tpar,it,dt=dt,ierr=error)
         
         if (error==1) then
            call output_error(s,sglobal,par,tpar)
         endif

         ! boundary conditions
         call wave_bc        (sglobal,s,par)
         if (par%gwflow==1)       call gw_bc          (s,par)
         if (par%flow+par%nonh>0) call flow_bc        (s,par)

         ! compute timestep
         if (par%ships==1)        call shipwave       (s,par,sh)
         if (par%swave==1)        call wave           (s,par)
         if (par%vegetation==1)   call vegatt         (s,par)
         if (par%gwflow==1)       call gwflow         (s,par)
         if (par%flow+par%nonh>0) call flow           (s,par)
         if (par%ndrifter>0)      call drifter        (s,par)
         if (par%sedtrans==1)     call transus        (s,par)
         ! Beach wizard
         if (par%bchwiz>0)        call assim          (s,par)
         ! Bed level update
         if ((par%morphology==1).and.(.not. par%bchwiz == 1).and.(.not. par%setbathy==1)) call bed_update(s,par)
         if (par%bchwiz>0)        call assim_update   (s, par)
         if (par%setbathy==1)     call setbathy_update(s, par)

      endif

      n = n + 1
      executestep = 0
   end function executestep
   !_____________________________________________________________________________


   integer(c_int) function final()


      !-----------------------------------------------------------------------------!
      ! Finalize simulation                                                         !
      !-----------------------------------------------------------------------------!

#ifdef USEMPI
      end_program = .true.
      call xmpi_send_sleep(xmpi_imaster,xmpi_omaster) ! wake up omaster
      call xmpi_bcast(end_program,toall)
      call xmpi_barrier(toall)
      call writelog_finalize(tbegin,n,par%t,par%nx,par%ny,t0,t01)
      call xmpi_finalize
#else
      call writelog_finalize(tbegin,n,par%t,par%nx,par%ny)
#endif
      final = 0
   end function final

   subroutine getversion(version)
      character(kind=c_char,len=*),intent(inout) :: version

      include 'version.def'
      include 'version.dat'

      version = trim(Build_Revision)
   end subroutine

end module libxbeach_module
