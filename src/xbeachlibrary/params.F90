module params
   use typesandkinds
   use mnemmodule
   use xmpi_module
   use paramsconst
   implicit none
   save

   ! before using any routine from this module, take care of the value of params_inio
   ! .true. : values will be broadcasted to all processes, including output process
   ! otherwise: values will be broadcasted to compute processes only
   !
   logical :: params_inio = .false.
   type parameters
      ! These parameters are constants, variables read from params.txt, or are scalars derived directly from read input
      !
      !  Please maintain current setup of type declaration in order to allow parsing of this file by script readers,
      !  such as the function "xb_get_params" in OpenEarth (http://public.deltares.nl/display/OET/OpenEarth).
      !
      !  Rules are:
      !  - [Section] markers indicate new sets of parameters, followed by the name of the set
      !  - Fortran declaration always as "kind  ::  name   = initial value"
      !  - After the declaration of a variable add exclamation mark, followed by [unit] and decription.
      !    If the parameter is essentially only for advanced users, follow the unit declation by "(advanced)"
      !    If the parameter is deprecated, but still used for backwards-compatibility, follow the unit declation by "(deprecated)"
      !  - Description of a variable may continue on a new line as long as the first character is "!" and the
      !    position of the first "!" is greater than 50 characters from the start of the line. Best practice
      !    is to keep in line with the start of the description on the line above.
      !  - Please keep the declaration of "globalvars", "meanvars" and "pointvars" on the line directly after their respective related size
      !    parameter declaration (i.e. declaration of "nglobalvars", "nmeanvar" and "npointvar"). This is needed for the autogeneration of
      !    parameters.inc and subsequent params.dat file.
      !  - To enable parsing of "all_input" subroutine please only use "if () then, ... , endif/else/elseif" blocks,
      !    rather than one line "if () ..." commands in the "all_input" subroutine.
      !
      !  Type             name                   initialize    !  [unit] (advanced/deprecated) description
      ! [Section] Physical processes
      integer*4     :: cyclic                     = -123    !  [-] Turn on cyclic boundary conditions
      integer*4     :: swave                      = -123    !  [-] Turn on short waves
      integer*4     :: lwave                      = -123    !  [-] Turn on short wave forcing on NLSW equations and boundary conditions
      integer*4     :: flow                       = -123    !  [-] Turn on flow calculation
      integer*4     :: sedtrans                   = -123    !  [-] Turn on sediment transport
      integer*4     :: morphology                 = -123    !  [-] Turn on morphology
      integer*4     :: avalanching                = -123    !  [-] Turn on avalanching
      integer*4     :: nonh                       = -123    !  [-] (advanced) Turn on non-hydrostatic pressure: 0 = NSWE, 1 = NSW + non-hydrostatic pressure compensation Stelling & Zijlema, 2003
      integer*4     :: gwflow                     = -123    !  [-] (advanced) Turn on groundwater flow
      integer*4     :: q3d                        = -123    !  [-] (advanced,silent) Turn on quasi-3D sediment transport
      integer*4     :: swrunup                    = -123    !  [-] (advanced,silent) Turn on short wave runup
      integer*4     :: ships                      = -123    !  [-] (advanced) Turn on ship waves
      integer*4     :: vegetation                 = -123    !  [-] (advanced) Turn on interaction of waves and flow with vegetation
      integer*4     :: snells                     = -123    !  [-] (advanced) Turn on Snell's law for wave refraction
      integer*4     :: single_dir                 = -123    !  [-] (advanced) Turn on stationary model for refraction, surfbeat based on mean direction
      integer*4     :: bchwiz                     = -123    !  [-] (advanced,silent) Turn on beachwizard
      integer*4     :: setbathy                   = -123    !  [-] Turn on timeseries of prescribed bathy input

      ! [Section] Grid parameters
      character(slen):: depfile                   = 'abc'   !  [file] Name of the input bathymetry file
      real*8        :: posdwn                     = -123    !  [-] Bathymetry is specified positive down (1) or positive up (-1)
      integer*4     :: nx                         = -123    !  [-] Number of computational cell corners in x-direction
      integer*4     :: ny                         = -123    !  [-] Number of computational cell corners in y-direction
      real*8        :: alfa                       = -123    !  [deg] Angle of x-axis from East
      integer*4     :: vardx                      = -123    !  [-] Switch for variable grid spacing
      real*8        :: dx                         = -123    !  [m] Regular grid spacing in x-direction
      real*8        :: dy                         = -123    !  [m] Regular grid spacing in y-direction
      character(slen) :: xfile                    = 'abc'   !  [file] Name of the file containing x-coordinates of the calculation grid
      character(slen) :: yfile                    = 'abc'   !  [file] Name of the file containing y-coordinates of the calculation grid
      real*8        :: xori                       = -123    !  [m] X-coordinate of origin of axis
      real*8        :: yori                       = -123    !  [m] Y-coordinate of origin of axis
      real*8        :: thetamin                   = -123    !  [deg] Lower directional limit (angle w.r.t computational x-axis)
      real*8        :: thetamax                   = -123    !  [deg] Higher directional limit (angle w.r.t computational x-axis)
      real*8        :: dtheta                     = -123    !  [deg] Directional resolution
      real*8        :: dtheta_s                   = -123    !  [deg] Directional resolution in case of stationary refraction
      integer*4     :: thetanaut                  = -123    !  [-] Switch to specify thetamin and thetamax in nautical convention rather than cartesian
      integer       :: gridform                   = -123    !  [name] Grid definition format
      character(slen) :: xyfile                   = 'abc'   !  [file] Name of the file containing Delft3D xy-coordinates of the calculation grid
      character(slen):: gridform_str              = ' '     !  [name] Grid definition format

      ! [Section] Model time
      real*8        :: tstop                      = -123    !  [s] Stop time of simulation, in morphological time
      real*8        :: CFL                        = -123    !  [-] Maximum Courant-Friedrichs-Lewy number
      character(slen)  :: tunits                  = 's'     !  [-] (advanced) Time units in udunits format (seconds since 1970-01-01 00:00:00.00 +1:00)

      ! [Section] Physical constants
      real*8        :: g                          = -123    !  [ms^-2] Gravitational acceleration
      real*8        :: rho                        = -123    !  [kgm^-3] Density of water
      real*8        :: depthscale                 = -123    !  [-] (advanced)  depthscale of (lab)test simulated, affects eps, hmin, hswitch and dzmax

      ! [Section] Initial conditions
      real*8        :: zs0                        = -123    !  [m] Inital water level
      character(slen):: zsinitfile                = 'abc'   !  [file] Name of inital water level file
      integer*4     :: hotstartflow               = -123    !  [-] (advanced) Switch for hotstart flow conditions with pressure gradient balanced by wind and bed stress

      ! [Section] Wave boundary condition parameters
      integer       :: instat                     = -123    !  [name] Wave boundary condition type
      character(slen):: instat_str                = ' '     !  [-] Wave boundary condition type
      real*8        :: taper                      = -123    !  [s] Spin-up time of wave boundary conditions, in morphological time
      real*8        :: Hrms                       = -123    !  [m] Hrms wave height for instat = stat, bichrom, ts_1 or ts_2
      real*8        :: Tm01                       = -123    !  [s] (deprecated) Old name for Trep
      real*8        :: Trep                       = -123    !  [s] Representative wave period for instat = stat, bichrom, ts_1 or ts_2
      real*8        :: Tlong                      = -123    !  [s] Wave group period for case instat = bichrom
      real*8        :: dir0                       = -123    !  [deg] Mean wave direction for instat = stat, bichrom, ts_1 or ts_2 (nautical convention)
      real*8        :: nmax                       = -123    !  [-] (advanced) maximum ratio of cg/c for computing long wave boundary conditions
      integer*4     :: m                          = -123    !  [-] Power in cos^m directional distribution for instat = stat, bichrom, ts_1 or ts_2
      integer       :: lateralwave                = -123    !  [name] Switch for lateral boundary at left
      character(slen):: lateralwave_str           = ' '
      integer       :: leftwave                   = -123    !  [-] (deprecated) old name for lateralwave
      character(slen) :: leftwave_str             =  ' '    !  [-] old name for lateralwave
      integer       :: rightwave                  = -123    !  [-] (deprecated) old name for lateralwave
      character(slen) :: rightwave_str            = ' '     !  [-] old name for lateralwave

      ! [Section] Wave-spectrum boundary condition parameters
      character(slen):: bcfile                    = 'abc'   !  [file] Name of spectrum file
      integer*4     :: random                     = -123    !  [-] (advanced) Switch to enable random seed for instat = jons, swan or vardens boundary conditions
      real*8        :: fcutoff                    = -123    !  [Hz] (advanced) Low-freq cutoff frequency for instat = jons, swan or vardens boundary conditions
      integer*4     :: nspr                       = -123    !  [-] (advanced) Switch to enable long wave direction forced into centres of short wave bins
      real*8        :: trepfac                    = -123    !  [-] (advanced) Compute mean wave period over energy band: par%trepfac*maxval(Sf) for instat jons, swan or vardens; converges to Tm01 for trepfac = 0.0 and
      real*8        :: sprdthr                    = -123    !  [-] (advanced) Threshold ratio to maximum value of S above which spectrum densities are read in
      integer*4     :: oldwbc                     = -123    !  [-] (deprecated,silent) (1) Use old version wave boundary conditions for instat jons, swan or vardens
      integer*4     :: correctHm0                 = -123    !  [-] (advanced,silent) Switch to enable Hm0 correction
      integer*4     :: oldnyq                     = -123    !  [-] (advanced,silent) Switch to enable old nyquist switch
      integer*4     :: Tm01switch                 = -123    !  [-] (advanced) Switch to enable Tm01 rather than Tm-10
      real*8        :: rt                         = -123    !  [s] Duration of wave spectrum at offshore boundary, in morphological time
      real*8        :: dtbc                       = -123    !  [s] (advanced) Timestep used to describe time series of wave energy and long wave flux at offshore boundary (not affected by morfac)
      real*8        :: dthetaS_XB                 = -123    !  [deg] (advanced) The (counter-clockwise) angle in the degrees needed to rotate from the x-axis in SWAN to the x-axis pointing East
      integer*4     :: nspectrumloc               = -123    !  [-] (advanced) Number of input spectrum locations
      integer*4     :: wbcversion                 = -123    !  [-] (advanced) Version of wave boundary conditions
      integer*4     :: nonhspectrum               = -123    !  [-] (advanced) Spectrum format for wave action balance of nonhydrostatic waves

      ! [Section] Flow boundary condition parameters
      integer         :: front                    = -123    !  [name] Switch for seaward flow boundary
      character(slen) :: front_str                =  ' '    !
      integer         :: left                     = -123    !  [name] Switch for lateral boundary at ny+1
      character(slen) :: left_str                 =  ' '    !
      integer         :: right                    = -123    !  [name] Switch for lateral boundary at 0
      character(slen) :: right_str                =  ' '    !
      integer         :: back                     = -123    !  [name] Switch for boundary at bay side
      character(slen) :: back_str                 =  ' '    !
      integer*4     :: ARC                        = -123    !  [-] (advanced) Switch for active reflection compensation at seaward boundary
      real*8        :: order                      = -123    !  [-] (advanced) Switch for order of wave steering, 1 = first order wave steering (short wave energy only), 2 = second oder wave steering (bound long wave corresponding to short wave forcing is added)
      integer*4     :: freewave                   = -123    !  [-] (advanced) Switch for free wave propagation 0 = use cg (default); 1 = use sqrt(gh) in instat = ts_2
      integer*4     :: carspan                    = -123    !  [-] (deprecated) Switch for Carrier-Greenspan test 0 = use cg (default); 1 = use sqrt(gh) in instat = ts_2 for c&g tests
      real*8        :: epsi                       = -123    !  [-] (advanced) Ratio of mean current to time varying current through offshore boundary
      integer*4     :: nc                         = -123    !  [-] (advanced) Smoothing distance for estimating umean (defined as nr of cells)
      integer       :: tidetype                   = -123    !  [name] (advanced) Switch for offfshore boundary, velocity boundary or instant water level boundary
      character(slen):: tidetype_str              =   ' '   !  [-] (advanced) Switch for offfshore boundary, velocity boundary or instant water level boundary (default)


      ! [Section] Tide boundary conditions
      character(slen):: zs0file                   = 'abc'   !  [file] Name of tide boundary condition series
      integer*4     :: tideloc                    = -123    !  [-] Number of corner points on which a tide time series is specified
      integer       :: paulrevere                 = -123    !  [name] Specifies tide on sea and land or two sea points if tideloc = 2

      character(slen):: paulrevere_str            =  ' '    !
      !      if tideloc =>2, then this indicates where the time series are to be
      !      applied. Input for tidal information to xbeach options (3):
      !      1. one tidal record --> specify tidal record everywhere
      !      2. two tidal records --> Need to specify keyword 'paulrevere'
      !      paulrevere=='land' implies to apply one tidal record to
      !      both sea corners and one tidal record to both land corners
      !      paulrevere=='sea' implies to apply the first tidal record
      !      (column 2 in zs0input.dat) to the (x=1,y=1) sea corner and
      !      the second tidal record (third column) to the (x=1,y=N) sea corner
      !      3. four tidal records --> Need to list tidal records in
      !      zs0input.dat in order of:
      !         (x=1,y=1)
      !         (x=1,y=N)
      !         (x=N,y=N)
      !         (x=N,y=1)
      !      NOTE:  clockwise from (1,1) corner

      ! [Section] Flow discharge boundary conditions
      integer*4     :: ndischarge                 = -123   !  [-] (advanced) Number of discharge locations
      integer*4     :: ntdischarge                = -123   !  [-] (advanced) Length of discharge time series
      character(slen):: disch_loc_file            = 'abc'  !  [file] (advanced) Name of discharge locations file
      character(slen):: disch_timeseries_file     = 'abc'  !  [file] (advanced) Name of discharge timeseries file

      ! [Section] Wave breaking parameters
      integer       :: break                      = -123    !  [name] Type of breaker formulation
      character(slen):: break_str                 =  ' '    !  [-] Type of breaker formulation (1=roelvink, 2=baldock, 3=roelvink adapted, 4=roelvink on/off breaking)
      real*8        :: gamma                      = -123    !  [-] Breaker parameter in Baldock or Roelvink formulation
      real*8        :: gamma2                     = -123    !  [-] End of breaking parameter in Roelvink Daly formulation
      real*8        :: alpha                      = -123    !  [-] (advanced) Wave dissipation coefficient in Roelvink formulation
      real*8        :: n                          = -123    !  [-] (advanced) Power in Roelvink dissipation model
      real*8        :: gammax                     = -123    !  [-] (advanced) Maximum ratio wave height to water depth
      real*8        :: delta                      = -123    !  [-] (advanced) Fraction of wave height to add to water depth
      real*8        :: fw                         = -123    !  [-] (advanced) Bed friction factor
      real*8        :: fwcutoff                   = -123    !  [-] Depth greater than which the bed friction factor is not applied
      integer*4     :: breakerdelay               = -123    !  [-] (advanced) Switch to enable breaker delay model
      integer*4     :: shoaldelay                 = -123    !  [-] (advanced,silent) Switch to enable shoaling delay
      real*8        :: facsd                      = -123    !  [-] (advanced,silent) fraction of the local wave length to use for shoaling delay depth
      real*8        :: facrun                     = -123    !  [-] (advanced,silent) calibration coefficient for short wave runup

      ! [Section] Roller parameters
      integer*4     :: roller                     = -123    !  [-] (advanced) Switch to enable roller model
      real*8        :: beta                       = -123    !  [-] (advanced) Breaker slope coefficient in roller model
      integer*4     :: rfb                        = -123    !  [-] (advanced) Switch to feed back maximum wave surface slope in roller energy balance, otherwise rfb = par%Beta

      ! [Section] Wave-current interaction parameters
      integer*4     :: wci                        = -123    !  [-] Turns on wave-current interaction
      real*8        :: hwci                       = -123    !  [m] (advanced) Minimum depth until which wave-current interaction is used
      real*8        :: cats                       = -123    !  [Trep] (advanced) Current averaging time scale for wci, in terms of mean wave periods

      ! [Section] Flow parameters
      integer       :: bedfriction                = -123    !  [name] Bed friction formulation
      character(slen):: bedfriction_str           =  ' '    !
      real*8        :: bedfriccoef                = -123    !  [-] Bed friction coefficient
      character(slen):: bedfricfile               = 'abc'   !  [file] Bed friction file (only valid with values of C)
      real*8        :: C                          = -123    !  [m^0.5s^-1] Chezy coefficient
      real*8        :: cf                         = -123    !  [-] (advanced) Friction coefficient flow
      real*8        :: nuh                        = -123    !  [m^2s^-1] Horizontal background viscosity
      real*8        :: nuhfac                     = -123    !  [-] (advanced) Viscosity switch for roller induced turbulent horizontal viscosity
      real*8        :: nuhv                       = -123    !  [-] (advanced,silent) Longshore viscosity enhancement factor, following Svendsen (?)
      integer*4     :: smag                       = -123    !  [-] (advanced) Switch for smagorinsky subgrid model for viscocity

      ! [Section] Coriolis force parameters
      real*8        :: wearth                     = -123    !  [hour^-1] (advanced) Angular velocity of earth calculated as: 1/rotation_time (in hours)
      real*8        :: lat                        = -123    !  [deg] (advanced) Latitude at model location  for computing Coriolis

      ! [Section] Wind parameters
      real*8        :: rhoa                       = -123    !  [kgm^-3] (advanced) Air density
      real*8        :: Cd                         = -123    !  [-] (advanced) Wind drag coefficient
      real*8        :: windv                      = -123    !  [ms^-1] Wind velocity, in case of stationary wind
      real*8        :: windth                     = -123    !  [deg] Nautical wind direction, in case of stationary wind
      character(slen) :: windfile                 = 'abc'   !  [file] Name of file with non-stationary wind data

      ! [Section] Groundwater parameters
      real*8        :: kx                         = -123    !  [ms^-1] (advanced) Darcy-flow permeability coefficient in x-direction
      real*8        :: ky                         = -123    !  [ms^-1] (advanced) Darcy-flow permeability coefficient in y-direction
      real*8        :: kz                         = -123    !  [ms^-1] (advanced) Darcy-flow permeability coefficient in z-direction
      real*8        :: dwetlayer                  = -123    !  [m] (advanced) Thickness of the top soil layer interacting more freely with the surface water
      real*8        :: aquiferbot                 = -123    !  [m] (advanced) Level of uniform aquifer bottom
      character(slen):: aquiferbotfile            = 'abc'   !  [file] (advanced) Name of the aquifer bottom file
      real*8        :: gw0                        = -123    !  [m] (advanced) Level initial groundwater level
      character(slen):: gw0file                   = 'abc'   !  [file] (advanced) Name of initial groundwater level file
      integer*4     :: gwnonh                     = -123    !  [-] (advanced) Switch to turn on or off non-hydrostatic pressure for groundwater
      integer*4     :: gwfastsolve                = -123    !  [-] (advanced,silent) Reduce full 2D non-hydrostatic solution to quasi-explicit in longshore direction
      integer       :: gwscheme                   = -123    !  [name] (advanced) Scheme for momentum equation
      character(slen):: gwscheme_str              =  ' '    !
      real*8        :: gwReturb                   = -123    !  [-] (advanced) Reynolds number for start of turbulent flow in case of gwscheme = turbulent
      integer       :: gwheadmodel                = -123    !  [name] (advanced) Model to use for vertical groundwater head
      character(slen):: gwheadmodel_str           =  ' '    !
      integer*4     :: gwhorinfil                 = -123    !  [-] (advanced) switch to include horizontal infiltration from surface water to groundwater

      ! [Section] Q3D sediment transport parameters
      real*8        :: vonkar                     = -123    !  [-] (advanced,silent) von Karman constant
      real*8        :: vicmol                     = -123    !  [-] (advanced,silent) molecular viscosity
      integer*4     :: kmax                       = -123    !  [-] (advanced,silent) Number of sigma layers in Quasi-3D model; kmax = 1 is without vertical structure of flow and suspensions
      real*8        :: sigfac                     = -123    !  [-] (advanced,silent) dsig scales with log(sigfac)

      ! [Section] Non-hydrostatic correction parameters
      integer       :: solver                     = -123    ! [name] (advanced) Solver used to solve the linear system
      character(slen):: solver_str                =  ' '    !
      integer*4     :: solver_maxit               = -123    ! [-] (advanced) Maximum number of iterations in the linear sip solver
      real*8        :: solver_acc                 = -123    ! [-] (advanced) accuracy with respect to the right-hand side used
      !     in the following termination criterion:
      !         ||b-Ax || < acc*||b||
      real*8        :: solver_urelax              = -123    ! [-] (advanced) Underrelaxation parameter
      real*8        :: kdmin                      = -123    ! [-] (advanced) Minimum value of kd (pi/dx > min(kd))
      real*8        :: dispc                      = -123    ! [?] (advanced) Coefficient in front of the vertical pressure gradient
      real*8        :: Topt                       = -123    ! [s] (advanced) Absolute period to optimize coefficient
      integer*4     :: nhbreaker                  = -123    ! [-] (advanced) Non-hydrostatic breaker model
      real*8        :: breakviscfac               = -123    ! [-] (advanced) Factor to increase viscosity during breaking
      real*8        :: maxbrsteep                 = -123    ! [-] (advanced) Maximum wave steepness criterium
      real*8        :: secbrsteep                 = -123    ! [-] (advanced) Secondary maximum wave steepness criterium
      real*8        :: reformsteep                = -123    ! [-] (advanced) Wave steepness criterium to reform after breaking
      real*8        :: breakvisclen               = -123    ! [-] (advanced) Ratio between local depth and length scale in extra breaking viscosity

      ! [Section] Bed composition parameters
      real*8        :: rhos                       = -123    !  [kgm^-3] Solid sediment density (no pores)
      integer*4     :: ngd                        = -123    !  [-] Number of sediment classes
      integer*4     :: nd                         = -123    !  [-] (advanced) Number of computational layers in the bed
      real*8        :: dzg1                       = -123    !  [m] (advanced) Thickness of top sediment class layers
      real*8        :: dzg2                       = -123    !  [m] (advanced) Nominal thickness of variable sediment class layer
      real*8        :: dzg3                       = -123    !  [m] (advanced) Thickness of bottom sediment class layers
      real*8        :: por                        = -123    !  [-] Porosity
      real*8,dimension(99)  :: D15                = -123    !  [m] D15 grain size per grain type
      real*8,dimension(99)  :: D50                = -123    !  [m] D50 grain size per grain type
      real*8,dimension(99)  :: D90                = -123    !  [m] D90 grain size per grain type
      real*8,dimension(99)  :: sedcal             = -123    !  [-] (advanced) Sediment transport calibration coefficient per grain type
      real*8,dimension(99)  :: ucrcal             = -123    !  [-] (advanced) Critical velocity calibration coefficient per grain type


      ! [Section] Sediment transport parameters
      integer       :: waveform                   = -123    !  [name] Wave shape model
      character(slen):: waveform_str              =  ' '    !
      integer       :: form                       = -123    !  [name] Equilibrium sediment concentration formulation
      character(slen) :: form_str                 =  ' '    !
      integer*4     :: sws                        = -123    !  [-] (advanced) Switch to enable short wave and roller stirring and undertow
      integer*4     :: lws                        = -123    !  [-] (advanced) Switch to enable long wave stirring
      real*8        :: BRfac                      = -123    !  [-] (advanced) Calibration factor surface slope
      real*8        :: facsl                      = -123    !  [-] (advanced) Factor bedslope effect
      real*8        :: z0                         = -123    !  [m] (advanced) Zero flow velocity level in Soulsby and van Rijn (1997) sediment concentration
      real*8        :: smax                       = -123    !  [-] (advanced) Maximum Shields parameter for equillibrium sediment concentration acc. Diane Foster
      real*8        :: tsfac                      = -123    !  [-] (advanced) Coefficient determining Ts = tsfac * h/ws in sediment source term
      real*8        :: facua                      = -123    !  [-] (advanced) Calibration factor time averaged flows due to wave skewness and asymmetry
      real*8        :: facSk                      = -123    !  [-] (advanced) Calibration factor time averaged flows due to wave skewness
      real*8        :: facAs                      = -123    !  [-] (advanced) Calibration factor time averaged flows due to wave asymmetry
      integer       :: turbadv                    = -123    !  [name] (advanced) Switch to activate turbulence advection model for short and or long wave turbulence
      character(slen):: turbadv_str               =  ' '    !
      integer       :: turb                       = -123    !  [name] (advanced) Switch to include short wave turbulence
      character(slen):: turb_str                  =  ' '    !
      real*8        :: Tbfac                      = -123    !  [-] (advanced) Calibration factor for bore interval Tbore: Tbore = Tbfac*Tbore
      real*8        :: Tsmin                      = -123    !  [s] (advanced) Minimum adaptation time scale in advection diffusion equation sediment
      integer*4     :: lwt                        = -123    !  [-] (advanced) Switch to enable long wave turbulence

      real*8        :: betad                      = -123    !  [-] (advanced) Dissipation parameter long wave breaking turbulence
      character(slen):: swtable                   = 'abc'   !  [-] (deprecated)Name of intra short wave assymetry and skewness table
      integer*4     :: sus                        = -123    !  [-] (advanced) Calibration factor for suspensions transports
      integer*4     :: bed                        = -123    !  [-] (advanced) Calibration factor for bed transports
      integer*4     :: bulk                       = -123    !  [-] (advanced) Switch to compute bulk transport rather than bed and suspended load separately
      real*8        :: facDc                      = -123    !  [-] (advanced) Option to control sediment diffusion coefficient
      real*8        :: jetfac                     = -123    !  [-] (advanced,silent) Option to mimic turbulence production near revetments
      integer*4     :: fallvelred                = -123     !  [-] Switch to reduce fall velocity for high concentrations
      integer*4     :: dilatancy                 = -123     !  [-] Switch to reduce critical shields number due dilatancy
      real*8        :: rheeA                     = -123     !  [-] A parameter in the Van Rhee expression
      real*8        :: pormax                    = -123     !  [-] Max porosity used in the experession of Van Rhee
      real*8        :: reposeangle               = -123     !  [deg] Angle of internal friction
      integer*4     :: bdslpeffmag               = -123     !  [name] Modify the magnitude of the sediment transport based on the bed slope, uses facsl
      integer*4     :: bdslpeffini               = -123     !  [name] Modify the critical shields parameter based on the bed slope
      integer*4     :: bdslpeffdir               = -123     !  [name] Modify the direction of the sediment transport based on the bed slope
      real*8        :: bdslpeffdirfac            = -123     !  [-] Calibration factor in the modification of the direction

      ! [Section] Morphology parameters
      real*8        :: morfac                     = -123    !  [-] Morphological acceleration factor
      integer*4     :: morfacopt                  = -123    !  [-] (advanced) Switch to adjusting output times for morfac
      real*8        :: morstart                   = -123    !  [s] Start time morphology, in morphological time
      real*8        :: morstop                    = -123    !  [s] Stop time morphology, in morphological time
      real*8        :: wetslp                     = -123    !  [-] Critical avalanching slope under water (dz/dx and dz/dy)
      real*8        :: dryslp                     = -123    !  [-] Critical avalanching slope above water (dz/dx and dz/dy)
      real*8        :: hswitch                    = -123    !  [m] (advanced) Water depth at which is switched from wetslp to dryslp
      real*8        :: dzmax                      = -123    !  [m/s/m] (advanced) Maximum bed level change due to avalanching
      integer*4     :: struct                     = -123    !  [-] Switch for enabling hard structures
      character(slen):: ne_layer                  = 'abc'   !  [file] Name of file containing depth of hard structure

      ! [Section] Output variables
      integer*4     :: timings                    = -123    !  [-] (advanced) Switch enable progress output to screen
      real*8        :: tstart                     = -123    !  [s] Start time of output, in morphological time
      real*8        :: tint                       = -123    !  [s] (deprecated) Interval time of global output (replaced by tintg)
      real*8        :: tintg                      = -123    !  [s] Interval time of global output
      real*8        :: tintp                      = -123    !  [s] Interval time of point and runup gauge output
      real*8        :: tintc                      = -123    !  [s] (advanced) Interval time of cross section output
      real*8        :: tintm                      = -123    !  [s] Interval time of mean, var, max, min output
      character(slen):: tsglobal                  = 'abc'   !  [-] (advanced) Name of file containing timings of global output
      character(slen):: tspoints                  = 'abc'   !  [-] (advanced) Name of file containing timings of point output
      character(slen):: tscross                   = 'abc'   !  [-] (advanced,silent) Name of file containing timings of cross section output
      character(slen):: tsmean                    = 'abc'   !  [-] (advanced) Name of file containing timings of mean, max, min and var output
      integer*4     :: nglobalvar                 = -123    !  [-] Number of global output variables (as specified by user)
      character(maxnamelen), dimension(numvars)   :: globalvars = 'abc' !  [-] (advanced) Mnems of global output variables, not per se the same size as nglobalvar (invalid variables, defaults)
      integer*4     :: nmeanvar                   = -123    !  [-] Number of mean, min, max, var output variables
      character(maxnamelen), dimension(numvars)   :: meanvars  = 'abc'  !  [-] (advanced) Mnems of mean output variables (by variables)
      integer*4     :: npointvar                  = -123    !  [-] Number of point output variables
      character(maxnamelen), dimension(numvars)   :: pointvars  = 'abc'  !  [-] (advanced) Mnems of point output variables (by variables)
      integer*4     :: npoints                    = -123    !  [-] Number of output point locations
      integer*4     :: nrugauge                   = -123    !  [-] Number of output runup gauge locations
      integer, dimension(:), pointer                     :: pointtypes => NULL()  !  [-] (advanced) Point types (0 = point, 1 = rugauge)
      real*8 ,dimension(:), pointer                      :: xpointsw => NULL()  ! (advanced) world x-coordinate of output points
      real*8 ,dimension(:), pointer                      :: ypointsw => NULL()  ! (advanced) world y-coordinate of output points

      integer*4     :: nrugdepth                  = -123    !  [-] (advanced) Number of depths to compute runup in runup gauge
      real*8,dimension(99) :: rugdepth            = -123    !  [m] (advanced) Minimum depth for determination of last wet point in runup gauge
      integer*4     :: ncross                     = -123    !  [-] (advanced,silent) Number of output cross sections
      integer        :: outputformat              = OUTPUTFORMAT_DEBUG  !  [name] (advanced) Output file format
      character(slen):: outputformat_str          = 'debug'    !
      character(slen):: ncfilename                = 'xboutput.nc' ! [file] (advanced) xbeach netcdf output file name
      integer        :: outputprecision           = -123    ! [name] switch between single and double precision output in NetCDF
      character(slen):: outputprecision_str       =  ' '    !

      ! Projection units (not to be used, only pass to output, this limit is too short for WKT....)
      ! This could be the proj4 string +init=epsg:28992
      ! [Section] Output projection
      character(slen)  :: projection              = ' '     !  [-] (advanced) projection string
      integer*4        :: rotate                  = -123    !  [-] Rotate output as postprocessing with given angle

      ! [Section] Drifters parameters
      integer*4     :: ndrifter                   = -123    !  [-] Number of drifers
      character(slen) :: drifterfile              = 'abc'   !  [file] Name of drifter data file

      ! [Section] Ship parameters
      character(slen) :: shipfile                 = 'abc'   !  [file] Name of ship data file
      integer*4       :: nship                    = -123    !  [-] (advanced) Number of ships

      ! [Section] Vegetation parameters
      character(slen) :: veggiefile                = 'abc'   !  [file] Name of vegetation species list file
      character(slen) :: veggiemapfile             = 'abc'   !  [file] Name of vegetation species map file
      integer*4       :: nveg                     = -123    !  [-] Number of vegetation species

      ! [Section] Wave numerics parameters
      integer       :: scheme                     = -123   !  [name] (advanced) Numerical scheme for wave propagation
      character(slen):: scheme_str                =  ' '    !  [-] (advanced) Use first-order upwind (upwind_1), second order upwind (upwind_2) or Lax-Wendroff (lax_wendroff)
      real*8        :: wavint                     = -123    !  [s] Interval between wave module calls (only in stationary wave mode)
      real*8        :: maxerror                   = -123    !  [m] (advanced) Maximum wave height error in wave stationary iteration
      integer*4     :: maxiter                    = -123    !  [-] (advanced) Maximum number of iterations in wave stationary

      ! [Section] Flow numerics parameters
      real*8        :: eps                        = -123    !  [m] Threshold water depth above which cells are considered wet
      real*8        :: eps_sd                     = -123    !  [m/s] Threshold velocity difference to determine conservation of energy head versus momentum
      real*8        :: umin                       = -123    !  [m/s] Threshold velocity for upwind velocity detection and for vmag2 in equilibrium sediment concentration
      real*8        :: hmin                       = -123    !  [m] Threshold water depth above which Stokes drift is included
      integer*4     :: secorder                   = -123    !  [-] (advanced) Use second order corrections to advection/non-linear terms based on MacCormack scheme
      integer*4     :: oldhu                      = -123    !  [-] (advanced,silent) Switch to enable old hu calculation

      ! [Section] Sediment transport numerics parameters
      real*8        :: thetanum                   = -123    !  [-] (advanced) Coefficient determining whether upwind (1) or central scheme (0.5) is used.
      integer*4     :: sourcesink                 = -123    !  [-] (advanced) Switch to enable source-sink terms to calculate bed level change rather than suspended transport gradients
      real*8        :: cmax                       = -123    !  [-] (advanced) Maximum allowed sediment concentration

      ! [Section] Bed update numerics parameters
      real*8        :: frac_dz                    = -123    !  [-] (advanced) Relative thickness to split time step for bed updating
      integer*4     :: nd_var                     = -123    !  [-] (advanced) Index of layer with variable thickness
      real*8        :: split                      = -123    !  [-] (advanced) Split threshold for variable sediment layer (ratio to nominal thickness)
      real*8        :: merge                      = -123    !  [-] (advanced) Merge threshold for variable sediment layer (ratio to nominal thickness)
      integer*4     :: nsetbathy                  = -123    !  [-] (advanced) Number of prescribed bed updates
      character(slen) :: setbathyfile             = 'abc'   !  [file] (advanced) Name of prescribed bed update file

      ! [Section] MPI parameters
      integer       :: mpiboundary                = -123   ! [name] (advanced) Fix mpi boundaries along y-lines, x-lines, use manual defined domains or find shortest boundary automatically
      character(slen)      :: mpiboundary_str     =  ' '    !
      integer*4     :: mmpi                       = -123    ! [-] (advanced) Number of domains in cross-shore direction when manually specifying mpi domains
      integer*4     :: nmpi                       = -123    ! [-] (advanced) Number of domains in alongshore direction when manually specifying mpi domains

      ! [Section] Constants, not read in params.txt
      real*8               :: px                  = -123    !  [-] Pi
      complex(kind(0.0d0)) :: compi               = -123    !  [-] Imaginary unit
      real*8               :: rhog8               = -123    !  [Nm^-3] 1/8*rho*g

      ! [Section] Variables, not read in params.txt
      real*8               :: dt                  = -123    !  [s] Computational time step, in hydrodynamic time
      real*8               :: t                   = -123    !  [s] Computational time, in hydrodynamic time
      real*8               :: tnext               = -123    !  [s] Next time point for output or wave stationary calculation, in hydrodynamic time

   end type parameters

contains

   subroutine all_input(par)

      use readkey_module
      use xmpi_module
      use filefunctions
      use logging_module

      implicit none

      type(parameters), intent(inout)                     :: par

      character(slen)                                     :: testc,line
      character(slen)                                     :: dummystring

      integer                                             :: filetype,mmax,nmax,ier,ic
      logical                                             :: comment
      logical                                             :: fe1,fe2,fe3

      logical, parameter :: toall = .true.
      readkey_inio = toall

      call writelog('sl','','Reading input parameters: ')
      !
      ! Check params.txt exists
      !
      call check_file_exist('params.txt')
      !
      !
      ! Physical processes
      call writelog('l','','--------------------------------')
      call writelog('l','','Physical processes: ')
      par%cyclic      = readkey_int ('params.txt','cyclic',        0,        0,     1)
      par%swave       = readkey_int ('params.txt','swave',         1,        0,     1)
      par%single_dir  = readkey_int ('params.txt','single_dir',    0,        0,     1)
      par%lwave       = readkey_int ('params.txt','lwave',         1,        0,     1)
      par%flow        = readkey_int ('params.txt','flow',          1,        0,     1)
      par%sedtrans    = readkey_int ('params.txt','sedtrans',      1,        0,     1)
      par%morphology  = readkey_int ('params.txt','morphology',    1,        0,     1)
      par%avalanching = readkey_int ('params.txt','avalanching',   par%morphology,0,1)
      par%nonh        = readkey_int ('params.txt','nonh',          0,        0,     1)
      par%gwflow      = readkey_int ('params.txt','gwflow',        0,        0,     1)
      par%q3d         = readkey_int ('params.txt','q3d',           0,        0,     1,silent=.true.)
      par%swrunup     = readkey_int ('params.txt','swrunup',       0,        0,     1,silent=.true.)
      par%ships       = readkey_int ('params.txt','ships',         0,        0,     1)
      if (par%ships == 0) par%nship = 0
      ! nship defined by shipfile
      par%bchwiz      = readkey_int ('params.txt','bchwiz',        0,        0,     1,silent=.true.)
      par%vegetation  = readkey_int ('params.txt','vegetation',    0,        0,     1)
      par%setbathy    = readkey_int ('params.txt','setbathy',      0,        0,     1)
      !
      !
      ! Grid parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Grid parameters: ')

      ! check gridform
      call setallowednames('xbeach',GRIDFORM_XBEACH,'delft3d',GRIDFORM_DELFT3D)
      call setoldnames('0','1')
      call parmapply('gridform',1,par%gridform,par%gridform_str)

      ! gridform switch
      if (par%gridform==GRIDFORM_XBEACH) then
         ! XBeach grid parameters
         par%xori  = readkey_dbl('params.txt','xori',  0.d0,   -1d9,      1d9)
         par%yori  = readkey_dbl('params.txt','yori',  0.d0,   -1d9,      1d9)
         par%alfa  = readkey_dbl('params.txt','alfa',  0.d0,   0.d0,   360.d0)
         par%nx    = readkey_int('params.txt','nx',     50,      2,     10000,required=.true.)
         par%ny    = readkey_int('params.txt','ny',      2,      0,     10000,required=.true.)
         par%posdwn= readkey_dbl('params.txt','posdwn', 1.d0,     -1.d0,     1.d0)
         if (par%setbathy .ne. 1) then
            par%depfile = readkey_name('params.txt','depfile',required=.true.)
            call check_file_exist(par%depfile)
            call check_file_length(par%depfile,par%nx+1,par%ny+1)
         else
            par%depfile = readkey_name('params.txt','depfile')
         endif
         par%vardx = readkey_int('params.txt','vardx',   0,      0,         1)

         if (par%vardx==0) then
            par%dx    = readkey_dbl('params.txt','dx',    -1.d0,   0.d0,      1d9,required=.true.)
            par%dy    = readkey_dbl('params.txt','dy',    -1.d0,   0.d0,      1d9,required=.true.)
         else
            par%dx    = readkey_dbl('params.txt','dx',    -1.d0,   0.d0,      1d9)       ! bas: not required, but can be used in superfast 1D (smagorinsky and timestep)
            par%dy    = readkey_dbl('params.txt','dy',    -1.d0,   0.d0,      1d9)
            par%xfile = readkey_name('params.txt','xfile')
            call check_file_exist(par%xfile)
            call check_file_length(par%xfile,par%nx+1,par%ny+1)

            par%yfile = readkey_name('params.txt','yfile')
            if (par%ny>0) then
               call check_file_exist(par%yfile)
               call check_file_length(par%yfile,par%nx+1,par%ny+1)
            end if

         endif
      elseif (par%gridform==GRIDFORM_DELFT3D) then
         par%depfile = readkey_name('params.txt','depfile',required=.true.)
         call check_file_exist(par%depfile)
         par%dx = -1.d0   ! Why?
         par%dy = -1.d0   ! Why?
         par%xyfile = readkey_name('params.txt','xyfile',required=.true.)
         call check_file_exist(par%xyfile)
         ! read grid properties from xyfile
         if (xmaster) then
            open(31,file=par%xyfile,status='old',iostat=ier)
            ! skip comment text in file...
            comment=.true.
            do while (comment .eqv. .true.)
               read(31,'(a)',iostat=ier)line
               if (ier .ne. 0) then
                  call report_file_read_error(par%xyfile)
               endif
               if (line(1:1)/='*') then
                  comment=.false.
               endif
            enddo
            ! Check if grid coordinates are Cartesian
            ic=scan(line,'Cartesian')
            if (ic<=0) then
               call writelog('ewsl','','Delft3D grid is not Cartesian')
               call halt_program
            endif
            ! read grid dimensions
            read(31,*,iostat=ier) mmax,nmax
            ! catch new grid format that  specifies missing value
            if (ier .ne. 0) then ! try reading the next line
               read(31,*,iostat=ier) mmax,nmax
            endif
            ! if still error, then XBeach cannot read this file
            if (ier .ne. 0) then
               call report_file_read_error(par%xyfile)
            endif
            close (31)
         endif
#ifdef USEMPI
         call xmpi_bcast(mmax,toall)
         call xmpi_bcast(nmax,toall)
#endif
         par%nx = mmax-1
         par%ny = nmax-1
         ! should we allow this input for gridform == 'Delft3D'
         par%xori  = readkey_dbl('params.txt','xori',  0.d0,   -1d9,      1d9)
         par%yori  = readkey_dbl('params.txt','yori',  0.d0,   -1d9,      1d9)
         par%alfa  = readkey_dbl('params.txt','alfa',  0.d0,   0.d0,   360.d0)
         par%posdwn= readkey_dbl('params.txt','posdwn', 1.d0,     -1.d0,     1.d0)
      endif
      par%thetamin = readkey_dbl ('params.txt','thetamin', -90.d0,    -180.d0,  180.d0,required=(par%swave==1))
      par%thetamax = readkey_dbl ('params.txt','thetamax',  90.d0,    -180.d0,  180.d0,required=(par%swave==1))
      par%dtheta   = readkey_dbl ('params.txt','dtheta',    10.d0,      0.1d0,   20.d0,required=(par%swave==1))
      par%thetanaut= readkey_int ('params.txt','thetanaut',    0,        0,     1)
      if (par%single_dir==1) then
         par%dtheta_s   = readkey_dbl ('params.txt','dtheta_s',    10.d0,      0.1d0,   20.d0,required=(par%swave==1))
      endif
      !
      !
      ! Model time parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Model time parameters: ')
      par%CFL     = readkey_dbl ('params.txt','CFL',     0.7d0,     0.1d0,      0.9d0)
      par%tstop   = readkey_dbl ('params.txt','tstop', 2000.d0,      1.d0, 1000000.d0,required=.true.)
      !
      !
      ! Physical constants
      call writelog('l','','--------------------------------')
      call writelog('l','','Physical constants: ')
      par%rho        = readkey_dbl ('params.txt','rho',       1025.0d0,  1000.0d0,  1040.0d0)
      par%g          = readkey_dbl ('params.txt','g',         9.81d0,    9.7d0,     9.9d0)
      par%depthscale = readkey_dbl ('params.txt','depthscale',1.0d0,     1.0d0,     200.d0)
      !
      !
      ! Initial conditions
      call writelog('l','','--------------------------------')
      call writelog('l','','Initial conditions: ')
      par%zsinitfile = readkey_name('params.txt','zsinitfile')
      if (par%zsinitfile==' ') then
         ! do nothing
      else
         call check_file_exist(par%zsinitfile)
         if (par%gridform==GRIDFORM_XBEACH) then ! nx and ny not known in case of Delft3D
            call check_file_length(par%zsinitfile,par%nx+1,par%ny+1)
         endif
      endif
      par%hotstartflow = readkey_int ('params.txt','hotstartflow',    0,        0,     1)
      !
      !
      ! Wave boundary condition parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Wave boundary condition parameters: ')
      call setallowednames('stat',       INSTAT_STAT,       &
      'bichrom',    INSTAT_BICHROM,    &
      'ts_1',       INSTAT_TS_1,       &
      'ts_2',       INSTAT_TS_2,       &
      'jons',       INSTAT_JONS,       &
      'swan',       INSTAT_SWAN,       &
      'vardens',    INSTAT_VARDENS,    &
      'reuse',      INSTAT_REUSE,      &
      'ts_nonh',    INSTAT_TS_NONH,    &
      'off',        INSTAT_OFF,        &
      'stat_table', INSTAT_STAT_TABLE, &
      'jons_table', INSTAT_JONS_TABLE)
      call setoldnames('0','1','2','3','4','5','6','7','8','9','40','41')
      call parmapply('instat',2,par%instat,par%instat_str, &
      required=(par%swave==1))

      if ( par%instat==INSTAT_JONS .or. &
      par%instat==INSTAT_SWAN .or. &
      par%instat==INSTAT_VARDENS.or. &
      par%instat==INSTAT_STAT_TABLE .or. &
      par%instat==INSTAT_JONS_TABLE &
      )then
         par%bcfile = readkey_name('params.txt','bcfile')
         call check_file_exist(par%bcfile)
         call checkbcfilelength(par%tstop,par%instat,par%bcfile,filetype)
         ! Only carried out on xmaster so:
#ifdef USEMPI
         call xmpi_bcast(filetype,toall)
#endif
      elseif (par%instat==INSTAT_REUSE) then
         ! See if this is reusing nonhydrostatic, or hydrostatic boundary conditions
         ! Note, check file length is done after recomputation of tstop due to morfacopt
         ! at the end of this subroutine.
         if (xmaster) then
            inquire(file='ebcflist.bcf',exist=fe1)
            inquire(file='qbcflist.bcf',exist=fe2)
            inquire(file='nhbcflist.bcf',exist=fe3)
         endif
#ifdef USEMPI
         call xmpi_bcast(fe1,toall)
         call xmpi_bcast(fe2,toall)
         call xmpi_bcast(fe3,toall)
#endif
         if (fe3 .and. .not. (fe1 .or. fe2)) then
            par%nonhspectrum = 1
            dummystring='nhbcflist.bcf'
            call checkbcfilelength(par%tstop,par%instat,dummystring,filetype,nonh=.true.)
         elseif (.not. fe3 .and. (fe1 .and. fe2)) then
            par%nonhspectrum = 0
            dummystring='ebcflist.bcf'
            call checkbcfilelength(par%tstop,par%instat,dummystring,filetype)
            dummystring='qbcflist.bcf'
            call checkbcfilelength(par%tstop,par%instat,dummystring,filetype)
         elseif (fe3 .and. (fe1 .or. fe2)) then
            call writelog('lswe','', &
            'If ''instat=reuse'' the model directory may not contain multiple boundary definition files.')
            call writelog('lswe','','Use either ebcflist.bcf/qbcflist.bcf, or nhbcflist.bcf')
            call halt_program
         elseif (.not. fe3 .and. .not. (fe1 .and. fe2)) then
            call writelog('lswe','', &
            'If ''instat=reuse'' the model directory may not contain sufficient boundary definition files.')
            if (.not. fe1) then
               call writelog('lswe','','Model currently missing ebcflist.bcf')
            elseif (.not. fe2) then
               call writelog('lswe','','Model currently missing qbcflist.bcf')
            endif
            call halt_program
         else
            call writelog('lswe','','If ''instat=reuse'' the model directory must contain boundary definition files.')
            call writelog('lswe','','Use either ebcflist.bcf/qbcflist.bcf, or nhbcflist.bcf')
            call halt_program
         endif
      else
         filetype=-1
      endif
      par%taper    = readkey_dbl ('params.txt','taper',   100.d0,      0.0d0, 1000.d0)
      par%nmax     = readkey_dbl ('params.txt','nmax',    0.8d0,       0.5d0, 1.d0)
      if (par%instat == INSTAT_STAT .or. par%single_dir==1) then
         par%nonhspectrum    = readkey_int ('params.txt','nonhspectrum', par%nonh,          0,       1 )
         par%Hrms  = readkey_dbl ('params.txt','Hrms',      1.d0,      0.d0,    10.d0)
         par%Tm01  = readkey_dbl ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
         par%Trep  = readkey_dbl ('params.txt','Trep',     par%Tm01,   1.d0,    20.d0)
         par%dir0  = readkey_dbl ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
         par%m     = readkey_int ('params.txt','m',        10,         2,      128)
      elseif (par%instat == INSTAT_BICHROM) then
         par%Hrms  = readkey_dbl ('params.txt','Hrms',      1.d0,      0.d0,    10.d0)
         par%Tm01  = readkey_dbl ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
         par%Trep  = readkey_dbl ('params.txt','Trep',     par%Tm01,   1.d0,    20.d0)
         par%Tlong = readkey_dbl ('params.txt','Tlong',    80.d0,     20.d0,   300.d0)
         par%dir0  = readkey_dbl ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
         par%m     = readkey_int ('params.txt','m',        10,         2,      128)
      elseif (par%instat == INSTAT_TS_1 .or. par%instat == INSTAT_TS_2) then
         par%Hrms  = readkey_dbl ('params.txt','Hrms',      1.d0,      0.d0,    10.d0)
         par%Tm01  = readkey_dbl ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
         par%Trep  = readkey_dbl ('params.txt','Trep',     par%Tm01,   1.d0,    20.d0)
         par%dir0  = readkey_dbl ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
         par%m     = readkey_int ('params.txt','m',        10,         2,      128)
         call check_file_exist('bc/gen.ezs')
      elseif (par%instat == INSTAT_TS_NONH) then
         par%Tm01  = readkey_dbl ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
         par%Trep  = readkey_dbl ('params.txt','Trep',     par%Tm01,   1.d0,    20.d0)
         call check_file_exist('boun_U.bcf')
      endif

      call setallowednames('neumann',LATERALWAVE_NEUMANN,'wavecrest',LATERALWAVE_WAVECREST,'cyclic',LATERALWAVE_CYCLIC)
      call setoldnames('0','1')
      call parmapply('lateralwave',1,par%lateralwave,par%lateralwave_str)



      ! TODO: fix
      !if (isSetParameter('params.txt','lateralwave')) then
      !   call setallowednames('neumann',LEFTWAVE_NEUMANN,'wavecrest',LEFTWAVE_WAVECREST)
      !   call setoldnames('0','1')
      !   call parmapply('leftwave',1,par%leftwave)
      !else
      !   if (isSetParameter('params.txt','leftwave') .or. isSetParameter('params.txt','rightwave')) then
      !      par%leftwave  = readkey_str('params.txt','leftwave','neumann',2,2,allowednames,oldnames)
      !      par%rightwave  = readkey_str('params.txt','rightwave','neumann',2,2,allowednames,oldnames)
      !      if (par%leftwave==par%rightwave) then
      !         par%lateralwave = par%leftwave
      !         call writelog('lsw','','LEFTWAVE and RIGHTWAVE parameters are deprecated.')
      !         call writelog('lsw','','Setting LATERALWAVE to ', trim(par%leftwave))
      !      else
      !         call writelog('lswe','','LEFTWAVE and RIGHTWAVE parameters are deprecated.')
      !         call writelog('lswe','','Left and Right wave boundary conditions must be equal')
      !         call writelog('lswe','','Use LATERALWAVE parameter to set wave boundary conditions')
      !         call halt_program
      !      endif
      !   else
      !      par%lateralwave  = readkey_str('params.txt','lateralwave','neumann',3,3,allowednames,oldnames)
      !   endif
      !endif
      !deallocate(allowednames,oldnames)
      !
      !
      ! Wave-spectrum boundary condition parameters
      if ( par%instat == INSTAT_JONS          .or.    &
      par%instat == INSTAT_SWAN          .or.    &
      par%instat == INSTAT_VARDENS       .or.    &
      par%instat == INSTAT_JONS_TABLE                ) then

         call writelog('l','','--------------------------------')
         call writelog('l','','Wave-spectrum boundary condition parameters: ')

         par%nonhspectrum    = readkey_int ('params.txt','nonhspectrum', par%nonh,          0,       1 )
         par%random          = readkey_int ('params.txt','random',       1,          0,          1       )
         par%fcutoff         = readkey_dbl ('params.txt','fcutoff',      0.d0,       0.d0,       40.d0   )
         par%nspr            = readkey_int ('params.txt','nspr',         0,          0,          1       )
         par%trepfac         = readkey_dbl ('params.txt','trepfac',      0.01d0,     0.d0,       1.d0    )
         if (par%nonhspectrum==1) then
            par%sprdthr         = readkey_dbl ('params.txt','sprdthr',      0.00d0,     0.d0,       1.d0    )
         else
            par%sprdthr         = readkey_dbl ('params.txt','sprdthr',      0.08d0,     0.d0,       1.d0    )
         endif
         par%oldwbc          = readkey_int ('params.txt','oldwbc',       0,          0,          1       ,silent=.true.)
         par%correctHm0      = readkey_int ('params.txt','correctHm0',   1,          0,          1       ,silent=.true.)
         par%oldnyq          = readkey_int ('params.txt','oldnyq',       0,          0,          1       ,silent=.true.)
         par%Tm01switch      = readkey_int ('params.txt','Tm01switch',   0,          0,          1       )

         if (filetype==0) then
            par%rt          = readkey_dbl('params.txt','rt',   min(3600.d0,par%tstop),    1200.d0,    7200.d0 )
            par%dtbc        = readkey_dbl('params.txt','dtbc',          1.0d0,      0.1d0,      2.0d0   )
         endif

         if (par%instat==INSTAT_SWAN) then
            par%dthetaS_XB  = readkey_dbl ('params.txt','dthetaS_XB',   0.0d0,      -360.d0,    360.0d0 )
         endif
         !                               wbcversion defaults to 3
         if (par%oldwbc==1) then
            par%wbcversion = 1
         else
            par%wbcversion = readkey_int ('params.txt','wbcversion', 3, 1, 3)
         endif
         if (par%wbcversion>2) then
            par%nspectrumloc    = readkey_int ('params.txt','nspectrumloc',   1,          1,       par%ny+1 )
         else
            par%nspectrumloc = 1
         endif
      endif
      !
      !
      ! Flow boundary condition parameters
      ! front
      call writelog('l','','--------------------------------')
      call writelog('l','','Flow boundary condition parameters: ')
      call setallowednames('abs_1d',    FRONT_ABS_1D,  &
      'abs_2d',    FRONT_ABS_2D,  &
      'wall',      FRONT_WALL,    &
      'wlevel',    FRONT_WLEVEL,  &
      'nonh_1d',   FRONT_NONH_1D, &
      'waveflume', FRONT_WAVEFLUME)
      call setoldnames('0','1','2','3','4','5')
      call parmapply('front',2,par%front,par%front_str)

      ! left and right
      call setallowednames('neumann',   LR_NEUMANN,    &
      'wall'   ,   LR_WALL,       &
      'no_advec',  LR_NO_ADVEC,   &
      'neumann_v', LR_NEUMANN_V)
      call setoldnames('0','1')
      call parmapply('left',1,par%left,par%left_str)

      call parmapply('right',1,par%right,par%right_str)

      ! back
      call setallowednames('wall',    BACK_WALL,     &
      'abs_1d',  BACK_ABS_1D,   &
      'abs_2d',  BACK_ABS_2D,   &
      'wlevel',  BACK_WLEVEL)
      call setoldnames('0','1','2','3')
      call parmapply('back',3,par%back,par%back_str)

      ! others
      par%ARC         = readkey_int ('params.txt','ARC',      1,              0,       1       )
      par%order       = readkey_dbl ('params.txt','order',    2.d0,           1.d0,    2.d0    )
      par%carspan     = readkey_int ('params.txt','carspan',  0,              0,       1       )      ! deprecated, added for backwards compatibility
      par%freewave    = readkey_int ('params.txt','freewave', par%carspan,    0,       1       )
      par%epsi        = readkey_dbl ('params.txt','epsi',     -1.d0,          -1.d0,   0.2d0   )
      par%nc          = readkey_int ('params.txt','nc',       par%ny+1,       1,       par%ny+1)

      call setallowednames('instant',   TIDETYPE_INSTANT,  &
      'velocity',  TIDETYPE_VELOCITY)
      call parmapply('tidetype',2,par%tidetype,par%tidetype_str)

      !
      !
      ! Tide boundary conditions
      call writelog('l','','--------------------------------')
      call writelog('l','','Tide boundary conditions: ')
      par%tideloc    = readkey_int ('params.txt','tideloc', 0,             0,      4)
      if (par%tideloc>0) then
         if (par%tideloc==2) then
            call setallowednames('land',    PAULREVERE_LAND,  &
            'sea',     PAULREVERE_SEA)
            call setoldnames('0','1')
            call parmapply('paulrevere',1,par%paulrevere,par%paulrevere_str)
         endif
         par%zs0file = readkey_name('params.txt','zs0file')
         call check_file_exist(par%zs0file)
      else
         par%zs0        = readkey_dbl ('params.txt','zs0',     0.0d0,     -5.d0,      5.d0)
      endif
      !
      !
      ! Discharge boundary conditions
      call writelog('l','','--------------------------------')
      call writelog('l','','Discharge boundary conditions: ')

      par%disch_loc_file          = readkey_name  ('params.txt','disch_loc_file'                       )
      par%disch_timeseries_file   = readkey_name  ('params.txt','disch_timeseries_file'                )

      par%ndischarge              = get_file_length(par%disch_loc_file                                 )
      par%ntdischarge             = get_file_length(par%disch_timeseries_file                          )

      par%ndischarge              = readkey_int   ('params.txt', 'ndischarge',  par%ndischarge,  0, 100)
      par%ntdischarge             = readkey_int   ('params.txt', 'ntdischarge', par%ntdischarge, 0, 100)

      if (par%ndischarge>0) then
         call check_file_exist(par%disch_loc_file)

         if (par%ntdischarge>0) then
            call check_file_exist(par%disch_timeseries_file)
         endif
      endif
      !
      !
      ! Wave breaking parameters
      par%beta     = readkey_dbl ('params.txt','beta',    0.10d0,     0.05d0,   0.3d0)
      if (par%swave==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Wave breaking parameters: ')
         call setallowednames('roelvink1',     BREAK_ROELVINK1,  &
         'baldock',       BREAK_BALDOCK,    &
         'roelvink2',     BREAK_ROELVINK2,  &
         'roelvink_daly', BREAK_ROELVINK_DALY,  &
         'janssen',       BREAK_JANSSEN)
         call setoldnames('1','2','3','4','5')
         if (par%instat == INSTAT_STAT .or. par%instat == INSTAT_STAT_TABLE) then
            call parmapply('break',2,par%break,par%break_str) ! default: baldock
         else
            call parmapply('break',3,par%break,par%break_str) ! default: roelvink2
         endif

         par%gamma    = readkey_dbl ('params.txt','gamma',   0.55d0,     0.4d0,     0.9d0)   !changed 28/11
         if (par%break == BREAK_ROELVINK_DALY) then
            par%gamma2   = readkey_dbl ('params.txt','gamma2',   0.3d0,     0.0d0,     0.5d0)
         endif
         par%alpha    = readkey_dbl ('params.txt','alpha',   1.0d0,     0.5d0,     2.0d0)
         par%n        = readkey_dbl ('params.txt','n',       10.0d0,     5.0d0,    20.0d0)   !changed 28/11
         par%gammax   = readkey_dbl ('params.txt','gammax',   2.d0,      .4d0,      5.d0)    !changed 28/11
         par%delta    = readkey_dbl ('params.txt','delta',   0.0d0,     0.0d0,     1.0d0)
         par%fw       = readkey_dbl ('params.txt','fw',       0.d0,   0d0,      1.0d0)
         par%fwcutoff = readkey_dbl ('params.txt','fwcutoff',  1000.d0,   0d0,      1000.d0)
         par%breakerdelay = readkey_int ('params.txt','breakerdelay',    1,   0,      1)
         par%shoaldelay = readkey_int ('params.txt','shoaldelay',    0,   0,      1,silent=.true.)
         if (par%shoaldelay==1) then
            par%facsd      = readkey_dbl ('params.txt','facsd',       1.d0,   0d0,      2.0d0)
         endif
         if (par%swrunup==1) then
            par%facrun     = readkey_dbl ('params.txt','facrun',      1.d0,   0d0,      2.0d0)
         endif
         !
         !
         ! Roller parameters
         call writelog('l','','--------------------------------')
         call writelog('l','','Roller parameters: ')
         par%roller   = readkey_int ('params.txt','roller',     1,        0,     1)
         par%rfb      = readkey_int ('params.txt','rfb',        0,        0,     1)
         !
         !
         ! Wave-current interaction parameters
         call writelog('l','','--------------------------------')
         call writelog('l','','Wave-current interaction parameters: ')
         par%wci      = readkey_int ('params.txt','wci',        0,        0,     1)
         par%hwci     = readkey_dbl ('params.txt','hwci',   0.1d0,   0.001d0,      1.d0)
         par%cats     = readkey_dbl ('params.txt','cats',   4.d0,     1.d0,      50.d0)
      endif
      !
      !
      ! Flow parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Flow parameters: ')

      call setallowednames('chezy',                     BEDFRICTION_CHEZY,  &
      'cf',                        BEDFRICTION_CF, &
      'white-colebrook',           BEDFRICTION_WHITE_COLEBROOK, &
      'manning',                   BEDFRICTION_MANNING, &
      'white-colebrook-grainsize', BEDFRICTION_WHITE_COLEBROOK_GRAINSIZE)
      call parmapply('bedfriction',1,par%bedfriction,par%bedfriction_str)
      ! Catch old exceptions
      if (.not. isSetParameter('params.txt','bedfriction')) then
         ! Catch really wierd old combinations
         if (isSetParameter('params.txt','bedfricfile')) then
            call writelog('lswe','','The use of keyword BEDFRICFILE without keyword BEDFRICTION is not allowed')
            call writelog('lswe','','Terminating simulation')
            call halt_program
         endif
         ! Else, continue as normal
         if (isSetParameter('params.txt','cf') .and. .not. isSetParameter('params.txt','C')) then
            par%bedfriction = BEDFRICTION_CF
            par%bedfriccoef = readkey_dbl ('params.txt','cf',      3.d-3, 1.d-3, 0.1d0)
            call writelog('lws','(a,a)','Warning: Specification of bed friction using parameter ''cf''', &
            ' will not be supported in future versions of XBeach')
            call writelog('lws','(a)','Use parameters ''bedfriction'' and ''bedfriccoef'' instead')
         elseif (isSetParameter('params.txt','C') .and. .not. isSetParameter('params.txt','cf')) then
            par%bedfriction = BEDFRICTION_CHEZY
            par%bedfriccoef = readkey_dbl ('params.txt','C',   55.d0 ,     20.d0,    100.d0)
            call writelog('lws','(a,a)','Warning: Specification of bed friction using parameter ''C''', &
            ' will not be supported in future versions of XBeach')
            call writelog('lws','(a)','Use parameters ''bedfriction'' and ''bedfriccoef'' instead')
         elseif (isSetParameter('params.txt','C') .and. isSetParameter('params.txt','cf')) then
            par%bedfriction = BEDFRICTION_CHEZY
            par%bedfriccoef = readkey_dbl ('params.txt','C',   55.d0 ,     20.d0,    100.d0)
            call writelog('lws','(a,a)','Warning: Specification of bed friction using parameters ''C'' and ''cf''', &
            ' will not be supported in future versions of XBeach')
            call writelog('lws','(a)','Use parameters ''bedfriction'' and ''bedfriccoef'' instead')
            call writelog('lws','(a)','Warning: C and cf both specified. C will take precedence')
         else
            ! Default to Chezy value
            par%bedfriction = BEDFRICTION_CHEZY
            par%bedfriccoef = readkey_dbl('params.txt','bedfriccoef', 55.d0, 20.d0, 100.d0)
         endif
      else
         if (par%bedfriction==BEDFRICTION_CHEZY .or. &
         par%bedfriction==BEDFRICTION_CF .or. &
         par%bedfriction==BEDFRICTION_MANNING .or. &
         par%bedfriction==BEDFRICTION_WHITE_COLEBROOK) then
            ! try to read a bed friction file
            par%bedfricfile = readkey_name('params.txt','bedfricfile')
            if (par%bedfricfile .ne. ' ') then
               call check_file_exist(par%bedfricfile)
               if (par%gridform==GRIDFORM_XBEACH) then
                  call check_file_length(par%bedfricfile,par%nx+1,par%ny+1)
               endif
               call writelog('lws','(a,a,a)','Warning: bed friction coefficient values from file ''',&
               trim(par%bedfricfile), &
               ''' will be used in computation')
            else
               if (par%bedfriction==BEDFRICTION_CHEZY) then
                  par%bedfriccoef = readkey_dbl('params.txt','bedfriccoef', 55.d0, 20.d0, 100.d0)
               elseif (par%bedfriction==BEDFRICTION_CF) then
                  par%bedfriccoef = readkey_dbl('params.txt','bedfriccoef', 3.d-3, 1.d-3,  0.1d0)
               elseif (par%bedfriction==BEDFRICTION_MANNING) then
                  par%bedfriccoef = readkey_dbl('params.txt','bedfriccoef', 0.02d0, 0.01d0, 0.05d0)
               elseif (par%bedfriction==BEDFRICTION_WHITE_COLEBROOK) then
                  par%bedfriccoef = readkey_dbl('params.txt','bedfriccoef', 0.01d0, 3.5d-5, 0.9d0)
               endif
            endif
         else
            ! other formulations do not require a bed friction coefficient
            par%bedfricfile = ''  ! empty string so doesn't go searching for file called 'abc'
            par%bedfriccoef = -999.d0
         endif
      endif
      par%nuh     = readkey_dbl ('params.txt','nuh',       0.1d0,     0.0d0,   1.0d0)
      par%nuhfac  = readkey_dbl ('params.txt','nuhfac',    1.0d0,     0.0d0,   1.0d0)
      par%nuhv    = readkey_dbl ('params.txt','nuhv',      1.d0,      1.d0,    20.d0,silent=.true.)
      par%smag    = readkey_int ('params.txt','smag',      1,         0,       1)
      !
      !
      ! Coriolis force parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Coriolis force parameters: ')
      par%wearth  = readkey_dbl ('params.txt','wearth',   1.d0/24.d0, 0.d0,    1.d0)
      par%lat     = readkey_dbl ('params.txt','lat',     0.d0,      -90.d0,   90.d0)
      !
      !
      ! Wind parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Wind parameters: ')
      par%rhoa    = readkey_dbl ('params.txt','rhoa',   1.25d0,     1.0d0,   2.0d0)
      par%Cd      = readkey_dbl ('params.txt','Cd',    0.002d0,  0.0001d0,  0.01d0)
      par%windfile = readkey_name('params.txt','windfile')
      if (par%windfile==' ') then
         par%windv   = readkey_dbl ('params.txt','windv',   0.0d0,     0.0d0, 200.0d0)
         par%windth  = readkey_dbl ('params.txt','windth', 270.0d0,  -360.0d0, 360.0d0)
      else
         call check_file_exist(par%windfile)
      endif
      !
      !
      ! Groundwater parameters
      if (par%gwflow==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Groundwater parameters: ')
         par%kx         = readkey_dbl ('params.txt','kx'        , 0.0001d0 , 0.00001d0, 0.1d0)
         par%ky         = readkey_dbl ('params.txt','ky'        , par%kx   , 0.00001d0, 0.1d0)
         par%kz         = readkey_dbl ('params.txt','kz'        , par%kx   , 0.00001d0, 0.1d0)
         par%dwetlayer  = readkey_dbl ('params.txt','dwetlayer' , 0.1d0    , 0.01d0     , 1.d0)
         par%aquiferbotfile = readkey_name('params.txt','aquiferbotfile')
         if (par%aquiferbotfile==' ') then
            !also read in groundwater.f90 which determines value
            par%aquiferbot = readkey_dbl('params.txt','aquiferbot',-10.d0,-100.d0,100.d0)
         else
            call check_file_exist(par%aquiferbotfile)
         endif
         par%gw0file = readkey_name('params.txt','gw0file')
         if (par%gw0file==' ') then
            par%gw0 = readkey_dbl('params.txt','gw0',0.d0,-5.d0,5.d0)
         else
            call check_file_exist(par%gw0file)
         endif
         par%gwnonh     = readkey_int ('params.txt','gwnonh',      0,           0,        1)
         if (par%gwnonh==1) then
            if (par%ny>2) then
               par%gwfastsolve = readkey_int ('params.txt','gwfastsolve',      0,    0,      1,silent=.true.)
            endif
         endif

         ! Type of momentum equation
         call setallowednames('laminar',    GWSCHEME_LAMINAR,  &
         'turbulent',  GWSCHEME_TURBULENT)
         call setoldnames('darcy','modflow')
         call parmapply('gwscheme',1,par%gwscheme,par%gwscheme_str)

         if (par%gwscheme==GWSCHEME_TURBULENT) then
            par%gwReturb    = readkey_dbl ('params.txt','gwReturb'   , 100.d0    , 1.d0     , 600.d0)
         endif
         call setallowednames('parabolic',       GWHEADMODEL_PARABOLIC,   &
         'exponential',     GWHEADMODEL_EXPONENTIAL)
         call parmapply('gwheadmodel',1,par%gwheadmodel,par%gwheadmodel_str)

         par%gwhorinfil = readkey_int ('params.txt','gwhorinfil',      0,           0,        1)
      endif
      !
      !
      ! Q3D sediment transport parameters
      if (par%q3d==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Q3D sediment transport parameters: ')
         par%vonkar  = readkey_dbl ('params.txt','vonkar',   0.4d0,     0.01d0,  1.d0)
         par%vicmol  = readkey_dbl ('params.txt','vicmol',   0.000001d0,   0.d0,    0.001d0)
         par%kmax    = readkey_int ('params.txt','kmax ',      1,           1,        1000)
         par%sigfac  = readkey_dbl ('params.txt','sigfac ',1.3d0,     0.00d0,   10.d0)
      endif
      !
      !
      ! Non-hydrostatic correction parameters
      if (par%nonh==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Non-hydrostatic correction parameters: ')
         call setallowednames('sip',       SOLVER_SIPP,  &
         'tridiag',   SOLVER_TRIDIAGG)
         call setoldnames('1','2')
         if (par%ny>2) then
            call parmapply('solver',1,par%solver,par%solver_str)  ! default: sip
         else
            call parmapply('solver',2,par%solver,par%solver_str)  ! default: tridiag
         endif
         if (par%solver==SOLVER_SIPP) then
            par%solver_maxit = readkey_int('params.txt','solver_maxit' ,30,1,1000)
            par%solver_acc   = readkey_dbl('params.txt','solver_acc' ,0.005d0,0.00001d0,0.1d0)
            par%solver_urelax= readkey_dbl('params.txt','solver_urelax' ,0.92d0,0.5d0,0.99d0)
         endif
         par%kdmin        = readkey_dbl('params.txt','kdmin' ,0.0d0,0.0d0,0.05d0)
         par%dispc        = readkey_dbl('params.txt','dispc' ,1.0d0,0.1d0,2.0d0)
         par%Topt         = readkey_dbl('params.txt','Topt',  10.d0, 1.d0, 20.d0)
         par%nhbreaker    = readkey_int('params.txt','nhbreaker' ,2,0,3)
         if (par%nhbreaker==1) then
            par%breakviscfac = readkey_dbl('params.txt','breakviscfac',1.5d0, 1.d0, 3.d0)
            par%maxbrsteep   = readkey_dbl('params.txt','maxbrsteep',0.6d0, 0.3d0, 0.8d0)
            par%reformsteep  = readkey_dbl('params.txt','reformsteep',0.25d0*par%maxbrsteep,0.d0,0.95d0*par%maxbrsteep)
         elseif (par%nhbreaker==2) then
            par%breakvisclen = readkey_dbl('params.txt','breakvisclen',1.0d0, 0.75d0, 3.d0)
            par%maxbrsteep   = readkey_dbl('params.txt','maxbrsteep',0.6d0, 0.3d0, 0.8d0)
            par%secbrsteep   = readkey_dbl('params.txt','secbrsteep',0.5d0*par%maxbrsteep,0.d0,0.95d0*par%maxbrsteep)
         elseif (par%nhbreaker==3) then
            !par%avis    = readkey_dbl ('params.txt','avis',       0.3d0,     0.0d0,   1.0d0)
            !par%maxbrsteep   = readkey_dbl('params.txt','maxbrsteep',0.6d0, 0.3d0, 0.8d0)
            !par%secbrsteep   = readkey_dbl('params.txt','secbrsteep',0.5d0*par%maxbrsteep,0.d0,0.95d0*par%maxbrsteep)
         endif
      endif
      !
      !
      ! Sediment transport parameters
      if (par%sedtrans==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Sediment transport parameters: ')

         call setallowednames('soulsby_vanrijn',    FORM_SOULSBY_VANRIJN,  &
         'vanthiel_vanrijn',   FORM_VANTHIEL_VANRIJN)
         call setoldnames('1','2')
         call parmapply('form',2,par%form, par%form_str)

         call setallowednames('ruessink_vanrijn',  WAVEFORM_RUESSINK_VANRIJN,  &
         'vanthiel',          WAVEFORM_VANTHIEL)
         call setoldnames('1','2')
         call parmapply('waveform',2,par%waveform,par%waveform_str)

         par%sws      = readkey_int ('params.txt','sws',           1,        0,     1)
         par%lws      = readkey_int ('params.txt','lws',           1,        0,     1)
         par%BRfac    = readkey_dbl ('params.txt','BRfac',    1.0d0,       0.d0, 1.d0)
         par%facsl    = readkey_dbl ('params.txt','facsl  ',  1.6d0,       0.d0, 1.6d0)
         par%z0       = readkey_dbl ('params.txt','z0     ',0.006d0,    0.0001d0,   0.05d0)
         par%smax     = readkey_dbl ('params.txt','smax',   -1.d0,    -1.d0,   3.d0)       !changed 28/11 and back 10/2
         par%tsfac    = readkey_dbl ('params.txt','tsfac',   0.1d0,    0.01d0,   1.d0)
         par%facua    = readkey_dbl ('params.txt','facua  ',0.10d0,    0.00d0,   1.0d0)
         par%facSk    = readkey_dbl ('params.txt','facSk  ',par%facua,    0.00d0,   1.0d0)
         par%facAs    = readkey_dbl ('params.txt','facAs  ',par%facua,    0.00d0,   1.0d0)
         call setallowednames('none',       TURBADV_NONE,        &
         'lagrangian', TURBADV_LAGRANGIAN,  &
         'eulerian',   TURBADV_EULERIAN)
         call parmapply('turbadv',1,par%turbadv,par%turbadv_str)

         call setallowednames('none',              TURB_NONE,           &
         'wave_averaged',     TURB_WAVE_AVERAGED,  &
         'bore_averaged',     TURB_BORE_AVERAGED)
         call setoldnames('0','1','2')
         call parmapply('turb',3,par%turb,par%turb_str)

         par%Tbfac    = readkey_dbl ('params.txt','Tbfac  ',1.0d0,     0.00d0,   1.0d0)
         par%Tsmin    = readkey_dbl ('params.txt','Tsmin  ',0.5d0,     0.01d0,   10.d0)
         par%lwt      = readkey_int ('params.txt','lwt    ',0,           0,            1)
         par%betad    = readkey_dbl ('params.txt','betad  ',1.0d0,     0.00d0,   10.0d0)
         par%sus      = readkey_int ('params.txt','sus    ',1,           0,            1)
         par%bed      = readkey_int ('params.txt','bed    ',1,           0,            1)
         par%bulk     = readkey_int ('params.txt','bulk   ',0,           0,            1)
         par%facDc    = readkey_dbl ('params.txt','facDc  ',1.0d0,     0.00d0,   1.0d0)
         par%jetfac   = readkey_dbl ('params.txt','jetfac ',0.0d0,     0.00d0,   1.0d0,silent=.true.)
         par%fallvelred          = readkey_int ('params.txt','fallvelred',    0,        0,     1)
         par%dilatancy           = readkey_int ('params.txt','dilatancy', 0,    0,     1)
         if (par%dilatancy==1) then
            par%rheeA               = readkey_dbl ('params.txt','rheeA',         0.75d0,   0.75d0, 2.d0) ! Between 3/4 and 1/(1-n), see paper Van Rhee (2010)
            par%pormax              = readkey_dbl ('params.txt','pormax',        0.5d0,    0.3d0, 0.6d0)
         endif
         par%reposeangle         = readkey_dbl ('params.txt','reposeangle',  30.d0,     0.d0,     45.d0)
         call setallowednames('none',              BDSLPEFFMAG_NONE,           &
         'roelvink_total',    BDSLPEFFMAG_ROELV_TOTAL,  &
         'roelvink_bed',      BDSLPEFFMAG_ROELV_BED, &
         'soulsby_total',     BDSLPEFFMAG_SOULS_TOTAL, &
         'soulsby_bed',       BDSLPEFFMAG_SOULS_BED)
         call setoldnames('0','1','2','3','4')
         call parmapply('bdslpeffmag',2,par%bdslpeffmag)

         call setallowednames('none',     BDSLPEFFINI_NONE,   &
         'total',    BDSLPEFFINI_TOTAL,  &
         'bed',      BDSLPEFFINI_BED)
         call setoldnames('0','1','2')
         call parmapply('bdslpeffini',1,par%bdslpeffini)

         call setallowednames('none',     BDSLPEFFDIR_NONE,   &
         'talmon',   BDSLPEFFDIR_TALMON)
         call setoldnames('0','1')
         call parmapply('bdslpeffdir',1,par%bdslpeffdir)
         if (par%bdslpeffdir>0) then
            par%bdslpeffdirfac      = readkey_dbl ('params.txt','bdslpeffdirfac',   1.d0,    0.d0,  2.d0)
         endif
      endif
      !
      !
      ! Bed composition parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Bed composition parameters: ')
      par%ngd      = readkey_int ('params.txt','ngd',        1,           1,        20)
      par%nd       = readkey_int ('params.txt','nd ',        3,           3,        1000)
      par%por      = readkey_dbl ('params.txt','por',    0.4d0,       0.3d0,    0.5d0)
      if (par%dilatancy==1) then
         par%D15      = readkey_dblvec('params.txt','D15',par%ngd,size(par%D15),0.00015d0,0.00001d0,0.0008d0) ! Lodewijk
      endif
      par%D50      = readkey_dblvec('params.txt','D50',par%ngd,size(par%D50),0.0002d0,0.00005d0,0.0008d0)
      par%D90      = readkey_dblvec('params.txt','D90',par%ngd,size(par%D90),0.0003d0,0.00010d0,0.0015d0)
      if (par%sedtrans==1) then
         par%rhos     = readkey_dbl ('params.txt','rhos',  2650d0,     2400.d0,  2800.d0)
         par%dzg1     = readkey_dbl ('params.txt','dzg',    0.1d0,      0.01d0,     1.d0)
         par%dzg1     = readkey_dbl ('params.txt','dzg1', par%dzg1,     0.01d0,     1.d0)
         par%dzg2     = readkey_dbl ('params.txt','dzg2', par%dzg1,     0.01d0,     1.d0)
         par%dzg3     = readkey_dbl ('params.txt','dzg3', par%dzg1,     0.01d0,     1.d0)
         !                              file,keyword,size read vector,size vector in par,default,min,max
         par%sedcal   = readkey_dblvec('params.txt','sedcal',par%ngd,size(par%sedcal),1.d0,0.d0,2.d0)
         par%ucrcal   = readkey_dblvec('params.txt','ucrcal',par%ngd,size(par%ucrcal),1.d0,0.d0,2.d0)
      endif
      !
      !
      ! Morphology parameters
      if (par%morphology==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Morphology parameters: ')
         par%morfac   = readkey_dbl ('params.txt','morfac', 1.0d0,        0.d0,  1000.d0)
         par%morfacopt= readkey_int ('params.txt','morfacopt', 1,        0,        1)
         par%morstart = readkey_dbl ('params.txt','morstart',120.d0,      0.d0, 10000000.d0)
         par%morstop  = readkey_dbl ('params.txt','morstop', par%tstop,      0.d0, 10000000.d0)
         par%wetslp   = readkey_dbl ('params.txt','wetslp', 0.3d0,       0.1d0,     1.d0)
         par%dryslp   = readkey_dbl ('params.txt','dryslp', 1.0d0,       0.1d0,     2.d0)
         par%hswitch  = readkey_dbl ('params.txt','hswitch',0.1d0,      0.01d0,    1.0d0)
         par%dzmax    = readkey_dbl ('params.txt','dzmax  ',0.05d0,    0.00d0,   1.0d0)
         par%struct   = readkey_int ('params.txt','struct ',0    ,      0,             1)
         if (par%struct==1) then
            par%ne_layer = readkey_name('params.txt','ne_layer')
            call check_file_exist(par%ne_layer)
            if (par%gridform==GRIDFORM_XBEACH) then
               call check_file_length(par%ne_layer,par%nx+1,par%ny+1)
            endif
         endif
         if (par%swrunup==1 .and. par%struct==0) then
            call writelog('lws','(a)','Warning: swrunup can only be used in combination with struct=1. swrunup will be turned off.')
            par%swrunup = 0;
         endif
      endif
      !
      !
      ! Output variables
      call writelog('l','','--------------------------------')
      call writelog('l','','Output variables: ')
      par%timings  = readkey_int ('params.txt','timings',      1,       0,      1)
      testc = readkey_name('params.txt','tunits')
      if (len(trim(testc)) .gt. 0) par%tunits = trim(testc)
      par%tstart  = readkey_dbl ('params.txt','tstart',   1.d0,      0.d0,1000000.d0)
      par%tint    = readkey_dbl ('params.txt','tint',     1.d0,     .01d0, 100000.d0)  ! Robert
      par%tsglobal = readkey_name('params.txt','tsglobal')
      if (par%tsglobal==' ') then
         par%tintg   = readkey_dbl ('params.txt','tintg', par%tint,     .01d0, 100000.d0)  ! Robert
      endif
      par%tspoints = readkey_name('params.txt','tspoints')
      if (par%tspoints==' ') then
         par%tintp   = readkey_dbl ('params.txt','tintp', par%tint,     .01d0, 100000.d0)  ! Robert
      endif
      par%tscross = readkey_name('params.txt','tscross',silent=.true.)
      if (par%tscross==' ') then
         par%tintc   = readkey_dbl ('params.txt','tintc', par%tint,     .01d0, 100000.d0,silent=.true.)  ! Robert
      endif
      par%tsmean = readkey_name('params.txt','tsmean')
      if (par%tsmean==' ') then
         par%tintm   = readkey_dbl ('params.txt','tintm', par%tstop-par%tstart,     1.d0, par%tstop-par%tstart)  ! Robert
      endif

      ! global output
      par%nglobalvar  = readkey_int ('params.txt','nglobalvar', -1, -1, 20)
      call readglobalvars(par)

      ! point output
      par%npoints     = readkey_int ('params.txt','npoints',     0,  0, 50)
      par%nrugauge    = readkey_int ('params.txt','nrugauge',    0,  0, 50)
      ! update the pointvariables
      !call readpointvars(par, par%xpointsw, par%ypointsw, par%pointtypes, par%pointvars)
      ! Robert: to deal with MPI some changes here
      par%npointvar   = readkey_int ('params.txt','npointvar',   0,  0, 50)
      call readpointvars(par)
      !    par%rugdepth    = readkey_dbl ('params.txt','rugdepth', 0.0d0,0.d0,0.05d0)
      par%nrugdepth   = readkey_int('params.txt','nrugdepth',1,1,10)
      par%rugdepth    = readkey_dblvec('params.txt','rugdepth',par%nrugdepth,size(par%rugdepth),0.0d0,0.0d0,0.1d0)

      ! mean output
      par%nmeanvar    = readkey_int ('params.txt','nmeanvar'  ,  0,  0, 15)
      call readmeans(par)
      par%ncross      = readkey_int ('params.txt','ncross',      0,  0, 50)

      call setallowednames('fortran',             OUTPUTFORMAT_FORTRAN,  &
      'netcdf ',             OUTPUTFORMAT_NETCDF,   &
      'debug  ',             OUTPUTFORMAT_DEBUG)
#ifdef USENETCDF 
      ! Default to NetCDF output in case of NetCDF-enabled executable
      call parmapply('outputformat',2,par%outputformat,par%outputformat_str, &
      required = .false.) ! wwvv-todo
#else     
      call parmapply('outputformat',1,par%outputformat,par%outputformat_str, &
      required = .false.) ! wwvv-todo
#endif
      if(par%outputformat==OUTPUTFORMAT_NETCDF .or. &
         par%outputformat==OUTPUTFORMAT_DEBUG) then
         ! type of output precision
         call setallowednames('single',  OUTPUTPRECISION_SINGLE,  &
                              'double',  OUTPUTPRECISION_DOUBLE)
         call parmapply('outputprecision',2,par%outputprecision,par%outputprecision_str,required = .false.)
         ! get the nc output file name from the parameter file
         par%ncfilename = readkey_name('params.txt','ncfilename')
         if (len(trim(par%ncfilename)) .eq. 0) par%ncfilename = 'xboutput.nc'
         call writelog('ls','','netcdf output to:' // par%ncfilename)
      endif         
      !
      !
      ! Output projection
      call writelog('l','','--------------------------------')
      call writelog('l','','Output projection: ')
      testc = readkey_name('params.txt','projection')
      if (len(trim(testc)) .gt. 0) par%projection = testc
      par%rotate = readkey_int ('params.txt','rotate',     1,  0, 1)
      !
      !
      ! Drifters parameters
      if (isSetParameter('params.txt','drifterfile')) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Drifters parameters: ')
         par%drifterfile = readkey_name  ('params.txt', 'drifterfile'                    )
         call check_file_exist(par%drifterfile)
         par%ndrifter    = get_file_length(par%drifterfile                               )
         par%ndrifter    = readkey_int   ('params.txt', 'ndrifter', par%ndrifter, 0, 50  )
      else
         par%ndrifter = 0
      endif
      !
      !
      ! Shipwaves parameters
      if (par%ships==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Shipwaves parameters: ')
         par%shipfile = readkey_name  ('params.txt', 'shipfile')
         call check_file_exist(par%shipfile)
         ! shipfile routine should set nship
      endif
      !
      !
      ! Vegetation parameters
      if (par%vegetation==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Vegetation parameters: ')
         par%veggiefile    = readkey_name  ('params.txt', 'veggiefile'                       )
         par%veggiemapfile = readkey_name  ('params.txt', 'veggiemapfile'                    )
         ! veggiefile routine should set nveg
      endif
      !
      !
      ! Wave numerics parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Wave numerics parameters: ')

      call setallowednames('upwind_1',      SCHEME_UPWIND_1,      &
      'lax_wendroff',  SCHEME_LAX_WENDROFF,  &
      'upwind_2',      SCHEME_UPWIND_2)
      call setoldnames('1','2','3')
      call parmapply('scheme',3,par%scheme,par%scheme_str)

      if (par%instat == INSTAT_STAT .or. par%instat == INSTAT_STAT_TABLE .or. par%single_dir==1) then
         par%wavint     = readkey_dbl ('params.txt','wavint',    60.d0,      1.d0,  3600.d0)
         par%maxerror   = readkey_dbl ('params.txt','maxerror', 0.00005d0, 0.00001d0, 0.001d0)
         par%maxiter    = readkey_int ('params.txt','maxiter',    500,         2,      1000)
      endif
      ! only default to Snell's Law if in 1D and only one directional bin
      if (par%ny==0 .and. nint(abs(par%thetamax-par%thetamin)/par%dtheta)<2) then
         par%snells      = readkey_int ('params.txt','snells',    1,        0,     1)
      else
         par%snells      = readkey_int ('params.txt','snells',    0,        0,     1)
      endif
      !
      !
      ! Flow numerics parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Flow numerics parameters: ')
      par%eps     = readkey_dbl ('params.txt','eps',     0.005d0,   0.001d0,      0.1d0)
      par%eps_sd  = readkey_dbl ('params.txt','eps_sd',  0.5d0,     0.000d0,      1.0d0)
      par%umin    = readkey_dbl ('params.txt','umin',    0.0d0,     0.0d0,        0.2d0)
      par%hmin    = readkey_dbl ('params.txt','hmin',    0.2d0,     0.001d0,      1.d0)
      par%secorder = readkey_int('params.txt','secorder' ,0,0,1)
      par%oldhu    = readkey_int('params.txt','oldhu' ,0,0,1,silent=.true.)
      !
      !
      ! Sediment transport numerics parameters
      if (par%sedtrans==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Sediment transport numerics parameters: ')
         par%thetanum   = readkey_dbl ('params.txt','thetanum',   1.d0,    0.5d0,   1.d0)
         par%sourcesink = readkey_int ('params.txt','sourcesink    ',0,     0,         1)
         par%cmax       = readkey_dbl ('params.txt','cmax',      0.1d0,    0.0d0,   1.d0)
      endif
      !
      !
      ! Bed update numerics parameters
      if (par%morphology==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Bed update numerics parameters: ')
         par%frac_dz   = readkey_dbl ('params.txt','frac_dz',   0.7d0,    0.5d0,   0.98d0)
         par%nd_var    = readkey_int ('params.txt','nd_var',     2,           2,        par%nd)
         par%split     = readkey_dbl ('params.txt','split',    1.01d0,  1.005d0,   1.10d0)
         par%merge     = readkey_dbl ('params.txt','merge',    0.01d0,  0.005d0,   0.10d0)
      endif
      !
      ! Prescribed bathy update
      if (par%setbathy==1) then
         par%nsetbathy    = readkey_int ('params.txt','nsetbathy',1,1,1000)
         par%setbathyfile = readkey_name  ('params.txt', 'setbathyfile',required=.true. )
         call check_file_exist(par%setbathyfile)
      endif
      !
      !
      ! MPI parameters
#ifdef USEMPI
      call writelog('l','','--------------------------------')
      call writelog('l','','MPI parameters: ')

      call setallowednames('auto',   MPIBOUNDARY_AUTO,   &
      'x',      MPIBOUNDARY_X,      &
      'y',      MPIBOUNDARY_Y,      &
      'man',    MPIBOUNDARY_MAN)
      call parmapply('mpiboundary',1,par%mpiboundary,par%mpiboundary_str)
      if (par%instat == INSTAT_STAT .or. par%instat == INSTAT_STAT_TABLE .or. par%single_dir==1) then
         par%mpiboundary=MPIBOUNDARY_X
         par%mpiboundary_str='x'
         call writelog('l','','mpiboundary set to x for stationary wave model')
      endif

      par%mmpi= readkey_int('params.txt','mmpi',2,1,100)
      par%nmpi= readkey_int('params.txt','nmpi',4,1,100)

#endif

      !
      !
      ! Finish
      call writelog('l','','--------------------------------')
      call writelog('sl','','Finished reading input parameters')
      call writelog('l','','--------------------------------')
      !
      !
      ! -------------------   Post-input processing -------------------------
      !
      !
      ! Fix input parameters for choosen depthscale
      if (par%depthscale .ne. 1.d0) then
         par%eps     = par%eps/par%depthscale
         par%hmin    = par%hmin/par%depthscale
         par%hswitch = par%hswitch/par%depthscale
         par%dzmax   = par%dzmax/par%depthscale**1.5d0

         call writelog('lws','(a)','Warning: input parameters eps, hmin, hswitch and dzmax are scaled with')
         call writelog('lws','(a)','         depthscale to:')
         call writelog('lws','(a,f0.4)','eps = ',    par%eps)
         call writelog('lws','(a,f0.4)','hmin = ',   par%hmin)
         call writelog('lws','(a,f0.4)','hswitch = ',par%hswitch)
         call writelog('lws','(a,f0.4)','dzmax = ',  par%dzmax)
      endif
      !
      !
      ! Constants
      par%px    = 4.d0*atan(1.d0)
      par%compi = (0.0d0,1.0d0)
      par%rhog8 = 1.0d0/8.0d0*par%rho*par%g
      !
      !
      if (par%posdwn<0.1d0) then
         par%posdwn=-1.d0  ! Backward compatibility, now posdwn = 0 also works in input (i.e. posdwn = false)
      endif
      !
      !
      ! Stop useless physical processes
      if (par%sedtrans==0 .and. par%morphology==1) then
         call writelog('lse','(a)','Error: Morphology cannot be computed without sediment transport.')
         call writelog('lse','(a)','       Set sedtrans=1 or morphology=0')
         call halt_program
      endif
      if (par%morphology==0 .and. par%avalanching==1) then
         call writelog('lsw','(a)','Warning: Avalanching cannot be computed without morphology.')
         call writelog('lsw','(a)','         Avalanching has been turned off')
         par%avalanching=0
      endif
      if (par%setbathy == 1 .and. par%morphology==1) then
         call writelog('lsw','(a)','Morphology and avalanching has been turned off for prescibed bed update (setbathy=1)')
         par%morphology=0
         par%avalanching=0
      endif
      !
      !
      ! Cyclic boundary conditions
#ifdef USEMPI
      if(par%cyclic==1) then
         call writelog('lsw','(a)','Warning: Cyclic boundary conditions only work in combination with MPI. Re-running this')
         call writelog('lsw','(a)','         simulation without MPI will not be possible')
      endif
#else
      if(par%cyclic==1) then
         call writelog('lswe','(a)','Error: Cyclic boundary conditions only work in combination with MPI.')
         call writelog('lsw','(a)','        Choose an MPI-compatable version of XBeach for this option')
         call halt_program
      endif
#endif
      !
      !
      ! Set taper to non-zero
      par%taper    = max(par%taper,1.d-6)
      !
      !
      ! Compute Coriolis
      par%lat = par%lat*par%px/180.d0
      par%wearth = par%px*par%wearth/1800.d0
      !
      !
      ! Only allow Baldock in stationary mode and Roelvink in non-stationary
      if (par%swave==1) then
         if (par%instat == INSTAT_STAT .or. par%instat == INSTAT_STAT_TABLE) then
            if (par%break .ne. BREAK_BALDOCK) then
               if(par%break == BREAK_ROELVINK_DALY) then
                  call writelog('lwse','','Error: Roelvink-Daly formulations not implemented in stationary wave mode, use Baldock')
                  call writelog('lwse','','         formulation.')
                  call halt_program
               else
                  call writelog('lwse','','Error: Roelvink formulations not allowed in stationary, use Baldock')
                  call writelog('lwse','','         formulation.')
                  call halt_program
               endif
            endif
         else
            if (par%break == BREAK_BALDOCK) then
               call writelog('lwse','','Error: Baldock formulation not allowed in non-stationary, use a Roelvink')
               call writelog('lwse','','         formulation.')
               call halt_program
            endif
         endif
      endif 
      !
      !
      ! Convert cf from C
      !par%cf      = par%g/par%C**2
      !
      !
      ! Only allow bore-averaged turbulence in combination with vanthiel waveform
      if ((par%waveform .ne. WAVEFORM_VANTHIEL) .and. (par%turb .eq. TURB_BORE_AVERAGED)) then
         call writelog('lse','','Error: Cannot compute bore-averaged turbulence without vanthiel wave form.')
         call writelog('lse','','       Please set waveform=vanthiel in params.txt, or choose another')
         call writelog('lse','','       turbulence model')
         call halt_program
      endif
      !
      !
      ! Set smax to huge if default is specified, but allow for computations to be
      !  done with smax  wwvv
      if (par%smax<0) par%smax=huge(0.d0)*1.0d-20
      !
      !
      ! Source-sink check
      if (par%morfac>1.d0) then
         if (par%sourcesink==1) then
            call writelog('lws','','Warning: Using source-sink terms for bed level change with morfac can lead to')
            call writelog('lws','','         loss of sediment mass conservation.')
         endif
      endif
      !
      !
      ! Check maximum grain size Soulsby-Van Rijn
      if ((par%form==FORM_SOULSBY_VANRIJN) .and. (maxval(par%D50(1:par%ngd))>= 0.05d0)) then
         call writelog('lews','','Error: Soulsby - Van Rijn cannot be used for sediment D50 > 0.05m')
         call halt_program
      endif
      !
      !
      ! If using tide, epsi should be on
      if (par%tideloc>0) then
         if (par%epsi<-1.d0) then
            call writelog('lws','','Automatically computing epsi using offshore boundary conditions')
            ! par%epsi = 0.05d0 --> Jaap do this in boundary conditions to account for vary wave conditions during simulation
         endif
      endif
      !
      !
      ! If using nonh, secorder should always be on
      if (par%nonh==1) then
         if (par%secorder==0) then
            call writelog('lws','','Warning: Automatically turning on 2nd order correction in flow for')
            call writelog('lws','','         non-hydrostatic module [secorder=1]')
            par%secorder = 1
         endif
      endif
      !
      !
      ! If using nonh, then the solver type is set by the grid size
      if (par%nonh==1) then
         if (par%ny>2 .and. par%solver==SOLVER_TRIDIAGG) then
            call writelog('lswe','','Tri-diagonal solver cannot be used if ny>2')
            call halt_program
         endif
         if (par%ny==0 .and. par%solver==SOLVER_SIPP) then
            call writelog('lse','','SIP solver cannot be used if ny==0')
            call halt_program
         endif
      endif
      !
      !
      ! If generating spectral time series for nonhydrostatic waves, you need at least wbcversion 3
      if (par%nonhspectrum==1 .and. par%wbcversion<3) then
         call writelog('lws','','Warning: Automatically changing to wbcversion 3 for')
         call writelog('lws','','         non-hydrostatic spectral boundary condition [nonhspectrum=1]')
         par%wbcversion=3
      endif
      if (par%wbcversion==1 .or. par%wbcversion==2) then
         call writelog('lwse','','************************** ERROR ******************************')
         call writelog('lwse','','wbcversion 1 and 2 are no longer supported, from v1.21 onwards')
         call writelog('lwse','','The current default wbcversion is 3')
         call writelog('lwse','','***************************************************************')
         call halt_program
      endif
      !
      !
      ! If using nhbreaker then maxbrsteep should be larger than reformsteep and reformsteep should be
      ! greater than zero
      if (par%nhbreaker==1) then
         if (par%reformsteep>0.95d0*par%maxbrsteep) then
            par%reformsteep=0.95d0*par%maxbrsteep
            call writelog('lws','(a)','Warning: Setting reformsteep to maximum of 95% of  maxbrsteep')
         elseif (par%reformsteep<0.0d0) then
            par%reformsteep=0.0d0
            call writelog('lws','(a)','Warning: Setting reformsteep to minimum of zero')
         endif
      endif
      !
      !
      ! The number of layers in the bed should be at least 3 (par%nd>=3)
      if (par%nd<3) then
         call writelog('lwse','','The number of sediment layers in the bed (nd) may not be smaller than 3')
         call halt_program
      endif
      !
      !
      ! fix minimum runup depth
      if (par%nrugauge>0) then
         ! Fill up remaining part of the array with minimum value
         par%rugdepth(par%nrugdepth+1:) = (1.0d0 + epsilon(0.d0))*par%eps
         if (minval(par%rugdepth)<=par%eps) then
            where(par%rugdepth<=par%eps)
               par%rugdepth = (1.0d0 + epsilon(0.d0))*par%eps
            endwhere
            call writelog('lws','(a,f0.5,a)','Warning: Setting rugdepth to minimum value greater than eps (', &
            (1.0d0 + epsilon(0.d0))*par%eps,')')
         endif
      endif
      !
      !
      ! Give an error if you ask for netcdf output if you don't have a netcdf executable
#ifndef USENETCDF
      if (par%outputformat .eq. OUTPUTFORMAT_NETCDF) then
         call writelog('lse', '', 'Error: You have asked for netcdf output [outputformat=netcdf] but this')
         call writelog('lse', '', '       executable is built without netcdf support. Use a netcdf enabled')
         call writelog('lse', '', '       executable or outputformat=fortran')
         call halt_program
      endif
#endif
      !
      !
      ! Lax-Wendroff not yet supported in curvilinear
      if (par%scheme==SCHEME_LAX_WENDROFF) then
         par%scheme=SCHEME_UPWIND_2
         call writelog('lws','','Warning: Lax Wendroff [scheme=lax_wendroff] scheme is not supported, changed')
         call writelog('lws','','         to 2nd order upwind [scheme=upwind_2]')
      endif
      !
      !
      ! Wave-current interaction with non-stationary waves still experimental
      if (par%instat/=INSTAT_STAT .and. par%instat/=INSTAT_STAT_TABLE .and. par%wci==1) then
         call writelog('lws','','Warning: Wave-current interaction with non-stationary waves is still')
         call writelog('lws','','         experimental, continue with computation nevertheless')
      endif
      !
      !
      ! 2D absorbing boundary limits to 1D absorbing boundary with 1D
      if (par%front==FRONT_ABS_2D .and. par%ny<3) then
         call writelog('lws','','Warning: 2D absorbing boundary condition [front=abs_2d] reduces to a')
         call writelog('lws','','         1D absorbing boundary condition [front=abs_1d] in')
         call writelog('lws','','         1D mode [ny=0]')
         par%front = FRONT_ABS_1D
      endif
      if (par%back==BACK_ABS_2D .and. par%ny<3) then
         call writelog('lws','','Warning: 2D absorbing boundary condition [back=abs_2d] reduces to a')
         call writelog('lws','','         1D absorbing boundary condition [back=abs_1d] in')
         call writelog('lws','','         1D mode [ny=0]')
         par%back = BACK_ABS_1D
      endif
      !
      !
      ! Mean output time fix
      if(par%tintm>(par%tstop-par%tstart) .and. par%nmeanvar>0) then
         call writelog('lws','','Warning: ''tintm'' is larger than output duration in the simulation.')
         call writelog('lws','','         Settting ''tintm'' = tstop-tstart = ',par%tstop-par%tstart,'s')
         par%tintm=par%tstop-par%tstart
      endif
      !
      !
      ! fix tint
      par%tint    = min(par%tintg,par%tintp,par%tintm,par%tintc)
      !
      !
      ! All input time frames converted to XBeach hydrodynamic time
      if (par%morfacopt==1) then
         par%tstart  = par%tstart / max(par%morfac,1.d0)
         par%tint    = par%tint   / max(par%morfac,1.d0)
         par%tintg   = par%tintg  / max(par%morfac,1.d0)
         par%tintp   = par%tintp  / max(par%morfac,1.d0)
         par%tintc   = par%tintc  / max(par%morfac,1.d0)
         par%tintm   = par%tintm  / max(par%morfac,1.d0)
         par%wavint  = par%wavint / max(par%morfac,1.d0)
         par%tstop   = par%tstop  / max(par%morfac,1.d0)
         par%morstart= par%morstart / max(par%morfac,1.d0)
         par%morstop = par%morstop / max(par%morfac,1.d0)
         par%rt      = par%rt / max(par%morfac,1.d0)
      endif
      !
      !
      ! Check bc file length in case of instat = reuse. In this case time is defined by
      ! morfacopt, which is not known at earlier stage
      if (par%instat == INSTAT_REUSE) then
         if (par%nonhspectrum==1) then
            dummystring='nhbcflist.bcf'
            call checkbcfilelength(par%tstop,par%instat,dummystring,filetype,nonh=.true.)
         else
            dummystring='ebcflist.bcf'
            call checkbcfilelength(par%tstop,par%instat,dummystring,filetype)
            dummystring='qbcflist.bcf'
            call checkbcfilelength(par%tstop,par%instat,dummystring,filetype)
         endif
      endif
      !
      !
      ! Check for unknown parameters
      if (xmaster) call readkey('params.txt','checkparams',dummystring)
      !
      !
      ! Distribute over MPI processes
#ifdef USEMPI
      call distribute_par(par)
#endif
      !
      !
      ! Write settings to file
      !include 'parameters.inc'
      !call outputparameters(par)

   end subroutine all_input


#ifdef USEMPI
   subroutine distribute_par(par)
      use mpi
      use xmpi_module
      implicit none
      type(parameters)        :: par
      integer                 :: ierror,i,npoints

      ! We're sending these over by hand, because intel fortran + vs2008 breaks things...
      integer, dimension(:), allocatable                     :: pointtypes !  [-] Point types (0 = point, 1=rugauge)
      real*8 ,dimension(:), allocatable                      :: xpointsw ! world x-coordinate of output points
      real*8 ,dimension(:), allocatable                      :: ypointsw ! world y-coordinate of output points

      logical, parameter                  :: toall = .true.
      integer                             :: parlen


      !
      ! distribute parameters

      ! This distributes all of the properties of par, including pointers. These point to memory adresses on the master
      ! We need to reset these on the non masters

      call xmpi_bcast(par%swave,toall)
      parlen = int(sizeof(par))

      if (toall) then
         call MPI_Bcast(par,parlen,MPI_BYTE,xmpi_imaster,xmpi_ocomm,ierror)
      else
         call MPI_Bcast(par,parlen,MPI_BYTE,xmpi_master,xmpi_comm,ierror)
      endif

      ! Ok now for the manual stuff to circumvent a bug in the intel compiler, which doesn't allow to send over arrays in derived types
      ! The only way to do it on all 3 compilers (gfortran, CVF, ifort) is with pointers.
      ! First let's store the number of variables, we need this to reserve some memory on all nodes

      do i=1,size(par%globalvars)
         call xmpi_bcast(par%globalvars(i),toall)
      enddo

      do i=1,size(par%pointvars)
         call xmpi_bcast(par%pointvars(i),toall)
      enddo

      do i=1,size(par%meanvars)
         call xmpi_bcast(par%meanvars(i),toall)
      enddo

      if (xmaster) npoints = size(par%pointtypes)
      ! send it over
      call xmpi_bcast(npoints,toall)

      ! now on all nodes allocate a array outside the par structure
      allocate(pointtypes(npoints))

      ! Par is only filled on the master, so use that one and put it in the seperate array
      if (xmaster) pointtypes = par%pointtypes

      ! Now for another ugly step, we can't broadcast the whole array but have to do it per variable.
      call xmpi_bcast(pointtypes,toall)

      ! so now everyone has the pointtypes, let's put it back in par
      ! first dereference the old one, otherwise we get a nasty error on ifort....
      if (.not. xmaster) par%pointtypes => NULL()

      ! now we can allocate the memory again
      if (.not. xmaster) allocate(par%pointtypes(npoints))

      ! and store the values from the local array
      if (.not. xmaster) par%pointtypes = pointtypes

      ! and clean up the local one.
      deallocate(pointtypes)

      if (xmaster) npoints = size(par%xpointsw)
      ! send it over
      call xmpi_bcast(npoints,toall)
      ! now on all nodes allocate a array outside the par structure
      allocate(xpointsw(npoints))
      ! Par is only filled on the master, so use that one and put it in the seperate array
      if (xmaster) xpointsw = par%xpointsw
      ! Now for another ugly step, we can't broadcast the whole array but have to do it per variable.
      call xmpi_bcast(xpointsw,toall)
      ! so now everyone has the xpointsw, let's put it back in par
      ! first dereference the old one, otherwise we get a nasty error on ifort....
      if (.not. xmaster) par%xpointsw => NULL()
      ! now we can allocate the memory again
      if (.not. xmaster) allocate(par%xpointsw(npoints))
      ! and store the values from the local array
      if (.not. xmaster) par%xpointsw = xpointsw
      ! and clean up the local one.
      deallocate(xpointsw)

      if (xmaster) npoints = size(par%ypointsw)
      ! send it over
      call xmpi_bcast(npoints,toall)
      ! now on all nodes allocate a array outside the par structure
      allocate(ypointsw(npoints))
      ! Par is only filled on the master, so use that one and put it in the seperate array
      if (xmaster) ypointsw = par%ypointsw
      ! Now for another ugly step, we can't broadcast the whole array but have to do it per variable.
      call xmpi_bcast(ypointsw,toall)
      ! so now everyone has the ypointsw, let's put it back in par
      ! first dereference the old one, otherwise we get a nasty error on ifort....
      if (.not. xmaster) par%ypointsw => NULL()
      ! now we can allocate the memory again
      if (.not. xmaster) allocate(par%ypointsw(npoints))
      ! and store the values from the local array
      if (.not. xmaster) par%ypointsw = ypointsw
      ! and clean up the local one.
      deallocate(ypointsw)




      return

      ! so, the following code is NOT used anymore. I left this here
      ! maybe method above does not work everywhere. wwvv

      ! For efficiency, this subroutine should use MPI_Pack and
      ! MPI_Unpack. However, since this subroutine is only called a
      ! few times, a more simple approach is used.
      !

      !call xmpi_bcast(par%px)
      !call xmpi_bcast(par%Hrms)
      !....

   end subroutine distribute_par
#endif

   !
   ! Some extra functions to make reading the output variables possible
   !
   subroutine readglobalvars(par)
      use logging_module
      use mnemmodule
      implicit none
      type(parameters), intent(inout)            :: par
      integer ::  i

      if (xmaster) then
         if (par%nglobalvar == -1) then
            par%globalvars(1:21) =  (/'H    ', 'zs   ', 'zs0  ', 'zb   ', 'hh   ', 'u    ', 'v    ', 'ue   ',&
            've   ', 'urms ', 'Fx   ', 'Fy   ', 'ccg  ', 'ceqsg', 'ceqbg', 'Susg ',&
            'Svsg ', 'E    ', 'R    ', 'D    ', 'DR   ' /)
            par%nglobalvar = 21
         elseif (par%nglobalvar == 999) then ! Output all
            par%nglobalvar = 0
            do i=1,numvars
               ! Don't output variables that don't work (misalignment)
               if (mnemonics(i) .eq. 'xyzs01') then
                  cycle
               elseif (mnemonics(i) .eq. 'xyzs02') then
                  cycle
               elseif (mnemonics(i) .eq. 'xyzs03') then
                  cycle
               elseif (mnemonics(i) .eq. 'xyzs04') then
                  cycle
               elseif (mnemonics(i) .eq. 'tideinpz') then
                  cycle
                  !elseif (mnemonics(i) .eq. 'umean') then
                  !   cycle
                  !elseif (mnemonics(i) .eq. 'vmean') then
                  !   cycle
               elseif (mnemonics(i) .eq. 'gw0back') then
                  cycle
               elseif (mnemonics(i) .eq. 'zi') then
                  cycle
               elseif (mnemonics(i) .eq. 'wi') then
                  cycle
               elseif (mnemonics(i) .eq. 'tideinpz') then
                  cycle
               end if
               par%globalvars(par%nglobalvar+1) = mnemonics(i)   ! list of all s% variables
               par%nglobalvar = par%nglobalvar + 1
            end do
         else
            ! User specified output
            ! Find nglobalvar keyword in params.txt
            call readOutputStrings(par,'global')
         end if ! globalvar
      end if ! xmaster
   end subroutine readglobalvars

   !
   ! FB:
   ! Now for a long one, reading the point and rugauges output options
   ! Just moved this from varoutput, so it can be used in ncoutput also
   ! Keeping as much as possible local to the subroutine...
   ! It can be reduced a bit by combining points and rugauges
   ! also instead of nvarpoints it can use a ragged array:
   ! For example: http://coding.derkeiler.com/Archive/Fortran/comp.lang.fortran/2004-05/0774.html
   !
   ! The outputformat for points is xworld, yworld, nvars, var1#var2#
   ! Split for nrugauges and npoints
   ! This is a bit inconvenient because you have different sets of points
   ! I think it's more logical to get value per variable than per point....
   ! So I read the data as follows:
   ! make a collection of all points
   ! store the types per point
   ! store all found variables in a combined list (not per point)
   ! So for each point all variables are stored
   !

   ! TODO after some discussion with Robert we will change this to:
   ! Rugauges don't need any output variables other then zs, x, y, t
   ! Pointvars need a different input specification:
   ! npointvars=3
   ! zs
   ! H
   ! u
   ! npoints=2
   ! 3.0 2.0
   ! 10.0 3.1
   ! nrugauges=1
   ! 4.1 3.2
   ! Using the old notation should give a warning and an explanation what will be output.

   ! Tasks:
   ! Create tests: FB+RMC
   ! Implement npointvars, stop using vars in #: RMC
   ! Give error if # present in points: RMC
   ! Use fixed outputvars for rugauges: RMC
   ! Fix matlab read routines (toolbox + zelt): FB
   ! Reimplement rugauges in ncoutput:  FB
   ! Update documentation: FB, check RMC


   subroutine readpointvars(par)
      use logging_module
      use mnemmodule
      use readkey_module
      implicit none
      type(parameters), intent(inout)          :: par

      real*8,dimension(:),allocatable          :: xpointsw,ypointsw
      integer, dimension(:), allocatable       :: pointtypes
      integer                                  :: i
      logical                                  :: xzfound, yzfound, zsfound

      ! These "targets" must be allocated by all processes
      allocate(pointtypes(par%npoints+par%nrugauge))
      allocate(xpointsw(par%npoints+par%nrugauge))
      allocate(ypointsw(par%npoints+par%nrugauge))
      if (xmaster) then
         ! set the point types to the indicator
         pointtypes(1:par%npoints)=0
         pointtypes(par%npoints+1:par%npoints+par%nrugauge)=1
         par%pointvars=''
         if (par%npoints>0) then
            call readPointPosition(par,'point',xpointsw,ypointsw)
            ! Is the keyword npointvar specified?
            if (isSetParameter('params.txt','npointvar',bcast=.false.)) then
               ! Look for keyword npointvar in params.txt
               call readOutputStrings(par,'point')
            else
               ! This branch of the else will change to a halt_program statement in later versions (written on 13 January 2011)
               call writelog('ls','','')
               call writelog('lws','','************************** WARNING  ***************************')
               call writelog('lws','','In future versions the keyword ''npointvar'' must be specified if ''npoints''>0')
               call writelog('lws','','Current order of output in all point output files is:')
               do i=1,par%npointvar
                  call writelog('lws','',trim(par%pointvars(i)))
               enddo
               call writelog('lws','','Order of point output variables stored in ''pointvars.idx''')
               call writelog('lws','','***************************************************************')
               call writelog('ls','','')
            endif  ! isSetParameter
         endif ! par%npoints>0

         if (par%nrugauge>0) then
            call readPointPosition(par,'rugauge',xpointsw,ypointsw)
            ! Make sure xz, yz and zs are in pointvars (default)
            xzfound = .false.
            yzfound = .false.
            zsfound = .false.
            do i=1,par%npointvar
               if (par%pointvars(i) .eq. 'xz') xzfound = .true.
               if (par%pointvars(i) .eq. 'yz') yzfound = .true.
               if (par%pointvars(i) .eq. 'zs') zsfound = .true.
            end do
            if (.not.xzfound) then
               par%pointvars(par%npointvar+1)='xz'
               par%npointvar = par%npointvar + 1
            end if
            if (.not.yzfound) then
               par%pointvars(par%npointvar+1)='yz'
               par%npointvar = par%npointvar + 1
            end if
            if (.not.zsfound) then
               par%pointvars(par%npointvar+1)='zs'
               par%npointvar = par%npointvar + 1
            end if
         end if

         ! Set pointers correct
         allocate(par%pointtypes(size(pointtypes)))
         allocate(par%xpointsw(size(xpointsw)))
         allocate(par%ypointsw(size(ypointsw)))
         par%pointtypes=pointtypes
         par%xpointsw= xpointsw
         par%ypointsw=ypointsw
      endif ! xmaster
      !
   end subroutine readpointvars


   subroutine readmeans(par)
      use logging_module
      use mnemmodule

      implicit none

      type(parameters),intent(inout)    :: par

      if (par%nmeanvar>0) then
         call readOutputStrings(par,'mean')
      endif
   end subroutine readmeans

   subroutine readOutputStrings(par,readtype)
      use logging_module
      use mnemmodule
      use readkey_module
      implicit none
      type(parameters), intent(inout)          :: par
      character(*), intent(in)                 :: readtype

      character(slen)                           :: okline,errline
      character(slen)                            :: line,keyword,keyread
      character(slen+slen)                     :: tempout
      integer                                  :: i,imax,id,ic,index,ier
      character(maxnamelen),dimension(numvars) :: tempnames

      imax = -123
      select case (readtype)
       case ('global')
         imax = par%nglobalvar
         keyword = 'nglobalvar'
         okline =  'nglobalvar: Will generate global output for variable: '
         errline = ' Unknown global output variable: '''
       case ('mean')
         imax = par%nmeanvar
         keyword = 'nmeanvar'
         okline =  'nmeanvar: Will generate mean, min, max and variance output for variable: '
         errline = ' Unknown mean output variable: '''
       case ('point')
         imax = par%npointvar
         keyword = 'npointvar'
         okline =  'npointvar: Will generate point output for variable: '
         errline = ' Unknown global output variable: '''
       case ('rugauge')
         imax = 0
         keyword = ''
         okline =  ''
         errline = ''
       case default
         write(*,*)'Programming error calling readOutputStrings'
         write(*,*)'Unknown calling type ''',trim(readtype),''''
         call halt_program
      end select
      tempnames=''

      if (xmaster) then
         id=0
         open(10,file='params.txt')   ! (this is done by xmaster only)
         do while (id == 0)
            read(10,'(a)',iostat=ier)line
            if (ier .ne. 0) then
               tempout = 'params.txt (looking for '//trim(keyword)//')'
               call report_file_read_error(tempout)
            endif
            ic=scan(line,'=')
            if (ic>0) then
               keyread=adjustl(line(1:ic-1))
               if (keyread == keyword) id=1
            endif
         enddo
         ! Read through the variables lines,
         do i=1,imax
            read(10,'(a)',iostat=ier)line
            if (ier .ne. 0) then
               tempout = 'params.txt (reading '//trim(keyword)//')'
               call report_file_read_error(tempout)
            endif
            line = line
            ! Check if this is a valid variable name
            index = chartoindex(line)
            if (index/=-1) then
               tempnames(i)=trim(line) ! wwvv use trim() to avoid compiler warning
               call writelog('ls','',trim(okline),trim(tempnames(i)))
            else
               call writelog('sle','',trim(errline),trim(line),'''')
               call halt_program
            endif
         end do
         close(10)
      endif

      ! only useful information on xmaster, but distributed later by distribute_pars
      select case (readtype)
       case ('global')
         par%globalvars=tempnames
       case ('mean')
         par%meanvars=tempnames
       case ('point')
         par%pointvars=tempnames
         call writelog('ls','','Order of point output variables stored in ''pointvars.idx''')
       case ('rugauge')
         ! no variables to store
      end select

   end subroutine readOutputStrings

   subroutine readPointPosition(par,readtype,xpoints,ypoints)
      use logging_module
      use mnemmodule
      use readkey_module
      implicit none
      type(parameters), intent(inout)          :: par
      character(*), intent(in)                 :: readtype
      real*8,dimension(:),intent(inout)        :: xpoints,ypoints

      character(slen)                          :: line,keyword,keyread,varline
      character(maxnamelen)                    :: varstr
      character(slen)                           :: fullline,errmes1,errmes2,okaymes
      integer                                  :: i,imax,id,ic,imark,imarkold,imin,nvar,ivar,index,j
      logical                                  :: varfound

      imin = 0
      imax = 0
      select case (readtype)
       case('point')
         imin = 0
         imax = par%npoints
         keyword = 'npoints'
         errmes1 = ' points '
         errmes2 = 'State point output variables using the "npointvar" keyword'
         okaymes = ' Output point '
       case('rugauge')
         imin = par%npoints
         imax = par%nrugauge
         keyword = 'nrugauge'
         errmes1 = ' runup gauge '
         errmes2 = 'Runup gauge automatically returns t,xz,yz and zs only'
         okaymes = ' Output runup gauge '
       case default
         write(*,*)'Programming error calling readOutputStrings'
         write(*,*)'Unknown calling type ''',trim(readtype),''''
         call halt_program
      end select

      if (xmaster) then
         id=0
         open(10,file='params.txt')   ! (this is done by xmaster only)
         do while (id == 0)
            read(10,'(a)')line
            ic=scan(line,'=')
            if (ic>0) then
               keyread=adjustl(line(1:ic-1))
               if (keyread == keyword) id=1
            endif
         enddo
         ! Read positions
         do i=1,imax
            read(10,'(a)')fullline
            ! Check params.txt has old (unsupported) method of defining points
            ic=scan(fullline,'#')
            if (ic .ne. 0) then ! This branch of the if/else will change to a halt_program statement in later versions
               ! (written on 13 January 2011)
               if (.not. isSetParameter('params.txt','npointvar',bcast=.false.)) then
                  read(fullline,*)xpoints(i+imin),ypoints(i+imin),nvar,varline
                  if (readtype=='point') then   ! no use at all for rugauge output, which is fixed at x,y,zs
                     imarkold = 0
                     do ivar=1,nvar
                        imark=scan(varline(imarkold+1:80),'#')
                        imark=imark+imarkold
                        varstr=varline(imarkold+1:imark-1)
                        index = chartoindex(varstr)
                        if (index/=-1) then
                           ! see if we already found this variable....
                           varfound = .false.
                           do j=1,numvars
                              if (par%pointvars(j) == varstr) then
                                 varfound = .true.
                              end if
                           end do
                           ! we have a new variable, store it
                           if (varfound .eqv. .false.) then
                              par%npointvar=par%npointvar+1
                              par%pointvars(par%npointvar)=varstr
                           end if
                        else
                           call writelog('sle','',' Unknown point output variable: ''',trim(varstr),'''')
                           call halt_program
                        endif
                        imarkold=imark
                     enddo
                  endif ! type point
               else   ! isset npointvar
                  read(fullline,*)xpoints(i+imin),ypoints(i+imin),nvar,varline
                  call writelog('lws','','Point output variables specified by ''npointvar'' will be selected over')
                  call writelog('lws','','variables specified on the point location line in params.txt')
               endif  ! isset npointvar
               call writelog('ls','','')
               call writelog('lws','','************************** WARNING  ***************************')
               call writelog('lws','','Unsupported method of defining output',trim(errmes1),'variables')
               call writelog('lws','',trim(fullline))
               call writelog('lws','','Please remove "nvar var1#Var2#..." from definition of',trim(errmes1))
               call writelog('lws','',trim(errmes2))
               call writelog('lws','','This warning will become and error in future versions of XBeach')
               call writelog('lws','','Refer to manual for complete documentation')
               call writelog('lws','','***************************************************************')
               call writelog('ls','','')
            else ! not old method of setting point output
               read(fullline,*)xpoints(i+imin),ypoints(i+imin)
            endif ! old method of point output
            call writelog('ls','(a,i0,a,f0.2,a,f0.2)',&
            trim(okaymes),i,' xpoint: ',xpoints(i+imin),'   ypoint: ',ypoints(i+imin))
         enddo
         close(10)
      endif

   end subroutine readPointPosition

end module params
