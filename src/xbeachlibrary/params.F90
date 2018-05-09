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
      !  type                                name                     initialize           !  [unit] (advanced/deprecated) description
      ! [Section] Use of default sets
      integer                           :: useXBeachGSettings       = -123                 !  [-] (advanced) Use XBeach-G default settings for all non-specified parameters
      
      ! [Section] Physical processes
      integer                             :: wavemodel              = -123                  !  [-] Stationary (0), surfbeat (1) or non-hydrostatic (2)
      character(slen)                     :: wavemodel_str          = ' '                  !  
      integer                             :: cyclic                 = -123                 !  [-] Turn on cyclic boundary conditions
      integer                             :: swave                  = -123                 !  [-] Turn on short waves
      integer                             :: lwave                  = -123                 !  [-] Turn on short wave forcing on NLSW equations and boundary conditions
      integer                             :: flow                   = -123                 !  [-] Turn on flow calculation
      integer                             :: sedtrans               = -123                 !  [-] Turn on sediment transport
      integer                             :: morphology             = -123                 !  [-] Turn on morphology
      integer                             :: avalanching            = -123                 !  [-] Turn on avalanching
      integer                             :: nonh                   = -123                 !  [-] (advanced) Turn on non-hydrostatic pressure: 0 = NSWE, 1 = NSW + non-hydrostatic pressure compensation Stelling & Zijlema, 2003
      integer                             :: gwflow                 = -123                 !  [-] (advanced) Turn on groundwater flow
      integer                             :: q3d                    = -123                 !  [-] (advanced,silent) Turn on quasi-3D sediment transport
      integer                             :: swrunup                = -123                 !  [-] (advanced,silent) Turn on short wave runup
      integer                             :: ships                  = -123                 !  [-] (advanced) Turn on ship waves
      integer                             :: vegetation             = -123                 !  [-] (advanced) Turn on interaction of waves and flow with vegetation
      integer                             :: snells                 = -123                 !  [-] (advanced) Turn on Snell's law for wave refraction
      integer                             :: single_dir             = -123                 !  [-] (advanced) Turn on stationary model for refraction, surfbeat based on mean direction
      integer                             :: bchwiz                 = -123                 !  [-] (advanced,silent) Turn on beachwizard
      integer                             :: setbathy               = -123                 !  [-] Turn on timeseries of prescribed bathy input
      integer                             :: viscosity              = -123                 !  [-] Include viscosity in flow solver
      integer                             :: advection              = -123                 !  [-] Include advection in flow solver
      integer                             :: wind                   = -123                 !  [-] Include wind in flow solver

      ! [Section] Grid parameters
      character(slen)                     :: depfile                = 'abc'                !  [file] Name of the input bathymetry file
      double precision                    :: posdwn                 = -123                 !  [-] Bathymetry is specified positive down (1) or positive up (-1)
      integer                             :: nx                     = -123                 !  [-] Number of computational cell corners in x-direction
      integer                             :: ny                     = -123                 !  [-] Number of computational cell corners in y-direction
      integer                             :: nz                     = -123                 !  [-] Number of computational cells in z-direction
      double precision                    :: alfa                   = -123                 !  [deg] Angle of x-axis from East
      integer                             :: vardx                  = -123                 !  [-] Switch for variable grid spacing
      double precision                    :: dx                     = -123                 !  [m] Regular grid spacing in x-direction
      double precision                    :: dy                     = -123                 !  [m] Regular grid spacing in y-direction
      character(slen)                     :: xfile                  = 'abc'                !  [file] Name of the file containing x-coordinates of the calculation grid
      character(slen)                     :: yfile                  = 'abc'                !  [file] Name of the file containing y-coordinates of the calculation grid
      double precision                    :: xori                   = -123                 !  [m] X-coordinate of origin of axis
      double precision                    :: yori                   = -123                 !  [m] Y-coordinate of origin of axis
      double precision                    :: thetamin               = -123                 !  [deg] Lower directional limit (angle w.r.t computational x-axis)
      double precision                    :: thetamax               = -123                 !  [deg] Higher directional limit (angle w.r.t computational x-axis)
      double precision                    :: dtheta                 = -123                 !  [deg] Directional resolution
      double precision                    :: dtheta_s               = -123                 !  [deg] Directional resolution in case of stationary refraction
      integer                             :: thetanaut              = -123                 !  [-] Switch to specify thetamin and thetamax in nautical convention rather than cartesian
      integer                             :: gridform               = -123                 !  [name] Grid definition format
      character(slen)                     :: xyfile                 = 'abc'                !  [file] Name of the file containing Delft3D xy-coordinates of the calculation grid
      character(slen)                     :: gridform_str           = ' '                  !  [name] Grid definition format

      ! [Section] Model time
      double precision                    :: tstop                  = -123                 !  [s] Stop time of simulation, in morphological time
      double precision                    :: CFL                    = -123                 !  [-] Maximum Courant-Friedrichs-Lewy number
      double precision                    :: dtset                  = -123                 !  [-] Fixed timestep, overrides use of CFL
      integer                             :: defuse                 = -123                 !  [-] (advanced,silent) Turn on timestep explosion prevention mechanism
      double precision                    :: maxdtfac               = -123                 !  [-] (advanced,silent) Maximum increase/decrease in time stp in explosion prevention mechanism
      character(slen)                     :: tunits                 = 's'                  !  [-] (advanced) Time units in udunits format (seconds since 1970-01-01 00:00:00.00 +1:00)

      ! [Section] Physical constants
      double precision                    :: g                      = -123                 !  [ms^-2] Gravitational acceleration
      double precision                    :: rho                    = -123                 !  [kgm^-3] Density of water
      double precision                    :: depthscale             = -123                 !  [-] (advanced)  depthscale of (lab)test simulated, affects eps, hmin, hswitch and dzmax
      
      ! [Section] Initial conditions
      double precision                  :: zs0                      = -123                 !  [m] Inital water level
      character(slen)                   :: zsinitfile               = 'abc'                !  [file] Name of inital water level file
      integer                           :: hotstartflow             = -123                 !  [-] (advanced) Switch for hotstart flow conditions with pressure gradient balanced by wind and bed stress
      
      ! [Section] Wave boundary condition parameters
      integer                           :: wbctype                  = -123                 !  [name] New wave boundary condition type
      character(slen)                   :: wbctype_str              = ' '                  !  [-] 
      integer                           :: instat                   = -123                 !  [name] Old wave boundary condition type
      character(slen)                   :: instat_str               = ' '                  !  [-] Wave boundary condition type
      double precision                  :: taper                    = -123                 !  [s] Spin-up time of wave boundary conditions, in morphological time
      double precision                  :: Hrms                     = -123                 !  [m] Hrms wave height for instat = stat, bichrom, ts_1 or ts_2
      double precision                  :: Tm01                     = -123                 !  [s] (deprecated) Old name for Trep
      double precision                  :: Trep                     = -123                 !  [s] Representative wave period for instat = stat, bichrom, ts_1 or ts_2
      double precision                  :: Tlong                    = -123                 !  [s] Wave group period for case instat = bichrom
      double precision                  :: dir0                     = -123                 !  [deg] Mean wave direction for instat = stat, bichrom, ts_1 or ts_2 (nautical convention)
      double precision                  :: nmax                     = -123                 !  [-] (advanced) maximum ratio of cg/c for computing long wave boundary conditions
      integer                           :: m                        = -123                 !  [-] Power in cos^m directional distribution for instat = stat, bichrom, ts_1 or ts_2
      integer                           :: lateralwave              = -123                 !  [name] Switch for lateral boundary at left
      character(slen)                   :: lateralwave_str          = ' '
      integer                           :: leftwave                 = -123                 !  [-] (deprecated) old name for lateralwave
      character(slen)                   :: leftwave_str             =  ' '                 !  [-] old name for lateralwave
      integer                           :: rightwave                = -123                 !  [-] (deprecated) old name for lateralwave
      character(slen)                   :: rightwave_str            = ' '                  !  [-] old name for lateralwav
      integer                           :: bclwonly                 = -123                 !  [-] (advanced,silent) switch to run boundary conditions with long waves only
      integer                           :: Sfold                    = -123                 !  [-] (advanced,silent) switch to run secoond order boundary conditions with Sf instead of Sfp
      
      ! [Section] Wave-spectrum boundary condition parameters
      character(slen)                   :: bcfile                   = 'abc'                !  [file] Name of spectrum file
      integer                           :: random                   = -123                 !  [-] (advanced) Switch to enable random seed for instat = jons, swan or vardens boundary conditions
      double precision                  :: fcutoff                  = -123                 !  [Hz] (advanced) Low-freq cutoff frequency for instat = jons, swan or vardens boundary conditions
      integer                           :: nspr                     = -123                 !  [-] (advanced,silent) Switch to enable long wave direction forced into centres of short wave bins
      double precision                  :: trepfac                  = -123                 !  [-] (advanced) Compute mean wave period over energy band: par%trepfac*maxval(Sf) for instat jons, swan or vardens; converges to Tm01 for trepfac = 0.0 and
      double precision                  :: sprdthr                  = -123                 !  [-] (advanced) Threshold ratio to maximum value of S above which spectrum densities are read in
      integer                           :: correctHm0               = -123                 !  [-] (advanced,silent) Switch to enable Hm0 correction
      integer                           :: Tm01switch               = -123                 !  [-] (advanced) Switch to enable Tm01 rather than Tm-10
      double precision                  :: rt                       = -123                 !  [s] Duration of wave spectrum at offshore boundary, in morphological time
      double precision                  :: dtbc                     = -123                 !  [s] (advanced) Timestep used to describe time series of wave energy and long wave flux at offshore boundary (not affected by morfac)
      double precision                  :: dthetaS_XB               = -123                 !  [deg] (advanced) The (counter-clockwise) angle in the degrees needed to rotate from the x-axis in SWAN to the x-axis pointing East
      integer                           :: nspectrumloc             = -123                 !  [-] (advanced) Number of input spectrum locations
      integer                           :: wbcversion               = -123                 !  [-] (advanced,silent) Version of wave boundary conditions
      integer                           :: nonhspectrum             = -123                 !  [-] (advanced) Spectrum format for wave action balance of nonhydrostatic waves

      ! [Section] Flow boundary condition parameters
      integer                           :: front                    = -123                 !  [name] Switch for seaward flow boundary
      character(slen)                   :: front_str                =  ' '                 !
      integer                           :: left                     = -123                 !  [name] Switch for lateral boundary at ny+1
      character(slen)                   :: left_str                 =  ' '                 !
      integer                           :: right                    = -123                 !  [name] Switch for lateral boundary at 0
      character(slen)                   :: right_str                =  ' '                 !
      integer                           :: back                     = -123                 !  [name] Switch for boundary at bay side
      character(slen)                   :: back_str                 =  ' '                 !
      integer                           :: ARC                      = -123                 !  [-] (advanced) Switch for active reflection compensation at seaward boundary
      double precision                  :: order                    = -123                 !  [-] (advanced) Switch for order of wave steering, 
      ! first order wave steering (short wave energy only), second oder wave steering (bound long wave corresponding to short wave forcing is added)
      integer                           :: highcomp                 = -123                 !  [-] (advanced) Switch for including the bound super harmonics in the boundary conditions, 
      integer                           :: freewave                 = -123                 !  [-] (advanced) Switch for free wave propagation 0 = use cg (default); 1 = use sqrt(gh) in instat = ts_2
      double precision                  :: epsi                     = -123                 !  [-] (advanced) Ratio of mean current to time varying current through offshore boundary
      integer                           :: nc                       = -123                 !  [-] (advanced,silent) Smoothing distance for estimating umean (defined as nr of cells)
      integer                           :: tidetype                 = -123                 !  [name] (advanced) Switch for offfshore boundary, velocity boundary or instant water level boundary
      character(slen)                   :: tidetype_str             =   ' '                !  [-] (advanced) Switch for offfshore boundary, velocity boundary or instant water level boundary (default)

      ! [Section] Tide boundary conditions
      character(slen)                   :: zs0file                  = 'abc'                !  [file] Name of tide boundary condition series
      integer                           :: tideloc                  = -123                 !  [-] Number of corner points on which a tide time series is specified
      integer                           :: paulrevere               = -123                 !  [name] Specifies tide on sea and land or two sea points if tideloc = 2
      character(slen)                   :: paulrevere_str           =  ' '                 !
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
      integer                           :: ndischarge               = -123                 !  [-] (advanced) Number of discharge locations
      integer                           :: ntdischarge              = -123                 !  [-] (advanced) Length of discharge time series
      character(slen)                   :: disch_loc_file           = 'abc'                !  [file] (advanced) Name of discharge locations file
      character(slen)                   :: disch_timeseries_file    = 'abc'                !  [file] (advanced) Name of discharge timeseries file

      ! [Section] Wave breaking parameters
      integer                           :: break                    = -123                 !  [name] Type of breaker formulation
      character(slen)                   :: break_str                =  ' '                 !  [-] Type of breaker formulation (1=roelvink, 2=baldock, 3=roelvink adapted, 4=roelvink on/off breaking)
      double precision                  :: gamma                    = -123                 !  [-] Breaker parameter in Baldock or Roelvink formulation
      double precision                  :: gamma2                   = -123                 !  [-] End of breaking parameter in Roelvink Daly formulation
      double precision                  :: alpha                    = -123                 !  [-] (advanced) Wave dissipation coefficient in Roelvink formulation
      double precision                  :: n                        = -123                 !  [-] (advanced) Power in Roelvink dissipation model
      double precision                  :: gammax                   = -123                 !  [-] (advanced) Maximum ratio wave height to water depth
      double precision                  :: delta                    = -123                 !  [-] (advanced) Fraction of wave height to add to water depth
      double precision                  :: wavfriccoef              = -123                 !  [-] Wave friction coefficient
      character(slen)                   :: wavfricfile              = 'abc'                !  [file] Wave friction file
      double precision                  :: fwcutoff                 = -123                 !  [-] Depth greater than which the bed friction factor is not applied
      double precision                  :: breakerdelay             = -123                 !  [-] (advanced) Switch to enable breaker delay model
      integer                           :: shoaldelay               = -123                 !  [-] (advanced,silent) Switch to enable shoaling delay
      double precision                  :: facsd                    = -123                 !  [-] (advanced,silent) fraction of the local wave length to use for shoaling delay depth
      double precision                  :: facrun                   = -123                 !  [-] (advanced,silent) calibration coefficient for short wave runup

      ! [Section] Roller parameters
      integer                           :: roller                   = -123                 !  [-] (advanced) Switch to enable roller model
      double precision                  :: beta                     = -123                 !  [-] (advanced) Breaker slope coefficient in roller model
      integer                           :: rfb                      = -123                 !  [-] (advanced) Switch to feed back maximum wave surface slope in roller energy balance, otherwise rfb = par%Beta

      ! [Section] Wave-current interaction parameters
      integer                           :: wci                      = -123                 !  [-] Turns on wave-current interaction
      double precision                  :: hwci                     = -123                 !  [m] (advanced) Minimum depth until which wave-current interaction is used
      double precision                  :: hwcimax                  = -123                 !  [m] (advanced) Maximum depth until which wave-current interaction is used
      double precision                  :: cats                     = -123                 !  [Trep] (advanced) Current averaging time scale for wci, in terms of mean wave periods

      ! [Section] Flow parameters
      integer                           :: bedfriction              = -123                 !  [name] Bed friction formulation
      character(slen)                   :: bedfriction_str          =  ' '                 !
      double precision                  :: bedfriccoef              = -123                 !  [-] Bed friction coefficient
      character(slen)                   :: bedfricfile              = 'abc'                !  [file] Bed friction file (valid for all bed friction coefficients)
      double precision                  :: maxcf                    = -123                 !  [-] max dimensionless friction coefficient (only for Manning and White-Colebrook)
      double precision                  :: nuh                      = -123                 !  [m^2s^-1] Horizontal background viscosity
      double precision                  :: nuhfac                   = -123                 !  [-] (advanced) Viscosity switch for roller induced turbulent horizontal viscosity
      double precision                  :: nuhv                     = -123                 !  [-] (advanced,silent) Longshore viscosity enhancement factor, following Svendsen (?)
      integer                           :: smag                     = -123                 !  [-] (advanced) Switch for smagorinsky subgrid model for viscocity
      integer*4                         :: friction_infiltration    = -123                 !  [-] turn on or off the effect of infiltration on bed roughness (Conley and Inman)
      integer*4                         :: friction_turbulence      = -123                 !  [-] turn on or off the effect of turbulence on bed roughness (Reniers Van Thiel)
      integer*4                         :: friction_acceleration    = -123                 !  [-] turn on or off the effect of acceleration on bed roughness (Morrison)
      character(slen)                   :: friction_acceleration_str=  ' '                 !
      double precision                  :: gamma_turb               = -123                 !  [-] calibration factor for turbulence contribution to bed roughness

      ! [Section] Coriolis force parameters
      double precision                  :: wearth                   = -123                 !  [hour^-1] (advanced) Angular velocity of earth calculated as: 1/rotation_time (in hours)
      double precision                  :: lat                      = -123                 !  [deg] (advanced) Latitude at model location  for computing Coriolis

      ! [Section] Wind parameters
      double precision                  :: rhoa                     = -123                 !  [kgm^-3] (advanced) Air density
      double precision                  :: Cd                       = -123                 !  [-] (advanced) Wind drag coefficient
      double precision                  :: windv                    = -123                 !  [ms^-1] Wind velocity, in case of stationary wind
      double precision                  :: windth                   = -123                 !  [deg] Nautical wind direction, in case of stationary wind
      character(slen)                   :: windfile                 = 'abc'                !  [file] Name of file with non-stationary wind data

      ! [Section] Groundwater parameters
      double precision                  :: kx                       = -123                 !  [ms^-1] (advanced) Darcy-flow permeability coefficient in x-direction
      double precision                  :: ky                       = -123                 !  [ms^-1] (advanced) Darcy-flow permeability coefficient in y-direction
      double precision                  :: kz                       = -123                 !  [ms^-1] (advanced) Darcy-flow permeability coefficient in z-direction
      double precision                  :: dwetlayer                = -123                 !  [m] (advanced) Thickness of the top soil layer interacting more freely with the surface water
      double precision                  :: aquiferbot               = -123                 !  [m] (advanced) Level of uniform aquifer bottom
      character(slen)                   :: aquiferbotfile           = 'abc'                !  [file] (advanced) Name of the aquifer bottom file
      double precision                  :: gw0                      = -123                 !  [m] (advanced) Level initial groundwater level
      character(slen)                   :: gw0file                  = 'abc'                !  [file] (advanced) Name of initial groundwater level file
      integer                           :: gwnonh                   = -123                 !  [-] (advanced) Switch to turn on or off non-hydrostatic pressure for groundwater
      integer                           :: gwfastsolve              = -123                 !  [-] (advanced,silent) Reduce full 2D non-hydrostatic solution to quasi-explicit in longshore direction
      integer                           :: gwscheme                 = -123                 !  [name] (advanced) Scheme for momentum equation
      character(slen)                   :: gwscheme_str             =  ' '                 !
      double precision                  :: gwReturb                 = -123                 !  [-] (advanced) Reynolds number for start of turbulent flow in case of gwscheme = turbulent
      integer                           :: gwheadmodel              = -123                 !  [name] (advanced) Model to use for vertical groundwater head
      character(slen)                   :: gwheadmodel_str          =  ' '                 !
      integer                           :: gwhorinfil               = -123                 !  [-] (advanced) switch to include horizontal infiltration from surface water to groundwater

      ! [Section] Q3D sediment transport parameters
      double precision                  :: vonkar                   = -123                 !  [-] (advanced,silent) von Karman constant
      double precision                  :: vicmol                   = -123                 !  [-] (advanced,silent) molecular viscosity
      integer                           :: kmax                     = -123                 !  [-] (advanced,silent) Number of sigma layers in Quasi-3D model; kmax = 1 is without vertical structure of flow and suspensions
      double precision                  :: sigfac                   = -123                 !  [-] (advanced,silent) dsig scales with log(sigfac)
      double precision                  :: deltar                   = -123                 !  [-] (advanced,silent) estimated ripple height
      double precision                  :: rwave                   = -123                  !  [-] (advanced,silent) user-defined wave roughness adjustment factor

      ! [Section] Non-hydrostatic correction parameters
      integer                           :: solver                   = -123                 !  [name] (advanced) Solver used to solve the linear system
      character(slen)                   :: solver_str               =  ' '                 ! 
      integer                           :: solver_maxit             = -123                 !  [-] (advanced) Maximum number of iterations in the linear sip solver
      double precision                  :: solver_acc               = -123                 !  [-] (advanced) accuracy with respect to the right-hand side used
                                                                                           !                 in the following termination criterion:
                                                                                           !                     ||b-Ax || < acc*||b||
      double precision                  :: solver_urelax            = -123                 !  [-] (advanced) Underrelaxation parameter
      double precision                  :: kdmin                    = -123                 !  [-] (advanced) Minimum value of kd (pi/dx > min(kd))
      double precision                  :: dispc                    = -123                 !  [?] (advanced) Coefficient in front of the vertical pressure gradient
      double precision                  :: Topt                     = -123                 !  [s] (advanced) Absolute period to optimize coefficient
      integer                           :: nhbreaker                = -123                 !  [-] (advanced) Non-hydrostatic breaker model
      double precision                  :: breakviscfac             = -123                 !  [-] (advanced) Factor to increase viscosity during breaking
      double precision                  :: maxbrsteep               = -123                 !  [-] (advanced) Maximum wave steepness criterium
      double precision                  :: secbrsteep               = -123                 !  [-] (advanced) Secondary maximum wave steepness criterium
      double precision                  :: reformsteep              = -123                 !  [-] (advanced) Wave steepness criterium to reform after breaking
      double precision                  :: breakvisclen             = -123                 !  [-] (advanced) Ratio between local depth and length scale in extra breaking viscosity
      integer                           :: nonhq3d                  = -123                 !  [-] (advanced,silent) Turn on or off the the reduced two-layer nonhydrostatic model, default = 0
      double precision                  :: nhlay                    = -123                 !  [-] (advanced) Layer distribution in the nonhydrostatic model, default = 0.33

      ! [Section] Bed composition parameters
      double precision                  :: rhos                     = -123                 !  [kgm^-3] Solid sediment density (no pores)
      integer                           :: ngd                      = -123                 !  [-] Number of sediment classes
      integer                           :: nd                       = -123                 !  [-] (advanced) Number of computational layers in the bed
      double precision                  :: dzg1                     = -123                 !  [m] (advanced) Thickness of top sediment class layers
      double precision                  :: dzg2                     = -123                 !  [m] (advanced) Nominal thickness of variable sediment class layer
      double precision                  :: dzg3                     = -123                 !  [m] (advanced) Thickness of bottom sediment class layers
      double precision                  :: por                      = -123                 !  [-] Porosity
      double precision                  :: D15(99)                  = -123                 !  [m] D15 grain size per grain type
      double precision                  :: D50(99)                  = -123                 !  [m] D50 grain size per grain type
      double precision                  :: D90(99)                  = -123                 !  [m] D90 grain size per grain type
      double precision                  :: ws                       = 0.02d0               !  [m/s] average fall velocity (is computed in morphevolution)
      double precision                  :: sedcal(99)               = -123                 !  [-] (advanced) Sediment transport calibration coefficient per grain type
      double precision                  :: ucrcal(99)               = -123                 !  [-] (advanced) Critical velocity calibration coefficient per grain type

      ! [Section] Sediment transport parameters
      integer                           :: waveform                 = -123                 !  [name] Wave shape model
      character(slen)                   :: waveform_str             =  ' '                 !
      integer                           :: form                     = -123                 !  [name] Equilibrium sediment concentration formulation
      character(slen)                   :: form_str                 =  ' '                 !
      integer                           :: sws                      = -123                 !  [-] (advanced) Switch to enable short wave and roller stirring and undertow
      integer                           :: lws                      = -123                 !  [-] (advanced) Switch to enable long wave stirring
      double precision                  :: BRfac                    = -123                 !  [-] (advanced) Calibration factor surface slope
      double precision                  :: facsl                    = -123                 !  [-] (advanced) Factor bedslope effect
      double precision                  :: z0                       = -123                 !  [m] (advanced) Zero flow velocity level in Soulsby and van Rijn (1997) sediment concentration
      double precision                  :: smax                     = -123                 !  [-] (advanced) Maximum Shields parameter for equilibrium sediment concentration acc. Diane Foster
      double precision                  :: tsfac                    = -123                 !  [-] (advanced) Coefficient determining Ts = tsfac * h/ws in sediment source term
      double precision                  :: facua                    = -123                 !  [-] (advanced) Calibration factor time averaged flows due to wave skewness and asymmetry
      double precision                  :: facSk                    = -123                 !  [-] (advanced) Calibration factor time averaged flows due to wave skewness
      double precision                  :: facAs                    = -123                 !  [-] (advanced) Calibration factor time averaged flows due to wave asymmetry
      integer                           :: turbadv                  = -123                 !  [name] (advanced) Switch to activate turbulence advection model for short and or long wave turbulence
      character(slen)                   :: turbadv_str              =  ' '                 !
      integer                           :: turb                     = -123                 !  [name] (advanced) Switch to include short wave turbulence
      character(slen)                   :: turb_str                 =  ' '                 !
      double precision                  :: Tbfac                    = -123                 !  [-] (advanced) Calibration factor for bore interval Tbore: Tbore = Tbfac*Tbore
      double precision                  :: Tsmin                    = -123                 !  [s] (advanced) Minimum adaptation time scale in advection diffusion equation sediment
      integer                           :: lwt                      = -123                 !  [-] (advanced) Switch to enable long wave turbulence

      double precision                  :: betad                    = -123                 !  [-] (advanced) Dissipation parameter long wave breaking turbulence
      character(slen)                   :: swtable                  = 'abc'                !  [-] (deprecated)Name of intra short wave assymetry and skewness table
      integer                           :: sus                      = -123                 !  [-] (advanced) Calibration factor for suspensions transports
      integer                           :: bed                      = -123                 !  [-] (advanced) Calibration factor for bed transports
      integer                           :: bulk                     = -123                 !  [-] (advanced) Switch to compute bulk transport rather than bed and suspended load separately
      double precision                  :: facDc                    = -123                 !  [-] (advanced) Option to control sediment diffusion coefficient
      double precision                  :: jetfac                   = -123                 !  [-] (advanced,silent) Option to mimic turbulence production near revetments
      integer                           :: fallvelred               = -123                 !  [-] Switch to reduce fall velocity for high concentrations
      integer                           :: dilatancy                = -123                 !  [-] Switch to reduce critical shields number due dilatancy
      double precision                  :: rheeA                    = -123                 !  [-] A parameter in the Van Rhee expression
      double precision                  :: pormax                   = -123                 !  [-] Max porosity used in the experession of Van Rhee
      double precision                  :: reposeangle              = -123                 !  [deg] Angle of internal friction
      integer                           :: bdslpeffmag              = -123                 !  [name] Modify the magnitude of the sediment transport based on the bed slope, uses facsl
      integer                           :: bdslpeffini              = -123                 !  [name] Modify the critical shields parameter based on the bed slope
      integer                           :: bdslpeffdir              = -123                 !  [name] Modify the direction of the sediment transport based on the bed slope
      double precision                  :: bdslpeffdirfac           = -123                 !  [-] Calibration factor in the modification of the direction
      double precision                  :: ci                       = -123                 !  [-] (advanced) Mass coefficient in Shields inertia term
      double precision                  :: phit                     = -123                 !  [-] (advanced) Phase lag angle in Nielsen transport equation 
      integer*4                         :: incldzdx                 = -123                 !  [-] (advanced,silent) Turn on or off dzsdx term in Shields
      integer*4                         :: inclrelweight            = -123                 !  [-] (advanced,silent) Turn on or off infilitration/exfiltration effect on particle weight
      integer*4                         :: streaming                = -123                 !  [-] (advanced,silent) Turn on or off streaming contribution in Nielsen 2006 
      real*8                            :: uprushfac                = -123                 !  [-] (advanced,silent) Factor to increase uprush transport
      real*8                            :: backwashfac              = -123                 !  [-] (advanced,silent) Factor to increase backwash transport
      real*8                            :: yturb                    = -123                 !  [-] (advanced,silent) factor for distribution of near-bed turbulence into bed load and suspended load transport
      real*8                            :: facthr                   = -123                 !  [-] (advanced,silent) multiplication factor for numerically estimated long wave roller thickness
      integer                           :: sedfricfac               = -123                 !  [name] (advanced,silent) Wave shape model
      character(slen)                   :: sedfricfac_str           =  ' '                 !
      real*8                            :: Arms                     = -123                 !  [m] (advanced,silent) swash excursion for Nielsen expression
      real*8                            :: Ctrans                   = -123                 !  [-] (advanced,silent) Constant in Nielsen expression (default=12, could be 20)
      integer*4                         :: slopecorr                = -123                 !  [name] (advanced,silent) which slope correction used in nielsen formula: 'nielsen' (2002) or hughes&masselink (xxxx) ('hughes_masselink')
      character(slen)                   :: slopecorr_str            =  ' '                 !
      real*8                            :: fsed                     = -123                 !  [-] (advanced,silent) constant sediment friction factor
      integer*4                         :: phaselag                 = -123                 !  [-] (advanced,silent) 1 = phase lag, 0 = no phase lag
      real*8                            :: thetcr                   = -123                 !  [-] (advanced,silent) critical shields param
      integer*4                         :: bermslopetransport       = -123                 !  [-] (advanced,silent) Turn on or off bermslope swash transport model
      integer*4                         :: bermslopebed             = -123                 !  [-] (advanced,silent) Turn on or off bermslope swash transport model for bed load
      integer*4                         :: bermslopesus             = -123                 !  [-] (advanced,silent) Turn on or off bermslope swash transport model for suspended load
      double precision                  :: bermslope                = -123                 !  [-] (advanced,silent) Swash zone slope for (semi-) reflective beaches
      double precision                  :: bermslopefac             = -123                 !  [-] (advanced,silent) Bed slope transport factor for bermslope model
      double precision                  :: bermslopegamma           = -123                 !  [-] (advanced,silent) Wave height - water depth ratio to turn on bermslope model in surf-beat
      double precision                  :: bermslopedepth           = -123                 !  [-] (advanced,silent) Water depth to turn on on bermslope model in stationary and nonh


      ! [Section] Morphology parameters
      double precision                  :: morfac                   = -123                 !  [-] Morphological acceleration factor
      integer                           :: morfacopt                = -123                 !  [-] (advanced) Switch to adjusting output times for morfac
      double precision                  :: morstart                 = -123                 !  [s] Start time morphology, in morphological time
      double precision                  :: morstop                  = -123                 !  [s] Stop time morphology, in morphological time
      double precision                  :: wetslp                   = -123                 !  [-] Critical avalanching slope under water (dz/dx and dz/dy)
      double precision                  :: dryslp                   = -123                 !  [-] Critical avalanching slope above water (dz/dx and dz/dy)
      double precision                  :: lsgrad                   = -123                 !  [1/m] Factor to include longshore transport gradient in 1D simulations
                                                                                           !        dSy/dy=lsgrad*Sy; dimension 1/length scale of longshore gradients
      double precision                  :: hswitch                  = -123                 !  [m] (advanced) Water depth at which is switched from wetslp to dryslp
      double precision                  :: dzmax                    = -123                 !  [m/s/m] (advanced) Maximum bed level change due to avalanching
      integer                           :: struct                   = -123                 !  [-] Switch for enabling hard structures
      character(slen)                   :: ne_layer                 = 'abc'                !  [file] Name of file containing thickness of the erodible layer

      ! [Section] Output variables
      integer                           :: timings                  = -123                 !  [-] (advanced) Switch enable progress output to screen
      double precision                  :: tstart                   = -123                 !  [s] Start time of output, in morphological time
      double precision                  :: tint                     = -123                 !  [s] (deprecated) Interval time of global output (replaced by tintg)
      double precision                  :: tintg                    = -123                 !  [s] Interval time of global output
      double precision                  :: tintp                    = -123                 !  [s] Interval time of point and runup gauge output
      double precision                  :: tintc                    = -123                 !  [s] (advanced) Interval time of cross section output
      double precision                  :: tintm                    = -123                 !  [s] Interval time of mean, var, max, min output
      character(slen)                   :: tsglobal                 = 'abc'                !  [-] (advanced) Name of file containing timings of global output
      character(slen)                   :: tspoints                 = 'abc'                !  [-] (advanced) Name of file containing timings of point output
      character(slen)                   :: tsmean                   = 'abc'                !  [-] (advanced) Name of file containing timings of mean, max, min and var output
      integer                           :: nglobalvar               = -123                 !  [-] Number of global output variables (as specified by user)
      character(maxnamelen)             :: globalvars(numvars)      = 'abc'                !  [-] (advanced) Mnems of global output variables, 
      ! not per se the same size as nglobalvar (invalid variables, defaults)
      integer                           :: nmeanvar                 = -123                 !  [-] Number of mean, min, max, var output variables
      character(maxnamelen)             :: meanvars(numvars)        = 'abc'                !  [-] (advanced) Mnems of mean output variables (by variables)
      integer                           :: npointvar                = -123                 !  [-] Number of point output variables
      character(maxnamelen)             :: pointvars(numvars)       = 'abc'                !  [-] (advanced) Mnems of point output variables (by variables)
      integer                           :: npoints                  = -123                 !  [-] Number of output point locations
      integer                           :: nrugauge                 = -123                 !  [-] Number of output runup gauge locations
      integer, pointer              :: pointtypes(:) => NULL()                                   !  [-] (advanced) Point types (0 = point, 1 = rugauge)
      double precision, pointer     :: xpointsw(:) => NULL()                                     !  (advanced) world x-coordinate of output points
      double precision, pointer     :: ypointsw(:) => NULL()                                     !  (advanced) world y-coordinate of output points
           
      integer                           :: nrugdepth                = -123                 !  [-] (advanced) Number of depths to compute runup in runup gauge
      double precision                  :: rugdepth(9999)             = -123               !  [m] (advanced) Minimum depth for determination of last wet point in runup gauge
      integer                           :: outputformat             = OUTPUTFORMAT_DEBUG   !  [name] (advanced) Output file format
      character(slen)                   :: outputformat_str         = 'debug'              !
      character(slen)                   :: ncfilename               = 'xboutput.nc'        !  [file] (advanced) xbeach netcdf output file name
      integer                           :: outputprecision          = -123                 !  [name] switch between single and double precision output in NetCDF
      character(slen)                   :: outputprecision_str      =  ' '                 !
      character(64)                     :: stationid(9999)            = 'abc'              !  [-] (advanced,silent) Station id names of output points

      ! Projection units (not to be used, only pass to output, this limit is too short for WKT....)
      ! This could be the proj4 string +init=epsg:28992
      ! [Section] Output projection
      character(slen)                   :: projection               = ' '                  !  [-] (advanced) projection string
      integer                           :: rotate                   = -123                 !  [-] Rotate output as postprocessing with given angle
      integer                           :: remdryoutput             = -123                 !  [-] Remove dry output points from output data of zs etc.

      ! [Section] Drifters parameters
      integer                           :: ndrifter                 = -123                 !  [-] Number of drifers
      character(slen)                   :: drifterfile              = 'abc'                !  [file] Name of drifter data file

      ! [Section] Ship parameters
      character(slen)                   :: shipfile                 = 'abc'                !  [file] Name of ship data file
      integer                           :: nship                    = -123                 !  [-] (advanced) Number of ships

      ! [Section] Vegetation parameters
      character(slen)                   :: veggiefile               = 'abc'                !  [-] Name of veggie species list file
      character(slen)                   :: veggiemapfile            = 'abc'                !  [-] Name of veggie species map file
      integer                           :: nveg                     = -123                 !  [-] Number of vegetation species
      integer                           :: vegnonlin                = -123                 !  [-] include non-linear wave effect [1] or not [0]
      integer                           :: vegcanflo                = -123                 !  [-] include incanopy flow [1] or not [0]
      integer                           :: veguntow                 = -123                 !  [-] include undertow in phase-averaged vegetati
      integer                           :: porcanflow               = -123                 !  [-] Compute in-canopy flow
      double precision                  :: Kp                       = -123                 !  [-] Laminar resistance factor (in-canopy flow)
      double precision                  :: Cm                       = -123                 !  [-] Inertia coefficient (in-canopy flow)

      ! [Section] Wave numerics parameters
      integer                           :: scheme                   = -123                 !  [name] (advanced) Numerical scheme for wave propagation
      character(slen)                   :: scheme_str               =  ' '                 !  [-] (advanced) Use first-order upwind (upwind_1), second order upwind (upwind_2) or Lax-Wendroff (lax_wendroff)
      double precision                  :: wavint                   = -123                 !  [s] Interval between wave module calls (only in stationary wave mode)
      double precision                  :: maxerror                 = -123                 !  [m] (advanced) Maximum wave height error in wave stationary iteration
      integer                           :: maxiter                  = -123                 !  [-] (advanced) Maximum number of iterations in wave stationary
      real*8                            :: swkhmin                  = -123                 !  [-] (advanced,silent) Minimum kh value to include in wave action balance, lower included in NLSWE (default -1.d0)

      ! [Section] Flow numerics parameters
      double precision                  :: eps                      = -123                 !  [m] Threshold water depth above which cells are considered wet
      double precision                  :: eps_sd                   = -123                 !  [m/s] Threshold velocity difference to determine conservation of energy head versus momentum
      double precision                  :: umin                     = -123                 !  [m/s] Threshold velocity for upwind velocity detection and for vmag2 in equilibrium sediment concentration
      double precision                  :: hmin                     = -123                 !  [m] Threshold water depth above which Stokes drift is included
      integer                           :: secorder                 = -123                 !  [-] (advanced) Use second order corrections to advection/non-linear terms based on MacCormack scheme
      integer                           :: oldhu                    = -123                 !  [-] (advanced,silent) Switch to enable old hu calculation

      ! [Section] Sediment transport numerics parameters
      double precision                  :: thetanum                 = -123                 !  [-] (advanced) Coefficient determining whether upwind (1) or central scheme (0.5) is used.
      integer                           :: sourcesink               = -123                 !  [-] (advanced) Switch to enable source-sink terms to calculate bed level change rather than suspended transport gradients
      double precision                  :: cmax                     = -123                 !  [-] (advanced) Maximum allowed sediment concentration

      ! [Section] Bed update numerics parameters
      double precision                  :: frac_dz                  = -123                 !  [-] (advanced) Relative thickness to split time step for bed updating
      integer                           :: nd_var                   = -123                 !  [-] (advanced) Index of layer with variable thickness
      double precision                  :: split                    = -123                 !  [-] (advanced) Split threshold for variable sediment layer (ratio to nominal thickness)
      double precision                  :: merge                    = -123                 !  [-] (advanced) Merge threshold for variable sediment layer (ratio to nominal thickness)
      integer                           :: nsetbathy                = -123                 !  [-] (advanced) Number of prescribed bed updates
      character(slen)                   :: setbathyfile             = 'abc'                !  [file] (advanced) Name of prescribed bed update file

      ! [Section] MPI parameters
      integer                           :: mpiboundary              = -123                 !  [name] (advanced) Fix mpi boundaries along y-lines, x-lines, use manual defined domains or find shortest boundary automatically
      character(slen)                   :: mpiboundary_str          =  ' '                 !
      integer                           :: mmpi                     = -123                 !  [-] (advanced) Number of domains in cross-shore direction when manually specifying mpi domains
      integer                           :: nmpi                     = -123                 !  [-] (advanced) Number of domains in alongshore direction when manually specifying mpi domains
      
      ! [Section] Constants, not read in params.txt
      double precision                  :: px                       = 4.d0*atan(1.d0)      !  [-] Pi
      complex(kind(0.0d0))              :: compi                    = -123                 !  [-] Imaginary unit
      double precision                  :: rhog8                    = -123                 !  [Nm^-3] 1/8*rho*g
      double precision                  :: irhog8                   = -123                 !  [N^-1m^3] (1/8*rho*g)^-1

      ! [Section] Variables, not read in params.txt
      double precision                  :: dt                       = -123                 !  [s] Computational time step, in hydrodynamic time
      double precision                  :: t                        = -123                 !  [s] Computational time, in hydrodynamic time
      double precision                  :: tnext                    = -123                 !  [s] Next time point for output or wave stationary calculation, in hydrodynamic time

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
      !
      call writelog('sl','','Reading input parameters: ')
      !
      ! Check params.txt exists
      call check_file_exist('params.txt')
      !
      !
      ! Collections of default sets
      par%useXBeachGSettings    = readkey_int ('params.txt','useXBeachGSettings',     0,0,1,strict=.true.,silent=.true.)
      ! par%useXBeachGSettings = 1
      !      
      !
      ! Backward compatibility
      call check_instat_backward_compatibility(par)
      !
      !
      ! Physical processes
      call writelog('l','','--------------------------------')
      call writelog('l','','Physical processes: ')
      !
      ! Wavemodel
      if (isSetParameter('params.txt','wavemodel')) then
          ! if defined in the params.txt
          call setallowednames('stationary', WAVEMODEL_STATIONARY, 'surfbeat', WAVEMODEL_SURFBEAT, 'nonh', WAVEMODEL_NONH)
          call setoldnames('0','1','2')
          if(par%useXBeachGSettings==0) then
             call parmapply('wavemodel',2, par%wavemodel, par%wavemodel_str)
          else
             call parmapply('wavemodel',3, par%wavemodel, par%wavemodel_str)
          endif
      elseif (par%wavemodel /= -123) then 
         ! if already determined with backward compatibility (instat)
         if (par%wavemodel == 0) then
             call writelog('l','','               wavemodel =stationary')
         elseif (par%wavemodel == 1) then  
             call writelog('l','','               wavemodel =surfbeat')
         else
             call writelog('l','','               wavemodel =nonh')
        endif
      elseif (par%wavemodel == -123) then   
         ! no wavemodel determined -> error
         call writelog('lse','(a)','Error: XBeach cannot run without a hydrodynamic type')
         call writelog('lse','(a)','       Set wavemodel to stationary, surfbeat or nonh')
         call halt_program
      endif
      !
      ! Rest of the processes
      par%cyclic      = readkey_int ('params.txt','cyclic',        0,        0,     1,strict=.true.)
      if (par%wavemodel==WAVEMODEL_NONH) then
         par%swave       = readkey_int ('params.txt','swave',         0, 0, 1,strict=.true.)
      else
         par%swave       = readkey_int ('params.txt','swave',      1, 0, 1,strict=.true.)
      endif
      par%single_dir  = readkey_int ('params.txt','single_dir',    0,        0,     1,strict=.true.)
      par%lwave       = readkey_int ('params.txt','lwave',         1,        0,     1,strict=.true.)
      par%flow        = readkey_int ('params.txt','flow',          1,        0,     1,strict=.true.)
      par%sedtrans    = readkey_int ('params.txt','sedtrans',      1,        0,     1,strict=.true.)
      par%morphology  = readkey_int ('params.txt','morphology',    par%sedtrans,  0,1,strict=.true.)
      par%avalanching = readkey_int ('params.txt','avalanching',   par%morphology,0,1,strict=.true.)
      par%gwflow      = readkey_int ('params.txt','gwflow',        par%useXBeachGSettings,        0,     1,strict=.true.)
      par%q3d         = readkey_int ('params.txt','q3d',           0,        0,     1,silent=.true.,strict=.true.)
      par%swrunup     = readkey_int ('params.txt','swrunup',       0,        0,     1,silent=.true.,strict=.true.)
      par%ships       = readkey_int ('params.txt','ships',         0,        0,     1,strict=.true.)
      ! nship defined by shipfile
      if (par%ships == 0) par%nship = 0
      par%bchwiz      = readkey_int ('params.txt','bchwiz',        0,        0,     1,silent=.true.,strict=.true.)
      par%vegetation  = readkey_int ('params.txt','vegetation',    0,        0,     1,strict=.true.)
      par%setbathy    = readkey_int ('params.txt','setbathy',      0,        0,     1,strict=.true.)
      par%viscosity   = readkey_int ('params.txt','viscosity',     1,        0,     1,strict=.true.)
      par%advection   = readkey_int ('params.txt','advection',     1,        0,     1,strict=.true.)
      par%wind        = readkey_int ('params.txt','wind',          1,        0,     1,strict=.true.)
      !
      ! Grid parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Grid parameters: ')
      !
      ! check gridform
      call setallowednames('xbeach',GRIDFORM_XBEACH,'delft3d',GRIDFORM_DELFT3D)
      call setoldnames('0','1')
      call parmapply('gridform',1,par%gridform,par%gridform_str)
      !
      ! gridform switch
      if (par%gridform==GRIDFORM_XBEACH) then
         ! XBeach grid parameters
         par%xori  = readkey_dbl('params.txt','xori',  0.d0,   -1d9,      1d9)
         par%yori  = readkey_dbl('params.txt','yori',  0.d0,   -1d9,      1d9)
         par%alfa  = readkey_dbl('params.txt','alfa',  0.d0,   0.d0,   360.d0)
         par%nx    = readkey_int('params.txt','nx',     50,      2,     10000,required=.true.)
         par%ny    = readkey_int('params.txt','ny',      2,      0,     10000,required=.true.)
         if(par%useXBeachGSettings==1) then
            if(par%ny>0) then
               call writelog('lws','(a)','XBeach-G settings cannot be applied in 2DH models')
               call writelog('lws','(a)','Set model to ny=0')
               call halt_program
            endif
         endif  
         par%posdwn= readkey_dbl('params.txt','posdwn', 1.d0,     -1.d0,     1.d0)
         if (par%setbathy .ne. 1) then
            par%depfile = readkey_name('params.txt','depfile',required=.true.)
            call check_file_exist(par%depfile)
            call check_file_length(par%depfile,par%nx+1,par%ny+1)
         else
            par%depfile = readkey_name('params.txt','depfile')
         endif
         par%vardx = readkey_int('params.txt','vardx',   0,      0,         1,strict=.true.)

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
        ! Bermslope only in 1D
        if(par%bermslope>0) then
        if(par%ny>0) then
            call writelog('lws','(a)','bermslope cannot be applied in 2DH models')
            call writelog('lws','(a)','Set model to ny=0')
            call halt_program
        endif
        endif  
      ! Q3d grid
         par%nz = readkey_int ('params.txt','nz',    1,        1,     100)
      ! Wave directional grid
      if(par%swave==1) then
         par%thetamin = readkey_dbl ('params.txt','thetamin', -90.d0,    -360.d0,  360.d0,required=.true.)
         par%thetamax = readkey_dbl ('params.txt','thetamax',  90.d0,    -360.d0,  360.d0,required=.true.)
         par%thetanaut= readkey_int ('params.txt','thetanaut',    0,        0,     1,strict=.true.)
         if (par%single_dir==1) then
            call writelog('ls','','dtheta will automatically be computed from thetamin and thetamax for single_dir = 1')
            par%dtheta_s   = readkey_dbl ('params.txt','dtheta_s',    10.d0,      0.1d0,   20.d0,required=.true.)
         else
            par%dtheta   = readkey_dbl ('params.txt','dtheta',    10.d0,      0.1d0,   180.d0,required=.true.)
         endif
      endif
      !
      ! Model time parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Model time parameters: ')
      if(par%useXBeachGSettings==0) then
         par%CFL     = readkey_dbl ('params.txt','CFL',     0.7d0,     0.1d0,      0.9d0)
      else
         par%CFL     = readkey_dbl ('params.txt','CFL',     0.6d0,     0.1d0,      0.9d0)
      endif
      par%dtset   = readkey_dbl ('params.txt','dtset',   0.0d0,     0.001d0,   100.d0)
      par%tstop   = readkey_dbl ('params.txt','tstop', 2000.d0,      1.d0, 1000000.d0,required=.true.)
      par%defuse  = readkey_int ('params.txt','defuse',        1,        0,     1,strict=.true.,silent=.true.)
      if (par%wavemodel==WAVEMODEL_STATIONARY .or. par%wavemodel==WAVEMODEL_SURFBEAT) then
         par%maxdtfac  = readkey_dbl ('params.txt','maxdtfac', 50.d0,      10.d0, 200.d0)
      else
         par%maxdtfac  = readkey_dbl ('params.txt','maxdtfac', 500.d0,      100.d0, 1000.d0)
      endif
      !
      ! Physical constants
      call writelog('l','','--------------------------------')
      call writelog('l','','Physical constants: ')
      par%rho        = readkey_dbl ('params.txt','rho',       1025.0d0,  1000.0d0,  1040.0d0)
      par%g          = readkey_dbl ('params.txt','g',         9.81d0,    9.7d0,     9.9d0)
      par%depthscale = readkey_dbl ('params.txt','depthscale',1.0d0,     1.0d0,     200.d0)
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
      par%hotstartflow = readkey_int ('params.txt','hotstartflow',    0,        0,     1,strict=.true.,silent=.true.)
      !
      ! Wave boundary condition parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Wave boundary condition parameters: ')
      !      
      ! New method: wbctype
      if (isSetParameter('params.txt','wbctype')) then
        call setallowednames('params',       WBCTYPE_PARAMS,       &
        'parametric', WBCTYPE_PARAMETRIC, &
        'swan',       WBCTYPE_SWAN,       &
        'vardens',    WBCTYPE_VARDENS,    &
        'off',        WBCTYPE_OFF,        &
        'jonstable',  WBCTYPE_JONS_TABLE, &
        'reuse',      WBCTYPE_REUSE,      &
        'ts_1',       WBCTYPE_TS_1,       &
        'ts_2',       WBCTYPE_TS_2,       &
        'ts_nonh',    WBCTYPE_TS_NONH)
        call setoldnames('0','1','2','3','4','5','6','7', '8', '9')
        call parmapply('wbctype', 1, par%wbctype, par%wbctype_str)
      elseif (par%wbctype /= -123) then 
          ! if already determined with backward compatibility (instat)
         if (par%wbctype == WBCTYPE_PARAMS) then
             call writelog('l','','                 wbctype =params')
         elseif (par%wbctype == WBCTYPE_PARAMETRIC) then
             call writelog('l','','                 wbctype =parametric')
         elseif (par%wbctype == WBCTYPE_SWAN) then  
             call writelog('l','','                 wbctype =swan')
         elseif (par%wbctype == WBCTYPE_VARDENS) then  
             call writelog('l','','                 wbctype =vardens')
         elseif (par%wbctype == WBCTYPE_OFF) then  
             call writelog('l','','                 wbctype =off')
         elseif (par%wbctype == WBCTYPE_JONS_TABLE) then  
             call writelog('l','','                 wbctype =jonstable')
         elseif (par%wbctype == WBCTYPE_REUSE) then  
             call writelog('l','','                 wbctype =reuse')
         elseif (par%wbctype == WBCTYPE_TS_1) then  
             call writelog('l','','                 wbctype =ts_1')
         elseif (par%wbctype == WBCTYPE_TS_2) then  
             call writelog('l','','                 wbctype =ts_2')
         elseif (par%wbctype == WBCTYPE_TS_NONH) then  
             call writelog('l','','                 wbctype =ts_nonh')
        endif
      else  
         ! no wbctype determined -> error
         call writelog('lse','(a)','Error: XBeach cannot run without a type of boundary condition')
         call writelog('lse','(a)','       Set wbctype')
         call halt_program
      endif
      !
      if ( par%wbctype==WBCTYPE_PARAMETRIC .or. &
      par%wbctype==WBCTYPE_SWAN .or. &
      par%wbctype==WBCTYPE_VARDENS .or. &
      par%wbctype==WBCTYPE_JONS_TABLE &
      )then
         par%bcfile = readkey_name('params.txt','bcfile')
         call check_file_exist(par%bcfile)
         call checkbcfilelength(par%tstop,par%wbctype,par%bcfile,filetype)
         ! Only carried out on xmaster so:
#ifdef USEMPI
         call xmpi_bcast(filetype,toall)
#endif
      elseif (par%wbctype==WBCTYPE_REUSE .or. par%instat == INSTAT_REUSE) then
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
            ! Check for file length is done later, after tstop is adjusted for morfac
         elseif (.not. fe3 .and. (fe1 .and. fe2)) then
            par%nonhspectrum = 0
            ! Check for file length is done later, after tstop is adjusted for morfac
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
      if (par%wbctype == WBCTYPE_PARAMS.or. par%single_dir==1) then
         par%nonhspectrum    = readkey_int ('params.txt','nonhspectrum', 0,          0,       1 ,strict=.true.)
         par%Hrms  = readkey_dbl ('params.txt','Hrms',      1.d0,      0.d0,    10.d0)
         par%Tm01  = readkey_dbl ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
         par%Trep  = readkey_dbl ('params.txt','Trep',     par%Tm01,   1.d0,    20.d0)
         par%dir0  = readkey_dbl ('params.txt','dir0',    270.d0,    -360.d0,   360.d0)
         par%m     = readkey_int ('params.txt','m',        10,         2,      128)
      elseif (par%wbctype == WBCTYPE_TS_1 .or. par%wbctype == WBCTYPE_TS_2) then
         call check_file_exist('bc/gen.ezs')
         par%Tm01  = readkey_dbl ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
         par%Trep  = readkey_dbl ('params.txt','Trep',     par%Tm01,   1.d0,    20.d0)
         par%dir0  = readkey_dbl ('params.txt','dir0',    270.d0,    -360.d0,   360.d0)
         par%m     = readkey_int ('params.txt','m',        10,         2,      128)
      elseif (par%wbctype == WBCTYPE_TS_NONH) then
         call check_file_exist('boun_U.bcf')
      endif
      !
      ! If Tlong is defined: than bichromatic waves        
      if (isSetParameter('params.txt','Tlong')) then 
         par%Tlong = readkey_dbl ('params.txt','Tlong',    80.d0,     20.d0,   300.d0)
      endif
      !
      call setallowednames('neumann',LATERALWAVE_NEUMANN,'wavecrest',LATERALWAVE_WAVECREST,'cyclic',LATERALWAVE_CYCLIC)
      call setoldnames('0','1')
      call parmapply('lateralwave',1,par%lateralwave,par%lateralwave_str)
      par%bclwonly   = readkey_int ('params.txt','bclwonly',  0,0,1,strict=.true.,silent=.true.)
      par%Sfold   = readkey_int ('params.txt','Sfold',  0,0,1,strict=.true.,silent=.true.)
      !
      !
      ! Wave-spectrum boundary condition parameters
      if ( par%wbctype==WBCTYPE_PARAMETRIC .or. &
      par%wbctype==WBCTYPE_SWAN .or. &
      par%wbctype==WBCTYPE_VARDENS .or. &
      par%wbctype==WBCTYPE_JONS_TABLE &
      )then
      !
         call writelog('l','','--------------------------------')
         call writelog('l','','Wave-spectrum boundary condition parameters: ')
         if (par%wavemodel==WAVEMODEL_NONH) then
             par%nonhspectrum= 1
         else
            par%nonhspectrum = 0
         endif
         par%random          = readkey_int ('params.txt','random',       1,          0,          1     ,strict=.true.)
         par%fcutoff         = readkey_dbl ('params.txt','fcutoff',      0.d0,       0.d0,       40.d0   )
         par%nspr            = readkey_int ('params.txt','nspr',         0,          0,          1     ,silent=.true.,strict=.true.)
         par%trepfac         = readkey_dbl ('params.txt','trepfac',      0.01d0,     0.d0,       1.d0    )
         if (par%nonhspectrum==1) then
            par%sprdthr         = readkey_dbl ('params.txt','sprdthr',      0.00d0,     0.d0,       1.d0    )
         else
            par%sprdthr         = readkey_dbl ('params.txt','sprdthr',      0.08d0,     0.d0,       1.d0    )
         endif
         par%correctHm0      = readkey_int ('params.txt','correctHm0',   1,          0,          1  ,silent=.true.,strict=.true.)
         par%Tm01switch      = readkey_int ('params.txt','Tm01switch',   0,          0,          1       ,strict=.true.)
         !
         if (filetype==0) then
            par%rt          = readkey_dbl('params.txt','rt',   min(3600.d0,par%tstop),    1200.d0,    7200.d0 )
            par%dtbc        = readkey_dbl('params.txt','dtbc',          1.0d0,      0.1d0,      2.0d0   )
         endif
         !
         if (par%wbctype==WBCTYPE_SWAN) then
            par%dthetaS_XB  = readkey_dbl ('params.txt','dthetaS_XB',   0.0d0,      -360.d0,    360.0d0,strict=.true. )
         endif
         ! 
         par%wbcversion = readkey_int ('params.txt','wbcversion', 3, 1, 3,strict=.true.,silent=.true.)  ! wbcversion defaults to 3
         !
         if (par%wbcversion>2) then
            par%nspectrumloc    = readkey_int ('params.txt','nspectrumloc',   1,          1,       par%ny+1 )
         else
            par%nspectrumloc = 1
         endif
      endif
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
      if (par%nonhspectrum==1) then
         call parmapply('front',5,par%front,par%front_str)
      else
         if (par%ny>2) then
            call parmapply('front',2,par%front,par%front_str)
         else
            call parmapply('front',1,par%front,par%front_str)
         endif
      endif
      !
      ! left and right
      call setallowednames('neumann',   LR_NEUMANN,    &
      'wall'   ,   LR_WALL,       &
      'no_advec',  LR_NO_ADVEC,   &
      'neumann_v', LR_NEUMANN_V,  &
      'abs_1d',    LR_ABS_1D)
      call setoldnames('0','1')
      call parmapply('left',1,par%left,par%left_str)

      call parmapply('right',1,par%right,par%right_str)

      ! back
      call setallowednames('wall',    BACK_WALL,     &
      'abs_1d',  BACK_ABS_1D,   &
      'abs_2d',  BACK_ABS_2D,   &
      'wlevel',  BACK_WLEVEL)
      call setoldnames('0','1','2','3')
      if (par%ny>2) then
         call parmapply('back',3,par%back,par%back_str)
      else
         call parmapply('back',2,par%back,par%back_str)
      endif
      
      ! Error for wlevel
      if ((par%front == FRONT_WLEVEL) .or. (par%back == BACK_WLEVEL)) then
         call writelog('lse','','Error: wlevel no longer supported. Change front and/or back boundary condition.')
         call halt_program
      endif
      
      ! others
      par%ARC         = readkey_int ('params.txt','ARC',      1,              0,       1       ,strict=.true.)
      par%order       = readkey_dbl ('params.txt','order',    2.d0,           1.d0,    2.d0    ,strict=.true.)
      par%highcomp    = readkey_int ('params.txt','highcomp',    0,           0,       1                     )
      par%freewave    = readkey_int ('params.txt','freewave', 0,              0,       1       ,strict=.true.)
      par%epsi        = readkey_dbl ('params.txt','epsi',     -1.d0,          -1.d0,   0.2d0   )
      par%nc          = readkey_int ('params.txt','nc',       par%ny+1,       1,       par%ny+1,strict=.true.,silent=.true.)

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
         if (par%wavemodel == WAVEMODEL_STATIONARY) then
            call parmapply('break',2,par%break,par%break_str) ! default: baldock
            par%gamma    = readkey_dbl ('params.txt','gamma',   0.78d0,     0.4d0,     0.9d0)   
            par%gammax   = readkey_dbl ('params.txt','gammax',   0.6d0,      .4d0,      5.d0)   
         elseif (par%wavemodel == WAVEMODEL_SURFBEAT) then
            call parmapply('break',3,par%break,par%break_str) ! default: roelvink2
            par%gamma    = readkey_dbl ('params.txt','gamma',   0.55d0,     0.4d0,     0.9d0)   
            par%gammax   = readkey_dbl ('params.txt','gammax',   2.d0,      .4d0,      5.d0)  
         else
         endif
         if (par%break == BREAK_ROELVINK_DALY) then
            par%gamma2   = readkey_dbl ('params.txt','gamma2',   0.3d0,     0.0d0,     0.5d0)
         endif
         par%alpha        = readkey_dbl ('params.txt','alpha',   1.0d0,     0.5d0,     2.0d0)
         par%n            = readkey_dbl ('params.txt','n',       10.0d0,     5.0d0,    20.0d0) 
         par%delta        = readkey_dbl ('params.txt','delta',   0.0d0,     0.0d0,     1.0d0,strict=.true.)
         par%wavfriccoef  = readkey_dbl ('params.txt','fw',       0.d0,   0d0,      1.0d0)
         ! try to read a wave friction file
         par%wavfricfile  = readkey_name('params.txt','fwfile')
         if (par%wavfricfile .ne. ' ') then
            call check_file_exist(par%wavfricfile)
            if (par%gridform==GRIDFORM_XBEACH) then
               call check_file_length(par%wavfricfile,par%nx+1,par%ny+1)
            endif
            call writelog('lws','(a,a,a)','Warning: wave friction coefficient values from file ''',&
            trim(par%wavfricfile), &
            ''' will be used in computation')
         end if
         par%fwcutoff = readkey_dbl ('params.txt','fwcutoff',  1000.d0,   0d0,      1000.d0)
         par%breakerdelay = readkey_dbl ('params.txt','breakerdelay',    1.d0,   0.d0,      3.d0,strict=.true.)
         par%shoaldelay = readkey_int ('params.txt','shoaldelay',    0,   0,      1,silent=.true.,strict=.true.)
         if (par%shoaldelay==1) then
            par%facsd      = readkey_dbl ('params.txt','facsd',       1.d0,   0d0,      2.0d0)
         endif
         if (par%swrunup==1) then
            par%facrun     = readkey_dbl ('params.txt','facrun',      1.d0,   0d0,      2.0d0)
         endif
         !
         ! Roller parameters
         call writelog('l','','--------------------------------')
         call writelog('l','','Roller parameters: ')
         par%roller   = readkey_int ('params.txt','roller',     1,        0,     1,strict=.true.)
         par%rfb      = readkey_int ('params.txt','rfb',        0,        0,     1,strict=.true.)
         !
         ! Wave-current interaction parameters
         call writelog('l','','--------------------------------')
         call writelog('l','','Wave-current interaction parameters: ')
         par%wci      = readkey_int ('params.txt','wci',        0,        0,     1,strict=.true.)
         par%hwci     = readkey_dbl ('params.txt','hwci',   0.1d0,   0.001d0,      1.d0)
         par%hwcimax  = readkey_dbl ('params.txt','hwcimax',   100.d0,   0.01d0,      100.d0)
         par%cats     = readkey_dbl ('params.txt','cats',   4.d0,     1.d0,      50.d0)
      endif
      !
      ! Flow parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Flow parameters: ')

      call setallowednames('chezy',     BEDFRICTION_CHEZY,  &
           'cf',                        BEDFRICTION_CF, &
           'white-colebrook',           BEDFRICTION_WHITE_COLEBROOK, &
           'manning',                   BEDFRICTION_MANNING, &
           'white-colebrook-grainsize', BEDFRICTION_WHITE_COLEBROOK_GRAINSIZE)
      if(par%useXBeachGSettings==0) then
         call parmapply('bedfriction',1,par%bedfriction,par%bedfriction_str)
      else
         call parmapply('bedfriction',5,par%bedfriction,par%bedfriction_str)
      endif

      ! prevent using bed friction files without explicitly setting bedfriction type
      if (.not. isSetParameter('params.txt','bedfriction') .and. &
           isSetParameter('params.txt','bedfricfile')) then
         call writelog('lswe','','The use of keyword BEDFRICFILE without keyword BEDFRICTION is not allowed')
         call writelog('lswe','','Terminating simulation')
         call halt_program
      endif

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
      par%maxcf   = readkey_dbl ('params.txt','maxcf',     0.04d0,     0.0d0,   1.0d0)  ! max cf, only used for Manning and White Colebrook (Chezy: 15)
      par%nuh     = readkey_dbl ('params.txt','nuh',       0.1d0,     0.0d0,   1.0d0)
      par%nuhfac  = readkey_dbl ('params.txt','nuhfac',    1.0d0,     0.0d0,   1.0d0)
      par%nuhv    = readkey_dbl ('params.txt','nuhv',      1.d0,      1.d0,    20.d0,silent=.true.)
      par%smag    = readkey_int ('params.txt','smag',      1,         0,       1,strict=.true.)
      
      if (par%gwflow ==1) then
         if(par%useXBeachGSettings==1) then
            ! default on in XBeach-G!
            par%friction_infiltration  = readkey_int('params.txt','friction_infiltration',1,0,1, strict=.true.) 
         else
            par%friction_infiltration  = readkey_int('params.txt','friction_infiltration',0,0,1, silent=.true.,strict=.true.) 
         endif
      else
         par%friction_infiltration  = 0
      endif
      par%friction_turbulence    = readkey_int('params.txt','friction_turbulence'  ,0,0,1, silent=.true.,strict=.true.)
      if (par%friction_turbulence==1) then 
         par%gamma_turb             = readkey_dbl ('params.txt','gamma_turb', 1.d0, 0.d0, 2.d0)
      endif
      
      call setallowednames('none', CF_ACC_NONE, &
                           'mccall', CF_ACC_MCCALL,  &
                           'nielsen',CF_ACC_NIELSEN)
      if(par%useXBeachGSettings==1) then
         ! Default should be MCCALL in XBeach-G
         call parmapply('friction_acceleration',2,par%friction_acceleration,par%friction_acceleration_str)  
      else
         call parmapply('friction_acceleration',1,par%friction_acceleration,par%friction_acceleration_str,silent=.true.)
      endif
      
      if(par%friction_acceleration==CF_ACC_MCCALL) then
         par%ci          = readkey_dbl ('params.txt','ci',1.0d0,    0.5d0,   1.5d0) 
      elseif (par%friction_acceleration==CF_ACC_NIELSEN) then
         par%phit        = readkey_dbl ('params.txt','phit',25.0d0,    0.00d0,  90.0d0) 
      endif
      
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
         if(par%useXBeachGSettings==1) then
            par%kx         = readkey_dbl ('params.txt','kx'        , 0.01d0 , 0.001d0, 0.9d0)
            par%ky         = readkey_dbl ('params.txt','ky'        , par%kx , 0.001d0, 0.9d0)
            par%kz         = readkey_dbl ('params.txt','kz'        , par%kx , 0.001d0, 0.9d0)
         else
            par%kx         = readkey_dbl ('params.txt','kx'        , 0.0001d0 , 0.00001d0, 0.1d0)
            par%ky         = readkey_dbl ('params.txt','ky'        , par%kx   , 0.00001d0, 0.1d0)
            par%kz         = readkey_dbl ('params.txt','kz'        , par%kx   , 0.00001d0, 0.1d0)
         endif
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
         par%gwnonh     = readkey_int ('params.txt','gwnonh',par%useXBeachGSettings,   0,    1,strict=.true.)
         if (par%gwnonh==1) then
            if (par%ny>2) then
               par%gwfastsolve = readkey_int ('params.txt','gwfastsolve',      0,    0,      1,silent=.true.,strict=.true.)
            endif
         endif

         ! Type of momentum equation
         call setallowednames('laminar',    GWSCHEME_LAMINAR,  &
         'turbulent',  GWSCHEME_TURBULENT)
         call setoldnames('darcy','modflow')
         if(par%useXBeachGSettings==1) then
            call parmapply('gwscheme',2,par%gwscheme,par%gwscheme_str)
         else
            call parmapply('gwscheme',1,par%gwscheme,par%gwscheme_str)
         endif

         if (par%gwscheme==GWSCHEME_TURBULENT) then
            par%gwReturb    = readkey_dbl ('params.txt','gwReturb'   , 100.d0    , 1.d0     , 600.d0)
         endif
         call setallowednames('parabolic',       GWHEADMODEL_PARABOLIC,   &
         'exponential',     GWHEADMODEL_EXPONENTIAL)
         call parmapply('gwheadmodel',1,par%gwheadmodel,par%gwheadmodel_str)

         par%gwhorinfil = readkey_int ('params.txt','gwhorinfil',      0,           0,        1,strict=.true.,silent=.true.)
      endif
      !
      ! Non-hydrostatic correction parameters
      if (par%wavemodel==WAVEMODEL_NONH) then
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
         par%Topt         = readkey_dbl('params.txt','Topt',  10.d0, 1.d0, 20.d0)
         par%nonhq3d      = readkey_int('params.txt','nonhq3d' ,0,0,1,silent=.true.,strict=.true.)
         if (par%nonhq3d==1) then
            par%nhbreaker    = readkey_int('params.txt','nhbreaker' ,1,0,1,strict=.true.)
            par%nhlay        = readkey_dbl('params.txt','nhlay' ,0.33d0,0.d0,1.d0)
         else
            if(par%useXBeachGSettings==0) then
               par%nhbreaker    = readkey_int('params.txt','nhbreaker' ,1,0,1,strict=.true.)
            else
               par%nhbreaker    = readkey_int('params.txt','nhbreaker' ,1,0,3,strict=.true.)
            endif
            par%dispc        = readkey_dbl('params.txt','dispc' ,-1.0d0,0.1d0,2.0d0)
         endif

         if (par%nhbreaker==1) then
            par%breakviscfac = readkey_dbl('params.txt','breakviscfac',1.5d0, 1.d0, 3.d0)
            par%maxbrsteep   = readkey_dbl('params.txt','maxbrsteep',0.4d0, 0.3d0, 0.8d0)
            par%reformsteep  = readkey_dbl('params.txt','reformsteep',0.25d0*par%maxbrsteep,0.d0,0.95d0*par%maxbrsteep)
         elseif (par%nhbreaker==2) then
            par%breakvisclen = readkey_dbl('params.txt','breakvisclen',1.0d0, 0.75d0, 3.d0)
            par%maxbrsteep   = readkey_dbl('params.txt','maxbrsteep',0.4d0, 0.3d0, 0.8d0)
            par%secbrsteep   = readkey_dbl('params.txt','secbrsteep',0.5d0*par%maxbrsteep,0.d0,0.95d0*par%maxbrsteep)
         elseif (par%nhbreaker==3) then
            !par%avis    = readkey_dbl ('params.txt','avis',       0.3d0,     0.0d0,   1.0d0)
            par%breakvisclen = readkey_dbl('params.txt','breakvisclen',1.0d0, 0.75d0, 3.d0)
            par%maxbrsteep   = readkey_dbl('params.txt','maxbrsteep',0.6d0, 0.3d0, 0.8d0)
            par%secbrsteep   = readkey_dbl('params.txt','secbrsteep',0.5d0*par%maxbrsteep,0.d0,0.95d0*par%maxbrsteep)
         endif
      endif
      !
      !
      ! Sediment transport parameters
      if (par%sedtrans==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Sediment transport parameters: ')

         call setallowednames('soulsby_vanrijn',    FORM_SOULSBY_VANRIJN,  &
         'vanthiel_vanrijn',   FORM_VANTHIEL_VANRIJN, &
         'vanrijn1993', FORM_VANRIJN1993, &
         'nielsen2006', FORM_NIELSEN2006, &
         'mccall_vanrijn', FORM_MCCALL_VANRIJN, &
         'wilcock_crow', FORM_WILCOCK_CROW, &
         'engelund_fredsoe', FORM_ENGELUND_FREDSOE, &
         'mpm', FORM_MPM, &
         'wong_parker', FORM_WONG_PARKER, &
         'fl_vb', FORM_FL_VB, &
         'fredsoe_deigaard', FORM_FREDSOE_DEIGAARD)
         call setoldnames('1','2','3','4','5','6','7','8','9','10','11')
         if(par%useXBeachGSettings==0) then
            call parmapply('form',2,par%form, par%form_str)
         else
            call parmapply('form',5,par%form, par%form_str)
         endif
         
         if (par%wavemodel==WAVEMODEL_STATIONARY .or. par%wavemodel==WAVEMODEL_SURFBEAT) then
            call setallowednames('ruessink_vanrijn',  WAVEFORM_RUESSINK_VANRIJN,  &
                                 'vanthiel',          WAVEFORM_VANTHIEL)
            call setoldnames('1','2')
            call parmapply('waveform',2,par%waveform,par%waveform_str)
            par%sws      = readkey_int ('params.txt','sws',           1,        0,     1,strict=.true.)
            par%lws      = readkey_int ('params.txt','lws',           1,        0,     1,strict=.true.)
            par%BRfac    = readkey_dbl ('params.txt','BRfac',    1.0d0,       0.d0, 1.d0)
            par%facua    = readkey_dbl ('params.txt','facua  ',0.10d0,    0.00d0,   1.0d0)
            par%facSk    = readkey_dbl ('params.txt','facSk  ',par%facua,    0.00d0,   1.0d0)
            par%facAs    = readkey_dbl ('params.txt','facAs  ',par%facua,    0.00d0,   1.0d0)
            if (par%waveform == WAVEFORM_VANTHIEL) then
               par%Tbfac    = readkey_dbl ('params.txt','Tbfac  ',1.0d0,     0.00d0,   1.0d0)
            endif
            call setallowednames('none',              TURB_NONE,           &
                                 'wave_averaged',     TURB_WAVE_AVERAGED,  &
                                 'bore_averaged',     TURB_BORE_AVERAGED)
            call setoldnames('0','1','2')
            call parmapply('turb',3,par%turb,par%turb_str)
            if (par%waveform == WAVEFORM_RUESSINK_VANRIJN) then
               call parmapply('turb',2,par%turb,par%turb_str)
            endif
            call setallowednames('none',       TURBADV_NONE,        &
                                 'lagrangian', TURBADV_LAGRANGIAN,  &
                                 'eulerian',   TURBADV_EULERIAN)
            call parmapply('turbadv',1,par%turbadv,par%turbadv_str)
            
         else
            par%sws      = 0
            par%lws      = 1
            par%BRfac    = 0.d0
            par%facua    = 0.d0
            par%facSk    = 0.d0
            par%facAs    = 0.d0
            par%turb     = TURB_WAVE_AVERAGED     ! note, needed in case we want to use lwt ... Better solution still waiting
            par%turbadv  = TURBADV_LAGRANGIAN
         endif
         
         ! Parameters for gravel-type equations
         if(par%form==FORM_NIELSEN2006 .or. &
            par%form==FORM_MCCALL_VANRIJN .or. &
            par%form==FORM_WILCOCK_CROW .or. &
            par%form==FORM_ENGELUND_FREDSOE .or. &
            par%form==FORM_MPM .or. &
            par%form==FORM_WONG_PARKER .or. &
            par%form==FORM_FL_VB .or. &
            par%form==FORM_FREDSOE_DEIGAARD) then
            
            if(par%form==FORM_NIELSEN2006) then
               par%Ctrans   = readkey_dbl ('params.txt','Ctrans',12d0,  0d0,  30d0)
               par%thetcr   = readkey_dbl ('params.txt','thetcr ',0.05d0,    0.00d0,   1.0d0)
               par%phaselag = readkey_int ('params.txt','phaselag',    0,    0,   1)
               if (par%phaselag==1) then
                  par%phit     = readkey_dbl ('params.txt','phit   ',25.00d0,    0.00d0,  90.0d0) 
               endif
               call setallowednames('constant',          SEDFRICFAC_CONSTANT,     &
                                    'flowfricfac',       SEDFRICFAC_FLOWFRIC,  &
                                    'nielsen',           SEDFRICFAC_NIELSEN, &
                                    'swart',             SEDFRICFAC_SWART, &
                                    'wilson',            SEDFRICFAC_WILSON)
               call setoldnames('0','1','2','3','4')
               call parmapply('sedfricfac',1,par%sedfricfac,par%sedfricfac_str)
               if (par%sedfricfac==SEDFRICFAC_SWART) then
                  par%Arms     = readkey_dbl ('params.txt','Arms', 10.d0, 0.d0, 100.d0) 
               endif
               if (par%sedfricfac==SEDFRICFAC_CONSTANT) then
                  par%fsed     = readkey_dbl ('params.txt','fsed', 0.025d0, 0.005d0, 0.05d0) 
               endif
               par%streaming  = readkey_int ('params.txt','streaming', 0,   0,   1,silent=.true.,strict=.true.)
            else ! all other transport equations
               par%thetcr   = readkey_dbl ('params.txt','thetcr ',-10.d0, 0.00d0, 1.0d0,silent=.true.) ! negative = self-compute
               par%sedfricfac = SEDFRICFAC_FLOWFRIC ! not actually used in the equations, all comes from taubx
               par%uprushfac    = readkey_dbl ('params.txt','uprushfac', 1.d0, 0.d0, 3.d0)
               par%backwashfac  = readkey_dbl ('params.txt','backwashfac', 1.d0, 0.d0, 3.d0)
               par%incldzdx     = readkey_int('params.txt','incldzdx',0,0,1,strict=.true.,silent=.true.)
               if (par%gwflow==1) then 
                  par%inclrelweight = readkey_int('params.txt','inclrelweight',1,0,1,strict=.true.) 
               else
                  par%inclrelweight = 0
               endif
            endif
            
            call setallowednames('none',             SLOPECORR_NONE,     &
                                 'nielsen',          SLOPECORR_NIELSEN,  &
                                 'fredsoe_deigaard', SLOPECORR_FREDSOE_DEIGAARD)
            call setoldnames('0','1','2')
            call parmapply('slopecorr',3,par%slopecorr,par%sedfricfac_str)
            
            par%sus      = readkey_int ('params.txt','sus    ',0,           0,            0,strict=.true.)
            par%bed      = readkey_int ('params.txt','bed    ',1,           1,            1,strict=.true.)
            par%bulk     = readkey_int ('params.txt','bulk   ',0,           0,            0,strict=.true.)
            if(par%bulk==0) then
               par%facsl = 0.d0
            else
               par%facsl    = readkey_dbl ('params.txt','facsl  ',0.0d0,       0.d0, 1.6d0)
            endif
            par%bdslpeffmag = BDSLPEFFMAG_NONE
            par%bdslpeffini = BDSLPEFFINI_NONE
            par%bdslpeffdir = BDSLPEFFDIR_NONE
            par%smax = -1.d0
         else ! non-gravel transport equations
            par%sus      = readkey_int ('params.txt','sus    ',1,           0,            1,strict=.true.)
            par%bed      = readkey_int ('params.txt','bed    ',1,           0,            1,strict=.true.)
            par%bulk     = readkey_int ('params.txt','bulk   ',0,           0,            1,strict=.true.)
            par%facsl    = readkey_dbl ('params.txt','facsl  ',  0.15d0,       0.d0, 1.6d0)   
            par%z0       = readkey_dbl ('params.txt','z0     ',0.006d0,    0.0001d0,   0.05d0)
            par%smax     = readkey_dbl ('params.txt','smax',   -1.d0,    -1.d0,   3.d0)       !changed 28/11 and back 10/2
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
         if(par%useXBeachGSettings==0) then
            par%reposeangle         = readkey_dbl ('params.txt','reposeangle',  30.d0,     0.d0,     45.d0)
         else
            par%reposeangle         = readkey_dbl ('params.txt','reposeangle',  35.d0,     20.d0,     60.d0)
         endif
         par%tsfac    = readkey_dbl ('params.txt','tsfac',   0.1d0,    0.01d0,   1.d0)
         par%Tsmin    = readkey_dbl ('params.txt','Tsmin  ',0.5d0,     0.01d0,   10.d0)
         par%facDc    = readkey_dbl ('params.txt','facDc  ',1.0d0,     0.00d0,   1.0d0)
         par%jetfac   = readkey_dbl ('params.txt','jetfac ',0.0d0,     0.00d0,   1.0d0,silent=.true.)
         par%lwt      = readkey_int ('params.txt','lwt    ',0,           0,            1,strict=.true.)
         if(par%lwt==1 .and. par%form==FORM_NIELSEN2006) then
            par%yturb    = readkey_dbl ('params.txt','yturb', 1.d0, 0.d0, 1.d0)
            par%facthr   = readkey_dbl ('params.txt','facthr',  1.0d0,    1.d0,    100.d0)
         endif
         par%betad    = readkey_dbl ('params.txt','betad  ',1.0d0,     0.00d0,   10.0d0)
         par%fallvelred          = readkey_int ('params.txt','fallvelred',    0,        0,     1,strict=.true.)
         par%dilatancy           = readkey_int ('params.txt','dilatancy', 0,    0,     1,strict=.true.)
         if (par%dilatancy==1) then
            par%rheeA               = readkey_dbl ('params.txt','rheeA',         0.75d0,   0.75d0, 2.d0) ! Between 3/4 and 1/(1-n), see paper Van Rhee (2010)
            par%pormax              = readkey_dbl ('params.txt','pormax',        0.5d0,    0.3d0, 0.6d0)
         endif
         par%bermslopetransport = readkey_int ('params.txt','bermslopetransport', 0, 0, 1,strict=.true.,silent=.true.)
         if (par%bermslopetransport==1) then
            par%bermslopebed = readkey_int ('params.txt','bermslopebed', 0, 0,  1,strict=.true.)
            par%bermslopesus = readkey_int ('params.txt','bermslopesus', 0, 0,  1,strict=.true.)        
            par%bermslope   = readkey_dbl ('params.txt','bermslope ',0.1d0,     0.00d0,   1.0d0)
            par%bermslopefac = readkey_dbl ('params.txt','bermslopefac ',15.0d0, 0.00d0, 30.0d0)
            if (par%wavemodel==WAVEMODEL_SURFBEAT) then 
               par%bermslopegamma = readkey_dbl ('params.txt','bermslopegamma ',1.0d0, 0.75d0, 1.5d0)
               par%bermslopedepth = readkey_dbl ('params.txt','bermslopedepth ',0.0d0, 0.0d0, 0.5d0)
            else
               par%bermslopedepth = readkey_dbl ('params.txt','bermslopedepth ',1.0d0, 0.5d0, 1.5d0)
            endif
         else
            par%bermslopebed = 0
            par%bermslopesus = 0
         endif
         
      endif
      !
      ! Q3D sediment transport parameters
      ! -> also make this part active if van Rijn, 1993 is used as ceq solver
      if (par%q3d==1 .or. par%form==FORM_VANRIJN1993) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Q3D sediment transport parameters: ')
         par%vonkar  = readkey_dbl ('params.txt','vonkar',   0.4d0,     0.01d0,  1.d0)
         par%vicmol  = readkey_dbl ('params.txt','vicmol',   0.000001d0,   0.d0,    0.001d0)
         par%kmax    = readkey_int ('params.txt','kmax ',      100,        1,        1000)
         par%nz      = par%kmax     ! work-around; discuss with Dano later difference nz / kmax
         par%sigfac  = readkey_dbl ('params.txt','sigfac ',1.3d0,     0.00d0,   10.d0)
         par%deltar  = readkey_dbl ('params.txt','deltar ',0.025d0,     0.001d0,   1.0d0)
         par%rwave  = readkey_dbl ('params.txt','rwave ',2.0d0,     0.1d0,   10.0d0)
      else ! workaround, needed for array allocation
         par%kmax = 1
      endif
      !
      ! Bed composition parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Bed composition parameters: ')
      par%ngd      = readkey_int ('params.txt','ngd',        1,           1,        20,strict=.true.)
      par%nd       = readkey_int ('params.txt','nd ',        3,           3,        1000,strict = .true.)
      par%por      = readkey_dbl ('params.txt','por',    0.4d0,       0.3d0,    0.5d0)
      if (par%dilatancy==1) then
         par%D15      = readkey_dblvec('params.txt','D15',par%ngd,size(par%D15),0.00015d0,0.00001d0,0.0008d0) ! Lodewijk
      endif
      if(par%useXBeachGSettings==1) then
         par%D50      = readkey_dblvec('params.txt','D50',par%ngd,size(par%D50),0.010d0,0.002d0,0.08d0)
         par%D90      = readkey_dblvec('params.txt','D90',par%ngd,size(par%D90),0.015d0,0.003d0,0.12d0)
      else
         par%D50      = readkey_dblvec('params.txt','D50',par%ngd,size(par%D50),0.0002d0,0.00005d0,0.0008d0)
         par%D90      = readkey_dblvec('params.txt','D90',par%ngd,size(par%D90),0.0003d0,0.00010d0,0.0015d0)
      endif
         
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
      ! Morphology parameters
      if (par%morphology==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Morphology parameters: ')
         par%morfac   = readkey_dbl ('params.txt','morfac', 1.0d0,        0.d0,  1000.d0)
         par%morfacopt= readkey_int ('params.txt','morfacopt', 1,        0,        1,strict=.true.)
         par%morstart = readkey_dbl ('params.txt','morstart',0.0d0,      0.d0, par%tstop)
         par%morstop  = readkey_dbl ('params.txt','morstop', par%tstop,      0.d0, 10000000.d0)
         if(par%useXBeachGSettings==1) then
            par%wetslp   = readkey_dbl ('params.txt','wetslp', tan(par%reposeangle/180.d0*4.d0*atan(1.d0)),0.1d0,1.d0)
            par%dryslp   = readkey_dbl ('params.txt','dryslp', par%wetslp,       0.1d0,     2.d0)
         else
            par%wetslp   = readkey_dbl ('params.txt','wetslp', 0.3d0,       0.1d0,     1.d0)
            par%dryslp   = readkey_dbl ('params.txt','dryslp', 1.0d0,       0.1d0,     2.d0)
         endif
         if (par%ny == 0) then
            par%lsgrad   = readkey_dbl ('params.txt','lsgrad', 0.0d0,       -0.1d0,     .1d0,silent=.true.)
         else
            par%lsgrad = 0.0d0 
         endif
         par%hswitch  = readkey_dbl ('params.txt','hswitch',0.1d0,      0.01d0,    1.0d0)
         par%dzmax    = readkey_dbl ('params.txt','dzmax  ',0.05d0,    0.00d0,   1.0d0)
         par%struct   = readkey_int ('params.txt','struct ',0    ,      0,             1,strict=.true.)
         if (par%struct==1) then
            par%ne_layer = readkey_name('params.txt','ne_layer')
            call check_file_exist(par%ne_layer)
            if (par%gridform==GRIDFORM_XBEACH) then
               call check_file_length(par%ne_layer,par%nx+1,par%ny+1)
            endif
         endif
         if (par%swrunup==1 .and. par%struct==0) then
            call writelog('lws','(a)', &
            'Warning: swrunup can only be used in combination with struct=1. swrunup will be turned off.')
            par%swrunup = 0;
         endif
      endif
      !
      !
      ! Output variables
      call writelog('l','','--------------------------------')
      call writelog('l','','Output variables: ')
      par%timings  = readkey_int ('params.txt','timings',      1,       0,      1,strict=.true.)
      testc = readkey_name('params.txt','tunits')
      if (len(trim(testc)) .gt. 0) par%tunits = trim(testc)
      par%tstart  = readkey_dbl ('params.txt','tstart',   0.d0,      0.d0,par%tstop)
      par%tint    = readkey_dbl ('params.txt','tint',     1.d0,     .01d0, par%tstop-par%tstart)  ! Robert
      par%tsglobal = readkey_name('params.txt','tsglobal')
      if (par%tsglobal==' ') then
         par%tintg   = readkey_dbl ('params.txt','tintg', par%tint,     .01d0, par%tstop-par%tstart)  ! Robert
      endif
      par%tspoints = readkey_name('params.txt','tspoints')
      if (par%tspoints==' ') then
         par%tintp   = readkey_dbl ('params.txt','tintp', par%tint,     .01d0, par%tstop-par%tstart)  ! Robert
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
      par%npointvar   = readkey_int ('params.txt','npointvar',   0,  0, 50)
      call readpointvars(par)
      par%nrugdepth   = readkey_int('params.txt','nrugdepth',1,1,10)
      par%rugdepth    = readkey_dblvec('params.txt','rugdepth',par%nrugdepth,size(par%rugdepth),0.0d0,0.0d0,0.1d0)

      ! mean output
      par%nmeanvar    = readkey_int ('params.txt','nmeanvar'  ,  0,  0, 15)
      call readmeans(par)

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
      if(par%outputformat==OUTPUTFORMAT_NETCDF .and. par%useXBeachGSettings==0) then
         par%remdryoutput = readkey_int ('params.txt','remdryoutput',     1,  0, 1,strict=.true.)
      else
         par%remdryoutput = readkey_int ('params.txt','remdryoutput',     0,  0, 1,strict=.true.)
      endif
      !
      !
      ! Output projection
      call writelog('l','','--------------------------------')
      call writelog('l','','Output projection: ')
      testc = readkey_name('params.txt','projection')
      if (len(trim(testc)) .gt. 0) par%projection = testc
      par%rotate = readkey_int ('params.txt','rotate',     1,  0, 1,strict=.true.)
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
         call check_file_exist(par%veggiefile)
         par%veggiemapfile = readkey_name  ('params.txt', 'veggiemapfile'                    )
         call check_file_exist(par%veggiemapfile)
         par%Trep          = readkey_dbl   ('params.txt','Trep',     1.d0,   0.01d0,    20.d0)
         par%vegnonlin     = readkey_int   ('params.txt', 'vegnonlin',0,0,1,silent=.true.)
         par%vegcanflo     = readkey_int   ('params.txt', 'vegcanflo',0,0,1,silent=.true.)
         par%veguntow      = readkey_int   ('params.txt', 'veguntow', 1,0,1,silent=.true.)
         par%porcanflow    = readkey_int   ('params.txt', 'porcanflow', 0,0,1)
         if (par%porcanflow ==1) then
             par%Kp             = readkey_dbl ('params.txt','Kp',    0.0001d0,      0.d0,     1.d0)
             par%Cm             = readkey_dbl ('params.txt','Cm',    1.0d0,         0.d0,     2.d0)
         endif
         ! veggiefile routine should set nveg
      endif
      !
      !
      ! Wave numerics parameters
      call writelog('l','','--------------------------------')
      call writelog('l','','Wave numerics parameters: ')

      call setallowednames('upwind_1',      SCHEME_UPWIND_1,      &
      'lax_wendroff',  SCHEME_LAX_WENDROFF,  &
      'upwind_2',      SCHEME_UPWIND_2,      &
      'warmbeam',      SCHEME_WARMBEAM)
      call setoldnames('1','2','3','4')
      call parmapply('scheme',4,par%scheme,par%scheme_str)

      if (par%wavemodel == WAVEMODEL_STATIONARY .or. par%single_dir==1) then
         par%wavint     = readkey_dbl ('params.txt','wavint',    60.d0,      1.d0,  3600.d0)
         if (par%single_dir==1) then
            par%maxerror   = readkey_dbl ('params.txt','maxerror', 0.005d0, 0.00001d0, 0.001d0)
         else
            par%maxerror   = readkey_dbl ('params.txt','maxerror', 0.0005d0, 0.00001d0, 0.001d0)
         endif
         par%maxiter    = readkey_int ('params.txt','maxiter',    500,         2,      1000)
      endif
      ! only default to Snell's Law if in 1D and only one directional bin
      if (par%single_dir == 0) then
         if (par%ny==0 .and. nint(abs(par%thetamax-par%thetamin)/par%dtheta)<2) then
            par%snells      = readkey_int ('params.txt','snells',    1,        0,     1,strict=.true.)
         else
            par%snells      = readkey_int ('params.txt','snells',    0,        0,     1,strict=.true.)
         endif
      else
         par%snells      = readkey_int ('params.txt','snells',    0,        0,     1,silent=.true.,strict=.true.)
      endif
      if(par%nonhspectrum == 0) then
         par%swkhmin = readkey_dbl ('params.txt','swkhmin',    -0.01d0, -0.01d0, 0.35d0,silent=.true.) ! Robert for Daan
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
      par%secorder = readkey_int('params.txt','secorder' ,par%nonh,0,1,strict=.true.)
      par%oldhu    = readkey_int('params.txt','oldhu' ,0,0,1,silent=.true.,strict=.true.)
      !
      !
      ! Sediment transport numerics parameters
      if (par%sedtrans==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Sediment transport numerics parameters: ')
         par%thetanum   = readkey_dbl ('params.txt','thetanum',   1.d0,    0.5d0,   1.d0)
         par%sourcesink = readkey_int ('params.txt','sourcesink    ',0,     0,         1,silent=.true.,strict=.true.)
         par%cmax       = readkey_dbl ('params.txt','cmax',      0.1d0,    0.0d0,   1.d0)
      endif
      !
      !
      ! Bed update numerics parameters
      if (par%morphology==1) then
         call writelog('l','','--------------------------------')
         call writelog('l','','Bed update numerics parameters: ')
         par%frac_dz   = readkey_dbl ('params.txt','frac_dz',   0.7d0,    0.5d0,   0.98d0)
         par%nd_var    = readkey_int ('params.txt','nd_var',     2,           1,        par%nd,strict=.true.)
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
      if (par%mpiboundary ==  MPIBOUNDARY_MAN) then
         par%mmpi= readkey_int('params.txt','mmpi',2,1,100)
         par%nmpi= readkey_int('params.txt','nmpi',4,1,100)
      endif
#endif
      !
      !
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
      ! Transport equations (gravel) that don't work in 2DH
      if(par%form==FORM_NIELSEN2006 .or. &
         par%form==FORM_MCCALL_VANRIJN .or. &
         par%form==FORM_WILCOCK_CROW .or. &
         par%form==FORM_ENGELUND_FREDSOE .or. &
         par%form==FORM_MPM .or. &
         par%form==FORM_WONG_PARKER .or. &
         par%form==FORM_FL_VB .or. &
         par%form==FORM_FREDSOE_DEIGAARD) then
            
         ! These transport equations only work in 1D
         if (par%ny>0) then
            call writelog('lswe','(a)','Error: The following sediment transport equations cannot be used in 2DH simulations:')
            call writelog('lswe','(a)','nielsen2006')
            call writelog('lswe','(a)','mccall_vanrijn')
            call writelog('lswe','(a)','wilcock_crow')
            call writelog('lswe','(a)','engelund_fredsoe')
            call writelog('lswe','(a)','mpm')
            call writelog('lswe','(a)','wong_parker')
            call writelog('lswe','(a)','fl_vb')
            call writelog('lswe','(a)','fredsoe_deigaard')
            call writelog('lswe','(a)','Change to other transport equation or 1D model') 
            call halt_program
         endif
      endif
      !
      !
      ! Fix input parameters for choosen depthscale
      if (par%depthscale .ne. 1.d0) then
         par%eps     = par%eps/par%depthscale
         par%hmin    = par%hmin/par%depthscale
         par%hswitch = par%hswitch/par%depthscale
         par%dzmax   = par%dzmax/par%depthscale**1.5d0
         par%maxerror= par%maxerror/par%depthscale

         call writelog('lws','(a)','Warning: input parameters eps, hmin, hswitch and dzmax are scaled with')
         call writelog('lws','(a)','         depthscale to:')
         call writelog('lws','(a,f0.4)','eps = ',    par%eps)
         call writelog('lws','(a,f0.4)','hmin = ',   par%hmin)
         call writelog('lws','(a,f0.4)','hswitch = ',par%hswitch)
         call writelog('lws','(a,f0.4)','dzmax = ',  par%dzmax)
         call writelog('lws','(a,f0.4)','maxerror = ',  par%maxerror)

      endif
      !
      !
      ! Constants
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
      ! Wavemodel and wbctypes: warnings
      ! swave and nonh
      if (par%swave == 1 .and. par%wavemodel == WAVEMODEL_NONH) then
            call writelog('lsw','(a)','Warning: swave with wavemodel = nonh not possible')
            par%swave     = 0 
      endif
      !
      ! Wave model and wbctypes: errors (halt program)
      ! Bichromatic and not surfbeat
      if (par%wbctype == WBCTYPE_PARAMS .and. par%Tlong /= -123 .and. par%wavemodel /= WAVEMODEL_SURFBEAT) then
            call writelog('lwse','(a)','Error: bichromatic waves can only be used in surfbeat mode')
            call halt_program
      endif
      ! ts_nonh and not nonh
      if (par%wbctype == WBCTYPE_TS_NONH .and. par%wavemodel /= WAVEMODEL_NONH) then
            call writelog('lwse','(a)','Error: nonh time series can only be used in nonh mode')
            call halt_program
      endif
      ! ts_1 or ts_2
      if (par%wbctype == WBCTYPE_TS_1 .or. par%wbctype == WBCTYPE_TS_2) then
          if (par%wavemodel /= WAVEMODEL_SURFBEAT) then
                call writelog('lwse','(a)','Error: ts_1 and ts_2 time series can only be used in surfbeat mode')
                call halt_program
          endif
      endif
      ! errors for nonh = 1 and a wavemodel which is not wavemodel nonh
      if (par%nonh == 1 .and. par%wavemodel /= WAVEMODEL_NONH) then
            call writelog('lwse','(a)','Error: found both wavemodel and nonh=1. Define hydrodynamic mode with wavemodel')
            call halt_program
      endif
      ! errors for wavemodel = stationary and wbctype = reuse
      if (par%wbctype == WBCTYPE_REUSE .and. par%wavemodel == WAVEMODEL_STATIONARY) then
            call writelog('lwse','(a)','Error: wbctype = reuse cannot be used in stationary mode')
            call halt_program
      endif
      ! errors for wavemodel = stationary and wbctype = swan
      if (par%wbctype == WBCTYPE_SWAN .and. par%wavemodel == WAVEMODEL_STATIONARY) then
            call writelog('lwse','(a)','Error: wbctype = swan cannot be used in stationary mode')
            call halt_program
      endif
      ! errors for wavemodel = vardens and wbctype = swan
      if (par%wbctype == WBCTYPE_VARDENS .and. par%wavemodel == WAVEMODEL_STATIONARY) then
            call writelog('lwse','(a)','Error: wbctype = vardens cannot be used in stationary mode')
            call halt_program
      endif
      ! errors for wavemodel = vardens and wbctype = swan
      if (par%wbctype == WBCTYPE_PARAMETRIC .and. par%wavemodel == WAVEMODEL_STATIONARY) then
            call writelog('lwse','(a)','Error: wbctype = parametric cannot be used in stationary mode')
            call halt_program
      endif
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
      ! Only allow Baldock or Janssen in stationary mode and Roelvink in non-stationary
      if (par%swave==1) then
         if (par%wavemodel == WAVEMODEL_STATIONARY) then
            if (par%break .ne. BREAK_BALDOCK .and. par%break .ne. BREAK_JANSSEN) then
               if(par%break == BREAK_ROELVINK_DALY) then
                  call writelog('lwse','','Error: Roelvink-Daly formulations not implemented in stationary wave mode,')
                  call writelog('lwse','','         use Baldock or Janssen formulation.')
                  call halt_program
               else
                  call writelog('lwse','','Error: Roelvink formulations not allowed in stationary,')
                  call writelog('lwse','','         use Baldock or Janssen formulation.')
                  call halt_program
               endif
            endif
         else
            if (par%break == BREAK_BALDOCK) then
               call writelog('lws','','Warning: Baldock formulation not allowed in non-stationary, use break=roelvink1')
               call writelog('lws','','         or roelvink2 or roelvink_daly formulation.')
            endif
            if (par%break == BREAK_JANSSEN) then
               call writelog('lws','','Warning: Janssen formulation not allowed in non-stationary, use break=roelvink1 ')
               call writelog('lws','','         or roelvink2 or roelvink_daly formulation.')
            endif
         endif
      endif
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
      ! Check on pormax value
      if (par%dilatancy==1) then
         if (par%pormax.le.par%por) then
            call writelog('lws','','Warning: Maximum porosity for dilatancy effect smaller than actual porosity.')
            call writelog('lws','','         Setting pormax equal to por, effectively switching off dilatancy.')
            par%pormax = par%por
         end if
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
      ! is not using nonh with front = nonh_1d
      if (par%front==FRONT_NONH_1D) then  
        if (par%wavemodel /= WAVEMODEL_NONH) then
            call writelog('lws','','Warning: not recommended to use front = nonh_1d without wavemodel = nonh')
        endif
      endif
      !
      ! is using nonh without front = nonh_1d
      if (par%front/=FRONT_NONH_1D) then  
        if (par%wavemodel == WAVEMODEL_NONH) then
            call writelog('lws','','Warning: automatically changing front to nonh_1d with wavemodel = nonh')
            par%front = FRONT_NONH_1D
        endif
      endif
      ! is using wavemodel = stat without front = abs_1d
      if (par%front/=FRONT_ABS_1D) then  
        if (par%wavemodel == WAVEMODEL_STATIONARY ) then
            call writelog('lws','','Warning: automatically changing front to abs_1d with wavemodel = stationary')
            par%front = FRONT_ABS_1D
        endif
      endif
      !
      ! If using nonh, secorder should always be on
      if (par%wavemodel==WAVEMODEL_NONH) then
         if (par%secorder==0) then
            call writelog('lws','','Warning: Automatically turning on 2nd order correction in flow for')
            call writelog('lws','','         non-hydrostatic module [secorder=1]')
            par%secorder = 1
         endif
      endif

      !
      ! If using nonh, then the solver type is set by the grid size
      if (par%wavemodel==WAVEMODEL_NONH) then
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
      ! Lax-Wendroff not yet supported in curvilinear
      if (par%scheme==SCHEME_LAX_WENDROFF) then
         par%scheme=SCHEME_WARMBEAM
         call writelog('lws','','Warning: Lax Wendroff [scheme=lax_wendroff] scheme is not supported, changed')
         call writelog('lws','','         to Warming and Beam [scheme=SCHEME_WARMBEAM]')
      endif
      !
      ! Wave-current interaction with non-stationary waves still experimental
      if (par%wavemodel == WAVEMODEL_STATIONARY .and. par%wci==1) then
         call writelog('lws','','Warning: Wave-current interaction with non-stationary waves is still')
         call writelog('lws','','         experimental, continue with computation nevertheless')
      endif
      !
      ! Check for setting Snells law and single_dir
      if (par%single_dir == 1 .and. par%snells==1) then
         call writelog('lse', '', 'The options ''single_dir = 1'' and ''snells = 1'' are not compatible')
         call writelog('lse', '', 'Terminating simulation')
         call halt_program
      endif
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
      ! Mean output time fix
      if(par%tintm>(par%tstop-par%tstart) .and. par%nmeanvar>0) then
         call writelog('lws','','Warning: ''tintm'' is larger than output duration in the simulation.')
         call writelog('lws','','         Setting ''tintm'' = tstop-tstart = ',par%tstop-par%tstart,'s')
         par%tintm=par%tstop-par%tstart
      endif
      !
      !
      ! MPI domains
#ifdef USEMPI
      if (par%swave==1) then
         if ((par%wavemodel == WAVEMODEL_STATIONARY .or. par%single_dir==1) .and. par%ny > 0) then
            ! We need to set to mpiboundary = x to solve the stationary wave model.
            ! However, this requires ny>3*xmpi_osize
            if (par%mpiboundary .ne.  MPIBOUNDARY_MAN) then
               if(par%ny<=2*xmpi_size) then
                  call writelog('ewsl','','This simulation cannot be run in current MPI mode:')
                  call writelog('ewsl','','The stationary wave solver requires MPI subdivision by "x" (split ny).')
                  call writelog('ewsl','','The number of subdomains selected to run the model is ',xmpi_size,'.')
                  call writelog('ewsl','','The total number of grid cells in y (ny) is ',par%ny,'.')
                  call writelog('ewsl','(a,f0.2,a)','The number of cells per domain is ',dble(par%ny)/xmpi_size, &
                  ', which is less than the minimum value of 3')
                  call writelog('ewsl','','If you really (!) know what you''re doing, use "mpiboundary = man" and deal')
                  call writelog('ewsl','','with any unsatisfactory results')
                  call halt_program
               else
                  par%mpiboundary=MPIBOUNDARY_X
                  par%mpiboundary_str='x'
                  call writelog('wsl','','Changing mpiboundary to "x" for stationary wave model')
               endif
            else
               call writelog('wsl','','Warning: the stationary wave model only works with "mpiboundary=x"!')
            endif
         endif
      endif
#endif
      !
      !
      ! fix tint
      par%tint    = min(par%tintg,par%tintp,par%tintm)
      !
      !
      ! All input time frames converted to XBeach hydrodynamic time
      if (par%morfacopt==1) then
         par%tstart  = par%tstart / max(par%morfac,1.d0)
         par%tint    = par%tint   / max(par%morfac,1.d0)
         par%tintg   = par%tintg  / max(par%morfac,1.d0)
         par%tintp   = par%tintp  / max(par%morfac,1.d0)
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
      if (par%wbctype == WBCTYPE_REUSE) then
         if (par%wavemodel==WAVEMODEL_NONH) then
            dummystring='nhbcflist.bcf'
            call checkbcfilelength(par%tstop,par%wbctype,dummystring,filetype,nonh=.true.)
         else
            dummystring='ebcflist.bcf'
            call checkbcfilelength(par%tstop,par%wbctype,dummystring,filetype)
            dummystring='qbcflist.bcf'
            call checkbcfilelength(par%tstop,par%wbctype,dummystring,filetype)
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
      integer, allocatable                :: pointtypes(:) !  [-] Point types (0 = point, 1=rugauge)
      double precision, allocatable       :: xpointsw(:) ! world x-coordinate of output points
      double precision, allocatable       :: ypointsw(:) ! world y-coordinate of output points

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

      double precision, dimension(:),allocatable          :: xpointsw,ypointsw
      integer, dimension(:), allocatable       :: pointtypes
      integer                                  :: i,j
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
               if(par%nrugauge<=0) then
                  ! This branch of the else will change to a halt_program statement in later versions (written on 13 January 2011)
                  call writelog('lswe','','Point output must be specified using keyword ''npointvar''')
                  call writelog('lswe','','Stopping simulation')
                  call halt_program
               else
                  call writelog('lsw','','All point output will contain same data as rugauge output (x,y,zs).')
                  call writelog('lsw','','Other point output must be specified using keyword ''npointvar''')
               endif
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

         ! Now check all point and runup gauge names
         if(par%nrugauge+par%npoints>0) then
            do i=1,par%nrugauge+par%npoints-1
               do j=i+1,par%nrugauge+par%npoints
                  if(par%stationid(i)==par%stationid(j)) then
                     call writelog('lswe','','Duplicate names used for point station ID:')
                     call writelog('lswe','',par%stationid(i))
                     call writelog('lswe','','Stopping simulation')
                     call halt_program
                  endif
               enddo
            enddo
         endif

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
      double precision, dimension(:),intent(inout)        :: xpoints,ypoints

      character(slen)                          :: line,keyword,keyread,varline
      character(maxnamelen)                    :: varstr
      character(slen)                           :: fullline,errmes1,errmes2,okaymes
      integer                                  :: i,imax,id,ic,imark,imarkold,imin,nvar,ivar,index,j,ier,ier2
      logical                                  :: varfound,readerror

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
            ! remove tab characters
            fullline = strippedline(fullline)
            ! Check params.txt has old (unsupported) method of defining points
            ic=scan(fullline,'#')
            if (ic .ne. 0) then ! This branch of the if/else will change to a halt_program statement in later versions
               ! (written on 13 January 2011)
               call writelog('lswe','','Error in definition of point output.')
               call writelog('lswe','','Use of #var1#var2#etc. is no longer valid')
               call writelog('lswe','','Stopping simulation')
               call halt_program
            else ! not old method of setting point output
               read(fullline,*,iostat=ier)xpoints(i+imin),ypoints(i+imin),par%stationid(i+imin)
               ! error checking
               if(ier==0) then
                  ! all fine
                  readerror=.false.
               elseif(ier==-1) then
                  ! line is too short, probably missing name of station
                  ! try reading just the coordinates
                  read(fullline,*,iostat=ier2)xpoints(i+imin),ypoints(i+imin)
                  if (ier2==0) then
                     ! coordinates okay, just name missing
                     select case (readtype)
                      case('point')
                        write(par%stationid(i+imin),'("point",i0.3)') i
                      case('rugauge')
                        write(par%stationid(i+imin),'("rugau",i0.3)') i
                     end select
                     readerror=.false.
                  else
                     readerror=.true.
                  endif
               else
                  readerror=.true.
               endif
               if(readerror) then
                  ! Error reading point/rugauge input. Stop
                  select case (readtype)
                   case('point')
                     call writelog('lswe','','Error reading output point location/name in the following line in params.txt:')
                   case('rugauge')
                     call writelog('lswe','','Error reading runup gauge location/name in the following line in params.txt:')
                  end select
                  call writelog('lswe','',trim(fullline))
                  call writelog('lswe','','Stopping simulation')
                  call halt_program
               else
                  call writelog('ls','(a,a,a,f0.2,a,f0.2)',&
                  trim(okaymes),trim(par%stationid(i+imin)),' xpoint: ',&
                  xpoints(i+imin),'   ypoint: ',ypoints(i+imin))
               endif
            endif ! old method of point output
         enddo
         close(10)
      endif

   end subroutine readPointPosition
   
   subroutine check_instat_backward_compatibility(par)
      use logging_module
      use mnemmodule
      use readkey_module
      
      implicit none
      
      type(parameters), intent(inout)          :: par
      
      if (isSetParameter('params.txt','nonh') .or. isSetParameter('params.txt','instat')) then
            call writelog('l','','--------------------------------')
            call writelog('l','','Backward compatibility:')
      endif
      !
      ! Part 1: nonh still needs to work
      if (isSetParameter('params.txt','nonh') .or. par%useXBeachGSettings==1) then
            call writelog('ws','(a,a,a)','Warning: Specification of nonh using parameter ''wavemodel''')
            if (par%useXBeachGSettings==1) then
               par%nonh      = readkey_int ('params.txt','nonh',        1,        0,     1)
            else
               par%nonh      = readkey_int ('params.txt','nonh',        0,        0,     1)
            endif
            if (par%nonh==1) then
                par%wavemodel = WAVEMODEL_NONH
                par%wavemodel_str = 'nonh'
            endif 
      else
      endif
      !
      ! Part 2:instat still needs to work and defines both wavemodel and wbctypef
      if (isSetParameter('params.txt','instat')) then
            !
            ! Write warning message
            call writelog('ws','(a,a,a)','Warning: Specification of instat using parameter ''wbctype''')
            !
            ! Read instat anyway
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
            call parmapply('instat',2,par%instat,par%instat_str)
            !
            ! Conditions without spectra (e.g. stat, ts_1, etc.)
            ! 1) INSTAT = STAT
            if (par%instat==INSTAT_STAT) then
                if (par%nonh==1) then
                   par%wavemodel       = WAVEMODEL_NONH     
                   par%wavemodel_str   = 'nonh' 
                else
                   par%wavemodel       = WAVEMODEL_STATIONARY     
                   par%wavemodel_str   = 'stationary' 
                endif
                par%wbctype         = WBCTYPE_PARAMS    
                par%wbctype_str    = 'params'
            endif
            ! 2) INSTAT = TS_1
            if (par%instat==INSTAT_TS_1) then
                par%wavemodel       = WAVEMODEL_SURFBEAT     
                par%wavemodel_str   = 'surfbeat' 
                par%wbctype         = WBCTYPE_TS_1    
                par%wbctype_str    = 'ts_1'
            endif
            ! 3) INSTAT = TS_2            
            if (par%instat==INSTAT_TS_2) then
                par%wavemodel       = WAVEMODEL_SURFBEAT     
                par%wavemodel_str   = 'surfbeat' 
                par%wbctype         = WBCTYPE_TS_2    
                par%wbctype_str    = 'ts_2'
            endif
            ! 4) INSTAT = OFF
            if (par%instat==INSTAT_OFF .and. par%wavemodel == WAVEMODEL_NONH)   then  
                par%wavemodel       = WAVEMODEL_NONH    
                par%wavemodel_str   = 'nonh'       
                par%wbctype         = WBCTYPE_OFF    
                par%wbctype_str    = 'off'
            elseif (par%instat==INSTAT_OFF) then
                par%wavemodel       = WAVEMODEL_SURFBEAT    
                par%wavemodel_str   = 'surfbeat'   
                par%wbctype         = WBCTYPE_OFF    
                par%wbctype_str    = 'off'
            endif
            ! 5) INSTAT = BICHROM
            if (par%instat==INSTAT_BICHROM)   then  
                par%wavemodel       = WAVEMODEL_SURFBEAT    
                par%wavemodel_str   = 'surfbeat'   
                par%wbctype         = WBCTYPE_PARAMS    
                par%wbctype_str    = 'params'
                ! But we should define the Tlong
                if (isSetParameter('params.txt','Tlong')) then 
                    ! we will read it in later
                else
                    par%Tlong = 80
                endif
            endif
            ! 6) INSTAT = STAT_TABLE
            if (par%instat==INSTAT_STAT_TABLE)   then  
                par%wavemodel       = WAVEMODEL_STATIONARY
                par%wavemodel_str   = 'stationary'   
                par%wbctype         = WBCTYPE_JONS_TABLE    
                par%wbctype_str    = 'jons_table'
            endif
            ! 7) INSTAT = ts_nonh
            if (par%instat==INSTAT_TS_NONH)   then  
                par%wavemodel       = WAVEMODEL_NONH
                par%wavemodel_str   = 'nonh'   
                par%wbctype         = WBCTYPE_TS_NONH    
                par%wbctype_str    = 'ts_nonh'
            endif
            !
            ! Conditions with spectra
            ! 1) INSTAT_JONS
            if (par%instat==INSTAT_JONS .and. par%wavemodel == WAVEMODEL_NONH)   then  
                par%wavemodel       = WAVEMODEL_NONH    
                par%wavemodel_str   = 'nonh'       
                par%wbctype         = WBCTYPE_PARAMETRIC    
                par%wbctype_str    = 'parametric'
            elseif (par%instat==INSTAT_JONS) then
                par%wavemodel       = WAVEMODEL_SURFBEAT    
                par%wavemodel_str   = 'surfbeat'   
                par%wbctype         = WBCTYPE_PARAMETRIC    
                par%wbctype_str    = 'parametric'
            endif
            !
            ! 2) INSTAT_JONS_TABLE
            if (par%instat==INSTAT_JONS_TABLE .and. par%wavemodel == WAVEMODEL_NONH)   then  
                par%wavemodel       = WAVEMODEL_NONH    
                par%wavemodel_str   = 'nonh'
                par%wbctype         = WBCTYPE_JONS_TABLE    
                par%wbctype_str    = 'jons_table'
            elseif (par%instat==INSTAT_JONS_TABLE) then
                par%wavemodel       = WAVEMODEL_SURFBEAT    
                par%wavemodel_str   = 'surfbeat'   
                par%wbctype         = WBCTYPE_JONS_TABLE    
                par%wbctype_str    = 'jons_table'
            endif
            !
            ! 3) INSTAT = SWAN
            if (par%instat==INSTAT_SWAN .and. par%wavemodel == WAVEMODEL_NONH)   then  
                par%wavemodel       = WAVEMODEL_NONH    
                par%wavemodel_str   = 'nonh'
                par%wbctype         = WBCTYPE_SWAN    
                par%wbctype_str    = 'swan'
            elseif (par%instat==INSTAT_SWAN) then
                par%wavemodel       = WAVEMODEL_SURFBEAT    
                par%wavemodel_str   = 'surfbeat'   
                par%wbctype         = WBCTYPE_SWAN    
                par%wbctype_str    = 'swan'
            endif
            !
            ! 4) INSTAT = VARDENS
            if (par%instat==INSTAT_VARDENS .and. par%wavemodel == WAVEMODEL_NONH)   then  
                par%wavemodel       = WAVEMODEL_NONH    
                par%wavemodel_str   = 'nonh'
                par%wbctype         = WBCTYPE_VARDENS    
                par%wbctype_str    = 'vardens'
            elseif (par%instat==INSTAT_VARDENS) then
                par%wavemodel       = WAVEMODEL_SURFBEAT    
                par%wavemodel_str   = 'surfbeat'   
                par%wbctype         = WBCTYPE_VARDENS    
                par%wbctype_str    = 'vardens'
            endif
            ! 
            ! Other: re-use
            if (par%instat==INSTAT_REUSE .and. par%wavemodel == WAVEMODEL_NONH)   then  
                par%wavemodel       = WAVEMODEL_NONH    
                par%wavemodel_str   = 'nonh'
                par%wbctype         = WBCTYPE_REUSE    
                par%wbctype_str    = 'reuse'
            elseif (par%instat==INSTAT_REUSE) then
                par%wavemodel       = WAVEMODEL_SURFBEAT    
                par%wavemodel_str   = 'surfbeat'   
                par%wbctype         = WBCTYPE_REUSE    
                par%wbctype_str    = 'reuse'
            endif
      endif
   end subroutine

end module params
