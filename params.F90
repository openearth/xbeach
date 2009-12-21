module params

type parameters
! These parameters are constants, variables read from params.txt, or are scalars derived directly from read input
!

   ! Grid parameters                                                                                                               
!  Type             name                   initialize    !  [unit] description                                                                                                                 
   character(256):: depfile                    = 'abc'   !  [-] Name of the input bathymetry file
   real*8        :: posdwn                     = -123    !  [-] Bathymetry is specified positive down (1) or positive up (-1)
   integer*4     :: nx                         = -123    !  [-] Number of computiation cell corners in x-direction
   integer*4     :: ny                         = -123    !  [-] Number of computiation cell corners in y-direction
   real*8        :: alfa                       = -123    !  [deg] Angle of x-axis from East
   integer*4     :: vardx                      = -123    !  [-] Switch for variable grid spacing: 1 = irregular spacing, 0 = regular grid spacing
   real*8        :: dx                         = -123    !  [m] Regular grid spacing in x-direction
   real*8        :: dy                         = -123    !  [m] Regular grid spacing in y-direction
   character(256) :: xfile                      = 'abc'   !  [name] Name of the file containing x-coordinates of the calculation grid
   character(256) :: yfile                      = 'abc'   !  [name] Name of the file containing y-coordinates of the calculation grid
   real*8        :: xori                       = -123    !  [m] X-coordinate of origin of axis
   real*8        :: yori                       = -123    !  [m] Y-coordinate of origin of axis
   real*8        :: thetamin                   = -123    !  [deg] Lower directional limit (angle w.r.t computational x-axis)
   real*8        :: thetamax                   = -123    !  [deg] Higher directional limit (angle w.r.t computational x-axis)
   real*8        :: dtheta                     = -123    !  [deg] Directional resolution
   integer*4     :: thetanaut                  = -123    !  [-] Thetamin,thetamax in cartesian (0) or nautical convention (1)
 
   ! Model time                                                                                                                    
   real*8        :: tstop                      = -123    !  [s] Stop time of simulation, in morphological time
   real*8        :: CFL                        = -123    !  [-] Maximum Courant-Friedrichs-Lewy number
    
   ! Physical constants                                                                                                            
   real*8        :: g                          = -123    !  [ms^-2] Gravitational acceleration
   real*8        :: rho                        = -123    !  [kgm^-3] Density of water
 
   ! Physical processes                                                                                                            
   integer*4     :: swave                      = -123    !  [-] Include short waves (1), exclude short waves (0)
   integer*4     :: lwave                      = -123    !  [-] Include short wave forcing on NLSW equations and boundary conditions (1), or exclude (0)
   integer*4     :: sedtrans                   = -123    !  [-] Include sediment transport (1) or exclude (0)
   integer*4     :: morphology                 = -123    !  [-] Include morphology (1) or exclude (0)
   integer*4     :: nonh                       = -123    !  [-] Non-hydrostatic pressure option: 0 = NSWE, 1 = NSW + non-hydrostatic pressure compensation Stelling & Zijlema, 2003 
   integer*4     :: gwflow                     = -123    !  [-] Turn on (1) or off (0) groundwater flow module
   integer*4     :: q3d                        = -123    !  [-] Turn on (1) or off (0) quasi-3D sediment transport module

   ! Wave numerics parameters                                                                                                      
   integer*4     :: scheme                     = -123    !  [-] Use first-order upwind (1) or Lax-Wendroff (2) for wave propagation
   real*8        :: wavint                     = -123    !  [s] Interval between wave module calls (only in stationary wave mode)
   real*8        :: maxerror                   = -123    !  [m] Maximum wave height error in wave stationary iteration
   integer*4     :: maxiter                    = -123    !  [-] Maximum number of iterations in wave stationary
 
   ! Flow numerics parameters                                                                                                      
   real*8        :: eps                        = -123    !  [m] Threshold water depth above which cells are considered wet
   real*8        :: hmin                       = -123    !  [m] Threshold water depth above which Stokes drift is included
   real*8        :: umin                       = -123    !  [m/s] Threshold velocity for upwind velocity detection and for vmag2 in eq. sediment concentration

   ! Sediment transport numerics parameters                                                                                  
   real*8        :: thetanum                   = -123    !  [-] Coefficient determining whether upwind (1) or central scheme (0.5) is used.
   integer*4     :: sourcesink                 = -123    !  [-] In suspended transport use source-sink terms to calculate bed level change (1) or sus transport gradients (0)
    
   ! Wave-current interaction parameters                                                                                           
   integer*4     :: wci                        = -123    !  [-] Turns on (1) or off (0) wave-current interaction
   real*8        :: hwci                       = -123    !  [m] Minimum depth until which wave-current interaction is used
   real*8        :: cats                       = -123    !  [Trep] Current averaging time scale for wci, in terms of mean wave periods
 
   ! Wave boundary condition parameters                                                                                            
   integer*4     :: instat                     = -123    !  [-] Wave boundary condtion type
   real*8        :: taper                      = -123    !  [s] Spin-up time of wave boundary conditions, in morphological time  
   real*8        :: Hrms                       = -123    !  [-] Hrms wave height for instat = 0,1,2,3
   real*8        :: Tm01                       = -123    !  [s] Old name for Trep
   real*8        :: Trep                       = -123    !  [s] Representative wave period for instat = 0,1,2,3
   real*8        :: Tlong                      = -123    !  [s] Wave group period for case instat = 1
   real*8        :: dir0                       = -123    !  [deg] Mean wave direction (Nautical convention) for instat = 0,1,2,3
   integer*4     :: m                          = -123    !  [-] Power in cos^m directional distribution for instat = 0,1,2,3
   character(24) :: lateralwave                = 'neumann'   !  [name] Switch for lateral boundary at left, 'neumann' = E Neumann, 'wavefront' = along wave front
   integer*4     :: leftwave                   = -123    !  [-] old name for lateralwave
   integer*4     :: rightwave                  = -123    !  [-] old name for lateralwave

   ! Wave-spectrum boundary condition parameters
   character(256):: bcfile                     = 'abc'   !  Note, will replace current lookup in boundary conditions [name] Name of spectrum file  
   integer*4     :: random                     = -123    !  [-] Random seed on (1) or off (0) for instat = 4,5,6 boundary conditions
   real*8        :: fcutoff                    = -123    !  [Hz] Low-freq cutoff frequency for instat = 4,5,6 boundary conditions
   integer*4     :: nspr                       = -123    !  [-] nspr = 1 bin all wave components for generation of qin (instat 4,5,6) in one direction, nspr = 0 regular long wave spreadin
   real*8        :: trepfac                    = -123    !  [-] Compute mean wave period over energy band: par%trepfac*maxval(Sf) for instat 4,5,6; converges to Tm01 for trepfac = 0.0 and
   real*8        :: sprdthr                    = -123    !  [-] Threshold ratio to maxval of S above which spec dens are read in (default 0.08*maxval)
   integer*4     :: oldwbc                     = -123    !  [-] (1) Use old version wave boundary conditions for instat 4,5,6
   real*8        :: rt                         = -123    !  [s] Duration of wave spectrum at offshore boundary, in morphological time 
   real*8        :: dtbc                       = -123    !  [s] Timestep used to describe time series of wave energy and long wave flux at offshore boundary (not affected by morfac)
   real*8        :: dthetaS_XB                 = -123    !  [deg] If SWAN input is not in nautical degrees, dthetaS_XB is the angle from SWAN x-axis to XBeach x-axis in cathesian degrees
      
   ! Flow boundary condition parameters
   integer*4     :: front                      = -123    !  [-] Switch for seaward flow boundary: 0 = radiating boundary(Ad), 1 = Van Dongeren, 1997
   integer*4     :: left                       = -123    !  [name] Switch for lateral boundary at ny+1, 'neumann' = vv computed from NSWE, 'wall' = reflective wall; vv=0
   integer*4     :: right                      = -123    !  [-] Switch for lateral boundary at right, 0 = vv computed from NSWE, 1 = reflective wall; vv=0 
   integer*4     :: back                       = -123    !  [-] Switch for boundary at bay side, 0 = radiating boundary (Ad), 1 = reflective boundary; uu=0
   integer*4     :: ARC                        = -123    !  [-] Switch for active reflection compensation at seaward boundary: 0 = reflective, 1 = weakly (non) reflective
   real*4        :: order                      = -123    !  [-] Switch for order of wave steering, 1 = first order wave steering (short wave energy only), 2 = second oder wave steering (bound long wave corresponding to short wave forcing is added)
   integer*4     :: carspan                    = -123    !  [-] Switch for Carrier-Greenspan test 0 = use cg (default); 1 = use sqrt(gh) in instat = 3 for c&g tests
   character(256):: zsinitfile                 = 'abc'   !  [name] Name of inital condition file zs

   ! Tide boundary conditions                                                                                                      
   real*8        :: zs0                        = -123    !  [m] Inital water level
   real*8        :: zs01                       = -123    !  [m] Initial water level first sea boundary
   real*8        :: zs02                       = -123    !  [m] Initial water level second sea boundary
   real*8        :: zs03                       = -123    !  [m] Initial water level first land boundary
   real*8        :: zs04                       = -123    !  [m] Initial water level second land boundary
   character(256):: zs0file                    = 'abc'   !  Note, will replace lookup in readtide [name] Name of tide boundary condition series
   integer*4     :: tideloc                    = -123    !  [-] Number of corner points on which a tide time series is specified
   integer*4     :: paulrevere                 = -123    !  [-] Specifies tide on sea and land (0) or two sea points (1) if tideloc = 2
                                                         !      if tideloc =>2, then this indicates where the time series are to be 
                                                         !      applied. Input for tidal information to xbeach options (3):
                                                         !      1. one tidal record --> specify tidal record everywhere
                                                         !      2. two tidal records --> Need to specify keyword 'paulrevere'
                                                         !      paulrevere==0 implies to apply one tidal record to
                                                         !      both sea corners and one tidal record to both land corners
                                                         !      paulrevere==1 implies to apply the first tidal record
                                                         !      (column 2 in zs0input.dat) to the (x=1,y=1) sea corner and
                                                         !      the second tidal record (third column) to the (x=1,y=N) sea corner
                                                         !      3. four tidal records --> Need to list tidal records in  
                                                         !      zs0input.dat in order of:
                                                         !         (x=1,y=1)
                                                         !         (x=1,y=N)
                                                         !         (x=N,y=N)
                                                         !         (x=N,y=1)
                                                         !      NOTE:  clockwise from (1,1) corner

   real*8,dimension(2) :: xyzs01               = -123    ! global xy coordinates of corner (x=1,y=1) 
   real*8,dimension(2) :: xyzs02               = -123    ! global xy coordinates of corner (x=1,y=N)
   real*8,dimension(2) :: xyzs03               = -123    ! global xy coordinates of corner (x=N,y=N)
   real*8,dimension(2) :: xyzs04               = -123    ! global xy coordinates of corner (x=N,y=1)

   ! Discharge boundary conditions 
   character(256):: disch_loc_file             = 'abc'   !  Note: will replace lookup in boundary conditions [name] Name of discharge locations file
   character(256):: disch_timeseries_file      = 'abc'   !  Note: will replace lookup in boundary conditions [name] Name of discharge timeseries file

   ! Wave breaking parameters                                                                                                      
   integer*4     :: break                      = -123    !  [-] Type of breaker formulation (1=roelvink, 2=baldock, 3=roelvink adapted, 4=roelvink on/off breaking)
   real*8        :: gamma                      = -123    !  [-] Breaker parameter in Baldock or Roelvink formulation
   real*8        :: gamma2                     = -123    !  [-] End of breaking parameter in break = 4 formulation
   real*8        :: alpha                      = -123    !  [-] Wave dissipation coefficient in Roelvink formulation
   real*8        :: n                          = -123    !  [-] Power in Roelvink dissipation model
   real*8        :: gammax                     = -123    !  [-] Maximum ratio wave height to water depth
   real*8        :: delta                      = -123    !  [-] Fraction of wave height to add to water depth
   real*8        :: waverr                     = -123    !  [-] max. absolute wave height difference between time steps
   real*8        :: fw                         = -123    !  [-] Bed friction factor
 
   ! Roller parameters                                                                                                             
   integer*4     :: roller                     = -123    !  [-] Turn on (1) or off(0) roller model
   real*8        :: beta                       = -123    !  [-] Breaker slope coefficient in roller model
   integer*4     :: rfb                        = -123    !  [-] If rfb = 1 then maximum wave surface slope is feeded back in roller energy balance; else rfb = par%Beta
 
   ! Flow parameters                                                                                                               
   real*8        :: C                          = -123    !  [m^0.5s^-1] Chezy coefficient
   real*8        :: cf                         = -123    !  [-] Friction coefficient flow
   real*8        :: nuh                        = -123    !  [m^2s^-1] Horizontal background viscosity 
   real*8        :: nuhfac                     = -123    !  [-] Viscosity switch for roller induced turbulent horizontal viscosity
   real*8        :: nuhv                       = -123    !  [-] Longshore viscosity enhancement factor, following Svendsen (?)
   integer*4     :: smag                       = -123    !  [-] (=1) Use smagorinsky subgrid model for viscocity
   
   ! Coriolis force parameters
   real*8        :: wearth                     = -123    !  [hour^-1] Angular velocity of earth calculated as: 1/rotation_time (in hours), later changed in calculation code to rad/s
   real*8        :: lat                        = -123    !  [deg] Latitude at model location  for computing Coriolis
   
   ! Wind parameters
   real*8        :: rhoa                       = -123    !  [kgm^-3] Air density
   real*8        :: Cd                         = -123    !  [-] Wind drag coefficient
   real*8        :: windv                      = -123    !  [ms^-1] Wind velocity, in case of stationary wind
   real*8        :: windth                     = -123    !  [deg] Nautical wind direction, in case of stationary wind
   character(256) :: windfile                   = 'abc'   !  Note, will replace lookup in readwind [name] Name of file with non-stationary wind data
   
   ! Groundwater parameters 
   real*8        :: kx                         = -123    !  [ms^-1] Darcy-flow permeability coefficient in x-direction [m/s]
   real*8        :: ky                         = -123    !  [ms^-1] Darcy-flow permeability coefficient in y-direction [m/s]
   real*8        :: kz                         = -123    !  [ms^-1] Darcy-flow permeability coefficient in z-direction [m/s]
   real*8        :: dwetlayer                  = -123    !  [m] Thickness of the top soil layer interacting more freely with the surface water
   real*8        :: aquiferbot                 = -123    !  Note, will replace lookup in groundwater module [m] Level of uniform aquifer bottom
   character(256):: aquiferbotfile             = 'abc'   !  Note, will replace lookup in groundwater module [name] Name of the aquifer bottom file
   real*8        :: gw0                        = -123    !  Note, will replace lookup in groundwater module [m] Level initial groundwater level
   character(256):: gw0file                    = 'abc'   !  Note, will replace lookup in groundwater module [name] Name of initial groundwater level file

   
   ! Q3D sediment transport parameters
   real*8        :: vonkar                     = -123    !   von Karman constant
   real*8        :: vicmol                     = -123    !   molecular viscosity
   integer*4     :: kmax                       = -123    !  [-] Number of sigma layers in Quasi-3D model; kmax = 1 (default) is without vertical structure of flow and suspensions
   real*8        :: sigfac                     = -123    !  [-] dsig scales with log(sigfac). Default = 1.3
   
   ! Non-hydrostatic correction parameters
   integer*4     :: secorder                   = -123    ! [-] Use second order corrections to advection/non-linear terms based on mcCormack scheme
   integer*4     :: solver_maxit               = -123    ! [-] Maximum number of iterations in the linear SIP solver
   real*8        :: solver_acc                 = -123    ! [?] accuracy with respect to the right-hand side used
                                                         !     in the following termination criterion:
                                                         !         ||b-Ax || < acc*||b||
   real*8        :: solver_alpha               = -123    ! [?] Underrelaxation parameter 
   integer*4     :: solver                     = -123    ! [-] Solver used to solve the linear system, 1=SIP, 2=TRIDIAG (only for 1d)
   real*8        :: kdmin                      = -123    ! Minimum value of kd ( pi/dx > minkd )
   real*8        :: dispc                      = -123    ! Coefficient in front of the vertical pressure gradient, Default = 1.
   real*8        :: Topt                       = -123    ! Absolute period to optimize coefficient
   
   ! Bed composition parameters
   real*8        :: rhos                       = -123    !  [kgm^-3] Solid sediment density (no pores)
   integer*4     :: ngd                        = -123    !  [-] Number of sediment classes
   integer*4     :: nd                         = -123    !  [-] Number of computational layers in the bed
   real*8        :: dzg1                       = -123    !  [m] Thickness of top sediment class layers
   real*8        :: dzg2                       = -123    !  [m] Nominal thickness of variable sediment class layer
   real*8        :: dzg3                       = -123    !  [m] Thickness of bottom sediment class layers
   real*8        :: por                        = -123    !  [-] Porosity
   !real*8,dimension(:),allocatable  :: D50               !  Note, will replace lookup in initialize [m] D50 grain size per grain type
   !real*8,dimension(:),allocatable  :: D90               !  Note, will replace lookup in initialize [m] D90 grain size per grain type
   !real*8,dimension(:),allocatable  :: sedcal            !  Note, will replace lookup in initialize [-] Sediment transport calibration coefficient per grain type
   !real*8,dimension(:),allocatable  :: ucrcal            !  Note, will replace lookup in initialize [-] Critical velocity calibration coefficient per grain type

   ! Bed update numerics parameters
   real*8        :: frac_dz                    = -123    !  [-] Relative thickness to split time step for bed updating
   integer*4     :: nd_var                     = -123    !  [-] Index of layer with variable thickness 
   real*8        :: split                      = -123    !  [-] Split threshold for variable sediment layer (ratio to nominal thickness)
   real*8        :: merge                      = -123    !  [-] Merge threshold for variable sediment layer (ratio to nominal thickness)


   ! Sediment transport parameters                                                                                                 
   integer*4     :: form                       = -123    !  [-] Equilibrium sed. conc. formulation: 1 = Soulsby van Rijn, 1997, 2 = Van Rijn 2008 with modifications by Van Thiel
   integer*4     :: sws                        = -123    !  [-] 1 = short wave & roller stirring and undertow, 0 = no short wave & roller stirring and undertow
   integer*4     :: lws                        = -123    !  [-] 1 = long wave stirring, 0 = no long wave stirring
   real*8        :: BRfac                      = -123    !  [-] Calibration factor surface slope
   real*8        :: facsl                      = -123    !  [-] Factor bedslope effect
   real*8        :: z0                         = -123    !  [m] Zero flow velocity level in Soulsby van Rijn (1997) sed.conc. expression
   real*8        :: smax                       = -123    !  [-] Being tested: maximum Shields parameter for ceq Diane Foster
   real*8        :: tsfac                      = -123    !  [-] Coefficient determining Ts = tsfac * h/ws in sediment source term
   real*8        :: facua                      = -123    !  [-] Calibration factor time averaged flows due to wave asymmetry
   integer*4     :: turb                       = -123    !  [-] Equlibrium sediment concentration is computed as function of:
   real*8        :: Tbfac                      = -123    !  [-] Calibration factor for bore interval Tbore: Tbore = Tbfac*Tbore
   real*8        :: Tsmin                      = -123    !  [s] Minimum adaptation time scale in advection diffusion equation sediment                                                      !      0 = no turbulence, 1 = wave averaged turbulence, 2 = maximum turbulence
   integer*4     :: lwt                        = -123    !  [-] Switch 0/1 long wave turbulence model
   real*8        :: betad                      = -123    !  [-] Dissipation parameter long wave breaking turbulence

   ! Morphology parameters                                                                                                         
   real*8        :: morfac                     = -123    !  [-] Morphological acceleration factor
   integer*4     :: morfacopt                  = -123    !  [-] Option indicating whether times should be adjusted (1) or not(0) for morfac
   real*8        :: morstart                   = -123    !  [s] Start time morphology, in morphological time
   real*8        :: wetslp                     = -123    !  [-] Critical avalanching slope under water (dz/dx and dz/dy)
   real*8        :: dryslp                     = -123    !  [-] Critical avalanching slope above water (dz/dx and dz/dy)
   real*8        :: hswitch                    = -123    !  [m] Water depth at which is switched from wetslp to dryslp
   real*8        :: dzmax                      = -123    !  [m/s/m] Maximum bedlevel change due to avalanching
   integer*4     :: struct                     = -123    !  [-] Switch for hard structures
   character(256):: ne_layer                   = 'abc'   !  [name] Name of file containing depth of hard structure
 
   ! Output variables                                                                                                              
   integer*4     :: timings                    = -123    !  [-] Switch to turn on (1) or off (0) progress output to screen
   real*8        :: tstart                     = -123    !  [s] Start time of output, in morphological time
   real*8        :: tint                       = -123    !  [s] Interval time of global output (replaced by tintg)
   real*8        :: tintg                      = -123    !  [s] Interval time of global output
   real*8        :: tintp                      = -123    !  [s] Interval time of point and runup gauge output
   real*8        :: tintc                      = -123    !  [s] Interval time of cross section output
   real*8        :: tintm                      = -123    !  [s] Interval time of mean,var,max,min output
   character(256):: tsglobal                   = 'abc'   !  [name] Name of file containing timings of global output
   character(256):: tspoints                   = 'abc'   !  [name] Name of file containing timings of point output
   character(256):: tscross                    = 'abc'   !  [name] Name of file containing timings of cross section output
   character(256):: tsmean                     = 'abc'   !  [name] Name of file containing timings of mean, max, min and var output
   integer*4     :: nglobalvar                 = -123    !  [-] Number of global output variables
   integer*4     :: nmeanvar                   = -123    !  [-] Number of mean,min,max,var output variables
   integer*4     :: npoints                    = -123    !  [-] Number of output point locations
   integer*4     :: nrugauge                   = -123    !  [-] Number of output runup gauge locations
   integer*4     :: ncross                     = -123    !  [-] Number of output cross sections
  
   ! Drifters parameters
   integer*4     :: ndrifter                   = -123    !  Note: will replace lookup in drifters module [-] Number of drifers
   character(256) :: drifterfile                = 'abc'   !  Note: will replace lookup in drifters module [name] Name of drifter data file

   ! MPI parameters
   character(4)   :: mpiboundary               = 'abc'   ! Fix mpi boundaries along y-lines ('y'), x-lines ('x'), or find shortest boundary ('auto')
  
   ! Constants, not read in params.txt
   real*8               :: px                  = -123    !  [-] Pi
   complex(kind(0.0d0)) :: compi               = -123    !  [-] Imaginary unit
   real*8               :: rhog8               = -123    !  [Nm^-3] 1/8*rho*g

   ! Variables, not read in params.txt
   real*8               :: dt                  = -123    !  [s] Computational time step, in hydrodynamic time
   real*8               :: t                   = -123    !  [s] Computational time, in hydrodynamic time
   real*8               :: tnext               = -123    !  [s] Next time point for output or wave stationary calculation, in hydrodynamic time
   real*8               :: Emean               = -123    !  Note: can be made local [Jm^-2] Mean wave energy at boundary
   real*8               :: w                   = -123    !  Note: can be made local [ms^-1] Fall velocity sediment
   integer*4            :: listline            = -123    !  Note: can be made local [-] Keeps rack of the record line in bcf-files in case of instat 4,5,6
   real*8               :: fc                  = -123    !  Note: can be made local [s^-1] Coriolis forcing parameter
   integer*4            :: tidelen             = -123    !  [-] Length of input tidal time series 
   integer*4            :: windlen             = -123    !  [-] Length of input wind time series 
   real*8               :: Llong               = -123    !  Note: can be made local [m] Alongshore wave group length for case instat=1


end type parameters

contains

subroutine all_input(par)

use readkey_module
use xmpi_module
use general_fileio
implicit none
type(parameters)            :: par

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!        Grid parameters       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
par%xori  = readkey_dbl('params.txt','xori',  0.d0,   -1d9,      1d9)
par%yori  = readkey_dbl('params.txt','yori',  0.d0,   -1d9,      1d9)
par%alfa  = readkey_dbl('params.txt','alfa',  0.d0,   -360.d0,   360.d0)
par%nx    = readkey_int('params.txt','nx',     50,      2,     10000)
par%ny    = readkey_int('params.txt','ny',      2,      2,     10000)
par%posdwn= readkey_dbl('params.txt','posdwn', 1.d0,     -1.d0,     1.d0)
call readkey('params.txt','depfile',par%depfile)  ! Bathymetry file name
call check_file_exist(par%depfile)
call check_file_length(par%depfile,par%nx+1,par%ny+1)
par%vardx = readkey_int('params.txt','vardx',   0,      0,         1) 
if (par%vardx==0) then
  par%dx    = readkey_dbl('params.txt','dx',    0.d0,   -1d9,      1d9)
  par%dy    = readkey_dbl('params.txt','dy',    0.d0,   -1d9,      1d9)
else
  call readkey('params.txt','xfile',par%xfile)    ! X-grid file
  call check_file_exist(par%xfile)
  call check_file_length(par%xfile,par%nx+1,par%ny+1)
  call readkey('params.txt','xfile',par%xfile)    ! Y-grid file
  call check_file_exist(par%xfile)
  call check_file_length(par%xfile,par%nx+1,par%ny+1)
endif
par%thetamin = readkey_dbl ('params.txt','thetamin', -80.d0,    -180.d0,  180.d0)
par%thetamax = readkey_dbl ('params.txt','thetamax',  80.d0,    -180.d0,  180.d0)
par%dtheta   = readkey_dbl ('params.txt','dtheta',    10.d0,      0.1d0,   20.d0)
par%thetanaut= readkey_int ('params.txt','thetanaut',    0,        0,     1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!        Model time parameters       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
par%CFL     = readkey_dbl ('params.txt','CFL',     0.2d0,     0.1d0,      0.9d0)
par%tstop   = readkey_dbl ('params.txt','tstop', 2000.d0,      1.d0, 1000000.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!        Physical constants       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
par%rho   = readkey_dbl ('params.txt','rho',  1025.0d0,  1000.0d0,  1040.0d0)
par%g     = readkey_dbl ('params.txt','g',      9.81d0,     9.7d0,     9.9d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!        Physical processes       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
par%swave      = readkey_int ('params.txt','swave',         1,        0,     1)
par%lwave      = readkey_int ('params.txt','lwave',         1,        0,     1)
par%sedtrans   = readkey_int ('params.txt','nonh',          0,        0,     1)
par%morphology = readkey_int ('params.txt','sedtrans',      1,        0,     1)
par%nonh       = readkey_int ('params.txt','morphology',    1,        0,     1)
par%gwflow     = readkey_int ('params.txt','gwflow',        0,        0,     1)
par%q3d        = readkey_int ('params.txt','q3d',           0,        0,     1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!        Post-input processing      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (par%alfa.lt.0) then 
   par%alfa = 360.d0+par%alfa
endif
par%alfa  = par%alfa*atan(1.0d0)/45.d0   ! All input converted directly to cathesian XBeach grid direction

if (par%posdwn<0.1d0) then 
   par%posdwn=-1.d0  ! Backward compatibility, now posdwn = 0 also works in input (i.e. posdwn = false)
endif

! All input time frames converted to XBeach hydrodynamic time
if (par%morfacopt==1) par%tstop   = par%tstop  / max(par%morfac,1.d0)


end subroutine all_input


subroutine wave_input(par)
use readkey_module
use xmpi_module
implicit none
type(parameters)            :: par

character(len=80)          :: dummystring

par%px    = 4.d0*atan(1.d0)
par%compi = (0.0d0,1.0d0)

par%instat   = readkey_int     ('params.txt','instat',    1,         0,        7)
par%fcutoff  = readkey_dbl     ('params.txt','fcutoff',   0.d0,      0.d0,     40.d0)
par%random   = readkey_int     ('params.txt','random',    0,         0,        1)
if (par%instat == 0) then
    par%dir0  = readkey_dbl    ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
    par%Hrms  = readkey_dbl    ('params.txt','Hrms',      1.d0,      0.d0,    10.d0)
    par%wavint   = readkey_dbl ('params.txt','wavint',    1.d0,      1.d0,  3600.d0)
    par%m     = readkey_int    ('params.txt','m',        10,         2,      128)
    par%Trep  = readkey_dbl    ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
    par%Trep  = readkey_dbl    ('params.txt','Trep',     par%Trep,   1.d0,    20.d0)
	 par%maxiter = readkey_int  ('params.txt','maxiter', 200,         2,      100)
    par%maxerror= readkey_dbl  ('params.txt','maxerror', 0.0001d0, 0.00001d0, 0.01d0)
!    par%omega    = 2.d0*par%px/par%Trep;
elseif (par%instat==40) then
    par%wavint   = readkey_dbl ('params.txt','wavint',    1.d0,      1.d0,  3600.d0)
	 par%maxiter = readkey_int  ('params.txt','maxiter',  15,         2,      100)
    par%maxerror= readkey_dbl  ('params.txt','maxerror', 0.0001d0, 0.00001d0, 0.01d0)
elseif (par%instat==1) then
    par%dir0  = readkey_dbl    ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
    par%Hrms  = readkey_dbl    ('params.txt','Hrms',      1.d0,      0.d0,    10.d0)
    par%Tlong = readkey_dbl    ('params.txt','Tlong',    80.d0,     20.d0,   300.d0)
    par%m     = readkey_int    ('params.txt','m',        10,         2,      128)
    par%Trep  = readkey_dbl    ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
    par%Trep  = readkey_dbl    ('params.txt','Trep',     par%Trep,   1.d0,    20.d0)
!    par%omega    = 2.d0*par%px/par%Trep;
elseif (par%instat==2 .or. par%instat==3) then
    par%dir0  = readkey_dbl    ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
    par%Hrms  = readkey_dbl    ('params.txt','Hrms',      1.d0,      0.d0,    10.d0)
    par%m     = readkey_int    ('params.txt','m',        10,         2,      128)
    par%Trep  = readkey_dbl    ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
    par%Trep  = readkey_dbl    ('params.txt','Trep',     par%Trep,   1.d0,    20.d0)
!    par%omega    = 2.d0*par%px/par%Trep;
elseif (par%instat==4 .or. par%instat==41 .or. par%instat==5 .or. par%instat==6) then
        ! Just a check .....
        if (xmaster) then
          call readkey('params.txt','bcfile',dummystring)
          call readkey('params.txt','rt',dummystring)
          call readkey('params.txt','dtbc',dummystring)
          call readkey('params.txt','dthetaS_XB',dummystring)
        endif
elseif (par%instat == 8) then        
    par%Topt  = readkey_dbl    ('params.txt','Topt',     10.d0,      1.d0,    20.d0)
    par%zs0   = readkey_dbl    ('params.txt','zs0',     0.0d0,      -100.0d0,    100.0d0)             		
elseif (par%instat > 8.and. (par%instat>41.or.par%instat<40)) then
    if(xmaster) then
      write(*,*)'Instat invalid option'
    endif
    call halt_program
end if
!                         Input file  Keyword Default  Minimum  Maximum
!par%dir0  = readkey_dbl    ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
par%hmin  = readkey_dbl ('params.txt','hmin',   0.01d0,   0.001d0,      1.d0)
par%gammax= readkey_dbl ('params.txt','gammax',   2.d0,      .4d0,      5.d0)    !changed 28/11
par%gamma = readkey_dbl ('params.txt','gamma',   0.55d0,     0.4d0,     0.9d0)   !changed 28/11
par%gamma2= readkey_dbl ('params.txt','gamma2',   0.3d0,     0.0d0,     0.5d0)   !added 16/3/09
par%alpha = readkey_dbl ('params.txt','alpha',   1.0d0,     0.5d0,     2.0d0)
par%delta = readkey_dbl ('params.txt','delta',   0.0d0,     0.0d0,     1.0d0)
par%n    =  readkey_dbl ('params.txt','n',       10.0d0,     5.0d0,    20.0d0)   !changed 28/11
par%rho   = readkey_dbl ('params.txt','rho',  1025.0d0,  1000.0d0,  1040.0d0)
par%g     = readkey_dbl ('params.txt','g',      9.81d0,     9.7d0,     9.9d0)
par%rhog8 = 1.0d0/8.0d0*par%rho*par%g
par%Emean = par%rhog8*par%Hrms**2
par%thetamin = readkey_dbl ('params.txt','thetamin', -80.d0,    -180.d0,  180.d0)
par%thetamax = readkey_dbl ('params.txt','thetamax',  80.d0,    -180.d0,  180.d0)
par%dtheta   = readkey_dbl ('params.txt','dtheta',    10.d0,      0.1d0,   20.d0)
par%thetanaut= readkey_int ('params.txt','thetanaut',    0,        0,     1)
par%wci      = readkey_int ('params.txt','wci',        0,        0,     1)
par%hwci     = readkey_dbl ('params.txt','hwci',   0.1d0,   0.001d0,      1.d0)
par%cats     = readkey_dbl ('params.txt','cats',  7.d0,   1.d0,      50.d0)
par%fw       = readkey_dbl ('params.txt','fw',  0.d0,   0d0,      0.3d0)
par%break    = readkey_int ('params.txt','break',      3,        1,     3)
! Only allow Baldock in stationary mode and Roelvink in non-stationary
if (par%instat==0) then
   if (par%break .ne. 2) then
        write(*,*)'Error: Roelvink formulation not allowed in stationary calculation, use Baldock formulation. Stopping.'
       	call halt_program
   endif
else
   if (par%break==2) then 
        write(*,*)'Warning: Baldock formulation not allowed in non-stationary calculation, use Roelvink formulation. Stopping.'
   endif
endif
par%roller   = readkey_int ('params.txt','roller',     1,        0,     1)
par%beta     = readkey_dbl ('params.txt','beta',    0.15d0,     0.05d0,   0.3d0)
par%rfb      = readkey_int ('params.txt','rfb',        1,        0,     1)
par%taper    = readkey_dbl ('params.txt','taper',   100.d0,      0.0d0, 1000.d0)
par%taper    = max(par%taper,1.d-6)
!par%refl     = readkey_int ('params.txt','refl',       0,        0,     1) 
par%nspr     = readkey_int ('params.txt','nspr',       0,        0,     1) 
par%scheme   = readkey_int ('params.txt','scheme',     2,        1,     2)
par%trepfac  = readkey_dbl  ('params.txt','trepfac', 0.8d0,       0.d0,  1.d0) 
par%sprdthr  = readkey_dbl ('params.txt','sprdthr', 0.08d0,      0.d0,  1.d0) 
par%lwave    = readkey_int ('params.txt','lwave',         1,        0,     1)
par%swave    = readkey_int ('params.txt','swave',         1,        0,     1)
par%sws      = readkey_int ('params.txt','sws',           1,        0,     1)
par%lws      = readkey_int ('params.txt','lws',           1,        0,     1)
!par%ut       = readkey_int ('params.txt','ut',            1,        0,     1)
par%BRfac    = readkey_dbl ('params.txt','BRfac',    1.0d0,       0.d0, 1.d0)
par%oldwbc   = readkey_int ('params.txt','oldwbc',       1,        0,     1)

if (xmaster) call readkey('params.txt','mpiboundary',par%mpiboundary)    ! X-grid file
if (par%mpiboundary==' ') par%mpiboundary='auto'   ! Default

end subroutine wave_input

subroutine flow_input(par)
use readkey_module
use xmpi_module
implicit none
type(parameters)            :: par


par%cf      = readkey_dbl ('params.txt','cf',      3.d-3,     0.d0,     0.1d0)
par%C       = readkey_dbl ('params.txt','C',       sqrt(par%g/par%cf),     20.d0,    100.d0)
par%cf      = par%g/par%C**2
par%eps     = readkey_dbl ('params.txt','eps',     0.1d0,   0.001d0,      1.d0)
par%umin    = readkey_dbl ('params.txt','umin',    0.1d0,   0.001d0,      5.d0)
par%zs01    = readkey_dbl ('params.txt','zs0',     0.0d0,     -5.d0,      5.d0)
par%tideloc = readkey_int ('params.txt','tideloc', 0,             0,      4)
par%paulrevere = readkey_int ('params.txt','paulrevere', 0,       0,      1)
par%tidelen = readkey_int ('params.txt','tidelen',       0,       0,      1000000)
par%tstart  = readkey_dbl ('params.txt','tstart',   1.d0,      0.d0,1000000.d0)
par%tstop   = readkey_dbl ('params.txt','tstop', 2000.d0,      1.d0,1000000.d0)
par%tint    = readkey_dbl ('params.txt','tint',     1.d0,     .01d0, 100000.d0)  ! Robert
par%tintg   = readkey_dbl ('params.txt','tintg', par%tint,     .01d0, 100000.d0)  ! Robert
par%tintp   = readkey_dbl ('params.txt','tintp',par%tintg,    .01d0, 100000.d0)  ! Robert
par%tintc   = readkey_dbl ('params.txt','tintc',par%tintg,    .01d0, 100000.d0)  ! Robert
par%tintm   = readkey_dbl ('params.txt','tintm',par%tintg,     1.d0, par%tstop)  ! Robert
par%tint    = min(par%tintg,par%tintp,par%tintm,par%tintc)                       ! Robert     
! adapt flow times to morfac if par%morfacopt=1
par%morfac   = readkey_dbl ('params.txt','morfac', 0.0d0,        0.d0,  1000.d0)
par%morfacopt= readkey_int ('params.txt','morfacopt', 1,        0,        1)
if (par%morfacopt==1) then
   par%tstart  = par%tstart / max(par%morfac,1.d0)
   par%tint    = par%tint   / max(par%morfac,1.d0)
   par%tintg   = par%tintg  / max(par%morfac,1.d0)
   par%tintp   = par%tintp  / max(par%morfac,1.d0)
   par%tintc   = par%tintc  / max(par%morfac,1.d0)
   par%tintm   = par%tintm  / max(par%morfac,1.d0)
   par%wavint  = par%wavint / max(par%morfac,1.d0)
   par%tstop   = par%tstop  / max(par%morfac,1.d0)
endif
if (min(par%tintg,par%tintm,par%tintp,par%tintc)<=0.d0) then
   write(*,*)'Error: Output time steps of 0 or less are not allowed.'
   call halt_program
endif
par%fcutoff  = readkey_dbl ('params.txt','fcutoff', 0.d0,      0.d0,  40.d0)
par%sprdthr = readkey_dbl ('params.txt','sprdthr', 0.08d0,      0.d0,  1.d0)
par%carspan = readkey_int    ('params.txt','carspan',0,         0,      1)
!par%rugauge = readkey_int    ('params.txt','rugauge',0,         0,      1)
!par%ntout=nint((par%tstop-par%tstart)/par%tint)+1
par%CFL     = readkey_dbl ('params.txt','CFL',     0.2d0,     0.1d0,      0.9d0)
par%front   = readkey_int ('params.txt','front',       1,         0,      1) 
par%ARC     = readkey_int ('params.txt','ARC',         1,         0,      1)
par%order   = readkey_dbl ('params.txt','order',       2.d0,         1.d0,      2.d0)
par%left    = readkey_int ('params.txt','left',        0,         0,      1)
par%right   = readkey_int ('params.txt','right',       0,         0,      1)
par%leftwave= readkey_int ('params.txt','leftwave',        0,         0,      1)
par%rightwave= readkey_int ('params.txt','rightwave',        0,         0,      1)
par%back    = readkey_int ('params.txt','back',        2,         0,      2)
par%nuh     = readkey_dbl ('params.txt','nuh',     0.5d0,     0.0d0,      1.0d0)
par%nuhfac  = readkey_dbl ('params.txt','nuhfac',      0.0d0,     0.0d0,  1.0d0)

par%nonh    = readkey_int ('params.txt','nonh',        0,         0,      1)
par%smag    = readkey_int ('params.txt','smag',        0,         0,      1)
par%nuhv    = readkey_dbl ('params.txt','nuhv',     1.d0,      1.d0,    20.d0)
par%lat     = readkey_dbl ('params.txt','lat',     0.d0,      0.d0,   90.d0)
par%wearth  = readkey_dbl ('params.txt','omega',   1.d0/24.d0, 0.d0,    1.d0)
par%vonkar  = readkey_dbl ('params.txt','vonkar',   0.4d0,     0.01d0,  1.d0)
par%vicmol  = readkey_dbl ('params.txt','vicmol',   0.000001d0,   0.d0,    0.001d0)

par%lat = par%lat*par%px/180.d0
par%wearth = par%px*par%wearth/1800.d0
par%fc = 2.d0*par%wearth*sin(par%lat)

! Ground water
par%gwflow  = readkey_int ('params.txt','gwflow',         0,           0,       1)
if (par%gwflow==1) then
   par%kx         = readkey_dbl ('params.txt','kx'        , 0.0001d0 , 0.00001d0, 0.01d0)
   par%ky         = readkey_dbl ('params.txt','ky'        , par%kx   , 0.00001d0, 0.01d0)
   par%kz         = readkey_dbl ('params.txt','kz'        , par%kx   , 0.00001d0, 0.01d0)
   par%dwetlayer  = readkey_dbl ('params.txt','dwetlayer' , 0.2d0    , 0.01d0     , 1.d0)
endif

! Initial condition water level from file
if (xmaster) then
   call readkey('params.txt','zsinitfile',par%zsinitfile)
endif

par%secorder     = readkey_int('params.txt','secorder' ,0,0,1)
par%solver_maxit = readkey_int('params.txt','solver_maxit' ,30,1,1000)
par%kdmin        = readkey_dbl('params.txt','kdmin' ,0.0d0,0.0d0,0.05d0)  
par%dispc        = readkey_dbl('params.txt','dispc' ,1.0d0,0.1d0,2.0d0)  
par%solver_acc   = readkey_dbl('params.txt','solver_acc' ,0.005d0,0.00001d0,0.1d0)  
par%solver_alpha = readkey_dbl('params.txt','solver_urelax' ,0.92d0,0.5d0,0.99d0)
par%solver       = readkey_int('params.txt','solver' ,1,0,2)



end subroutine flow_input

subroutine sed_input(par)
use readkey_module
use xmpi_module
use general_fileio
implicit none
type(parameters)            :: par
character(len=80)           :: dummystring


!par%dico     = readkey_dbl ('params.txt','dico',    1.d0,        0.d0,    10.d0)
par%kmax     = readkey_int ('params.txt','kmax ',      1,           1,        1000)
par%ngd      = readkey_int ('params.txt','ngd',        1,           1,        20)
par%nd       = readkey_int ('params.txt','nd ',        1,           1,        1000)
par%nd_var   = readkey_int ('params.txt','nd_var', min(par%nd,2),   1,        par%nd)
par%dzg1     = readkey_dbl ('params.txt','dzg',    0.1d0,      0.01d0,     1.d0)
par%dzg1     = readkey_dbl ('params.txt','dzg1', par%dzg1,     0.01d0,     1.d0)
par%dzg2     = readkey_dbl ('params.txt','dzg2', par%dzg1,     0.01d0,     1.d0)
par%dzg3     = readkey_dbl ('params.txt','dzg3', par%dzg1,     0.01d0,     1.d0)
par%rhos     = readkey_dbl ('params.txt','rhos',  2650d0,     2400.d0,  2800.d0)
par%morfac   = readkey_dbl ('params.txt','morfac', 0.0d0,        0.d0,  1000.d0)
par%morfacopt= readkey_int ('params.txt','morfacopt', 1,        0,        1)
par%morstart = readkey_dbl ('params.txt','morstart',120.d0,      0.d0, 10000.d0)
if (par%morfacopt==1) par%morstart = par%morstart / max(par%morfac,1.d0)
par%wetslp   = readkey_dbl ('params.txt','wetslp', 0.3d0,       0.1d0,     1.d0)
par%dryslp   = readkey_dbl ('params.txt','dryslp', 1.0d0,       0.1d0,     2.d0)
par%por      = readkey_dbl ('params.txt','por',    0.4d0,       0.3d0,    0.5d0)
par%hswitch  = readkey_dbl ('params.txt','hswitch',0.1d0,      0.01d0,    1.0d0)  
par%z0       = readkey_dbl ('params.txt','z0     ',0.006d0,    0.0001d0,   0.05d0)  
par%facsl    = readkey_dbl ('params.txt','facsl  ',0.00d0,    0.00d0,   1.6d0)  
!par%struct   = readkey_int ('params.txt','struct ',0,           0,            1) 
par%form     = readkey_int ('params.txt','form',   1,           1,            3)
par%smax     = readkey_dbl ('params.txt','smax',   -1.d0,    -1.d0,   3.d0)       !changed 28/11 and back 10/2
if (par%smax<0) par%smax=huge(0.d0)
par%thetanum = readkey_dbl ('params.txt','thetanum',   1.d0,    0.5d0,   1.d0) 
par%tsfac    = readkey_dbl ('params.txt','tsfac',   0.1d0,    0.01d0,   1.d0) 
par%facua    = readkey_dbl ('params.txt','facua  ',0.00d0,    0.00d0,   1.0d0) 
par%dzmax    = readkey_dbl ('params.txt','dzmax  ',0.05d0,    0.00d0,   1.0d0) 
par%turb     = readkey_int ('params.txt','turb',   2,           0,            2)
par%Tbfac    = readkey_dbl ('params.txt','Tbfac  ',1.0d0,     0.00d0,   1.0d0) 
par%Tsmin    = readkey_dbl ('params.txt','Tsmin  ',0.2d0,     0.01d0,   10.d0) 
!par%impact   = readkey_int ('params.txt','impact ',0,           0,            1)  
!par%CE       = readkey_dbl ('params.txt','CE     ',0.2d0,     0.00d0,   100.d0) 
par%betad    = readkey_dbl ('params.txt','betad  ',1.0d0,     0.00d0,   10.0d0) 
par%lwt      = readkey_int ('params.txt','lwt    ',0,           0,            1)
par%sigfac   = readkey_dbl ('params.txt','sigfac ',1.3d0,     0.00d0,   10.d0) 
par%sourcesink      = readkey_int ('params.txt','sourcesink    ',0,           0,            1)
if (par%morfac>1.d0) then
   if (par%sourcesink==1) then
       write(*,*)'Warning: Using source-sink terms for bed level change with morfac can lead to loss of sediment mass conservation.'
   endif
endif
par%struct   = readkey_int ('params.txt','struct ',0    ,      0,             1)
if (par%struct==1) then
   call readkey('params.txt','ne_layer',par%ne_layer)  ! Bathymetry file name
   call check_file_exist(par%ne_layer)
   call check_file_length(par%ne_layer,par%nx+1,par%ny+1)
endif


! Just a check to see if they are there.....
if (xmaster) then
 call readkey('params.txt','D50',dummystring)
 call readkey('params.txt','D90',dummystring)
 call readkey('params.txt','sedcal',dummystring)
 call readkey('params.txt','ucrcal',dummystring)
endif
end subroutine sed_input

#ifdef USEMPI
subroutine distribute_par(par)
use xmpi_module
implicit none
type(parameters)        :: par
integer                 :: parlen,w,ierror

! 
! distribute parameters 

! new method: distribute at once

!inquire(iolength=parlen) par
!inquire(iolength=w) 1.0d0

! convert parlen to number of bytes, assuming that 
! 1.0d0 takes 8 bytes

!parlen = (8/w)*parlen

!call MPI_Bcast(par,parlen,MPI_BYTE,xmpi_master,xmpi_comm,ierror)
call MPI_Bcast(par,sizeof(par),MPI_BYTE,xmpi_master,xmpi_comm,ierror)
return

! so, the following code is NOT used anymore. I left this here
! maybe method above does not work everywhere. wwvv

! For efficiency, this subroutine should use MPI_Pack and 
! MPI_Unpack. However, since this subroutine is only called a 
! few times, a more simple approach is used.
!

!call xmpi_bcast(par%px)
!call xmpi_bcast(par%Hrms)
!call xmpi_bcast(par%Trep)
!call xmpi_bcast(par%dir0)
!call xmpi_bcast(par%m)
!call xmpi_bcast(par%nt)
!call xmpi_bcast(par%hmin)
!call xmpi_bcast(par%gammax)
!call xmpi_bcast(par%Tlong)
!call xmpi_bcast(par%Llong)
!call xmpi_bcast(par%gamma)
!call xmpi_bcast(par%delta)
!call xmpi_bcast(par%rho)
!call xmpi_bcast(par%g)
!call xmpi_bcast(par%rhog8)
!!call xmpi_bcast(par%omega)
!call xmpi_bcast(par%thetamin)
!call xmpi_bcast(par%thetamax)
!call xmpi_bcast(par%dtheta)
!call xmpi_bcast(par%thetanaut)
!call xmpi_bcast(par%wci)
!call xmpi_bcast(par%hwci)
!call xmpi_bcast(par%dt)
!call xmpi_bcast(par%break)
!call xmpi_bcast(par%instat)
!call xmpi_bcast(par%wavint)
!call xmpi_bcast(par%alpha)
!call xmpi_bcast(par%n)
!call xmpi_bcast(par%roller)
!call xmpi_bcast(par%beta)
!call xmpi_bcast(par%taper)
!call xmpi_bcast(par%t)
!call xmpi_bcast(par%tnext)
!call xmpi_bcast(par%it)
!call xmpi_bcast(par%tstart)
!call xmpi_bcast(par%tint)
!call xmpi_bcast(par%tintp)
!call xmpi_bcast(par%tintg)
!call xmpi_bcast(par%tintm)
!call xmpi_bcast(par%tstop)
!call xmpi_bcast(par%ntout)
!call xmpi_bcast(par%C)
!call xmpi_bcast(par%cf)
!call xmpi_bcast(par%eps)
!call xmpi_bcast(par%umin)
!call xmpi_bcast(par%zs01)
!call xmpi_bcast(par%zs02)
!call xmpi_bcast(par%zs03)
!call xmpi_bcast(par%zs04)
!call xmpi_bcast(par%tideloc)
!call xmpi_bcast(par%paulrevere)
!call xmpi_bcast(par%tidelen)
!
!call xmpi_bcast(par%facsl)
!call xmpi_bcast(par%nuh)
!call xmpi_bcast(par%nuhfac)
!call xmpi_bcast(par%rhos)
!call xmpi_bcast(par%morfac)
!call xmpi_bcast(par%morstart)
!call xmpi_bcast(par%Emean)
!call xmpi_bcast(par%CFL)
!call xmpi_bcast(par%ngd)
!call xmpi_bcast(par%nd)
!call xmpi_bcast(par%por)
!call xmpi_bcast(par%wetslp)
!call xmpi_bcast(par%dryslp)
!call xmpi_bcast(par%sw)
!call xmpi_bcast(par%front)
!call xmpi_bcast(par%ARC)
!call xmpi_bcast(par%order)
!call xmpi_bcast(par%left)
!call xmpi_bcast(par%right)
!call xmpi_bcast(par%back)
!call xmpi_bcast(par%refl)
!call xmpi_bcast(par%hswitch)
!call xmpi_bcast(par%z0)
!call xmpi_bcast(par%w)
!call xmpi_bcast(par%compi)
!call xmpi_bcast(par%listline)
!call xmpi_bcast(par%rhoa)
!call xmpi_bcast(par%Cd)
!call xmpi_bcast(par%windv)
!call xmpi_bcast(par%windth)
!call xmpi_bcast(par%nonh)
!call xmpi_bcast(par%nuhv)
!call xmpi_bcast(par%wearth)
!call xmpi_bcast(par%lat)
!call xmpi_bcast(par%fc)
!call xmpi_bcast(par%fcutoff)
!call xmpi_bcast(par%sprdthr)
!call xmpi_bcast(par%smax)
!call xmpi_bcast(par%form)
!call xmpi_bcast(par%carspan)
!call xmpi_bcast(par%nspr)
!call xmpi_bcast(par%thetanum)
!call xmpi_bcast(par%tsfac)
!call xmpi_bcast(par%scheme)
!call xmpi_bcast(par%random)
!call xmpi_bcast(par%trepfac)
!call xmpi_bcast(par%facua)
!call xmpi_bcast(par%dzmax)
!call xmpi_bcast(par%turb)
!call xmpi_bcast(par%rfb)
!call xmpi_bcast(par%lwave)
!call xmpi_bcast(par%swave)
!call xmpi_bcast(par%sws)
!call xmpi_bcast(par%ut)
!call xmpi_bcast(par%Tbfac)
!call xmpi_bcast(par%Tsmin)
!call xmpi_bcast(par%impact)
!call xmpi_bcast(par%CE)
!call xmpi_bcast(par%BRfac)
!call xmpi_bcast(par%betad)
!call xmpi_bcast(par%lwt)
end subroutine distribute_par
#endif
!
! printparams for debugging only
!
subroutine printparams(par,str)
#ifdef USEMPI
  use xmpi_module
#endif
  type (parameters) :: par
  character(*), intent(in) :: str
  integer :: f
  integer :: id
  id = 0
#ifdef USEMPI
  id = xmpi_rank
#endif
  f = 70+id
  if (id .gt. 0) then
    return
  endif

!  write(f,*) 'printparams ', id, ' ', str
!  write(f,*) 'printpar ',id,' ','px:',par%px
!  write(f,*) 'printpar ',id,' ','Hrms:',par%Hrms
!  write(f,*) 'printpar ',id,' ','Trep:',par%Trep
!  write(f,*) 'printpar ',id,' ','dir0:',par%dir0
!  write(f,*) 'printpar ',id,' ','m:',par%m
!!  write(f,*) 'printpar ',id,' ','nt:',par%nt
!  write(f,*) 'printpar ',id,' ','hmin:',par%hmin
!  write(f,*) 'printpar ',id,' ','gammax:',par%gammax
!  write(f,*) 'printpar ',id,' ','Tlong:',par%Tlong
!  write(f,*) 'printpar ',id,' ','Llong:',par%Llong
!  write(f,*) 'printpar ',id,' ','gamma:',par%gamma
!  write(f,*) 'printpar ',id,' ','delta:',par%delta
!  write(f,*) 'printpar ',id,' ','rho:',par%rho
!  write(f,*) 'printpar ',id,' ','g:',par%g
!  write(f,*) 'printpar ',id,' ','rhog8:',par%rhog8
!!  write(f,*) 'printpar ',id,' ','omega:',par%omega
!  write(f,*) 'printpar ',id,' ','thetamin:',par%thetamin
!  write(f,*) 'printpar ',id,' ','thetamax:',par%thetamax
!  write(f,*) 'printpar ',id,' ','dtheta:',par%dtheta
!  write(f,*) 'printpar ',id,' ','thetanaut:',par%thetanaut
!  write(f,*) 'printpar ',id,' ','wci:',par%wci
!  write(f,*) 'printpar ',id,' ','hwci:',par%hwci
!  write(f,*) 'printpar ',id,' ','dt:',par%dt
!  write(f,*) 'printpar ',id,' ','break:',par%break
!  write(f,*) 'printpar ',id,' ','instat:',par%instat
!  write(f,*) 'printpar ',id,' ','wavint:',par%wavint
!  write(f,*) 'printpar ',id,' ','alpha:',par%alpha
!  write(f,*) 'printpar ',id,' ','n:',par%n
!  write(f,*) 'printpar ',id,' ','roller:',par%roller
!  write(f,*) 'printpar ',id,' ','beta:',par%beta
!  write(f,*) 'printpar ',id,' ','taper:',par%taper
!  write(f,*) 'printpar ',id,' ','t:',par%t
!  write(f,*) 'printpar ',id,' ','tnext:',par%tnext
!!  write(f,*) 'printpar ',id,' ','it:',par%it
!  write(f,*) 'printpar ',id,' ','tstart:',par%tstart
!  write(f,*) 'printpar ',id,' ','tint:',par%tint
!  write(f,*) 'printpar ',id,' ','tintp:',par%tintp
!  write(f,*) 'printpar ',id,' ','tintg:',par%tintg
!  write(f,*) 'printpar ',id,' ','tintm:',par%tintm
!  write(f,*) 'printpar ',id,' ','tstop:',par%tstop
!!  write(f,*) 'printpar ',id,' ','ntout:',par%ntout
!  write(f,*) 'printpar ',id,' ','C:',par%C
!  write(f,*) 'printpar ',id,' ','cf:',par%cf
!  write(f,*) 'printpar ',id,' ','eps:',par%eps
!  write(f,*) 'printpar ',id,' ','umin:',par%umin
!  write(f,*) 'printpar ',id,' ','zs01:',par%zs01
!  write(f,*) 'printpar ',id,' ','zs02:',par%zs02
!  write(f,*) 'printpar ',id,' ','zs03:',par%zs03
!  write(f,*) 'printpar ',id,' ','zs04:',par%zs04
!  write(f,*) 'printpar ',id,' ','tideloc:',par%tideloc
!  write(f,*) 'printpar ',id,' ','paulrevere:',par%paulrevere
!  write(f,*) 'printpar ',id,' ','tidelen:',par%tidelen
!!  write(f,*) 'printpar ',id,' ','dico:',par%dico
!  write(f,*) 'printpar ',id,' ','facsl:',par%facsl
!  write(f,*) 'printpar ',id,' ','nuh:',par%nuh
!  write(f,*) 'printpar ',id,' ','nuhfac:',par%nuhfac
!  write(f,*) 'printpar ',id,' ','rhos:',par%rhos
!  write(f,*) 'printpar ',id,' ','morfac:',par%morfac
!  write(f,*) 'printpar ',id,' ','morstart:',par%morstart
!  write(f,*) 'printpar ',id,' ','Emean:',par%Emean
!  write(f,*) 'printpar ',id,' ','CFL:',par%CFL
!  write(f,*) 'printpar ',id,' ','ngd:',par%ngd
!  write(f,*) 'printpar ',id,' ','por:',par%por
!  write(f,*) 'printpar ',id,' ','wetslp:',par%wetslp
!  write(f,*) 'printpar ',id,' ','dryslp:',par%dryslp
!!  write(f,*) 'printpar ',id,' ','sw:',par%sw
!  write(f,*) 'printpar ',id,' ','front:',par%front
!  write(f,*) 'printpar ',id,' ','ARC:',par%ARC
!  write(f,*) 'printpar ',id,' ','order:',par%order
!  write(f,*) 'printpar ',id,' ','left:',par%left
!  write(f,*) 'printpar ',id,' ','right:',par%right
!  write(f,*) 'printpar ',id,' ','back:',par%back
!!  write(f,*) 'printpar ',id,' ','refl:',par%refl
!  write(f,*) 'printpar ',id,' ','hswitch:',par%hswitch
!  write(f,*) 'printpar ',id,' ','z0:',par%z0
!  write(f,*) 'printpar ',id,' ','w:',par%w
!  write(f,*) 'printpar ',id,' ','compi:',par%compi
!  write(f,*) 'printpar ',id,' ','listline:',par%listline
!  write(f,*) 'printpar ',id,' ','rhoa:',par%rhoa
!  write(f,*) 'printpar ',id,' ','Cd:',par%Cd
!  write(f,*) 'printpar ',id,' ','windv:',par%windv
!  write(f,*) 'printpar ',id,' ','windth:',par%windth
!  write(f,*) 'printpar ',id,' ','nonh:',par%nonh
!  write(f,*) 'printpar ',id,' ','nuhv:',par%nuhv
!  write(f,*) 'printpar ',id,' ','wearth:',par%wearth
!  write(f,*) 'printpar ',id,' ','lat:',par%lat
!  write(f,*) 'printpar ',id,' ','fc:',par%fc
!  write(f,*) 'printpar ',id,' ','fcutoff:',par%fcutoff
!  write(f,*) 'printpar ',id,' ','sprdthr:',par%sprdthr
!
!  write(f,*) 'printpar ',id,' ','smax:',par%smax
!  write(f,*) 'printpar ',id,' ','form:',par%form
!  write(f,*) 'printpar ',id,' ','carspan:',par%carspan
!  write(f,*) 'printpar ',id,' ','nspr:',par%nspr
!  write(f,*) 'printpar ',id,' ','thetanum:',par%thetanum
!  write(f,*) 'printpar ',id,' ','tsfac:',par%tsfac
!  write(f,*) 'printpar ',id,' ','scheme:',par%scheme
!  write(f,*) 'printpar ',id,' ','random:',par%random
!  write(f,*) 'printpar ',id,' ','trepfac:',par%trepfac
!  write(f,*) 'printpar ',id,' ','facua:',par%facua
!  write(f,*) 'printpar ',id,' ','dzmax:',par%dzmax
!  write(f,*) 'printpar ',id,' ','turb:',par%turb
!  write(f,*) 'printpar ',id,' ','rfb:',par%rfb
!  write(f,*) 'printpar ',id,' ','lwave:',par%lwave
!  write(f,*) 'printpar ',id,' ','swave:',par%swave
!  write(f,*) 'printpar ',id,' ','sws:',par%sws
!!  write(f,*) 'printpar ',id,' ','ut:',par%ut
!  write(f,*) 'printpar ',id,' ','Tbfac:',par%Tbfac
!  write(f,*) 'printpar ',id,' ','Tsmin:',par%Tsmin
!!  write(f,*) 'printpar ',id,' ','impact:',par%impact
!!  write(f,*) 'printpar ',id,' ','CE:',par%CE
!  write(f,*) 'printpar ',id,' ','BRfac:',par%BRfac
!  write(f,*) 'printpar ',id,' ','betad:',par%betad
!  write(f,*) 'printpar ',id,' ','lwt:',par%lwt
  
end subroutine printparams


end module params
