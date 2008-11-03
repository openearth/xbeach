module params
type parameters
! for debugging, I found it very useful to initialize
! the variables to something like -123   wwvv
real*8     :: px        = -123 ! pi
real*8     :: Hrms      = -123 ! Hrms wave height
real*8     :: Trep      = -123 ! representative wave period 
real*8     :: dir0      = -123 ! mean wave direction (Nautical convention)
integer*4  :: m         = -123 ! power in cos^m directional distribution
integer*4  :: nt        = -123 ! max. number of time steps
real*8     :: hmin      = -123 ! threshold water depth 
real*8     :: gammax    = -123 ! maximum ratio Hrms/hh
real*8     :: Tlong     = -123 ! wave group period for case instat=1
real*8     :: Llong     = -123 ! alongshore wave group length for case instat=1
real*8     :: gamma     = -123 ! breaker parameter in Baldock or Roelvink formulation
real*8     :: delta     = -123 ! fraction of wave height to add to depth in computation of celerity
real*8     :: rho       = -123 ! water density
real*8     :: g         = -123 ! acceleration of gravity
real*8     :: rhog8     = -123 ! 1/8*rho*g
real*8     :: omega     = -123 ! angular wave frequency
real*8     :: thetamin  = -123 ! lower directional limit (angle w.r.t computational x-axis)
real*8     :: thetamax  = -123 ! upper directional limit (angle w.r.t computational x-axis)
real*8     :: dtheta    = -123 ! directional resolution (deg)
integer*4  :: thetanaut = -123 ! option to enter thetamin,thetamax in nautical convention
real*8     :: wci       = -123 ! option wave/current interaction 0/1
real*8     :: hwci      = -123 ! min depth fro wci
real*8     :: dt        = -123 ! time step
integer*4  :: break     = -123 ! option breaker model (1=roelvink, 2=baldock, 3=roelvink adapted)
integer*4  :: instat    = -123 ! option time-varying wave b.c. 
                               ! (0=stationary, 1=regular wave groups, 3=long-crested random wave groups)
integer*4  :: wavint    = -123 ! (only in stationary mode) interval between wave module calls in tint
real*8     :: alpha     = -123 ! wave dissipation coefficient
real*8     :: n         = -123 ! power in roelvink dissipation model
integer*4  :: roller    = -123 ! option to turn off/on roller model (0/1) (not implemented yet)
real*8     :: beta      = -123 ! breaker slope coefficient in roller model
real*8     :: taper     = -123 ! time to spin up wave b.c. in case of stationary waves
real*8     :: t         = -123 ! time (s)
real*8     :: tnext     = -123 ! next time point for output
integer*4  :: it        = -123 ! output time step number
real*8     :: tstart    = -123 ! start time of simulation output
real*8     :: tint      = -123 ! time interval output
real*8     :: tintp     = -123 ! time interval output points
real*8     :: tintg     = -123 ! time interval output global variables
real*8     :: tintm     = -123 ! time interval output mean global variables
real*8     :: tstop     = -123 ! stop time simulation
integer*4  :: ntout     = -123 ! number of output time steps
real*8     :: C         = -123 ! Chezy value
real*8     :: cf        = -123 ! friction coefficient flow [-]
real*8     :: eps       = -123 ! threshold depth
real*8     :: umin      = -123 ! threshold velocity upwind scheme
real*8     :: zs01      = -123 ! initial water level first sea boundary
real*8     :: zs02      = -123 ! initial water level second sea boundary
real*8     :: zs03      = -123 ! initial water level first land boundary
real*8     :: zs04      = -123 ! initial water level second land boundary
integer*4  :: tideloc   = -123 ! number of input tidal time series
real*8     :: paulrevere= -123 ! if tideloc =>2, then this indicates where the time series are to be 
                               ! applied. Input for tidal information to xbeach options (3):
                               ! 1. one tidal record --> specify tidal record everywhere
                               ! 2. two tidal records --> Need to specify keyword 'paulrevere'
                               ! paulrevere==0 implies to apply one tidal record to
                               ! both sea corners and one tidal record to both land corners
                               ! paulrevere==1 implies to apply the first tidal record
                               ! (column 2 in zs0input.dat) to the (x=1,y=1) sea corner and
                               ! the second tidal record (third column) to the (x=1,y=N) sea corner
                               ! 3. four tidal records --> Need to list tidal records in  
                               ! zs0input.dat in order of:
                               !      (x=1,y=1)
                               !      (x=1,y=N)
                               !      (x=N,y=N)
                               !      (x=N,y=1)
                               !      NOTE:  clockwise from (1,1) corner 
integer*4  :: tidelen   = -123 ! length of input tidal time series 
real*8     :: A         = -123 ! obsolete
real*8     :: dico      = -123 ! diffusion coefficient
real*8     :: facsl     = -123 ! factor bedslope effect
real*8     :: nuh       = -123 ! horizontal background viscosity 
real*8     :: nuhfac    = -123 ! viscosity coefficient for roller induced turbulent horizontal viscosity
real*8     :: rhos      = -123 ! sediment density
real*8     :: morfac    = -123 ! morphological factor
real*8     :: morstart  = -123 ! start time morphology
real*8     :: Emean     = -123 ! mean wave energy at boundary
real*8     :: CFL       = -123 ! maximum courant number
integer*4  :: ngd       = -123 ! number of sediment classes
integer*4  :: nd        = -123 ! number of sediment class layers
real*8     :: dzg       = -123 ! thickness of sediment class layers
real*8     :: D501      = -123 ! D50 grain diameter second class of sediment
real*8     :: D901      = -123 ! D90 grain diameter second class of sediment
real*8     :: D502      = -123 ! D50 grain diameter second class of sediment
real*8     :: D902      = -123 ! D90 grain diameter second class of sediment
real*8     :: D503      = -123 ! D50 grain diameter third class of sediment
real*8     :: D903      = -123 ! D90 grain diameter third class of sediment
real*8     :: sedcal1   = -123 ! calibration factor for sediment class 1 
real*8     :: sedcal2   = -123 ! calibration factor for sediment class 2 
real*8     :: sedcal3   = -123 ! calibration factor for sediment class 2 
real*8     :: por       = -123 ! porosity
real*8     :: wetslp    = -123 ! critical avalanching slope under water
real*8     :: dryslp    = -123 ! critical avalanching slope above water
integer*4  :: sw        = -123 ! short wave contribution: 0 = urms=0 & ust =0, 1 = default model, 2 = urms=1 and ust=0, 3 = urms=0 & ust=1 :: Not used
integer*4  :: front     = -123 ! switch for seaward flow boundary: 0 = radiating boundary(Ad), 1 = Van Dongeren, 1997
integer*4  :: ARC       = -123 ! switch for active reflection compensation at seaward boundary: 0 = reflective, 1 = weakly (non) reflective
real*4     :: order     = -123 ! switch for order of wave steering, 1 = first order wave steering (short wave energy only), 2 = second oder wave steering (bound long wave corresponding to short wave forcing is added)
! wwvv: why is order real*4 and not integer?
integer*4  :: left      = -123 ! switch for lateral boundary at left, 0 = vv computed from NSWE, 1 = reflective wall; vv=0
integer*4  :: right     = -123 ! switch for lateral boundary at right, 0 = vv computed from NSWE, 1 = reflective wall; vv=0
integer*4  :: back      = -123 ! switch for boundary at bay side, 0 = radiating boundary (Ad), 1 = reflective boundary; uu=0
integer*4  :: refl      = -123 ! 1 = compensate for reflected wave and roller massflux, 0 = no compensation
real*8     :: hswitch   = -123 ! is the water depth at which is switched from wetslp to dryslp
real*8     :: z0        = -123 ! zero flow velocity level in Soulsby van Rijn (1997) sed.conc. expression
real*8     :: w         = -123 ! fall velocity sediment
complex(kind(0.0d0)):: compi = -123       ! complex i, sqrt(-1)
integer*4  :: listline  = -123 ! keeps rack of the record line in bcf-files 
real*8     :: rhoa      = -123 ! air density
real*8     :: Cd        = -123 ! wind drag coefficient
real*8     :: windv     = -123 ! wind velocity
real*8     :: windth    = -123 ! wind direction (nautical input)
real*8     :: epsi      = -123 ! weighting factor for actual flow in computing time avergaed flow at seawrd boundary 1>=epsi>=0
integer*4  :: nonh      = -123 ! 0 = NSWE, 1 = NSW + non-hydrostatic pressure compensation Stelling & Zijlema, 2003 
real*8     :: nuhv      = -123 ! longshore viscosity enhancement factor
real*8     :: wearth    = -123 ! angular velocity of earth for computing Coriolis forces
real*8     :: lat       = -123 ! estimated latitude at model location  for computing Coriolis
real*8     :: fc        = -123 !
real*8     :: fcutoff   = -123 ! lo freq cutoff frequency for boundary conditions
real*8     :: sprdthr   = -123 ! threshold above which spec dens are read in (default 0.08*maxval)
real*8     :: struct    = -123 ! 0 = no revetment, 1 = multiple sediment classes with non-erodable fractions
real*8     :: smax      = -123 ! Being tested: maximum Shields parameter for ceq   Diane Foster
integer*4  :: form      = -123 ! equilibrium sed. conc. formulation: 1 = Soulsby van rijn, 1997, 2 = Van Rijn 2008
integer*4  :: carspan   = -123 ! 0 = use cg (default); 1 = use sqrt(gh) in instat = 3 for c&g tests
!integer*4  :: rugauge  = -123  ! 0 = normal obs. point (default) ; 1 = runupgauge obs. point moving with the shoreline.
integer*4  :: nspr      = -123 ! Expert tool: nspr = 1 bin all wave components for generation of qin (instat 4+) in one direction
                               !              nspr = 0 regular long wave spreading (default)
real*8     :: thetanum  = -123 ! Coefficient determining whether upwind (1) or central scheme (0.5) is used.
real*8     :: tsfac     = -123 ! Coefficient determining Ts = tsfac * h/ws in sediment source term
integer*4  :: scheme    = -123 ! Numerical scheme for wave and roller energy : 1=upwind, 2=Lax-Wendroff
integer*4  :: random    = -123 ! 1 = random seed, 0 = seed values are zero
real*8     :: trepfac   = -123 ! compute mean wave period over energy band par%trepfac*maxval(Sf); converges to Tm01 for trepfac = 0.0 and to Tp for trepfac = 1.0 
real*8     :: facua     = -123 ! calibration factor time averaged flows due to wave asymmetry
real*8     :: dzmax     = -123 ! maximum bedlevel change due to avalanching [m/s/m]
integer*4  :: turb      = -123 ! equlibrium sediment concentration is computed as function of:
                               ! 0 = no turbulence, 1 = wave averaged turbulence, 2 = maximum turbulence
integer*4  :: rfb      = -123  ! if rfb = 1 then maximum wave surface slope is feeded back in roller energy balance; else rfb = par%Beta
integer*4  :: lwave    = -123  ! 1 = long waves, 0 = no long waves
integer*4  :: swave    = -123  ! 1 = short waves, 0 = no short waves
integer*4  :: sws      = -123  ! 1 = short wave & roller undertow, 0 = no short wave & roller undertow
integer*4  :: ut       = -123  ! 1 = short wave up-stirring, 0 = no short wave up-stirring
real*8     :: Tbfac    = -123  ! Calibration factor for bore interval Tbore: Tbore = Tbfac*Tbore
real*8     :: Tsmin    = -123  ! Minimum adaptation time scale in advection diffusion equation sediment
real*8     :: impact   = -123  ! Include Fisher Overton approach in avalanching
real*8     :: CE       = -123  ! Dune face erosion coefficient for avalanching
real*8     :: BRfac    = -123  ! calibration factor surface slope


end type parameters

contains

subroutine wave_input(par)
use readkey_module
use xmpi_module
implicit none
type(parameters)            :: par

character(len=80)          :: dummystring

par%px    = 3.14159265358979d0
par%compi = (0.0d0,1.0d0)

par%instat   = readkey_int     ('params.txt','instat',    1,         0,        7)
par%fcutoff  = readkey_dbl     ('params.txt','fcutoff',   0.d0,      0.d0,     40.d0)
par%random   = readkey_int     ('params.txt','random',    0,         0,        1)
if (par%instat == 0) then
    par%dir0  = readkey_dbl    ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
    par%Hrms  = readkey_dbl    ('params.txt','Hrms',      1.d0,      0.d0,    10.d0)
    par%wavint   = readkey_int ('params.txt','wavint',    1,         1,     3600)
    par%m     = readkey_int    ('params.txt','m',        10,         2,      128)
    par%Trep  = readkey_dbl    ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
    par%Trep  = readkey_dbl    ('params.txt','Trep',     par%Trep,   1.d0,    20.d0)
    par%omega    = 2.d0*par%px/par%Trep;
elseif (par%instat==1) then
    par%dir0  = readkey_dbl    ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
    par%Hrms  = readkey_dbl    ('params.txt','Hrms',      1.d0,      0.d0,    10.d0)
    par%Tlong = readkey_dbl    ('params.txt','Tlong',    80.d0,     20.d0,   300.d0)
    par%m     = readkey_int    ('params.txt','m',        10,         2,      128)
    par%Trep  = readkey_dbl    ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
    par%Trep  = readkey_dbl    ('params.txt','Trep',     par%Trep,   1.d0,    20.d0)
    par%omega    = 2.d0*par%px/par%Trep;
elseif (par%instat==2 .or. par%instat==3) then
    par%dir0  = readkey_dbl    ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
    par%Hrms  = readkey_dbl    ('params.txt','Hrms',      1.d0,      0.d0,    10.d0)
    par%m     = readkey_int    ('params.txt','m',        10,         2,      128)
    par%Trep  = readkey_dbl    ('params.txt','Tm01',     10.d0,      1.d0,    20.d0)
    par%Trep  = readkey_dbl    ('params.txt','Trep',     par%Trep,   1.d0,    20.d0)
    par%omega    = 2.d0*par%px/par%Trep;
elseif (par%instat==4 .or. par%instat==5 .or. par%instat==6) then
        ! Just a check .....
        if (xmaster) then
          call readkey('params.txt','bcfile',dummystring)
          call readkey('params.txt','rt',dummystring)
          call readkey('params.txt','dtbc',dummystring)
          call readkey('params.txt','dthetaS_XB',dummystring)
        endif
elseif (par%instat > 8.and. par%instat/=41) then
    if(xmaster) then
      write(*,*)'Instat invalid option'
    endif
    call halt_program
end if
!                         Input file  Keyword Default  Minimum  Maximum
!par%dir0  = readkey_dbl    ('params.txt','dir0',    270.d0,    180.d0,   360.d0)
par%hmin  = readkey_dbl ('params.txt','hmin',   0.01d0,   0.001d0,      1.d0)
par%gammax= readkey_dbl ('params.txt','gammax',   5.d0,      .4d0,      5.d0)
par%gamma = readkey_dbl ('params.txt','gamma',   0.6d0,     0.4d0,     0.9d0)
par%alpha = readkey_dbl ('params.txt','alpha',   1.0d0,     0.5d0,     2.0d0)
par%delta = readkey_dbl ('params.txt','delta',   0.0d0,     0.0d0,     1.0d0)
par%n    =  readkey_dbl ('params.txt','n',       5.0d0,     5.0d0,    20.0d0)
par%rho   = readkey_dbl ('params.txt','rho',  1025.0d0,  1000.0d0,  1040.0d0)
par%g     = readkey_dbl ('params.txt','g',      9.81d0,     9.7d0,     9.9d0)
par%rhog8 = 1.0d0/8.0d0*par%rho*par%g
par%Emean = par%rhog8*par%Hrms**2
par%thetamin = readkey_dbl ('params.txt','thetamin', -80.d0,    -180.d0,  180.d0)
par%thetamax = readkey_dbl ('params.txt','thetamax',  80.d0,    -180.d0,  180.d0)
par%dtheta   = readkey_dbl ('params.txt','dtheta',    10.d0,      0.1d0,   20.d0)
par%thetanaut= readkey_int ('params.txt','thetanaut',    0,        0,     1)
par%wci      = readkey_int ('params.txt','wci',        0,        0,     1)
par%hwci  = readkey_dbl ('params.txt','hwci',   0.01d0,   0.001d0,      1.d0)
par%break    = readkey_int ('params.txt','break',      3,        1,     3)
par%roller   = readkey_int ('params.txt','roller',     1,        0,     1)
par%beta     = readkey_dbl ('params.txt','beta',    0.15d0,     0.05d0,   0.3d0)
par%rfb      = readkey_int ('params.txt','rfb',        1,        0,     1)
par%taper    = readkey_dbl ('params.txt','taper',   100.d0,      0.0d0, 1000.d0)
par%taper=max(par%taper,1.d-6)
par%refl     = readkey_int ('params.txt','refl',       0,        0,     1) 
par%nspr     = readkey_int ('params.txt','nspr',       0,        0,     1) 
par%scheme   = readkey_int ('params.txt','scheme',     1,        1,     2)
par%trepfac  = readkey_dbl  ('params.txt','sprdthr', 0.8d0,       0.d0,  1.d0) 
par%sprdthr  = readkey_dbl ('params.txt','sprdthr', 0.08d0,      0.d0,  1.d0) 
par%lwave    = readkey_int ('params.txt','lwave',         1,        0,     1)
par%swave    = readkey_int ('params.txt','swave',         1,        0,     1)
par%sws      = readkey_int ('params.txt','sws',           1,        0,     1)
par%ut       = readkey_int ('params.txt','ut',            1,        0,     1)
par%BRfac    = readkey_dbl ('params.txt','BRfac',    1.0d0,       0.d0, 1.d0)
end subroutine wave_input

subroutine flow_input(par)
use readkey_module
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
par%tint    = readkey_dbl ('params.txt','tint',     1.d0,     .01d0, 100000.d0)  ! Robert
par%tintg   = readkey_dbl ('params.txt','tintg', par%tint,     .01d0, 100000.d0)  ! Robert
par%tintp   = readkey_dbl ('params.txt','tintp',par%tintg,    .01d0, 100000.d0)  ! Robert
par%tintm   = readkey_dbl ('params.txt','tintm',par%tintg,    1.d0, par%tstop)   ! Robert
par%tint    = min(par%tintg,par%tintp,par%tintm)                                 ! Robert                                                                                                ! Robert
par%tstop   = readkey_dbl ('params.txt','tstop', 2000.d0,      1.d0,1000000.d0)
! adapt flow times to morfac
par%morfac   = readkey_dbl ('params.txt','morfac', 0.0d0,        0.d0,  1000.d0)
par%tstart  = par%tstart / max(par%morfac,1.d0)
par%tint    = par%tint   / max(par%morfac,1.d0)
par%tintg   = par%tintg  / max(par%morfac,1.d0)
par%tintp   = par%tintp  / max(par%morfac,1.d0)
par%tintm   = par%tintm  / max(par%morfac,1.d0)
par%tint    = par%tint   / max(par%morfac,1.d0)
par%tstop   = par%tstop  / max(par%morfac,1.d0)
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
par%back    = readkey_int ('params.txt','back',        2,         0,      2)
par%nuh     = readkey_dbl ('params.txt','nuh',     0.5d0,     0.0d0,      1.0d0)
par%nuhfac  = readkey_dbl ('params.txt','nuhfac',      0.0d0,     0.0d0,  1.0d0)
par%rhoa    = readkey_dbl ('params.txt','rhoa',   1.25d0,     1.0d0,   2.0d0)
par%Cd      = readkey_dbl ('params.txt','Cd',    0.002d0,  0.0001d0,  0.01d0)
par%windv   = readkey_dbl ('params.txt','windv',   0.0d0,     0.0d0, 200.0d0)
par%windth  = readkey_dbl ('params.txt','windth', 90.0d0,  -180.0d0, 180.0d0)
par%epsi    = readkey_dbl ('params.txt','epsi',     0.d0,      0.d0,    1.d0)
par%nonh    = readkey_int ('params.txt','nonh',        0,         0,      1)
par%nuhv    = readkey_dbl ('params.txt','nuhv',     1.d0,      1.d0,    20.d0)
par%lat     = readkey_dbl ('params.txt','lat',     0.d0,      0.d0,   90.d0)
par%wearth  = readkey_dbl ('params.txt','omega',   1.d0/24.d0, 0.d0,    1.d0)
! Convert from Nautical to cartesian convention
par%windth=(270.d0-par%windth)*par%px/180.d0
par%lat = par%lat*par%px/180.d0
par%wearth = par%px*par%wearth/1800.d0
par%fc = 2.d0*par%wearth*sin(par%lat)

end subroutine flow_input

subroutine sed_input(par)
use readkey_module
implicit none
type(parameters)            :: par


par%dico     = readkey_dbl ('params.txt','dico',    1.d0,        0.d0,    10.d0)
par%ngd      = readkey_int ('params.txt','ngd',        1,           1,        2)
par%nd       = readkey_int ('params.txt','nd',        1,           1,        20)
par%dzg      = readkey_dbl ('params.txt','dzg',    0.1d0,      0.01d0,     1.d0)
par%D501     = readkey_dbl ('params.txt','D50', 0.0002d0,   0.00005d0,  0.001d0)
par%D901     = readkey_dbl ('params.txt','D90', 0.0003d0,   0.00005d0,  0.001d0)
par%D502     = readkey_dbl ('params.txt','D502', 0.0000d0,   0.0000d0,  0.001d0)
par%D902     = readkey_dbl ('params.txt','D902', 0.0000d0,   0.0000d0,  0.001d0)
par%D503     = readkey_dbl ('params.txt','D503', 0.0000d0,   0.0000d0,  0.001d0)
par%D903     = readkey_dbl ('params.txt','D903', 0.0000d0,   0.0000d0,  0.001d0)
par%sedcal1  = readkey_dbl ('params.txt','sedcal1', 1.0000d0,   0.0000d0, 10.00d0)
par%sedcal2  = readkey_dbl ('params.txt','sedcal2', 1.0000d0,   0.0000d0, 10.00d0)
par%sedcal3  = readkey_dbl ('params.txt','sedcal3', 1.0000d0,   0.0000d0, 10.00d0)
par%rhos     = readkey_dbl ('params.txt','rhos',  2650d0,     2400.d0,  2800.d0)
par%morfac   = readkey_dbl ('params.txt','morfac', 0.0d0,        0.d0,  1000.d0)
par%morstart = readkey_dbl ('params.txt','morstart',120.d0,      0.d0, 10000.d0)
par%morstart = par%morstart / max(par%morfac,1.d0)
par%wetslp   = readkey_dbl ('params.txt','wetslp', 0.3d0,       0.1d0,     1.d0)
par%dryslp   = readkey_dbl ('params.txt','dryslp', 1.0d0,       0.1d0,     2.d0)
par%por      = readkey_dbl ('params.txt','por',    0.4d0,       0.3d0,    0.5d0)
par%hswitch  = readkey_dbl ('params.txt','hswitch',0.1d0,      0.01d0,    1.0d0)  
par%z0       = readkey_dbl ('params.txt','z0     ',0.006d0,    0.0001d0,   0.05d0)  
par%facsl    = readkey_dbl ('params.txt','facsl  ',0.00d0,    0.00d0,   1.6d0)  
par%struct   = readkey_int ('params.txt','struct ',0,           0,            1) 
par%form     = readkey_int ('params.txt','form',   1,           1,            3)
par%smax     = readkey_dbl ('params.txt','smax',   -1.d0,    0.d0,   3.d0)
if (par%smax<0) par%smax=huge(0.d0)
par%thetanum = readkey_dbl ('params.txt','thetanum',   1.d0,    0.5d0,   1.d0) 
par%tsfac    = readkey_dbl ('params.txt','tsfac',   0.1d0,    0.01d0,   1.d0) 
par%facua    = readkey_dbl ('params.txt','facua  ',0.00d0,    0.00d0,   1.0d0) 
par%dzmax    = readkey_dbl ('params.txt','dzmax  ',0.05d0,    0.00d0,   1.0d0) 
par%turb     = readkey_int ('params.txt','turb',   2,           0,            2)
par%Tbfac    = readkey_dbl ('params.txt','Tbfac  ',1.0d0,     0.00d0,   1.0d0) 
par%Tsmin    = readkey_dbl ('params.txt','Tsmin  ',0.2d0,     0.01d0,   10.d0) 
par%impact   = readkey_int ('params.txt','impact ',0,           0,            1)  
par%CE       = readkey_dbl ('params.txt','CE     ',0.2d0,     0.00d0,   100.d0) 

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

inquire(iolength=parlen) par
inquire(iolength=w) 1.0d0

! convert parlen to number of bytes, assuming that 
! 1.0d0 takes 8 bytes

parlen = (8/w)*parlen

call MPI_Bcast(par,parlen,MPI_BYTE,xmpi_master,xmpi_comm,ierror)
return

! so, the following code is NOT used anymore. I left this here
! maybe method above does not work everywhere. wwvv

! For efficiency, this subroutine should use MPI_Pack and 
! MPI_Unpack. However, since this subroutine is only called a 
! few times, a more simple approach is used.
!

call xmpi_bcast(par%px)
call xmpi_bcast(par%Hrms)
call xmpi_bcast(par%Trep)
call xmpi_bcast(par%dir0)
call xmpi_bcast(par%m)
call xmpi_bcast(par%nt)
call xmpi_bcast(par%hmin)
call xmpi_bcast(par%gammax)
call xmpi_bcast(par%Tlong)
call xmpi_bcast(par%Llong)
call xmpi_bcast(par%gamma)
call xmpi_bcast(par%delta)
call xmpi_bcast(par%rho)
call xmpi_bcast(par%g)
call xmpi_bcast(par%rhog8)
call xmpi_bcast(par%omega)
call xmpi_bcast(par%thetamin)
call xmpi_bcast(par%thetamax)
call xmpi_bcast(par%dtheta)
call xmpi_bcast(par%thetanaut)
call xmpi_bcast(par%wci)
call xmpi_bcast(par%hwci)
call xmpi_bcast(par%dt)
call xmpi_bcast(par%break)
call xmpi_bcast(par%instat)
call xmpi_bcast(par%wavint)
call xmpi_bcast(par%alpha)
call xmpi_bcast(par%n)
call xmpi_bcast(par%roller)
call xmpi_bcast(par%beta)
call xmpi_bcast(par%taper)
call xmpi_bcast(par%t)
call xmpi_bcast(par%tnext)
call xmpi_bcast(par%it)
call xmpi_bcast(par%tstart)
call xmpi_bcast(par%tint)
call xmpi_bcast(par%tintp)
call xmpi_bcast(par%tintg)
call xmpi_bcast(par%tintm)
call xmpi_bcast(par%tstop)
call xmpi_bcast(par%ntout)
call xmpi_bcast(par%C)
call xmpi_bcast(par%cf)
call xmpi_bcast(par%eps)
call xmpi_bcast(par%umin)
call xmpi_bcast(par%zs01)
call xmpi_bcast(par%zs02)
call xmpi_bcast(par%zs03)
call xmpi_bcast(par%zs04)
call xmpi_bcast(par%tideloc)
call xmpi_bcast(par%paulrevere)
call xmpi_bcast(par%tidelen)
call xmpi_bcast(par%A)
call xmpi_bcast(par%dico)
call xmpi_bcast(par%facsl)
call xmpi_bcast(par%nuh)
call xmpi_bcast(par%nuhfac)
call xmpi_bcast(par%rhos)
call xmpi_bcast(par%morfac)
call xmpi_bcast(par%morstart)
call xmpi_bcast(par%Emean)
call xmpi_bcast(par%CFL)
call xmpi_bcast(par%ngd)
call xmpi_bcast(par%nd)
call xmpi_bcast(par%dzg)
call xmpi_bcast(par%D501)
call xmpi_bcast(par%D901)
call xmpi_bcast(par%D502)
call xmpi_bcast(par%D902)
call xmpi_bcast(par%D503)
call xmpi_bcast(par%D903)
call xmpi_bcast(par%sedcal1)
call xmpi_bcast(par%sedcal2)
call xmpi_bcast(par%sedcal3)
call xmpi_bcast(par%por)
call xmpi_bcast(par%wetslp)
call xmpi_bcast(par%dryslp)
call xmpi_bcast(par%sw)
call xmpi_bcast(par%front)
call xmpi_bcast(par%ARC)
call xmpi_bcast(par%order)
call xmpi_bcast(par%left)
call xmpi_bcast(par%right)
call xmpi_bcast(par%back)
call xmpi_bcast(par%refl)
call xmpi_bcast(par%hswitch)
call xmpi_bcast(par%z0)
call xmpi_bcast(par%w)
call xmpi_bcast(par%compi)
call xmpi_bcast(par%listline)
call xmpi_bcast(par%rhoa)
call xmpi_bcast(par%Cd)
call xmpi_bcast(par%windv)
call xmpi_bcast(par%windth)
call xmpi_bcast(par%epsi)
call xmpi_bcast(par%nonh)
call xmpi_bcast(par%nuhv)
call xmpi_bcast(par%wearth)
call xmpi_bcast(par%lat)
call xmpi_bcast(par%fc)
call xmpi_bcast(par%fcutoff)
call xmpi_bcast(par%sprdthr)
call xmpi_bcast(par%struct)
call xmpi_bcast(par%smax)
call xmpi_bcast(par%form)
call xmpi_bcast(par%carspan)
call xmpi_bcast(par%nspr)
call xmpi_bcast(par%thetanum)
call xmpi_bcast(par%tsfac)
call xmpi_bcast(par%scheme)
call xmpi_bcast(par%random)
call xmpi_bcast(par%trepfac)
call xmpi_bcast(par%facua)
call xmpi_bcast(par%dzmax)
call xmpi_bcast(par%turb)
call xmpi_bcast(par%rfb)
call xmpi_bcast(par%lwave)
call xmpi_bcast(par%swave)
call xmpi_bcast(par%sws)
call xmpi_bcast(par%ut)
call xmpi_bcast(par%Tbfac)
call xmpi_bcast(par%Tsmin)
call xmpi_bcast(par%impact)
call xmpi_bcast(par%CE)
call xmpi_bcast(par%BRfac)
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

  write(f,*) 'printparams ', id, ' ', str
  write(f,*) 'printpar ',id,' ','px:',par%px
  write(f,*) 'printpar ',id,' ','Hrms:',par%Hrms
  write(f,*) 'printpar ',id,' ','Trep:',par%Trep
  write(f,*) 'printpar ',id,' ','dir0:',par%dir0
  write(f,*) 'printpar ',id,' ','m:',par%m
  write(f,*) 'printpar ',id,' ','nt:',par%nt
  write(f,*) 'printpar ',id,' ','hmin:',par%hmin
  write(f,*) 'printpar ',id,' ','gammax:',par%gammax
  write(f,*) 'printpar ',id,' ','Tlong:',par%Tlong
  write(f,*) 'printpar ',id,' ','Llong:',par%Llong
  write(f,*) 'printpar ',id,' ','gamma:',par%gamma
  write(f,*) 'printpar ',id,' ','delta:',par%delta
  write(f,*) 'printpar ',id,' ','rho:',par%rho
  write(f,*) 'printpar ',id,' ','g:',par%g
  write(f,*) 'printpar ',id,' ','rhog8:',par%rhog8
  write(f,*) 'printpar ',id,' ','omega:',par%omega
  write(f,*) 'printpar ',id,' ','thetamin:',par%thetamin
  write(f,*) 'printpar ',id,' ','thetamax:',par%thetamax
  write(f,*) 'printpar ',id,' ','dtheta:',par%dtheta
  write(f,*) 'printpar ',id,' ','thetanaut:',par%thetanaut
  write(f,*) 'printpar ',id,' ','wci:',par%wci
  write(f,*) 'printpar ',id,' ','hwci:',par%hwci
  write(f,*) 'printpar ',id,' ','dt:',par%dt
  write(f,*) 'printpar ',id,' ','break:',par%break
  write(f,*) 'printpar ',id,' ','instat:',par%instat
  write(f,*) 'printpar ',id,' ','wavint:',par%wavint
  write(f,*) 'printpar ',id,' ','alpha:',par%alpha
  write(f,*) 'printpar ',id,' ','n:',par%n
  write(f,*) 'printpar ',id,' ','roller:',par%roller
  write(f,*) 'printpar ',id,' ','beta:',par%beta
  write(f,*) 'printpar ',id,' ','taper:',par%taper
  write(f,*) 'printpar ',id,' ','t:',par%t
  write(f,*) 'printpar ',id,' ','tnext:',par%tnext
  write(f,*) 'printpar ',id,' ','it:',par%it
  write(f,*) 'printpar ',id,' ','tstart:',par%tstart
  write(f,*) 'printpar ',id,' ','tint:',par%tint
  write(f,*) 'printpar ',id,' ','tintp:',par%tintp
  write(f,*) 'printpar ',id,' ','tintg:',par%tintg
  write(f,*) 'printpar ',id,' ','tintm:',par%tintm
  write(f,*) 'printpar ',id,' ','tstop:',par%tstop
  write(f,*) 'printpar ',id,' ','ntout:',par%ntout
  write(f,*) 'printpar ',id,' ','C:',par%C
  write(f,*) 'printpar ',id,' ','cf:',par%cf
  write(f,*) 'printpar ',id,' ','eps:',par%eps
  write(f,*) 'printpar ',id,' ','umin:',par%umin
  write(f,*) 'printpar ',id,' ','zs01:',par%zs01
  write(f,*) 'printpar ',id,' ','zs02:',par%zs02
  write(f,*) 'printpar ',id,' ','zs03:',par%zs03
  write(f,*) 'printpar ',id,' ','zs04:',par%zs04
  write(f,*) 'printpar ',id,' ','tideloc:',par%tideloc
  write(f,*) 'printpar ',id,' ','paulrevere:',par%paulrevere
  write(f,*) 'printpar ',id,' ','tidelen:',par%tidelen
  write(f,*) 'printpar ',id,' ','A:',par%A
  write(f,*) 'printpar ',id,' ','dico:',par%dico
  write(f,*) 'printpar ',id,' ','facsl:',par%facsl
  write(f,*) 'printpar ',id,' ','nuh:',par%nuh
  write(f,*) 'printpar ',id,' ','nuhfac:',par%nuhfac
  write(f,*) 'printpar ',id,' ','rhos:',par%rhos
  write(f,*) 'printpar ',id,' ','morfac:',par%morfac
  write(f,*) 'printpar ',id,' ','morstart:',par%morstart
  write(f,*) 'printpar ',id,' ','Emean:',par%Emean
  write(f,*) 'printpar ',id,' ','CFL:',par%CFL
  write(f,*) 'printpar ',id,' ','ngd:',par%ngd
  write(f,*) 'printpar ',id,' ','nd:',par%nd
  write(f,*) 'printpar ',id,' ','dzg:',par%dzg
  write(f,*) 'printpar ',id,' ','D501:',par%D501
  write(f,*) 'printpar ',id,' ','D901:',par%D901
  write(f,*) 'printpar ',id,' ','D502:',par%D502
  write(f,*) 'printpar ',id,' ','D902:',par%D902
  write(f,*) 'printpar ',id,' ','D503:',par%D503
  write(f,*) 'printpar ',id,' ','D903:',par%D903
  write(f,*) 'printpar ',id,' ','sedcal1:',par%sedcal1
  write(f,*) 'printpar ',id,' ','sedcal2:',par%sedcal2
  write(f,*) 'printpar ',id,' ','sedcal3:',par%sedcal3
  write(f,*) 'printpar ',id,' ','por:',par%por
  write(f,*) 'printpar ',id,' ','wetslp:',par%wetslp
  write(f,*) 'printpar ',id,' ','dryslp:',par%dryslp
  write(f,*) 'printpar ',id,' ','sw:',par%sw
  write(f,*) 'printpar ',id,' ','front:',par%front
  write(f,*) 'printpar ',id,' ','ARC:',par%ARC
  write(f,*) 'printpar ',id,' ','order:',par%order
  write(f,*) 'printpar ',id,' ','left:',par%left
  write(f,*) 'printpar ',id,' ','right:',par%right
  write(f,*) 'printpar ',id,' ','back:',par%back
  write(f,*) 'printpar ',id,' ','refl:',par%refl
  write(f,*) 'printpar ',id,' ','hswitch:',par%hswitch
  write(f,*) 'printpar ',id,' ','z0:',par%z0
  write(f,*) 'printpar ',id,' ','w:',par%w
  write(f,*) 'printpar ',id,' ','compi:',par%compi
  write(f,*) 'printpar ',id,' ','listline:',par%listline
  write(f,*) 'printpar ',id,' ','rhoa:',par%rhoa
  write(f,*) 'printpar ',id,' ','Cd:',par%Cd
  write(f,*) 'printpar ',id,' ','windv:',par%windv
  write(f,*) 'printpar ',id,' ','windth:',par%windth
  write(f,*) 'printpar ',id,' ','epsi:',par%epsi
  write(f,*) 'printpar ',id,' ','nonh:',par%nonh
  write(f,*) 'printpar ',id,' ','nuhv:',par%nuhv
  write(f,*) 'printpar ',id,' ','wearth:',par%wearth
  write(f,*) 'printpar ',id,' ','lat:',par%lat
  write(f,*) 'printpar ',id,' ','fc:',par%fc
  write(f,*) 'printpar ',id,' ','fcutoff:',par%fcutoff
  write(f,*) 'printpar ',id,' ','sprdthr:',par%sprdthr
  write(f,*) 'printpar ',id,' ','struct:',par%struct
  write(f,*) 'printpar ',id,' ','smax:',par%smax
  write(f,*) 'printpar ',id,' ','form:',par%form
  write(f,*) 'printpar ',id,' ','carspan:',par%carspan
  write(f,*) 'printpar ',id,' ','nspr:',par%nspr
  write(f,*) 'printpar ',id,' ','thetanum:',par%thetanum
  write(f,*) 'printpar ',id,' ','tsfac:',par%tsfac
  write(f,*) 'printpar ',id,' ','scheme:',par%scheme
  write(f,*) 'printpar ',id,' ','random:',par%random
  write(f,*) 'printpar ',id,' ','trepfac:',par%trepfac
  write(f,*) 'printpar ',id,' ','facua:',par%facua
  write(f,*) 'printpar ',id,' ','dzmax:',par%dzmax
  write(f,*) 'printpar ',id,' ','turb:',par%turb
  write(f,*) 'printpar ',id,' ','rfb:',par%rfb
  write(f,*) 'printpar ',id,' ','lwave:',par%lwave
  write(f,*) 'printpar ',id,' ','swave:',par%swave
  write(f,*) 'printpar ',id,' ','sws:',par%sws
  write(f,*) 'printpar ',id,' ','ut:',par%ut
  write(f,*) 'printpar ',id,' ','Tbfac:',par%Tbfac
  write(f,*) 'printpar ',id,' ','Tsmin:',par%Tsmin
  write(f,*) 'printpar ',id,' ','impact:',par%impact
  write(f,*) 'printpar ',id,' ','CE:',par%CE
  write(f,*) 'printpar ',id,' ','BRfac:',par%BRfac
  
end subroutine printparams
end module params
