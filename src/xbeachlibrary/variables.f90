module variables
  double precision, allocatable, target :: x(:,:)           !< [m] x-coord. original cmp. grid {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "projection_x_coordinate", "broadcast": "d"}
  double precision, allocatable, target :: y(:,:)           !< [m] y-coord. original comp. grid {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "projection_y_coordinate", "broadcast": "d"}
  double precision, allocatable, target :: xz(:,:)          !< [m] x-coord. comp. grid (positive shoreward, perp. to coastline) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "projection_x_coordinate", "broadcast": "d"}
  double precision, allocatable, target :: yz(:,:)          !< [m] y-coord. comp. grid {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "projection_y_coordinate", "broadcast": "d"}
  double precision, allocatable, target :: xu(:,:)          !< [m] x-coord. comp. grid u-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "projection_x_coordinate", "broadcast": "d"}
  double precision, allocatable, target :: yu(:,:)          !< [m] y-coord. comp. grid u-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "projection_y_coordinate", "broadcast": "d"}
  double precision, allocatable, target :: xv(:,:)          !< [m] x-coord. comp. grid v-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "projection_x_coordinate", "broadcast": "d"}
  double precision, allocatable, target :: yv(:,:)          !< [m] y-coord. comp. grid v-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "projection_y_coordinate", "broadcast": "d"}
  double precision, allocatable, target :: dsu(:,:)         !< [m] grid distance in s-direction, centered around u-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dsv(:,:)         !< [m] grid distance in s-direction, centered around v-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dsz(:,:)         !< [m] grid distance in s-direction, centered around z-point (=eta-point) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dsc(:,:)         !< [m] grid distance in s-direction, centered around c-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dnu(:,:)         !< [m] grid distance in n-direction, centered around u-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dnv(:,:)         !< [m] grid distance in n-direction, centered around v-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dnz(:,:)         !< [m] grid distance in n-direction, centered around z-point (=eta-point) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dnc(:,:)         !< [m] grid distance in n-direction, centered around c-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dsdnui(:,:)      !< [1/m2] inverse of grid cell surface, centered around u-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dsdnvi(:,:)      !< [1/m2] inverse of grid cell surface, centered around v-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dsdnzi(:,:)      !< [1/m2] inverse of grid cell surface, centered around z-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: alfaz(:,:)       !< [rad] grid orientation at z-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "angle_of_rotation_from_east_to_x", "broadcast": "d"}
  double precision, allocatable, target :: alfau(:,:)       !< [rad] grid orientation at u-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "angle_of_rotation_from_east_to_x", "broadcast": "d"}
  double precision, allocatable, target :: alfav(:,:)       !< [rad] grid orientation at v-point {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "angle_of_rotation_from_east_to_x", "broadcast": "d"}
  double precision, allocatable, target :: sdist(:,:)       !< [m] cum. distance from offshore boundary along s-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ndist(:,:)       !< [m] cum. distance from right boundary along n-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision,              target :: dx               !< [m] grid size x-direction {"shape": [], "standard_name": "", "broadcast": "b"}
  double precision,              target :: dy               !< [m] grid size y-direction {"shape": [], "standard_name": "", "broadcast": "b"}
  double precision,              target :: xori             !< [m] x-origin of grid in world coordinate {"shape": [], "standard_name": "projection_x_coordinate", "broadcast": "b"}
  double precision,              target :: yori             !< [m] y-origin of grid in world coordinates {"shape": [], "standard_name": "projection_y_coordinate", "broadcast": "b"}
  double precision,              target :: alfa             !< [rad] (deg on input) angle of grid w.r.t. East {"shape": [], "standard_name": "angle_of_rotation_from_east_to_x", "broadcast": "b"}
  double precision,              target :: posdwn           !< [-] depths defined positive downwards (1) or upwards(-1) {"shape": [], "standard_name": "", "broadcast": "b"}
  double precision,              target :: thetamin         !< [rad] minimum angle of computational wave grid (cart. in rad) {"shape": [], "standard_name": "angle_of_rotation_from_east_to_x", "broadcast": "b"}
  double precision,              target :: thetamax         !< [rad] minimum angle of computational wave grid (cart. in rad) {"shape": [], "standard_name": "angle_of_rotation_from_east_to_x", "broadcast": "b"}
  integer,                       target :: nx               !< [-] local number of grid cells x-direction {"shape": [], "standard_name": "", "broadcast": "b"}
  integer,                       target :: ny               !< [-] local number of grid cells y-direction {"shape": [], "standard_name": "", "broadcast": "b"}
  double precision,              target :: zs01             !< [m] Initial water level first sea boundary {"shape": [], "standard_name": "sea_surface_height_above_sea_level", "broadcast": "b"}
  double precision,              target :: zs02             !< [m] Initial water level second sea boundary {"shape": [], "standard_name": "sea_surface_height_above_sea_level", "broadcast": "b"}
  double precision,              target :: zs03             !< [m] Initial water level first land boundary {"shape": [], "standard_name": "sea_surface_height_above_sea_level", "broadcast": "b"}
  double precision,              target :: zs04             !< [m] Initial water level second land boundary {"shape": [], "standard_name": "sea_surface_height_above_sea_level", "broadcast": "b"}
  double precision,              target :: xyzs01(:)        !< [-] global xy coordinates of corner (x=1,y=1) {"shape": [2], "standard_name": "", "broadcast": "b"}
  double precision,              target :: xyzs02(:)        !< [-] global xy coordinates of corner (x=1,y=N) {"shape": [2], "standard_name": "", "broadcast": "b"}
  double precision,              target :: xyzs03(:)        !< [-] global xy coordinates of corner (x=N,y=N) {"shape": [2], "standard_name": "", "broadcast": "b"}
  double precision,              target :: xyzs04(:)        !< [-] global xy coordinates of corner (x=N,y=1) {"shape": [2], "standard_name": "", "broadcast": "b"}
  integer,                       target :: tidelen          !< [-] length of tide time series {"shape": [], "standard_name": "", "broadcast": "b"}
  integer,                       target :: windlen          !< [-] length of tide time series {"shape": [], "standard_name": "", "broadcast": "b"}
  double precision, allocatable, target :: zb(:,:)          !< [m] bed level {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "altitude", "broadcast": "d"}
  double precision, allocatable, target :: zb0(:,:)         !< [m] initial bed level {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "altitude", "broadcast": "d"}
  double precision, allocatable, target :: theta(:)         !< [rad] wave angles directional distribution w.r.t. comp. x-axis {"shape": ["s%ntheta"], "standard_name": "sea_surface_wind_wave_to_direction", "broadcast": "b"}
  double precision, allocatable, target :: theta_s(:)       !< [rad] wave angles directional distribution w.r.t. comp. x-axis {"shape": ["s%ntheta_s"], "standard_name": "sea_surface_wind_wave_to_direction", "broadcast": "b"}  
  integer,                       target :: ntheta           !< [-] number of wave direction bins {"shape": [], "standard_name": "", "broadcast": "b"}
  integer,                       target :: ntheta_s         !< [-] number of wave direction bins {"shape": [], "standard_name": "", "broadcast": "b"}
  double precision,              target :: dtheta           !< [rad] wave direction bin size {"shape": [], "standard_name": "", "broadcast": "b"}
  double precision,              target :: dtheta_s         !< [rad] wave direction bin size {"shape": [], "standard_name": "", "broadcast": "b"}
  double precision,              target :: theta0           !< [rad] mean incident wave angle {"shape": [], "standard_name": "sea_surface_wind_wave_to_direction", "broadcast": "b"}
  double precision, allocatable, target :: thetamean(:,:)   !< [rad] mean wave angle {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_surface_wind_wave_to_direction", "broadcast": "d"}
  double precision, allocatable, target :: Fx(:,:)          !< [N/m2] wave force, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Fy(:,:)          !< [N/m2] wave force, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Sxy(:,:)         !< [N/m] radiation stress, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Syy(:,:)         !< [N/m] radiation stress, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Sxx(:,:)         !< [N/m] radiation stress, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: n(:,:)           !< [-] ratio group velocity/wave celerity {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: H(:,:)           !< [m] Hrms wave height based on instantaneous wave energy {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cgx(:,:,:)       !< [m/s] group velocity, x-component {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cgy(:,:,:)       !< [m/s] group velocity, y-component {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cx(:,:,:)        !< [m/s] wave celerity, x-component {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cy(:,:,:)        !< [m/s] wave celerity, y-component {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ctheta(:,:,:)    !< [rad/s] wave celerity theta-direction (refraction) {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ee(:,:,:)        !< [J/m2/rad] directionally distributed wave energy {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: thet(:,:,:)      !< [rad] wave angles {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "sea_surface_wind_wave_to_direction", "broadcast": "d"}
  double precision, allocatable, target :: costh(:,:,:)     !< [-] cos of wave angles relative to grid direction {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: sinth(:,:,:)     !< [-] sin of wave angles relative to grid direction {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: sigt(:,:,:)      !< [rad/s] relative frequency {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: rr(:,:,:)        !< [J/m2/rad] directionally distributed roller energy {"shape": ["s%nx+1", "s%ny+1", "s%ntheta"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cgx_s(:,:,:)     !< [m/s] group velocity, x-component {"shape": ["s%nx+1", "s%ny+1", "s%ntheta_s"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cgy_s(:,:,:)     !< [m/s] group velocity, y-component {"shape": ["s%nx+1", "s%ny+1", "s%ntheta_s"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ctheta_s(:,:,:)  !< [rad/s] wave celerity theta-direction (refraction) {"shape": ["s%nx+1", "s%ny+1", "s%ntheta_s"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ee_s(:,:,:)      !< [J/m2/rad] directionally distributed wave energy {"shape": ["s%nx+1", "s%ny+1", "s%ntheta_s"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: thet_s(:,:,:)    !< [rad] wave angles {"shape": ["s%nx+1", "s%ny+1", "s%ntheta_s"], "standard_name": "sea_surface_wind_wave_to_direction", "broadcast": "d"}
  double precision, allocatable, target :: costh_s(:,:,:)   !< [-] cos of wave angles relative to grid direction {"shape": ["s%nx+1", "s%ny+1", "s%ntheta_s"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: sinth_s(:,:,:)   !< [-] sin of wave angles relative to grid direction {"shape": ["s%nx+1", "s%ny+1", "s%ntheta_s"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: k(:,:)           !< [rad/m] wave number {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: c(:,:)           !< [m/s] wave celerity {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cg(:,:)          !< [m/s] group velocity {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: sigm(:,:)        !< [rad/s] mean frequency {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: wm(:,:)          !< [rad/s] mean abs frequency {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: hh(:,:)          !< [m] water depth {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: zs(:,:)          !< [m] water level {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_surface_height_above_sea_level", "broadcast": "d"}
  double precision, allocatable, target :: zs0(:,:)         !< [m] water level due to tide alone {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_surface_height_above_sea_level", "broadcast": "d"}
  double precision, allocatable, target :: tideinpt(:)      !< [s] input time of input tidal signal {"shape": ["s%tidelen"], "standard_name": "time", "broadcast": "b"}
  double precision, allocatable, target :: tideinpz(:,:)    !< [m] input tidal signal {"shape": ["s%tidelen", "par%tideloc"], "standard_name": "sea_surface_height_above_sea_level", "broadcast": "b"}
  double precision, allocatable, target :: windinpt(:)      !< [s] input time of input wind signal {"shape": ["s%windlen"], "standard_name": "time", "broadcast": "b"}
  double precision, allocatable, target :: windvelts(:)     !< [m/s] input wind velocity {"shape": ["s%windlen"], "standard_name": "", "broadcast": "b"}
  double precision, allocatable, target :: winddirts(:)     !< [deg_nautical] input wind direction {"shape": ["s%windlen"], "standard_name": "", "broadcast": "b"}
  double precision, allocatable, target :: windxts(:)       !< [m/s] time series of input wind velocity (not S direction), x-component {"shape": ["s%windlen"], "standard_name": "", "broadcast": "b"}
  double precision, allocatable, target :: windyts(:)       !< [m/s] time series of input wind velocity (not N direction), y-component {"shape": ["s%windlen"], "standard_name": "", "broadcast": "b"}
  double precision, allocatable, target :: windsu(:,:)      !< [m/s] wind velocity in S direction in u point at current time step {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: windnv(:,:)      !< [m/s] wind velocity in N direction in v point at current time step {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dzsdt(:,:)       !< [m/s] rate of change water level {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dzsdx(:,:)       !< [m/s] water surface gradient in x-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dzsdy(:,:)       !< [m/s] water surface gradient in y-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dzbdx(:,:)       !< [-] bed level gradient in x-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dzbdy(:,:)       !< [-] bed level gradient in y-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dzbdt(:,:)       !< [m/s] rate of change bed level {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dzbnow(:,:)      !< [m] bed level change in current time step {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: uu(:,:)          !< [m/s] GLM velocity in u-points, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: vu(:,:)          !< [m/s] GLM velocity in u-points, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_y_velocity", "broadcast": "d"}
  double precision, allocatable, target :: uv(:,:)          !< [m/s] GLM velocity in v-points, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: vv(:,:)          !< [m/s] GLM velocity in v-points, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_y_velocity", "broadcast": "d"}
  double precision, allocatable, target :: qx(:,:)          !< [m2/s] discharge in u-points, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: qy(:,:)          !< [m2/s] discharge in u-points, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: sedero(:,:)      !< [m] cum. sedimentation/erosion {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dcbdx(:,:)       !< [kg/m3/m] bed concentration gradient x-dir. {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dcbdy(:,:)       !< [kg/m3/m] bed concentration gradient y-dir. {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dcsdx(:,:)       !< [kg/m3/m] suspended concentration gradient x-dir. {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dcsdy(:,:)       !< [kg/m3/m] suspended concentration gradient y-dir. {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ui(:,:)          !< [m/s] incident bound wave velocity in, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vi(:,:)          !< [m/s] incident bound wave velocity in, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: E(:,:)           !< [Nm/m2] wave energy {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: R(:,:)           !< [Nm/m2] roller energy {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: urms(:,:)        !< [m/s] orbital velocity {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: D(:,:)           !< [W/m2] dissipation {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Qb(:,:)          !< [-] fraction breaking waves {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ust(:,:)         !< [m/s] Stokes drift {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ueu(:,:)         !< [m/s] Eulerian velocity in u-points, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: vev(:,:)         !< [m/s] Eulerian velocity in u-points, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_y_velocity", "broadcast": "d"}
  double precision, allocatable, target :: vmagu(:,:)       !< [m/s] GLM velocity magnitude u-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vmageu(:,:)      !< [m/s] Eulerian velocity magnitude u-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vmagv(:,:)       !< [m/s] GLM velocity magnitude v-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vmagev(:,:)      !< [m/s] Eulerian velocity magnitude v-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: u(:,:)           !< [m/s] GLM velocity in cell centre, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: v(:,:)           !< [m/s] GLM velocity in cell centre, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_y_velocity", "broadcast": "d"}
  double precision, allocatable, target :: ue(:,:)          !< [m/s] Eulerian velocity in cell centre, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: ve(:,:)          !< [m/s] Eulerian velocity in cell centre, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_y_velocity", "broadcast": "d"}
  double precision, allocatable, target :: ue_sed(:,:)      !< [m/s] Advection velocity sediment in cell centre, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: ve_sed(:,:)      !< [m/s] Advection velocity sediment in cell centre, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: hold(:,:)        !< [m] water depth previous time step {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  integer,          allocatable, target :: wetu(:,:)        !< [-] mask wet/dry u-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  integer,          allocatable, target :: wetv(:,:)        !< [-] mask wet/dry v-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  integer,          allocatable, target :: wetz(:,:)        !< [-] mask wet/dry eta-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  integer,          allocatable, target :: wete(:,:)        !< [-] mask wet/dry wave-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: hu(:,:)          !< [m] water depth in u-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: hv(:,:)          !< [m] water depth in v-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: hum(:,:)         !< [m] water depth in u-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: hvm(:,:)         !< [m] water depth in v-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vmag(:,:)        !< [m/s] velocity magnitude in cell centre {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ccg(:,:,:)       !< [m3/m3] depth-averaged suspended concentration for each sediment fraction {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: uwf(:,:)         !< [m/s] Stokes drift, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vwf(:,:)         !< [m/s] Stokes drift, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ustr(:,:)        !< [m/s] return flow due to roller {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: usd(:,:)         !< [m/s] return flow due to roller after breaker delay {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision,              target :: bi(:)            !< [m] incoming bound long wave {"shape": ["s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: DR(:,:)          !< [W/m2] roller energy dissipation {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: umean(:,:)       !< [m/s] longterm mean velocity at bnds in u-points, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: vmean(:,:)       !< [m/s] longterm mean velocity at bnds in u-points, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_y_velocity", "broadcast": "d"}
  double precision, allocatable, target :: ur(:,:)          !< [m/s] reflected velocity at bnds in u-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  integer,                       target :: vardx            !< [-] 0 = uniform grid size, 1 = variable grid size {"shape": [], "standard_name": "", "broadcast": "b"}
  double precision,              target :: D15(:)           !< [m] D15 grain diameters for all sediment classses {"shape": ["par%ngd"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: D50(:)           !< [m] D50 grain diameters for all sediment classses {"shape": ["par%ngd"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: D90(:)           !< [m] D90 grain diameters for all sediment classses {"shape": ["par%ngd"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: sedcal(:)        !< [-] equilibrium sediment concentartion factor for each sediment class {"shape": ["par%ngd"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: ucrcal(:)        !< [-] calibration factor for u critical for each sediment class {"shape": ["par%ngd"], "standard_name": "", "broadcast": "b"}
  double precision, allocatable, target :: Tsg(:,:,:)       !< [s] sediment response time for each sediment class {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Susg(:,:,:)      !< [m2/s] suspended sediment transport for each sediment class (excluding pores), x-component {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Svsg(:,:,:)      !< [m2/s] suspended sediment transport for each sediment class (excluding pores), y-component {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Subg(:,:,:)      !< [m2/s] bed sediment transport for each sediment class (excluding pores), x-component {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Svbg(:,:,:)      !< [m2/s] bed sediment transport for each sediment class (excluding pores), y-component {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ceqbg(:,:,:)     !< [m3/m3] depth-averaged bed equilibrium concentration for each sediment class {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ceqsg(:,:,:)     !< [m3/m3] depth-averaged suspended equilibrium concentration for each sediment class {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ua(:,:)          !< [m/s] time averaged flow velocity due to wave assymetry {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: BR(:,:)          !< [-] maximum wave surface slope used in roller dissipation formulation {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: kb(:,:)          !< [m2/s2] near bed turbulence intensity due to depth induces breaking {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Tbore(:,:)       !< [s] wave period interval associated with breaking induced turbulence {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dzav(:,:)        !< [m] total bed level change due to avalanching {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: maxzs(:,:)       !< [m] maximum elevation in simulation {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: minzs(:,:)       !< [m] minimum elevation in simulation {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: L1(:,:)          !< [m] wave length (used in dispersion relation) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Sk(:,:)          !< [-] skewness of short waves {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: As(:,:)          !< [-] asymmetry of short waves {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gwhead(:,:)      !< [m] groundwater head (differs from gwlevel) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gwheadb(:,:)     !< [m] groundwater head at bottom (differs from gwlevel) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gwlevel(:,:)     !< [m] groundwater table (min(zb,gwhead)) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gwheight(:,:)    !< [m] vertical size of aquifer through which groundwater can flow {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gwbottom(:,:)    !< [m] level of the bottom of the aquifer {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gwu(:,:)         !< [m/s] groundwater flow in x-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gwv(:,:)         !< [m/s] groundwater flow in y-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gwqx(:,:)        !< [m/s] groundwater discharge in x-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gwqy(:,:)        !< [m/s] groundwater discharge in y-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gww(:,:)         !< [m/s] groundwater flow in z-direction (interaction between surface and ground water) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gwcurv(:,:)      !< [-] curvature coefficient of groundwater head function {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dinfil(:,:)      !< [m] Infiltration layer depth used in quasi-vertical flow model for groundwater {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: infil(:,:)       !< [m/s] Rate of exchange of water between surface and groundwater (positive from sea to groundwater) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: gw0back(:,:)     !< [m] boundary condition back boundary for groundwater head {"shape": [2, "s%ny+1"], "standard_name": "", "broadcast": "2"}
  double precision, allocatable, target :: Kx(:,:)          !< [m/s] (Turbulent) Hydraulic conductivity in x-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Ky(:,:)          !< [m/s] (Turbulent) Hydraulic conductivity in y-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Kz(:,:)          !< [m/s] (Turbulent) Hydraulic conductivity in z-direction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Kzinf(:,:)       !< [m/s] (Turbulent) Hydraulic conductivity in z-direction for infiltration {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: kturb(:,:)       !< [m2/s2] depth averaged turbulence intensity due to long wave breaking {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ero(:,:,:)       !< [m/s] bed erosion rate per fraction {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: depo_im(:,:,:)   !< [m/s] implicit bed deposition rate per fraction {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: depo_ex(:,:,:)   !< [m/s] explicit bed deposition rate per fraction {"shape": ["s%nx+1", "s%ny+1", "par%ngd"], "standard_name": "", "broadcast": "d"}
  integer,          allocatable, target :: nd(:,:)          !< [-] number of bed layers (can be different for each computational cell) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision,              target :: pbbed(:,:,:,:)   !< [-] NO DESCRIPTION {"shape": ["s%nx+1", "s%ny+1", "par%nd", "par%ngd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dzbed(:,:,:)     !< [-] NO DESCRIPTION {"shape": ["s%nx+1", "s%ny+1", "par%nd"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: z0bed(:,:)       !< [-] NO DESCRIPTION {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ureps(:,:)       !< [m/s] representative flow velocity for sediment advection and diffusion, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: vreps(:,:)       !< [m/s] representative flow velocity for sediment advection and diffusion, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_y_velocity", "broadcast": "d"}
  double precision, allocatable, target :: urepb(:,:)       !< [m/s] representative flow velocity for sediment advection and diffusion, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: vrepb(:,:)       !< [m/s] representative flow velocity for sediment advection and diffusion, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_y_velocity", "broadcast": "d"}
  double precision, allocatable, target :: umwci(:,:)       !< [m/s] velocity (time-averaged) for wci, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_x_velocity", "broadcast": "d"}
  double precision, allocatable, target :: vmwci(:,:)       !< [m/s] velocity (time-averaged) for wci, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_water_y_velocity", "broadcast": "d"}
  double precision, allocatable, target :: rolthick(:,:)    !< [m] long wave roller thickness {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: zswci(:,:)       !< [m] waterlevel (time-averaged) for wci {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: pres(:,:)        !< [m2/s2] normalized dynamic pressure {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dU(:,:)          !< [m2/s2] u-velocity difference between two vertical layers (reduced 2-layer non-hydrostatic model) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dV(:,:)          !< [m2/s2] v-velocity difference between two vertical layers (reduced 2-layer non-hydrostatic model) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: wb(:,:)          !< [m/s] vertical velocity at the bottom {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ws(:,:)          !< [m/s] vertical velocity at the free surface {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: wscrit(:,:)      !< [m/s] critial vertical velocity at the free surface for breaking {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: bedfriccoef(:,:) !< [-] dimensional/dimensionless input bed friction coefficient; depends on value of parameter bedfriction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: taubx(:,:)       !< [N/m2] bed shear stress, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: tauby(:,:)       !< [N/m2] bed shear stress, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Df(:,:)          !< [W/m2] dissipation rate due to bed friction {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Dp(:,:)          !< [W/m2] dissipation rate in the swash due to transformation of kinetic wave energy to potential wave energy {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Sutot(:,:)       !< [m2/s] Sediment transport integrated over bed load and suspended and for all sediment grains, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Svtot(:,:)       !< [m2/s] Sediment transport integrated over bed load and suspended and for all sediment grains, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cctot(:,:)       !< [m3/m3] Sediment concentration integrated over bed load and suspended and for all sediment grains {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: wi(:,:)          !< [m/s] Vertical velocity at boundary due to (short) waves {"shape": [2, "s%ny+1"], "standard_name": "", "broadcast": "2"}
  double precision, allocatable, target :: dUi(:,:)         !< [m/s] Velocity difference at boundary due to (short) waves {"shape": [2, "s%ny+1"], "standard_name": "", "broadcast": "2"}
  double precision, allocatable, target :: zi(:,:)          !< [m] Surface elevation at boundary due to (short) waves {"shape": [2, "s%ny+1"], "standard_name": "", "broadcast": "2"}
  double precision, allocatable, target :: nuh(:,:)         !< [m2/s] horizontal viscosity coefficient {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cf(:,:)          !< [-] Friction coefficient flow {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cfu(:,:)         !< [-] Friction coefficient flow in u-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cfv(:,:)         !< [-] Friction coefficient flow in v-points {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: D50top(:,:)      !< [-] Friction coefficient flow {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: D90top(:,:)      !< [-] Friction coefficient flow {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: structdepth(:,:) !< [m] Depth of structure in relation to instantaneous bed level {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: zs0fac(:,:,:)    !< [-] relative weight of offshore boundary and bay boundary for each grid point is stored in zs0fac {"shape": ["s%nx+1", "s%ny+1", 2], "standard_name": "", "broadcast": "d"}
  double precision,              target :: tdisch(:)        !< [-] Discharge time series {"shape": ["par%ntdischarge"], "standard_name": "", "broadcast": "b"}
  double precision, allocatable, target :: pdisch(:,:)      !< [-] Discharge locations {"shape": ["par%ndischarge", 4], "standard_name": "", "broadcast": "b"}
  integer,          allocatable, target :: pntdisch(:)      !< [-] Point discharge locations (no momentum) {"shape": ["par%ndischarge"], "standard_name": "", "broadcast": "b"}
  double precision, allocatable, target :: qdisch(:,:)      !< [m2/s] Discharges {"shape": ["par%ntdischarge", "par%ndischarge"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: idrift(:)        !< [-] Drifter x-coordinate in grid space {"shape": ["par%ndrifter"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: jdrift(:)        !< [-] Drifter y-coordinate in grid space {"shape": ["par%ndrifter"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: tdriftb(:)       !< [s] Drifter release time {"shape": ["par%ndrifter"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: tdrifte(:)       !< [s] Drifter retrieval time {"shape": ["par%ndrifter"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: runup(:)         !< [m] Short wave runup height {"shape": ["s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision,              target :: Hrunup(:)        !< [m] Short wave height used in runup formulation {"shape": ["s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision,              target :: xHrunup(:)       !< [m] Location at which short wave height for runup is taken {"shape": ["s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision,              target :: istruct(:)       !< [-] location of revetments toe {"shape": ["s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision,              target :: iwl(:)           !< [-] location of water line (including long wave runup) {"shape": ["s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision,              target :: strucslope(:)    !< [-] slope of structure {"shape": ["s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Dc(:,:)          !< [m2/s] diffusion coefficient {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ph(:,:)          !< [m] pressure head due to ship {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  integer,                       target :: newstatbc        !< [-] (1) Use new stationary boundary conditions for instat is stat or stat_table {"shape": [], "standard_name": "", "broadcast": "b"}
  double precision, allocatable, target :: dobs(:,:)        !< [W/m2] beachwizard observed dissipation {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: sig2prior(:,:)   !< [m2] beachwizard prior std squared {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: zbobs(:,:)       !< [m] beachwizard observed depth {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: shobs(:,:)       !< [m] beachwizard observed shoreline {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: bwalpha(:,:)     !< [-] beachwizard weighting factor {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dcmdo(:,:)       !< [W/m2] beachwizard computed minus observed dissipation {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dassim(:,:)      !< [m] beachwizard depth change {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: cobs(:,:)        !< [m/s] beachwizard observed wave celerity {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision,              target :: shipxCG(:)       !< [m] x-coordinate of ship center of gravity {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shipyCG(:)       !< [m] y-coordinate of ship center of gravity {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shipzCG(:)       !< [m] z-coordinate of ship center of gravity {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shipFx(:)        !< [N] force on ship in x-direction {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shipFy(:)        !< [N] force on ship in y-direction {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shipFz(:)        !< [N] force on ship in z-direction {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shipMx(:)        !< [Nm] moment on ship around x-axis {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shipMy(:)        !< [Nm] moment on ship around y-axis {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shipMz(:)        !< [Nm] moment on ship around z-axis {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shipphi(:)       !< [deg] turning angle arround x-axis {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shipchi(:)       !< [deg] turning angle arround y-axis {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  double precision,              target :: shippsi(:)       !< [deg] turning angle arround z-axis {"shape": ["par%nship"], "standard_name": "", "broadcast": "b"}
  integer,          allocatable, target :: vegtype(:,:)     !< [-] vegetation type index {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Cdrag(:,:)       !< [-] Vegetation drag coefficient {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Dveg(:,:)        !< [W/m2] dissipation due to short wave attenuation by vegetation {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Fvegu(:,:)       !< [N/m2] x-forcing due to long wave attenuation by vegetation {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: Fvegv(:,:)       !< [N/m2] y-forcing due to long wave attenuation by vegetation {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ududx(:,:)       !< [m2/s2] advection {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vdvdy(:,:)       !< [m2/s2] advection {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: udvdx(:,:)       !< [m2/s2] advection {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vdudy(:,:)       !< [m2/s2] advection {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: viscu(:,:)       !< [m2/s2] viscosity {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: viscv(:,:)       !< [m2/s2] viscosity {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: setbathy(:,:,:)  !< [m] prescribed bed levels {"shape": ["s%nx+1", "s%ny+1", "par%nsetbathy"], "standard_name": "", "broadcast": "d"}
  double precision,              target :: tsetbathy(:)     !< [s] points in time of prescibed bed levels {"shape": ["par%nsetbathy"], "standard_name": "", "broadcast": "b"}
  integer,          allocatable, target :: breaking(:,:)    !< [-] indicator whether cell has breaking nonh waves {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: fw(:,:)          !< [-] wave friction coefficient {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: sigz(:)          !< [-] vertical distribution of sigma layers Q3D {"shape": ["par%nz"], "standard_name": "", "broadcast": "b"}
  double precision, allocatable, target :: uz(:,:,:)        !< [-] velocity (Q3D) ksi-comp {"shape": ["s%nx+1","s%ny+1","par%nz"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vz(:,:,:)        !< [-] velocity (Q3D) eta-comp {"shape": ["s%nx+1","s%ny+1","par%nz"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ustz(:,:,:)      !< [-] stokes velocity (Q3D) {"shape": ["s%nx+1","s%ny+1","par%nz"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: nutz(:,:,:)      !< [-] turbulence viscosity {"shape": ["s%nx+1","s%ny+1","par%nz"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: dzs0dn(:,:)      !< [-] alongshore water level gradient due to tide alone {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ccz(:,:,:)       !< [m3/m3] concentration profile {"shape": ["s%nx+1","s%ny+1","par%nz"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: refA(:,:)        !< [m] reference level {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: ca(:,:)          !< [m3/m3] reference concentration {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: zs1(:,:)         !< [m] water level minus tide {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "sea_surface_height_above_sea_level_without_tide", "broadcast": "d"}
  double precision, allocatable, target :: taubx_add(:,:)   !< [N/m2] additional bed shear stress due to boundary layer effects, x-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: tauby_add(:,:)   !< [N/m2] additional bed shear stress due to boundary layer effects, y-component {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: hhw(:,:)         !< [m] water depth used in all wave computations, includes H*par%delta {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: hhws(:,:)        !< [m] water depth used in wave stationary computation (and single_dir wave directions) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: uws(:,:)         !< [m/s] u-velocity used in wave stationary computation (and single_dir wave directions) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vws(:,:)         !< [m/s] v-velocity used in wave stationary computation (and single_dir wave directions) {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: hhwcins(:,:)     !< [m] water depth used in wave instationary computation in case of wci {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: uwcins(:,:)      !< [m/s] u-velocity used in wave stationary computation in case of wci {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
  double precision, allocatable, target :: vwcins(:,:)      !< [m/s] v-velocity used in wave stationary computation in case of wci {"shape": ["s%nx+1", "s%ny+1"], "standard_name": "", "broadcast": "d"}
end module variables
