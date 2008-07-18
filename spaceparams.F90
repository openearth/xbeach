module spaceparams
type spacepars
real*8,dimension(:,:),pointer   :: x   => NULL()      ! [m]        x-coord. comp. grid (positive shoreward, perp. to coastline)
real*8,dimension(:,:),pointer   :: y   =>NULL()      ! [m]        y-coord. comp. grid
real*8,dimension(:),pointer     :: xz  =>NULL()      ! [m]        x-coord. comp. grid (positive shoreward, perp. to coastline)
real*8,dimension(:),pointer     :: yz  =>NULL()      ! [m]        y-coord. comp. grid
real*8,dimension(:),pointer     :: xu  =>NULL()      ! [m]        x-coord. in u points
real*8,dimension(:),pointer     :: yv  =>NULL()      ! [m]        y-coord. in v points
real*8,dimension(:,:),pointer   :: xw  =>NULL()      ! [m]        world x-coordinates
real*8,dimension(:,:),pointer   :: yw  =>NULL()      ! [m]        world y-coordinates
real*8                          :: dx        !`[m]        grid size x-direction
real*8                          :: dy        ! [m]        grid size y-direction
real*8                          :: xori      ! [m]        x-origin of grid in world coordinates
real*8                          :: yori      ! [m]        y-origin of grid in world coordinates
real*8                          :: alfa      ! [rad]      (deg on input) angle of grid w.r.t. East
real*8                          :: posdwn    ! [-]        depths defined positive downwards (1) or upwards(-1)
integer                         :: nx        ! [-]        local number of grid cells x-direction
integer                         :: ny        ! [-]        local number of grid cells y-direction
real*8,dimension(:,:),pointer   :: zb      =>NULL()      ! [m]        bed level
real*8,dimension(:,:),pointer   :: zb0     =>NULL()      ! [m]        initial bed level
real*8,dimension(:),pointer     :: theta   =>NULL()    ! [rad]      wave angles directional distribution
                                                 !            w.r.t. comp. x-axis
integer                         :: ntheta    ! [-]        number of wave direction bins
real*8                          :: dtheta    ! [rad]      wave direction bin size
real*8                          :: theta0    ! [rad]      mean incident wave angle
real*8,dimension(:),pointer     :: cxsth    =>NULL()   ! [-]        cos(theta)
real*8,dimension(:),pointer     :: sxnth    =>NULL()   ! [-]        sin(theta)
real*8,dimension(:,:),pointer   :: thetamean   =>NULL()! [rad]      mean wave angle
real*8,dimension(:,:),pointer   :: Fx     =>NULL()     ! [N/m2]     wave force x-direction
real*8,dimension(:,:),pointer   :: Fy     =>NULL()     ! [N/m2]     wave force y-direction
real*8,dimension(:,:),pointer   :: Sxy    =>NULL()     ! [N/m]      radiation stress
real*8,dimension(:,:),pointer   :: Syy    =>NULL()     ! [N/m]      radiation stress
real*8,dimension(:,:),pointer   :: Sxx     =>NULL()    ! [N/m]      radiation stress
real*8,dimension(:,:),pointer   :: n       =>NULL()    ! [-]        ratio group velocity/wave celerity
real*8,dimension(:,:),pointer   :: H      =>NULL()     ! [m]        wave height 
real*8,dimension(:,:,:),pointer :: cgx   =>NULL()    ! [m/s]      group velocity x-direction
real*8,dimension(:,:,:),pointer :: cgy   =>NULL()    ! [m/s]      group velocity y-direction
real*8,dimension(:,:,:),pointer :: cx   =>NULL()     ! [m/s]      wave celerity x-direction
real*8,dimension(:,:,:),pointer :: cy   =>NULL()     ! [m/s]      wave celerity y-direction
real*8,dimension(:,:,:),pointer :: ctheta =>NULL()   ! [rad/s]    wave celerity theta-direction (refraction)
real*8,dimension(:,:,:),pointer :: ee     =>NULL()   ! [J/m2/rad] directionally distributed wave energy
real*8,dimension(:,:,:),pointer :: thet   =>NULL()   ! [rad]      wave angles  
real*8,dimension(:,:,:),pointer :: costhet =>NULL()  ! [-]        cos of wave angles  
real*8,dimension(:,:,:),pointer :: sinthet =>NULL()  ! [-]        sin of wave angles  
real*8,dimension(:,:,:),pointer :: sigt   =>NULL()   ! [rad/s]    relative frequency 
real*8,dimension(:,:,:),pointer :: rr    =>NULL()    ! [J/m2/rad] directionally distributed roller energy
real*8,dimension(:,:),pointer   :: k     =>NULL()    ! [rad/m]    wave number
real*8,dimension(:,:),pointer   :: c     =>NULL()    ! [m/s]      wave celerity
real*8,dimension(:,:),pointer   :: cg    =>NULL()    ! [m/s]      group velocity
real*8,dimension(:,:),pointer   :: sigm  =>NULL()    ! [rad/s]    mean frequency
real*8,dimension(:,:),pointer   :: hh    =>NULL()    ! [m]        water depth
real*8,dimension(:,:),pointer   :: zs    =>NULL()    ! [m]        water level
real*8,dimension(:,:),pointer   :: zs0   =>NULL()    ! [m]        water level due to tide alone
real*8,dimension(:),pointer     :: tideinpt =>NULL() ! [s]        input time of input tidal signal
real*8,dimension(:,:),pointer   :: tideinpz=>NULL()  ! [m]        input tidal signal
real*8,dimension(:,:),pointer   :: dzsdt   =>NULL()  ! [m/s]      rate of change water level
real*8,dimension(:,:),pointer   :: dzbdt   =>NULL()  ! [m/s]      rate of change bed level 
real*8,dimension(:,:),pointer   :: uu     =>NULL()   ! [m/s]      (GLM) x-velocity in u-points
real*8,dimension(:,:),pointer   :: vv     =>NULL()   ! [m/s]      (GLM) y-velocity in v-points
real*8,dimension(:,:),pointer   :: qx     =>NULL()   ! [m2/s]     x-discharge in u-points
real*8,dimension(:,:),pointer   :: qy     =>NULL()   ! [m2/s]     y-discharge in u-points
real*8,dimension(:,:),pointer   :: sedero  =>NULL()  ! [m]        cum. sedimentation/erosion
real*8,dimension(:,:),pointer   :: dcdx   =>NULL()   ! [kg/m3/m]  concentration gradient x-dir.
real*8,dimension(:,:),pointer   :: dcdy   =>NULL()   ! [kg/m3/m]  concentration gradient y-dir.
real*8,dimension(:,:),pointer   :: ui     =>NULL()   ! [m/s]      incident bound wave velocity
real*8,dimension(:,:),pointer   :: E      =>NULL()   ! [Nm/m2]    wave energy
real*8,dimension(:,:),pointer   :: R      =>NULL()   ! [Nm/m2]    roller energy
real*8,dimension(:,:),pointer   :: urms   =>NULL()   ! [m/s]      orbital velocity
real*8,dimension(:,:),pointer   :: D      =>NULL()   ! [W/m2]     dissipation
real*8,dimension(:,:),pointer   :: ust     =>NULL()  ! [m/s]      Stokes drift
real*8,dimension(:,:),pointer   :: tm      =>NULL()  ! [rad]      mean wave direction
real*8,dimension(:,:),pointer   :: ueu     =>NULL()  ! [m/s]      Eulerian mean velocity x-dir.
real*8,dimension(:,:),pointer   :: vev    =>NULL()   ! [m/s]      Eulerian mean velocity y-dir.
real*8,dimension(:,:),pointer   :: vmagu  =>NULL()   ! [m/s]      (GLM) velocity magnitude u-points
real*8,dimension(:,:),pointer   :: vmageu =>NULL()   ! [m/s]      (GLM) velocity magnitude u-points
real*8,dimension(:,:),pointer   :: vmagv  =>NULL()   ! [m/s]      (GLM) velocity magnitude v-points
real*8,dimension(:,:),pointer   :: vmagev =>NULL()   ! [m/s]      (GLM) velocity magnitude v-points
real*8,dimension(:,:),pointer   :: u      =>NULL()   ! [m/s]      (GLM) x-velocity cell centre (for output)
real*8,dimension(:,:),pointer   :: v      =>NULL()   ! [m/s]      (GLM) y-velocity cell centre (for output)
real*8,dimension(:,:),pointer   :: ue     =>NULL()   ! [m/s]      Eulerian mean x-velocity cell centre (for output)
real*8,dimension(:,:),pointer   :: ve     =>NULL()   ! [m/s]      Eulerian mean y-velocity cell centre (for output)
real*8,dimension(:,:),pointer   :: hold   =>NULL()   ! [m]        water depth previous time step
integer,dimension(:,:),pointer  :: wetu   =>NULL()   ! [-]        mask wet/dry u-points
integer,dimension(:,:),pointer  :: wetv   =>NULL()   ! [-]        mask wet/dry v-points
integer,dimension(:,:),pointer  :: wetz   =>NULL()   ! [-]        mask wet/dry eta-points
real*8,dimension(:,:),pointer   :: hu     =>NULL()   ! [m]        water depth in u-points
real*8,dimension(:,:),pointer   :: hv     =>NULL()   ! [m]        water depth in v-points
real*8,dimension(:,:),pointer   :: hum    =>NULL()   ! [m]        water depth in u-points
real*8,dimension(:,:),pointer   :: hvm     =>NULL()  ! [m]        water depth in v-points
!real*8,dimension(:,:),pointer  :: ceq     =>NULL()  ! [m3/m3]    depth-averaged equilibrium concentration
real*8,dimension(:,:),pointer   :: vmag    =>NULL()  ! [m/s]      velocity magnitude in cell centre
!real*8,dimension(:,:),pointer   :: Su     =>NULL()   ! [m2/s]     sediment transport x-dir. (excluding pores)
!real*8,dimension(:,:),pointer   :: Sv     =>NULL()   ! [m2/s]     sediment transport y-dir. (excluding pores)
!real*8,dimension(:,:),pointer   :: Ts     =>NULL()   ! [s]        adaptation time scale
!real*8,dimension(:,:),pointer   :: cc     =>NULL()   ! [m3/m3]    depth-averaged concentration
real*8,dimension(:,:,:),pointer :: ccg     =>NULL()  ! [m3/m3]    depth-averaged concentration for each sediment fraction
real*8,dimension(:,:),pointer   :: uwf     =>NULL()  ! [m/s]      x-comp. Stokes drift
real*8,dimension(:,:),pointer   :: vwf     =>NULL()  ! [m/s]      y-comp. Stokes drift
real*8,dimension(:,:),pointer   :: ustr    =>NULL()  ! [m/s]      return flow due to roller
real*8,dimension(:,:),pointer   :: usd     =>NULL()  ! [m/s]      return flow due to roller after breaker delay
real*8,dimension(:),pointer     :: bi     =>NULL()   ! [m]        incoming bound long wave
real*8,dimension(:,:),pointer   :: DR     =>NULL()   ! [W/m2]     roller energy dissipation
real*8,dimension(:,:),pointer   :: umean  =>NULL()   ! [m/s]      longterm mean velocity at bnds
integer                         :: vardx     ! [-]        0 = uniform grid size, 1 = variable grid size
real*8,dimension(:,:),pointer   :: vu     =>NULL()   ! [m/s]      y velocity in u points          
real*8,dimension(:,:),pointer   :: uv            ! [m/s]          x velocity in v points
real*8,dimension(:,:,:,:),pointer :: graindistr ! [-]     fractional graindistribution for sediment classes
real*8,dimension(:),pointer      :: D50      ! [m]        D50 grain diameters for all sediment classses
real*8,dimension(:),pointer      :: D90      ! [m]        D90 grain diameters for all sediment classses
real*8,dimension(:),pointer      :: sedcal   ! [-]        equilibrium sediment concentartion factor for each sediment class
real*8,dimension(:,:,:),pointer  :: Tsg      ! [s]        sediment response time for each sediment class
real*8,dimension(:,:,:),pointer  :: Sug      ! [m2/s]     sediment transport x-dir. for each sediment class (excluding pores)
real*8,dimension(:,:,:),pointer  :: Svg      ! [m2/s]     sediment transport y-dir. for each sediment class (excluding pores)
real*8,dimension(:,:,:),pointer  :: ceqg     ! [m3/m3]    depth-averaged equilibrium concentration for each sediment class
real*8,dimension(:,:),pointer    :: ua       ! [m/s]      time averaged flow velocity due to wave assymetry
real*8,dimension(:,:),pointer    :: BR       ! [-]        maximum wave surface slope used in roller dissipation formulation
real*8,dimension(:,:),pointer    :: kb       ! [m^2/s^2]  near bed turbulence intensity due to depth induces breaking
real*8,dimension(:,:),pointer    :: Tbore    ! [s]        wave period interval associated with breaking induced turbulence
real*8,dimension(:,:),pointer    :: uon      ! [m/s]      onshore directed peak orbital velocity
real*8,dimension(:,:),pointer    :: uoff     ! [m/s]      offshore directed peak orbital velocity
real*8,dimension(:,:),pointer    :: dzav     ! [m]        total bed level change due to avalanching
real*8,dimension(:,:),pointer    :: maxzs    ! [m]        maximum elevation in simulation  
real*8,dimension(:,:),pointer    :: minzs    ! [m]        minimum elevation in simulation 
#ifdef USEMPI
!
! administration of the lay-out of the distributed matrices
! In the following: p is the MPI-process number: p = 0 .. numprocs-1
! The array indices start with 1, so, for example, is(2) describes
! the situation for process rank 1.
! a is the global matrix, b is a local matrix on process p.
!
! is(:) and js(:) describe the location of the submatrix in process p:
!    b(1,1) coincides with a(is(p+1),js(p+1))
! 
! lm(:) and ln(:) describe the extend of b:
!    the dimensions of b on p are (lm(p+1),ln(p+1))
!    b concides with 
!       a(is(p+1):is(p+1)+lm(p+1)-1,js(p+1):js(p+1)+ln(p+1)-1)
!
! isleft(:), isright(:), istop(:), isbot(:) tell if 
!    matrix b respectively map from a:
!    the first column
!    the last  column
!    the first row
!    lhe last  row
! 
integer, dimension(:), pointer  :: is      => NULL()   
integer, dimension(:), pointer  :: js      => NULL()
integer, dimension(:), pointer  :: lm      => NULL()
integer, dimension(:), pointer  :: ln      => NULL()
logical, dimension(:), pointer  :: isleft  => NULL()
logical, dimension(:), pointer  :: isright => NULL()
logical, dimension(:), pointer  :: istop   => NULL()
logical, dimension(:), pointer  :: isbot   => NULL()
#endif
end type                                         
                                                 
#ifdef USEMPI

interface space_distribute
module procedure space_distribute_matrix_real8
module procedure space_distribute_matrix_integer
module procedure space_distribute_block_real8
module procedure space_distribute_block4_real8
module procedure space_distribute_vector
module procedure space_distribute_block_vector
end interface space_distribute

interface space_distribute_p
module procedure space_distribute_matrix_real8_p
module procedure space_distribute_matrix_integer_p
module procedure space_distribute_block_real8_p
module procedure space_distribute_block4_real8_p
module procedure space_distribute_vector_p
end interface space_distribute_p

interface space_shift_borders
module procedure space_shift_borders_matrix_real8
module procedure space_shift_borders_block_real8
end interface space_shift_borders

interface space_collect
module procedure space_collect_matrix_real8
end interface space_collect

#endif

interface printsum
module procedure printsum0
module procedure printsum1
module procedure printsum2
module procedure printsum3
module procedure printsum4
module procedure printsumi0
module procedure printsumi1
module procedure printsumi2
module procedure printsumi3
end interface printsum


contains                                         
                                                 
subroutine grid_bathy(s,par)                         
          
use params                
use xmpi_module
use general_mpi_module
use readkey_module
                                                         
implicit none                                    
                                                 
type(spacepars)                     :: s         
type(parameters)                    :: par         

character*80                        :: fnameh,fnamex,fnamey
integer                             :: i
integer                             :: j
integer                             :: itheta
real*8                              :: degrad,thetamin,thetamax

!                     Input file  Keyword Default  Minimum  Maximum
s%nx    = readkey_int('params.txt','nx',     50,      2,     10000)
s%ny    = readkey_int('params.txt','ny',      2,      2,     10000)
s%dx    = readkey_dbl('params.txt','dx',    0.d0,   -1d9,      1d9)
s%dy    = readkey_dbl('params.txt','dy',    0.d0,   -1d9,      1d9)
s%xori  = readkey_dbl('params.txt','xori',  0.d0,   -1d9,      1d9)
s%yori  = readkey_dbl('params.txt','yori',  0.d0,   -1d9,      1d9)
s%alfa  = readkey_dbl('params.txt','alfa',  0.d0,   -360.d0,   360.d0)
s%alfa  = s%alfa*atan(1.0d0)/45.d0
s%posdwn= readkey_dbl('params.txt','posdwn',1.d0,   -1.d0,     1.d0)
s%vardx = readkey_int('params.txt','vardx',   0,      0,         1)     !Jaap

if (xreader) then
  allocate(s%x(1:s%nx+1,1:s%ny+1))
  allocate(s%y(1:s%nx+1,1:s%ny+1))
  allocate(s%xz(1:s%nx+1))
  allocate(s%yz(1:s%ny+1))
  allocate(s%xu(1:s%nx+1))
  allocate(s%yv(1:s%ny+1))
  allocate(s%xw(1:s%nx+1,1:s%ny+1))
  allocate(s%yw(1:s%nx+1,1:s%ny+1))
  allocate(s%zb(1:s%nx+1,1:s%ny+1))
  allocate(s%zb0(1:s%nx+1,1:s%ny+1))
endif
!
! Create grid
!
if(s%vardx==0)then

  call readkey('params.txt','depfile',fnameh)

  if (xreader) then
    open(31,file=fnameh)
    do j=1,s%ny+1
        read(31,*)(s%zb(i,j),i=1,s%nx+1)
    end do
    close(31)
  endif
  if (xreader) then
    s%zb=-s%zb*s%posdwn
    do j=1,s%ny+1
       do i=1,s%nx+1
         s%x(i,j)=(i-1)*s%dx
         s%y(i,j)=(j-1)*s%dy
       end do
    end do
  endif

elseif(s%vardx==1)then

  call readkey('params.txt','depfile',fnameh)
  call readkey('params.txt','xfile',fnamex)
  call readkey('params.txt','yfile',fnamey)

  if (xreader) then
    open(31,file=fnameh)
    open(32,file=fnamex)
    open(33,file=fnamey)
  endif
!  do j = 1,s%ny+1
!    if (xreader) then
!          read(31,*)(s%zb(i,j),i=1,s%nx+1)
!          read(32,*)(s%x(i,j),i=1,s%nx+1)
!          read(33,*)(s%y(i,j),i=1,s%nx+1)
!    endif
!  enddo
  if (xreader) then
    read(31,*)((s%zb(i,j),i=1,s%nx+1),j=1,s%ny+1)
    read(32,*)((s%x(i,j),i=1,s%nx+1),j=1,s%ny+1)
    read(33,*)((s%y(i,j),i=1,s%nx+1),j=1,s%ny+1)
    close(31)
    close(32)
    close(33)
    if (abs(s%x(1,2)-s%x(1,1))>par%eps) then
       ! Apparently input grid is at an angle, therefore defined in world coordinates
       ! Find out grid orientation
       s%alfa=atan2(s%y(2,1)-s%y(1,1),s%x(2,1)-s%x(1,1))
       s%xori=s%x(1,1)
       s%yori=s%y(1,1)
       s%xw=s%x
       s%yw=s%y
       s%x= cos(s%alfa)*(s%xw-s%xori)+sin(s%alfa)*(s%yw-s%yori)
       s%y=-sin(s%alfa)*(s%xw-s%xori)+cos(s%alfa)*(s%yw-s%yori)
    endif
  s%zb=-s%zb*s%posdwn
  endif

endif

if(xreader) then
  s%xz = s%x(:,1)
  s%yz = s%y(1,:)
  s%xu(1:s%nx) = 0.5d0*(s%xz(1:s%nx)+s%xz(2:s%nx+1))
  s%xu(s%nx+1) = s%xz(s%nx+1)+0.5d0*(s%xz(s%nx+1)-s%xz(s%nx))
  s%yv(1:s%ny) = 0.5d0*(s%yz(1:s%ny)+s%yz(2:s%ny+1))
  s%yv(s%ny+1) = s%yz(s%ny+1)+0.5d0*(s%yz(s%ny+1)-s%yz(s%ny))

  s%xw=s%xori+s%x*cos(s%alfa)-s%y*sin(s%alfa)
  s%yw=s%yori+s%x*sin(s%alfa)+s%y*cos(s%alfa)

  s%zb0  = s%zb
  !
  ! Specify theta-grid
  !
!
! from Nautical wave directions in degrees to Cartesian wave directions in radian !!!
!
  s%theta0=(1.5d0*par%px-s%alfa)-par%dir0*atan(1.d0)/45 ! Updated in waveparams.f90 for instat 4,5,6,7
  degrad=par%px/180.d0
  if (par%thetanaut==1) then  
     thetamin=(270-par%thetamax)*degrad-s%alfa
     thetamax=(270-par%thetamin)*degrad-s%alfa
     if (thetamax>par%px) then
        thetamax=thetamax-2*par%px
        thetamin=thetamin-2*par%px
     endif
     if (thetamin<-par%px) then
        thetamax=thetamax+2*par%px
        thetamin=thetamin+2*par%px
     endif
  else
  thetamin=par%thetamin*degrad
  thetamax=par%thetamax*degrad
  endif
  s%dtheta=par%dtheta*degrad
  s%ntheta=(thetamax-thetamin)/s%dtheta

  allocate(s%theta(1:s%ntheta))
  allocate(s%cxsth(1:s%ntheta))
  allocate(s%sxnth(1:s%ntheta))

  do itheta=1,s%ntheta
      s%theta(itheta)=thetamin+s%dtheta/2+s%dtheta*(itheta-1)
  end do

  s%cxsth=cos(s%theta)
  s%sxnth=sin(s%theta)
endif
end subroutine grid_bathy                         

#ifdef USEMPI

subroutine space_distribute_matrix_real8(sg,sl,a,b)
use xmpi_module
use general_mpi_module
implicit none
type (spacepars), intent(in)        :: sg
type (spacepars), intent(inout)     :: sl
real*8, dimension(:,:), intent(in)  :: a
real*8, dimension(:,:), intent(out) :: b

!wwvv to quiet the compiler: sl is indeed not used
sl%nx = sl%nx

call matrix_distr(a,b,sg%is,sg%lm,sg%js,sg%ln,xmpi_master,xmpi_comm)

end subroutine space_distribute_matrix_real8

subroutine space_distribute_matrix_real8_p(sg,sl,a,b)
use xmpi_module
use general_mpi_module
implicit none
type (spacepars), intent(in)        :: sg
type (spacepars), intent(inout)     :: sl
real*8, dimension(:,:), intent(in)  :: a
real*8, pointer, dimension(:,:)     :: b

if (associated(b)) then
  deallocate(b)
endif

allocate(b(sl%nx+1,sl%ny+1))

call space_distribute(sg,sl,a,b)

end subroutine space_distribute_matrix_real8_p

subroutine space_distribute_matrix_integer(sg,sl,a,b)
use xmpi_module
use general_mpi_module
implicit none
type (spacepars), intent(in)         :: sg
type (spacepars), intent(inout)      :: sl
integer, dimension(:,:), intent(in)  :: a
integer, dimension(:,:), intent(out) :: b

!wwvv to quiet the compiler: sl is indeed not used
sl%nx = sl%nx

call matrix_distr(a,b,sg%is,sg%lm,sg%js,sg%ln,xmpi_master,xmpi_comm)

end subroutine space_distribute_matrix_integer

subroutine space_distribute_matrix_integer_p(sg,sl,a,b)
use xmpi_module
use general_mpi_module
implicit none
type (spacepars), intent(in)        :: sg
type (spacepars), intent(inout)     :: sl
integer, dimension(:,:), intent(in)  :: a
integer, pointer, dimension(:,:)     :: b

if (associated(b)) then
  deallocate(b)
endif

allocate(b(sl%nx+1,sl%ny+1))

call space_distribute(sg,sl,a,b)

end subroutine space_distribute_matrix_integer_p

subroutine space_distribute_block_real8(sg,a,b,d3)
use xmpi_module
use general_mpi_module
implicit none
type (spacepars), intent(in)                   :: sg
real*8, dimension(:,:,:), intent(in)           :: a
real*8, dimension(:,:,:), intent(out)          :: b
integer, intent(in)                            :: d3

integer                                        :: i

do i=1,d3
  call matrix_distr(a(:,:,i),b(:,:,i),sg%is,sg%lm,sg%js,sg%ln,xmpi_master,xmpi_comm)
enddo

end subroutine space_distribute_block_real8

subroutine space_distribute_block4_real8(sg,a,b,d3,d4)
use xmpi_module
use general_mpi_module
implicit none
type (spacepars), intent(in)                   :: sg
real*8, dimension(:,:,:,:), intent(in)         :: a
real*8, dimension(:,:,:,:), intent(out)        :: b
integer, intent(in)                            :: d3,d4

integer                                        :: i,j

do i=1,d3
  do j=1,d4
    call matrix_distr(a(:,:,i,j),b(:,:,i,j),sg%is,sg%lm,sg%js,sg%ln,xmpi_master,xmpi_comm)
  enddo
enddo

end subroutine space_distribute_block4_real8

subroutine space_distribute_block_real8_p(sg,sl,a,b)
use xmpi_module
use general_mpi_module
implicit none
type (spacepars), intent(in)                   :: sg
type (spacepars), intent(inout)                :: sl
real*8, dimension(:,:,:), intent(in)           :: a
real*8, pointer, dimension(:,:,:)              :: b
integer                                        :: d3

if (associated(b)) then
  deallocate(b)
endif

if (xmpi_rank .eq. xmpi_master) then
  d3=size(a,3)
endif
call xmpi_bcast(d3)

allocate(b(sl%nx+1,sl%ny+1,d3))

call space_distribute(sg,a,b,d3)

end subroutine space_distribute_block_real8_p

subroutine space_distribute_block4_real8_p(sg,sl,a,b)
use xmpi_module
use general_mpi_module
implicit none
type (spacepars), intent(in)                   :: sg
type (spacepars), intent(inout)                :: sl
real*8, dimension(:,:,:,:), intent(in)         :: a
real*8, pointer, dimension(:,:,:,:)            :: b
integer                                        :: d3,d4

if (associated(b)) then
  deallocate(b)
endif

if (xmpi_rank .eq. xmpi_master) then
  d3=size(a,3)
  d4=size(a,4)
endif
call xmpi_bcast(d3)
call xmpi_bcast(d4)

allocate(b(sl%nx+1,sl%ny+1,d3,d4))

call space_distribute(sg,a,b,d3,d4)

end subroutine space_distribute_block4_real8_p
 
subroutine space_distribute_vector(xy,sg,sl,a,b)
use xmpi_module
use general_mpi_module
implicit none
character, intent(in)             :: xy
type (spacepars), intent(in)      :: sg
type (spacepars), intent(inout)   :: sl
real*8, dimension(:), intent(in)  :: a
real*8, dimension(:), intent(out) :: b

integer, dimension(:), pointer    :: ijs,lmn
integer                           :: nxy

if(xy .eq.'x') then
  ijs => sg%is
  lmn => sg%lm
  nxy =  sl%nx
else
  ijs => sg%js
  lmn => sg%ln
  nxy =  sl%ny
endif

call vector_distr_send(a,b,ijs,lmn,xmpi_master,xmpi_comm)

end subroutine space_distribute_vector

subroutine space_distribute_vector_p(xy,sg,sl,a,b)
use xmpi_module
use general_mpi_module
implicit none
character, intent(in)             :: xy
type (spacepars), intent(in)      :: sg
type (spacepars), intent(inout)   :: sl
real*8, dimension(:), intent(in)  :: a
real*8, dimension(:), pointer     :: b

integer                           :: nxy

if(associated(b)) then
  deallocate(b)
endif

if(xy .eq.'x') then
  nxy =  sl%nx
else
  nxy =  sl%ny
endif

allocate(b(nxy+1))

call space_distribute(xy,sg,sl,a,b)

end subroutine space_distribute_vector_p

subroutine space_distribute_block_vector(xy,sg,sl,a,b)
use xmpi_module
use general_mpi_module
implicit none
character, intent(in)               :: xy
type (spacepars), intent(in)        :: sg
type (spacepars), intent(inout)     :: sl
real*8, dimension(:,:), intent(in)  :: a
real*8, dimension(:,:), intent(out) :: b

integer                             :: i

do i=1,sl%ntheta
  call space_distribute(xy,sg,sl,a(:,i),b(:,i))
enddo

end subroutine space_distribute_block_vector

subroutine space_distribute_space(sg,sl,par)
use xmpi_module
use general_mpi_module
use params
implicit none
type(spacepars), intent(inout)  :: sg
type(spacepars), intent(out)    :: sl 
type(parameters)                :: par

integer                         :: i
real*8, pointer, dimension(:)   :: umean1, umean2

!
! This subroutine takes care that all contents of the global
! space is distributed to the local space.
!
! Also, the isleft, isright, istop and isbot logicals from
! the xmpi module are put in sg and sl.
!

!
! first the scalars, except nx and ny
!

if (xmpi_rank .eq. xmpi_master) then
  sl%dx         = sg%dx
  sl%dy         = sg%dy
  sl%xori       = sg%xori      
  sl%yori       = sg%yori      
  sl%alfa       = sg%alfa      
  sl%posdwn     = sg%posdwn    
! sl%nx         = sg%nx        
! sl%ny         = sg%ny        
  sl%ntheta     = sg%ntheta    
  sl%dtheta     = sg%dtheta    
  sl%theta0     = sg%theta0    
  sl%vardx      = sg%vardx
endif

  call xmpi_bcast(sl%dx)
  call xmpi_bcast(sl%dy)
  call xmpi_bcast(sl%xori)
  call xmpi_bcast(sl%yori)
  call xmpi_bcast(sl%alfa)
  call xmpi_bcast(sl%posdwn)
! call xmpi_bcast(sl%nx)
! call xmpi_bcast(sl%ny)
  call xmpi_bcast(sl%ntheta)
  call xmpi_bcast(sl%dtheta)
  call xmpi_bcast(sl%theta0)
  call xmpi_bcast(sl%vardx)

! 
! Distribute is,js,ln,lm,isleft,isright,istop,isbot
!

  if (associated(sg%is)) then
    deallocate(sg%is)
    deallocate(sg%js)
    deallocate(sg%lm)
    deallocate(sg%ln)
    deallocate(sg%isleft)
    deallocate(sg%isright)
    deallocate(sg%istop)
    deallocate(sg%isbot)
  endif


  allocate(sg%is(xmpi_size))
  allocate(sg%js(xmpi_size))
  allocate(sg%lm(xmpi_size))
  allocate(sg%ln(xmpi_size))
  allocate(sg%isleft(xmpi_size))
  allocate(sg%isright(xmpi_size))
  allocate(sg%istop(xmpi_size))
  allocate(sg%isbot(xmpi_size))

  if (associated(sl%is)) then
    deallocate(sl%is)
    deallocate(sl%js)
    deallocate(sl%lm)
    deallocate(sl%ln)
    deallocate(sl%isleft)
    deallocate(sl%isright)
    deallocate(sl%istop)
    deallocate(sl%isbot)
  endif
  
  allocate(sl%is(xmpi_size))
  allocate(sl%js(xmpi_size))
  allocate(sl%lm(xmpi_size))
  allocate(sl%ln(xmpi_size))
  allocate(sl%isleft(xmpi_size))
  allocate(sl%isright(xmpi_size))
  allocate(sl%istop(xmpi_size))
  allocate(sl%isbot(xmpi_size))

  call det_submatrices(sg%nx+1, sg%ny+1, xmpi_m, xmpi_n, &
                             sg%is, sg%lm, sg%js, sg%ln, &
                             sg%isleft, sg%isright, sg%istop, sg%isbot)
  if (xreader) then
    print *,'Distribution of matrix on processors'
    print *,' proc   is   lm   js   ln'
    do i=1,xmpi_size
       print "(1x,5i5)",i-1,sg%is(i),sg%lm(i),sg%js(i),sg%ln(i)
    enddo
    print *,' proc   left right top bot'
    do i=1,xmpi_size
       print "(1x,i5,4l5)",i-1,sg%isleft(i),sg%isright(i),sg%istop(i),sg%isbot(i)
    enddo
  endif


  if (xmpi_rank .eq. xmpi_master) then
    sl%is       = sg%is
    sl%js       = sg%js
    sl%lm       = sg%lm
    sl%ln       = sg%ln
    sl%isleft   = sg%isleft
    sl%isright  = sg%isright
    sl%istop    = sg%istop
    sl%isbot    = sg%isbot
  endif

  call xmpi_bcast(sl%is)
  call xmpi_bcast(sl%js)
  call xmpi_bcast(sl%lm)
  call xmpi_bcast(sl%ln)
  call xmpi_bcast(sl%isleft)
  call xmpi_bcast(sl%isright)
  call xmpi_bcast(sl%istop)
  call xmpi_bcast(sl%isbot)

!
! compute the values for local nx and ny
!

  sl%nx = sl%lm(xmpi_rank+1) - 1 
  sl%ny = sl%ln(xmpi_rank+1) - 1

! 
! distribute all arrays and matrices
!

! dimension s%ntheta:
  if (associated(sl%theta)) then
    deallocate(sl%theta)
    deallocate(sl%cxsth)
    deallocate(sl%sxnth)
  endif

  allocate(sl%theta(sl%ntheta))
  allocate(sl%cxsth(sl%ntheta))
  allocate(sl%sxnth(sl%ntheta))

  if (xmpi_rank .eq. xmpi_master) then
    sl%theta = sg%theta
    sl%cxsth = sg%cxsth
    sl%sxnth = sg%sxnth
  endif

  call xmpi_bcast(sl%theta) 
  call xmpi_bcast(sl%cxsth) 
  call xmpi_bcast(sl%sxnth) 
  
! dimension par%ngd

  if (associated(sl%D50)) then
    deallocate(sl%D50)
    deallocate(sl%D90)
    deallocate(sl%sedcal)
  endif
  allocate(sl%D50(par%ngd))
  allocate(sl%D90(par%ngd))
  allocate(sl%sedcal(par%ngd))
  if (xmpi_rank .eq. xmpi_master) then
    sl%D50    = sg%D50
    sl%D90    = sg%D90
    sl%sedcal = sg%sedcal
  endif
  
! dimension s%nx+1:
  if (associated(sl%xz)) then
    deallocate(sl%xz)
    deallocate(sl%xu)
  endif
  allocate(sl%xz(sl%nx+1))
  allocate(sl%xu(sl%nx+1))
  call space_distribute("x", sg, sl, sg%xz, sl%xz)
  call space_distribute("x", sg, sl, sg%xu, sl%xu)

! dimension s%ny+1
  
  if (associated(sl%yz)) then
    deallocate(sl%yz)
    deallocate(sl%yv)
    deallocate(sl%bi)
  endif
  allocate(sl%yz(sl%ny+1))
  allocate(sl%yv(sl%ny+1))
  allocate(sl%bi(sl%ny+1))

  call space_distribute("y", sg, sl, sg%yz, sl%yz)
  call space_distribute("y", sg, sl, sg%yv, sl%yv)
  call space_distribute("y", sg, sl, sg%bi, sl%bi)
! dimension par%tidelen
  if (associated(sl%tideinpt)) then
    deallocate(sl%tideinpt)
  endif
  allocate(sl%tideinpt(par%tidelen))
  if (xmpi_rank .eq. xmpi_master) then
    sl%tideinpt = sg%tideinpt
  endif
  call xmpi_bcast(sl%tideinpt) 
  
! dimension s%nx+1, s%ny+1
#if 0
  if(associated(sl%x)) then
    deallocate(sl%x)
    deallocate(sl%y)
    deallocate(sl%xw)
    deallocate(sl%yw)
    deallocate(sl%zb)
    deallocate(sl%zb0)
    deallocate(sl%thetamean)
    deallocate(sl%Fx)
    deallocate(sl%Fy)
    deallocate(sl%Sxy)
    deallocate(sl%Syy)
    deallocate(sl%Sxx)
    deallocate(sl%n)
    deallocate(sl%H)
    deallocate(sl%k)
    deallocate(sl%c)
    deallocate(sl%cg)
    deallocate(sl%sigm)
    deallocate(sl%hh)
    deallocate(sl%zs)
    deallocate(sl%zs0)
    deallocate(sl%dzsdt)
    deallocate(sl%dzbdt)
    deallocate(sl%uu)
    deallocate(sl%vv)
    deallocate(sl%qx)
    deallocate(sl%qy)
    deallocate(sl%sedero)
    deallocate(sl%dcdx)
    deallocate(sl%dcdy)
    deallocate(sl%ui)
    deallocate(sl%E)
    deallocate(sl%R)
    deallocate(sl%urms)
    deallocate(sl%D)
    deallocate(sl%ust)
    deallocate(sl%tm)
    deallocate(sl%ueu)
    deallocate(sl%vev)
    deallocate(sl%vmagu)
    deallocate(vmageu)
    deallocate(sl%vmagv)
    deallocate(vmagev)
    deallocate(sl%u)
    deallocate(sl%v)
    deallocate(sl%ue)
    deallocate(sl%ve)
    deallocate(sl%hold)
    deallocate(sl%wetu)
    deallocate(sl%wetv)
    deallocate(sl%wetz)
    deallocate(sl%hu)
    deallocate(sl%hv)
    deallocate(sl%hum)
    deallocate(sl%hvm)
    deallocate(sl%vmag)
    deallocate(sl%ccg)
    deallocate(sl%uwf)
    deallocate(sl%vwf)
    deallocate(sl%ustr)
    deallocate(sl%usd)
    deallocate(sl%DR)
    deallocate(sl%vu)
    deallocate(sl%uv)
    deallocate(sl%graindistr)
    deallocate(sl%Tsg)
    deallocate(sl%Sug)
    deallocate(sl%Svg)
    deallocate(sl%ua)
    deallocate(sl%BR)
    deallocate(sl%kb)
    deallocate(sl%Tbore)
    deallocate(sl%uon)
    deallocate(sl%uoff)
    deallocate(sl%dzav)
    deallocate(sl%maxzs)
    deallocate(sl%minzs)
  endif
  allocate(sl%x(sl%nx+1,sl%ny+1))
  allocate(sl%y(sl%nx+1,sl%ny+1))
  allocate(sl%xw(sl%nx+1,sl%ny+1))
  allocate(sl%yw(sl%nx+1,sl%ny+1))
  allocate(sl%zb(sl%nx+1,sl%ny+1))
  allocate(sl%zb0(sl%nx+1,sl%ny+1))
  allocate(sl%thetamean(sl%nx+1,sl%ny+1))
  allocate(sl%Fx(sl%nx+1,sl%ny+1))
  allocate(sl%Fy(sl%nx+1,sl%ny+1))
  allocate(sl%Sxy(sl%nx+1,sl%ny+1))
  allocate(sl%Syy(sl%nx+1,sl%ny+1))
  allocate(sl%Sxx(sl%nx+1,sl%ny+1))
  allocate(sl%n(sl%nx+1,sl%ny+1))
  allocate(sl%H(sl%nx+1,sl%ny+1))
  allocate(sl%k(sl%nx+1,sl%ny+1))
  allocate(sl%c(sl%nx+1,sl%ny+1))
  allocate(sl%cg(sl%nx+1,sl%ny+1))
  allocate(sl%sigm(sl%nx+1,sl%ny+1))
  allocate(sl%hh(sl%nx+1,sl%ny+1))
  allocate(sl%zs(sl%nx+1,sl%ny+1))
  allocate(sl%zs0(sl%nx+1,sl%ny+1))
  allocate(sl%dzsdt(sl%nx+1,sl%ny+1))
  allocate(sl%dzbdt(sl%nx+1,sl%ny+1))
  allocate(sl%uu(sl%nx+1,sl%ny+1))
  allocate(sl%vv(sl%nx+1,sl%ny+1))
  allocate(sl%qx(sl%nx+1,sl%ny+1))
  allocate(sl%qy(sl%nx+1,sl%ny+1))
  allocate(sl%sedero(sl%nx+1,sl%ny+1))
  allocate(sl%dcdx(sl%nx+1,sl%ny+1))
  allocate(sl%dcdy(sl%nx+1,sl%ny+1))
  allocate(sl%ui(sl%nx+1,sl%ny+1))
  allocate(sl%E(sl%nx+1,sl%ny+1))
  allocate(sl%R(sl%nx+1,sl%ny+1))
  allocate(sl%urms(sl%nx+1,sl%ny+1))
  allocate(sl%D(sl%nx+1,sl%ny+1))
  allocate(sl%ust(sl%nx+1,sl%ny+1))
  allocate(sl%tm(sl%nx+1,sl%ny+1))
  allocate(sl%ueu(sl%nx+1,sl%ny+1))
  allocate(sl%vev(sl%nx+1,sl%ny+1))
  allocate(sl%vmagu(sl%nx+1,sl%ny+1))
  allocate(sl%vmageu(sl%nx+1,sl%ny+1))
  allocate(sl%vmagv(sl%nx+1,sl%ny+1))
  allocate(sl%vmagev(sl%nx+1,sl%ny+1))
  allocate(sl%u(sl%nx+1,sl%ny+1))
  allocate(sl%v(sl%nx+1,sl%ny+1))
  allocate(sl%ue(sl%nx+1,sl%ny+1))
  allocate(sl%ve(sl%nx+1,sl%ny+1))
  allocate(sl%hold(sl%nx+1,sl%ny+1))
  allocate(sl%wetu(sl%nx+1,sl%ny+1))
  allocate(sl%wetv(sl%nx+1,sl%ny+1))
  allocate(sl%wetz(sl%nx+1,sl%ny+1))
  allocate(sl%hu(sl%nx+1,sl%ny+1))
  allocate(sl%hv(sl%nx+1,sl%ny+1))
  allocate(sl%hum(sl%nx+1,sl%ny+1))
  allocate(sl%hvm(sl%nx+1,sl%ny+1))
  allocate(sl%vmag(sl%nx+1,sl%ny+1))
  allocate(sl%ccg(sl%nx+1,sl%ny+1))
  allocate(sl%uwf(sl%nx+1,sl%ny+1))
  allocate(sl%vwf(sl%nx+1,sl%ny+1))
  allocate(sl%ustr(sl%nx+1,sl%ny+1))
  allocate(sl%usd(sl%nx+1,sl%ny+1))
  allocate(sl%DR(sl%nx+1,sl%ny+1))
  allocate(sl%vu(sl%nx+1,sl%ny+1))
  allocate(sl%uv(sl%nx+1,sl%ny+1))
  allocate(sl%graindistr(sl%nx+1,sl%ny+1))
  allocate(sl%Tsg(sl%nx+1,sl%ny+1))
  allocate(sl%Sug(sl%nx+1,sl%ny+1))
  allocate(sl%Svg(sl%nx+1,sl%ny+1))
  allocate(sl%ceqg(sl%nx+1,sl%ny+1))
  allocate(sl%ua(sl%nx+1,sl%ny+1))
  allocate(sl%BR(sl%nx+1,sl%ny+1))
  allocate(sl%kb(sl%nx+1,sl%ny+1))
  allocate(sl%Tbore(sl%nx+1,sl%ny+1))
  allocate(sl%uon(sl%nx+1,sl%ny+1))
  allocate(sl%uoff(sl%nx+1,sl%ny+1))
  allocate(sl%dzav(sl%nx+1,sl%ny+1))
  allocate(sl%maxzs(sl%nx+1,sl%ny+1))
  allocate(sl%minzs(sl%nx+1,sl%ny+1))

#endif

  call space_distribute_p(sg,sl,sg%x,sl%x)
  call space_distribute_p(sg,sl,sg%y,sl%y)
  call space_distribute_p(sg,sl,sg%xw,sl%xw)
  call space_distribute_p(sg,sl,sg%yw,sl%yw)
  call space_distribute_p(sg,sl,sg%zb,sl%zb)
  call space_distribute_p(sg,sl,sg%zb0,sl%zb0)
  call space_distribute_p(sg,sl,sg%thetamean,sl%thetamean)
  call space_distribute_p(sg,sl,sg%Fx,sl%Fx)
  call space_distribute_p(sg,sl,sg%Fy,sl%Fy)
  call space_distribute_p(sg,sl,sg%Sxy,sl%Sxy)
  call space_distribute_p(sg,sl,sg%Syy,sl%Syy)
  call space_distribute_p(sg,sl,sg%Sxx,sl%Sxx)
  call space_distribute_p(sg,sl,sg%n,sl%n)
  call space_distribute_p(sg,sl,sg%H,sl%H)
  call space_distribute_p(sg,sl,sg%k,sl%k)
  call space_distribute_p(sg,sl,sg%c,sl%c)
  call space_distribute_p(sg,sl,sg%cg,sl%cg)
  call space_distribute_p(sg,sl,sg%sigm,sl%sigm)
  call space_distribute_p(sg,sl,sg%hh,sl%hh)
  call space_distribute_p(sg,sl,sg%zs,sl%zs)
  call space_distribute_p(sg,sl,sg%zs0,sl%zs0)
  call space_distribute_p(sg,sl,sg%dzsdt,sl%dzsdt)
  call space_distribute_p(sg,sl,sg%dzbdt,sl%dzbdt)
  call space_distribute_p(sg,sl,sg%uu,sl%uu)
  call space_distribute_p(sg,sl,sg%vv,sl%vv)
  call space_distribute_p(sg,sl,sg%qx,sl%qx)
  call space_distribute_p(sg,sl,sg%qy,sl%qy)
  call space_distribute_p(sg,sl,sg%sedero,sl%sedero)
  call space_distribute_p(sg,sl,sg%dcdx,sl%dcdx)
  call space_distribute_p(sg,sl,sg%dcdy,sl%dcdy)
  call space_distribute_p(sg,sl,sg%ui,sl%ui)
  call space_distribute_p(sg,sl,sg%E,sl%E)
  call space_distribute_p(sg,sl,sg%R,sl%R)
  call space_distribute_p(sg,sl,sg%urms,sl%urms)
  call space_distribute_p(sg,sl,sg%D,sl%D)
  call space_distribute_p(sg,sl,sg%ust,sl%ust)
  call space_distribute_p(sg,sl,sg%tm,sl%tm)
  call space_distribute_p(sg,sl,sg%ueu,sl%ueu)
  call space_distribute_p(sg,sl,sg%vev,sl%vev)
  call space_distribute_p(sg,sl,sg%vmagu,sl%vmagu)
  call space_distribute_p(sg,sl,sg%vmageu,sl%vmageu)
  call space_distribute_p(sg,sl,sg%vmagv,sl%vmagv)
  call space_distribute_p(sg,sl,sg%vmagev,sl%vmagev)
  call space_distribute_p(sg,sl,sg%u,sl%u)
  call space_distribute_p(sg,sl,sg%v,sl%v)
  call space_distribute_p(sg,sl,sg%ue,sl%ue)
  call space_distribute_p(sg,sl,sg%ve,sl%ve)
  call space_distribute_p(sg,sl,sg%hold,sl%hold)
  call space_distribute_p(sg,sl,sg%wetu,sl%wetu)
  call space_distribute_p(sg,sl,sg%wetv,sl%wetv)
  call space_distribute_p(sg,sl,sg%wetz,sl%wetz)
  call space_distribute_p(sg,sl,sg%hu,sl%hu)
  call space_distribute_p(sg,sl,sg%hv,sl%hv)
  call space_distribute_p(sg,sl,sg%hum,sl%hum)
  call space_distribute_p(sg,sl,sg%hvm,sl%hvm)
  call space_distribute_p(sg,sl,sg%vmag,sl%vmag)
  call space_distribute_p(sg,sl,sg%ccg,sl%ccg)
  call space_distribute_p(sg,sl,sg%uwf,sl%uwf)
  call space_distribute_p(sg,sl,sg%vwf,sl%vwf)
  call space_distribute_p(sg,sl,sg%ustr,sl%ustr)
  call space_distribute_p(sg,sl,sg%usd,sl%usd)
  call space_distribute_p(sg,sl,sg%DR,sl%DR)
  call space_distribute_p(sg,sl,sg%vu,sl%vu)
  call space_distribute_p(sg,sl,sg%uv,sl%uv)
  call space_distribute_p(sg,sl,sg%Tsg,sl%Tsg)
  call space_distribute_p(sg,sl,sg%Sug,sl%Sug)
  call space_distribute_p(sg,sl,sg%Svg,sl%Svg)
  call space_distribute_p(sg,sl,sg%ceqg,sl%ceqg)
  call space_distribute_p(sg,sl,sg%ua,sl%ua)
  call space_distribute_p(sg,sl,sg%BR,sl%BR)
  call space_distribute_p(sg,sl,sg%kb,sl%kb)
  call space_distribute_p(sg,sl,sg%Tbore,sl%Tbore)
  call space_distribute_p(sg,sl,sg%uon,sl%uon)
  call space_distribute_p(sg,sl,sg%uoff,sl%uoff)
  call space_distribute_p(sg,sl,sg%dzav,sl%dzav)
  call space_distribute_p(sg,sl,sg%maxzs,sl%maxzs)
  call space_distribute_p(sg,sl,sg%minzs,sl%minzs)


#if 0
 x         
 y         
 xw        
 yw        
 zb        
 zb0       
 thetamean 
 Fx        
 Fy        
 Sxy       
 Syy       
 Sxx       
 n         
 H         
 k         
 c         
 cg        
 sigm      
 hh        
 zs        
 zs0       
 dzsdt     
 dzbdt     
 uu        
 vv        
 qx        
 qy        
 sedero    
 dcdx      
 dcdy      
 ui        
 E         
 R         
 urms      
 D         
 ust       
 tm        
 ueu       
 vev       
 vmagu     
 vmageu     
 vmagv     
 vmagev     
 u         
 v         
 ue        
 ve        
 hold      
 wetu      
 wetv      
 wetz      
 hu        
 hv        
 hum       
 hvm       
 vmag      
 ccg       
 uwf       
 vwf       
 ustr      
 usd       
 DR        
 vu        
 uv        
  graindistr
  D50
  D90
  sedcal
  Tsg
  Sug
  Svg
  ceqg
  ua
  BR
  kb
  Tbore
  uon
  uoff
  dzav
  maxzs
  minzs
#endif
! dimension s%nx+1, s%ny+1, s%ntheta
  call space_distribute_p(sg,sl,sg%cgx,sl%cgx)
  call space_distribute_p(sg,sl,sg%cgy,sl%cgy)
  call space_distribute_p(sg,sl,sg%cx,sl%cx)
  call space_distribute_p(sg,sl,sg%cy,sl%cy)
  call space_distribute_p(sg,sl,sg%ctheta,sl%ctheta)
  call space_distribute_p(sg,sl,sg%ee,sl%ee)
  call space_distribute_p(sg,sl,sg%thet,sl%thet)
  call space_distribute_p(sg,sl,sg%costhet,sl%costhet)
  call space_distribute_p(sg,sl,sg%sinthet,sl%sinthet)
  call space_distribute_p(sg,sl,sg%sigt,sl%sigt)
  call space_distribute_p(sg,sl,sg%rr,sl%rr)

! dimension s%nx+1,s%ny+1,par%nd,par%ngd
  call space_distribute_p(sg,sl,sg%graindistr,sl%graindistr)

! dimension 2,s%ny+1
  allocate(umean1(sg%ny+1))
  allocate(umean2(sl%ny+1))
  if (associated(sl%umean)) then
    deallocate(sl%umean)
  endif
  allocate(sl%umean(2,sl%ny+1))
  do i=1,2
    if (xmpi_rank .eq. xmpi_master) then
      umean1 = sg%umean(i,:)
    endif
    call space_distribute("y", sg,sl,umean1,umean2)
    sl%umean(i,:) = umean2
  enddo
  deallocate(umean1,umean2)

! dimension par%tidelen, par%tideloc
  if (associated(sl%tideinpz)) then
    deallocate(sl%tideinpz)
  endif
  allocate(sl%tideinpz(par%tidelen,par%tideloc))
  if (xmpi_rank .eq. xmpi_master) then
    sl%tideinpz = sg%tideinpz
  endif
  call xmpi_bcast(sl%tideinpz) 
end subroutine space_distribute_space

subroutine space_shift_borders_matrix_real8(s,a)
  use general_mpi_module
  use xmpi_module
  implicit none
  type(spacepars), intent(in)          :: s
  real*8, intent(inout),dimension(:,:) :: a
!
! s is not used here, but for symmetry, compared with the block
! variety, we have it as first parameter
! to quiet the compiler, do something with it
!
  integer i
  do i=1, s%ntheta - s%ntheta + 1
    call shift_borders_matrix_real8(a,xmpi_left,xmpi_right, &
             xmpi_top,xmpi_bot,xmpi_comm)
  enddo
end subroutine space_shift_borders_matrix_real8

subroutine space_shift_borders_block_real8(s,a)
  use general_mpi_module
  use xmpi_module
  implicit none
  type (spacepars), intent(in)           :: s
  real*8, intent(inout),dimension(:,:,:) :: a

  integer i

  do i=1,s%ntheta
    call shift_borders_matrix_real8(a(:,:,i),xmpi_left,xmpi_right,&
              xmpi_top,xmpi_bot,xmpi_comm)
  enddo
end subroutine space_shift_borders_block_real8

subroutine space_collect_matrix_real8(s,a,b)
use general_mpi_module
use xmpi_module
implicit none
type(spacepars), intent(in)          :: s
real*8, dimension(:,:), intent(out)  :: a
real*8, dimension(:,:), intent(in)   :: b

call matrix_coll(a,b,s%is,s%lm,s%js,s%ln, &
                 s%isleft,s%isright,s%istop,s%isbot, &
                 xmpi_master,xmpi_comm)

end subroutine space_collect_matrix_real8
#endif

subroutine printsum0(f,str,id,val)
implicit none
integer, intent(in) :: f
character(*),intent(in) :: str
integer, intent(in) :: id
real*8, intent(in) :: val
write(f,*) 'printsum ',id,' ',str,':',val
end subroutine printsum0

subroutine printsum1(f,str,id,val)
implicit none
integer, intent(in) :: f
character(*), intent(in) :: str
integer, intent(in) :: id
real*8, pointer, dimension(:) :: val
if (associated(val) ) then
  write(f,*) 'printsum ',id,' ',str,':',sum(val)
else
  write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
endif
end subroutine printsum1

subroutine printsum2(f,str,id,val)
implicit none
integer, intent(in) :: f
character(*), intent(in) :: str
integer, intent(in) :: id
real*8, pointer, dimension(:,:) :: val
if (associated(val) ) then
  write(f,*) 'printsum ',id,' ',str,':',sum(val)
else
  write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
endif

end subroutine printsum2

subroutine printsum3(f,str,id,val)
implicit none
integer, intent(in) :: f
character(*), intent(in) :: str
integer, intent(in) :: id
real*8, pointer, dimension(:,:,:) :: val
if (associated(val) ) then
  write(f,*) 'printsum ',id,' ',str,':',sum(val)
else
  write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
endif
end subroutine printsum3

subroutine printsum4(f,str,id,val)
implicit none
integer, intent(in) :: f
character(*), intent(in) :: str
integer, intent(in) :: id
real*8, pointer, dimension(:,:,:,:) :: val
if (associated(val) ) then
  write(f,*) 'printsum ',id,' ',str,':',sum(val)
else
  write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
endif
end subroutine printsum4

subroutine printsumi0(f,str,id,val)
implicit none
integer, intent(in) :: f
character(*),intent(in) :: str
integer, intent(in) :: id
integer, intent(in) :: val
write(f,*) 'printsum ',id,' ',str,':',val
end subroutine printsumi0

subroutine printsumi1(f,str,id,val)
implicit none
integer, intent(in) :: f
character(*), intent(in) :: str
integer, intent(in) :: id
integer, pointer, dimension(:) :: val
if (associated(val) ) then
  write(f,*) 'printsum ',id,' ',str,':',sum(val)
else
  write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
endif
end subroutine printsumi1

subroutine printsumi2(f,str,id,val)
implicit none
integer, intent(in) :: f
character(*), intent(in) :: str
integer, intent(in) :: id
integer, pointer, dimension(:,:) :: val
if (associated(val) ) then
  write(f,*) 'printsum ',id,' ',str,':',sum(val)
else
  write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
endif

end subroutine printsumi2

subroutine printsumi3(f,str,id,val)
implicit none
integer, intent(in) :: f
character(*), intent(in) :: str
integer, intent(in) :: id
integer, pointer, dimension(:,:,:) :: val
if (associated(val) ) then
  write(f,*) 'printsum ',id,' ',str,':',sum(val)
else
  write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
endif
end subroutine printsumi3

subroutine printssums(s,str)
#ifdef USEMPI
use xmpi_module
#endif
type (spacepars), intent(in) :: s
character(*), intent(in) :: str
integer :: id, f

#ifdef USEMPI
id = xmpi_rank
if (id .gt. 0 ) then
  return
endif
#else
id=0
#endif

f = 170+id

write(f, *) 'printsum: ',id,'Start of printssums ',str

call printsum(f,'s%x',id,s%x)
call printsum(f,'s%y',id,s%y)
call printsum(f,'s%xw',id,s%y)
call printsum(f,'s%yw',id,s%y)
call printsum(f,'s%dx',id,s%dx)
call printsum(f,'s%dy',id,s%dy)
call printsum(f,'s%xz',id,s%xz)
call printsum(f,'s%yz',id,s%yz)
call printsum(f,'s%xu',id,s% xu)
call printsum(f,'s%yv',id,s%yv)
call printsum(f,'s%nx',id,s%nx)
call printsum(f,'s%ny',id,s%ny)
call printsum(f,'s%zb',id,s%zb)
call printsum(f,'s%zb0',id,s%zb0)
call printsum(f,'s%theta',id,s%theta)
call printsum(f,'s%ntheta',id,s%ntheta)
call printsum(f,'s%dtheta',id,s%dtheta)
call printsum(f,'s%theta0',id,s%theta0)
call printsum(f,'s%cxsth',id,s%cxsth)
call printsum(f,'s%sxnth',id,s%sxnth)
call printsum(f,'s%thetamean',id,s%thetamean)
call printsum(f,'s%Fx',id,s%Fx)
call printsum(f,'s%Fy',id,s%Fy)
call printsum(f,'s%Sxy',id,s%Sxy)
call printsum(f,'s%Syy',id,s%Syy)
call printsum(f,'s%Sxx',id,s%Sxx)
call printsum(f,'s%n',id,s%n)
call printsum(f,'s%H',id,s%H)
call printsum(f,'s%cgx',id,s%cgx)
call printsum(f,'s%cgy',id,s%cgy)
call printsum(f,'s%cx',id,s%cx)
call printsum(f,'s%cy',id,s%cy)
call printsum(f,'s%ctheta',id,s%ctheta)
call printsum(f,'s%ee',id,s%ee)
call printsum(f,'s%rr',id,s%rr)
call printsum(f,'s%thet',id,s%thet)
call printsum(f,'s%costhet',id,s%costhet)
call printsum(f,'s%sinthet',id,s%sinthet)
call printsum(f,'s%sigt',id,s%sigt)
call printsum(f,'s%k',id,s%k)
call printsum(f,'s%c',id,s%c)
call printsum(f,'s%cg',id,s%cg)
call printsum(f,'s%sigm',id,s%sigm)
call printsum(f,'s%hh',id,s%hh)
call printsum(f,'s%zs',id,s%zs)
call printsum(f,'s%zs0',id,s%zs0)
call printsum(f,'s%tideinpt',id,s%tideinpt)
call printsum(f,'s%tideinpz',id,s%tideinpz)
call printsum(f,'s%dzsdt',id,s%dzsdt)
call printsum(f,'s%dzbdt',id,s%dzbdt)
call printsum(f,'s%uu',id,s%uu)
call printsum(f,'s%vv',id,s%vv)
call printsum(f,'s%qx',id,s%qx)
call printsum(f,'s%qy',id,s%qy)
call printsum(f,'s%sedero',id,s%sedero)
call printsum(f,'s%dcdx',id,s%dcdx)
call printsum(f,'s%dcdy',id,s%dcdy)
call printsum(f,'s%ui',id,s%ui)
call printsum(f,'s%E',id,s%E)
call printsum(f,'s%R',id,s%R)
call printsum(f,'s%urms',id,s%urms)
call printsum(f,'s%D',id,s%D)
call printsum(f,'s%ust',id,s%ust)
call printsum(f,'s%tm',id,s%tm)
call printsum(f,'s%ueu',id,s%ueu)
call printsum(f,'s%vev',id,s%vev)
call printsum(f,'s%vmagu',id,s%vmagu)
call printsum(f,'s%vmageu',id,s%vmagu)
call printsum(f,'s%vmagv',id,s%vmagv)
call printsum(f,'s%vmagev',id,s%vmagev)
call printsum(f,'s%u',id,s%u)
call printsum(f,'s%v',id,s%v)
call printsum(f,'s%ue',id,s%ue)
call printsum(f,'s%ve',id,s%ve)
call printsum(f,'s%hold',id,s%hold)
call printsum(f,'s%wetu',id,s%wetu)
call printsum(f,'s%wetv',id,s%wetv)
call printsum(f,'s%wetz',id,s%wetz)
call printsum(f,'s%hu',id,s%hu)
call printsum(f,'s%hv',id,s%hv)
call printsum(f,'s%hum',id,s%hum)
call printsum(f,'s%hvm',id,s%hvm)
call printsum(f,'s%vmag',id,s%vmag)
call printsum(f,'s%ccg',id,s%ccg)
call printsum(f,'s%uwf',id,s%uwf)
call printsum(f,'s%vwf',id,s%vwf)
call printsum(f,'s%ustr',id,s%ustr)
call printsum(f,'s%usd',id,s%usd)
call printsum(f,'s%bi',id,s%bi)
call printsum(f,'s%DR ',id,s%DR )
call printsum(f,'s%umean ',id,s%umean )
call printsum(f,'s%vardx',id,s%vardx)
call printsum(f,'s%vu',id,s%vu)
call printsum(f,'s%uv',id,s%uv)
call printsum(f,'s%graindistr',id,s%graindistr)
call printsum(f,'s%D50',id,s%D50)
call printsum(f,'s%D90',id,s%D90)
call printsum(f,'s%sedcal',id,s%sedcal)
call printsum(f,'s%Tsg',id,s%Tsg)
call printsum(f,'s%Sug',id,s%Sug)
call printsum(f,'s%Svg',id,s%Svg)
call printsum(f,'s%ceqg',id,s%ceqg)
call printsum(f,'s%ua',id,s%ua)
call printsum(f,'s%BR',id,s%BR)
call printsum(f,'s%kb',id,s%kb)
call printsum(f,'s%Tbore',id,s%Tbore)
call printsum(f,'s%uon',id,s%uon)
call printsum(f,'s%uoff',id,s%uoff)
call printsum(f,'s%dzav',id,s%dzav)
call printsum(f,'s%maxzs',id,s%maxzs)
call printsum(f,'s%minzs',id,s%minzs)
#ifdef USEMPI
call printsum(f,'s%is',id,s%is) 
call printsum(f,'s%js',id,s%js) 
call printsum(f,'s%lm',id,s%lm) 
call printsum(f,'s%ln',id,s%ln) 
#endif
end subroutine printssums

subroutine printssumso(s)
type (spacepars), intent(in) :: s

print *,'Start of printssumso'

print *,'s%xw',sum(s%xw)
print *,'s%yw',sum(s%yw)
print *,'s%x',sum(s%x)
print *,'s%y',sum(s%y)
print *,'s%zs',sum(s%zs)
print *,'s%u',sum(s%u)
print *,'s%v',sum(s%v)
print *,'s%ue',sum(s%ue)
print *,'s%ve',sum(s%ve)
print *,'s%H',sum(s%H)
print *,'s%urms',sum(s%urms)
print *,'s%zb',sum(s%zb)
print *,'s%hh',sum(s%hh)
print *,'s%Fx',sum(s%Fx)
print *,'s%Fy',sum(s%Fy)
print *,'s%E',sum(s%E)
print *,'s%R',sum(s%R)
print *,'s%D',sum(s%D)
end subroutine printssumso
end module spaceparams
