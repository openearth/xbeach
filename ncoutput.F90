!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   MODULE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this module will provide for netcdf output in the way varoutput provides for binary output. 
! It should be called optionally through one of the params.txt settings
! It will also be compiled conditionally like mpi. Only if it proves usefull will it be added by default. 
! it will add dependencies on the netcdf fortran library (http://www.unidata.ucar.edu/software/netcdf/)
! 
module ncoutput_module
#ifdef USENETCDF

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

  use netcdf

implicit none
private
public ncoutput_init, nc_output 

integer*4                           :: nglobalvar  ! number of global output variables
integer*4                           :: npoints     ! number of output points
integer*4                           :: nrugauge    ! number of runup gauges
integer*4,dimension(:),allocatable  :: nassocvar   ! vector with number of output variable per output point
integer*4                           :: nmeanvar    ! number of time-average variables

integer, save :: ncid
character(80), save                 :: ncfilename          


! grid
integer, save :: xdimid, ydimid
integer, save :: xvarid, yvarid 

! Wave angle
integer, save :: thetadimid
! Sediment
integer, save :: sedimentclassesdimid, bedlayersdimid

! global
integer, dimension(:), allocatable, save :: globalvarids
! default output (fixed length)

! points 
integer, save :: pointsdimid
integer, save :: xpointsvarid, ypointsvarid, pointtypesvarid, xpointindexvarid, ypointindexvarid
integer, dimension(:), allocatable, save :: pointsvarids
integer, dimension(:),allocatable, save  :: xpoints     ! model x-coordinate of output points
integer, dimension(:),allocatable, save  :: ypoints     ! model y-coordinate of output points

! global
integer, dimension(:), allocatable, save :: meanvarids

! time 
integer, save :: globaltimedimid, pointtimedimid
integer, save :: globaltimevarid, pointtimevarid

! TODO: check out why these are sometimes used....
integer, save :: tidetimedimid, windtimedimid
integer, save :: inoutdimid, tidecornersdimid

! local variables
integer, save :: npointstotal
logical, save :: pointoutput




contains

  
  ! Error handling of netcdf errors 
  subroutine handle_err(status) 
    use netcdf

    integer, intent ( in) :: status
    integer :: status2
    
    if(status /= nf90_noerr) then
       !UNIT=6 for stdout and UNIT=0 for stderr.
       write(0,*) trim(nf90_strerror(status))
       write(0,*) 'closing file'
       status2 = nf90_close(ncid)
       if (status2 /= nf90_noerr) then
          write(0,*) trim(nf90_strerror(status2))
       end if
       stop 1
    end if
  end subroutine handle_err
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!   INITIALISE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ncoutput_init(s, sl, par, tpar)
  use params
  use spaceparams
  use readkey_module
  use timestep_module
  use mnemmodule
  use postprocessmod
  use logging_module
  implicit none
  integer :: status ! file id and status returned from a file operation

  type(spacepars), intent(in)                  :: s,sl ! why do we get 2 space params??
  type(parameters), intent(in)                 :: par 
  type(timepars), intent(in)                   :: tpar
  
  type(arraytype)                                     :: t
  integer                                      :: i,j,k,l,m,n
  character(len=10)                             :: mnem
  
  integer                                      :: npointstotal
  logical                                      :: outputp, outputg, outputw, outputm, outputc
  integer, dimension(:), allocatable           :: dimids ! store the dimids in a vector
  
  ! variables for point grid snapping
  real*8,dimension(s%nx+1,s%ny+1)	  :: mindist
  integer,dimension(2)                :: minlocation


  ! initialize values
  ! global
  ! store netcdf variable ids for each variable
  allocate(globalvarids(par%nglobalvar))
  allocate(meanvarids(par%nmeanvar))
  globalvarids = -1 ! initialize to -1, so an error is raised when we miss something... 
  meanvarids = -1
  outputg = .true.

  npointstotal = par%npoints+par%nrugauge
  outputp = (npointstotal .gt. 0) .and. (size(tpar%tpp) .gt. 0)
  allocate(pointsvarids(par%npoints))

  ncfilename = 'xboutput.nc' 
  ! create a file
  status = nf90_create(path = ncfilename, cmode=NF90_CLOBBER, ncid = ncid)
  if (status /= nf90_noerr) call handle_err(status)


  ! dimensions TODO: only output dimensions that are used
  ! grid
  status = nf90_def_dim(ncid, 'x', s%nx+1, xdimid) 
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, 'y', s%ny+1, ydimid)
  if (status /= nf90_noerr) call handle_err(status)
  ! wave angles
  status = nf90_def_dim(ncid, 'wave_angle', s%ntheta, thetadimid)
  if (status /= nf90_noerr) call handle_err(status)
  ! computational layers in bed ...
  ! TODO: Clean this up, why max(par%nd,2)???
  status = nf90_def_dim(ncid, 'bed_layers', max(par%nd,2), bedlayersdimid)
  if (status /= nf90_noerr) call handle_err(status)
  ! sediment classes
  status = nf90_def_dim(ncid, 'sediment_classes', par%ngd, sedimentclassesdimid)
  if (status /= nf90_noerr) call handle_err(status)

  ! dimensions of length 2.... what is this.... TODO: find out what this is 
  status = nf90_def_dim(ncid, 'inout', 2, inoutdimid)
  if (status /= nf90_noerr) call handle_err(status)

  ! time dimensions are fixed, only defined if there are points
  if (outputg) then
     status = nf90_def_dim(ncid, 'globaltime', size(tpar%tpg), globaltimedimid)
     if (status /= nf90_noerr) call handle_err(status)
  end if
  if (outputp) then
     ! points
     status = nf90_def_dim(ncid, 'points', npointstotal, pointsdimid)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_def_dim(ncid, 'pointtime', size(tpar%tpp), pointtimedimid)
     if (status /= nf90_noerr) call handle_err(status)
  end if

  ! TODO: par%tidelen par%tideloc par%windlen 
  if (par%tidelen > 0) then
     status = nf90_def_dim(ncid, 'tidetime', par%tidelen, tidetimedimid)
     if (status /= nf90_noerr) call handle_err(status)
  endif

  if (par%tideloc > 0) then
     status = nf90_def_dim(ncid, 'tidecorners', par%tideloc, tidecornersdimid)
     if (status /= nf90_noerr) call handle_err(status)
  endif

  if (par%windlen > 0) then
     status = nf90_def_dim(ncid, 'windtime', par%windlen, windtimedimid)
     if (status /= nf90_noerr) call handle_err(status)
  endif

  ! define space & time variables
  ! grid
  status = nf90_def_var(ncid, 'x', NF90_DOUBLE, (/ xdimid /), xvarid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, xvarid, 'units', 'm')
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, xvarid, 'long_name', 'local x coordinate')
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_var(ncid, 'y', NF90_DOUBLE, (/ ydimid /), yvarid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, yvarid, 'units', 'm')
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, yvarid, 'long_name', 'local y coordinate')
  if (status /= nf90_noerr) call handle_err(status)

  

  ! global
  if (outputg) then
     status = nf90_def_var(ncid, 'globaltime', NF90_DOUBLE, (/ globaltimedimid /), globaltimevarid)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_att(ncid, globaltimevarid, 'units', trim(par%tunits))
     if (status /= nf90_noerr) call handle_err(status)
     
     ! default global output variables
     do i=1,par%nglobalvar
        mnem = trim(par%globalvars(i))
        j = chartoindex(mnem)
        call indextos(s,j,t)
        select case(t%type)
        case('r')
           ! Build the array with dimension ids
           call writelog('ls', '', 'Creating netcdf variable: ', trim(mnem) )
           allocate(dimids(t%rank+1))
           select case(t%rank)
           case(2)
              dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), globaltimedimid /)
           case(3)
              dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), &
                   dimensionid(t%dimensions(3)), globaltimedimid /)
           case(4)
              call writelog('ls', '', 'Variable ' // trim(mnem) // ' is of rank 4. This may not work due to an' // &
                   ' unresolved issue. If so, remove the variable or use the fortran outputformat option.')
              dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), &
                   dimensionid(t%dimensions(3)), dimensionid(t%dimensions(4)), globaltimedimid /)
           case default
              call writelog('lse', '', 'mnem: ' // mnem // ' not supported, rank:', t%rank)
              stop 1
           end select
           status = nf90_def_var(ncid, trim(mnem), NF90_DOUBLE, &
                dimids, globalvarids(i))
           if (status /= nf90_noerr) call handle_err(status)
           deallocate(dimids)
        case default
           write(0,*) 'mnem', mnem, ' not supported, type:', t%type
        end select
        status = nf90_put_att(ncid, globalvarids(i), 'units', trim(t%units))
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_put_att(ncid, globalvarids(i), 'long_name', trim(t%description))
        if (status /= nf90_noerr) call handle_err(status)
     end do
  end if

!  ! points
  ! default global output variables
  if (outputp) then
     status = nf90_def_var(ncid, 'pointtime', NF90_DOUBLE, (/ pointtimedimid /), pointtimevarid)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_att(ncid, pointtimevarid, 'units', trim(par%tunits))
     if (status /= nf90_noerr) call handle_err(status)

     ! points
     status = nf90_def_var(ncid, 'xpoint', NF90_DOUBLE, (/ pointsdimid /), xpointsvarid)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_att(ncid, xpointsvarid, 'units', 'm')
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_att(ncid, xpointsvarid, 'long_name', 'local x coordinate')
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_def_var(ncid, 'ypoint', NF90_DOUBLE, (/ pointsdimid /), ypointsvarid)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_att(ncid, ypointsvarid, 'units', 'm')
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_att(ncid, ypointsvarid, 'long_name', 'local y coordinate')
     if (status /= nf90_noerr) call handle_err(status)
     
     status = nf90_def_var(ncid, 'xpointindex', NF90_DOUBLE, (/ pointsdimid /), xpointindexvarid)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_att(ncid, xpointindexvarid, 'long_name', 'nearest x grid cell')
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_def_var(ncid, 'ypointindex', NF90_DOUBLE, (/ pointsdimid /), ypointindexvarid)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_att(ncid, ypointindexvarid, 'long_name', 'nearest y grid cell')
     if (status /= nf90_noerr) call handle_err(status)
     
     status = nf90_def_var(ncid, 'pointtypes', NF90_DOUBLE, (/ pointsdimid /), pointtypesvarid)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_att(ncid, pointtypesvarid, 'long_name', 'type of point (0=point, 1=rugauge)')
     if (status /= nf90_noerr) call handle_err(status)

     do i=1,par%npointvar
        mnem = trim(par%pointvars(i))
        j = chartoindex(mnem)
        call indextos(s,j,t)
        select case(t%type)
        case('r')
           ! Build the array with dimension ids
           call writelog('ls', '', 'Creating netcdf variable: ', 'point_'// trim(mnem) )
           allocate(dimids(t%rank))
           ! Make sure the variable has x and y as the first 2 dimensions
           if ((dimensionid(t%dimensions(1)) .ne. xdimid) .or. (dimensionid(t%dimensions(2)) .ne. ydimid)) then
              call writelog('lse','', 'Tried to store variable ' // trim(mnem) // ', but this variable is not a function of x,y')
           endif
           select case(t%rank)
           case(2)
              dimids = (/ pointsdimid, pointtimedimid /)
           case(3)
              dimids = (/ pointsdimid, dimensionid(t%dimensions(3)), pointtimedimid /)
           case(4)
              call writelog('ls', '', 'Variable ' // trim(mnem) // ' is of rank 4. This may not work due to an' // &
                   ' unresolved issue. If so, remove the variable or use the fortran outputformat option.')
              dimids = (/ pointsdimid, dimensionid(t%dimensions(3)), dimensionid(t%dimensions(4)), pointtimedimid /)
           case default
              call writelog('lse', '', 'mnem: ' // mnem // ' not supported, rank:', t%rank)
              stop 1
           end select
           status = nf90_def_var(ncid, 'point_' // trim(mnem), NF90_DOUBLE, &
                dimids, pointsvarids(i))
           if (status /= nf90_noerr) call handle_err(status)
           deallocate(dimids)
        case default
           write(0,*) 'mnem', mnem, ' not supported, type:', t%type
        end select
        status = nf90_put_att(ncid, pointsvarids(i), 'units', trim(t%units))
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_put_att(ncid, pointsvarids(i), 'long_name', trim(t%description))
        if (status /= nf90_noerr) call handle_err(status)
     end do
  end if
  
  ! done defining variables
  call writelog('ls', '', 'Writing file definition.')
  status = nf90_enddef(ncid)
  if (status /= nf90_noerr) call handle_err(status)

  ! Fill meta variables
  ! Grid
  j = chartoindex('xz')
  call indextos(s,j,t)
  status = nf90_put_var(ncid, xvarid, t%r1)
  if (status /= nf90_noerr) call handle_err(status)

  j = chartoindex('yz')
  call indextos(s,j,t)
  status = nf90_put_var(ncid, yvarid, t%r1)
  if (status /= nf90_noerr) call handle_err(status)
  
  if (outputp) then
     call writelog('ls', '', 'Writing point vars.')
     status = nf90_put_var(ncid, xpointsvarid, par%xpointsw)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_var(ncid, ypointsvarid, par%ypointsw)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_var(ncid, pointtypesvarid, par%pointtypes)
     if (status /= nf90_noerr) call handle_err(status)

     ! Convert world coordinates of points to nearest (lsm) grid point
     ! This could be done in some postprocessing function
     allocate(xpoints(par%nrugauge+par%npoints))
     allocate(ypoints(par%nrugauge+par%npoints))
     call snappointstogrid(par, s, xpoints, ypoints)

     status = nf90_put_var(ncid, xpointindexvarid, xpoints)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_var(ncid, ypointindexvarid, ypoints)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_var(ncid, pointtypesvarid, par%pointtypes)
     if (status /= nf90_noerr) call handle_err(status)
  end if

  status = nf90_close(ncid)
  if (status /= nf90_noerr) call handle_err(status)
  
end subroutine ncoutput_init

subroutine nc_output(s,sl,par, tpar)
  use params
  use spaceparams
  use timestep_module
  use mnemmodule
  use postprocessmod

  implicit none

  type(spacepars), intent(in)            :: s,sl ! s-> spaceparams, what is isl?
  type(parameters), intent(in)           :: par
  type(timepars), intent(in)             :: tpar
  
  type(arraytype)                        :: t
  integer                                :: i,j,k,l,m,n, ii
  character(len=10)                      :: mnem

  ! some local variables to pass the data through the postprocessing function.
  integer :: i0
  integer, dimension(:), allocatable :: i1
  integer, dimension(:,:), allocatable :: i2
  integer, dimension(:,:,:), allocatable :: i3
  integer, dimension(:,:,:,:), allocatable :: i4
  real*8 :: r0
  real*8, dimension(:), allocatable :: r1
  real*8, dimension(:,:), allocatable :: r2
  real*8, dimension(:,:,:), allocatable :: r3
  real*8, dimension(:,:,:,:), allocatable :: r4

  integer :: status

  ! not a output step

  status = nf90_open(ncid=ncid, path=ncfilename, mode=nf90_write)
  if (status /= nf90_noerr) call handle_err(status)  


  if (tpar%outputg) then
     status = nf90_put_var(ncid, globaltimevarid, par%t, (/tpar%itg/))
     if (status /= nf90_noerr) call handle_err(status) 
     ! write global output variables
     do i=1,par%nglobalvar
        mnem = trim(par%globalvars(i))
        j = chartoindex(mnem)
        ! lookup the proper array
        call indextos(s,j,t)
        select case(t%type)
        case('r')
           select case(t%rank)
           case(2)
              allocate(r2(size(t%r2,1),size(t%r2,2)))
              call gridrotate(s, t, r2)
              status = nf90_put_var(ncid, globalvarids(i), r2, start=(/1,1,tpar%itg/) )
              deallocate(r2)
              if (status /= nf90_noerr) call handle_err(status)
           case(3)
              allocate(r3(size(t%r3,1),size(t%r3,2),size(t%r3,3)))
              call gridrotate(s, t, r3)
              status = nf90_put_var(ncid, globalvarids(i), r3, start=(/1,1,1, tpar%itg/) )
              deallocate(r3)
              if (status /= nf90_noerr) call handle_err(status)
           case(4)
              allocate(r4(size(t%r4,1),size(t%r4,2),size(t%r4,3),size(t%r4,4)))
              call gridrotate(s, t, r4)
              status = nf90_put_var(ncid, globalvarids(i), r4, start=(/1,1,1,1, tpar%itg/) )
              deallocate(r4)
              if (status /= nf90_noerr) call handle_err(status)
           case default
              write(0,*) 'Can''t handle rank: ', t%rank, ' of mnemonic', mnem
           end select
        case default
           write(0,*) 'Can''t handle type: ', t%type, ' of mnemonic', mnem
        end select
     end do
!  tpar%outputg = .false. ! not sure if this is required
  end if
  if (tpar%outputp) then

     status = nf90_put_var(ncid, pointtimevarid, par%t, (/tpar%itp/))
     if (status /= nf90_noerr) call handle_err(status) 

     do i=1,par%npointvar
        mnem = trim(par%pointvars(i))
        j = chartoindex(mnem)
        ! lookup the proper array
        call indextos(s,j,t)
        ! get the proper output points ....
        ! I have no idea what is happening in varouput so I'll try it in a different way
        select case(t%type)
           !TODO This is not very efficient because we are using the outer counters, reorder dimensions....
        case('r')
           do ii = 1, (par%npoints + par%nrugauge)           
              select case(t%rank)
              case(2)
                 ! This postprocessing creates an ugly dependency.
                 ! it would be nice if we could call gridrotate as a function
                 ! or if we could just have the postprocessing insert some reference processing routines to call
                 ! or if we could split this out of the case statement (dry)
                 ! or if we could defer this to a postprocessing routine (for example ncks)
                 allocate(r2(size(t%r2,1),size(t%r2,2)))
                 call gridrotate(s, t, r2)
                 status = nf90_put_var(ncid, pointsvarids(i), r2(xpoints(ii), ypoints(ii)), start=(/ii,tpar%itp/) )
                 deallocate(r2)
                 if (status /= nf90_noerr) call handle_err(status)
              case(3)
                 allocate(r3(size(t%r3,1),size(t%r3,2),size(t%r3,3)))
                 call gridrotate(s, t, r3)
                 status = nf90_put_var(ncid, pointsvarids(i), r3(xpoints(ii), ypoints(ii),:), start=(/ii,1,tpar%itp/) )
                 deallocate(r3)
                 if (status /= nf90_noerr) call handle_err(status)
              case(4)
                 allocate(r4(size(t%r4,1),size(t%r4,2),size(t%r4,3),size(t%r4,4)))
                 call gridrotate(s, t, r4)
                 status = nf90_put_var(ncid, pointsvarids(i), r4(xpoints(ii), ypoints(ii),:,:), start=(/ii,1,1,tpar%itp/) )
                 deallocate(r4)
                 if (status /= nf90_noerr) call handle_err(status)
              case default
                 write(0,*) 'Can''t handle rank: ', t%rank, ' of mnemonic', mnem
              end select
           end do
           case default
              write(0,*) 'Can''t handle type: ', t%type, ' of mnemonic', mnem
           end select
        end do
  end if
  status = nf90_close(ncid=ncid)
  if (status /= nf90_noerr) call handle_err(status) 
 

end subroutine nc_output


integer function dimensionid(expression)
  ! Function to transform the expression in spaceparams.tmpl to an id, we might want this in the 
  ! makeincludes module
  use logging_module
  implicit none
  character(len=20),intent(in)                       :: expression
  integer :: i, ic
  select case(trim(expression))
  case('s%nx+1')
     dimensionid = xdimid
  case('s%ny+1')
     dimensionid = ydimid
  case('s%ntheta')
     dimensionid = thetadimid
  case('par%tidelen')
     dimensionid = tidetimedimid
  case('par%tideloc')
     dimensionid = tidecornersdimid
  case('par%windlen')
     dimensionid = windtimedimid
  case('par%ngd')
     dimensionid = sedimentclassesdimid
  case('2')
     dimensionid = inoutdimid
  case('max(par%nd,2)')
     dimensionid = bedlayersdimid
  case default
     call writelog('els','','Unknown dimension expression:'  // expression)
     stop 1
  end select
end function
#endif
end module ncoutput_module
