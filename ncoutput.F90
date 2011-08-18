!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   MODULE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this module will provide for netcdf output in the way varoutput provides for binary output. 
! It should be called optionally through one of the params.txt settings
! It will also be compiled conditionally like mpi. Only if it proves usefull will it be added by default. 
! it will add dependencies on the netcdf fortran library (http://www.unidata.ucar.edu/software/netcdf/)
! 
! With contributions from Uwe Rosebrock (CSIRO).
module ncoutput_module
#ifdef USENETCDF

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
  use xmpi_module
  use netcdf

  implicit none
  private
  public ncoutput_init, ncoutput 

  integer, save :: ncid

  ! parameters
  integer, save :: parvarid

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

  ! mean
  ! number of variables by number of parameters per variable (mean, sigma^2, min, max)
  integer, dimension(:,:), allocatable, save       :: meanvarids
  character(8), dimension(:), allocatable, save    :: meanvartypes
  integer*4                           :: nmeanvartypes  = 4   ! number of time-average variable types



  ! time 
  integer, save :: globaltimedimid, pointtimedimid, meantimedimid
  integer, save :: globaltimevarid, pointtimevarid, meantimevarid

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
    use xmpi_module
    use params
    use spaceparams
    use timestep_module
    use mnemmodule
    use means_module
    use postprocessmod
    use logging_module
    ! This module is awaiting comments from Robert.
    use getkey_module
    
    implicit none
    integer :: status ! file id and status returned from a file operation

    type(spacepars), intent(in)                  :: s,sl ! Use s => global data and sl => local data
    type(parameters), intent(in)                 :: par 
    type(timepars), intent(in)                   :: tpar

    ! Part of the getkey
    type(parameter)                              :: val
    type(arraytype)                              :: t
    type(meanspars)                              :: meanvar
    integer                                      :: i,j
    character(len=maxnamelen)                    :: mnem

    integer                                      :: npointstotal
    logical                                      :: outputp, outputg, outputm
    integer, dimension(:), allocatable           :: dimids ! store the dimids in a vector
    character(256)                               :: coordinates
    character(8)                                 :: cellmethod

    character(len=maxnamelen), dimension(:), allocatable       :: keys
    ! subversion information
    include 'version.def'
    include 'version.dat'

    if (xmaster) then

       ! initialize values
       ! global

       ! store netcdf variable ids for each variable
       allocate(globalvarids(par%nglobalvar))
       globalvarids = -1 ! initialize to -1, so an error is raised when we miss something... 
       outputg = .true.

       npointstotal = par%npoints+par%nrugauge
       outputp = (npointstotal .gt. 0) .and. (size(tpar%tpp) .gt. 0)
       allocate(pointsvarids(par%npointvar))


       allocate(meanvarids(par%nmeanvar,nmeanvartypes))
       meanvarids = -1
       outputm = (par%nmeanvar .gt. 0)
       allocate(meanvartypes(nmeanvartypes))
       meanvartypes = (/ 'mean    ', 'var     ', 'min     ', 'max     '  /)


       ! create a file
       status = nf90_create(path = par%ncfilename, cmode=ior(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid = ncid)
       if (status /= nf90_noerr) call handle_err(status)

       ! dimensions TODO: only output dimensions that are used
       ! grid
       status = nf90_def_dim(ncid, 'globalx', s%nx+1, xdimid) 
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_def_dim(ncid, 'globaly', s%ny+1, ydimid)
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
          status = nf90_def_dim(ncid, 'globaltime', NF90_unlimited, globaltimedimid)
          if (status /= nf90_noerr) call handle_err(status)
       end if
       if (outputp) then
          ! points
          status = nf90_def_dim(ncid, 'points', npointstotal, pointsdimid)
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_def_dim(ncid, 'pointtime', size(tpar%tpp), pointtimedimid)
          if (status /= nf90_noerr) call handle_err(status)
       end if
       if (outputm) then
          status = nf90_def_dim(ncid, 'meantime', size(tpar%tpm)-1, meantimedimid)
          if (status /= nf90_noerr) call handle_err(status)
       end if

       if (s%tidelen > 0) then
          status = nf90_def_dim(ncid, 'tidetime', s%tidelen, tidetimedimid)
          if (status /= nf90_noerr) call handle_err(status)
       endif

       if (par%tideloc > 0) then
          status = nf90_def_dim(ncid, 'tidecorners', par%tideloc, tidecornersdimid)
          if (status /= nf90_noerr) call handle_err(status)
       endif

       if (s%windlen > 0) then
          status = nf90_def_dim(ncid, 'windtime', s%windlen, windtimedimid)
          if (status /= nf90_noerr) call handle_err(status)
       endif

       ! define empty parameter variable
       status = nf90_def_var(ncid, 'parameter', NF90_DOUBLE, varid=parvarid)
       if (status /= nf90_noerr) call handle_err(status)
       ! define space & time variables
       ! grid
       status = nf90_def_var(ncid, 'globalx', NF90_DOUBLE, (/ xdimid, ydimid /), xvarid)
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid, xvarid, 'units', 'm')
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid, xvarid, 'long_name', 'local x coordinate')
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid, xvarid, 'standard_name', 'projection_x_coordinate')
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid, xvarid, 'axis', 'X')
       if (status /= nf90_noerr) call handle_err(status)
       ! For compatibility with CSIRO Dive software
       if (len(trim(par%projection)) .ne. 0)  then
          status = nf90_put_att(ncid, xvarid, 'projection', par%projection)
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, xvarid, 'rotation',( s%alfa/atan(1.0d0)*45.d0))
          if (status /= nf90_noerr) call handle_err(status)
       end if

       status = nf90_def_var(ncid, 'globaly', NF90_DOUBLE, (/ xdimid, ydimid /), yvarid)
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid, yvarid, 'units', 'm')
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid, yvarid, 'long_name', 'local y coordinate')
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid, yvarid, 'standard_name', 'projection_y_coordinate')
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid, yvarid, 'axis', 'Y')
       if (status /= nf90_noerr) call handle_err(status)
       ! For compatibility with CSIRO Dive software
       if (len(trim(par%projection)) .ne. 0)  then
          status = nf90_put_att(ncid, yvarid, 'projection', par%projection)
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, yvarid, 'rotation',( s%alfa/atan(1.0d0)*45.d0))
          if (status /= nf90_noerr) call handle_err(status)
       end if

       ! Some metadata attributes
       status = nf90_put_att(ncid,nf90_global, "Conventions", "CF-1.4")
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid,nf90_global, "Producer", "XBeach littoral zone wave model (http://www.xbeach.org)")
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid,nf90_global, "Build-Revision", trim(Build_Revision))
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid,nf90_global, "Build-Date", trim(Build_Date))
       if (status /= nf90_noerr) call handle_err(status)
       status = nf90_put_att(ncid,nf90_global, "URL", trim(Build_URL))
       if (status /= nf90_noerr) call handle_err(status)

       ! Store all the parameters
       ! This part is awaiting comments from Robert McCall
       call getkeys(par, keys)
       do i=1,size(keys)
          call getkey(par, keys(i), val)
          if (val%type == 'i') then
             status = nf90_put_att(ncid, parvarid, keys(i), val%i0 )
          elseif (val%type == 'c') then
             status = nf90_put_att(ncid, parvarid, keys(i), val%c0 )
          elseif (val%type == 'r') then
             status = nf90_put_att(ncid, parvarid, keys(i), val%r0 )
          end if
          if (status /= nf90_noerr) call handle_err(status)
       end do

       ! global
       if (outputg) then
          status = nf90_def_var(ncid, 'globaltime', NF90_DOUBLE, (/ globaltimedimid /), globaltimevarid)
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, globaltimevarid, 'units', trim(par%tunits))
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, globaltimevarid, 'axis', 'T')
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, globaltimevarid, 'standard_name', 'time')
          if (status /= nf90_noerr) call handle_err(status)

          ! default global output variables
          do i=1,par%nglobalvar
             mnem = trim(par%globalvars(i))
             coordinates = ''
             j = chartoindex(mnem)
             call indextos(s,j,t)
             call writelog('ls', '', 'Creating netcdf variable: ', trim(mnem) )

             ! Build the array with dimension ids
             allocate(dimids(t%rank+1))
             select case(t%rank)
             case(0)
                dimids = (/ globaltimedimid /)
                coordinates = ''
             case(1)
                dimids = (/ dimensionid(t%dimensions(1)), globaltimedimid /)
                coordinates = ''
                if (dimids(1) .eq. xdimid) coordinates = 'globalx'
                if (dimids(1) .eq. ydimid) coordinates = 'globaly'
             case(2)
                dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), globaltimedimid /)
                coordinates = 'globalx globaly'
             case(3)
                dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), &
                     dimensionid(t%dimensions(3)), globaltimedimid /)
                coordinates = 'globalx globaly' 
                ! Do we have a vertical level?
                if (dimids(3) .eq. bedlayersdimid) coordinates = trim(coordinates) // ' bed_layers'
             case(4)
                call writelog('ls', '', 'Variable ' // trim(mnem) // ' is of rank 4. This may not work due to an' // &
                     ' unresolved issue. If so, remove the variable or use the fortran outputformat option.')
                dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), &
                     dimensionid(t%dimensions(3)), dimensionid(t%dimensions(4)), globaltimedimid /)
                coordinates = 'globalx globaly' 
                ! Do we have a vertical level?
                if ((dimids(3) .eq. bedlayersdimid) .or. (dimids(4) .eq. bedlayersdimid)) then
                   coordinates = trim(coordinates) // ' bed_layers'
                end if
             case default
                call writelog('lse', '', 'mnem: ' // mnem // ' not supported, rank:', t%rank)
                stop 1
             end select
             select case(t%type)
             case('i')
                status = nf90_def_var(ncid, trim(mnem), NF90_INT, &
                     dimids, globalvarids(i))
                if (status /= nf90_noerr) call handle_err(status)
             case('r')
                status = nf90_def_var(ncid, trim(mnem), NF90_DOUBLE, &
                     dimids, globalvarids(i))
                if (status /= nf90_noerr) call handle_err(status)
             case default
                write(0,*) 'mnem', mnem, ' not supported, type:', t%type
             end select
             status = nf90_put_att(ncid, globalvarids(i), 'coordinates', trim(coordinates))
             if (status /= nf90_noerr) call handle_err(status)
             deallocate(dimids)
             status = nf90_put_att(ncid, globalvarids(i), 'units', trim(t%units))
             if (status /= nf90_noerr) call handle_err(status)
             if (.not.(trim(t%standardname) .eq. '')) then
               status = nf90_put_att(ncid, globalvarids(i), 'standard_name', trim(t%standardname))
               if (status /= nf90_noerr) call handle_err(status)
             endif
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
          status = nf90_put_att(ncid, pointtimevarid, 'axis', 'T')
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, pointtimevarid, 'standard_name', 'time')
          if (status /= nf90_noerr) call handle_err(status)

          ! points
          status = nf90_def_var(ncid, 'pointx', NF90_DOUBLE, (/ pointsdimid /), xpointsvarid)
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, xpointsvarid, 'units', 'm')
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, xpointsvarid, 'long_name', 'local x coordinate')
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, xpointsvarid, 'standard_name', 'projection_x_coordinate')
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, xpointsvarid, 'axis', 'X')
          if (status /= nf90_noerr) call handle_err(status)


          status = nf90_def_var(ncid, 'pointy', NF90_DOUBLE, (/ pointsdimid /), ypointsvarid)
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, ypointsvarid, 'units', 'm')
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, ypointsvarid, 'long_name', 'local y coordinate')
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, ypointsvarid, 'standard_name', 'projection_y_coordinate')
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, ypointsvarid, 'axis', 'Y')
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
             coordinates = ''
             call indextos(s,j,t)
             select case(t%type)
             case('r')
                ! Build the array with dimension ids
                call writelog('ls', '', 'Creating netcdf variable: ', 'point_'// trim(mnem) )
                allocate(dimids(t%rank))
                ! Make sure the variable has x and y as the first 2 dimensions
                if ((dimensionid(t%dimensions(1)) .ne. xdimid) .or. (dimensionid(t%dimensions(2)) .ne. ydimid)) then
                   call writelog('lse','', 'Tried to store variable ' // trim(mnem) // ', but it is not a function of x,y')
                endif
                select case(t%rank)
                case(2)
                   dimids = (/ pointsdimid, pointtimedimid /)
                   coordinates = 'pointx pointy'
                case(3)
                   dimids = (/ pointsdimid, dimensionid(t%dimensions(3)), pointtimedimid /)
                   coordinates = 'pointx pointy'
                   ! Do we have a vertical level?
                   if (dimids(3) .eq. bedlayersdimid) coordinates = trim(coordinates) // ' bed_layers'
                case(4)
                   call writelog('ls', '', 'Variable ' // trim(mnem) // ' is of rank 4. This may not work due to an' // &
                        ' unresolved issue. If so, remove the variable or use the fortran outputformat option.')
                   dimids = (/ pointsdimid, dimensionid(t%dimensions(3)), dimensionid(t%dimensions(4)), pointtimedimid /)
                   coordinates = 'pointx pointy'
                   ! Do we have a vertical level?
                   if ((dimids(3) .eq. bedlayersdimid) .or. (dimids(4) .eq. bedlayersdimid)) then
                      coordinates = trim(coordinates) // ' bed_layers'
                   end if
                case default
                   call writelog('lse', '', 'mnem: ' // mnem // ' not supported, rank:', t%rank)
                   stop 1
                end select
                status = nf90_def_var(ncid, 'point_' // trim(mnem), NF90_DOUBLE, &
                     dimids, pointsvarids(i))
                if (status /= nf90_noerr) call handle_err(status)
                status = nf90_put_att(ncid, pointsvarids(i), 'coordinates', trim(coordinates))
                if (status /= nf90_noerr) call handle_err(status)
                deallocate(dimids)
             case default
                write(0,*) 'mnem', mnem, ' not supported, type:', t%type
             end select
             status = nf90_put_att(ncid, pointsvarids(i), 'units', trim(t%units))
             if (status /= nf90_noerr) call handle_err(status)
             if (.not.(trim(t%standardname) .eq. '')) then
               status = nf90_put_att(ncid, pointsvarids(i), 'standard_name', trim(t%standardname))
               if (status /= nf90_noerr) call handle_err(status)
             endif
             status = nf90_put_att(ncid, pointsvarids(i), 'long_name', trim(t%description))
             if (status /= nf90_noerr) call handle_err(status)
          end do
       end if

       if (outputm) then
          status = nf90_def_var(ncid, 'meantime', NF90_DOUBLE, (/ meantimedimid /), meantimevarid)
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, meantimevarid, 'units', trim(par%tunits))
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, meantimevarid, 'axis', 'T')
          if (status /= nf90_noerr) call handle_err(status)
          status = nf90_put_att(ncid, meantimevarid, 'standard_name', 'time')
          if (status /= nf90_noerr) call handle_err(status)
          ! default global output variables
          do i=1,par%nmeanvar
             ! Not sure if this is required here, but it is used in varoutput
! #ifdef USEMPI
!              ! No need to collect here, we're just using the types
!              ! call means_collect(sl,meansparsglobal(i),meansparslocal(i))
! #else
!              meansparsglobal(i)=meansparslocal(i)
! #endif
             coordinates = ''
             meanvar = meansparsglobal(i)
             t = meanvar%t
             select case(t%type)
             case('r')
                ! Build the array with dimension ids
                allocate(dimids(t%rank+1))
                select case(t%rank)
                case(2)
                   dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), meantimedimid /)
                   coordinates = 'globalx globaly'
                case(3)
                   dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), &
                        dimensionid(t%dimensions(3)), meantimedimid /)
                   coordinates = 'globalx globaly' 
                   ! Do we have a vertical level?
                   if (dimids(3) .eq. bedlayersdimid) coordinates = trim(coordinates) // ' bed_layers'
                case(4)
                   call writelog('ls', '', 'Variable ' // trim(mnem) // ' is of rank 4. This may not work due to an' // &
                        ' unresolved issue. If so, remove the variable or use the fortran outputformat option.')
                   dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), &
                        dimensionid(t%dimensions(3)), dimensionid(t%dimensions(4)), meantimedimid /)
                   coordinates = 'globalx globaly' 
                   ! Do we have a vertical level?
                   if ((dimids(3) .eq. bedlayersdimid) .or. (dimids(4) .eq. bedlayersdimid)) then
                      coordinates = trim(coordinates) // ' bed_layers'
                   end if
                case default
                   call writelog('lse', '', 'mnem: ' // mnem // ' not supported, rank:', t%rank)
                   stop 1
                end select

                ! Create a variable for all types of meanvars (mean, var, min, max)
                do j = 1,nmeanvartypes
                   cellmethod = trim(meanvartypes(j))
                   call writelog('ls', '', 'Creating netcdf variable: ',  trim(t%name) // '_' // cellmethod)
                   status = nf90_def_var(ncid, trim(t%name) // '_' // trim(cellmethod), NF90_DOUBLE, &
                        dimids, meanvarids(i,j))
                   if (status /= nf90_noerr) call handle_err(status)
                   status = nf90_put_att(ncid, meanvarids(i,j), 'coordinates', trim(coordinates))
                   if (status /= nf90_noerr) call handle_err(status)
                   status = nf90_put_att(ncid, meanvarids(i,j), 'units', trim(t%units))
                   if (status /= nf90_noerr) call handle_err(status)
                   if (.not.(trim(t%standardname) .eq. '')) then
                     status = nf90_put_att(ncid, meanvarids(i,j), 'standard_name', trim(t%standardname))
                     if (status /= nf90_noerr) call handle_err(status)
                   endif
                   status = nf90_put_att(ncid, meanvarids(i,j), 'long_name', trim(t%description))
                   if (status /= nf90_noerr) call handle_err(status)
                   ! For H and urms we don't compute the mean but the rms of the rms.....
                   if (cellmethod .eq. 'mean' .and. ((t%name .eq. 'H') .or. (t%name .eq. 'urms')))  then
                      cellmethod = 'rms'
                   end if
                   status = nf90_put_att(ncid, meanvarids(i,j), 'cell_methods', 'meantime: ' // trim(cellmethod))
                   if (status /= nf90_noerr) call handle_err(status)
                end do

                deallocate(dimids)
             case default
                write(0,*) 'mnem', mnem, ' not supported, type:', t%type
             end select
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
       
       status = nf90_put_var(ncid, xvarid, t%r2)
       if (status /= nf90_noerr) call handle_err(status)

       j = chartoindex('yz')
       call indextos(s,j,t)
       status = nf90_put_var(ncid, yvarid, t%r2)
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
          allocate(xpoints(npointstotal))
          allocate(ypoints(npointstotal))
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
    end if
  end subroutine ncoutput_init

  subroutine ncoutput(s,sl,par, tpar)
    use logging_module
    use xmpi_module
    use params
    use spaceparams
    use timestep_module
    use mnemmodule
    use means_module
    use postprocessmod

    implicit none

    type(spacepars), intent(in)            :: s,sl ! s-> spaceparams, what is isl?
    type(parameters), intent(in)           :: par
    type(timepars), intent(in)             :: tpar

    type(arraytype)                        :: t
    integer                                :: i,j,ii
    character(len=maxnamelen)              :: mnem

    ! some local variables to pass the data through the postprocessing function.
    integer :: i0
    integer, dimension(:,:), allocatable :: i2
    integer, dimension(:,:,:), allocatable :: i3
    real*8 :: r0
    real*8, dimension(:), allocatable :: r1
    real*8, dimension(:,:), allocatable :: r2
    real*8, dimension(:,:,:), allocatable :: r3
    real*8, dimension(:,:,:,:), allocatable :: r4

    integer :: status

    ! Open the output file
    if (xmaster) then
       status = nf90_open(ncid=ncid, path=par%ncfilename, mode=nf90_write)
       if (status /= nf90_noerr) call handle_err(status)  
    end if


    ! If we're gonna write some global output
    if (tpar%outputg) then
#ifdef USEMPI
       ! we'll need to collect the information from all nodes.
       do i=1,par%nglobalvar
          mnem = trim(par%globalvars(i))
          j = chartoindex(mnem)
          call space_collect_index(s,sl,j)
       end do
#endif
       ! only write the information on the xmaster node
       if (xmaster) then
          ! Store the time (in morphological time)
          status = nf90_put_var(ncid, globaltimevarid, par%t*max(par%morfac,1.d0), (/tpar%itg/))
          if (status /= nf90_noerr) call handle_err(status) 
          ! write global output variables
          do i=1,par%nglobalvar
             mnem = trim(par%globalvars(i))
             j = chartoindex(mnem)
             ! lookup the proper array (should have been collected already)
             call indextos(s,j,t)

             select case(t%type)
             case('i')
                select case(t%rank)
                case(0)
                   ! no need to allocate here
                   status = nf90_put_var(ncid, globalvarids(i), i0, start=(/1,tpar%itg/) )
                   if (status /= nf90_noerr) call handle_err(status)
                case(2)
                   allocate(i2(size(t%i2,1),size(t%i2,2)))
                   call gridrotate(s, t, i2)
                   status = nf90_put_var(ncid, globalvarids(i), i2, start=(/1,1,tpar%itg/) )
                   deallocate(i2)
                   if (status /= nf90_noerr) call handle_err(status)
                case(3)
                   allocate(i3(size(t%i3,1),size(t%i3,2),size(t%i3,3)))
                   call gridrotate(s, t, i3)
                   status = nf90_put_var(ncid, globalvarids(i), i3, start=(/1,1,1,tpar%itg/) )
                   deallocate(i3)
                   if (status /= nf90_noerr) call handle_err(status)
                case default
                   write(0,*) 'Can''t handle rank: ', t%rank, ' of mnemonic', mnem
                end select
             case('r')
                select case(t%rank)
                case(0)
                   ! no need to allocate here
                   status = nf90_put_var(ncid, globalvarids(i), r0, start=(/1,tpar%itg/) )
                   if (status /= nf90_noerr) call handle_err(status)
                case(1)
                   allocate(r1(size(t%r1,1)))
                   ! no need to rotate here
                   r1 = t%r1
                   status = nf90_put_var(ncid, globalvarids(i), r1, start=(/1,tpar%itg/) )
                   deallocate(r1)
                   if (status /= nf90_noerr) call handle_err(status)
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
       end if
    end if


    if (tpar%outputp) then
#ifdef USEMPI
       ! Collect all data for which we store the points.
       ! TODO: This will be a lot faster if nodes write their own point. Use the parallel netcdf for that.
       ! Let's wait till someone needs it...
       do i=1,par%npointvar
          mnem = trim(par%pointvars(i))
          j = chartoindex(mnem)
          call space_collect_index(s,sl,j)
       end do
#endif
       if (xmaster) then
          status = nf90_put_var(ncid, pointtimevarid, par%t*max(par%morfac,1.d0), (/tpar%itp/))
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
    end if



    ! If we're gonna write some mean output
    if (tpar%outputm .and. tpar%itm>1) then
       ! only write the information on the xmaster node

       do i=1,par%nmeanvar
#ifdef USEMPI
          call means_collect(sl,meansparsglobal(i),meansparslocal(i))
#else
          meansparsglobal(i)=meansparslocal(i)
#endif
       end do
       if (xmaster) then
          ! Store the time (in morphological time)
          status = nf90_put_var(ncid, meantimevarid, par%t*max(par%morfac,1.d0), (/tpar%itm-1/))
          if (status /= nf90_noerr) call handle_err(status) 
          ! write global output variables
          do i=1,par%nmeanvar
             t = meansparsglobal(i)%t
             do j=1,nmeanvartypes

                select case(t%type)
                case('r')
                   select case(t%rank)
                   case(2)
                      select case(meanvartypes(j))
                      case('mean')
                         if ((t%name .eq. 'H') .or. (t%name .eq. 'urms'))  then
                            status = nf90_put_var(ncid, meanvarids(i,j), sqrt(meansparsglobal(i)%variancesquareterm2d), &
                                 start=(/1,1,tpar%itm-1/) )
                         else
                            status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%mean2d, start=(/1,1,tpar%itm-1/) )
                         end if
                         if (status /= nf90_noerr) call handle_err(status)
                      case('var')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%variance2d, start=(/1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case('min')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%min2d, start=(/1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case('max')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%max2d, start=(/1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case default
                         write(0,*) 'Can''t handle cell method: ', trim(meanvartypes(j)), ' of mnemonic', trim(t%name)
                      end select
                   case(3)
                      select case(meanvartypes(j))
                      case('mean')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%mean3d, start=(/1,1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case('var')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%variance3d, start=(/1,1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case('min')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%min3d, start=(/1,1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case('max')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%max3d, start=(/1,1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case default
                         write(0,*) 'Can''t handle cell method: ', trim(meanvartypes(j)), ' of mnemonic', trim(t%name)
                      end select
                   case(4)
                      select case(meanvartypes(j))
                      case('mean')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%mean4d, start=(/1,1,1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case('var')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%variance4d, start=(/1,1,1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case('min')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%min4d, start=(/1,1,1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case('max')
                         status = nf90_put_var(ncid, meanvarids(i,j), meansparsglobal(i)%max4d, start=(/1,1,1,1,tpar%itm-1/) )
                         if (status /= nf90_noerr) call handle_err(status)
                      case default
                         write(0,*) 'Can''t handle cell method: ', trim(meanvartypes(j)), ' of mnemonic', trim(t%name)
                      end select
                   case default
                      write(0,*) 'Can''t handle rank: ', t%rank, ' of mnemonic', mnem
                   end select
                case default
                   write(0,*) 'Can''t handle type: ', t%type, ' of mnemonic', mnem
                end select
             end do
          end do
       end if
    end if


    if (xmaster) then
       status = nf90_close(ncid=ncid)
       if (status /= nf90_noerr) call handle_err(status) 
    end if
  end subroutine ncoutput

  character(256) function dimensionnames(dimids)
    implicit none
    integer, dimension(:), intent(in)           :: dimids ! store the dimids in a vector

    integer :: i, status
    character(80)  :: dimensionname 
    ! combine all the dimensionnames
    ! assumes all dimensions have an accompanying variable that should be used for coordinates.
    ! ",".join would have been nice here....
    dimensionnames = ''
    ! Fortran array dimensions are in reverse order
    do i=size(dimids),2,-1
       status = nf90_inquire_dimension(ncid, dimids(i), name=dimensionname)
       if (status /= nf90_noerr) call handle_err(status) 
       dimensionnames = trim(dimensionnames) // trim(dimensionname) // ',' 
    end do
    status = nf90_inquire_dimension(ncid, dimids(1), name=dimensionname)
    if (status /= nf90_noerr) call handle_err(status) 
    dimensionnames = trim(dimensionnames) // trim(dimensionname)
  end function dimensionnames

  integer function dimensionid(expression)
    ! Function to transform the expression in spaceparams.tmpl to an id, we might want this in the 
    ! makeincludes module
    use logging_module
    implicit none
    character(len=20),intent(in)                       :: expression
    
    select case(trim(expression))
    case('s%nx+1')
       dimensionid = xdimid
    case('s%ny+1')
       dimensionid = ydimid
    case('s%ntheta')
       dimensionid = thetadimid
    case('s%tidelen')
       dimensionid = tidetimedimid
    case('par%tideloc')
       dimensionid = tidecornersdimid
    case('s%windlen')
       dimensionid = windtimedimid
    case('par%ngd')
       dimensionid = sedimentclassesdimid
    case('s%ntdisch')
       dimensionid = inoutdimid
    case('2')
       dimensionid = inoutdimid
    case('max(par%nd,2)')
       dimensionid = bedlayersdimid
    case default
       call writelog('els','','Unknown dimension expression:'  // expression)
       stop 1
    end select
  end function dimensionid
#endif
end module ncoutput_module
