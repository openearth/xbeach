!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   MODULE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this module will provide for netcdf output in the way varoutput provides for binary output.
! It should be called optionally through one of the params.txt settings
! It will also be compiled conditionally like mpi. Only if it proves usefull will it be added by default.
! it will add dependencies on the netcdf fortran library (http://www.unidata.ucar.edu/software/netcdf/)
!
! With contributions from Uwe Rosebrock (CSIRO).
!
! if NCSINGLE is defined, output of all real variables in single precision
!                         else double precision
!
! P+R: tijdelijk uit voor testbed
!#define NCSINGLE
!#ifdef NCSINGLE
!#define NCREAL NF90_REAL
!#define CONVREAL sngl
!#define CONVREALTYPE real*4
!#else
#define NCREAL NF90_DOUBLE
#define CONVREAL
#define CONVREALTYPE real*8
!#endif

module ncoutput_module

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
   use xmpi_module
#ifdef USENETCDF
   use netcdf
#endif
   use typesandkinds
   use mnemmodule
   implicit none
   save
   private
   public ncoutput, fortoutput_init
#ifdef USENETCDF
   public ncoutput_init
#endif

   ! wwvv todo why the save's?
   !        see http://stackoverflow.com/questions/2893097/fortran-save-statement
   !        so, we should use save everywhere in modules, or use them in main program
   !        we could also add a line
   !        save
   !        after 'implicit none'

   integer :: ncid

   ! parameters
   integer :: parvarid

   ! grid
   integer :: xdimid, ydimid
   integer :: xvarid, yvarid

   ! Wave angle
   integer :: thetadimid
   ! Sediment
   integer :: sedimentclassesdimid, bedlayersdimid
   ! Drifters
   integer :: drifterdimid
   ! Ships
   integer :: shipdimid

   ! global
   integer, dimension(:), allocatable :: globalvarids
   ! default output (fixed length)

   ! points
   integer :: pointsdimid
   integer :: xpointsvarid, ypointsvarid, pointtypesvarid, xpointindexvarid, ypointindexvarid
   integer, dimension(:), allocatable :: pointsvarids
   integer, dimension(:),allocatable  :: xpoints     ! model x-coordinate of output points
   integer, dimension(:),allocatable  :: ypoints     ! model y-coordinate of output points

   ! mean
   ! number of variables by number of parameters per variable (mean, sigma^2, min, max)
   integer, dimension(:,:), allocatable       :: meanvarids
   character(slen), dimension(:), allocatable :: meanvartypes
   integer*4                           :: nmeanvartypes  = 4   ! number of time-average variable types
   integer*4,dimension(:),allocatable         :: rugrowindex ! Array with row index where runup gauge can be found



   ! time
   integer :: globaltimedimid, pointtimedimid, meantimedimid
   integer :: globaltimevarid, pointtimevarid, meantimevarid

   ! TODO: check out why these are sometimes used....
   integer :: tidetimedimid, windtimedimid
   integer :: inoutdimid, tidecornersdimid

   ! local variables
   integer :: npointstotal
   logical :: pointoutput
   integer :: itg,itp,itc,itm,itd
   ! Only alive at xmaster
   integer :: stpm        ! size of tpm
   integer :: wordsize    ! size of word in bytes

   integer                     :: noutnumbers = 0  ! the number of outnumbers
   integer, dimension(numvars) :: outnumbers  ! numbers, corrsponding to mnemonics, which are to be output

contains


#ifdef USENETCDF
   ! Error handling of netcdf errors
   subroutine handle_err(status,file,line)
      use netcdf

      integer, intent ( in) :: status
      character(*), intent(in) :: file
      integer, intent ( in)    :: line
      integer :: status2

      if(status /= nf90_noerr) then
         !UNIT=6 for stdout and UNIT=0 for stderr.
         write(0,'(a,i6,":",a)') file,line,trim(nf90_strerror(status))
         write(0,*) 'closing file'
         status2 = nf90_close(ncid)
         if (status2 /= nf90_noerr) then
            write(0,*) trim(nf90_strerror(status2))
         end if
         call halt_program
      end if
   end subroutine handle_err
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!   INITIALISE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USENETCDF
   subroutine ncoutput_init(s, sl, par, tpar)
      use xmpi_module
      use indextos_module
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

      type(spacepars), intent(inout)               :: s ! Use s => global data and sl => local data
      type(spacepars), intent(in)                  :: sl ! Use s => global data and sl => local data
      type(parameters), intent(in)                 :: par
      type(timepars), intent(in)                   :: tpar

      ! Part of the getkey
      type(parameter)                              :: val
      type(arraytype)                              :: t
      type(meanspars)                              :: meanvar
      integer                                      :: i,j
      integer                                      :: rc ! return code
      character(slen)                              :: mnem

      integer                                      :: npointstotal
      logical                                      :: outputp, outputg, outputm
      integer, dimension(:), allocatable           :: dimids ! store the dimids in a vector
      character(slen)                              :: coordinates
      character(slen)                              :: cellmethod

      character(slen), dimension(:), allocatable       :: keys
      logical :: dofortran, donetcdf
      ! subversion information
      include 'version.def'
      include 'version.dat'


      if (.not. xomaster) return

      dofortran = par%outputformat .eq. OUTPUTFORMAT_FORTRAN .or. &
      &par%outputformat .eq. OUTPUTFORMAT_DEBUG
      donetcdf = par%outputformat .eq. OUTPUTFORMAT_NETCDF .or. &
      &par%outputformat .eq. OUTPUTFORMAT_DEBUG


      outputp = .false.


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
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

      ! dimensions TODO: only output dimensions that are used
      ! grid
      status = nf90_def_dim(ncid, 'globalx', s%nx+1, xdimid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_def_dim(ncid, 'globaly', s%ny+1, ydimid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

      ! wave angles
      status = nf90_def_dim(ncid, 'wave_angle', s%ntheta, thetadimid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

      ! computational layers in bed ...
      ! TODO: Clean this up, why max(par%nd,2)???
      status = nf90_def_dim(ncid, 'bed_layers', max(par%nd,2), bedlayersdimid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

      ! sediment classes
      status = nf90_def_dim(ncid, 'sediment_classes', par%ngd, sedimentclassesdimid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

      ! dimensions of length 2.... what is this.... TODO: find out what this is
      status = nf90_def_dim(ncid, 'inout', 2, inoutdimid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

      write(*,*) 'Writing ndrifter', par%ndrifter
      if (par%ndrifter .gt. 0) then
         status = nf90_def_dim(ncid, 'ndrifter', par%ndrifter, drifterdimid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      end if

      !write(*,*) 'Writing nship', par%nship
      if (par%nship .gt. 0) then
         status = nf90_def_dim(ncid, 'nship', par%nship, shipdimid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      end if

      ! time dimensions are fixed, only defined if there are points
      if (outputg) then
         status = nf90_def_dim(ncid, 'globaltime', NF90_unlimited, globaltimedimid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      end if
      if (outputp) then
         ! points
         status = nf90_def_dim(ncid, 'points', npointstotal, pointsdimid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_def_dim(ncid, 'pointtime', size(tpar%tpp), pointtimedimid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      end if
      if (outputm) then
         status = nf90_def_dim(ncid, 'meantime', size(tpar%tpm)-1, meantimedimid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      end if

      if (s%tidelen > 0) then
         status = nf90_def_dim(ncid, 'tidetime', s%tidelen, tidetimedimid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      endif

      if (par%tideloc > 0) then
         status = nf90_def_dim(ncid, 'tidecorners', par%tideloc, tidecornersdimid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      endif

      if (s%windlen > 0) then
         status = nf90_def_dim(ncid, 'windtime', s%windlen, windtimedimid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      endif

      ! define empty parameter variable
      status = nf90_def_var(ncid, 'parameter', NCREAL, varid=parvarid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      ! define space & time variables
      ! grid
      status = nf90_def_var(ncid, 'globalx', NCREAL, (/ xdimid, ydimid /), xvarid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid, xvarid, 'units', 'm')
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid, xvarid, 'long_name', 'local x coordinate')
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid, xvarid, 'standard_name', 'projection_x_coordinate')
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid, xvarid, 'axis', 'X')
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      ! For compatibility with CSIRO Dive software
      if (len(trim(par%projection)) .ne. 0)  then
         status = nf90_put_att(ncid, xvarid, 'projection', par%projection)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, xvarid, 'rotation',( s%alfa/atan(1.0d0)*45.d0))
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      end if

      status = nf90_def_var(ncid, 'globaly', NCREAL, (/ xdimid, ydimid /), yvarid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid, yvarid, 'units', 'm')
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid, yvarid, 'long_name', 'local y coordinate')
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid, yvarid, 'standard_name', 'projection_y_coordinate')
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid, yvarid, 'axis', 'Y')
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      ! For compatibility with CSIRO Dive software
      if (len(trim(par%projection)) .ne. 0)  then
         status = nf90_put_att(ncid, yvarid, 'projection', par%projection)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, yvarid, 'rotation',( s%alfa/atan(1.0d0)*45.d0))
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      end if

      ! Some metadata attributes
      status = nf90_put_att(ncid,nf90_global, "Conventions", "CF-1.4")
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid,nf90_global, "Producer", "XBeach littoral zone wave model (http://www.xbeach.org)")
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid,nf90_global, "Build-Revision", trim(Build_Revision))
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid,nf90_global, "Build-Date", trim(Build_Date))
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      status = nf90_put_att(ncid,nf90_global, "URL", trim(Build_URL))
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

      ! Store all the parameters
      ! This part is awaiting comments from Robert McCall
      call getkeys(par, keys)
      do i=1,size(keys)
         rc = getkey(par, keys(i), val)
         if (val%type == 'i') then
            status = nf90_put_att(ncid, parvarid, keys(i), val%i0 )
         elseif (val%type == 'c') then
            status = nf90_put_att(ncid, parvarid, keys(i), val%c0 )
         elseif (val%type == 'r') then
            status = nf90_put_att(ncid, parvarid, keys(i), val%r0 )
         end if
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      end do

      ! global
      if (outputg) then
         status = nf90_def_var(ncid, 'globaltime', NCREAL, (/ globaltimedimid /), globaltimevarid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, globaltimevarid, 'units', trim(par%tunits))
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, globaltimevarid, 'axis', 'T')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, globaltimevarid, 'standard_name', 'time')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

         ! default global output variables
         do i=1,par%nglobalvar
            mnem = par%globalvars(i)
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
               &dimensionid(t%dimensions(3)), globaltimedimid /)
               coordinates = 'globalx globaly'
               ! Do we have a vertical level?
               if (dimids(3) .eq. bedlayersdimid) coordinates = trim(coordinates) // ' bed_layers'
             case(4)
               call writelog('ls', '', 'Variable ' // trim(mnem) // ' is of rank 4. This may not work due to an' // &
               &' unresolved issue. If so, remove the variable or use the fortran outputformat option.')
               dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), &
               &dimensionid(t%dimensions(3)), dimensionid(t%dimensions(4)), globaltimedimid /)
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
               &dimids, globalvarids(i))
               if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
             case('r')
               status = nf90_def_var(ncid, trim(mnem), NCREAL, &
               &dimids, globalvarids(i))
               if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
             case default
               write(0,*) 'mnem', mnem, ' not supported, type:', t%type
            end select
            status = nf90_put_att(ncid, globalvarids(i), 'coordinates', trim(coordinates))
            if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
            deallocate(dimids)
            status = nf90_put_att(ncid, globalvarids(i), 'units', trim(t%units))
            if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
            if (.not.(trim(t%standardname) .eq. '')) then
               status = nf90_put_att(ncid, globalvarids(i), 'standard_name', trim(t%standardname))
               if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
            endif
            status = nf90_put_att(ncid, globalvarids(i), 'long_name', trim(t%description))
            if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         end do
      end if

      !  ! points
      ! default global output variables
      if (outputp) then
         status = nf90_def_var(ncid, 'pointtime', NCREAL, (/ pointtimedimid /), pointtimevarid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, pointtimevarid, 'units', trim(par%tunits))
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, pointtimevarid, 'axis', 'T')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, pointtimevarid, 'standard_name', 'time')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

         ! points
         status = nf90_def_var(ncid, 'pointx', NCREAL, (/ pointsdimid /), xpointsvarid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, xpointsvarid, 'units', 'm')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, xpointsvarid, 'long_name', 'local x coordinate')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, xpointsvarid, 'standard_name', 'projection_x_coordinate')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, xpointsvarid, 'axis', 'X')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)


         status = nf90_def_var(ncid, 'pointy', NCREAL, (/ pointsdimid /), ypointsvarid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, ypointsvarid, 'units', 'm')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, ypointsvarid, 'long_name', 'local y coordinate')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, ypointsvarid, 'standard_name', 'projection_y_coordinate')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, ypointsvarid, 'axis', 'Y')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

         status = nf90_def_var(ncid, 'xpointindex', NF90_INT, (/ pointsdimid /), xpointindexvarid)
         ! wwvv above was NF90_DOUBLE
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, xpointindexvarid, 'long_name', 'nearest x grid cell')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_def_var(ncid, 'ypointindex', NF90_INT, (/ pointsdimid /), ypointindexvarid)
         ! wwvv above was NF90_DOUBLE
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, ypointindexvarid, 'long_name', 'nearest y grid cell')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

         status = nf90_def_var(ncid, 'pointtypes', NF90_INT, (/ pointsdimid /), pointtypesvarid)
         ! wwvv above was NF90_DOUBLE
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, pointtypesvarid, 'long_name', 'type of point (0=point, 1=rugauge)')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

         do i=1,par%npointvar
            mnem = par%pointvars(i)
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
                  &' unresolved issue. If so, remove the variable or use the fortran outputformat option.')
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
               status = nf90_def_var(ncid, 'point_' // trim(mnem), NCREAL, &
               &dimids, pointsvarids(i))
               if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
               status = nf90_put_att(ncid, pointsvarids(i), 'coordinates', trim(coordinates))
               if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
               deallocate(dimids)
             case default
               write(0,*) 'mnem', mnem, ' not supported, type:', t%type
            end select
            status = nf90_put_att(ncid, pointsvarids(i), 'units', trim(t%units))
            if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
            if (.not.(trim(t%standardname) .eq. '')) then
               status = nf90_put_att(ncid, pointsvarids(i), 'standard_name', trim(t%standardname))
               if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
            endif
            status = nf90_put_att(ncid, pointsvarids(i), 'long_name', trim(t%description))
            if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         end do
      endif ! outputp

      if (outputm) then
         status = nf90_def_var(ncid, 'meantime', NCREAL, (/ meantimedimid /), meantimevarid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, meantimevarid, 'units', trim(par%tunits))
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, meantimevarid, 'axis', 'T')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_att(ncid, meantimevarid, 'standard_name', 'time')
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
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
                  &dimensionid(t%dimensions(3)), meantimedimid /)
                  coordinates = 'globalx globaly'
                  ! Do we have a vertical level?
                  if (dimids(3) .eq. bedlayersdimid) coordinates = trim(coordinates) // ' bed_layers'
                case(4)
                  call writelog('ls', '', 'Variable ' // trim(mnem) // ' is of rank 4. This may not work due to an' // &
                  &' unresolved issue. If so, remove the variable or use the fortran outputformat option.')
                  dimids = (/ dimensionid(t%dimensions(1)), dimensionid(t%dimensions(2)), &
                  &dimensionid(t%dimensions(3)), dimensionid(t%dimensions(4)), meantimedimid /)
                  coordinates = 'globalx globaly'
                  ! Do we have a vertical level?
                  if ((dimids(3) .eq. bedlayersdimid) .or. (dimids(4) .eq. bedlayersdimid)) then
                     coordinates = trim(coordinates) // ' bed_layers'
                  end if
                case default
                  call writelog('lse', '', 'mnem: ' // mnem // ' not supported, rank:', t%rank)
                  call halt_program
                  stop 1
               end select

               ! Create a variable for all types of meanvars (mean, var, min, max)
               do j = 1,nmeanvartypes
                  cellmethod = meanvartypes(j)
                  call writelog('ls', '', 'Creating netcdf variable: ',  trim(t%name) // '_' // cellmethod)
                  status = nf90_def_var(ncid, trim(t%name) // '_' // trim(cellmethod), NCREAL, &
                  &dimids, meanvarids(i,j))
                  if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  status = nf90_put_att(ncid, meanvarids(i,j), 'coordinates', trim(coordinates))
                  if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  if (cellmethod .eq. 'var') then
                     status = nf90_put_att(ncid, meanvarids(i,j), 'units', '(' // trim(t%units) // ')^2')
                  else
                     status = nf90_put_att(ncid, meanvarids(i,j), 'units', trim(t%units))
                  endif
                  if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  if (.not.(trim(t%standardname) .eq. '')) then
                     status = nf90_put_att(ncid, meanvarids(i,j), 'standard_name', trim(t%standardname))
                     if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  endif
                  status = nf90_put_att(ncid, meanvarids(i,j), 'long_name', trim(t%description))
                  if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  ! For H and urms we don't compute the mean but the rms of the rms.....
                  if (cellmethod .eq. 'mean' .and. ((t%name .eq. 'H') .or. (t%name .eq. 'urms')))  then
                     cellmethod = 'rms'
                  end if
                  status = nf90_put_att(ncid, meanvarids(i,j), 'cell_methods', 'meantime: ' // trim(cellmethod))
                  if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
               end do

               deallocate(dimids)
             case default
               write(0,*) 'mnem', mnem, ' not supported, type:', t%type
            end select
         end do
      endif  ! outputm

      ! done defining variables
      call writelog('ls', '', 'Writing file definition.')
      status = nf90_enddef(ncid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

      ! Fill meta variables
      ! Grid
      j = chartoindex('xz')
      call indextos(s,j,t)

      status = nf90_put_var(ncid, xvarid, CONVREAL(t%r2))
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

      j = chartoindex('yz')
      call indextos(s,j,t)
      status = nf90_put_var(ncid, yvarid, CONVREAL(t%r2))
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

      if (outputp) then
         call writelog('ls', '', 'Writing point vars.')
         status = nf90_put_var(ncid, xpointsvarid, CONVREAL(par%xpointsw))
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_var(ncid, ypointsvarid, CONVREAL(par%ypointsw))
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_var(ncid, pointtypesvarid, par%pointtypes)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)

         ! Convert world coordinates of points to nearest (lsm) grid point
         ! This could be done in some postprocessing function

         if(.not. allocated(xpoints)) then
            allocate(xpoints(npointstotal))
            allocate(ypoints(npointstotal))
            call snappointstogrid(par, s, xpoints, ypoints)
         endif

         status = nf90_put_var(ncid, xpointindexvarid, xpoints)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_var(ncid, ypointindexvarid, ypoints)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         status = nf90_put_var(ncid, pointtypesvarid, par%pointtypes)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      end if

      status = nf90_close(ncid)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      ! wwvv sl is not used so it should be removed.
      ! wwvv for now, to avoid warning:
      if (sl%nx .ne. -1) return
   end subroutine ncoutput_init
#endif
   ! USENETCDF


   !___________________________________________________________________________________

   subroutine ncoutput(s,sl,par, tpar)
      use logging_module
      use indextos_module
#ifdef USEMPI
      use xmpi_module
#endif
      use params
      use spaceparams
      use timestep_module
      use mnemmodule
      use means_module
      use postprocessmod

      implicit none

      type(spacepars), intent(inout)         :: s ! s-> spaceparams, what is isl?
      !                                         s describes the global system, sl the local system
      !                                           in this MPI process
      type(spacepars), intent(inout)         :: sl ! s-> spaceparams, what is isl?
      type(parameters), intent(inout)        :: par
      type(timepars), intent(in)             :: tpar

      type(arraytype)                        :: t
      integer                                :: i,j,ii
#ifdef USEMPI
      integer                                :: index
#endif
      character(slen)                        :: mnem
      real*8, dimension(:,:), allocatable    :: points

      ! some local variables to pass the data through the postprocessing function.
      integer :: i0
      integer, dimension(:,:), allocatable :: i2
      integer, dimension(:,:,:), allocatable :: i3
      CONVREALTYPE                                  :: r0conv
      real*8, dimension(:), allocatable :: r1
      CONVREALTYPE, dimension(:), allocatable       :: r1conv
      real*8, dimension(:,:), allocatable :: r2
      CONVREALTYPE, dimension(:,:), allocatable     :: r2conv
      real*8, dimension(:,:,:), allocatable :: r3
      CONVREALTYPE, dimension(:,:,:), allocatable   :: r3conv
      real*8, dimension(:,:,:,:), allocatable :: r4
      CONVREALTYPE, dimension(:,:,:,:), allocatable :: r4conv
      real*8, allocatable                           :: tempvectorr(:)
      real*8,dimension(size(tpar%tpg)+size(tpar%tpp)+size(tpar%tpc)+size(tpar%tpm)) :: outputtimes

      integer                                       :: jtg,reclen,unit,idumhl,ird,xmax
      integer                                       :: iz,jz
      real*8                                        :: di,dj,dx,dy
#ifdef USEMPI
      real*8, dimension(par%ndrifter)               :: idriftlocal,jdriftlocal
#endif


#ifdef USENETCDF
      integer :: status
#endif

      logical :: dofortran, donetcdf, dofortran_compat
      logical :: dooutput_global, dooutput_mean, dooutput_point

      ! fortran output requested?
      dofortran = par%outputformat .eq. OUTPUTFORMAT_FORTRAN .or. &
      &par%outputformat .eq. OUTPUTFORMAT_DEBUG

      ! for compatibility with older version, at some places
      ! fortran output is only done when tpar%output = .true.

      dofortran_compat = dofortran .and. tpar%output

      !netcdf output requested?
      donetcdf = par%outputformat .eq. OUTPUTFORMAT_NETCDF .or. &
      &par%outputformat .eq. OUTPUTFORMAT_DEBUG

      ! time for global output?
      dooutput_global = tpar%outputg .and. par%nglobalvar .gt. 0

      ! time for mean output?
      dooutput_mean   = tpar%outputm .and. par%nmeanvar   .gt. 0 .and. tpar%itm .gt. 1

      ! time for point output?
      dooutput_point  = tpar%outputp .and. par%npointvar .gt. 0


#ifdef USEMPI
      ! clear collected items
      s%collected = s%precollected

      ! If we're gonna write some global output
      if (dooutput_global) then
         ! we'll need to collect the information from all nodes.
         do i=1,par%nglobalvar
            mnem = par%globalvars(i)
            index = chartoindex(mnem)
            call space_collect_index(s,sl,par,index)
            !  we have to make sure that the extra information needed
            !  is also collected
            !
            ! get extra needed s%vars
            !
            call indextos(s,index,t)
            !
            ! the following calls to gridrotate will ensure that
            ! data, required for rotation of arrays mentioned in
            ! par%globalvars will be collected.
            !  Note that this special call of gridrotate only collects
            !  missing data, no rotation is done.
            !
            select case(t%type)
             case('r')
               select case(t%rank)
                case(2)
                  allocate(r2(0,0))
                  call gridrotate(par, s, t, r2, sl)
                  deallocate(r2)
                case(3)
                  allocate(r3(0,0,0))
                  call gridrotate(par, s, t, r3, sl)
                  deallocate(r3)
               end select
            end select
         end do
      endif

      ! The same for mean output
      if (dooutput_mean) then
         do i=1,par%nmeanvar
            mnem = par%meanvars(i)
            index = chartoindex(mnem)
            call space_collect_index(s,sl,par,index)
            call indextos(s,index,t)
            select case(t%type)
             case('r')
               select case(t%rank)
                case(2)
                  allocate(r2(0,0))
                  call gridrotate(par, s, t, r2, sl)
                  deallocate(r2)
                case(3)
                  allocate(r3(0,0,0))
                  call gridrotate(par, s, t, r3, sl)
                  deallocate(r3)
               end select
            end select
         end do
      endif
#endif


#ifdef USEMPI
      if(dooutput_point) then
         ! Collect all data for which we store the points.
         ! TODO: This will be a lot faster if nodes write their own point. Use the parallel netcdf for that.
         ! Let's wait till someone needs it...
         ! wwvv Collecting the data could als be done using sends from the compute processes
         ! wwvv and receives on the xomaster process.
         ! wwvv the info to do this properly could be computed in _init
         ! wwvv in stead of sends and receives, probably more simple is the use of mpi_alltoallw
         !
         do i=1,par%npointvar
            mnem = par%pointvars(i)
            call space_collect_mnem(s,sl,par,mnem)
         end do
         ! wwvv temporary method to determine runup values
         ! wwvv we need hh and zs:
         if (par%nrugauge .ge. 0) then
            call space_collect_mnem(s,sl,par,mnem_hh)
            call space_collect_mnem(s,sl,par,mnem_zs)
            ! xz and yz should already be available, but just to be sure:
            call space_collect_mnem(s,sl,par,mnem_xz)
            call space_collect_mnem(s,sl,par,mnem_yz)
            ! if already present, they will not be recollected
         endif
      endif
#endif

      ! If we're gonna write some mean output
      if(dooutput_mean) then
         ! only write the information on the xomaster node

         do i=1,par%nmeanvar
#ifdef USEMPI
            call means_collect(sl,meansparsglobal(i),meansparslocal(i))
#else
            meansparsglobal(i)=meansparslocal(i)
#endif
         end do
      endif


#ifdef USEMPI
      if (par%ndrifter .gt. 0) then
         if (xmaster) then
            idriftlocal = sl%idrift
            jdriftlocal = sl%jdrift
         endif
         call xmpi_send(xmpi_imaster,xmpi_omaster,idriftlocal)
         call xmpi_send(xmpi_imaster,xmpi_omaster,jdriftlocal)
         if (xomaster) then
            s%idrift = idriftlocal
            s%jdrift = jdriftlocal
         endif
      endif
#endif


      ! writing is done by xomaster, the others processes go back to work

      if( .not. xomaster) return

      ! Open the output file

#ifdef USENETCDF
      if(donetcdf) then
         status = nf90_open(ncid=ncid, path=par%ncfilename, mode=nf90_write)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      endif
#endif

      !
      ! some variables can only be output if others are available, see the code
      ! for gridrotate.
      !
      if (dooutput_global) then
         itg = itg+1
         ! Store the time (in morphological time)
#ifdef USENETCDF
         if(donetcdf) then
            status = nf90_put_var(ncid, globaltimevarid, CONVREAL(par%t*max(par%morfac,1.d0)), (/tpar%itg/))
            if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         endif
#endif
         ! write global output variables
         do i=1,par%nglobalvar
            mnem = par%globalvars(i)
            j = chartoindex(mnem)
            ! lookup the proper array (should have been collected already)
            call indextos(s,j,t)

            select case(t%type)
             case('i')
               select case(t%rank)
                case(0)
                  ! no need to allocate here
                  call gridrotate(t,i0)
#ifdef USENETCDF
                  if(donetcdf) then
                     status = nf90_put_var(ncid, globalvarids(i), i0, start=(/1,tpar%itg/) )
                     if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  endif
#endif
                  if(dofortran) then
                     inquire(iolength=reclen) i0
                     call checkfile(i,unit,reclen,jtg)
                     write(unit,rec=jtg) i0
                     call flush(unit)
                  endif
                case(2)
                  allocate(i2(size(t%i2,1),size(t%i2,2)))
                  call gridrotate(t, i2)
#ifdef USENETCDF
                  if(donetcdf) then
                     status = nf90_put_var(ncid, globalvarids(i), i2, start=(/1,1,tpar%itg/) )
                     if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  endif
#endif
                  if(dofortran) then
                     inquire(iolength=reclen) i2
                     call checkfile(i,unit,reclen,jtg)
                     write(unit,rec=jtg) i2
                     call flush(unit)
                  endif
                  deallocate(i2)
                case(3)
                  allocate(i3(size(t%i3,1),size(t%i3,2),size(t%i3,3)))
                  call gridrotate(t, i3)
#ifdef USENETCDF
                  if(donetcdf) then
                     status = nf90_put_var(ncid, globalvarids(i), i3, start=(/1,1,1,tpar%itg/) )
                     if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  endif
#endif
                  if(dofortran) then
                     inquire(iolength=reclen) i3
                     call checkfile(i,unit,reclen,jtg)
                     write(unit,rec=jtg) i3
                     call flush(unit)
                  endif
                  deallocate(i3)
                case default
                  write(0,*) 'Can''t handle rank: ', t%rank, ' of mnemonic', mnem
               end select
             case('r')
               select case(t%rank)
                case(0)
                  ! no need to allocate here
                  r0conv = CONVREAL(t%r0)
#ifdef USENETCDF
                  if(donetcdf) then
                     status = nf90_put_var(ncid, globalvarids(i), r0conv, start=(/1,tpar%itg/) )
                     if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  endif
#endif
                  if(dofortran) then
                     inquire(iolength=reclen) r0conv
                     call checkfile(i,unit,reclen,jtg)
                     write(unit,rec=jtg) r0conv
                     call flush(unit)
                  endif
                case(1)
                  allocate(r1(size(t%r1,1)))
                  allocate(r1conv(size(t%r1,1)))
                  ! no need to rotate here
                  r1conv = CONVREAL(t%r1)
#ifdef USENETCDF
                  if(donetcdf) then
                     status = nf90_put_var(ncid, globalvarids(i), r1conv, start=(/1,tpar%itg/) )
                     if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  endif
#endif
                  if(dofortran) then
                     inquire(iolength=reclen) r1conv
                     call checkfile(i,unit,reclen,jtg)
                     write(unit,rec=jtg) r1conv
                     call flush(unit)
                  endif
                  deallocate(r1)
                  deallocate(r1conv)
                case(2)
                  allocate(r2(size(t%r2,1),size(t%r2,2)))
                  allocate(r2conv(size(t%r2,1),size(t%r2,2)))
                  call gridrotate(par, s, t, r2)
                  r2conv = CONVREAL(r2)
#ifdef USENETCDF
                  if(donetcdf) then
                     status = nf90_put_var(ncid, globalvarids(i), r2conv, start=(/1,1,tpar%itg/) )
                     if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  endif
#endif
                  if(dofortran) then
                     inquire(iolength=reclen) r2conv
                     call checkfile(i,unit,reclen,jtg)
                     write(unit,rec=jtg) r2conv
                     call flush(unit)
                  endif
                  deallocate(r2)
                  deallocate(r2conv)
                case(3)
                  allocate(r3(size(t%r3,1),size(t%r3,2),size(t%r3,3)))
                  allocate(r3conv(size(t%r3,1),size(t%r3,2),size(t%r3,3)))
                  call gridrotate(par, s, t, r3)
                  r3conv = CONVREAL(r3)
#ifdef USENETCDF
                  if(donetcdf) then
                     status = nf90_put_var(ncid, globalvarids(i), r3conv, start=(/1,1,1, tpar%itg/) )
                     if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  endif
#endif
                  if(dofortran) then
                     inquire(iolength=reclen) r3conv
                     call checkfile(i,unit,reclen,jtg)
                     write(unit,rec=jtg) r3conv
                     call flush(unit)
                  endif
                  deallocate(r3)
                  deallocate(r3conv)
                case(4)
                  allocate(r4(size(t%r4,1),size(t%r4,2),size(t%r4,3),size(t%r4,4)))
                  allocate(r4conv(size(t%r4,1),size(t%r4,2),size(t%r4,3),size(t%r4,4)))
                  call gridrotate(t, r4)
                  r4conv = CONVREAL(r4)
#ifdef USENETCDF
                  if(donetcdf) then
                     status = nf90_put_var(ncid, globalvarids(i), r4conv, start=(/1,1,1,1, tpar%itg/) )
                     if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                  endif
#endif
                  if(dofortran) then
                     inquire(iolength=reclen) r4conv
                     call checkfile(i,unit,reclen,jtg)
                     write(unit,rec=jtg) r4conv
                     call flush(unit)
                  endif
                  deallocate(r4)
                  deallocate(r4conv)
                case default
                  write(0,*) 'Can''t handle rank: ', t%rank, ' of mnemonic', mnem
               end select
             case default
               write(0,*) 'Can''t handle type: ', t%type, ' of mnemonic', mnem
            end select
         end do
      end if


      if(dooutput_point) then
         itp=itp+1

#ifdef USENETCDF
         if(donetcdf) then
            status = nf90_put_var(ncid, pointtimevarid, CONVREAL(par%t*max(par%morfac,1.d0)), (/tpar%itp/))
            if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         endif
#endif

         if(dofortran_compat) then
            if (par%npoints .gt. 0) then
               allocate(points(par%npoints,par%npointvar+1))
            else
               allocate(points(0,0))
            endif
         else
            allocate(points(0,0))
         endif

         do i=1,par%npointvar
            mnem = par%pointvars(i)
            j = chartoindex(mnem)
            ! lookup the proper array
            call indextos(s,j,t)
            ! get the proper output points ....
            ! I have no idea what is happening in varouput so I'll try it in a different way
            !TODO This is not very efficient because we are using the outer counters, reorder dimensions....
            select case(t%type)
             case('r')
               do ii = 1, par%npoints + par%nrugauge
                  select case(t%rank)
                   case(2)
                     ! This postprocessing creates an ugly dependency.
                     ! it would be nice if we could call gridrotate as a function
                     ! or if we could just have the postprocessing insert some reference processing routines to call
                     ! or if we could split this out of the case statement (dry)
                     ! or if we could defer this to a postprocessing routine (for example ncks)
                     allocate(r2(size(t%r2,1),size(t%r2,2)))
                     call gridrotate(par, s, t, r2)
#ifdef USENETCDF
                     if (donetcdf) then
                        status = nf90_put_var(ncid, pointsvarids(i), CONVREAL(r2(xpoints(ii), ypoints(ii))), start=(/ii,tpar%itp/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                     endif
#endif
                     if(dofortran_compat) then
                        if (ii .le. par%npoints) then
                           points(ii,i+1) = r2(xpoints(ii), ypoints(ii))
                        endif
                     endif
                     deallocate(r2)
                   case(3)
                     allocate(r3(size(t%r3,1),size(t%r3,2),size(t%r3,3)))
                     call gridrotate(par, s, t, r3)
#ifdef USENETCDF
                     if (donetcdf) then
                        status = nf90_put_var(ncid, pointsvarids(i), CONVREAL(r3(xpoints(ii), ypoints(ii),:)), start=(/ii,1,tpar%itp/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                     endif
#endif
                     if(dofortran_compat) then
                        if (ii .le. par%npoints) then
                           points(ii,i+1) = r3(xpoints(ii), ypoints(ii),1)    ! wwvv todo
                        endif
                     endif
                     deallocate(r3)
                   case(4)
                     allocate(r4(size(t%r4,1),size(t%r4,2),size(t%r4,3),size(t%r4,4)))
                     call gridrotate(t, r4)
#ifdef USENETCDF
                     if (donetcdf) then
                        status = nf90_put_var(ncid, pointsvarids(i), CONVREAL(r4(xpoints(ii), ypoints(ii),:,:)), &
                        &start=(/ii,1,1,tpar%itp/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                     endif
#endif
                     if(dofortran_compat) then
                        if (ii .le. par%npoints) then
                           points(ii,i+1) = r4(xpoints(ii), ypoints(ii),1,1)    ! wwvv todo
                        endif
                     endif
                     deallocate(r4)
                   case default
                     write(0,*) 'Can''t handle rank: ', t%rank, ' of mnemonic', mnem
                  end select
               end do
             case default
               write(0,*) 'Can''t handle type: ', t%type, ' of mnemonic', mnem
            end select
         enddo ! i=1,par%npointvar

         if(dofortran_compat) then
            do ii = 1,par%npoints
               if (par%morfacopt==1) then
                  points(ii,1)=par%t*max(par%morfac,1.d0)
               else
                  points(ii,1)=par%t
               endif
            end do
            do ii = 1,par%npoints
               write(indextopointsunit(ii),rec=tpar%itp)CONVREAL(points(ii,:))
               call flush(indextopointsunit(ii))
            enddo
            deallocate(points)

            if (par%nrugauge .gt. 0) then
               ! Set up runup gauge output vector
               allocate(tempvectorr(1+par%nrugdepth*3))
               tempvectorr=huge(0.d0)
               do i=1,par%nrugauge
                  do ird=1,par%nrugdepth
                     xmax = s%nx+1
                     idumhl = xmax        ! Set default
                     if (rugrowindex(i)>0) then  ! master domain always contains this runup gauge
                        ! local index of minimum location where hh<rugdepth
                        do ii=2,xmax
                           if ((s%hh(ii,rugrowindex(i))<=par%rugdepth(ird)) .and. &
                           &(s%hh(ii-1,rugrowindex(i))>par%rugdepth(ird)) ) then
                              idumhl=ii-1
                              exit
                           end if
                        enddo
                     end if
                     if (par%morfacopt==1) then
                        tempvectorr(1)=par%t*max(par%morfac,1.d0)
                     else
                        tempvectorr(1)=par%t
                     endif
                     tempvectorr((ird-1)*3+2)=s%xz(idumhl,rugrowindex(i))
                     tempvectorr((ird-1)*3+3)=s%yz(idumhl,rugrowindex(i))
                     tempvectorr((ird-1)*3+4)=s%zs(idumhl,rugrowindex(i))
                  enddo
                  write(indextopointsunit(i+par%npoints),rec=tpar%itp)CONVREAL(tempvectorr)
                  call flush(indextopointsunit(i+par%npoints))
               enddo  ! i=1,par%nrugauge
            endif !par%nrugauge
         endif ! dofortran
      endif   ! dooutput_point


#ifdef USENETCDF
      if(dooutput_mean .and. donetcdf) then
         ! Store the time (in morphological time)
         status = nf90_put_var(ncid, meantimevarid, CONVREAL(par%t*max(par%morfac,1.d0)), (/tpar%itm-1/))
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
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
                           status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(sqrt(meansparsglobal(i)%variancesquareterm2d)), &
                           &start=(/1,1,tpar%itm-1/) )
                        elseif (t%name .eq. 'thetamean') then
                           status = nf90_put_var(ncid, meanvarids(i,j), &
                           &CONVREAL( &
                           &mod(2.d0*par%px &
                           &+ atan2(nint(meansparsglobal(i)%mean2d)/1d7, &
                           &mod(meansparsglobal(i)%mean2d,1.d0)*1d1), 2.d0*par%px) / par%px * 180), &
                           &start=(/1,1,tpar%itm-1/) )
                        else
                           status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%mean2d), start=(/1,1,tpar%itm-1/) )
                        end if
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case('var')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%variance2d), start=(/1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case('min')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%min2d), start=(/1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case('max')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%max2d), start=(/1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case default
                        write(0,*) 'Can''t handle cell method: ', trim(meanvartypes(j)), ' of mnemonic', trim(t%name)
                     end select
                   case(3)
                     select case(meanvartypes(j))
                      case('mean')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%mean3d), start=(/1,1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case('var')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%variance3d), start=(/1,1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case('min')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%min3d), start=(/1,1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case('max')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%max3d), start=(/1,1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case default
                        write(0,*) 'Can''t handle cell method: ', trim(meanvartypes(j)), ' of mnemonic', trim(t%name)
                     end select
                   case(4)
                     select case(meanvartypes(j))
                      case('mean')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%mean4d), start=(/1,1,1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case('var')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%variance4d), start=(/1,1,1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case('min')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%min4d), start=(/1,1,1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                      case('max')
                        status = nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%max4d), start=(/1,1,1,1,tpar%itm-1/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
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
      endif ! dooutput_mean .and. donetcdf
#endif

      !if (par%nmeanvar>0 .and. dofortran .and. tpar%output) then
      ! Not at the first in tpm as this is the start of averaging. Only output after second in tpm
      if(dooutput_mean .and. dofortran_compat) then
         itm=itm+1  ! Note, this is a local counter, used to position in output file
         do i=1,par%nmeanvar
            select case (meansparsglobal(i)%rank)
             case (2)
               if (par%meanvars(i)=='H') then                ! Hrms changed to H
                  write(indextomeanunit(i),rec=itm)CONVREAL(sqrt(meansparsglobal(i)%variancesquareterm2d))
               elseif (par%meanvars(i)=='urms') then       ! urms
                  write(indextomeanunit(i),rec=itm)CONVREAL(sqrt(meansparsglobal(i)%variancesquareterm2d))
               elseif (par%meanvars(i)=='thetamean') then       ! thetamean
                  write(indextomeanunit(i),rec=itm) &
                  &CONVREAL( &
                  &mod(2.d0*par%px + atan2(nint(meansparsglobal(i)%mean2d)/1d7, &
                  &mod(meansparsglobal(i)%mean2d,1.d0)*1d1), 2.d0*par%px) / par%px * 180 &
                  &)
               else                                                    ! non-rms variables
                  write(indextomeanunit(i),rec=itm)CONVREAL(meansparsglobal(i)%mean2d)
               endif
               write(indextovarunit(i),rec=itm)CONVREAL(meansparsglobal(i)%variance2d)
               where(meansparsglobal(i)%min2d>0.99d0*huge(0.d0))
                  meansparsglobal(i)%min2d=-999.d0
               endwhere
               where(meansparsglobal(i)%max2d<-0.99d0*huge(0.d0))
                  meansparsglobal(i)%max2d=-999.d0
               endwhere
               write(indextominunit(i),rec=itm)CONVREAL(meansparsglobal(i)%min2d)
               write(indextomaxunit(i),rec=itm)CONVREAL(meansparsglobal(i)%max2d)
             case (3)
               write(indextomeanunit(i),rec=itm)CONVREAL(meansparsglobal(i)%mean3d)
               write(indextovarunit(i),rec=itm)CONVREAL(meansparsglobal(i)%variance3d)
               where(meansparsglobal(i)%min3d>0.99d0*huge(0.d0))
                  meansparsglobal(i)%min3d=-999.d0
               endwhere
               where(meansparsglobal(i)%max3d<-0.99d0*huge(0.d0))
                  meansparsglobal(i)%max3d=-999.d0
               endwhere
               write(indextominunit(i),rec=itm)CONVREAL(meansparsglobal(i)%min3d)
               write(indextomaxunit(i),rec=itm)CONVREAL(meansparsglobal(i)%max3d)
             case (4)
               write(indextomeanunit(i),rec=itm)CONVREAL(meansparsglobal(i)%mean4d)
               write(indextovarunit(i),rec=itm)CONVREAL(meansparsglobal(i)%variance4d)
               where(meansparsglobal(i)%min4d>0.99d0*huge(0.d0))
                  meansparsglobal(i)%min4d=-999.d0
               endwhere
               where(meansparsglobal(i)%max4d<-0.99d0*huge(0.d0))
                  meansparsglobal(i)%max4d=-999.d0
               endwhere
               write(indextominunit(i),rec=itm)CONVREAL(meansparsglobal(i)%min4d)
               write(indextomaxunit(i),rec=itm)CONVREAL(meansparsglobal(i)%max4d)
            end select
            call flush(indextomeanunit(i))
            call flush(indextovarunit(i))
            call flush(indextominunit(i))
            call flush(indextomaxunit(i))
         enddo
         par%tintm=tpar%tpm(min(itm+2,stpm))-tpar%tpm(itm+1)  ! Next averaging period (min to stop array out of bounds)
         par%tintm=max(par%tintm,tiny(0.d0))        ! to prevent par%tintm=0 after last output
      endif  ! dooutput_mean .and. dofortran_compat

#ifdef USENETCDF
      if(donetcdf) then
         status = nf90_close(ncid=ncid)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      endif
#endif

      if(dofortran_compat) then
         if (abs(mod(par%t,par%tintp))<1.d-6) then
            itd = itd+1
            do i=1,par%ndrifter
               if (  par%t>=s%tdriftb(i) .and. par%t<=s%tdrifte(i) .and. &
               &s%idrift(i)>1       .and. s%idrift(i)<=s%nx   .and. &
               &s%jdrift(i)>1       .and. s%jdrift(i)<=s%ny             ) then

                  iz = int(s%idrift(i))
                  jz = int(s%jdrift(i))

                  di = mod(s%idrift(i),1.d0)
                  dj = mod(s%jdrift(i),1.d0)

                  dx = di*s%dsu(iz,jz)*cos(s%alfaz(iz,jz)) - &
                  &dj*s%dnv(iz,jz)*sin(s%alfaz(iz,jz))
                  dy = di*s%dsu(iz,jz)*sin(s%alfaz(iz,jz)) + &
                  &dj*s%dnv(iz,jz)*cos(s%alfaz(iz,jz))

                  write(indextodrifterunit(i),rec=itd)    &
                  &CONVREAL(s%xz(iz,jz)+dx),                     &
                  &CONVREAL(s%yz(iz,jz)+dy),                     &
                  &CONVREAL(par%t)
               else
                  write(indextodrifterunit(i),rec=itd)    &
                  &CONVREAL(-999d0),                    &
                  &CONVREAL(-999d0),                    &
                  &CONVREAL(par%t)
               endif
               call flush(indextodrifterunit(i))
            enddo
         endif
      endif ! dofortran_compat

      if(dofortran_compat) then
         outputtimes=-999.d0
         outputtimes(1:itg)=tpar%tpg(1:itg)
         outputtimes(itg+1:itg+itp)=tpar%tpp(1:itp)
         outputtimes(itg+itp+1:itg+itp+itc)=tpar%tpc(1:itc)
         outputtimes(itg+itp+itc+1:itg+itp+itc+itm)=tpar%tpm(2:itm+1)          ! mean output always shifted by 1
         if (par%morfacopt==1) outputtimes=outputtimes*max(par%morfac,1.d0)
         open(1998,file='dims.dat',form='unformatted',access='direct',recl=wordsize*(10+size(outputtimes)))
         write(1998,rec=1) CONVREAL(itg*1.d0),&
         &CONVREAL(s%nx*1.d0),&
         &CONVREAL(s%ny*1.d0),&
         &CONVREAL(s%ntheta*1.d0),&
         &CONVREAL(par%kmax*1.d0),&
         &CONVREAL(par%ngd*1.d0),&
         &CONVREAL(par%nd*1.d0), &
         &CONVREAL(itp*1.d0),&
         &CONVREAL(itc*1.d0),&
         &CONVREAL(itm*1.d0),&
         &CONVREAL(outputtimes)
         call flush(1998)
      end if
      ! wwvv avoid warning about unused sl:
      if (sl%nx .eq. -1) return
   end subroutine ncoutput

#ifdef USENETCDF
   character(slen) function dimensionnames(dimids)
      implicit none
      integer, dimension(:), intent(in)           :: dimids ! store the dimids in a vector

      integer :: i, status
      character(slen)  :: dimensionname
      ! combine all the dimensionnames
      ! assumes all dimensions have an accompanying variable that should be used for coordinates.
      ! ",".join would have been nice here....
      dimensionnames = ''
      ! Fortran array dimensions are in reverse order
      do i=size(dimids),2,-1
         status = nf90_inquire_dimension(ncid, dimids(i), name=dimensionname)
         if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
         dimensionnames = trim(dimensionnames) // trim(dimensionname) // ','
      end do
      status = nf90_inquire_dimension(ncid, dimids(1), name=dimensionname)
      if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
      dimensionnames = trim(dimensionnames) // trim(dimensionname)
   end function dimensionnames

   integer function dimensionid(expression)
      ! Function to transform the expression in spaceparams.tmpl to an id, we might want this in the
      ! makeincludes module
      use logging_module
      implicit none
      character(len=*),intent(in) :: expression
      select case(expression)
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
       case('par%ndrifter')
         dimensionid = drifterdimid
       case('par%nship')
         dimensionid = shipdimid
       case default
         call writelog('els','','Unknown dimension expression:'  // expression)
         stop 1
      end select
   end function dimensionid
#endif
   ! USENETCDF

   subroutine fortoutput_init(s,par,tpar)
      use params
      use spaceparams
      use readkey_module
      use indextos_module
      use timestep_module
      use logging_module
      use postprocessmod
      use filefunctions
      use mnemmodule

      implicit none

      type(spacepars),intent(in)          :: s
      type(parameters),intent(in)         :: par
      type(timepars),intent(in)           :: tpar

      integer                             :: i,fid
      integer                             :: reclen,reclenm,reclenp
      character(100)                      :: fname,fnamemean,fnamevar,fnamemin,fnamemax
      type(arraytype)                     :: t

      if (.not. xomaster) return

      reclenm = -123

      ! Initialize places in output files
      itg = 0
      itm = 0
      itp = 0
      itc = 0
      itd = 0
      stpm = size(tpar%tpm)

      ! Record size for global and mean output
      inquire(iolength=wordsize) CONVREAL(1.0d0)
      reclen=wordsize*(s%nx+1)*(s%ny+1)

      open(100,file='xy.dat',form='unformatted',access='direct',recl=reclen,status='REPLACE')
      write(100,rec=1)CONVREAL(s%xz)
      write(100,rec=2)CONVREAL(s%yz)
      write(100,rec=3)CONVREAL(s%x)
      write(100,rec=4)CONVREAL(s%y)
      close(100)

      !     GLOBAL VARS

      noutnumbers = par%nglobalvar
      ! store all indices for the global variables
      do i= 1,noutnumbers
         outnumbers(i) = chartoindex(par%globalvars(i))
      enddo

      do i=1,par%npoints+par%nrugauge
         if (par%pointtypes(i)==0) then
            write(fname,'("point",i0.3,".dat")') i
         else
            write(fname,'("rugau",i0.3,".dat")') i-par%npoints
         endif
         if (par%pointtypes(i)==0) then
            reclenp=wordsize*(par%npointvar+1)
         else
            reclenp=wordsize*(1+par%nrugdepth*3)
         endif
         open(indextopointsunit(i),file=fname,&
         &form='unformatted',access='direct',recl=reclenp,status='REPLACE')
      enddo
      if (par%npoints>0) then
         ! write index file of point output variables
         fid=create_new_fid()
         open(fid,file='pointvars.idx',status='replace',action='write')
         do i=1,par%npointvar
            write(fid,*)trim(par%pointvars(i))
         enddo
         close(fid)
      endif

      if (par%npoints + par%nrugauge > 0) then
         if(.not. allocated(xpoints)) then
            allocate(xpoints(par%npoints+par%nrugauge))
            allocate(ypoints(par%npoints+par%nrugauge))
            xpoints=0
            ypoints=0

            ! Convert world coordinates of points to nearest (lsm) grid point
            call snappointstogrid(par, s, xpoints, ypoints)
         endif
      endif

      if (par%nrugauge>0) then
         allocate(rugrowindex(par%nrugauge))
         do i=1,par%nrugauge
            rugrowindex(i)=ypoints(par%npoints+i)
         enddo
      endif

      ! TIME-AVERAGE, VARIANCE and MIN-MAX ARRAYS

      if (par%nmeanvar>0) then
         !! First time file opening for time-average output
         do i=1,par%nmeanvar
            call makeaveragenames(chartoindex(par%meanvars(i)),fnamemean,fnamevar,fnamemin,fnamemax)
            call indextos(s,chartoindex(par%meanvars(i)),t)
            reclenm = wordsize
            select case(t%rank)
             case (2)
               reclenm = wordsize*size(t%r2)
             case (3)
               reclenm = wordsize*size(t%r3)
             case (4)
               reclenm = wordsize*size(t%r4)
            end select
            open(indextomeanunit(i),file=fnamemean, form='unformatted',access='direct',recl=reclenm,status='REPLACE')
            open(indextovarunit(i), file=fnamevar,  form='unformatted',access='direct',recl=reclenm,status='REPLACE')
            open(indextominunit(i), file=fnamemin,  form='unformatted',access='direct',recl=reclenm,status='REPLACE')
            open(indextomaxunit(i), file=fnamemax,  form='unformatted',access='direct',recl=reclenm,status='REPLACE')
         enddo
      endif ! par%nmeanvar > 0

      !
      ! drifter output files
      !
      if (par%ndrifter>0) then
         reclen=wordsize*3
         do i=1,par%ndrifter
            write(fname,'("drifter",i0.3,".dat")') i
            open(indextodrifterunit(i),file=fname,form='unformatted',access='direct',recl=reclen,status='REPLACE')
         enddo
      endif ! par%ndrifter >0

   end subroutine fortoutput_init

   subroutine checkfile(index,unit,reclen,jtg)
      implicit none
      integer, intent(in)  :: index,reclen
      integer, intent(out) :: unit,jtg
      logical              :: lopen
      character(len=1000)  :: filename

      unit = indextoglobalunit(index)
      inquire(unit=unit, opened=lopen)
      if ( .not. lopen ) then
         filename = trim(mnemonics(outnumbers(index)))//'.dat'
         open(unit, file=filename,form='unformatted',&
         &access='direct',recl=reclen)
      endif
      inquire(unit=unit,nextrec=jtg)
   end subroutine checkfile

   integer function outunit(ind,s)
      !
      ! given the type of output file 's' and the index 'ind'
      ! returns the appopriate file unit number
      !
      implicit none
      integer, intent(in)      :: ind
      character(*), intent(in) :: s

      integer, parameter :: offset = 10000
      select case(s)
       case('points')
         outunit = offset +              ind
       case('global')
         outunit = offset + 10*numvars + ind
       case('mean')
         outunit = offset + 20*numvars + ind
       case('min')
         outunit = offset + 30*numvars + ind
       case('max')
         outunit = offset + 40*numvars + ind
       case('var')
         outunit = offset + 50*numvars + ind
       case('drifter')
         outunit = offset + 60*numvars + ind
       case default
         print *,'internal error in outunit: no such type: ',trim(s)
         outunit = -1
         call halt_program
      end select
   end function outunit

   integer function indextopointsunit(index)
      implicit none
      integer, intent(in) :: index
      indextopointsunit = outunit(index,'points')
   end function indextopointsunit


   integer function indextoglobalunit(index)
      implicit none
      integer, intent(in) :: index
      indextoglobalunit = outunit(index,'global')
   end function indextoglobalunit

   integer function indextomeanunit(index)
      implicit none
      integer, intent(in) :: index
      indextomeanunit = outunit(index,'mean')
   end function indextomeanunit

   integer function indextominunit(index)
      implicit none
      integer, intent(in) :: index
      indextominunit = outunit(index,'min')
   end function indextominunit

   integer function indextomaxunit(index)
      implicit none
      integer, intent(in) :: index
      indextomaxunit = outunit(index,'max')
   end function indextomaxunit

   integer function indextovarunit(index)
      implicit none
      integer, intent(in) :: index
      indextovarunit = outunit(index,'var')
   end function indextovarunit

   integer function indextodrifterunit(index)
      implicit none
      integer, intent(in) :: index
      indextodrifterunit = outunit(index,'drifter')
   end function indextodrifterunit

   subroutine makeaveragenames(counter,fnamemean,fnamevar,fnamemin&
   &,fnamemax)
      use mnemmodule

      implicit none

      character(*)       :: fnamemean,fnamevar,fnamemin,fnamemax
      integer            :: counter

      fnamemean = trim(mnemonics(counter))//'_mean.dat'
      fnamevar  = trim(mnemonics(counter))//'_var.dat'
      fnamemin  = trim(mnemonics(counter))//'_min.dat'
      fnamemax  = trim(mnemonics(counter))//'_max.dat'

   end subroutine makeaveragenames


end module ncoutput_module
