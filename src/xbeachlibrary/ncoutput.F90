!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   MODULE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!#define NCREAL NF90_DOUBLE
#define CONVREAL
#define CONVREALTYPE real*8
!#endif

! NF90: macro to call nf90 function: if return value .ne. nf90_noerr,
!   an error message is produced, including file and lineno of the error,
!     and the program is halted
! example:
!  replace these two lines:
!     status = nf90_create(path = par%ncfilename, cmode=ior(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid = ncid)
!     if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
!  with this one line:
!  NF90(nf90_create(path = par%ncfilename, cmode=ior(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid = ncid))
!
!  NOTE: the NF90 call must be on one line because of preprocessor restrictions
!
#ifdef USENETCDF
#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)
#endif

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
   public ncoutput, fortoutput_init, points_output_init
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
   integer              :: driftersdimid, driftersdimid2, drifterstimedimid, drifterstimevarid
   integer, allocatable :: driftersvarids(:)
   real*8               :: drift(3)

   ! Ships
   integer :: shipdimid
   
   ! Q3D
   integer :: Q3Ddimid

   ! global
   integer, dimension(:), allocatable :: globalvarids
   ! default output (fixed length)

   ! points
   integer :: pointsdimid, pointnamelengthdimid
   integer :: xpointsvarid, ypointsvarid, pointtypesvarid, xpointindexvarid, ypointindexvarid, stationidvarid
   integer, dimension(:), allocatable :: pointsvarids
   integer, dimension(:),allocatable  :: xpoints     ! model x-coordinate of output points
   integer, dimension(:),allocatable  :: ypoints     ! model y-coordinate of output points

   ! mean
   ! number of variables by number of parameters per variable (mean, sigma^2, min, max)
   integer, dimension(:,:), allocatable       :: meanvarids
   character(slen), dimension(:), allocatable :: meanvartypes
   integer*4                                  :: nmeanvartypes  = 4   ! number of time-average variable types
   integer*4,dimension(:),allocatable         :: rugrowindex ! Array with row index where runup gauge can be found


   ! time
   integer :: globaltimedimid, pointtimedimid, meantimedimid
   integer :: globaltimevarid, pointtimevarid, meantimevarid


   ! TODO: check out why these are sometimes used....
   integer :: tidetimedimid, windtimedimid
   integer :: inoutdimid, tidecornersdimid

   ! local variables
   integer :: npointstotal
   integer :: itg,itp,itc,itm,itd
   ! Only alive at xmaster
   integer :: stpm        ! size of tpm
   integer :: wordsize    ! size of word in bytes

   integer                     :: noutnumbers = 0  ! the number of outnumbers
   integer, dimension(numvars) :: outnumbers  ! numbers, corrsponding to mnemonics, which are to be output

   ! Output type
   integer                     :: NCREAL ! can be NF90_double or NF90_float
contains

#ifdef USENETCDF
   ! Error handling of netcdf errors
   subroutine handle_err(status,file,line)
      use netcdf
      implicit none

      integer, intent ( in)    :: status
      character(*), intent(in) :: file
      integer, intent ( in)    :: line
      integer :: status2

      if(status /= nf90_noerr) then
         !UNIT=6 for stdout and UNIT=0 for stderr.
         write(0,'("NETCDF ERROR: ",a,i6,":",a)') file,line,trim(nf90_strerror(status))
         write(0,*) 'closing file'
         status2 = nf90_close(ncid)
         if (status2 /= nf90_noerr) then
            write(0,*) 'NETCDF ERROR: ', __FILE__,__LINE__,trim(nf90_strerror(status2))
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
      use params
      use paramsconst
      use spaceparams
      use timestep_module
      use mnemmodule
      use means_module
      use postprocessmod
      use logging_module
      ! This module is awaiting comments from Robert.
      use getkey_module

      implicit none

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

      !integer                                      :: npointstotal
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
      &           par%outputformat .eq. OUTPUTFORMAT_DEBUG
      donetcdf  = par%outputformat .eq. OUTPUTFORMAT_NETCDF .or. &
      &           par%outputformat .eq. OUTPUTFORMAT_DEBUG

      outputp = .false.

      ! initialize values

      ! set output precision for NetCDF
      if(par%outputprecision == OUTPUTPRECISION_SINGLE) then
         NCREAL = NF90_REAL
      else
         NCREAL = NF90_DOUBLE
      endif

      ! global

      ! store netcdf variable ids for each variable
      allocate(globalvarids(par%nglobalvar))
      globalvarids = -1 ! initialize to -1, so an error is raised when we miss something...
      outputg = .true.

      npointstotal = par%npoints+par%nrugauge
      outputp = (npointstotal .gt. 0) .and. (size(tpar%tpp) .gt. 0)
      allocate(pointsvarids(par%npointvar))

      allocate(driftersvarids(par%ndrifter))

      allocate(meanvarids(par%nmeanvar,nmeanvartypes))
      meanvarids = -1
      outputm = (par%nmeanvar .gt. 0)
      allocate(meanvartypes(nmeanvartypes))
      meanvartypes = (/ 'mean    ', 'var     ', 'min     ', 'max     '  /)

      ! create a file

      NF90(nf90_create(path = par%ncfilename, cmode=ior(NF90_CLOBBER,NF90_64BIT_OFFSET), ncid = ncid))

      ! dimensions TODO: only output dimensions that are used
      ! grid
      NF90(nf90_def_dim(ncid, 'globalx', s%nx+1, xdimid))
      NF90(nf90_def_dim(ncid, 'globaly', s%ny+1, ydimid))

      ! wave angles
      NF90(nf90_def_dim(ncid, 'wave_angle', s%ntheta, thetadimid))

      ! computational layers in bed ...
      NF90(nf90_def_dim(ncid, 'bed_layers', par%nd, bedlayersdimid))

      ! sediment classes
      NF90(nf90_def_dim(ncid, 'sediment_classes', par%ngd, sedimentclassesdimid))

      ! dimensions of length 2.... what is this.... TODO: find out what this is
      NF90(nf90_def_dim(ncid, 'inout', 2, inoutdimid))

      ! write(*,*) 'Writing ndrifter', par%ndrifter
      if (par%ndrifter .gt. 0) then
         ! create dimensions for drifters and drifterstime:
         NF90(nf90_def_dim(ncid, 'ndrifter',     par%ndrifter,   driftersdimid))
         NF90(nf90_def_dim(ncid, 'drifterstime', size(tpar%tpp), drifterstimedimid))
      end if

      ! write(*,*) 'Writing nship', par%nship
      if (par%nship .gt. 0) then
         NF90(nf90_def_dim(ncid, 'nship', par%nship, shipdimid))
      end if

            ! write(*,*) 'Writing nz', par%nz
      if (par%nz .gt. 1) then
         NF90(nf90_def_dim(ncid, 'nz', par%nz, Q3Ddimid))
      end if

      ! time dimensions are fixed, only defined if there are points
      if (outputg) then
         NF90(nf90_def_dim(ncid, 'globaltime', NF90_unlimited, globaltimedimid))
      end if
      if (outputp) then
         ! points
         NF90(nf90_def_dim(ncid, 'points', npointstotal, pointsdimid))
         NF90(nf90_def_dim(ncid, 'pointtime', size(tpar%tpp), pointtimedimid))
         NF90(nf90_def_dim(ncid, 'pointnamelength', 64, pointnamelengthdimid))
      end if
      if (outputm) then
         NF90(nf90_def_dim(ncid, 'meantime', size(tpar%tpm)-1, meantimedimid))
      end if

      if (s%tidelen > 0) then
         NF90(nf90_def_dim(ncid, 'tidetime', s%tidelen, tidetimedimid))
      endif

      if (par%tideloc > 0) then
         NF90(nf90_def_dim(ncid, 'tidecorners', par%tideloc, tidecornersdimid))
      endif

      if (s%windlen > 0) then
         NF90(nf90_def_dim(ncid, 'windtime', s%windlen, windtimedimid))
      endif

      ! define empty parameter variable
      NF90(nf90_def_var(ncid, 'parameter', NCREAL, varid=parvarid))
      ! define space & time variables
      ! grid
      NF90(nf90_def_var(ncid, 'globalx', NCREAL, (/ xdimid, ydimid /), xvarid))
      NF90(nf90_put_att(ncid, xvarid, 'units', 'm'))
      NF90(nf90_put_att(ncid, xvarid, 'long_name', 'local x coordinate'))
      NF90(nf90_put_att(ncid, xvarid, 'standard_name', 'projection_x_coordinate'))
      NF90(nf90_put_att(ncid, xvarid, 'axis', 'X'))
      ! For compatibility with CSIRO Dive software
      if (len(trim(par%projection)) .ne. 0)  then
         NF90(nf90_put_att(ncid, xvarid, 'projection', par%projection))
         NF90(nf90_put_att(ncid, xvarid, 'rotation',( s%alfa/atan(1.0d0)*45.d0)))
      end if

      NF90(nf90_def_var(ncid, 'globaly', NCREAL, (/ xdimid, ydimid /), yvarid))
      NF90(nf90_put_att(ncid, yvarid, 'units', 'm'))
      NF90(nf90_put_att(ncid, yvarid, 'long_name', 'local y coordinate'))
      NF90(nf90_put_att(ncid, yvarid, 'standard_name', 'projection_y_coordinate'))
      NF90(nf90_put_att(ncid, yvarid, 'axis', 'Y'))
      ! For compatibility with CSIRO Dive software
      if (len(trim(par%projection)) .ne. 0)  then
         NF90(nf90_put_att(ncid, yvarid, 'projection', par%projection))
         NF90(nf90_put_att(ncid, yvarid, 'rotation',( s%alfa/atan(1.0d0)*45.d0)))
      end if

      ! Some metadata attributes
      NF90(nf90_put_att(ncid,nf90_global, "Conventions", "CF-1.4"))
      NF90(nf90_put_att(ncid,nf90_global, "Producer", "XBeach littoral zone wave model (http://www.xbeach.org)"))
      NF90(nf90_put_att(ncid,nf90_global, "Build-Revision", trim(Build_Revision)))
      NF90(nf90_put_att(ncid,nf90_global, "Build-Date", trim(Build_Date)))
      NF90(nf90_put_att(ncid,nf90_global, "URL", trim(Build_URL)))

      ! Store all the parameters
      ! This part is awaiting comments from Robert McCall
      call getkeys(par, keys)
      do i=1,size(keys)
         rc = getkey(par, keys(i), val)
         if (val%type == 'i') then
            NF90(nf90_put_att(ncid, parvarid, keys(i), val%i0 ))
         elseif (val%type == 'c') then
            NF90(nf90_put_att(ncid, parvarid, keys(i), val%c0 ))
         elseif (val%type == 'r') then
            NF90(nf90_put_att(ncid, parvarid, keys(i), val%r0 ))
         end if
      end do

      ! global
      if (outputg) then
         NF90(nf90_def_var(ncid, 'globaltime', NCREAL, (/ globaltimedimid /), globaltimevarid))
         NF90(nf90_put_att(ncid, globaltimevarid, 'units', trim(par%tunits)))
         NF90(nf90_put_att(ncid, globaltimevarid, 'axis', 'T'))
         NF90(nf90_put_att(ncid, globaltimevarid, 'standard_name', 'time'))

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
               NF90(nf90_def_var(ncid, trim(mnem), NF90_INT, dimids, globalvarids(i)))
             case('r')
               NF90(nf90_def_var(ncid, trim(mnem), NCREAL, dimids, globalvarids(i)))
             case default
               write(0,*) 'mnem', mnem, ' not supported, type:', t%type
            end select
            NF90(nf90_put_att(ncid, globalvarids(i), 'coordinates', trim(coordinates)))
            deallocate(dimids)
            NF90(nf90_put_att(ncid, globalvarids(i), 'units', trim(t%units)))
            if (.not.(trim(t%standardname) .eq. '')) then
               NF90(nf90_put_att(ncid, globalvarids(i), 'standard_name', trim(t%standardname)))
            endif
            NF90(nf90_put_att(ncid, globalvarids(i), 'long_name', trim(t%description)))
         end do
      end if

      !  ! points
      ! default global output variables
      if (outputp) then
         NF90(nf90_def_var(ncid, 'pointtime', NCREAL, (/ pointtimedimid /), pointtimevarid))
         NF90(nf90_put_att(ncid, pointtimevarid, 'units', trim(par%tunits)))
         NF90(nf90_put_att(ncid, pointtimevarid, 'axis', 'T'))
         NF90(nf90_put_att(ncid, pointtimevarid, 'standard_name', 'time'))

         ! points
         NF90(nf90_def_var(ncid, 'pointx', NCREAL, (/ pointsdimid /), xpointsvarid))
         NF90(nf90_put_att(ncid, xpointsvarid, 'units', 'm'))
         NF90(nf90_put_att(ncid, xpointsvarid, 'long_name', 'local x coordinate'))
         NF90(nf90_put_att(ncid, xpointsvarid, 'standard_name', 'projection_x_coordinate'))
         NF90(nf90_put_att(ncid, xpointsvarid, 'axis', 'X'))


         NF90(nf90_def_var(ncid, 'pointy', NCREAL, (/ pointsdimid /), ypointsvarid))
         NF90(nf90_put_att(ncid, ypointsvarid, 'units', 'm'))
         NF90(nf90_put_att(ncid, ypointsvarid, 'long_name', 'local y coordinate'))
         NF90(nf90_put_att(ncid, ypointsvarid, 'standard_name', 'projection_y_coordinate'))
         NF90(nf90_put_att(ncid, ypointsvarid, 'axis', 'Y'))

         NF90(nf90_def_var(ncid, 'station_id', NF90_CHAR, (/ pointnamelengthdimid, pointsdimid /), stationidvarid))
         NF90(nf90_put_att(ncid, stationidvarid, 'long_name', 'station identification code'))
         NF90(nf90_put_att(ncid, stationidvarid, 'standard_name', 'station_id'))

         NF90(nf90_def_var(ncid, 'xpointindex', NF90_INT, (/ pointsdimid /), xpointindexvarid))
         ! wwvv above was NF90_DOUBLE
         NF90(nf90_put_att(ncid, xpointindexvarid, 'long_name', 'nearest x grid cell'))
         NF90(nf90_def_var(ncid, 'ypointindex', NF90_INT, (/ pointsdimid /), ypointindexvarid))
         ! wwvv above was NF90_DOUBLE
         NF90(nf90_put_att(ncid, ypointindexvarid, 'long_name', 'nearest y grid cell'))

         NF90(nf90_def_var(ncid, 'pointtypes', NF90_INT, (/ pointsdimid /), pointtypesvarid))
         ! wwvv above was NF90_DOUBLE
         NF90(nf90_put_att(ncid, pointtypesvarid, 'long_name', 'type of point (0=point, 1=rugauge)'))

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
               NF90(nf90_def_var(ncid, 'point_' // trim(mnem), NCREAL, dimids, pointsvarids(i)))
               NF90(nf90_put_att(ncid, pointsvarids(i), 'coordinates', trim(coordinates)))
               deallocate(dimids)
             case default
               write(0,*) 'mnem', mnem, ' not supported, type:', t%type
            end select
            NF90(nf90_put_att(ncid, pointsvarids(i), 'units', trim(t%units)))
            if (.not.(trim(t%standardname) .eq. '')) then
               NF90(nf90_put_att(ncid, pointsvarids(i), 'standard_name', trim(t%standardname)))
            endif
            NF90(nf90_put_att(ncid, pointsvarids(i), 'long_name', trim(t%description)))
         end do
      endif ! outputp

      if (par%ndrifter .gt. 0) then
         allocate(dimids(2))

         NF90(nf90_def_var(ncid, 'driftertime', NCREAL, (/ drifterstimedimid /), drifterstimevarid))
         NF90(nf90_put_att(ncid, drifterstimevarid, 'units', trim(par%tunits)))
         NF90(nf90_put_att(ncid, drifterstimevarid, 'axis', 'T'))
         NF90(nf90_put_att(ncid, drifterstimevarid, 'standard_name', 'time'))

         ! create netcdf variables for drifters:
         !  each variable is a 2-d array, containing the x-y values
         !  the netcdf name of e.g. the 3rd array is drifter_003

         NF90(nf90_def_dim(ncid, 'drifterstime2', 2, driftersdimid2))
         do i=1,par%ndrifter
            ! we need the following dimensions per drifter:
            !    driftersdimid2    : 2
            !    drifterstimedimid : number of output time steps
            dimids(1) = driftersdimid2
            dimids(2) = drifterstimedimid
            write(mnem,'("drifter_",I0.3)') i
            NF90(nf90_def_var(ncid, trim(mnem), NCREAL, dimids, driftersvarids(i)))
            NF90(nf90_put_att(ncid, driftersvarids(i), 'coordinates', 'pointx pointy'))
            NF90(nf90_put_att(ncid, driftersvarids(i), 'units', 'm'))
         enddo
         deallocate(dimids)

      endif  ! par%ndrifter

      if(outputm) then
         NF90(nf90_def_var(ncid, 'meantime', NCREAL, (/ meantimedimid /), meantimevarid))
         NF90(nf90_put_att(ncid, meantimevarid, 'units', trim(par%tunits)))
         NF90(nf90_put_att(ncid, meantimevarid, 'axis', 'T'))
         NF90(nf90_put_att(ncid, meantimevarid, 'standard_name', 'time'))
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
                  NF90(nf90_def_var(ncid, trim(t%name) // '_' // trim(cellmethod), NCREAL, dimids, meanvarids(i,j)))
                  NF90(nf90_put_att(ncid, meanvarids(i,j), 'coordinates', trim(coordinates)))
                  if (cellmethod .eq. 'var') then
                     NF90(nf90_put_att(ncid, meanvarids(i,j), 'units', '(' // trim(t%units) // ')^2'))
                  else
                     NF90(nf90_put_att(ncid, meanvarids(i,j), 'units', trim(t%units)))
                  endif
                  if (.not.(trim(t%standardname) .eq. '')) then
                     NF90(nf90_put_att(ncid, meanvarids(i,j), 'standard_name', trim(t%standardname)))
                  endif
                  NF90(nf90_put_att(ncid, meanvarids(i,j), 'long_name', trim(t%description)))
                  ! For H and urms we don't compute the mean but the rms of the rms.....
                  if (cellmethod .eq. 'mean' .and. ((t%name .eq. 'H') .or. (t%name .eq. 'urms')))  then
                     cellmethod = 'rms'
                  end if
                  NF90(nf90_put_att(ncid, meanvarids(i,j), 'cell_methods', 'meantime: ' // trim(cellmethod)))
               end do

               deallocate(dimids)
             case default
               write(0,*) 'mnem', mnem, ' not supported, type:', t%type
            end select
         end do
      endif  ! outputm

      ! done defining variables
      call writelog('ls', '', 'Writing file definition.')
      NF90(nf90_enddef(ncid))

      ! Fill meta variables
      ! Grid
      j = chartoindex('xz')
      call indextos(s,j,t)
      NF90(nf90_put_var(ncid, xvarid, CONVREAL(t%r2)))

      j = chartoindex('yz')
      call indextos(s,j,t)
      NF90(nf90_put_var(ncid, yvarid, CONVREAL(t%r2)))

      if (outputp) then
         call writelog('ls', '', 'Writing point vars.')
         NF90(nf90_put_var(ncid, xpointsvarid, CONVREAL(par%xpointsw)))
         NF90(nf90_put_var(ncid, ypointsvarid, CONVREAL(par%ypointsw)))
         NF90(nf90_put_var(ncid, pointtypesvarid, par%pointtypes))
         NF90(nf90_put_var(ncid, xpointindexvarid, xpoints))
         NF90(nf90_put_var(ncid, ypointindexvarid, ypoints))
         NF90(nf90_put_var(ncid, pointtypesvarid, par%pointtypes))
         NF90(nf90_put_var(ncid, stationidvarid, par%stationid(1:npointstotal)))
      end if

      NF90(nf90_close(ncid))
      ! wwvv sl is not used so it should be removed.
      ! wwvv for now, to avoid warning:

      if (sl%nx .ne. -1) return
   end subroutine ncoutput_init
#endif
   ! USENETCDF


   !___________________________________________________________________________________

   subroutine ncoutput(s,sl,par, tpar)
      use logging_module
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
      character(maxnamelen)                  :: mnem,sistermnemalloc
      real*8, dimension(:,:), allocatable    :: points

      ! some local variables to pass the data through the postprocessing function.
      integer                                       :: i0
      integer, dimension(:,:),          allocatable :: i2
      integer, dimension(:,:,:),        allocatable :: i3
      CONVREALTYPE                                  :: r0conv
      real*8,       dimension(:),       allocatable :: r1
      CONVREALTYPE, dimension(:),       allocatable :: r1conv
      real*8,       dimension(:,:),     allocatable :: r2
      CONVREALTYPE, dimension(:,:),     allocatable :: r2conv
      real*8,       dimension(:,:,:),   allocatable :: r3
      CONVREALTYPE, dimension(:,:,:),   allocatable :: r3conv
      real*8,       dimension(:,:,:,:), allocatable :: r4
      CONVREALTYPE, dimension(:,:,:,:), allocatable :: r4conv
      real*8,                           allocatable :: tempvectorr(:)
      real*8,dimension(size(tpar%tpg)+size(tpar%tpp)+size(tpar%tpc)+size(tpar%tpm)) :: outputtimes

      integer                                       :: jtg,reclen,unit,idumhl,ird,iru,xmax,xmin
      integer                                       :: iz,jz,idum
      real*8                                        :: di,dj,dx,dy
      integer                                       :: ilocal, jlocal
#ifdef USEMPI
      real*8, dimension(par%ndrifter)               :: idriftlocal,jdriftlocal
#endif
      integer                                       :: pii

#ifdef USENETCDF
      integer :: status
#endif

      logical :: dofortran, donetcdf, dofortran_compat
      logical :: dooutput_global, dooutput_mean, dooutput_point, dooutput_drifter

      type pointoutput
         integer                                  :: rank   ! rank of the data
         character(len=maxnamelen)                :: name   ! name of the variable
         real*8, dimension(:,:), allocatable      :: r2     ! contains 2-d data  (rank = 2)
         real*8, dimension(:,:,:), allocatable    :: r3     ! contains 3-d data  (rank = 3)
         real*8, dimension(:,:,:,:), allocatable  :: r4     ! contains 4-d data  (rank = 4)
         !                                                    note that r2 will be allocated (1,1)
         !                                                              r3 will be allocated (1,1,nnn)
         !                                                              r4 will be allocated  (1,1,mmm,nnn)
         !                                                    why ? because gridrotate demands this
      end type pointoutput

      integer                                     :: xpii, ypii
      integer                                     :: rugx, rugy

      type (pointoutput), dimension(:,:), allocatable :: pointoutputs

      real*8, dimension(:,:), allocatable :: runups,runups1   ! runups(:,i) will contain the 1+par%nrugdepth*3 runup values
      !                                                       ! for runup number i(i=1 .. par%nrugauge)
      integer, dimension(:), allocatable:: xpoints1  ! for netcdf runup output

      ! fortran output requested?
      dofortran = par%outputformat .eq. OUTPUTFORMAT_FORTRAN .or. &
      &           par%outputformat .eq. OUTPUTFORMAT_DEBUG

      ! for compatibility with older version, at some places
      ! fortran output is only done when tpar%output = .true.

      dofortran_compat = dofortran .and. tpar%output

      !netcdf output requested?
      donetcdf = par%outputformat .eq. OUTPUTFORMAT_NETCDF .or. &
      &          par%outputformat .eq. OUTPUTFORMAT_DEBUG

      ! time for global output?
      dooutput_global = tpar%outputg .and. par%nglobalvar .gt. 0

      ! time for mean output?
      dooutput_mean   = tpar%outputm .and. par%nmeanvar   .gt. 0 .and. tpar%itm .gt. 1

      ! time for point output?
      dooutput_point  = tpar%outputp .and. par%npointvar .gt. 0  &
      .and. par%npoints+par%nrugauge .gt. 0

      ! time for drifter output?
      dooutput_drifter  = tpar%outputp .and. par%ndrifter .gt. 0

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
            if (par%rotate==1) then
               sistermnemalloc = get_sister_mnem(mnem)

               select case (sistermnemalloc)
                case ('none')
                  ! nothing
                case default
                  call space_collect_mnem(s,sl,par,sistermnemalloc)
               end select

               select case(mnem)
                case(mnem_Sutot,mnem_Svtot)
                  call space_collect_mnem(s,sl,par,mnem_Subg)
                  call space_collect_mnem(s,sl,par,mnem_Svbg)
                  call space_collect_mnem(s,sl,par,mnem_Susg)
                  call space_collect_mnem(s,sl,par,mnem_Svsg)
               end select

            endif
         end do
         if (par%rotate==1) then
            call space_collect_mnem(s,sl,par,mnem_alfaz)
         endif
      endif

#endif
      ! USEMPI


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
      endif  !dooutput_point
#endif
      ! USEMPI

      ! If we're gonna write some mean output
      if(dooutput_mean) then

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
         if(xmaster) then
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

      ! wwvv a more efficient implementation:
      !      step 1: determine which process contains the desired point
      !      step 2: gridrotate this point
      !      step 3: store this point
      !      step 4: send the value of this point to xomaster (xmpi_send) after
      !              if( .not. xomaster) return
      !
      !      the points are assembled in pointoutputs, which is an 2-d array
      !        of type(pointoutput), which is described above.
      !
      !        pointoutputs(i,j) contains the values for pointvariable i
      !                                              and coordinates   j

      if (dooutput_point) then
         if (xomaster) then
            allocate(pointoutputs(par%npointvar,par%npoints+par%nrugauge))
         endif
         allocate(tempvectorr(1+par%nrugdepth*3))
         allocate(runups (size(tempvectorr),par%nrugauge))
         allocate(runups1(size(tempvectorr),par%nrugauge))
         ! WD: new code
         ! first: compute run gauge
         if(xcompute) then
            allocate(xpoints1(par%nrugauge))
            xpoints1 = s%nx+1
            ! the compute processes will determine their own runup values
            ! if no runp is found, a value of huge(0.0d0) is used
            ! The runup values are collected in xomaster, using xmpi_reduce,
            ! taking the minimum values
            do ii = par%npoints+1, par%npoints+par%nrugauge
               iru         = ii-par%npoints   !  iru will run from 1 to par%nrugauge
               tempvectorr = huge(0.d0)
               ! for netcdf output:
               do ird=1,par%nrugdepth
#ifdef USEMPI
                  call space_global_to_local(sl,1,rugrowindex(iru),rugx,rugy)
#else
                  rugx = 1
                  rugy = rugrowindex(iru)
#endif
                  ! rugy now contains the local y-coordinate
                  ! rugx is not used
                  ! if rugy is within the computational domain of this process,
                  ! we determine the local runup values
                  ! the computational domain in this process:
                  ! 1st dimension: sl%icls(xmpi_rank+1) .. sl%icle(xmpi_rank+1)
                  ! 2nd dimension: sl%jcls(xmpi_rank+1) .. sl%jcle(xmpi_rank+1)
#ifdef USEMPI
                  if (rugy .ge. sl%jcls(xmpi_rank+1) .and. rugy .le. sl%jcle(xmpi_rank+1)) then
                     xmin   = sl%icls(xmpi_rank+1)
                     xmax   = sl%icle(xmpi_rank+1)
#else
                  if (.true.) then
                     xmin   = 2
                     xmax   = sl%nx+1
#endif
                     idumhl = -1        ! Set default
                     if (rugrowindex(iru)>0) then  ! master domain always contains this runup gauge
                        ! local index of minimum location where hh<rugdepth
                        do j=xmin+1,xmax
                           if ((sl%hh(j,  rugy)<=par%rugdepth(ird)) .and. &
                           &   (sl%hh(j-1,rugy)> par%rugdepth(ird)) ) then
                              idumhl=j-1
                              exit
                           endif
                        enddo
                     endif

                     if (par%morfacopt==1) then
                        tempvectorr(1) = par%t*max(par%morfac,1.d0)
                     else
                        tempvectorr(1)=par%t
                     endif

                     ! for netcdf output:
                     if(ird == 1) then
                        ! convert from local coordinates to global, y-coordinate is not relevant:
                        if (idumhl .gt. 0) then
#ifdef USEMPI
                           call space_local_to_global(sl,idumhl,10,xpoints1(iru),idum)
#else
                           xpoints1(iru) = idumhl
                           idum           = 10
#endif
                        else
                           xpoints1(iru) = s%nx+1
                        endif
                     endif

                     if (idumhl .gt. 0) then
                        tempvectorr((ird-1)*3+2) = sl%xz(idumhl,rugy)
                        tempvectorr((ird-1)*3+3) = sl%yz(idumhl,rugy)
                        tempvectorr((ird-1)*3+4) = sl%zs(idumhl,rugy)
                     endif

                  endif  ! rugy within computational domain
               enddo     ! ird=1,par%nrugdepth
               runups(:,iru) = tempvectorr
            enddo        ! ii = par%npoints+1, par%npoints+par%nrugauge
            !  reduce the runups to runup1 on xmaster:
#ifdef USEMPI
            call xmpi_reduce(runups,runups1,MPI_MIN)
#else
            runups1 = runups
#endif

            ! reduce xpoints1 to (par%npoints+1:) on xmaster:
#ifdef USEMPI
            call xmpi_reduce(xpoints1,xpoints(par%npoints+1:),MPI_MIN)
#else
            xpoints(par%npoints+1:) = xpoints1
#endif
         endif ! xcompute
         !  send runups1 to xomaster who will receive it in runups
#ifdef USEMPI
         if (xomaster) then
            call xmpi_send(xmpi_imaster,xmpi_omaster,runups)
         else
            call xmpi_send(xmpi_imaster,xmpi_omaster,runups1)
         endif
         ! same for xpoints:
         call xmpi_send(xmpi_imaster,xmpi_omaster,xpoints(par%npoints+1:))
#else
         runups = runups1
#endif
         ! WD: /new code
         !  end runup gauge computations
         do i=1,par%npointvar
            mnem = par%pointvars(i)
            j = chartoindex(mnem)
            ! lookup the proper array
            call indextos(sl,j,t)
            ! get the proper output points ....
            select case(t%type)
             case('r')
               do ii = 1, par%npoints + par%nrugauge
                  ! wwvv above line should probably be
                  !  do ii = 1, par%npoints + par%nrugauge
                  !  have to check this with trunk
                  if(xomaster) then
                     pointoutputs(i,ii)%name = mnem
                     pointoutputs(i,ii)%rank = t%rank
                  endif
                  if (xomaster) then
                     xpii = xpoints(ii)
                     ypii = ypoints(ii)
                  endif
#ifdef USEMPI
                  call xmpi_bcast(xpii,xmpi_omaster,xmpi_ocomm)
                  call xmpi_bcast(ypii,xmpi_omaster,xmpi_ocomm)
                  call space_who_has(sl,xpii,ypii,pii)
                  if (xmpi_orank .eq. pii) then
                     call space_global_to_local(sl, xpii, ypii, ilocal,jlocal)
                  endif
#else
                  ilocal = xpii
                  jlocal = ypii
                  pii    = xmpi_orank
#endif
                  select case(t%rank)
                   case(2)
#ifdef USEMPI
                     if(xmpi_orank .eq. pii) then
                        allocate(r2(1,1))
                        call gridrotate(par, sl, t, r2, ilocal, jlocal)
                        call xmpi_send(pii,xmpi_omaster,r2)
                     elseif(xomaster) then
                        allocate(pointoutputs(i,ii)%r2(1,1))
                        call xmpi_send(pii,xmpi_omaster,pointoutputs(i,ii)%r2)
                     endif
                     if (xmpi_orank .eq. pii) deallocate(r2)
#else
                     allocate(pointoutputs(i,ii)%r2(1,1))
                     call gridrotate(par, sl, t, pointoutputs(i,ii)%r2, ilocal, jlocal)
#endif
                   case(3)
#ifdef USEMPI
                     if (xmpi_orank .eq. pii) then
                        allocate(r3(1,1,size(t%r3,3)))
                        call gridrotate(par, sl, t, r3, ilocal, jlocal)
                        call xmpi_send(pii,xmpi_omaster,r3)
                     elseif(xomaster) then
                        allocate(pointoutputs(i,ii)%r3(1,1,size(t%r3,3)))
                        call xmpi_send(pii,xmpi_omaster,pointoutputs(i,ii)%r3)
                     endif
                     if (xmpi_orank .eq. pii) deallocate(r3)
#else
                     allocate(pointoutputs(i,ii)%r3(1,1,size(t%r3,3)))
                     call gridrotate(par, sl, t, pointoutputs(i,ii)%r3, ilocal, jlocal)
#endif
                   case(4)
#ifdef USEMPI
                     if (xmpi_orank .eq. pii) then
                        allocate(r4(1,1,size(t%r4,3),size(t%r4,4)))
                        call gridrotate(t, r4, ilocal, jlocal)
                        call xmpi_send(pii,xmpi_omaster,r4)
                     elseif(xomaster) then
                        allocate(pointoutputs(i,ii)%r4(1,1,size(t%r4,3),size(t%r4,4)))
                        call xmpi_send(pii,xmpi_omaster,pointoutputs(i,ii)%r4)
                     endif
                     if (xmpi_orank .eq. pii) deallocate(r4)
#else
                     allocate(pointoutputs(i,ii)%r4(1,1,size(t%r3,3),size(t%r4,4)))
                     call gridrotate(t, pointoutputs(i,ii)%r4, ilocal, jlocal)
#endif
                   case default
                     write(0,*) 'Can''t handle rank: ', t%rank, ' of mnemonic', mnem
                  end select
               end do
             case default
               write(0,*) 'Can''t handle type: ', t%type, ' of mnemonic', mnem
            end select
         enddo ! i=1,par%npointvar
         !  Here follows go the code to determine the runup values and put them
         !    in runups(:,1:par%nrugauge)
         !
         !   This code assumes that hh and zs are available on xomaster
         !     we will change this asap
      endif ! dooutput_point

      ! writing is done by xomaster, the others processes go back to work

      if( .not. xomaster) return

      ! Open the output file

#ifdef USENETCDF
      if(donetcdf) then
         NF90(nf90_open(ncid=ncid, path=par%ncfilename, mode=nf90_write))
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
            NF90(nf90_put_var(ncid, globaltimevarid, CONVREAL(par%t*max(par%morfac,1.d0)), (/tpar%itg/)))
         endif
#endif
         ! write global output variables
         do i=1,par%nglobalvar
            mnem = par%globalvars(i)
            j    = chartoindex(mnem)
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
                     NF90(nf90_put_var(ncid, globalvarids(i), i0, start=(/1,tpar%itg/) ))
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
                     NF90(nf90_put_var(ncid, globalvarids(i), i2, start=(/1,1,tpar%itg/) ))
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
                     NF90(nf90_put_var(ncid, globalvarids(i), i3, start=(/1,1,1,tpar%itg/) ))
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
                     NF90(nf90_put_var(ncid, globalvarids(i), r0conv, start=(/1,tpar%itg/) ))
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
                     NF90(nf90_put_var(ncid, globalvarids(i), r1conv, start=(/1,tpar%itg/) ))
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
                  allocate(r2    (size(t%r2,1),size(t%r2,2)))
                  allocate(r2conv(size(t%r2,1),size(t%r2,2)))
                  call gridrotate(par, s, t, r2)
                  r2conv = CONVREAL(r2)
#ifdef USENETCDF
                  if(donetcdf) then
                     NF90(nf90_put_var(ncid, globalvarids(i), r2conv, start=(/1,1,tpar%itg/) ))
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
                  allocate(r3    (size(t%r3,1),size(t%r3,2),size(t%r3,3)))
                  allocate(r3conv(size(t%r3,1),size(t%r3,2),size(t%r3,3)))
                  call gridrotate(par, s, t, r3)
                  r3conv = CONVREAL(r3)
#ifdef USENETCDF
                  if(donetcdf) then
                     NF90(nf90_put_var(ncid, globalvarids(i), r3conv, start=(/1,1,1, tpar%itg/) ))
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
                  allocate(r4    (size(t%r4,1),size(t%r4,2),size(t%r4,3),size(t%r4,4)))
                  allocate(r4conv(size(t%r4,1),size(t%r4,2),size(t%r4,3),size(t%r4,4)))
                  call gridrotate(t, r4)
                  r4conv = CONVREAL(r4)
#ifdef USENETCDF
                  if(donetcdf) then
                     NF90(nf90_put_var(ncid, globalvarids(i), r4conv, start=(/1,1,1,1, tpar%itg/) ))
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
      end if   ! dooutput_global


      if(dooutput_point) then
         itp = itp+1

#ifdef USENETCDF
         if(donetcdf) then
            NF90(nf90_put_var(ncid, pointtimevarid, CONVREAL(par%t*max(par%morfac,1.d0)), (/tpar%itp/)))
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
            j    = chartoindex(mnem)
            ! lookup the proper array
            call indextos(s,j,t)
            ! get the proper output points ....
            ! I have no idea what is happening in varouput so I'll try it in a different way
            !TODO This is not very efficient because we are using the outer counters, reorder dimensions....
            !
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

#ifdef USENETCDF
                     if (donetcdf) then
                        NF90(nf90_put_var(ncid, pointsvarids(i), CONVREAL(pointoutputs(i,ii)%r2(1,1)), start=(/ii,tpar%itp/) ))
                     endif
#endif
                     if(dofortran_compat) then
                        if (ii .le. par%npoints) then
                           !points(ii,i+1) = r2(xpoints(ii), ypoints(ii))
                           points(ii,i+1) = pointoutputs(i,ii)%r2(1,1)
                        endif
                     endif
                   case(3)
#ifdef USENETCDF
                     if (donetcdf) then
                        NF90(nf90_put_var(ncid, pointsvarids(i), CONVREAL(r3(xpoints(ii), ypoints(ii),:)), start=(/ii,1,tpar%itp/)))
                     endif
#endif
                     if(dofortran_compat) then
                        if (ii .le. par%npoints) then
                           points(ii,i+1) = r3(xpoints(ii), ypoints(ii),1)    ! wwvv todo
                        endif
                     endif
                   case(4)
#ifdef USENETCDF
                     if (donetcdf) then
                        status      = nf90_put_var(ncid, pointsvarids(i), &
                        &                          CONVREAL(r4(xpoints(ii), ypoints(ii),:,:)), &
                        &                          start=(/ii,1,1,tpar%itp/) )
                        if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                     endif
#endif
                     if(dofortran_compat) then
                        if (ii .le. par%npoints) then
                           points(ii,i+1) = r4(xpoints(ii), ypoints(ii),1,1)    ! wwvv todo
                        endif
                     endif
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
                  points(ii,1) = par%t*max(par%morfac,1.d0)
               else
                  points(ii,1) = par%t
               endif
            enddo
            do ii = 1,par%npoints
               write(indextopointsunit(ii),rec=tpar%itp)CONVREAL(points(ii,:))
               call flush(indextopointsunit(ii))
            enddo
            deallocate(points)

            ! WD: new code
            do ii=1,par%nrugauge
               write(indextopointsunit(ii+par%npoints),rec=tpar%itp)CONVREAL(runups(:,ii))
               call flush(indextopointsunit(i+par%npoints))
            enddo
            ! WD: /new code
         endif ! dofortran
      endif   ! dooutput_point


#ifdef USENETCDF
      if(dooutput_mean .and. donetcdf) then
         ! Store the time (in morphological time)
         NF90(nf90_put_var(ncid, meantimevarid, CONVREAL(par%t*max(par%morfac,1.d0)), (/tpar%itm-1/)))
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
                           status      = nf90_put_var(ncid, meanvarids(i,j), &
                           &                          CONVREAL(sqrt(meansparsglobal(i)%variancesquareterm2d)), &
                           &                          start=(/1,1,tpar%itm-1/) )
                           if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                        elseif (t%name .eq. 'thetamean') then
                           status = nf90_put_var(ncid, meanvarids(i,j), &
                           &                     CONVREAL( &
                           &                     mod(2.d0*par%px &
                           &                     + atan2(nint(meansparsglobal(i)%mean2d)/1d7, &
                           &                     mod(meansparsglobal(i)%mean2d,1.d0)*1d1), 2.d0*par%px) / par%px * 180), &
                           &                     start=(/1,1,tpar%itm-1/) )
                           if (status /= nf90_noerr) call handle_err(status,__FILE__,__LINE__)
                        else
                           NF90(nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%mean2d),start=(/1,1,tpar%itm-1/)))
                        end if
                      case('var')
                        NF90(nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%variance2d),start=(/1,1,tpar%itm-1/)))
                      case('min')
                        NF90(nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%min2d),start=(/1,1,tpar%itm-1/)))
                      case('max')
                        NF90(nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%max2d),start=(/1,1,tpar%itm-1/)))
                      case default
                        write(0,*) 'Can''t handle cell method: ', trim(meanvartypes(j)), ' of mnemonic', trim(t%name)
                     end select
                   case(3)
                     select case(meanvartypes(j))
                      case('mean')
                        NF90(nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%mean3d), start=(/1,1,1,tpar%itm-1/) ))
                      case('var')
                        NF90(nf90_put_var(ncid,meanvarids(i,j),CONVREAL(meansparsglobal(i)%variance3d),start=(/1,1,1,tpar%itm-1/)))
                      case('min')
                        NF90(nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%min3d), start=(/1,1,1,tpar%itm-1/) ))
                      case('max')
                        NF90(nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%max3d), start=(/1,1,1,tpar%itm-1/) ))
                      case default
                        write(0,*) 'Can''t handle cell method: ', trim(meanvartypes(j)), ' of mnemonic', trim(t%name)
                     end select
                   case(4)
                     select case(meanvartypes(j))
                      case('mean')
                        NF90(nf90_put_var(ncid, meanvarids(i,j),CONVREAL(meansparsglobal(i)%mean4d), start=(/1,1,1,1,tpar%itm-1/)))
                      case('var')
                        NF90(nf90_put_var(ncid,meanvarids(i,j),CONVREAL(meansparsglobal(i)%variance4d),start=(/1,1,1,1,tpar%itm-1/)))
                      case('min')
                        NF90(nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%min4d), start=(/1,1,1,1,tpar%itm-1/)))
                      case('max')
                        NF90(nf90_put_var(ncid, meanvarids(i,j), CONVREAL(meansparsglobal(i)%max4d), start=(/1,1,1,1,tpar%itm-1/)))
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
                  &  CONVREAL( &
                  &            mod(2.d0*par%px + atan2(nint(meansparsglobal(i)%mean2d)/1d7, &
                  &            mod(meansparsglobal(i)%mean2d,1.d0)*1d1), 2.d0*par%px) / par%px * 180 &
                  &  )
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
         par%tintm = tpar%tpm(min(itm+2,stpm))-tpar%tpm(itm+1)  ! Next averaging period (min to stop array out of bounds)
         par%tintm = max(par%tintm,tiny(0.d0))                  ! to prevent par%tintm=0 after last output
      endif  ! dooutput_mean .and. dofortran_compat


      if (dooutput_drifter) then
         itd = itd+1
#ifdef USENETCDF
         if (donetcdf) then
            ! output time:
            NF90(nf90_put_var(ncid, drifterstimevarid, CONVREAL(par%t), start=(/itd/)))
         endif
#endif
         do i=1,par%ndrifter
            if (  par%t>=s%tdriftb(i) .and. par%t<=s%tdrifte(i) .and. &
            &     s%idrift(i)>1       .and. s%idrift(i)<=s%nx   .and. &
            &     s%jdrift(i)>1       .and. s%jdrift(i)<=s%ny             ) then

               iz = int(s%idrift(i))
               jz = int(s%jdrift(i))

               di = mod(s%idrift(i),1.d0)
               dj = mod(s%jdrift(i),1.d0)

               dx = di*s%dsu(iz,jz)*cos(s%alfaz(iz,jz)) - &
               &    dj*s%dnv(iz,jz)*sin(s%alfaz(iz,jz))
               dy = di*s%dsu(iz,jz)*sin(s%alfaz(iz,jz)) + &
               &    dj*s%dnv(iz,jz)*cos(s%alfaz(iz,jz))

               drift(1) = s%xz(iz,jz)+dx
               drift(2) = s%yz(iz,jz)+dy
               drift(3) = par%t
            else
               drift(1) = -999
               drift(2) = -999
               drift(3) = par%t
            endif

            if (dofortran_compat) then
               write(indextodrifterunit(i),rec=itd) CONVREAL(drift)
               call flush(indextodrifterunit(i))
            endif
#ifdef USENETCDF
            if (donetcdf) then
               ! output x,y:
               NF90(nf90_put_var(ncid,driftersvarids(i),CONVREAL(drift(1:2)),start=(/1,itd/)))
            endif
#endif

         enddo
      endif
#ifdef USENETCDF
      if(donetcdf) then
         NF90(nf90_close(ncid=ncid))
      endif
#endif

      if(dofortran_compat) then
         outputtimes                                = -999.d0
         outputtimes(1:itg)                         = tpar%tpg(1:itg)
         outputtimes(itg+1:itg+itp)                 = tpar%tpp(1:itp)
         outputtimes(itg+itp+1:itg+itp+itc)         = tpar%tpc(1:itc)
         outputtimes(itg+itp+itc+1:itg+itp+itc+itm) = tpar%tpm(2:itm+1)          ! mean output always shifted by 1

         if (par%morfacopt==1) outputtimes=outputtimes*max(par%morfac,1.d0)

         open(1998,file='dims.dat',form='unformatted',access='direct',recl=wordsize*(10+size(outputtimes)))
         write(1998,rec=1) CONVREAL(itg*1.d0),&
         &                 CONVREAL(s%nx*1.d0),&
         &                 CONVREAL(s%ny*1.d0),&
         &                 CONVREAL(s%ntheta*1.d0),&
         &                 CONVREAL(par%kmax*1.d0),&
         &                 CONVREAL(par%ngd*1.d0),&
         &                 CONVREAL(par%nd*1.d0), &
         &                 CONVREAL(itp*1.d0),&
         &                 CONVREAL(itc*1.d0),&
         &                 CONVREAL(itm*1.d0),&
         &                 CONVREAL(outputtimes)
         call flush(1998)
      endif
      ! wwvv avoid warning about unused sl:
      if (sl%nx .eq. -1) return
   end subroutine ncoutput

#ifdef USENETCDF
   character(slen) function dimensionnames(dimids)
      implicit none
      integer, dimension(:), intent(in)           :: dimids ! store the dimids in a vector

      integer :: i
      character(slen)  :: dimensionname
      ! combine all the dimensionnames
      ! assumes all dimensions have an accompanying variable that should be used for coordinates.
      ! ",".join would have been nice here....
      dimensionnames = ''
      ! Fortran array dimensions are in reverse order
      do i=size(dimids),2,-1
         NF90(nf90_inquire_dimension(ncid, dimids(i), name=dimensionname))
         dimensionnames = trim(dimensionnames) // trim(dimensionname) // ','
      end do
      NF90(nf90_inquire_dimension(ncid, dimids(1), name=dimensionname))
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
       case('par%nd')
         dimensionid = bedlayersdimid
       case('par%ndrifter')
         dimensionid = driftersdimid
       case('par%nship')
         dimensionid = shipdimid
       case('par%nz')
         dimensionid = Q3Ddimid
       case default
         call writelog('els','','Unknown dimension expression:'  // expression)
         stop 1
      end select
   end function dimensionid
#endif
   ! USENETCDF

   subroutine points_output_init(s,par)
      ! this initialize things for point output
      ! has to be called before fourtoutput_init and ncoutput_init
      use spaceparams
      use params
      use postprocessmod

      type(spacepars), intent(in) :: s
      type(parameters),intent(in) :: par

      allocate(rugrowindex(par%nrugauge))

      allocate(xpoints(par%npoints+par%nrugauge))
      allocate(ypoints(par%npoints+par%nrugauge))

      ! Convert world coordinates of points to nearest (lsm) grid point
      if(xomaster) then
         call snappointstogrid(par, s, xpoints, ypoints)
      endif

      if(xomaster) then
         if (par%nrugauge>0) then
            rugrowindex = ypoints(par%npoints+1:)
            !do i=1,par%nrugauge
            !   rugrowindex(i)=ypoints(par%npoints+i)
            !enddo
         endif
      endif
#ifdef USEMPI
      call xmpi_bcast(rugrowindex,xmpi_omaster,xmpi_ocomm)
      call xmpi_bcast(xpoints,    xmpi_omaster,xmpi_ocomm)
      call xmpi_bcast(ypoints,    xmpi_omaster,xmpi_ocomm)
#endif
   end subroutine points_output_init

   subroutine fortoutput_init(s,par,tpar)
      use params
      use spaceparams
      use readkey_module
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
         fname = trim(par%stationid(i)) // '.dat'
         if (par%pointtypes(i)==0) then
            reclenp=wordsize*(par%npointvar+1)
         else
            reclenp=wordsize*(1+par%nrugdepth*3)
         endif
         open(indextopointsunit(i),file=fname,&
         &    form='unformatted',access='direct',recl=reclenp,status='REPLACE')
      enddo
      if (par%npoints>0) then
         ! write index file of point output variables
         fid = create_new_fid()
         open(fid,file='pointvars.idx',status='replace',action='write')
         do i=1,par%npointvar
            write(fid,*)trim(par%pointvars(i))
         enddo
         close(fid)
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
         &    access='direct',recl=reclen)
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
   &                           ,fnamemax)
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
