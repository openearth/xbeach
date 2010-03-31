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

! global
integer, dimension(:), allocatable, save :: globalvarids
! default output (fixed length)
character(len=5), parameter :: default_mnems(7) = (/'H    ','zs   ','zs0  ','zb   ',&
     &'hh   ','u    ','v    '/) 
! One of these doesn't work:
!      &'ue   ','ve   ', 'urms ','Fx   ','Fy   ',&
!      &'ccg  ','ceqsg','ceqbg', &
!      &'E    ','R    ','D    ','DR   '/)

! points 
integer, save :: pointsdimid
integer, save :: xpointsvarid, ypointsvarid

! time 
integer, save :: globaltimedimid, pointtimedimid
integer, save :: globaltimevarid, pointtimevarid

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

  ! initialize values
  ! global
  allocate(globalvarids(par%nglobalvar))
  globalvarids = -1 ! initialize to 0
  outputg = .true.

  npointstotal = par%npoints+par%nrugauge
  outputp = (npointstotal .gt. 0) .and. (size(tpar%tpp) .gt. 0)

  ncfilename = 'xboutput.nc' 
  ! create a file
  status = nf90_create(path = ncfilename, cmode=NF90_CLOBBER, ncid = ncid)
  if (status /= nf90_noerr) call handle_err(status)

  ! dimensions
  ! grid
  status = nf90_def_dim(ncid, 'x', s%nx+1, xdimid) 
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, 'y', s%ny+1, ydimid)
  if (status /= nf90_noerr) call handle_err(status)
  ! wave angles
  status = nf90_def_dim(ncid, 'wave_angle', s%ntheta, thetadimid)
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

  ! define space & time variables
  ! grid
  status = nf90_def_var(ncid, 'x', NF90_DOUBLE, (/ xdimid /), xvarid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, xvarid, 'units', 'm')
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, xvarid, 'description', 'local x coordinate')
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_var(ncid, 'y', NF90_DOUBLE, (/ ydimid /), yvarid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, yvarid, 'units', 'm')
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_put_att(ncid, yvarid, 'description', 'local y coordinate')
  if (status /= nf90_noerr) call handle_err(status)



  ! global
  if (outputg) then
     status = nf90_def_var(ncid, 'globaltime', NF90_DOUBLE, (/ globaltimedimid /), globaltimevarid)
     if (status /= nf90_noerr) call handle_err(status)
     status = nf90_put_att(ncid, globaltimevarid, 'units', trim(par%tunits))
     if (status /= nf90_noerr) call handle_err(status)
     
     ! default global output variables
     do i=1,size(default_mnems)
        mnem = default_mnems(i)
        j = chartoindex(mnem)
        call indextos(s,j,t)
        select case(t%type)
        case('r')
           select case(t%rank)
           case(2)
              status = nf90_def_var(ncid, trim(mnem), NF90_DOUBLE, &
                   (/xdimid, ydimid, globaltimedimid /), globalvarids(i))
              if (status /= nf90_noerr) call handle_err(status)
           case(3)
              status = nf90_def_var(ncid, trim(mnem), NF90_DOUBLE,  &
                   (/xdimid, ydimid, thetadimid, globaltimedimid /), globalvarids(i))
              if (status /= nf90_noerr) call handle_err(status)
           case default
              write(0,*) 'mnem', mnem, ' not supported, rank:', t%rank
           end select
        case default
           write(0,*) 'mnem', mnem, ' not supported, type:', t%type
        end select
        status = nf90_put_att(ncid, globalvarids(i), 'units', trim(t%units))
        if (status /= nf90_noerr) call handle_err(status)
        status = nf90_put_att(ncid, globalvarids(i), 'description', trim(t%description))
        if (status /= nf90_noerr) call handle_err(status)
     end do
  end if
!  ! points


  ! done defining variables
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



  status = nf90_close(ncid)
  if (status /= nf90_noerr) call handle_err(status)
  
end subroutine ncoutput_init

subroutine nc_output(s,sl,par, tpar)
  use params
  use spaceparams
  use timestep_module
  use mnemmodule

  implicit none

  type(spacepars), intent(in)                         :: s,sl ! s-> spaceparams, what is isl?
  type(parameters), intent(in)                        :: par
  type(timepars), intent(in)                        :: tpar
  
  type(arraytype)                                     :: t
  integer                                      :: i,j,k,l,m,n
  character(len=10)                             :: mnem

  character, dimension(:), allocatable                               :: a
  integer :: status
  

  ! not a output step

  status = nf90_open(ncid=ncid, path=ncfilename, mode=nf90_write)
  if (status /= nf90_noerr) call handle_err(status)  


  if (tpar%outputg) then
     status = nf90_put_var(ncid, globaltimevarid, par%t, (/tpar%itg/))
     if (status /= nf90_noerr) call handle_err(status) 
     ! default global output variables
     do i=1,size(default_mnems)
        mnem = default_mnems(i)
        j = chartoindex(mnem)
        call indextos(s,j,t)
        select case(t%type)
        case('r')
           select case(t%rank)
           case(2)
              status = nf90_put_var(ncid, globalvarids(i), t%r2, start=(/1,1,tpar%itg/) )
              if (status /= nf90_noerr) call handle_err(status)
           case(3)
              status = nf90_put_var(ncid, globalvarids(i), t%r3, start=(/1,1,1, tpar%itg/) )
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


  status = nf90_close(ncid=ncid)
  if (status /= nf90_noerr) call handle_err(status) 
 

end subroutine nc_output

#endif
end module ncoutput_module
