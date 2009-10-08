!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   MODULE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this module will provide for netcdf output in the way varoutput provides for binary output. 
! It should be called optionally through one of the params.txt settings
! It will also be compiled conditionally like mpi. Only if it proves usefull will it be added by default. 
! it will add dependencies on the netcdf fortran library (http://www.unidata.ucar.edu/software/netcdf/)
! 
module ncoutputmod
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


implicit none
private
public init_output, var_output !what's this...

integer*4                           :: nglobalvar  ! number of global output variables
integer*4                           :: npoints     ! number of output points
integer*4                           :: nrugauge    ! number of runup gauges
integer*4,dimension(:),allocatable  :: nassocvar   ! vector with number of output variable per output point
integer*4                           :: nmeanvar    ! number of time-average variables

contains

  
  ! Error handling of netcdf errors 
  subroutine handle_err(status) 
    use netcdf
    
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine handle_err
  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!   INITIALISE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_output(s, sl, par, it)
  use params
  use spaceparams
  use readkey_module
  use netcdf

  implicit NONE
  integer :: ncid, status ! file id and status returned from a file operation

  type(spacepars), intent(in)                  :: s,sl ! why do we get 2 space params??
  type(parameters), intent(in)                 :: par 
  integer                                      :: it
  character(80)                                :: fname ! xy? h? zs?.... output.nc



  integer :: xdimid, ydimid, timedimid
  integer :: xvarid, yvarid, timevarid


  fname = "foo.nc" 
  ! create a file
  status = nf90_create(path = fname, cmode = nf90_noclobber, ncid = ncid)
  if (status /= nf90_noerr) call handle_err(status)

  ! dimensions
  status = nf90_def_dim(ncid, "xdim", s%nx, xdimid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "ydim", s%ny, ydimid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_dim(ncid, "timedim", nf90_unlimited, timedimid)
  if (status /= nf90_noerr) call handle_err(status)
  
  ! define space & time variables
  ! nf90_def_var(ncid, name, xtype, dimids, varid)
  status = nf90_def_var(ncid, "x", NF90_DOUBLE, (/ xdimid /), xvarid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_var(ncid, "y", NF90_DOUBLE, (/ ydimid /), yvarid)
  if (status /= nf90_noerr) call handle_err(status)
  status = nf90_def_var(ncid, "time", nf90_unlimited, (/ timedimid /), timevarid)
  if (status /= nf90_noerr) call handle_err(status)

  ! variables are saved on nglobal??, meanvar (time averaged?) 
  ! define output variables, not sure what "global" means in this context
  nglobalvar  = readkey_int ('params.txt','nglobalvar',   -1,       -1,     20) ! variables stored on grid time
  nmeanvar = readkey_int ('params.txt','nmeanvar',0,0,15) ! time averaged variables stored on grid
  !nassocvar =  ! point variables stored on point time, different variables can be stored per point? skip for now....
  
  
end subroutine init_output

subroutine var_output(it,s,sl,par)
  use params
  use spaceparams
  
  IMPLICIT NONE

  type(spacepars), intent(in)                         :: s,sl ! s-> spaceparams, what is isl?
  type(parameters), intent(in)                        :: par
  integer, intent(in)                                 :: it ! what is this

end subroutine var_output

end module ncoutputmod
