  logical, optional, intent(in) :: toall

  integer src,comm
  src  = xmpi_master
  comm = xmpi_comm
  if (present(toall)) then
    if(toall) then
      if(bcast_backtrace) then
        call backtrace
      endif
      src  = xmpi_imaster
      comm = xmpi_ocomm
    endif
  endif
  call xmpi_bcast(x,src,comm)
!directions for vi vim: filetype=fortran : syntax=fortran
