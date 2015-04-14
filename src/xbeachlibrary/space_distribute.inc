    logical, optional, intent(in)       :: toall
    integer src,comm
    logical iamsrc
    src    = xmpi_master
    comm   = xmpi_comm
    iamsrc = xmaster
    if(present(toall)) then
      if(toall) then
        src    = xmpi_imaster
        comm   = xmpi_ocomm
        iamsrc = xomaster
      endif
    endif
      if (xomaster) return
!directions for vi vim: filetype=fortran : syntax=fortran
