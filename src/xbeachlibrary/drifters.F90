module drifter_module
  implicit none
  save
contains
  subroutine drifter(s,par)
    use params
    use spaceparams
    use readkey_module
    use xmpi_module
    use mnemmodule

    IMPLICIT NONE

    type(spacepars), target     :: s
    type(parameters)            :: par

    integer                     :: i
    integer                     :: ishift,jshift
    integer                     :: iu,ju,iv,jv
    real*8                      :: di,dj


#ifdef USEMPI
    ishift  = s%is(xmpi_rank+1)-1
    jshift  = s%js(xmpi_rank+1)-1
#else
    ishift  = 0
    jshift  = 0
#endif

    do i=1,par%ndrifter
       if (par%t>s%tdriftb(i) .and. par%t<s%tdrifte(i)) then

          ! determine closest u- and v-points
          iu          = nint(s%idrift(i)-ishift-0.d5)
          ju          = nint(s%jdrift(i)-jshift     )
          iv          = nint(s%idrift(i)-ishift     )
          jv          = nint(s%jdrift(i)-jshift-0.d5)

          ! update drifter if still inside domain
          if (    iu >= 1 .and. iu <= s%nx .and.  &
               ju >= 1 .and. ju <= s%ny .and.  &
               iv >= 1 .and. iv <= s%nx .and.  &
               jv >= 1 .and. jv <= s%ny            ) then

             ! determine movement of drifter relative to grid size
             di      = s%uu(iu,ju)/s%dsu(iu,ju)*par%dt
             dj      = s%vv(iv,jv)/s%dnv(iv,jv)*par%dt

             s%idrift(i) = s%idrift(i) + di
             s%jdrift(i) = s%jdrift(i) + dj

#ifdef USEMPI
          else
             s%idrift(i)  = huge(0.0d0)-10000.d0
             s%jdrift(i)  = huge(0.0d0)-10000.d0
#endif
          endif

#ifdef USEMPI
          call xmpi_allreduce(s%idrift(i),MPI_MIN)
          call xmpi_allreduce(s%jdrift(i),MPI_MIN)
#endif

       endif
    enddo
  end subroutine drifter

end module drifter_module
