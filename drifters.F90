module drifter_module
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

    include 's.ind'
    include 's.inp'
    
#ifdef USEMPI
    ishift  = s%is-1
    jshift  = s%js-1
#else
    ishift  = 0.d0
    jshift  = 0.d0
#endif

    do i=1,s%ndrift
        if (par%t>s%tdriftb(i) .and. par%t<s%tdrifte(i)) then
     
            ! determine closest u- and v-points
            iu          = nint(s%idrift(i)-ishift-0.d5)
            ju          = nint(s%jdrift(i)-jshift     )
            iv          = nint(s%idrift(i)-ishift     )
            jv          = nint(s%jdrift(i)-jshift-0.d5)
            
            ! determine movement of drifter relative to grid size
            di          = uu(iu,ju)/dsu(iu,ju)*par%dt
            dj          = vv(iv,jv)/dnv(iv,jv)*par%dt
            
            ! update drifter if still inside domain
            if (s%idrift(i)+di<=s%nx .and. s%jdrift(i)+dj<=s%ny) then
                s%idrift(i) = s%idrift(i) + di + ishift
                s%jdrift(i) = s%jdrift(i) + dj + jshift
#ifdef USEMPI
            else
               s%idrift(i)  = huge(0.0d0)
               s%jdrift(i)  = huge(0.0d0)
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
