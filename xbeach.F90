program xbeach
use libxbeach_module
use iso_c_binding
implicit none

integer                                             :: rc
real(c_double)                                      :: t, tstop

! Initialize program                                                          !
rc = 0
rc = init()

! Start simulation                                                            !
rc = getdoubleparameter("t", t, 1)
rc = getdoubleparameter("tstop", tstop, 5)
do while (t<tstop)
   rc = executestep()
   rc = getdoubleparameter("t", t, 1)
   ! output
   rc = outputext()
enddo

! Cleanup
rc = finalize()
end program
