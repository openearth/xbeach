program run_xbeach_tests

    use ftnunit
    use interp_tests

    implicit none

    call prepareTests

    call runtests_init

    ! Run tests here
    call runtests( allInterpTests )
    
    call runtests_final

    !!TODO: showResult will not be called, since runtests_final already stops the program
!!    call showResult

contains

!> Routine to start the testing
!! Note:
!!     This routine merely takes care that the unit tests
!!     are indeed run
subroutine prepareTests

    open( 10, file = 'ftnunit.run' )
    write( 10, '(a)' ) 'ALL'
    close( 10 )

end subroutine prepareTests

!> Start the browser to show the result
!!
subroutine showResult
    !character(len=1) :: answer
    !
    !write(*,*)     'Press ENTER ...'
    !read(*,'(a)' ) answer

!!    call system( 'ftnunit.html' )

end subroutine showResult

end program run_xbeach_tests
