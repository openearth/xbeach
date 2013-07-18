module process_input
use xmpi_module
use iso_c_binding


    contains
    ! Reads command line input
    integer(c_int) function readinput()
    character(len=100) :: arg
    
    include 'version.def'
    include 'version.dat'
    
        n = iargc ()
        if (n > 0) then
            do iarg=1,n
                
                call getarg(iarg,arg)
                
                if (arg=='-i') then
                     write(*,*)'**********************************************************'
                     write(*,*)'                   Welcome to XBeach                      '
                     write(*,*)'                                                          '
                     write(*,*)'            revision ',trim(Build_Revision)             
                     write(*,*)'            date ',trim(Build_Date)                        
                     write(*,*)' URL: ',trim(Build_URL)                                    
                     write(*,*)'**********************************************************'
                     call halt_program
                endif
                
                if (arg=='-v') then
                     write(*,*)'revision ',trim(Build_Revision)
                     call halt_program
                endif
                
                if (arg=='?') then
                     write(*,*)'**********************************************************'
                     write(*,*)'                   Welcome to XBeach                      '
                     write(*,*)'                                                          '
                     write(*,*)'Usage:'
                     write(*,*)'    xbeach.exe'
                     write(*,*)'    xbeach.exe [options]'
                     write(*,*)' '
                     write(*,*)'Options:'
                     write(*,*)'    -i Shows information about the xbeach version on hand'
                     write(*,*)'    -v SHows the version of this xbeach executable'
                     write(*,*)'**********************************************************'
                     call halt_program
                endif
            enddo            
        endif
        readinput = 0
    end function readinput
    
end module process_input
    