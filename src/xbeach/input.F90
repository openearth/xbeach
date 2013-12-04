module process_input
use iso_c_binding
use libxbeach_module

    contains
    ! Reads command line input
    integer(c_int) function readinput()
    
    character(len=100) :: arg
    character(len=500) :: version
    integer         :: narguments
    
    readinput = 0
    
        narguments = command_argument_count()
        if (narguments > 0) then
            do iarg=1,narguments
                
                call get_command_argument(iarg,arg)
                
                if (arg=='-V') then
                    call getversion(version)
                    write(*,*)'**********************************************************'
                    write(*,*)'You are using XBeach version 1.20.',trim(version),' December 2013 Sinterklaas release'
                    write(*,*)'**********************************************************'
                    readinput = 1
                endif
                
                if (arg.eq.'-h' .or. arg.eq.'--help') then
                    write(*,*)' '
                    write(*,*)'**********************************************************'
                    write(*,*)'                   Welcome to XBeach                      '
                    write(*,*)' '
                    write(*,*)'Usage:'
                    write(*,*)'    xbeach.exe'
                    write(*,*)'    xbeach.exe [options]'
                    write(*,*)' '
                    write(*,*)'Options:'
                    write(*,*)'    -V SHows the version of this xbeach executable'
                    write(*,*)'**********************************************************'
                    write(*,*)' '
                    readinput = 1
                endif
            enddo            
        endif
    end function readinput
    
end module process_input
    