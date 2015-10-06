

! Module to defer logging to a function that can be set using the set_logger function

module logging_module

  use iso_c_utils
  use typesandkinds
  use xmpi_module

  implicit none
  integer,save     :: logfileid
  integer,save     :: errorfileid
  integer,save     :: warningfileid
  !integer,save     :: pardatfileid

  abstract interface
     subroutine ilogger(level, msg)
       use iso_c_binding
       use iso_c_utils
       integer(c_int), value, intent(in) :: level !< severity
       character(c_char), intent(in) :: msg(MAXSTRINGLEN) !< c message null terminated
     end subroutine ilogger
  end interface


  procedure(ilogger), pointer :: logging_callback => null()

  ! Levels correspond to log4net and log4j

  integer, parameter, public :: LEVEL_ALL = 0
  integer, parameter, public :: LEVEL_DEBUG = 1
  integer, parameter, public :: LEVEL_INFO  = 2
  integer, parameter, public :: LEVEL_WARN  = 3
  integer, parameter, public :: LEVEL_ERROR = 4
  integer, parameter, public :: LEVEL_FATAL = 5
  integer, parameter, public :: LEVEL_OFF = 6


include 'writeloginterface.inc'

contains
  subroutine set_logger(c_callback) bind(C, name="set_logger")
    !DEC$ ATTRIBUTES DLLEXPORT::set_logger

    type(c_funptr), value :: c_callback

    ! Set a callback that will be cauled with new messages

    call c_f_procpointer(c_callback, logging_callback)
  end subroutine set_logger

  subroutine logmsg(level, msg)
    integer(c_int), intent(in) :: level
    character(len=*), intent(in) :: msg

    character(c_char)             :: c_string(MAXSTRINGLEN)


    if (associated(logging_callback)) then
       c_string = string_to_char_array(msg)
       call logging_callback(level, c_string)
    end if
  end subroutine logmsg

  subroutine start_logfiles(error)

    implicit none

    integer         :: error

    if (xmaster) then

       logfileid       = generate_logfileid()
       if (xmaster)    open(logfileid,     file='XBlog.txt',       status='replace')

       errorfileid     = generate_logfileid()
       if (xmaster)    open(errorfileid,   file='XBerror.txt',     status='replace')

       warningfileid   = generate_logfileid()
       if (xmaster)    open(warningfileid, file='XBwarning.txt',   status='replace')

       if (logfileid < 0 .or. errorfileid < 0 .or. warningfileid < 0) error = 1

    endif ! xmaster

    if (error==1) then
       write(*,*) 'Error: not able to open log file. Stopping simulation'
       stop
    endif

  end subroutine start_logfiles

  subroutine close_logfiles

    if (xmaster) then
       close(logfileid                         )
       close(errorfileid,      STATUS='DELETE' )
       close(warningfileid                     )
    endif

  end subroutine close_logfiles


  subroutine get_logfileid(lid,eid,wid)

    implicit none
    integer, intent(out)     :: lid,eid,wid

    lid = logfileid
    eid = errorfileid
    wid = warningfileid

  endsubroutine get_logfileid

  function generate_logfileid() result (tryunit)

    implicit none

    integer     :: tryunit,error
    logical     :: fileopen

    tryunit  = 98
    fileopen = .true.
    error    = 0

    do while (fileopen)
       inquire(tryunit,OPENED=fileopen)
       if (fileopen) then
          tryunit=tryunit-1
       endif
       if (tryunit<=10) then
          tryunit     = -1
          fileopen    = .false.
          return
       endif
    enddo

  end function generate_logfileid

  subroutine progress_indicator(initialize,curper,dper,dt)

    implicit none

    logical,intent(in)      :: initialize    ! initialize current progress indicator
    real*8,intent(in)       :: curper        ! current percentage done
    real*8,intent(in)       :: dper          ! steps in percentage between output
    real*8,intent(in)       :: dt            ! steps in time (s) between output
    ! whichever reached earlier (dper,dt) will determin output
    ! internal
    real*8,save             :: lastper,lastt
    real*8                  :: tnow
    integer*4               :: count,count_rate,count_max


    if (initialize) then
       lastper = 0.d0
       call system_clock (count,count_rate,count_max)
       lastt = dble(count)/count_rate
    else
       call system_clock (count,count_rate,count_max)
       tnow = dble(count)/count_rate
       if (curper>=lastper+dper .or. tnow>=lastt+dt) then
          call writelog('ls','(f0.1,a)',curper,'% done')
          if (curper>=lastper+dper) then
             lastper = curper-mod(curper,dper)
          else
             lastper = curper
          endif
          lastt = tnow
       endif
    endif


  end subroutine progress_indicator

  subroutine report_file_read_error(filename)

     implicit none

     character(*)    :: filename

     call writelog('lswe','','Error reading file ''',trim(filename),'''')
     call writelog('lswe','','Check file for incorrect decimal format, line breaks and tab characters')
     call halt_program

  end subroutine report_file_read_error

  subroutine writelog_startup()

    use xmpi_module
    implicit none

    character(len=8)                                :: date
    character(len=10)                               :: time
    character(len=5)                                :: zone

    ! subversion information
    include 'version.def'

    ! get current working directory (gcc only)
#ifdef HAVE_CONFIG_H
#include "config.h"
    character(slen)                              :: cwd
    call getcwd(cwd)
#endif

    include 'version.dat'

    call date_and_time(DATE=date, TIME=time, ZONE=zone)

    if (xmaster) then
       call writelog('ls','','**********************************************************')
       call writelog('ls','','                   Welcome to XBeach                      ')
       call writelog('ls','','                                                          ')
       call writelog('ls','','            version 1.22.',trim(Build_Revision)            )
       call writelog('ls','','            date ',trim(Build_Date)                        )
       call writelog('ls','','  URL: ',trim(Build_URL)                                   )
       call writelog('ls','','**********************************************************')
       call writelog('ls','','                                                          ')
       call writelog('ls','','Simulation started: YYYYMMDD    hh:mm:ss     time zone (UTC)')
       call writelog('ls','','                    '//date //'  '//time(1:2)//':'//time(3:4)//':'//time(5:6)//'     '//zone)
       call writelog('ls','','                                                          ')
#ifdef HAVE_CONFIG_H
       call writelog('ls','',' running in: ',cwd)
#endif
       call writelog('ls','','General Input Module')
#ifdef USEMPI
       call writelog('ls','','MPI version, running on ',xmpi_size,'processes')
#endif
    endif

  end subroutine writelog_startup


#ifdef USEMPI
   subroutine writelog_mpi(mpiboundary,error)
      use xmpi_module

      implicit none

      integer, intent(in)                             :: error
      integer, intent(in)                             :: mpiboundary

      if (xmaster) then
         if (error==1) then
            call writelog('elws','','Unknown mpi division ',mpiboundary)
            call halt_program
         elseif (error==2) then
            call writelog('elws','','Number of domains specified does not match available number of computation cores ')
            call halt_program
         else
            call writelog('ls','','processor grid: ',xmpi_m,' X ',xmpi_n)
         endif
      endif

   end subroutine writelog_mpi
#endif

  subroutine writelog_finalize(tbegin, n, t, nx, ny, t0, t01)

    use xmpi_module
    implicit none

    integer                                         :: n,nx,ny
    real*8                                          :: tbegin,tend
    real*8                                          :: t,duration,dt,performance
    real*8, optional                                :: t0,t01

#ifdef USEMPI
    real*8                                          :: t1
#endif

    if (xmaster) then

       call cpu_time(tend)

       duration    = tend-tbegin
       dt          = t/n
       performance = duration/(nx+1)/(ny+1)/n

       call writelog('ls','','Duration   : ',duration,' seconds'       )
       call writelog('ls','','Timesteps  : ',n                         )
       call writelog('ls','','Average dt : ',dt,' seconds'             )
       call writelog('ls','','Unit speed : ',performance,' seconds/1'  )

#ifdef USEMPI
       if (present(t0) .and. present(t01)) then
          t1 = MPI_Wtime()
          call writelog('ls','','MPI timing : procs      : ',xmpi_size                )
          call writelog('ls','','             seconds    : total : ',t1-t0, ' seconds')
          call writelog('ls','','                          loop  : ',t1-t01,' seconds')
       endif
#endif

       call writelog('ls','','End of program xbeach')
    endif
    
    call close_logfiles
    ! reset callback
    logging_callback => null()
  end subroutine writelog_finalize

  subroutine writelog_distribute(destination,display)

    implicit none

    character(*), intent(in) :: destination
    character(*), intent(in) :: display
    integer                  :: level
    
    logical                  :: has_logger
    
    has_logger = associated(logging_callback)
    
    if (xmaster) then
        level = 0
        if (scan(destination,'s')>0) then
            level = 1
        end if
        if (scan(destination,'l')>0) then
           level = 2
           if (.not. has_logger) then
              ! Don't log to screen if callback is associated
              write(6,*) trim(display)
           end if
           write(logfileid,*)     trim(display)
        end if
        if (scan(destination,'w')>0) then
           level = 3
           write(0,*) trim(display)
           write(warningfileid,*) trim(display)
        end if
        if (scan(destination,'e')>0) then
           level = 4
           write(0,*)   trim(display)
           write(errorfileid,*)   trim(display)
        end if
        call logmsg(level, trim(display))
    endif

  end subroutine writelog_distribute

  include 'writelog.inc'

end module logging_module

