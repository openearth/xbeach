module logging_module
   use typesandkinds
   use xmpi_module

   implicit none
   save

   integer     :: logfileid
   integer     :: errorfileid
   integer     :: warningfileid


   procedure(distributeloginterface), pointer :: distributelog => null()

   abstract interface

      subroutine distributeloginterface(code,message,len)
         implicit none

         integer, intent(in) :: code
         integer, intent(in) :: len
         character(*),intent(in) :: message

      end subroutine distributeloginterface

   end interface

   !
   ! Options for destiantion in writelog
   ! 's' = screen
   ! 'l' = log file
   ! 'e' = error file
   ! 'w' = warning file
   !
   ! Combinations also allowed, f.i.
   ! 'le' = log file and error file
   ! 'el' ditto
   ! 'sel' = screen, log file and error file
   !
   interface writelog
      module procedure writelog_a
      module procedure writelog_aa
      module procedure writelog_ai
      module procedure writelog_ia
      module procedure writelog_aaa
      module procedure writelog_aaaa
      module procedure writelog_aai
      module procedure writelog_aii
      module procedure writelog_aaai
      module procedure writelog_aaia
      module procedure writelog_aia
      module procedure writelog_aiaa
      module procedure writelog_aiaaa
      module procedure writelog_aiai
      module procedure writelog_aiaia
      module procedure writelog_aaiai
      module procedure writelog_aaaiai
      module procedure writelog_aiafa
      module procedure writelog_aiafaf
      module procedure writelog_aiaiai
      module procedure writelog_aiaiaia
      module procedure writelog_aiaiaf
      module procedure writelog_aiaiafa
      module procedure writelog_iiiii
      module procedure writelog_af
      module procedure writelog_aaf
      module procedure writelog_afa
      module procedure writelog_afaf
      module procedure writelog_afafa
      module procedure writelog_aaaf
      module procedure writelog_aafa
      module procedure writelog_afaaa
      module procedure writelog_aafaf
      module procedure writelog_aaafaf
      module procedure writelog_afafafaf
      module procedure writelog_illll
      module procedure writelog_fa
      module procedure writelog_afaiaaa
   end interface writelog



CONTAINS

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
         call writelog('ls','','            version 1.22.',trim(Build_Revision),' Kingsday release')
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

   subroutine writelog_finalize(tbegin, n, t, nx, ny &
#ifdef USEMPI
   , t0, t01 &
#endif
   )

      use xmpi_module
      implicit none

      integer                                         :: n,nx,ny
      real*8                                          :: tbegin,tend
      real*8                                          :: t,duration,dt,performance
#ifdef USEMPI
      real*8, optional                                :: t0,t01
#endif

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

   end subroutine writelog_finalize

   subroutine writelog_distribute(destination,display,lomaster)

      implicit none

      character(*), intent(in)      :: destination
      character(slen), intent(in)   :: display
      logical, intent(in), optional :: lomaster

      integer                 :: level
      logical                       :: doit

      ! wwvv: sometimes, this routine is called from xomaster, hence this extra code
      ! wwvv: see snappointstogrid
      ! wwvv: only implemented in writelog_a

      doit = xmaster

      if (present(lomaster)) then
         if (lomaster) then
            doit = xomaster
         endif
      endif

      if (doit) then
         level = -1
         if (scan(destination,'s')>0) then
            level = 3
         end if
         if (scan(destination,'l')>0) then
            level = 2
            write(6,*)     trim(display)
            write(logfileid,*)     trim(display)
         end if
         if (scan(destination,'w')>0) then
            level = 1
            write(0,*) trim(display)
            write(warningfileid,*) trim(display)
         end if
         if (scan(destination,'e')>0) then
            level = 0
            write(0,*)   trim(display)
            write(errorfileid,*)   trim(display)
         end if
         if (associated(distributelog)) then
            call distributelog(level,trim(display), len(trim(display)))
         endif
      endif

   end subroutine writelog_distribute

   subroutine writelog_a(destination,form,message_char,lomaster)
      implicit none
      character(*),intent(in)    ::  form,message_char
      character(*),intent(in)    ::  destination
      logical, intent(in), optional ::  lomaster
      character(slen)            ::  display

      if (form=='') then
         write(display,*)trim(message_char)
      else
         write(display,form)trim(message_char)
      endif

      call writelog_distribute(destination, display, lomaster)

   end subroutine writelog_a

   subroutine writelog_aa(destination,form,message_char1,message_char2)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2
      character(*),intent(in)       ::  destination
      character(slen)            ::  display

      if (form=='') then
         write(display,*) message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         trim(message_char2)
      else
         write(display,form) message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         trim(message_char2)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aa

   subroutine writelog_ai(destination,form,message_char,message_int)
      implicit none
      character(*),intent(in)    ::  form,message_char
      character(*),intent(in)    ::  destination
      integer,intent(in)         ::  message_int
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char(1:min(len(message_char),len_trim(message_char)+1)), &
         message_int
      else
         write(display,form)message_char(1:min(len(message_char),len_trim(message_char)+1)), &
         message_int
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_ai

   subroutine writelog_ia(destination,form,mint1,mchar1)
      implicit none
      character(*),intent(in)    ::  form,mchar1
      character(*),intent(in)    ::  destination
      integer,intent(in)         ::  mint1
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mint1,trim(mchar1)
      else
         write(display,form)mint1,trim(mchar1)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_ia

   subroutine writelog_aaa(destination,form,message_char1,message_char2,message_char3)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
      character(*),intent(in)       ::  destination
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         trim(message_char3)
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         trim(message_char3)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aaa

   subroutine writelog_aaaa(destination,form,message_char1,message_char2,message_char3,message_char4)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2,message_char3,message_char4
      character(*),intent(in)       ::  destination
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         trim(message_char4)
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         trim(message_char4)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aaaa


   subroutine writelog_aai(destination,form,message_char1,message_char2,message_int)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2
      character(*),intent(in)    ::  destination
      integer,intent(in)         ::  message_int
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_int
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_int
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aai

   subroutine writelog_aii(destination,form,message_char1,message_int1,message_int2)
      implicit none
      character(*),intent(in)    ::  form,message_char1
      character(*),intent(in)    ::  destination
      integer,intent(in)         ::  message_int1,message_int2
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_int1,message_int2
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_int1,message_int2
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aii

   subroutine writelog_aia(destination,form,message_char1b,message_intb,message_char2b)
      implicit none
      character(*),intent(in)    ::  form,message_char1b,message_char2b
      character(*),intent(in)    ::  destination
      integer,intent(in)         ::  message_intb
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1b(1:min(len(message_char1b),len_trim(message_char1b)+1)), &
         message_intb,trim(message_char2b)
      else
         write(display,form)message_char1b(1:min(len(message_char1b),len_trim(message_char1b)+1)), &
         message_intb,trim(message_char2b)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aia

   subroutine writelog_aaai(destination,form,message_char1,message_char2,message_char3,message_int)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
      character(*),intent(in)    ::  destination
      integer,intent(in)         ::  message_int
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_int
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_int
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aaai

   subroutine writelog_aaia(destination,formb,message_char1b,message_char2b,message_int,message_char3b)
      implicit none
      character(*),intent(in)    ::  formb,message_char1b,message_char2b,message_char3b
      character(*),intent(in)    ::  destination
      integer,intent(in)         ::  message_int
      character(slen)            ::  display

      if (formb=='') then
         write(display,*)message_char1b(1:min(len(message_char1b),len_trim(message_char1b)+1)), &
         message_char2b(1:min(len(message_char2b),len_trim(message_char2b)+1)), &
         message_int,trim(message_char3b)
      else
         write(display,formb)message_char1b(1:min(len(message_char1b),len_trim(message_char1b)+1)), &
         message_char2b(1:min(len(message_char2b),len_trim(message_char2b)+1)), &
         message_int,trim(message_char3b)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aaia

   subroutine writelog_aiaa(destination,form,message_char1b,message_intb,message_char2b,message_char3b)
      implicit none
      character(*),intent(in)    ::  form,message_char1b,message_char2b,message_char3b
      character(*),intent(in)       ::  destination
      integer,intent(in)         ::  message_intb
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1b(1:min(len(message_char1b),len_trim(message_char1b)+1)), &
         message_intb, &
         message_char2b(1:min(len(message_char2b),len_trim(message_char2b)+1)), &
         trim(message_char3b)
      else
         write(display,form)message_char1b(1:min(len(message_char1b),len_trim(message_char1b)+1)), &
         message_intb, &
         message_char2b(1:min(len(message_char2b),len_trim(message_char2b)+1)), &
         trim(message_char3b)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aiaa

   subroutine writelog_aiaaa(destination,form,message_char1b,message_intb,message_char2b,message_char3b,message_char4b)
      implicit none
      character(*),intent(in)    ::  form,message_char1b,message_char2b,message_char3b,message_char4b
      character(*),intent(in)       ::  destination
      integer,intent(in)         ::  message_intb
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1b(1:min(len(message_char1b),len_trim(message_char1b)+1)), &
         message_intb, &
         message_char2b(1:min(len(message_char2b),len_trim(message_char2b)+1)), &
         message_char3b(1:min(len(message_char3b),len_trim(message_char3b)+1)), &
         trim(message_char4b)
      else
         write(display,form)message_char1b(1:min(len(message_char1b),len_trim(message_char1b)+1)), &
         message_intb, &
         message_char2b(1:min(len(message_char2b),len_trim(message_char2b)+1)), &
         message_char3b(1:min(len(message_char3b),len_trim(message_char3b)+1)), &
         trim(message_char4b)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aiaaa


   subroutine writelog_aiai(destination,form,message_char1,message_int1,message_char2,message_int2)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2
      character(*),intent(in)       ::  destination
      integer,intent(in)         ::  message_int1,message_int2
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_int1, &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_int2
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_int1, &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_int2
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aiai

   subroutine writelog_aiaia(destination,form,mc1,mi1,mc2,mi2,mc3)
      implicit none
      character(*),intent(in)    ::  form,mc1,mc2,mc3
      character(*),intent(in)       ::  destination
      integer,intent(in)         ::  mi1,mi2
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mi1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mi2, &
         trim(mc3)
      else
         write(display,form)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mi1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mi2, &
         trim(mc3)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aiaia

   subroutine writelog_aaiai(destination,form,message_char1,message_char2,message_i1,message_char3,message_i2)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
      character(*),intent(in)       ::  destination
      integer*4,intent(in)          ::  message_i1,message_i2
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_i1, &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_i2
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_i1, &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_i2
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aaiai

   subroutine writelog_aaaiai(destination,form,message_char1,message_char2,message_char3,message_i1,message_char4,message_i2)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2,message_char3,message_char4
      character(*),intent(in)       ::  destination
      integer*4,intent(in)          ::  message_i1,message_i2
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_i1, &
         message_char4(1:min(len(message_char4),len_trim(message_char4)+1)), &
         message_i2
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_i1, &
         message_char4(1:min(len(message_char4),len_trim(message_char4)+1)), &
         message_i2
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aaaiai

   subroutine writelog_aiafa(destination,form,mc1,mi1,mc2,mf1,mc3)
      implicit none
      character(*),intent(in)    ::  form,mc1,mc2,mc3
      character(*),intent(in)       ::  destination
      integer*4,intent(in)          ::  mi1
      real*8,intent(in)             ::  mf1
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mi1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mf1,trim(mc3)
      else
         write(display,form)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mi1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mf1,trim(mc3)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aiafa

   subroutine writelog_aiafaf(destination,form,mc1,mi1,mc2,mf1,mc3,mf2)
      implicit none
      character(*),intent(in)    ::  form,mc1,mc2,mc3
      character(*),intent(in)       ::  destination
      integer*4,intent(in)          ::  mi1
      real*8,intent(in)             ::  mf1,mf2
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mi1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mf1, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         mf2
      else
         write(display,form)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mi1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mf1, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         mf2
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aiafaf

   subroutine writelog_aiaiai(destination,form,message_char1,message_i1,message_char2,message_i2,message_char3,message_i3)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
      character(*),intent(in)    ::  destination
      integer*4,intent(in)       ::  message_i1,message_i2,message_i3
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_i1, &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_i2, &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_i3
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_i1, &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_i2, &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_i3
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aiaiai

   subroutine writelog_aiaiaia(destination,form,mc1,message_i1,mc2,message_i2,mc3,message_i3, &
   mc4)
      implicit none
      character(*),intent(in)    ::  form,mc1,mc2,mc3,mc4
      character(*),intent(in)    ::  destination
      integer*4,intent(in)       ::  message_i1,message_i2,message_i3
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         message_i1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         message_i2, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         message_i3,trim(mc4)
      else
         write(display,form)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         message_i1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         message_i2, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         message_i3,trim(mc4)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aiaiaia

   subroutine writelog_aiaiaf(destination,form,mc1,mi1,mc2,mi2,mc3,mf1)
      implicit none
      character(*),intent(in)    ::  form,mc1,mc2,mc3
      character(*),intent(in)    ::  destination
      integer*4,intent(in)       ::  mi1,mi2
      real*8,intent(in)          ::  mf1
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mi1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mi2, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         mf1
      else
         write(display,form)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mi1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mi2, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         mf1
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aiaiaf

   subroutine writelog_aiaiafa(destination,form,mc1,mi1,mc2,mi2,mc3,mf1,mc4)
      implicit none
      character(*),intent(in)    ::  form,mc1,mc2,mc3,mc4
      character(*),intent(in)    ::  destination
      integer*4,intent(in)       ::  mi1,mi2
      real*8,intent(in)          ::  mf1
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mi1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mi2, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         mf1,trim(mc4)
      else
         write(display,form)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mi1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mi2, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         mf1,trim(mc4)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aiaiafa

   subroutine writelog_iiiii(destination,form,mi1,mi2,mi3,mi4,mi5)
      implicit none
      character(*),intent(in)    ::  form
      character(*),intent(in)    ::  destination
      integer*4,intent(in)       ::  mi1,mi2,mi3,mi4,mi5
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mi1,mi2,mi3,mi4,mi5
      else
         write(display,form)mi1,mi2,mi3,mi4,mi5
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_iiiii

   subroutine writelog_af(destination,form,message_char1,message_f1)
      implicit none
      character(*),intent(in)    ::  form, message_char1
      character(*),intent(in)    ::  destination
      real*8,intent(in)          ::  message_f1
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_f1
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_f1
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_af

   subroutine writelog_aaf(destination,form,message_char1,message_char2,message_f1)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2
      character(*),intent(in)       ::  destination
      real*8,intent(in)          ::  message_f1
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_f1
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_f1
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aaf

   subroutine writelog_afa(destination,form,mc1,mf1,mc2)
      implicit none
      character(*),intent(in)    ::  form,mc1,mc2
      character(*),intent(in)       ::  destination
      real*8,intent(in)          ::  mf1
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mf1,trim(mc2)
      else
         write(display,form)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mf1,trim(mc2)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_afa

   subroutine writelog_afaf(destination,form,message_char1,message_f1,message_char2,message_f2)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2
      character(*),intent(in)       ::  destination
      real*8,intent(in)          ::  message_f1,message_f2
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_f1, &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_f2
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_f1, &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_f2
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_afaf

   subroutine writelog_afafa(destination,form,mc1,mf1,mc2,mf2,mc3)
      implicit none
      character(*),intent(in)    ::  form,mc1,mc2,mc3
      character(*),intent(in)       ::  destination
      real*8,intent(in)          ::  mf1,mf2
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mf1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mf2,trim(mc3)
      else
         write(display,form)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mf1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mf2,trim(mc3)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_afafa

   subroutine writelog_aaaf(destination,form,message_char1,message_char2,message_char3,message_f1)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
      character(*),intent(in)       ::  destination
      real*8,intent(in)          ::  message_f1
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_f1
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_f1
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aaaf

   subroutine writelog_aafa(destination,form,message_char1b,message_char2b,message_f1b,message_char3b)
      implicit none
      character(*),intent(in)    ::  form,message_char1b,message_char2b,message_char3b
      character(*),intent(in)       ::  destination
      real*8,intent(in)          ::  message_f1b
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1b(1:min(len(message_char1b),len_trim(message_char1b)+1)), &
         message_char2b(1:min(len(message_char2b),len_trim(message_char2b)+1)), &
         message_f1b,trim(message_char3b)
      else
         write(display,form)message_char1b(1:min(len(message_char1b),len_trim(message_char1b)+1)), &
         message_char2b(1:min(len(message_char2b),len_trim(message_char2b)+1)), &
         message_f1b,trim(message_char3b)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aafa

   subroutine writelog_afaaa(destination,form,mc1a,mfa,mc2a,mc3a,mc4a)
      implicit none
      character(*),intent(in)    ::  form,mc1a,mc2a,mc3a,mc4a
      character(*),intent(in)       ::  destination
      real*8,intent(in)          ::  mfa
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1a(1:min(len(mc1a),len_trim(mc1a)+1)), &
         mfa, &
         mc2a(1:min(len(mc2a),len_trim(mc2a)+1)), &
         mc3a(1:min(len(mc3a),len_trim(mc3a)+1)), &
         trim(mc4a)
      else
         write(display,form)mc1a(1:min(len(mc1a),len_trim(mc1a)+1)), &
         mfa, &
         mc2a(1:min(len(mc2a),len_trim(mc2a)+1)), &
         mc3a(1:min(len(mc3a),len_trim(mc3a)+1)), &
         trim(mc4a)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_afaaa

   subroutine writelog_aafaf(destination,form,message_char1,message_char2,message_f1,message_char3,message_f2)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
      character(*),intent(in)       ::  destination
      real*8,intent(in)          ::  message_f1,message_f2
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_f1, &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_f2
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_f1, &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_f2
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aafaf

   subroutine writelog_aaafaf(destination,form,message_char1,message_char2,message_char3,message_f1,message_char4,message_f2)
      implicit none
      character(*),intent(in)    ::  form,message_char1,message_char2,message_char3,message_char4
      character(*),intent(in)       ::  destination
      real*8,intent(in)          ::  message_f1,message_f2
      character(slen)            ::  display

      if (form=='') then
         write(display,*)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_f1, &
         message_char4(1:min(len(message_char4),len_trim(message_char4)+1)), &
         message_f2
      else
         write(display,form)message_char1(1:min(len(message_char1),len_trim(message_char1)+1)), &
         message_char2(1:min(len(message_char2),len_trim(message_char2)+1)), &
         message_char3(1:min(len(message_char3),len_trim(message_char3)+1)), &
         message_f1, &
         message_char4(1:min(len(message_char4),len_trim(message_char4)+1)), &
         message_f2
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_aaafaf

   subroutine writelog_afafafaf(destination,form,mc1,mf1,mc2,mf2,mc3,mf3,mc4,mf4)
      implicit none
      character(*),intent(in)    ::  form,mc1,mc2,mc3,mc4
      character(*),intent(in)    ::  destination
      real*8,intent(in)          ::  mf1,mf2,mf3,mf4
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mf1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mf2, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         mf3, &
         mc4(1:min(len(mc4),len_trim(mc4)+1)), &
         mf4
      else
         write(display,form)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mf1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mf2, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         mf3, &
         mc4(1:min(len(mc4),len_trim(mc4)+1)), &
         mf4
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_afafafaf

   subroutine writelog_illll(destination,form,mi1,ml1,ml2,ml3,ml4)
      implicit none
      character(*),intent(in)    ::  form
      character(*),intent(in)    ::  destination
      integer*4,intent(in)       ::  mi1
      logical,intent(in)         ::  ml1,ml2,ml3,ml4
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mi1,ml1,ml2,ml3,ml4
      else
         write(display,form)mi1,ml1,ml2,ml3,ml4
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_illll

   subroutine writelog_fa(destination,form,mf1,mc1)
      implicit none
      character(*),intent(in)    ::  form, mc1
      character(*),intent(in)    ::  destination
      real*8,intent(in)          ::  mf1
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mf1,trim(mc1)
      else
         write(display,form)mf1,trim(mc1)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_fa

   subroutine writelog_afaiaaa(destination,form,mc1,mf1,mc2,mi1,mc3,mc4,mc5)
      implicit none
      character(*),intent(in)    ::  form, mc1,mc2,mc3,mc4,mc5
      character(*),intent(in)    ::  destination
      real*8,intent(in)          ::  mf1
      integer,intent(in)         ::  mi1
      character(slen)            ::  display

      if (form=='') then
         write(display,*)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mf1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mi1, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         mc4(1:min(len(mc4),len_trim(mc4)+1)), &
         trim(mc5)
      else
         write(display,form)mc1(1:min(len(mc1),len_trim(mc1)+1)), &
         mf1, &
         mc2(1:min(len(mc2),len_trim(mc2)+1)), &
         mi1, &
         mc3(1:min(len(mc3),len_trim(mc3)+1)), &
         mc4(1:min(len(mc4),len_trim(mc4)+1)), &
         trim(mc5)
      endif

      call writelog_distribute(destination, display)

   end subroutine writelog_afaiaaa

   subroutine assignlogdelegate_internal(fPtr)
      use iso_c_binding
      type(c_funptr), VALUE :: fPtr

      distributelog => null()
      if (c_associated(fPtr)) then
         call c_f_procpointer (fPtr, distributelog )
      endif

   end subroutine assignlogdelegate_internal
end module logging_module
