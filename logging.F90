MODULE logging_module

use xmpi_module

implicit none

integer,save     :: logfileid
integer,save     :: errorfileid
!integer,save     :: pardatfileid


!
! Options for destiantion in writelog
! 's' = screen
! 'l' = log file
! 'e' = error file
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
   module procedure writelog_aaaiai
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
   module procedure writelog_aaaf
   module procedure writelog_aafa
   module procedure writelog_afaaa
   module procedure writelog_aafaf
   module procedure writelog_aaafaf
   module procedure writelog_afafafaf
   module procedure writelog_illll
end interface writelog

CONTAINS

  subroutine start_logfiles(error)
    implicit none
    integer    :: error
    integer    :: tryunit = 98
    logical    :: fileopen

    if (xmaster) then 
       fileopen = .true.    
       error    = 0
       do while (fileopen)
          inquire(tryunit,OPENED=fileopen)
          if (fileopen) then
             tryunit=tryunit-1
          endif
          if (tryunit<=10) then 
             tryunit = -1
             error  = 1
             fileopen = .false.
             return
          endif
       enddo
       logfileid = tryunit
       if (xmaster) then
          open(logfileid,file='XBlog.txt',status='replace')
       endif

       fileopen = .true.    
       error    = 0
       do while (fileopen)
          inquire(tryunit,OPENED=fileopen)
          if (fileopen) then
             tryunit=tryunit-1
          endif
          if (tryunit<=10) then 
             tryunit = -1
             error  = 1
             fileopen = .false.
             return
          endif
       enddo
       errorfileid = tryunit
       if (xmaster) then
          open(errorfileid,file='XBerror.txt',status='replace')
       endif
    endif ! xmaster


  end subroutine start_logfiles

subroutine close_logfiles

if (xmaster) then 
close(logfileid)
close(errorfileid, STATUS ='DELETE')
endif

end subroutine close_logfiles


subroutine get_logfileid(lid,eid)
   implicit none
   integer    :: lid,eid

   lid = logfileid
   eid = errorfileid
endsubroutine


subroutine writelog_a(destination,form,message_char)
   implicit none
   character(*),intent(in)    ::  form,message_char
   character(*),intent(in)       ::  destination
   character(256)             ::  display

   if (form=='') then
      write(display,*)message_char
   else
      write(display,form)message_char
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_a

subroutine writelog_aa(destination,form,message_char1,message_char2)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2
   character(*),intent(in)       ::  destination
   character(256)             ::  display
   
   if (form=='') then
      write(display,*)message_char1,message_char2
   else
      write(display,form)message_char1,message_char2
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aa

subroutine writelog_ai(destination,form,message_char,message_int)
   implicit none
   character(*),intent(in)    ::  form,message_char
   character(*),intent(in)    ::  destination
   integer,intent(in)         ::  message_int
   character(256)             ::  display

   if (form=='') then
      write(display,*)message_char,message_int
   else
      write(display,form)message_char,message_int
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_ai

subroutine writelog_ia(destination,form,mint1,mchar1)
   implicit none
   character(*),intent(in)    ::  form,mchar1
   character(*),intent(in)    ::  destination
   integer,intent(in)         ::  mint1
   character(256)             ::  display

   if (form=='') then
      write(display,*)mint1,mchar1
   else
      write(display,form)mint1,mchar1
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_ia

subroutine writelog_aaa(destination,form,message_char1,message_char2,message_char3)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
   character(*),intent(in)       ::  destination
   character(256)             ::  display
   
   if (form=='') then
      write(display,*)message_char1,message_char2,message_char3
   else
      write(display,form)message_char1,message_char2,message_char3
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aaa

subroutine writelog_aaaa(destination,form,message_char1,message_char2,message_char3,message_char4)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2,message_char3,message_char4
   character(*),intent(in)       ::  destination
   character(256)             ::  display
   
   if (form=='') then
      write(display,*)message_char1,message_char2,message_char3,message_char4
   else
      write(display,form)message_char1,message_char2,message_char3,message_char4
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aaaa


subroutine writelog_aai(destination,form,message_char1,message_char2,message_int)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2
   character(*),intent(in)    ::  destination
   integer,intent(in)         ::  message_int
   character(256)             ::  display

   if (form=='') then
      write(display,*)message_char1,message_char2,message_int
   else
      write(display,form)message_char1,message_char2,message_int
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aai

subroutine writelog_aii(destination,form,message_char1,message_int1,message_int2)
implicit none
   character(*),intent(in)    ::  form,message_char1
   character(*),intent(in)    ::  destination
   integer,intent(in)         ::  message_int1,message_int2
   character(256)             ::  display

   if (form=='') then
      write(display,*)message_char1,message_int1,message_int2
   else
      write(display,form)message_char1,message_int1,message_int2
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aii

subroutine writelog_aia(destination,form,message_char1b,message_intb,message_char2b)
   implicit none
   character(*),intent(in)    ::  form,message_char1b,message_char2b
   character(*),intent(in)    ::  destination
   integer,intent(in)         ::  message_intb
   character(256)             ::  display

   if (form=='') then
      write(display,*)message_char1b,message_intb,message_char2b
   else
      write(display,form)message_char1b,message_intb,message_char2b
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aia

subroutine writelog_aaai(destination,form,message_char1,message_char2,message_char3,message_int)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
   character(*),intent(in)    ::  destination
   integer,intent(in)         ::  message_int
   character(256)             ::  display

   if (form=='') then
      write(display,*)message_char1,message_char2,message_char3,message_int
   else
      write(display,form)message_char1,message_char2,message_char3,message_int
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aaai

subroutine writelog_aaia(destinationb,formb,message_char1b,message_char2b,message_int,message_char3b)
   implicit none
   character(*),intent(in)    ::  formb,message_char1b,message_char2b,message_char3b
   character(*),intent(in)    ::  destinationb
   integer,intent(in)         ::  message_int
   character(256)             ::  display

   if (formb=='') then
      write(display,*)message_char1b,message_char2b,message_int,message_char3b
   else
      write(display,formb)message_char1b,message_char2b,message_int,message_char3b
   endif

   if (xmaster) then 
      if (scan(destinationb,'s')>0) write(*,*)trim(display)
      if (scan(destinationb,'l')>0) write(logfileid,*)trim(display)
      if (scan(destinationb,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aaia

subroutine writelog_aiaa(destination,form,message_char1b,message_intb,message_char2b,message_char3b)
   implicit none
   character(*),intent(in)    ::  form,message_char1b,message_char2b,message_char3b
   character(*),intent(in)       ::  destination
   integer,intent(in)         ::  message_intb
   character(256)             ::  display

   if (form=='') then
      write(display,*)message_char1b,message_intb,message_char2b,message_char3b
   else
      write(display,form)message_char1b,message_intb,message_char2b,message_char3b
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aiaa

subroutine writelog_aiaaa(destination,form,message_char1b,message_intb,message_char2b,message_char3b,message_char4b)
   implicit none
   character(*),intent(in)    ::  form,message_char1b,message_char2b,message_char3b,message_char4b
   character(*),intent(in)       ::  destination
   integer,intent(in)         ::  message_intb
   character(256)             ::  display

   if (form=='') then
      write(display,*)message_char1b,message_intb,message_char2b,message_char3b,message_char4b
   else
      write(display,form)message_char1b,message_intb,message_char2b,message_char3b,message_char4b
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aiaaa


subroutine writelog_aiai(destination,form,message_char1,message_int1,message_char2,message_int2)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2
   character(*),intent(in)       ::  destination
   integer,intent(in)         ::  message_int1,message_int2
   character(256)             ::  display

   if (form=='') then
      write(display,*)message_char1,message_int1,message_char2,message_int2
   else
      write(display,form)message_char1,message_int1,message_char2,message_int2
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aiai

subroutine writelog_aiaia(destination,form,mc1,mi1,mc2,mi2,mc3)
   implicit none
   character(*),intent(in)    ::  form,mc1,mc2,mc3
   character(*),intent(in)       ::  destination
   integer,intent(in)         ::  mi1,mi2
   character(256)             ::  display

   if (form=='') then
      write(display,*)mc1,mi1,mc2,mi2,mc3
   else
      write(display,form)mc1,mi1,mc2,mi2,mc3
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aiaia

subroutine writelog_aaaiai(destination,form,message_char1,message_char2,message_char3,message_i1,message_char4,message_i2)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2,message_char3,message_char4
   character(*),intent(in)       ::  destination
   integer*4,intent(in)          ::  message_i1,message_i2
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)message_char1,message_char2,message_char3,message_i1,message_char4,message_i2
   else
      write(display,form)message_char1,message_char2,message_char3,message_i1,message_char4,message_i2
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aaaiai

subroutine writelog_aiafaf(destination,form,mc1,mi1,mc2,mf1,mc3,mf2)
   implicit none
   character(*),intent(in)    ::  form,mc1,mc2,mc3
   character(*),intent(in)       ::  destination
   integer*4,intent(in)          ::  mi1
   real*8,intent(in)             ::  mf1,mf2
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)mc1,mi1,mc2,mf1,mc3,mf2
   else
      write(display,form)mc1,mi1,mc2,mf1,mc3,mf2
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aiafaf

subroutine writelog_aiaiai(destination,form,message_char1,message_i1,message_char2,message_i2,message_char3,message_i3)
implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
   character(*),intent(in)    ::  destination
   integer*4,intent(in)       ::  message_i1,message_i2,message_i3
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)message_char1,message_i1,message_char2,message_i2,message_char3,message_i3
   else
      write(display,form)message_char1,message_i1,message_char2,message_i2,message_char3,message_i3
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aiaiai

subroutine writelog_aiaiaia(destination,form,message_char1,message_i1,message_char2,message_i2,message_char3,message_i3, & 
                            message_char4)
implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2,message_char3,message_char4
   character(*),intent(in)    ::  destination
   integer*4,intent(in)       ::  message_i1,message_i2,message_i3
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)message_char1,message_i1,message_char2,message_i2,message_char3,message_i3,message_char4
   else
      write(display,form)message_char1,message_i1,message_char2,message_i2,message_char3,message_i3,message_char4
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aiaiaia

subroutine writelog_aiaiaf(destination,form,mc1,mi1,mc2,mi2,mc3,mf1)
implicit none
   character(*),intent(in)    ::  form,mc1,mc2,mc3
   character(*),intent(in)    ::  destination
   integer*4,intent(in)       ::  mi1,mi2
   real*8,intent(in)          ::  mf1
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)mc1,mi1,mc2,mi2,mc3,mf1
   else
      write(display,form)mc1,mi1,mc2,mi2,mc3,mf1
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aiaiaf

subroutine writelog_aiaiafa(destination,form,mc1,mi1,mc2,mi2,mc3,mf1,mc4)
implicit none
   character(*),intent(in)    ::  form,mc1,mc2,mc3,mc4
   character(*),intent(in)    ::  destination
   integer*4,intent(in)       ::  mi1,mi2
   real*8,intent(in)          ::  mf1
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)mc1,mi1,mc2,mi2,mc3,mf1,mc4
   else
      write(display,form)mc1,mi1,mc2,mi2,mc3,mf1,mc4
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aiaiafa

subroutine writelog_iiiii(destination,form,mi1,mi2,mi3,mi4,mi5)
implicit none
   character(*),intent(in)    ::  form
   character(*),intent(in)    ::  destination
   integer*4,intent(in)       ::  mi1,mi2,mi3,mi4,mi5
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)mi1,mi2,mi3,mi4,mi5
   else
      write(display,form)mi1,mi2,mi3,mi4,mi5
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif

end subroutine writelog_iiiii

subroutine writelog_af(destination,form,message_char1,message_f1)
implicit none
   character(*),intent(in)    ::  form, message_char1
   character(*),intent(in)    ::  destination
   real*8,intent(in)          ::  message_f1
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)message_char1,message_f1
   else
      write(display,form)message_char1,message_f1
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_af

subroutine writelog_aaf(destination,form,message_char1,message_char2,message_f1)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2
   character(*),intent(in)       ::  destination
   real*8,intent(in)          ::  message_f1
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)message_char1,message_char2,message_f1
   else
      write(display,form)message_char1,message_char2,message_f1
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aaf

subroutine writelog_afa(destination,form,mc1,mf1,mc2)
   implicit none
   character(*),intent(in)    ::  form,mc1,mc2
   character(*),intent(in)       ::  destination
   real*8,intent(in)          ::  mf1
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)mc1,mf1,mc2
   else
      write(display,form)mc1,mf1,mc2
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_afa

subroutine writelog_afaf(destination,form,message_char1,message_f1,message_char2,message_f2)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2
   character(*),intent(in)       ::  destination
   real*8,intent(in)          ::  message_f1,message_f2
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)message_char1,message_f1,message_char2,message_f2
   else
      write(display,form)message_char1,message_f1,message_char2,message_f2
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_afaf


subroutine writelog_aaaf(destination,form,message_char1,message_char2,message_char3,message_f1)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
   character(*),intent(in)       ::  destination
   real*8,intent(in)          ::  message_f1
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)message_char1,message_char2,message_char3,message_f1
   else
      write(display,form)message_char1,message_char2,message_char3,message_f1
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aaaf

subroutine writelog_aafa(destination,form,message_char1b,message_char2b,message_f1b,message_char3b)
   implicit none
   character(*),intent(in)    ::  form,message_char1b,message_char2b,message_char3b
   character(*),intent(in)       ::  destination
   real*8,intent(in)          ::  message_f1b
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)message_char1b,message_char2b,message_f1b,message_char3b
   else
      write(display,form)message_char1b,message_char2b,message_f1b,message_char3b
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aafa

subroutine writelog_afaaa(destination,form,mc1a,mfa,mc2a,mc3a,mc4a)
   implicit none
   character(*),intent(in)    ::  form,mc1a,mc2a,mc3a,mc4a
   character(*),intent(in)       ::  destination
   real*8,intent(in)          ::  mfa
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)mc1a,mfa,mc2a,mc3a,mc4a
   else
      write(display,form)mc1a,mfa,mc2a,mc3a,mc4a
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_afaaa

subroutine writelog_aafaf(destination,form,message_char1,message_char2,message_f1,message_char3,message_f2)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2,message_char3
   character(*),intent(in)       ::  destination
   real*8,intent(in)          ::  message_f1,message_f2
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)message_char1,message_char2,message_f1,message_char3,message_f2
   else
      write(display,form)message_char1,message_char2,message_f1,message_char3,message_f2
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aafaf

subroutine writelog_aaafaf(destination,form,message_char1,message_char2,message_char3,message_f1,message_char4,message_f2)
   implicit none
   character(*),intent(in)    ::  form,message_char1,message_char2,message_char3,message_char4
   character(*),intent(in)       ::  destination
   real*8,intent(in)          ::  message_f1,message_f2
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)message_char1,message_char2,message_char3,message_f1,message_char4,message_f2
   else
      write(display,form)message_char1,message_char2,message_char3,message_f1,message_char4,message_f2
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_aaafaf

subroutine writelog_afafafaf(destination,form,mc1,mf1,mc2,mf2,mc3,mf3,mc4,mf4)
   implicit none
   character(*),intent(in)    ::  form,mc1,mc2,mc3,mc4
   character(*),intent(in)    ::  destination
   real*8,intent(in)          ::  mf1,mf2,mf3,mf4
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)mc1,mf1,mc2,mf2,mc3,mf3,mc4,mf4
   else
      write(display,form)mc1,mf1,mc2,mf2,mc3,mf3,mc4,mf4
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif
end subroutine writelog_afafafaf

subroutine writelog_illll(destination,form,mi1,ml1,ml2,ml3,ml4)
implicit none
   character(*),intent(in)    ::  form
   character(*),intent(in)    ::  destination
   integer*4,intent(in)       ::  mi1
   logical,intent(in)         ::  ml1,ml2,ml3,ml4
   character(256)             ::  display
 
   if (form=='') then
      write(display,*)mi1,ml1,ml2,ml3,ml4
   else
      write(display,form)mi1,ml1,ml2,ml3,ml4
   endif

   if (xmaster) then 
      if (scan(destination,'s')>0) write(*,*)trim(display)
      if (scan(destination,'l')>0) write(logfileid,*)trim(display)
      if (scan(destination,'e')>0) write(errorfileid,*)trim(display)   
   endif

end subroutine writelog_illll

end module logging_module
