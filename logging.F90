MODULE logging_module

use xmpi_module

implicit none

integer,save     :: logfileid
integer,save     :: errorfileid


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
   module procedure writelog_aaa
   module procedure writelog_aaaa
   module procedure writelog_aai
   module procedure writelog_aaai
   module procedure writelog_aaia
   module procedure writelog_aia
   module procedure writelog_aiaa
   module procedure writelog_aiai
   module procedure writelog_aaaiai
   module procedure writelog_aaf
   module procedure writelog_afaf
   module procedure writelog_aaaf
   module procedure writelog_aafa
   module procedure writelog_aafaf
   module procedure writelog_aaafaf
end interface writelog

CONTAINS

subroutine start_logfiles(error)
   implicit none
   integer    :: error
   integer    :: tryunit = 98
   logical    :: fileopen
   
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
end subroutine start_logfiles

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

end module logging_module
