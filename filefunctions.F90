module filefunctions

interface check_file_length
   module procedure check_file_length_1D
   module procedure check_file_length_2D
   module procedure check_file_length_3D
end interface check_file_length

contains

integer function create_new_fid()
use general_fileio
use xmpi_module
use logging_module

   integer    :: fileunit

   if (xmaster) then 
      fileunit = create_new_fid_generic()
   endif
#ifdef USEMPI
   call xmpi_bcast(fileunit)
#endif   
   if (fileunit==-1) then
      call writelog('les','','Serious problem: not enough free unit ids to create new file')
      call halt_program
   endif
   create_new_fid = fileunit   
end function create_new_fid

subroutine check_file_exist(filename)
use general_fileio
use xmpi_module
use logging_module

   implicit none

   character(*)               :: filename
   integer                    :: error
   logical                    :: file_exists


   if (xmaster) call check_file_exist_generic(filename,error)
#ifdef USEMPI
   call xmpi_bcast(error)
#endif 

   if (error==1) then
      call writelog('sle','','File ''',trim(filename),''' not found. Terminating simulation')
      call halt_program
   endif
end subroutine check_file_exist



subroutine check_file_length_1D(fname,d1)
use xmpi_module
  
   implicit none
   character(*)                   ::  fname
   integer, intent(in)            ::  d1
   integer                        ::  fid,iost
   integer                        ::  i
   real,dimension(:),allocatable  ::  dat
  
   allocate(dat(d1))
   
   fid = create_new_fid()
   if (xmaster) then
      open(fid,file=fname)
      read(fid,*,iostat=iost)(dat(i),i=1,d1)
   endif
#ifdef USEMPI
   call xmpi_bcast(iost)
#endif 
   if (iost .ne. 0) then
      if (xmaster) then
         write(*,*)'Error processing file ''',trim(fname),'''.',' File may be too short or contains invalid values.', & 
				                                                 ' Terminating simulation' 
         call halt_program
      endif
   endif
   close(fid)
end subroutine check_file_length_1D

subroutine check_file_length_2D(fname,d1,d2)
use xmpi_module
  
   implicit none
   character(*)                     :: fname
   integer, intent(in)              :: d1,d2
   integer                          :: fid,iost
   integer                          :: i,j
   real,dimension(:,:),allocatable  :: dat

   allocate(dat(d1,d2))
   
   fid = create_new_fid()
   if (xmaster) then 
      open(fid,file=fname)
      read(fid,*,iostat=iost)((dat(i,j),i=1,d1),j=1,d2)
   endif
#ifdef USEMPI
   call xmpi_bcast(iost)
#endif 
   if (iost .ne. 0) then
      if (xmaster) then
         write(*,*)'Error processing file ''',trim(fname),'''.',' File may be too short or contains invalid values.', & 
				                                                 ' Terminating simulation'
         call halt_program
      endif
   endif
   close(fid)
end subroutine check_file_length_2D  

subroutine check_file_length_3D(fname,d1,d2,d3)
use xmpi_module
  
   implicit none
   character(*)                       ::  fname
   integer, intent(in)                ::  d1,d2,d3
   integer                            ::  fid,iost
   integer                            ::  i,j,k
   real,dimension(:,:,:),allocatable  ::  dat

   allocate(dat(d1,d2,d3))
   
   fid = create_new_fid()
   if (xmaster) then 
      open(fid,file=fname)
      read(fid,*,iostat=iost)(((dat(i,j,k),i=1,d1),j=1,d2),k=1,d3)
   endif
#ifdef USEMPI
   call xmpi_bcast(iost)
#endif 
   if (iost .ne. 0) then
      if (xmaster) then
         write(*,*)'Error processing file ''',trim(fname),'''.',' File may be too short or contains invalid values.', & 
				                                                 ' Terminating simulation'
         call halt_program
      endif
   endif
   close(fid)
end subroutine check_file_length_3D  

subroutine checkbcfilelength(tstop,instat,filename,filetype)
use logging_module

IMPLICIT NONE
real*8, intent(in) :: tstop
integer, intent(in):: instat
character*80      :: filename,dummy
character*8       :: testc
character*1       :: ch
integer           :: i,ier=0,nlines,filetype,fid
real*8            :: t,dt,total,d1,d2,d3,d4,d5

ier = 0
fid = create_new_fid()
open(fid,file=filename)
i=0
do while (ier==0)
   read(fid,'(a)',iostat=ier)ch
   if (ier==0)i=i+1
enddo
nlines=i
rewind(fid)    
	
if (instat==4 .or. instat==5 .or. instat==6) then 
	read(fid,*)testc
    if (testc=='FILELIST') then
		filetype = 1
		nlines=nlines-1
	else
		filetype = 0
	endif
elseif (instat==40 .or. instat==41) then
    filetype = 2
endif
    
total=0.d0
i=0
select case (filetype)
    case(0)
	   total=2.d0*tstop
	case(1)
	   do while (total<tstop .and. i<nlines)
	      read(fid,*)t,dt,dummy
		  total=total+t
		  i=i+1
	   enddo
	case(2)
	   do while (total<tstop .and. i<nlines)
		  read(fid,*)d1,d2,d3,d4,d5,t,dt
		  total=total+t
		  i=i+1
	   enddo
end select
close(fid)

if (total<tstop) then
	call writelog('sle',' ','Error: Wave boundary condition time series too short in ',trim(filename))
    call writelog('sle','(a,f0.2,a,f0.2)',' Total wave condition time series is ',total, ' but simulation length is ',tstop)
	call writelog('sle',' ','Stopping calculation')
	call halt_program
endif

end subroutine checkbcfilelength

end module filefunctions


module general_fileio

contains

subroutine check_file_exist_generic(filename,error)
   implicit none

   character(*)               :: filename
   integer                    :: error
   logical                    :: file_exists

   inquire(file=filename,exist=file_exists)
   
   error = 0

   if (.not. file_exists) then
	      error = 1
   endif

end subroutine check_file_exist_generic


integer function create_new_fid_generic()
   integer    :: tryunit = 98
   logical    :: fileopen
   
   fileopen = .true.    
   do while (fileopen)
         inquire(tryunit,OPENED=fileopen)
	     if (fileopen) then
	        tryunit=tryunit-1
	     endif
	     if (tryunit<=10) then 
		   tryunit = -1
		   fileopen = .false.
	     endif	      
   enddo
   create_new_fid_generic = tryunit   
end function create_new_fid_generic

end module general_fileio
