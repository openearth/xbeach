module general_fileio

interface check_file_length
   module procedure check_file_length_1D
   module procedure check_file_length_2D
   module procedure check_file_length_3D
!   module procedure check_file_length_auto
end interface check_file_length

contains

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
   open(fid,file=fname)
   read(fid,*,iostat=iost)(dat(i),i=1,d1)
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
   open(fid,file=fname)
   read(fid,*,iostat=iost)((dat(i,j),i=1,d1),j=1,d2)
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
   open(fid,file=fname)
   read(fid,*,iostat=iost)(((dat(i,j,k),i=1,d1),j=1,d2),k=1,d3)
   if (iost .ne. 0) then
      if (xmaster) then
         write(*,*)'Error processing file ''',trim(fname),'''.',' File may be too short or contains invalid values.', & 
				                                                 ' Terminating simulation'
         call halt_program
      endif
   endif
   close(fid)
end subroutine check_file_length_3D  

subroutine check_file_exist(filename)
use xmpi_module
   implicit none

   character(*)               :: filename
   logical                    :: file_exists

   inquire(file=filename,exist=file_exists)
   
   if (.not. file_exists) then
      if (xmaster) then
         write(*,*)'File ''',trim(filename),''' not found. Terminating simulation'
	     call halt_program
	  endif
   endif

end subroutine check_file_exist

integer function create_new_fid()
use xmpi_module
   integer    :: tryunit = 98
   logical    :: fileopen
   
   fileopen = .true.    
   if (xmaster) then
      do while (fileopen==.true.)
         inquire(tryunit,OPENED=fileopen)
	     if (fileopen==.true.) then
	        tryunit=tryunit-1
	     endif
	     if (tryunit<=10) then 
	       write(*,*)'Serious problem: not enough free unit ids to create new file'
	       call halt_program
	     endif	      
      enddo
   endif
   create_new_fid = tryunit   
end function create_new_fid

end module general_fileio
