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
      do while (fileopen)
         inquire(tryunit,OPENED=fileopen)
	     if (fileopen) then
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

function s_eqi ( s1, s2 )
!
!*******************************************************************************
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!
!  Examples:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  implicit none
!
  character c1
  character c2
  integer i
  integer len1
  integer len2
  integer lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2
!
  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end function s_eqi

function ch_eqi ( c1, c2 )
!
!*******************************************************************************
!
!! CH_EQI is a case insensitive comparison of two characters for equality.
!
!
!  Examples:
!
!    CH_EQI ( 'A', 'a' ) is .TRUE.
!
!  Modified:
!
!    28 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character C1, C2, the characters to compare.
!
!    Output, logical CH_EQI, the result of the comparison.
!
  implicit none
!
  logical ch_eqi
  character c1
  character c1_cap
  character c2
  character c2_cap
!
  c1_cap = c1
  c2_cap = c2

  call ch_cap ( c1_cap )
  call ch_cap ( c2_cap )

  if ( c1_cap == c2_cap ) then
    ch_eqi = .true.
  else
    ch_eqi = .false.
  end if

  return
end function ch_eqi


subroutine ch_cap ( c )
!
!*******************************************************************************
!
!! CH_CAP capitalizes a single character.
!
!
!  Modified:
!
!    30 November 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  implicit none
!
  character c
  integer itemp
!
  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end subroutine ch_cap

end module general_fileio
