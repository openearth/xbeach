module mnemiso_module
use iso_c_binding
use mnemmodule


type, bind(c) :: b_arraytype

  character(kind=c_char) type         ! 'i' or 'r': integer or real*8
  character(kind=c_char) btype        ! 'b' or 'd':
  integer(c_int) rank           ! 0,1,2,3,4
  character(kind=c_char, len=maxnamelen) :: name     ! 'v','ve', .....
  character(kind=c_char, len=20) :: units     ! m, following udunits convention
  character(kind=c_char, len=1024) :: description
  character(kind=c_char, len=20), dimension(maxrank) :: dimensions ! the dimensions of the variable, for example (s%nx, s%ny)

  type (c_ptr) :: array
  
end type b_arraytype

type, bind(c) :: carraytype

  integer(c_int) rank           ! 0,1,2,3,4
  character(kind=c_char) type         ! 'i' or 'r': integer or real*8
  character(kind=c_char) btype        ! 'b' or 'd' '2': broadcast or distribute or 
  type(c_ptr) :: array
  
end type carraytype

contains 
type(carraytype) function arrayf2c(farray)
  type(arraytype), intent(in) :: farray
  real(c_double), target :: a(3)
  a = (/ 1.0d0,  2.0d0, 3.0d0 /)
  arrayf2c%type = farray%type
  arrayf2c%btype = farray%btype
  arrayf2c%rank = farray%rank
  ! if (farray%type == 'd' .and. farray%rank == 2) then
  !    call c_f_pointer(arrayf2c%array, farray%r2, shape=shape(farray%r2))
  ! else 
  !    arrayf2c%array =C_NULL_PTR
  ! end if
  arrayf2c%array = c_loc(a)
end function arrayf2c



! Utility functions
integer(c_int) function stringlength(char_array)
    character(c_char), intent(in) :: char_array(:)
    integer :: inull, i 
    stringlength = 0
    do i = 1, size(char_array)
        if (char_array(i) .eq. C_NULL_CHAR) then
            stringlength = i
        end if
    end do
    stringlength = size(char_array)
    
end function

function char_array_to_string(char_array, length)
   integer(c_int) :: length
    character(c_char) :: char_array(length)
    character(len=length) :: char_array_to_string
    integer :: i
    do i = 1, length
        char_array_to_string(i:i) = char_array(i)
    enddo
end function
function string_to_char_array(string, length)
  character(len=length) :: string
  character(kind=c_char,len=1) :: string_to_char_array(length+1)
  integer(c_int) :: length
  integer :: i
  do i = 1, length
     string_to_char_array(i) = string(i:i)
  enddo
  string_to_char_array(length+1) = C_NULL_CHAR
end function


!   FUNCTION C_F_STRING(CPTR) RESULT(FPTR)
!      ! Convert a null-terminated C string into a Fortran character array pointer
!      TYPE(C_PTR), INTENT(IN) :: CPTR ! The C address
!      CHARACTER(KIND=C_CHAR), DIMENSION(:), POINTER :: FPTR
!      
!      INTERFACE ! strlen is a standard C function from <string.h>
!         ! int strlen(char *string)
!         FUNCTION strlen(string) RESULT(len) BIND(C,NAME="strlen")
!            USE ISO_C_BINDING
!            TYPE(C_PTR), VALUE :: string ! A C pointer
!         END FUNCTION
!      END INTERFACE   
!      
!      IF(C_ASSOCIATED(CPTR)) THEN
!         CALL C_F_POINTER(FPTR=FPTR, CPTR=CPTR, SHAPE=[strlen(CPTR)])
!      ELSE
!         ! To avoid segfaults, associate FPTR with a dummy target:
!         FPTR=>dummy_string
!      END IF
!            
!   END FUNCTION
end module mnemiso_module
