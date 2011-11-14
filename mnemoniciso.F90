module mnemiso_module
use iso_c_binding
use mnemmodule


type, bind(c) :: b_arraytype

  character(kind=c_char) type         ! 'i' or 'r': integer or real*8
  character(kind=c_char) btype        ! 'b' or 'd' '2': broadcast or distribute or 
                         !                 distribute umean
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
function char_array_to_string(char_array, length)
   integer(c_int) :: length
    character(c_char) :: char_array(length)
    character(len=length) :: char_array_to_string
    integer :: i
    do i = 1, length
        char_array_to_string(i:i) = char_array(i)
    enddo
end function


end module mnemiso_module
