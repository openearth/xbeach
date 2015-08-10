<%
def element0(var):
    """find the first element array"""
    parts = ["lbound(s%%%s, %d)" % (var["name"], i+1) for i in range(var["rank"])]
    if var["rank"] > 0:
        txt = "(" + ",".join(parts) + ")"
    else:
        txt = ""
    return txt
%>
subroutine get_var(c_var_name, x) bind(C, name="get_var")
  !DEC$ ATTRIBUTES DLLEXPORT :: get_var

  ! Return a pointer to the variable

  character(kind=c_char), intent(in) :: c_var_name(*)
  type(c_ptr), intent(inout) :: x

  character(len=strlen(c_var_name)) :: var_name
  ! Store the name
  ! The real array
  integer :: index
  type(arraytype) :: array

  var_name = char_array_to_string(c_var_name)

  index =  chartoindex(var_name)
  if (index .eq. -1) return
  call indextos(s,index,array)

  select case(var_name)
%for variable in variables:
   case ('${variable["name"]}')
      ! if you want to reference a pointer variable
      ! you have to point to the first element
      x = c_loc(s%${variable["name"]}${element0(variable)})
%endfor
  end select

end subroutine get_var
