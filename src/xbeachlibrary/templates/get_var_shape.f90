subroutine get_var_shape(c_var_name, shape) bind(C, name="get_var_shape")
  !DEC$ ATTRIBUTES DLLEXPORT :: get_var_shape

  character(kind=c_char), intent(in) :: c_var_name(*)
  integer(c_int), intent(inout) :: shape(MAXDIMS)

  character(len=strlen(c_var_name)) :: var_name

  var_name = char_array_to_string(c_var_name)
  shape = (/0, 0, 0, 0, 0, 0/)

  select case(var_name)

%for variable in variables:
    case("${variable['name']}")
%if variable['rank'] == 0:
       shape(1) = 0
%else:
%for dim in range(variable['rank']):
       ! return in c memory order
%if "shape" in variable:
       shape(${variable['rank'] - dim}) = ${variable["shape"][dim]}
%else:
       shape(${variable['rank'] - dim}) = size(${variable['name']}, ${dim+1})
%endif
%endfor
%endif
%endfor
  end select

end subroutine get_var_shape


