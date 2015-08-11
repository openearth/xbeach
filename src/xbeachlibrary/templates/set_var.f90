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

  subroutine set_var(c_var_name, xptr) bind(C, name="set_var")
    !DEC$ ATTRIBUTES DLLEXPORT :: set_var
    ! Return a pointer to the variable
    use iso_c_binding, only: c_double, c_char, c_loc, c_f_pointer

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), value, intent(in) :: xptr

    real(c_double), pointer :: x_0d_double_ptr
    real(c_double), pointer :: x_1d_double_ptr(:)
    real(c_double), pointer :: x_2d_double_ptr(:,:)
    real(c_double), pointer :: x_3d_double_ptr(:,:,:)
    real(c_double), pointer :: x_4d_double_ptr(:,:,:,:)
    integer(c_int), pointer :: x_0d_int_ptr
    integer(c_int), pointer :: x_1d_int_ptr(:)
    integer(c_int), pointer :: x_2d_int_ptr(:,:)
    integer(c_int), pointer :: x_3d_int_ptr(:,:,:)
    real(c_float), pointer  :: x_1d_float_ptr(:)
    real(c_float), pointer  :: x_2d_float_ptr(:,:)
    real(c_float), pointer  :: x_3d_float_ptr(:,:,:)
    ! The fortran name of the attribute name
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
%for var in variables:
    case("${var['name']}")
%if var['rank'] > 0:
       call c_f_pointer(xptr, x_${var['rank']}d_${var['type']}_ptr, shape(s%${var["name"]}))
%else:
       call c_f_pointer(xptr, x_${var['rank']}d_${var['type']}_ptr) ! skip shape in case of rank 0
%endif
       s%${var["name"]}${dimstr(":"*var['rank'])} = x_${var['rank']}d_${var['type']}_ptr
%endfor
    case default
       write(*,*) 'Setting unknown variable', var_name
    end select

  end subroutine set_var
