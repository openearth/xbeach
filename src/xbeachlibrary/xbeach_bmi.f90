module xbeach_bmi
  use iso_c_binding
  use iso_c_utils
  use libxbeach_module
  use timestep_module
  use params

  implicit none

  ! initialize  -> loadmodel/initmodel
  ! perform single timestep -> update
  ! finalize -> finalize



  ! This is assumed.....
  integer(c_int), parameter :: MAXDIMS = 6


contains



  integer(c_int) function finalize() result(ierr) bind(C, name="finalize")
    !DEC$ ATTRIBUTES DLLEXPORT::finalize
    ierr = 0
    call logmsg(LEVEL_INFO, 'Finalize')
    ierr = final()
  end function finalize


  integer(c_int) function initialize(c_configfile) result(ierr) bind(C, name="initialize")
    !DEC$ ATTRIBUTES DLLEXPORT::initialize

    implicit none

    ! Variables
    character(kind=c_char), intent(in) :: c_configfile(*)
    character(len=strlen(c_configfile)) :: configfile

    ! Convert c string to fortran string
    configfile = char_array_to_string(c_configfile)

    ierr = init()

  end function initialize


  !> Performs a single timestep with the current model.
  integer(c_int) function update(dt) result(ierr) bind(C,name="update")
    !DEC$ ATTRIBUTES DLLEXPORT::update

    !< Custom timestep size, use -1 to use model default.
    real(c_double), value, intent(in) :: dt

    if (dt >= 0) then
       ierr = executestep(dt)
    else
       ierr = executestep()
    end if

    ! enable output to be able to use aggregated variable values (min,
    ! max, mean, var) through BMI. user can still control output
    ! through params.txt file
    ierr = outputext()

  end function update


  ! Void function is a subroutine
  subroutine get_var_type(c_var_name, c_type_name)  bind(C, name="get_var_type")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_type

    character(kind=c_char), intent(in) :: c_var_name(*)
    character(kind=c_char), intent(out) :: c_type_name(MAXSTRINGLEN)

    character(len=strlen(c_var_name)) :: var_name
    character(len=MAXSTRINGLEN) :: type_name

    character :: typecode
    integer :: index
    type(arraytype) :: array


    var_name = char_array_to_string(c_var_name)


    index =  chartoindex(var_name)
    if (index .eq. -1) return
    call indextos(s,index,array)
    typecode = array%type

    select case(typecode)
    case("r")
       type_name = "double"
    case("i")
       type_name = "int"
    case("c")
       type_name = "character"
    case default
       continue
    end select
    c_type_name = string_to_char_array(trim(type_name))

  end subroutine get_var_type

  subroutine get_var_rank(c_var_name, rank) bind(C, name="get_var_rank")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_rank

    character(kind=c_char), intent(in) :: c_var_name(*)
    integer(c_int), intent(out) :: rank

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    type(arraytype) :: array
    integer :: index

    ! Store the name
    var_name = char_array_to_string(c_var_name)

    index =  chartoindex(var_name)
    if (index .eq. -1) return
    call indextos(s,index,array)
    rank = array%rank
  end subroutine get_var_rank

  include 'get_var_shape.inc'
  include 'get_var.inc'
  include 'set_var.inc'
  subroutine set_current_time(xptr) bind(C, name="set_current_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: set_current_time

    use iso_c_binding, only: c_ptr, c_double, c_f_pointer

    type(c_ptr), value, intent(in) :: xptr
    real(c_double), pointer :: x_0d_double_ptr

    call c_f_pointer(xptr, x_0d_double_ptr)
    if (par%morfacopt == 1) then
       par%t = x_0d_double_ptr / max(par%morfac, 1.d0)
    else
       par%t = x_0d_double_ptr
    endif

  end subroutine set_current_time

  subroutine get_current_time(time) bind(C, name="get_current_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_current_time

    real(c_double) :: time

    if (par%morfacopt == 1) then
       time = par%t * max(par%morfac, 1.d0)
    else
       time = par%t
    endif
  end subroutine get_current_time

  subroutine get_start_time(time) bind(C, name="get_start_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_start_time

    real(c_double) :: time

    time = 0
  end subroutine get_start_time

  subroutine get_time_step(timestep) bind(C, name="get_time_step")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_time_step

    real(c_double) :: timestep
    integer :: ilim = 0
    integer :: jlim = 0
    real*8 :: dtref = 0.d0
    
    call compute_dt(s,par, tpar, it, ilim, jlim, dtref)
    timestep = par%dt
  end subroutine get_time_step

  subroutine get_end_time(time) bind(C, name="get_end_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_end_time

    real(c_double) :: time
    if (par%morfacopt == 1) then
       time = par%tstop * max(par%morfac, 1.d0)
    else
       time = par%tstop
    endif
  end subroutine get_end_time


end module xbeach_bmi
