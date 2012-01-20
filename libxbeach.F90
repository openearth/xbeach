module libxbeach_module
  use iso_c_binding
  use mnemiso_module
  use mnemmodule
  use getkey_module
  use params
  use spaceparams
  use xmpi_module
  use initialize
  use boundaryconditions
  use drifter_module
  use flow_timestep_module
  use morphevolution
  use readtide_module
  use readwind_module
  use wave_timestep_module
  use timestep_module
  use readkey_module
  use groundwaterflow
  use logging_module
  use means_module
  use output_module

  implicit none

  type(parameters), save                              :: par
  type(timepars), save                                     :: tpar
  type(spacepars), pointer                            :: s
  type(spacepars), target, save                       :: sglobal

  integer                                             :: n,it,error
  real*8                                              :: tbegin

#ifdef USEMPI
  type(spacepars), target, save                       :: slocal
  real*8                                              :: t0,t01
#endif

  interface getnparameter
     module procedure getnparameter_fortran
  end interface getnparameter

  interface getparametername
     module procedure getparametername_fortran
  end interface getparametername
  
  interface getparametertype
     module procedure getparametertype_fortran
  end interface getparametertype

  interface getdoubleparameter
      module procedure getdoubleparameter_fortran
  end interface getdoubleparameter

  interface setdoubleparameter
      module procedure setdoubleparameter_fortran
  end interface setdoubleparameter

  interface getintparameter
      module procedure getintparameter_fortran
  end interface getintparameter
      
  interface getarraytype
     module procedure getarraytype_fortran
  end interface getarraytype

  interface get2ddoublearray
      module procedure get2ddoublearray_fortran
  end interface get2ddoublearray

  interface get2dintarray
      module procedure get2dintarray_fortran
   end interface get2dintarray

  interface set2ddoublearray
      module procedure set2ddoublearray_fortran
  end interface set2ddoublearray

  !startinit

  !-----------------------------------------------------------------------------!
  ! Initialize program                                                          !
  !-----------------------------------------------------------------------------!


contains
  integer(c_int) function init() bind(C, name="init")
    !DEC$ ATTRIBUTES DLLEXPORT::Init

    error   = 0

    ! setup of MPI 
#ifdef USEMPI
    s=>slocal
    call xmpi_initialize
    t0 = MPI_Wtime()
#endif

    ! create log files
    call start_logfiles(error)

    ! set starting time and date
    call cpu_time(tbegin)

    ! show statup message
    call writelog_startup()

    !-----------------------------------------------------------------------------!
    ! Initialize simulation                                                       !
    !-----------------------------------------------------------------------------!

    ! initialize time counter
    it      = 0

    ! read input from params.txt
    call all_input(par)

    ! allocate space scalars
    call space_alloc_scalars(sglobal)
    s => sglobal

    ! read grid and bathymetry
    call grid_bathy(s,par)

    ! distribute grid over processors
#ifdef USEMPI
    call xmpi_determine_processor_grid(s%nx,s%ny,par%mpiboundary,error)
    call writelog_mpi(par%mpiboundary,error)
#endif

    ! initialize timestep
    call timestep_init(par, tpar)

    if (xmaster) then

       call writelog('ls','','Initializing .....')

       ! initialize physics
       call readtide           (s,par)
       call readwind           (s,par)

       call flow_init          (s,par)
       call discharge_init     (s,par)
       call drifter_init       (s,par)
       call wave_init          (s,par)
       call gw_init            (s,par)
       call sed_init           (s,par)

    endif

#ifdef USEMPI
    call distribute_par(par)
    s => slocal
    call space_distribute_space (sglobal,s,par     )
#endif

    ! initialize output
    call means_init             (sglobal,s,par     )
    call output_init            (sglobal,s,par,tpar)


    ! store first timestep
    call output(sglobal,s,par,tpar)
    init = 0
  end function init

  integer(c_int) function outputext() bind(C, name="outputext")
    !DEC$ ATTRIBUTES DLLEXPORT::outputext
    ! store first timestep
    call output(sglobal,s,par,tpar)
    outputext = 0
  end function outputext
  !-----------------------------------------------------------------------------!
  ! Start simulation                                                            !
  !-----------------------------------------------------------------------------!

  integer(c_int) function executestep() bind(C, name="executestep")
    !DEC$ ATTRIBUTES DLLEXPORT::executestep
     

    ! This is now called in output, but timestep depends on it...
    ! call outputtimes_update(par, tpar)
    !n = 0
    executestep = -1
    !do while (par%t<par%tstop)
    ! determine timestep
    call timestep(s,par,tpar,it,ierr=error)
    if (error==1) call output_error(s,sglobal,par,tpar)

    ! boundary conditions
    call wave_bc        (sglobal,s,par)
    if (par%gwflow==1)       call gw_bc          (s,par)
    if (par%flow+par%nonh>0) call flow_bc        (s,par)

#ifdef USEMPI
    if (it==0) t01 = MPI_Wtime()
#endif

    ! compute timestep
    if (par%swave==1)        call wave           (s,par)
    if (par%gwflow==1)       call gwflow         (s,par)
    if (par%flow+par%nonh>0) call flow           (s,par)
    if (par%ndrifter>0)      call drifter        (s,par)
    if (par%sedtrans==1)     call transus        (s,par)
    if (par%morphology==1)   call bed_update     (s,par)

    ! output
    ! This has a nasty dependency on outputtimes_update....
    ! call output(sglobal,s,par,tpar)

    ! n = n + 1
    executestep = 0
    ! enddo
  end function executestep

  ! No C for this one, no chars to mess up...
  integer(c_int) function getnparameter_fortran(n) bind(C, name="getnparameter")
    integer(c_int), intent(inout) :: n
    character(len=slen), allocatable :: keys(:)
    getnparameter_fortran = -1
    call getkeys(par, keys)
    n = size(keys,1)
    getnparameter_fortran = 0
  end function getnparameter_fortran

  integer(c_int) function getparametertype_fortran(name, typecode)
    character(kind=c_char,len=*),intent(in) :: name
    character(kind=c_char, len=1), intent(out) :: typecode
    
    integer(c_int) :: length
    character(1), dimension(len(name)) :: cname 
    integer :: i
    length  = len(name)
    cname = string_to_char_array(name,length)
    
    getparametertype_fortran = getparametertype_c(cname, typecode, length)
  end function getparametertype_fortran

  integer(c_int) function getparametertype_c(name, typecode, length) bind(C,name="getparametertype")
    character(kind=c_char,len=1),intent(in) :: name(length)
    character(kind=c_char, len=1), intent(out) :: typecode
    integer(c_int) :: length
    
    character(len=length) :: key
    integer :: index
    key = char_array_to_string(name, length)
    call getkey_indextype(par, key, index, typecode)
  end function getparametertype_c

  integer(c_int) function getparametername_fortran(index, name)
    integer(c_int), intent(in) :: index
    character(kind=c_char, len=*), intent(out) :: name
    character(kind=c_char, len=1), pointer :: cname(:)
    integer :: length

    getparametername_fortran = getparametername_c(index, cname, length)
    name = char_array_to_string(cname, length)
  end function getparametername_fortran


  integer(c_int) function getparametername_c(index, name, length) bind(C,name="getparametername")
    integer(c_int), intent(in) :: index
    integer(c_int), intent(out) :: length
    character(kind=c_char, len=1), intent(out) :: name(slen)

    integer :: i,j
    character(len=slen), allocatable :: keys(:)
    character(kind=c_char,len=slen) :: key
    getparametername_c = -1
    ! These are the keys in fortran format.
    call getkeys(par, keys)
    ! We need to conver them to C format (char1's)
    key = keys(index)
    length = len(trim(key))
    name = string_to_char_array(key, length)
    getparametername_c = 0
  end function getparametername_c

  integer(c_int) function getdoubleparameter_fortran(name,value) 
    USE iso_c_binding
    ! use inout otherwise things break
    real(c_double), intent(inout) :: value

    ! String
    character(kind=c_char,len=*),intent(in) :: name

    ! Transform name to a fortran character... 
    character(1), dimension(len(name)) :: cname 
    integer :: i
    do i = 1,len(name)
        cname(i) = name(i:i) 
    enddo
    getdoubleparameter_fortran = getdoubleparameter_c(cname,value,len(name))
  end function getdoubleparameter_fortran


  integer(c_int) function getdoubleparameter_c(name,value, length) bind(C,name="getdoubleparameter")
    !DEC$ ATTRIBUTES DLLEXPORT::getdoubleparameter_c
    USE iso_c_binding
    use getkey_module
    ! use inout otherwise things break
    real(c_double), intent(inout) :: value
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    ! Transform name to a fortran character... 
    type(parameter) :: myparam
    character(length) :: myname 
    ! Return -1 for invalid parameters
    getdoubleparameter_c = -1
    myname = char_array_to_string(name, length)
    getdoubleparameter_c =  getkey(par, myname, myparam)
    if (getdoubleparameter_c .eq. -1) return
    value = myparam%r0
    getdoubleparameter_c = 0
  end function getdoubleparameter_c


  integer(c_int) function setdoubleparameter_fortran(name,value) 
    USE iso_c_binding
    ! use inout otherwise things break
    real(c_double), intent(inout) :: value

    ! String
    character(kind=c_char,len=*),intent(in) :: name

    ! Transform name to a fortran character... 
    character(1), dimension(len(name)) :: myname 
    integer :: i
    do i = 1,len(name)
        myname(i) = name(i:i) 
    enddo
    setdoubleparameter_fortran = setdoubleparameter_c(myname,value,len(name))
  end function setdoubleparameter_fortran

  integer(c_int) function setdoubleparameter_c(name,value, length) bind(C,name="setdoubleparameter")
    !DEC$ ATTRIBUTES DLLEXPORT::setdoubleparameter_c
    USE iso_c_binding
    ! use inout otherwise things break
    real(c_double), intent(in) :: value
    ! and we need the string length ....
    integer(c_int), value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    ! Transform name to a fortran character... 
    character(length) :: myname 
    myname = char_array_to_string(name, length)
    select case (myname)
    case ('t')
       par%t = value
    case ('tstop')
       par%tstop = value
    case ('tnext')
       tpar%tnext = value
    case default
       setdoubleparameter_c = -1
       return
    end select
    setdoubleparameter_c = 0
  end function setdoubleparameter_c

  integer(c_int) function getintparameter_fortran(name,value) 
    USE iso_c_binding
    ! use inout otherwise things break
    integer(c_int), intent(inout) :: value

    ! String
    character(kind=c_char,len=*),intent(in) :: name

    ! Transform name to a fortran character... 
    character(1), dimension(len(name)) :: myname 
    integer :: i
    do i = 1,len(name)
        myname(i) = name(i:i) 
    enddo
    getintparameter_fortran = getintparameter_c(myname,value,len(name))
  end function getintparameter_fortran

  integer(c_int) function getintparameter_c(name,value, length) bind(C,name="getintparameter")
    !DEC$ ATTRIBUTES DLLEXPORT::getintparameter_c

    USE iso_c_binding
    use getkey_module
    ! use inout otherwise things break
    integer(c_int), intent(inout) :: value
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    ! Transform name to a fortran character... 
    character(length) :: myname 
    type(parameter) :: myparam
    myname = char_array_to_string(name, length)
    ! Lookup the parameter by name
    getintparameter_c = getkey(par, myname, myparam)
    if (getintparameter_c .eq. -1) return
    value = myparam%i0
    getintparameter_c = 0
  end function getintparameter_c

  integer(c_int) function getcharparameter_fortran(name,value) 
    USE iso_c_binding

    ! String
    character(kind=c_char,len=*),intent(in) :: name
    character(kind=c_char, len=*), intent(out) :: value

    ! Transform name to a fortran character... 
    character(1), dimension(len(name)) :: cname 
    character(kind=c_char,len=1), pointer :: cvalue(:)
    integer(c_int) :: valuelength
    cname = string_to_char_array(name, len(name))
    getcharparameter_fortran = getcharparameter_c(cname,cvalue,len(name),valuelength)
  end function getcharparameter_fortran


  integer(c_int) function getcharparameter_c(name,value, namelength, valuelength) bind(C,name="getcharparameter")
    !DEC$ ATTRIBUTES DLLEXPORT::getcharparameter_c

    USE iso_c_binding
    use getkey_module
    ! String
    character(kind=c_char),intent(in)         :: name(namelength)
    character(kind=c_char,len=1), intent(out) :: value(slen)
    integer(c_int), intent(in)         :: namelength
    integer(c_int), intent(out)        :: valuelength

    ! Transform name to a fortran character... 
    character(namelength) :: fname 
    type(parameter) :: myparam
    fname = char_array_to_string(name, namelength)
    ! Lookup the parameter by name
    getcharparameter_c = getkey(par, fname, myparam)
    if (getcharparameter_c .eq. -1) return
    valuelength = len(trim(myparam%c0))
    value = string_to_char_array(trim(myparam%c0), valuelength)
    getcharparameter_c = 0
  end function getcharparameter_c

  integer(c_int) function getarraytype_fortran(name, typecode) 
    !DEC$ ATTRIBUTES DLLEXPORT::getarray_type
    character(kind=c_char,len=*),intent(in) :: name
    character(kind=c_char, len=1), intent(out) :: typecode

    ! and we need the string length ....
    integer(c_int) :: length
    character(1), dimension(len(name)) :: cname 
    integer :: i
    length  = len(name)
    cname = string_to_char_array(name,length)
    getarraytype_fortran = getarraytype_c(cname, typecode, length)
  end function getarraytype_fortran

  integer(c_int) function getarraytype_c(name, typecode, length) bind(C,name="getarraytype")
    character(kind=c_char,len=1),intent(in) :: name(length)
    character(kind=c_char, len=1), intent(out) :: typecode
    integer(c_int) :: length
    
    character(len=length) :: key
    integer :: index
    type(arraytype) :: array
    key = char_array_to_string(name, length)
    getarraytype_c = -1
    index =  chartoindex(key)
    if (index .eq. -1) return
    call indextos(s,index,array)
    typecode = array%type
  end function getarraytype_c

  integer(c_int) function getarray(name, x, length) bind(C, name="getarray")
    !DEC$ ATTRIBUTES DLLEXPORT::getarray

    ! use inout otherwise things break
    type(carraytype), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(inout) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    getarray = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    x = arrayf2c(array)
    getarray = 0
  end function getarray

  integer(c_int) function get1ddoublearray(name, x, length) bind(C, name="get1ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::get1ddoublearray

    ! use inout otherwise things break
    type (c_ptr), value :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(inout) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), target, allocatable, dimension(:)  :: r1

    get1ddoublearray = -1

    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    allocate(r1(size(array%r1,1)))
    r1 = array%r1
    ! array%r1 => r1
    x = c_loc(r1)
    get1ddoublearray = 0
  end function get1ddoublearray

  integer(c_int) function get2ddoublearray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    real(c_double), intent(inout) :: x(:,:)

    type(arraytype) :: array
    integer :: i, index

    get2ddoublearray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    x = array%r2
    get2ddoublearray_fortran = 0
  end function get2ddoublearray_fortran


  integer(c_int) function get2ddoublearray_c(name, x, length) bind(C, name="get2ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::get2ddoublearray_c

    ! use inout otherwise things break
    type(c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), target, allocatable, dimension(:,:)  :: r2

    get2ddoublearray_c = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    allocate(r2(size(array%r2,1), size(array%r2,2)))
    r2(:,:) = array%r2(:,:)
    ! array%r2 => r2
    x = c_loc(r2)
    get2ddoublearray_c = 0
  end function get2ddoublearray_c


  integer(c_int) function get2dintarray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    integer(c_int), intent(inout) :: x(:,:)

    type(arraytype) :: array
    integer :: i, index

    get2dintarray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    x = array%i2
    get2dintarray_fortran = 0
  end function get2dintarray_fortran

  integer(c_int) function get2dintarray_c(name, x, length) bind(C, name="get2dintarray")
    !DEC$ ATTRIBUTES DLLEXPORT::get2dintarray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    integer(c_int), target, allocatable, dimension(:,:)  :: i2

    get2dintarray_c = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    allocate(i2(size(array%i2,1), size(array%i2,2)))
    i2(:,:) = array%i2(:,:)
    ! array%r2 => r2
    x = c_loc(i2)
    get2dintarray_c = 0
  end function get2dintarray_c



  integer(c_int) function set2ddoublearray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    real(c_double), intent(inout) :: x(:,:)

    integer :: i, index
    type(arraytype) :: array

    set2ddoublearray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    array%r2 = x
    set2ddoublearray_fortran = 0
  end function set2ddoublearray_fortran

  integer(c_int) function set2ddoublearray_c(name, x, length) bind(C, name="set2ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::set2ddoublearray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), pointer, dimension(:,:)  :: r2

    set2ddoublearray_c = -1

    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    ! Transform the c pointer into a fortran pointer
    call c_f_pointer(x, r2, shape(array%r2))
    ! Copy the values, or the pointer... not sure.
    array%r2 = r2
    
    set2ddoublearray_c = 0
  end function set2ddoublearray_c


  integer(c_int) function finalize() bind(C, name="finalize")
    !DEC$ ATTRIBUTES DLLEXPORT::finalize


    !-----------------------------------------------------------------------------!
    ! Finalize simulation                                                         !
    !-----------------------------------------------------------------------------!

#ifdef USEMPI
    call writelog_finalize(tbegin,n,par%t,par%nx,par%ny,t0,t01)
    call xmpi_finalize
#else
    call writelog_finalize(tbegin,n,par%t,par%nx,par%ny)
#endif
    finalize = 0
  end function finalize
end module libxbeach_module
