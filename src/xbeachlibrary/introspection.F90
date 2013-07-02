module introspection_module
  use iso_c_binding
  use libxbeach_module

  implicit none

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

  interface getnarray
     module procedure getnarray_fortran
  end interface getnarray

  interface getarraytype
     module procedure getarraytype_fortran
  end interface getarraytype

  interface getarrayrank
     module procedure getarrayrank_fortran
  end interface getarrayrank

  interface getarrayname
     module procedure getarrayname_fortran
  end interface getarrayname

  interface get0ddoublearray
     module procedure get0ddoublearray_fortran
  end interface get0ddoublearray

  interface get1ddoublearray
     module procedure get1ddoublearray_fortran
  end interface get1ddoublearray

  interface get2ddoublearray
     module procedure get2ddoublearray_fortran
  end interface get2ddoublearray

  interface get3ddoublearray
     module procedure get3ddoublearray_fortran
  end interface get3ddoublearray

  interface get0dintarray
     module procedure get0dintarray_fortran
  end interface get0dintarray

  interface get1dintarray
     module procedure get1dintarray_fortran
  end interface get1dintarray

  interface get2dintarray
     module procedure get2dintarray_fortran
  end interface get2dintarray

  interface set0dintarray
     module procedure set0dintarray_fortran
  end interface set0dintarray

  interface set1dintarray
     module procedure set1dintarray_fortran
  end interface set1dintarray

  interface set0ddoublearray
     module procedure set0ddoublearray_fortran
  end interface set0ddoublearray

  interface set1ddoublearray
     module procedure set1ddoublearray_fortran
  end interface set1ddoublearray

  interface set2ddoublearray
     module procedure set2ddoublearray_fortran
  end interface set2ddoublearray

  interface set3ddoublearray
     module procedure set3ddoublearray_fortran
  end interface set3ddoublearray

  interface set4ddoublearray
     module procedure set4ddoublearray_fortran
  end interface set4ddoublearray
contains
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

  integer(c_int) function getnarray_fortran(n) bind(C, name="getnarray")
    integer(c_int), intent(inout) :: n
    character(len=slen), allocatable :: keys(:)
    getnarray_fortran = -1
    n = numvars
    getnarray_fortran = 0
  end function getnarray_fortran


  integer(c_int) function getarrayname_fortran(index, name)
    integer(c_int), intent(in) :: index
    character(kind=c_char, len=*), intent(out) :: name
    character(kind=c_char, len=1), pointer :: cname(:)
    integer :: length

    getarrayname_fortran = getarrayname_c(index, cname, length)
    name = char_array_to_string(cname, length)
  end function getarrayname_fortran


  integer(c_int) function getarrayname_c(index, name, length) bind(C,name="getarrayname")
    integer(c_int), intent(in) :: index
    integer(c_int), intent(out) :: length
    character(kind=c_char, len=1), intent(out) :: name(slen)


    integer :: i,j
    character(kind=c_char,len=slen) :: key
    type(arraytype) :: array

    getarrayname_c = -1
    ! This is the index in fortran format
    call indextos(s,index,array)
    ! We need to conver them to C format (char1's)
    key = array%name
    length = len(trim(key))
    name = string_to_char_array(key, length)
    getarrayname_c = 0
  end function getarrayname_c

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

  integer(c_int) function getarrayrank_fortran(name, rank) 
    !DEC$ ATTRIBUTES DLLEXPORT::getarray_rank
    character(kind=c_char,len=*),intent(in) :: name
    integer(c_int), intent(out) :: rank

    ! and we need the string length ....
    integer(c_int) :: length
    character(1), dimension(len(name)) :: cname 
    integer :: i
    length  = len(name)
    cname = string_to_char_array(name,length)
    getarrayrank_fortran = getarrayrank_c(cname, rank, length)
  end function getarrayrank_fortran

  integer(c_int) function getarrayrank_c(name, rank, length) bind(C,name="getarrayrank")
    character(kind=c_char,len=1),intent(in) :: name(length)
    integer(c_int), intent(out) :: rank
    integer(c_int) :: length

    character(len=length) :: key
    integer :: index
    type(arraytype) :: array
    key = char_array_to_string(name, length)
    getarrayrank_c = -1
    index =  chartoindex(key)
    if (index .eq. -1) return
    call indextos(s,index,array)
    rank = array%rank
  end function getarrayrank_c


  integer(c_int) function getarraydimsize_fortran(name, dim, size) 
    !DEC$ ATTRIBUTES DLLEXPORT::getarray_dimsize
    character(kind=c_char,len=*),intent(in) :: name
    integer(c_int), intent(in) :: dim ! dimension number
    integer(c_int), intent(out) :: size

    ! and we need the string length ....
    integer(c_int) :: length
    character(1), dimension(len(name)) :: cname 
    integer :: i
    length  = len(name)
    cname = string_to_char_array(name,length)
    getarraydimsize_fortran = getarraydimsize_c(cname, dim, size, length)
  end function getarraydimsize_fortran

  integer(c_int) function getarraydimsize_c(name, dim, dimsize, length) bind(C,name="getarraydimsize")
    character(kind=c_char,len=1),intent(in) :: name(length)
    integer(c_int), intent(in) :: dim ! dimension number
    integer(c_int), intent(out) :: dimsize
    integer(c_int), intent(in) :: length

    character(len=length) :: key
    integer :: index
    type(arraytype) :: array
    key = char_array_to_string(name, length)
    getarraydimsize_c = -1
    index =  chartoindex(key)
    if (index .eq. -1) return
    call indextos(s,index,array)
    if (array%rank < dim) return
    if (array%rank == 0) then
       dimsize = 0
    elseif (array%rank == 1) then
       if (array%type == 'i') dimsize = size(array%i1)
       if (array%type == 'r') dimsize = size(array%r1)
    elseif (array%rank == 2) then
       if (array%type == 'i') dimsize = size(array%i2, dim)
       if (array%type == 'r') dimsize = size(array%r2, dim)
    elseif (array%rank == 3) then
       if (array%type == 'i') dimsize = size(array%i3, dim)
       if (array%type == 'r') dimsize = size(array%r3, dim)
    elseif (array%rank == 4) then
       if (array%type == 'i') dimsize = size(array%i4, dim)
       if (array%type == 'r') dimsize = size(array%r4, dim)
    endif
  end function getarraydimsize_c

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


  integer(c_int) function get0ddoublearray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    real(c_double), intent(inout) :: x

    type(arraytype) :: array
    integer :: i, index

    get0ddoublearray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    x = array%r0
    get0ddoublearray_fortran = 0
  end function get0ddoublearray_fortran


  integer(c_int) function get0ddoublearray_c(name, x, length) bind(C, name="get0ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::get0ddoublearray_c

    ! use inout otherwise things break
    type(c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), target  :: r0

    get0ddoublearray_c = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    r0 = array%r0
    x = c_loc(r0)
    get0ddoublearray_c = 0
  end function get0ddoublearray_c


  integer(c_int) function get1ddoublearray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    real(c_double), intent(inout) :: x(:)

    type(arraytype) :: array
    integer :: i, index

    get1ddoublearray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    x = array%r1
    get1ddoublearray_fortran = 0
  end function get1ddoublearray_fortran


  integer(c_int) function get1ddoublearray_c(name, x, length) bind(C, name="get1ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::get1ddoublearray_c

    ! use inout otherwise things break
    type(c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), target, allocatable, save, dimension(:)  :: r1

    get1ddoublearray_c = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    
    if (allocated (r1)) deallocate (r1)
    allocate(r1(size(array%r1,1)))
    r1(:) = array%r1(:)

    x = c_loc(r1)
    get1ddoublearray_c = 0
  end function get1ddoublearray_c

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
    real(c_double), target, allocatable, save, dimension(:,:)  :: r2

    get2ddoublearray_c = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    if (allocated (r2)) deallocate (r2)
    allocate(r2(size(array%r2,1), size(array%r2,2)))
    r2(:,:) = array%r2(:,:)
    ! array%r2 => r2
    x = c_loc(r2)
    get2ddoublearray_c = 0
  end function get2ddoublearray_c

  integer(c_int) function get3ddoublearray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    real(c_double), intent(inout) :: x(:,:,:)

    type(arraytype) :: array
    integer :: i, index

    get3ddoublearray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    x = array%r3
    get3ddoublearray_fortran = 0
  end function get3ddoublearray_fortran


  integer(c_int) function get3ddoublearray_c(name, x, length) bind(C, name="get3ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::get3ddoublearray_c

    ! use inout otherwise things break
    type(c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), target, allocatable, save, dimension(:,:,:)  :: r3

    get3ddoublearray_c = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    if (allocated (r3)) deallocate (r3)
    allocate(r3(size(array%r3,1), size(array%r3,2), size(array%r3,3)))
    r3(:,:,:) = array%r3(:,:,:)
    ! array%r3 => r3
    x = c_loc(r3)
    get3ddoublearray_c = 0
  end function get3ddoublearray_c


  integer(c_int) function get4ddoublearray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    real(c_double), intent(inout) :: x(:,:,:,:)

    type(arraytype) :: array
    integer :: i, index

    get4ddoublearray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    x = array%r4
    get4ddoublearray_fortran = 0
  end function get4ddoublearray_fortran


  integer(c_int) function get4ddoublearray_c(name, x, length) bind(C, name="get4ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::get4ddoublearray_c

    ! use inout otherwise things break
    type(c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), target, allocatable, save, dimension(:,:,:,:)  :: r4

    get4ddoublearray_c = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    if (allocated (r4)) deallocate (r4)
    allocate(r4(size(array%r4,1), size(array%r4,2), size(array%r4,3), size(array%r4,4)))
    r4(:,:,:,:) = array%r4(:,:,:,:)
    ! array%r4 => r4
    x = c_loc(r4)
    get4ddoublearray_c = 0
  end function get4ddoublearray_c


  integer(c_int) function get0dintarray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    integer(c_int), intent(inout) :: x

    type(arraytype) :: array
    integer :: i, index

    get0dintarray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    x = array%i0
    get0dintarray_fortran = 0
  end function get0dintarray_fortran

  integer(c_int) function get0dintarray_c(name, x, length) bind(C, name="get0dintarray")
    !DEC$ ATTRIBUTES DLLEXPORT::get0dintarray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    integer(c_int), target  :: i0

    get0dintarray_c = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    i0 = array%i0
    x = c_loc(i0)
    get0dintarray_c = 0
  end function get0dintarray_c

  integer(c_int) function get1dintarray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    integer(c_int), intent(inout) :: x(:)

    type(arraytype) :: array
    integer :: i, index

    get1dintarray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    x = array%i1
    get1dintarray_fortran = 0
  end function get1dintarray_fortran

  integer(c_int) function get1dintarray_c(name, x, length) bind(C, name="get1dintarray")
    !DEC$ ATTRIBUTES DLLEXPORT::get1dintarray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    integer(c_int), target, allocatable, save, dimension(:)  :: i1

    get1dintarray_c = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    if (allocated (i1)) deallocate (i1)
    allocate(i1(size(array%i1,1)))
    i1(:) = array%i1(:)
    ! array%r2 => r2
    x = c_loc(i1)
    get1dintarray_c = 0
  end function get1dintarray_c



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
    integer(c_int), target, allocatable, save, dimension(:,:)  :: i2

    get2dintarray_c = -1
    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    if (allocated (i2)) deallocate (i2)
    allocate(i2(size(array%i2,1), size(array%i2,2)))
    i2(:,:) = array%i2(:,:)
    ! array%r2 => r2
    x = c_loc(i2)
    get2dintarray_c = 0
  end function get2dintarray_c


  integer(c_int) function set0dintarray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    integer(c_int), intent(inout) :: x

    integer :: i, index
    type(arraytype) :: array

    set0dintarray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    array%i0 = x
    set0dintarray_fortran = 0
  end function set0dintarray_fortran

  integer(c_int) function set0dintarray_c(name, x, length) bind(C, name="set0dintarray")
    !DEC$ ATTRIBUTES DLLEXPORT::set0dintarray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    integer(c_int), pointer :: i0

    set0dintarray_c = -1

    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    ! Transform the c pointer into a fortran pointer
    call c_f_pointer(x, i0, shape(array%i0))
    ! Copy the values, or the pointer... not sure.
    array%i0 = i0

    set0dintarray_c = 0
  end function set0dintarray_c

  integer(c_int) function set1dintarray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    integer(c_int), intent(inout) :: x(:)

    integer :: i, index
    type(arraytype) :: array

    set1dintarray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    array%i1 = x
    set1dintarray_fortran = 0
  end function set1dintarray_fortran

  integer(c_int) function set1dintarray_c(name, x, length) bind(C, name="set1dintarray")
    !DEC$ ATTRIBUTES DLLEXPORT::set1dintarray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    integer(c_int), pointer, dimension(:)  :: i1

    set1dintarray_c = -1

    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    ! Transform the c pointer into a fortran pointer
    call c_f_pointer(x, i1, shape(array%i1))
    ! Copy the values, or the pointer... not sure.
    array%i1 = i1

    set1dintarray_c = 0
  end function set1dintarray_c

  integer(c_int) function set2dintarray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    integer(c_int), intent(inout) :: x(:,:)

    integer :: i, index
    type(arraytype) :: array

    set2dintarray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    array%i2 = x
    set2dintarray_fortran = 0
  end function set2dintarray_fortran

  integer(c_int) function set2dintarray_c(name, x, length) bind(C, name="set2dintarray")
    !DEC$ ATTRIBUTES DLLEXPORT::set2dintarray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    integer(c_int), pointer, dimension(:,:)  :: i2

    set2dintarray_c = -1

    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    ! Transform the c pointer into a fortran pointer
    call c_f_pointer(x, i2, shape(array%i2))
    ! Copy the values, or the pointer... not sure.
    array%i2 = i2

    set2dintarray_c = 0
  end function set2dintarray_c

  integer(c_int) function set0ddoublearray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    real(c_double), intent(inout) :: x

    integer :: i, index
    type(arraytype) :: array

    set0ddoublearray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    array%r0 = x
    set0ddoublearray_fortran = 0
  end function set0ddoublearray_fortran

  integer(c_int) function set0ddoublearray_c(name, x, length) bind(C, name="set0ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::set0ddoublearray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), pointer  :: r0

    set0ddoublearray_c = -1

    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    ! Transform the c pointer into a fortran pointer
    call c_f_pointer(x, r0, shape(array%r0))
    ! Copy the values, or the pointer... not sure.
    array%r0 = r0

    set0ddoublearray_c = 0
  end function set0ddoublearray_c

  integer(c_int) function set1ddoublearray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    real(c_double), intent(inout) :: x(:)

    integer :: i, index
    type(arraytype) :: array

    set1ddoublearray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    array%r1 = x
    set1ddoublearray_fortran = 0
  end function set1ddoublearray_fortran

  integer(c_int) function set1ddoublearray_c(name, x, length) bind(C, name="set1ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::set1ddoublearray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), pointer, dimension(:)  :: r1

    set1ddoublearray_c = -1

    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    ! Transform the c pointer into a fortran pointer
    call c_f_pointer(x, r1, shape(array%r1))
    ! Copy the values, or the pointer... not sure.
    array%r1 = r1

    set1ddoublearray_c = 0
  end function set1ddoublearray_c


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


  integer(c_int) function set3ddoublearray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    real(c_double), intent(inout) :: x(:,:,:)

    integer :: i, index
    type(arraytype) :: array

    set3ddoublearray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    array%r3 = x
    set3ddoublearray_fortran = 0
  end function set3ddoublearray_fortran

  integer(c_int) function set3ddoublearray_c(name, x, length) bind(C, name="set3ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::set3ddoublearray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), pointer, dimension(:,:,:)  :: r3

    set3ddoublearray_c = -1

    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    ! Transform the c pointer into a fortran pointer
    call c_f_pointer(x, r3, shape(array%r3))
    ! Copy the values, or the pointer... not sure.
    array%r3 = r3

    set3ddoublearray_c = 0
  end function set3ddoublearray_c


  integer(c_int) function set4ddoublearray_fortran(name,x) 
    USE iso_c_binding
    ! String
    character(kind=c_char,len=*),intent(in) :: name
    ! use inout otherwise things break
    real(c_double), intent(inout) :: x(:,:,:,:)

    integer :: i, index
    type(arraytype) :: array

    set4ddoublearray_fortran = -1
    index =  chartoindex(trim(name))
    if (index .eq. -1) return
    call indextos(s,index,array)
    array%r4 = x
    set4ddoublearray_fortran = 0
  end function set4ddoublearray_fortran

  integer(c_int) function set4ddoublearray_c(name, x, length) bind(C, name="set4ddoublearray")
    !DEC$ ATTRIBUTES DLLEXPORT::set4ddoublearray_c

    ! use inout otherwise things break
    type (c_ptr), intent(inout) :: x
    ! and we need the string length ....
    integer(c_int),value  ,intent(in)    :: length
    ! String
    character(kind=c_char),intent(in) :: name(length)

    character(length) :: myname 
    integer :: index
    type(arraytype) :: array
    real(c_double), pointer, dimension(:,:,:,:)  :: r4

    set4ddoublearray_c = -1

    myname = char_array_to_string(name, length)
    index =  chartoindex(myname)
    if (index .eq. -1) return
    call indextos(s,index,array)
    ! Transform the c pointer into a fortran pointer
    call c_f_pointer(x, r4, shape(array%r4))
    ! Copy the values, or the pointer... not sure.
    array%r4 = r4

    set4ddoublearray_c = 0
  end function set4ddoublearray_c

end module introspection_module
