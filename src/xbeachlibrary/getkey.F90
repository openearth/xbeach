module getkey_module
  ! This module makes parameters available through introspection/reflection
  ! You can ask for all available parameters using getkeys 
  ! You can ask for the value of a specific parameter using getkey
  ! The naming is still a bit inconsistent so that might change in the future
  use typesandkinds
  use mnemmodule
  implicit none
  save
  ! analogue to the arraytype in mnemonic
  type parameter
     character type         ! 'c', 'i' or 'r': integer or real*8
     integer rank           ! 0,1,2,3,4
     character(slen) :: name     ! 'v','ve', .....
     character(slen) :: units     ! m, following udunits convention
     character(slen) :: description
     character(slen), dimension(maxrank) :: dimensions ! the dimensions of the variable, for example (s%nx, s%ny)

     character(slen), pointer :: c0 ! pointer to characters
     character(slen), dimension(:), pointer :: c1 ! pointer to array of characters

     real*8, pointer             :: r0  ! pointer to real8 scalar
     real*8, dimension(:), pointer :: r1  ! pointer to real8 (:)

     integer, pointer             :: i0  ! pointer to integer scalar
     integer, dimension(:), pointer :: i1  ! pointer to integer (:)
  end type parameter

contains
  ! Lookup the index and valuetype of a parameter, currently only works for rank 0
  subroutine getkey_indextype(par, key, index, type)
    use params
    implicit none
    type(parameters), intent(in) :: par 
    character(len=*), intent(in) :: key
    integer, intent(out)         :: index
    character, intent(out)       :: type
    integer                      :: i
    include 'getkey.gen'
    index = -1
    type = ''
    do i=1,nintegerkeys
       if (integerkeys(i) == key) then
          index = i
          type = 'i'
          return
       end if
    end do
    do i=1,nrealkeys
       if (realkeys(i) == key) then
          index = i
          type = 'r'
          return
       end if
    end do
    do i=1,ncharacterkeys
       if (characterkeys(i) == key) then
          index = i
          type = 'c'
          return
       end if
    end do
  end subroutine getkey_indextype

  integer function getkey(par, key, value)
    ! This subroutine returns a value (parameter), given a key (parameter name).
    use params
    use mnemmodule
    implicit none

    type(parameters), intent(in) :: par 
    character(len=*), intent(in) :: key
    type(parameter), intent(inout) :: value
    ! Ok this is ugly, I changed the targets to save, so the pointer is still available after the function ends
    ! But this means that if you call getkey twice, then your old value will now refer to the new value
    character(slen),target, save :: charvalue 
    integer, target, save :: intvalue
    real*8, target, save :: realvalue
    integer :: index
    character :: type
    include 'getkey.gen'
    ! This sets index and type
    call getkey_indextype(par, key, index, type)
    getkey = -1
    if (index .eq. -1 ) return
    value%name = key
    value%rank = 0
    value%type = type
    if (type == 'c') then
       charvalue = charactervalues(index)
       value%c0 => charvalue
    elseif (type == 'i') then
       intvalue = integervalues(index)
       value%i0 => intvalue
    elseif (type == 'r') then
       realvalue = realvalues(index)
       value%r0 => realvalue
    end if
    getkey = 0
  end function getkey

  subroutine getnkeys(par, n)
    use params
    implicit none
    type(parameters), intent(in) :: par 
    integer, intent(out) :: n 
    include 'getkey.gen'
    n = ncharacterkeys+nintegerkeys+nrealkeys
  end subroutine getnkeys

  subroutine getkeys(par, keys)
    ! This subroutine gives a list of all available keys (parameter names)
    use params
    implicit none
    type(parameters), intent(in) :: par 
    character(slen), dimension(:), allocatable, intent(out) :: keys

    include 'getkey.gen'
    allocate(keys(ncharacterkeys+nintegerkeys+nrealkeys))
    keys(1:ncharacterkeys) = characterkeys(1:ncharacterkeys)
    keys(ncharacterkeys+1:ncharacterkeys+nintegerkeys) = integerkeys(1:nintegerkeys)
    keys(ncharacterkeys+nintegerkeys+1:ncharacterkeys+nintegerkeys+nrealkeys) = realkeys(1:nrealkeys)


  end subroutine getkeys





end module getkey_module

