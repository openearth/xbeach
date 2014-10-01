module indextos_module
contains
subroutine indextos(s,index,t)
  use mnemmodule
  use logging_module
  use spaceparamsdef
  implicit none
  type (spacepars), intent(in),target :: s
  integer, intent(in)                 :: index
  type(arraytype), intent(out)        :: t

  ! wwvv
  ! it appears that compilers take long time to
  ! analyse pointer assignments in user defined types,
  !  resulting in long compile times.
  ! therefore, first a local pointer is assigned, and at
  !  the end the appropriate pointer in t is assigned.

  real*8, pointer :: pointr0 
  real*8, pointer :: pointr1(:)
  real*8, pointer :: pointr2(:,:)
  real*8, pointer :: pointr3(:,:,:)
  real*8, pointer :: pointr4(:,:,:,:)
  integer, pointer :: pointi0
  integer, pointer :: pointi1(:)
  integer, pointer :: pointi2(:,:)

  if (index .lt. 1 .or. index .gt. numvars) then
    call writelog('els','(a,i3,a)','invalid index ',index,' in indextos. Program will stop')
    call halt_program
  endif

  select case(index)
  include 'indextos.gen'
  end select

  select case(t%type)
  case("r")
    select case(t%rank)
    case(0)
      t%r0 => pointr0
    case(1)
      t%r1 => pointr1
    case(2)
      t%r2 => pointr2
    case(3)
      t%r3 => pointr3
    case(4)
      t%r4 => pointr4
    case default
      print *,'error'
    end select
  case("i")
    select case(t%rank)
    case(0)
      t%i0 => pointi0
    case(1)
      t%i1 => pointi1
    case(2)
      t%i2 => pointi2
    case default
      print *,'error'
    end select
  case default
    print *,'error'
  end select

end subroutine indextos
end module indextos_module
