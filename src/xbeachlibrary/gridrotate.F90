! just testing ...
module gridrotate_module
use params
use spaceparams
use mnemmodule
contains

real*8 elemental function gridrotate(par,s,t,x)
type(parameters), intent(in)        :: par
type(spacepars), intent(in)         :: s
type(arraytype), intent(in)         :: t
real*8, intent(in) :: x

real*8, parameter :: pi=4*atan(1.0d0)

if (par%rotate .eq. 1) then
  select case(t%name)
  case(mnem_theta)
    !gridrotate=270-(t%r1*(180/pi))
    gridrotate = sin(x)
  case(mnem_theta0)
  !  gridrotate=270-(t%r1*(180/pi))
    gridrotate = sin(x)
  case default
  !  gridrotate=t%r1
    gridrotate = sin(x)
  end select
else
  !gridrotate=t%r1
    gridrotate = sin(x)
endif



end function gridrotate
end module gridrotate_module
