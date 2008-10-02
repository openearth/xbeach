module mnemmodule
include 'mnemonic.gen'
type arraytype

  character type         ! 'i' or 'r': integer or real*8
  character btype        ! 'b' or 'd' '2': broadcast or distribute or 
                         !                 distribute umean
  integer rank           ! 0,1,2,3,4
  character(len=maxnamelen) :: name     ! 'v','ve', .....

  real*8, pointer             :: r0  ! pointer to real8 scalar
  real*8, dimension(:), pointer :: r1  ! pointer to real8 (:)
  real*8, dimension(:,:), pointer :: r2  ! pointer to real8 (:,:)
  real*8, dimension(:,:,:), pointer :: r3  ! pointer to real8 (:,:,:)
  real*8, dimension(:,:,:,:), pointer :: r4  ! pointer to real8 (:,:,:,:)

  integer, pointer             :: i0  ! pointer to integer scalar
  integer, dimension(:), pointer :: i1  ! pointer to integer (:)
  integer, dimension(:,:), pointer :: i2  ! pointer to integer (:,:)
  integer, dimension(:,:,:), pointer :: i3  ! pointer to integer (:,:,:)
  integer, dimension(:,:,:,:), pointer :: i4  ! pointer to integer (:,:,:,:)
  
end type arraytype

contains

integer function chartoindex(line)

  implicit none

  character(len=*), intent(in)  :: line

  chartoindex = -1

include 'chartoindex.gen'

end function chartoindex

subroutine printvar(t)
type (arraytype) :: t

write(*,*)'Name:    ',trim(t%name)
write(*,*)'Type:    ',t%type
write(*,*)'Btype:   ',t%btype
write(*,*)'Rank:    ',t%rank

end subroutine printvar
end module mnemmodule
