module mnemmodule
use typesandkinds
include 'mnemonic.gen'

type arraytype

  character type         ! 'i' or 'r': integer or real*8
  character btype        ! 'b' or 'd': broadcast or distribute
  integer rank           ! 0,1,2,3,4
  character(slen) :: name     ! 'v','ve', .....
  character(slen) :: units     ! m, following udunits convention
  character(slen) :: standardname
  character(slen) :: description
  character(slen), dimension(maxrank) :: dimensions ! the dimensions of the variable, for example (s%nx, s%ny)

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
  use logging_module
  implicit none

  character(len=*), intent(in)  :: line

  chartoindex = -1
  
include 'chartoindex.gen'
  if (chartoindex == -1) then
     call writelog('sl','','Could not find index for: ', line)
  end if
end function chartoindex

subroutine printvar(t,lid,eid,wid)
type (arraytype) :: t
integer, intent(in) :: lid,eid, wid

write(*,*)'Name:    ',trim(t%name)
write(*,*)'Type:    ',t%type
write(*,*)'Btype:   ',t%btype
write(*,*)'Rank:    ',t%rank

write(lid,*)'Name:    ',trim(t%name)
write(lid,*)'Type:    ',t%type
write(lid,*)'Btype:   ',t%btype
write(lid,*)'Rank:    ',t%rank

write(eid,*)'Name:    ',trim(t%name)
write(eid,*)'Type:    ',t%type
write(eid,*)'Btype:   ',t%btype
write(eid,*)'Rank:    ',t%rank

write(wid,*)'Name:    ',trim(t%name)
write(wid,*)'Type:    ',t%type
write(wid,*)'Btype:   ',t%btype
write(wid,*)'Rank:    ',t%rank

end subroutine printvar
end module mnemmodule
