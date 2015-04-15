!
! demonstration how to use the mnemonics stuff
!
program demo
   use spaceparams
   use typesandkinds
   use mnemmodule
   implicit none

   type (spacepars) :: s

   integer, parameter :: nnames = 4
   character(slen), dimension(nnames),parameter :: names=(/'xz','x','cx','nx'/)
   integer :: i
   !
   ! alternatively, to let the compiler catch typing errors:

   !character(slen), dimension(nnames), parameter :: names=(/mnem_xz,mnem_x,mnem_cx,mnem_nx/)

   ! allocate scalars:
   call space_alloc_scalars(s)

   ! allocate some arrays in s

   allocate(s%xz(100))

   allocate(s%x(120,10))
   allocate(s%cx(10,10,10))
   !
   ! note, there is a subroutine space_alloc_arrays in
   ! spaceparams.F90 to allocate all arrays in s.
   ! That subroutine needs filled-in variables i s and par
   !


   ! Give the variables a value, using their names in array names

   do i=1,nnames
      call setvar(s,names(i),dble(i))
   enddo

   ! and show that the values are there:

   write(*,*)'nx:',s%nx
   write(*,*)'nxz:',s%xz(1)
   write(*,*)'x:',s%x(1,1)
   write(*,*)'cx:',s%cx(1,1,1)

end program demo

! example code to demonstrate indextos and chartoindex
! The first element of the named variable will be set to
! zero
subroutine setvar(s,name,value)
   use spaceparams
   use mnemmodule
   implicit none
   type(spacepars) :: s
   character(len=*) :: name
   real*8 :: value

   integer :: index
   type (arraytype) :: t

   index = chartoindex(name)   !determine index from name
   call indextos(s,index,t)    !get info and pointer
   write(*,*)'setting '//trim(name)//' to ',value
   write(*,*)'some properties of '//trim(name)//':'
   call printvar(t)
   !
   ! print some
   ! depending on s%type and s%rank, we have now a pointer
   ! to the desired variable
   !
   ! Fortran90 pointer type checking is quite strict: no void pointer
   ! available. So we have to catch every type and rank of array.
   ! See also mnemonic.F90
   !

   select case (t%type)
    case ('r')      ! type is integer
      select case (t%rank)
       case(0)              ! scalar
         t%r0 = value
       case(1)              ! (:)
         t%r1 = value
       case(2)              ! (:,:)
         t%r2 = value
       case(3)              ! (:,:,:)
         t%r3 = value
       case(4)              ! (:,:,:,:)
         t%r4 = value
      end select
    case ('i')     ! type is real*8
      select case (t%rank)
       case(0)              ! scalar
         t%i0 = value
       case(1)              ! (:)
         t%i1 = value
       case(2)              ! (:,:)
         t%i2 = value
       case(3)              ! (:,:,:)
         t%i3 = value
       case(4)              ! (:,:,:,:)
         t%i4 = value
      end select
   end select
end subroutine setvar

