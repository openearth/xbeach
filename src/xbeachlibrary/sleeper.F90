! provides subroutine myusleep(n)
!          sleeps for n microseconds
module sleeper
#if __GFORTRAN__
   use, intrinsic :: iso_c_binding
   implicit none
   save
   interface
      subroutine usleep(i) bind(C,name='usleep')
         import
         integer(c_int), value :: i
      end subroutine usleep
   end interface
#else
   implicit none
   save
#endif
   private
   public:: myusleep

contains
#if __GFORTRAN__
   subroutine myusleep(i)
      integer, intent(in) :: i
      call usleep(i)
   end subroutine myusleep
#elif __INTEL_COMPILER
   subroutine myusleep(i)
      integer, intent(in) :: i
      integer j
      j = i/1000
      j = max(1,j)
      call sleepqq(j)
   end subroutine myusleep
#elif __PGI
   subroutine myusleep(i)
      integer, intent(in) :: i
      integer j
      j = i/1000
      j = max(1,j)
      call sleepqq(j)
   end subroutine myusleep
#else
   subroutine myusleep(i)
      ! this is the final resort: do not sleep at all
      integer, intent(in) :: i
      integer j
      j = i/1000000
      j = max(0,j)
      ! j = max(1,j) ! will sleep for one second
      call sleep(j)
   end subroutine myusleep
#endif

end module sleeper
