module timestep_module
contains
subroutine timestep(s,par,it)
use params
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars)                     :: s
type(parameters)                    :: par

integer                     :: it
integer                     :: i
integer                     :: j
integer                     :: itheta
real*8                                          :: mdx,mdy

! Robert new time step criterion
if (par%t==0.0d0) then          ! conservative estimate
  par%dt    = par%CFL*min(minval(s%xz(2:s%nx+1)-s%xz(1:s%nx)),   &
                          minval(s%yz(2:s%ny+1)-s%yz(1:s%ny)))   &
                          /(maxval(s%hh)*par%g)
else
  par%dt=huge(0.0d0)      ! Seed dt
  do j=2,s%ny
    do i=2,s%nx     
      ! u-points
      mdx=(s%xz(i+1)-s%xz(i))
      mdy=min((s%yz(j+1)-s%yz(j)),(s%yz(j)-s%yz(j-1)))
      ! x-component
      par%dt=min(par%dt,mdx/(sqrt(par%g*s%hu(i,j))+abs(s%uu(i,j))))
      ! y-component
      par%dt=min(par%dt,mdy/(sqrt(par%g*s%hu(i,j))+abs(s%vu(i,j))))
      
      ! v-points
      mdx=min((s%xz(i+1)-s%xz(i)),(s%xz(i)-s%xz(i-1)))
      mdy=(s%yz(j+1)-s%yz(j))
      ! x-component
      par%dt=min(par%dt,mdx/(sqrt(par%g*s%hv(i,j))+abs(s%uv(i,j))))
      ! y-component
      par%dt=min(par%dt,mdy/(sqrt(par%g*s%hv(i,j))+abs(s%vv(i,j))))
    enddo
  enddo
  par%dt=par%dt*par%CFL*0.5d0
  if (par%instat>0 .and. par%instat/=40) then
    par%dt=min(par%dt,par%CFL*s%dtheta/(maxval(maxval(abs(s%ctheta),3)*real(s%wetz))+tiny(0.0d0)))
!    !do j=2,s%ny
!    !  do i=2,s%nx     
!    !    do itheta=1,s%ntheta
!    !      ! Theta points
!    !      par%dt=min(par%dt,par%CFL*s%dtheta/max(abs(s%ctheta(i,j,itheta)),tiny(0.0d0)))
!    !    end do
!    !  end do
!    !end do
  end if
end if
! wwvv: In the mpi version par%dt will be calculated different
! on different processes. So, we take the minimum of all dt's
#ifdef USEMPI
    call xmpi_allreduce(par%dt,MPI_MIN)
#endif

par%t=par%t+par%dt
if(par%t>=par%tnext) then
    par%dt=par%dt-(par%t-par%tnext)
    par%t=par%tnext
    it=it+1
!    write(*,*)'time (s) = ',par%t,' dt (s) = ',par%dt,' output step ',it
!else
!    write(*,*)'time (s) = ',par%t,' dt (s) = ',par%dt
end if  

end subroutine timestep
end module timestep_module
