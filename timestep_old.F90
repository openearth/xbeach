subroutine timestep(s,par)
use params
use spaceparams

IMPLICIT NONE

type(spacepars)			    :: s
type(parameters)		    :: par

integer                     :: i
integer                     :: j
integer                     :: itheta

! Robert new time step criterion
if (par%t==0.0d0) then		! conservative estimate
	par%dt    = par%CFL*min(minval(s%xz(2:s%nx+1)-s%xz(1:s%nx)),minval(s%yz(2:s%ny+1)-s%yz(1:s%ny)))/(maxval(s%hh)*par%g)
else
	par%dt=huge(0.0d0)	! Seed dt
	do itheta=1,s%ntheta
		do j=2,s%ny
			do i=2,s%nx					
				! u direction
					par%dt=min(par%dt,par%CFL*min(s%xu(i+1)-s%xu(i),s%xu(i)-s%xu(i-1))&
												/(max(par%eps,sqrt(par%g*max(s%hu(i,j),s%hu(i+1,j)))+max(abs(s%advecu(i,j)),s%vmageu(i,j)))))
				! v direction
					par%dt=min(par%dt,par%CFL*min(s%yv(j+1)-s%yv(j),s%yv(j)-s%yv(j-1))&
												/(max(par%eps,sqrt(par%g*max(s%hv(i,j),s%hv(i,j+1)))+max(abs(s%advecv(i,j)),s%vmagev(i,j)))))
				! Theta direction
					par%dt=min(par%dt,par%CFL*s%dtheta/max(abs(s%ctheta(i,j,itheta)),tiny(0.0d0)))
				! Riemann equation
				!    par%dt=min(par%dt,par%CFL*min(s%xu(i+1)-s%xu(i),s%xu(i)-s%xu(i-1))&
				!	                            /(max(par%eps,2*sqrt(par%g*max(s%hh(i,j),s%hh(i+1,j)))+max(abs(s%advecu(i,j)),s%vmageu(i,j)))))
			end do
		end do
	end do
end if

par%t=par%t+par%dt;
if(par%t>=par%tnext) then
    par%dt=par%dt-(par%t-par%tnext)
    par%t=par%tnext
    par%tnext=par%tnext+par%tint
end if  

end subroutine