module varianceupdate_module
use xmpi_module
real*8,dimension(:,:,:), allocatable:: meanarrays  ! Keep time average variables 
real*8,dimension(:,:,:), allocatable:: variancearrays     ! Keep variance of variables
real*8,dimension(:,:,:), allocatable:: variancecrossterm  ! Needed for variance calculation
real*8,dimension(:,:,:), allocatable:: variancesquareterm ! Needed for variance calculation
real*8,dimension(:,:,:), allocatable,target:: minarrays   ! Keep time min variables
real*8,dimension(:,:,:), allocatable,target:: maxarrays   ! Keep time max variables





!!!! Routine to calculate time-average variables
contains

subroutine initialize_mean_arrays(nx,ny,nmeanvar)
  
  Implicit none

  integer, intent(in)             :: nx,ny,nmeanvar
  
  allocate(meanarrays(nx+1,ny+1,nmeanvar+2))   ! two extra needed for straight average of Hrms and urms
  allocate(variancearrays(nx+1,ny+1,nmeanvar))
  allocate(variancecrossterm(nx+1,ny+1,nmeanvar))
  allocate(variancesquareterm(nx+1,ny+1,nmeanvar))
  allocate(minarrays(nx+1,ny+1,nmeanvar))
  allocate(maxarrays(nx+1,ny+1,nmeanvar))

  meanarrays = 0.d0
  variancearrays = 0.d0
  variancecrossterm = 0.d0
  variancesquareterm = 0.d0
  minarrays = huge(0.d0)
  maxarrays = -1.d0*huge(0.d0)
end subroutine


subroutine makeaverage(s,sl,par, meanvec, nmeanvar)

  use params
  use spaceparams
  use mnemmodule

  IMPLICIT NONE

  type(parameters), intent(IN)    :: par
  type(spacepars), intent(IN)     :: s
  type(spacepars), intent(IN)     :: sl
  integer*4, intent(IN)           :: nmeanvar    ! number of time-average variables
  integer*4,dimension(:),allocatable, intent(IN)  :: meanvec     ! keep track of which mean variables are used


  integer                         :: i
  real*8                          :: mult
  type(arraytype)                 :: t
  real*8,dimension(:,:),pointer   :: MI,MA
  real*8,dimension(s%nx+1,s%ny+1) :: oldmean,tvar

  ! wwvv this subroutine needs all of the arrays in question

  ! to comfort the compiler
  i=sl%nx
  
  mult = max(par%dt/par%tintm,0.d0) ! Catch initialization at t=0
!Dano  if(xmaster) then
    do i=1,nmeanvar
#ifdef USEMPI
      call space_collect_index(s,sl,meanvec(i))
#endif
     if(xmaster) then
      call indextos(s,meanvec(i),t)
		MI=>minarrays(:,:,i)
		MA=>maxarrays(:,:,i)
      select case(t%name)
        case (mnem_Fx)
		      oldmean=meanarrays(:,:,i)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
            tvar=s%Fx*cos(s%alfa) - s%Fy*sin(s%alfa)
				meanarrays(:,:,i) = meanarrays(:,:,i) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
				                         mult*2.d0*tvar*meanarrays(:,:,i)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
            MA = max(MA,tvar)
			   MI = min(MI,tvar)
        case (mnem_Fy)
		      oldmean=meanarrays(:,:,i)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
            tvar=s%Fx*sin(s%alfa) + s%Fy*cos(s%alfa)
				meanarrays(:,:,i) = meanarrays(:,:,i) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
				                         mult*2.d0*tvar*meanarrays(:,:,i)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
            MA = max(MA,tvar)
			   MI = min(MI,tvar)       
        case (mnem_H)
            meanarrays(:,:,i)=meanarrays(:,:,i) + mult*s%H**2
				oldmean=meanarrays(:,:,nmeanvar+1)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
            tvar=s%H
				meanarrays(:,:,nmeanvar+1) = meanarrays(:,:,nmeanvar+1) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,nmeanvar+1) + &
				                         mult*2.d0*tvar*meanarrays(:,:,nmeanvar+1)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,nmeanvar+1)**2
				MA = max(MA,s%H)
			   MI = min(MI,s%H)
        case (mnem_urms)
            meanarrays(:,:,i)=meanarrays(:,:,i) + mult*s%urms**2
				oldmean=meanarrays(:,:,nmeanvar+2)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
            tvar=s%urms
				meanarrays(:,:,nmeanvar+2) = meanarrays(:,:,nmeanvar+2) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,nmeanvar+2) + &
				                         mult*2.d0*tvar*meanarrays(:,:,nmeanvar+2)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,nmeanvar+2)**2
			   MA = max(MA,s%urms)
			   MI = min(MI,s%urms)
        case (mnem_u)
		      oldmean=meanarrays(:,:,i)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
                tvar=s%u*cos(s%alfa) - s%v*sin(s%alfa)
				meanarrays(:,:,i) = meanarrays(:,:,i) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
				                         mult*2.d0*tvar*meanarrays(:,:,i)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
            MA = max(MA,tvar)
			   MI = min(MI,tvar)
	    case (mnem_gwu)
		        oldmean=meanarrays(:,:,i)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
                tvar=s%gwu*cos(s%alfa) - s%gwv*sin(s%alfa)
				meanarrays(:,:,i) = meanarrays(:,:,i) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
				                         mult*2.d0*tvar*meanarrays(:,:,i)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
               MA = max(MA,tvar)
			   MI = min(MI,tvar)
        case (mnem_v)
		        oldmean=meanarrays(:,:,i)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
                tvar=s%u*sin(s%alfa) + s%v*cos(s%alfa)
				meanarrays(:,:,i) = meanarrays(:,:,i) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
				                         mult*2.d0*tvar*meanarrays(:,:,i)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
                MA = max(MA,tvar)
			    MI = min(MI,tvar)
        case (mnem_gwv)
		        oldmean=meanarrays(:,:,i)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
                tvar=s%gwu*sin(s%alfa) + s%gwv*cos(s%alfa)
				meanarrays(:,:,i) = meanarrays(:,:,i) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
				                         mult*2.d0*tvar*meanarrays(:,:,i)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
                MA = max(MA,tvar)
			    MI = min(MI,tvar)
        case (mnem_ue)
		      oldmean=meanarrays(:,:,i)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
            tvar=s%ue*cos(s%alfa) - s%ve*sin(s%alfa)
				meanarrays(:,:,i) = meanarrays(:,:,i) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
				                         mult*2.d0*tvar*meanarrays(:,:,i)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
            MA = max(MA,tvar)
			   MI = min(MI,tvar)
        case (mnem_ve)
		      oldmean=meanarrays(:,:,i)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
            tvar=s%ue*sin(s%alfa) + s%ve*cos(s%alfa)
				meanarrays(:,:,i) = meanarrays(:,:,i) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
				                         mult*2.d0*tvar*meanarrays(:,:,i)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
            MA = max(MA,tvar)
			   MI = min(MI,tvar)
        case (mnem_uwf)
		      oldmean=meanarrays(:,:,i)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
            tvar=s%uwf*cos(s%alfa) - s%vwf*sin(s%alfa)
				meanarrays(:,:,i) = meanarrays(:,:,i) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
				                         mult*2.d0*tvar*meanarrays(:,:,i)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
            MA = max(MA,tvar)
			   MI = min(MI,tvar)
        case (mnem_vwf)
		      oldmean=meanarrays(:,:,i)
				! oldmean is initiated at 0, so:
				where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
				   oldmean=tiny(0.d0)
				elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
				   oldmean=-1.d0*tiny(0.d0)
				endwhere
            tvar=s%uwf*sin(s%alfa) + s%vwf*cos(s%alfa)
				meanarrays(:,:,i) = meanarrays(:,:,i) + mult*tvar
				variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
				                         mult*2.d0*tvar*meanarrays(:,:,i)
				variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(tvar)**2
				variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
            MA = max(MA,tvar)
			   MI = min(MI,tvar)
        case default
          select case (t%type)
            case ('i')
              select case(t%rank)
                case(2)
					   ! Variance = sum(X^2) - 2*sum(XXm) + sum((Xm)^2)
						! X^2 needs previous instantanious value and current value
						! 2*sum(XXm) must be updated by:
						!     2*sum(XXm) = (2*sum(XXm))(n)/Xm(n)*Xm(n+1)+2*X(n+1)*Xm(n+1)
						! sum((Xm)^2) is updated using Xm(n+1)
                        oldmean=meanarrays(:,:,i)
						! oldmean is initiated at 0, so:
						where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
						   oldmean=tiny(0.d0)
						elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
						   oldmean=-1.d0*tiny(0.d0)
						endwhere
                        meanarrays(:,:,i)=meanarrays(:,:,i) + mult*dble(t%i2)
						variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
						                         mult*2.d0*dble(t%i2)*meanarrays(:,:,i)
						variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*dble(t%i2)**2
						variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
		                MA = max(MA,dble(t%i2))
			            MI = min(MI,dble(t%i2))
                case default
                  write(*,*)'Error in makeaverage, variable"'//t%name// &
                          '" is of wrong rank'
                  call halt_program
              end select  ! rank
            case ('r')
              select case(t%rank)
                case(2)
					   oldmean=meanarrays(:,:,i)
						! oldmean is initiated at 0, so:
						where (oldmean<tiny(0.d0) .and. oldmean>=0.d0)
						   oldmean=tiny(0.d0)
						elsewhere (oldmean>-1.d0*tiny(0.d0) .and. oldmean<0.d0)
						   oldmean=-1.d0*tiny(0.d0)
						endwhere
                        meanarrays(:,:,i)=meanarrays(:,:,i) + mult*t%r2
                        variancecrossterm(:,:,i)=variancecrossterm(:,:,i)/oldmean*meanarrays(:,:,i) + &
						                         mult*2.d0*t%r2*meanarrays(:,:,i)
						variancesquareterm(:,:,i)=variancesquareterm(:,:,i)+mult*(t%r2)**2
						variancearrays(:,:,i)=variancesquareterm(:,:,i)-variancecrossterm(:,:,i)+meanarrays(:,:,i)**2
		                MA = max(MA,t%r2)
			            MI = min(MI,t%r2)
                case default
                  write(*,*)'Error in makeaverage, variable"'//t%name// &
                          '" is of wrong rank'
                  call halt_program
              end select  ! rank
          end select  !type
      end select  ! name
   endif ! xmaster
   enddo ! nmeanvar
 !Dano  endif ! xmaster

end subroutine makeaverage

! Subroutine make time series file name

subroutine makeaveragenames(counter,fnamemean,fnamevar,fnamemin,fnamemax)
  use mnemmodule

  implicit none

  character(99)      :: fnamemean,fnamevar,fnamemin,fnamemax
  integer            :: counter

  fnamemean = trim(mnemonics(counter))//'_mean.dat'
  fnamevar  = trim(mnemonics(counter))//'_var.dat'
  fnamemin  = trim(mnemonics(counter))//'_min.dat'
  fnamemax  = trim(mnemonics(counter))//'_max.dat'

end subroutine makeaveragenames

subroutine makecrossvector(s,sl,crossvararray,nvar,varindexvec,mg,cstype)
  use params
  use spaceparams
  use mnemmodule

  IMPLICIT NONE

  type(spacepars), intent(in)       :: s,sl    
  real*8,dimension(:,:)             :: crossvararray
  integer,intent(in)                :: mg,cstype,nvar
  integer,dimension(nvar),intent(in):: varindexvec
  type(arraytype)                   :: t

  integer                           :: i

  crossvararray=-999.d0
  do i=1,nvar
#ifdef USEMPI
    call space_collect_index(s,sl,varindexvec(i))
#endif  
    if (xmaster) then
	   call indextos(s,varindexvec(i),t)	       
      select case(t%type)
		  case('r')
		    if(cstype==0) then
			   crossvararray(i,:)=t%r2(:,mg)
			 else
			   crossvararray(i,:)=t%r2(mg,:)
			 endif
		  case('i')
			 if(cstype==0) then
			   crossvararray(i,:)=t%i2(:,mg)
			 else
			   crossvararray(i,:)=t%i2(mg,:)
			 endif
	   end select
    endif
  enddo
#ifdef USEMPI
  call xmpi_bcast(crossvararray)
#endif 
end subroutine makecrossvector



end module
