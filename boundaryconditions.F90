module boundaryconditions
contains
subroutine wave_bc(sg,sl,par,newstatbc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Copyright (C) 2007 UNESCO-IHE, WL|Delft Hydraulics and Delft University !
! Dano Roelvink, Ap van Dongeren, Ad Reniers, Jamie Lescinski,            !
! Jaap van Thiel de Vries, Robert McCall                                  !       
!                                                                         !
! d.roelvink@unesco-ihe.org                                               !
! UNESCO-IHE Institute for Water Education                                !
! P.O. Box 3015                                                           !
! 2601 DA Delft                                                           !
! The Netherlands                                                         !
!                                                                         !
! This library is free software; you can redistribute it and/or           !
! modify it under the terms of the GNU Lesser General Public              !
! License as published by the Free Software Foundation; either            !
! version 2.1 of the License, or (at your option) any later version.      !
!                                                                         !
! This library is distributed in the hope that it will be useful,         !
! but WITHOUT ANY WARRANTY; without even the implied warranty of          !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        !
! Lesser General Public License for more details.                         !
!                                                                         !
! You should have received a copy of the GNU Lesser General Public        !
! License along with this library; if not, write to the Free Software     !
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307     !
! USA                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use params
use waveparams
use spaceparams
use interp
use wave_timestep_module
use xmpi_module
use readkey_module

IMPLICIT NONE

type(spacepars), target                     :: sg, sl
type(spacepars), pointer                    :: s
type(parameters)                            :: par
type(waveparameters)                        :: wp

integer, save                               :: nt
integer                                     :: i,new,reclen,wordsize
integer, save                               :: old
integer, save                               :: recpos
integer                                     :: j
integer                                     :: itheta
integer                                     :: E_idx
integer                                     :: dummy
real*8                                      :: E1,ei,dum,Hm0, dum1, spreadpar, bcdur, dum2
real*8, save                                :: dtbcfile,rt,bcendtime
real*8                                      :: em,tshifted,tnew
real*8,dimension(:)     ,allocatable,save   :: e01       ! [J/m2/rad] directional distribution of wave energy at boundary
real*8,dimension(:)     ,allocatable,save   :: fac1,fac2
real*8,dimension(:)     ,allocatable,save   :: tE,dataE,databi
real*8,dimension(:,:)   ,allocatable,save   :: ht
real*8,dimension(:)     ,allocatable,save   :: q1,q2,q
real*8,dimension(:)     ,allocatable,save   :: wcrestpos
real*8,dimension(:,:)   ,allocatable,save   :: ee1, ee2
real*8,dimension(:)     ,allocatable,save   :: gq1,gq2,gq
real*8,dimension(:,:)   ,allocatable,save   :: gee1, gee2
character*1                                 :: bl
character*80                                :: ebcfname,qbcfname,fname
real*8                                      :: E0
real*8,dimension(:),allocatable,save        :: dist,factor
logical                                     :: startbcf,newstatbc


include 's.ind'
s=>sl
#ifndef USEMPI
s=>sg
#endif
include 's.inp'

if (.not. allocated(fac1)) then
   allocate(ht      (2,ny+1))
   allocate(fac1    (nx))
   allocate(fac2    (nx))
   allocate(e01(1:ntheta))
   allocate(dist(1:ntheta))
   allocate(factor(1:ntheta))
   allocate(wcrestpos(nx+1))
endif
!
!  GENERATE AND READ-IN WAVE BOUNDARY CONDITIONS
!
! added for bound long wave comp Ad 28 march 2006
dtheta = par%dtheta*par%px/180
startbcf=.false.
if(abs(par%t-par%dt)<1.d-6) then
    if (xmaster) then
      write(*,*)'Setting up boundary conditions'
    endif
    startbcf=.true.                     ! trigger read from bcf for instat 3,4,5,7
    bcendtime=huge(0.0d0)               ! initial assumption for instat 3,4,5,7
    if (par%instat==2) then
       if(xmaster) then
		  open( unit=7, file='bc/gen.ezs')
       endif
       if (xmaster) then
5        continue
         read(7,'(a)')bl
         if(bl.eq.'*') goto 5
         read(7,*)nt
       endif
#ifdef USEMPI
       call xmpi_bcast(nt)
#endif
       allocate(dataE  (nt))
       allocate(tE     (nt))
       do i=1,nt
          if(xmaster) then
            read(7,*) tE(i),dum,dataE(i)
          endif
#ifdef USEMPI
          call xmpi_bcast(tE(i))
          call xmpi_bcast(dataE(i))
#endif
       end do
       if (xmaster) then
         close(7)
       endif
       par%Emean=sum(dataE)/nt
    elseif (par%instat==3) then
       if (xmaster) then
		  open( unit=7, file='bc/gen.ezs')
       endif
       if (xmaster) then
6        continue
         read(7,'(a)')bl
         if(bl.eq.'*') goto 6
         read(7,*)nt
       endif
#ifdef USEMPI
       call xmpi_bcast(nt)
#endif
       allocate(dataE  (nt))
       allocate(databi (nt))
       allocate(tE     (nt))
       do i=1,nt
          if(xmaster) then
            read(7,*) tE(i),databi(i),dataE(i)
          endif
#ifdef USEMPI
          call xmpi_bcast(tE(i))
          call xmpi_bcast(databi(i))
          call xmpi_bcast(dataE(i))
#endif
       end do
       if (xmaster) then
         close(7)
       endif
       par%Emean=sum(dataE)/nt
	elseif (par%instat==40) then
       if (xmaster) then
         call readkey('params.txt','bcfile',fname)
		 call checkbcfilelength(par,fname)
	     open( unit=7, file=fname)
!	     open( unit=7, file='jonswap1.txt')
	     read(7,*) Hm0, par%Trep,par%dir0, dum1, spreadpar, bcendtime, dum2
	     par%Hrms = Hm0/sqrt(2.d0)
	     par%m = 0.5d0*spreadpar
         bcendtime=bcendtime/max(par%morfac,1.d0)
	     s%theta0=(1.5d0*par%px-s%alfa)-par%dir0*atan(1.d0)/45.d0
       endif
	   newstatbc=.true.
#ifdef USEMPI
	   call xmpi_bcast(bcendtime)
	   call xmpi_bcast(par%Hrms)
	   call xmpi_bcast(par%Trep)
	   call xmpi_bcast(par%m)
	   call xmpi_bcast(s%theta0)
#endif  
       do itheta=1,s%ntheta
           s%sigt(:,:,itheta) = 2.d0*par%px/par%Trep
       end do
       s%sigm = sum(s%sigt,3)/s%ntheta
       call dispersion(par,s)     
	   
    elseif ((par%instat==4.or.par%instat==41).and.xmaster) then
       call makebcf(par,sg,wp)
    elseif (par%instat==5.and.xmaster) then
       call makebcf(par,sg,wp)
    elseif (par%instat==6.and.xmaster) then
       call makebcf(par,sg,wp) 
    elseif (par%instat==7.and.xmaster) then
       par%listline=1
    endif
    !
    ! Directional distribution
    !
    if(par%instat==0 .or. par%instat==1 .or. par%instat==2 .or. par%instat==3 .or. par%instat==40)then
        dist=(cos(theta-theta0))**par%m
        do i=1,ntheta
            if(abs(theta(i)-theta0)>par%px/2) then
                dist(i)=0
            end if
        end do
        E0=par%rhog8*par%Hrms**2
    
        ! energy density distribution

        if (sum(dist)>0.d0) then 
           factor = (dist/sum(dist))/dtheta
        else
           factor=0.d0
        endif
        e01    = factor*E0;
        e01    = max(e01,0.0d0);

        if (abs(theta0)<1.d-9) theta0=1.d-9
        par%Llong=par%Tlong*cg(1,1)*cos(theta0)

    endif
    if (xmaster) then
        write(*,*)'Boundary conditions complete, starting computation'
    endif
end if

if (par%t .ge. bcendtime) then  ! Recalculate bcf-file 
    if (par%instat==40) then
       if (xmaster) then
         write(*,*) 'Reading new wave conditions'
	     read(7,*) Hm0, par%Trep,par%dir0, dum1, spreadpar, bcdur, dum2
         par%Hrms = Hm0/sqrt(2.d0)
	     par%m = 0.5d0*spreadpar
	     bcendtime=bcendtime+bcdur/max(par%morfac,1.d0)
	     s%theta0=(1.5d0*par%px-s%alfa)-par%dir0*atan(1.d0)/45.d0
       endif
	   newstatbc=.true.
#ifdef USEMPI
	   call xmpi_bcast(bcendtime)
	   call xmpi_bcast(par%Hrms)
	   call xmpi_bcast(par%Trep)
	   call xmpi_bcast(par%m)
	   call xmpi_bcast(s%theta0)
#endif
	          
       do itheta=1,s%ntheta
           s%sigt(:,:,itheta) = 2*par%px/par%Trep
       end do
       s%sigm = sum(s%sigt,3)/s%ntheta
       call dispersion(par,s)
	
       dist=(cos(theta-theta0))**par%m
       do i=1,ntheta
           if(abs(theta(i)-theta0)>par%px/2.d0) then
               dist(i)=0
           end if
       end do
       E0=par%rhog8*par%Hrms**2
    
       ! energy density distribution

       if (sum(dist)>0.d0) then 
          factor = (dist/sum(dist))/dtheta
       else
          factor=0.d0
       endif
       e01    = factor*E0;
       e01    = max(e01,0.0d0);
    elseif ((par%instat==4.or.par%instat==41).and.xmaster) then
!    if ((par%instat==4.or.par%instat==41).and.xmaster) then
        close(71)
        close(72)
        call makebcf(par,sg,wp)
        startbcf=.true.
    elseif (par%instat==5.and.xmaster) then 
        close(71)
        close(72)
        call makebcf(par,sg,wp)
        startbcf=.true.
    elseif (par%instat==6.and.xmaster) then 
        close(71)
        close(72)
        call makebcf(par,sg,wp)
        startbcf=.true.
    elseif (par%instat==7.and.xmaster) then
        close(71)
        close(72)
        startbcf=.true.
        if (par%t <= (par%tstop-par%dt)) then
            par%listline=par%listline+1
        end if
    end if
#ifdef USEMPI
       call xmpi_bcast(startbcf)
#endif
end if  

!
! COMPUTE WAVE BOUNDARY CONDITIONS CURRENT TIMESTEP
!
! instat = 0 or 40 => stationary wave conditions
! instat = 1 => wave energy fluctuations associated with Tlong in params.txt; bound long wave from LHS 1962
! instat = 2 => wave energy time series from Gen.ezs; bound long wave from LHS 1962
! instat = 3 => wave energy time series and (bound) long waves from Gen.ezs
! instat = 4 or 41 => directional wave energy time series for Jonswap spectrum from user specified file; bound long wave from van Dongeren, 19??
! instat = 5 => directional wave energy time series from SWAN 2D spectrum file; bound long wave from van Dongeren, 19??    
! instat = 6 => directional wave energy time series from spectrum file; bound long wave from van Dongeren, 19??
! instat = 7 => as instat = 4/5/6; reading from previously computed wave boundary condition file.
if (par%instat==0 .or. par%instat==40) then
   do j=1,ny+1
       ee(1,j,:)=e01*min(par%t/par%taper,1.0d0)
       bi(1) = 0.0d0
       ui(1,j) = 0.0d0
   end do
elseif (par%instat==1) then
   do j=1,ny+1
       ee(1,j,:)=e01*0.5d0*(1.d0+cos(2*par%px*(par%t/par%Tlong-sin(theta0)*y(1,j)/par%Llong))) *min(par%t/par%taper,1.d0)
       em = (sum(0.5d0*e01))*dtheta *min(par%t/par%taper,1.d0)
       ei =  sum(ee(1,j,1:ntheta))*dtheta
       bi(1) = -(2*cg(1,j)/c(1,j)-0.5d0)*(em-ei)/(cg(1,j)**2-par%g*hh(1,j))/par%rho
       ht=s%zs0(1:2,:)-zb(1:2,:)
       ui(1,j) = cg(1,j)*bi(1)/ht(1,j)*cos(theta0)
   end do
elseif (par%instat==2) then
   do j=1,ny+1
      if (abs(theta0)<1e-3) then
         call linear_interp(tE,dataE,nt,par%t,E1,E_idx)
      else
         tshifted=max(par%t-(y(1,j)-y(1,1))*sin(theta0)/cg(1,1),0.d0)
         call linear_interp(tE,dataE,nt,tshifted,E1,E_idx) 
      endif
      ee(1,j,:)=e01*E1/max(par%Emean,0.000001d0)*min(par%t/par%taper,1.d0)
      em = par%Emean *min(par%t/par%taper,1.d0)
      ei = sum(ee(1,j,:))*dtheta
      bi(1) = -(2*cg(1,j)/c(1,j)-0.5d0)*(em-ei)/(cg(1,j)**2-par%g*hh(1,j))/par%rho
      ht=s%zs0(1:2,:)-zb(1:2,:)
      ui(1,j) = cg(1,j)*bi(1)/ht(1,j)*cos(theta0)
   end do
elseif (par%instat==3) then
    ht=s%zs0(1:2,:)-zb(1:2,:)
    do j=1,ny+1
        if (abs(theta0)<1e-3) then
           call linear_interp(tE,dataE,nt,par%t,E1,E_idx) 
           call linear_interp(tE,databi,nt,par%t,bi(1),E_idx)
        else
           tshifted=max(par%t-(y(1,j)-y(1,1))*sin(theta0)/cg(1,1),0.d0)
           call linear_interp(tE,dataE,nt,tshifted,E1,E_idx) 
           call linear_interp(tE,databi,nt,tshifted,bi(1),E_idx)
        endif
        ee(1,j,:)=e01*E1/max(par%Emean,0.000001d0)*min(par%t/par%taper,1.d0)
        ui(1,j) = cg(1,j)*bi(1)/ht(1,j)*cos(theta0)*min(par%t/par%taper,1.d0)
        if (par%carspan==1) then
          ui(1,j) = sqrt(par%g/ht(1,j))*bi(1)! Carrier and Greenspan
        endif

    end do
elseif ((par%instat==4).or.(par%instat==41).or.(par%instat==5) .or. (par%instat==6) .or. (par%instat==7)) then  
    ! open file if first time
    if (startbcf) then
        if(xmaster) then
          open(53,file='ebcflist.bcf',form='formatted',position='rewind')
          open(54,file='qbcflist.bcf',form='formatted',position='rewind')
        endif
        if (xmaster) then
           do i=1,par%listline
              read(53,*)bcendtime,rt,dtbcfile,par%Trep,s%theta0,ebcfname
              read(54,*)bcendtime,rt,dtbcfile,par%Trep,s%theta0,qbcfname
           enddo  ! wwvv strange
        endif
#ifdef USEMPI
            call xmpi_bcast(bcendtime)
            call xmpi_bcast(rt)
            call xmpi_bcast(dtbcfile)
            call xmpi_bcast(par%Trep)
            call xmpi_bcast(s%theta0)
            call xmpi_bcast(ebcfname)
#endif
        if (xmaster) then
          close(53)
          close(54)
        endif
        ! Robert and Jaap : Initialize for new wave conditions
        ! par%Trep = par%Trep
        ! par%omega = 2*par%px/par%Trep            
        do itheta=1,s%ntheta
            s%sigt(:,:,itheta) = 2*par%px/par%Trep
        end do
        s%sigm = sum(s%sigt,3)/s%ntheta
        call dispersion(par,s)     
        ! End initialize
        if (xmaster) then
		  inquire(iolength=wordsize) 1.d0
          reclen=wordsize*(sg%ny+1)*(s%ntheta)
          open(71,file=ebcfname,status='old',form='unformatted',access='direct',recl=reclen)
          reclen=wordsize*(sg%ny+1)
          open(72,file=qbcfname,status='old',form='unformatted',access='direct',recl=reclen)
        endif
!
! wwvv note that we need the global value of ny here
!
! masterprocess reads and distributes
!
        if (xmaster) then
          if (.not. allocated(gq1) ) then
              allocate(gq1(sg%ny+1),gq2(sg%ny+1),gq(sg%ny+1))
              allocate(gee1(sg%ny+1,ntheta),gee2(sg%ny+1,ntheta))
          endif
        else
          if (.not. allocated(gq1) ) then ! to get valid addresses for
                                          ! gq1, gq2, gq, gee1, gee2
              allocate(gq1(1),gq2(1),gq(1))
              allocate(gee1(1,ntheta),gee2(1,ntheta))
          endif
        endif
        if (.not. allocated(q1) ) then
            allocate(q1(ny+1),q2(ny+1),q(ny+1))
            allocate(ee1(ny+1,ntheta),ee2(ny+1,ntheta))
        end if
        if (xmaster) then
          read(71,rec=1)gee1       ! Earlier in time
          read(71,rec=2)gee2       ! Later in time
          read(72,rec=1)gq1        ! Earlier in time
          read(72,rec=2)gq2        ! Later in time
        endif
#ifdef USEMPI
!
! wwvv todo
! would be cleaner if we had defined general subroutines 
! for this in spaceparams
!
        call space_distribute("y",sl,gee1,ee1)
        call space_distribute("y",sl,gee2,ee2)
        call space_distribute("y",sl,gq1,q1)
        call space_distribute("y",sl,gq2,q2)
#else
        ee1=gee1
        ee2=gee2
        q1=gq1
        q2=gq2
#endif
        old=floor((par%t/dtbcfile)+1)
		recpos=1
    end if 

    new=floor((par%t/dtbcfile)+1)

    ! Check for next level in boundary condition file
    if (new/=old) then
        recpos=recpos+(new-old)
		! Check for how many bcfile steps are jumped
        if (new-old>1) then  ! Many steps further in the bc file
           if(xmaster) then
              read(72,rec=recpos+1)gq2
              read(71,rec=recpos+1)gee2
              read(72,rec=recpos)gq1
              read(71,rec=recpos)gee1
           endif

#ifdef USEMPI
           call space_distribute("y",sl,gee2,ee2)
           call space_distribute("y",sl,gq2,q2)
		   call space_distribute("y",sl,gee1,ee1)
           call space_distribute("y",sl,gq1,q1)
#else
           ee1=gee1
		   ee2=gee2
		   q1=gq2
           q2=gq2
#endif 
		else  ! Only one step further in the bc file
		   ee1=ee2
		   q1=q2
		   if(xmaster) then
		      read(72,rec=recpos+1)gq2
              read(71,rec=recpos+1)gee2
		   endif
#ifdef USEMPI
           call space_distribute("y",sl,gee2,ee2)
           call space_distribute("y",sl,gq2,q2)
#else  
           ee2=gee2
           q2=gq2
#endif
        endif
        old=new
    end if
    ht=s%zs0(1:2,:)-zb(1:2,:)
    tnew = dble(new)*dtbcfile
    ee(1,:,:) = (dtbcfile-(tnew-par%t))/dtbcfile*ee2 + & !Jaap
                (tnew-par%t)/dtbcfile*ee1
    q = (dtbcfile-(tnew-par%t))/dtbcfile*q2 + &          !Jaap
        (tnew-par%t)/dtbcfile*q1
    ui(1,:) = q/ht(1,:)*min(par%t/par%taper,1.0d0)
    ee(1,:,:)=ee(1,:,:)*min(par%t/par%taper,1.0d0)
else
   if (xmaster) then
     write(*,*)' instat = ',par%instat, ' invalid option'
   endif
   call halt_program
endif
! Jaap: set incoming short wave energy to zero
ee(1,:,:) = par%swave*ee(1,:,:)
! Jaap set incoming long waves to zero
ui = par%lwave*(par%order-1.d0)*ui
! wwvv need to communicate ui here, it is used later on in this
! subroutine
#ifdef USEMPI
call xmpi_shift(ui,'1:')
! wwvv also fill in ee(1,:,:)
call xmpi_shift(ee,'1:')
#endif
if (par%t>0.0d0) then
    if (par%rightwave==0) then
		
		!
		! Lateral boundary at y=0;
		!
		if(par%instat>0 .and. par%instat/=40) then
			fac1=(y(2:nx+1,2)-y(2:nx+1,1))*abs(tan(thetamean(2:nx+1,2)))/(x(2:nx+1,2)-x(1:nx,2))
		else
			fac1=0
		end if
		! fac1=min(fac1,1.d0)
		! fac1=max(fac1,0.d0)\
		! the above approach causes a shoreline jet from outside to inside the domain
		! and mass gain in the domain, so it is turned off here. This means we extrapolate ee laterally in y
		! so
		fac1 = 0.d0 !Ap 28/11   
		fac2=1.d0-fac1
		do itheta=1,ntheta 

			  ee(2:nx+1,1,itheta)=ee(1:nx+1-1,2,itheta)*fac1+ee(2:nx+1,2,itheta)*fac2

		end do
	elseif (par%rightwave==1) then
		wcrestpos=xz+tan(thetamean(:,2))*(yz(2)-yz(1))
		do itheta=1,ntheta
		   do i=1,nx+1
		      call Linear_interp((/-huge(0.d0)   ,wcrestpos     ,huge(0.d0)/),&
			                     (/ee(1,2,itheta),ee(:,2,itheta),ee(nx+1,2,itheta)/),&
								  nx+1,xz(i),ee(i,1,itheta),dummy)
		   enddo
	    enddo
	endif
	
	if (par%leftwave==0) then
		!
		! lateral; boundary at y=ny*dy
		!
		if(par%instat>0 .and. par%instat/=40) then
			fac1=(y(2:nx+1,ny+1)-y(2:nx+1,ny))*abs(tan(thetamean(2:nx+1,ny)))/(x(2:nx+1,ny)-x(1:nx,ny))
		else
			fac1=0.0d0
		end if
		! fac1=min(fac1,1.d0)
		! fac1=max(fac1,0.d0)
		! the above approach causes a shoreline jet from outside to inside the domain
		! and mass gain in the domain, so it is turned off here. This means we extrapolate ee laterally in y
		! so   
		fac1=0.d0 !Ap 28/11

		fac2=1.d0-fac1
		 do itheta=1,ntheta

			   ee(2:nx+1,ny+1,itheta)= &
			   ee(1:nx+1-1,ny+1-1,itheta)*fac1+ee(2:nx+1,ny+1-1,itheta)*fac2

		 end do
	elseif (par%leftwave==1) then
		wcrestpos=xz-tan(thetamean(:,ny))*(yz(ny+1)-yz(ny))
		do itheta=1,ntheta
		   do i=1,nx+1
		      call Linear_interp((/-huge(0.d0)   ,wcrestpos     ,huge(0.d0)/),&
			                     (/ee(1,ny,itheta),ee(:,ny,itheta),ee(nx+1,ny,itheta)/),&
								  nx+1,xz(i),ee(i,ny+1,itheta),dummy)
		   enddo
	    enddo
	 endif
end if
! wwvv communicate ee(:,1,:)
#ifdef USEMPI
call xmpi_shift(ee,':1')
! wwvv and ee(:,ny+1,:)
call xmpi_shift(ee,':n')
! wwv and ee(1,:,:) again
call xmpi_shift(ee,'1:')
#endif

end subroutine wave_bc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! FLOW BOUNDARY CONDITIONS
!
subroutine flow_bc(s,par)
use params
use spaceparams
use interp
use xmpi_module

IMPLICIT NONE

type(spacepars), target                     :: s
type(parameters)                            :: par

integer                                     :: i,ig
integer                                     :: j,jj,indt
real*8                                      :: qxr,alphanew,vert,factime
real*8 , dimension(2)                       :: yzs0,szs0
real*8 , dimension(:,:)  ,allocatable,save  :: zs0old
real*8 , dimension(:,:)  ,allocatable,save  :: ht,beta,betanp1
real*8 , dimension(:)    ,allocatable,save  :: bn,alpha2
real*8 , dimension(:,:)  ,allocatable,save  :: dhdx,dhdy,dvdx,dvdy
real*8 , dimension(:,:)  ,allocatable,save  :: dbetadx,dbetady,dvudy
real*8 , dimension(:)    ,allocatable,save  :: inv_ht

include 's.ind'
include 's.inp'

if (.not. allocated(ht)) then
   allocate(ht      (2,ny+1))
   allocate(dhdx    (2,ny+1))
   allocate(dhdy    (2,ny+1))   ! wwvv not used
   allocate(dvdx    (2,ny+1))   ! wwvv not used
   allocate(dvdy    (2,ny+1))
   allocate(dvudy   (1,ny+1))   ! wwvv not used
   allocate(inv_ht  (ny+1))
   allocate(dbetadx (2,ny+1))
   allocate(dbetady (2,ny+1))
   allocate(beta    (2,ny+1))
   allocate(bn      (ny+1))
   allocate(alpha2  (ny+1))
   allocate(betanp1 (1,ny+1))
   allocate(zs0old(nx+1,ny+1))   ! wwvv not used
endif

factime=1.d0/par%cats/par%Trep*par%dt

!allocate(szs0(1:2))  ! wwvv changed this, now defined as an automatic array
!allocate(yzs0(1:2)) 
!
! Sea boundary at x=0;
!
!
! UPDATE TIDE AND SURGE
! 
! Need to interpolate input tidal signal to xbeach par%t to 
! compute proper tide contribution 
!write(*,*)  'made it into flow_bc and about to interp tide(t)'

if (par%tideloc>0) then

  call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,1), par%tidelen, par%t, par%zs01, indt)

  if(par%tideloc.eq.1) then 
    par%zs02=par%zs01
  end if

  if(par%tideloc.eq.2 .and. par%paulrevere.eq.0) then
    call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,2), par%tidelen, par%t, par%zs03, indt)
    par%zs02=par%zs01
    par%zs04=par%zs03
  endif

  if(par%tideloc.eq.2 .and. par%paulrevere.eq.1) then
    call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,2), par%tidelen, par%t, par%zs02, indt)
    par%zs03=0.d0
    par%zs04=0.d0
  endif

  if(par%tideloc.eq.4) then
    call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,2), par%tidelen, par%t, par%zs02, indt)
    call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,3), par%tidelen, par%t, par%zs03, indt)
    call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,4), par%tidelen, par%t, par%zs04, indt)
  endif

  if(par%tideloc.eq.1) s%zs0 = par%zs01

  if(par%tideloc.eq.2 .and. par%paulrevere.eq.1) then
    yzs0(1)=s%yz(1)
    yzs0(2)=s%yz(s%ny+1)
    szs0(1)=par%zs01
    szs0(2)=par%zs02
    do i = 1,s%ny+1
        call LINEAR_INTERP(yzs0, szs0, 2, s%yz(i), s%zs0(1,i), indt)
    enddo
    do j = 1,s%ny+1 
        do i = 1,s%nx+1
            s%zs0(i,j) = s%zs0(1,j)
        enddo
    enddo   
  endif

  if(par%tideloc.eq.2 .and. par%paulrevere.eq.0) then
    yzs0(1)=s%xz(1)
    yzs0(2)=s%xz(s%nx+1)
    szs0(1)=par%zs01
    szs0(2)=par%zs04
    s%zs0(1,:)=par%zs01
    s%zs0(s%nx+1,:)=par%zs03
! 
! wwvv the following j-loop needs special care in the parallel case
!      zs0(:,j) will get equal to one of
!      -> zs0(1,j)
!      -> zs0(nx+1)
!      -> some interpolation function
!      -> not changed
!      The first two cases require, that the 'real' zs0(1,:)
!      and zs0(nx+1,:) are present, otherwise the local
!      values would be used.
!      We fix this by first copying the first row of the
!      top matrices to the first row of the local matrix
!      and the last row of the bottom matrices to the last
!      row of the local matrix.
!      Then new values of zs0 are determined and afterwards
!      the first and last rows are copied from top and lower neighbours

#ifdef USEMPI
    call xmpi_getrow(s%zs0,ny+1,'1',1,s%zs0(1,:))
    call xmpi_getrow(s%zs0,ny+1,'m',1,s%zs0(nx+1,:))
#endif
    do j = 1,s%ny+1 
        do i = 1,s%nx+1
            if (s%zb(i,j).gt.s%zs0(1,j)+par%eps) goto 302 ! wwvv can use exit here
            s%zs0(i,j) = s%zs0(1,j)
        enddo
302     if (i.lt.s%nx+1) then
            do ig=i+1,s%nx+1      ! wwvv then s%zs0(i) is untouched, ok?
                s%zs0(ig,j) = s%zs0(s%nx+1,j)   
            enddo
        else
            do i = 1,s%nx+1
                call LINEAR_INTERP(yzs0, szs0, 2, s%xz(i), s%zs0(i,j), indt)
            enddo
        endif
    enddo
  endif

  if(par%tideloc.eq.4) then
    yzs0(1)=s%yz(1)
    yzs0(2)=s%yz(s%ny+1)
    szs0(1)=par%zs01
    szs0(2)=par%zs02
    do i = 1,s%ny+1
        call LINEAR_INTERP(yzs0, szs0, 2, s%yz(i), s%zs0(1,i), indt)
    enddo
    yzs0(1)=s%yz(1)
    yzs0(2)=s%yz(s%ny+1)
    szs0(1)=par%zs04
    szs0(2)=par%zs03
    do i = 1,s%ny+1
        call LINEAR_INTERP(yzs0, szs0, 2, s%yz(i), s%zs0(s%nx+1,i), indt)
    enddo

! wwvv see above
#ifdef USEMPI
    call xmpi_getrow(s%zs0,ny+1,'1',1,s%zs0(1,:))
    call xmpi_getrow(s%zs0,ny+1,'m',1,s%zs0(nx+1,:))
#endif
    do j = 1,s%ny+1 
        do i = 1,s%nx+1
            if (s%zb(i,j).gt.s%zs0(1,j)+par%eps) goto 303
            s%zs0(i,j) = s%zs0(1,j)
        enddo
303     if (i.lt.s%nx+1) then
            do ig=i+1,s%nx+1
                s%zs0(ig,j) = s%zs0(s%nx+1,j)   
            enddo
        else
            do i = 1,s%nx+1
                call LINEAR_INTERP(yzs0, szs0, 2, s%xz(i), s%zs0(i,j), indt)
            enddo
        endif
    enddo
  endif
! wwvv in the parallel case, we need to get valid values for the top
!  and bottom row 
#ifdef USEMPI
  call xmpi_shift(s%zs0,'1:')
  call xmpi_shift(s%zs0,'m:')
#endif

!endif
else ! ie if tideloc=0
s%zs0 = par%zs01
endif
!
! UPDATE (LONG) WAVES
!
if (par%instat/=8)then
! wwvv the following is probably only to do in the top processes, but take care for
! the mpi_shift calls in horizontal directions
  if(xmpi_istop) then
    if (par%front==0) then ! Ad's radiating boundary
       uu(1,:)=2*ui(1,:)-(sqrt(par%g/hh(1,:))*(zs(2,:)-s%zs0(2,:)))
       vv(1,:)=vv(2,:)
    elseif (par%front==1) then ! Van Dongeren (1997), weakly reflective boundary condition
       ht(1:2,:)=max(s%zs0(1:2,:)-zb(1:2,:),par%eps)
       beta=uu(1:2,:)-2.*dsqrt(par%g*hum(1:2,:)) !cjaap : replace hh with hum

       do j=2,ny
            ! compute gradients in u-points....
            dvdy(1,j)=(vu(1,min(j,ny)+1)-vu(1,max(j,2)-1))/(yz(min(j,ny)+1)-yz(max(j,2)-1))
            dhdx(1,j)=(ht(2,j)-ht(1,j))/(xz(2)-xz(1))
            dbetadx(1,j)=(beta(2,j)-beta(1,j))/(xu(2)-xu(1))
            dbetady(1,j)=(beta(1,min(j,ny)+1)-beta(1,max(j,2)-1))/(yz(min(j,ny)+1)-yz(max(j,2)-1))

            inv_ht(j) = 1.d0/hum(1,j)                                  !Jaap replaced hh with hum

            bn(j)=-(uu(1,j)-dsqrt(par%g*hum(1,j)))*dbetadx(1,j) &      !Jaap replaced hh with hum 
                  -vu(1,j)*dbetady(1,j)& !Ap vu
                  +dsqrt(par%g*hum(1,j))*dvdy(1,j)&                    !Jaap replaced hh with hum 
                  +Fx(1,j)*inv_ht(j)/par%rho-par%g/par%C**2.d0&        !Ap
                  *sqrt(uu(1,j)**2+vu(1,j)**2)*uu(1,j)/hum(1,j)&       !Jaap replaced hh with hum
                  +par%g*dhdx(1,j)
       end do

       do j=2,ny
          betanp1(1,j) = beta(1,j)+ bn(j)*par%dt
          alpha2(j)=-theta0
          alphanew = 0.d0
          s%umean(1,j) = (factime*uu(1,j)+(1-factime)*s%umean(1,j))  
          do jj=1,50
             !---------- Lower order bound. cond. ---
             qxr = dcos(alpha2(j))/(dcos(alpha2(j))+1.d0)&
               *(0.5d0*(ht(1,j)+ht(2,j))*(betanp1(1,j)-s%umean(1,j)+2.d0*DSQRT(par%g*0.5d0*(ht(1,j)+ht(2,j))))&  !Jaap replaced ht(1,j) with 0.5*(ht(1,j)+ht(2,j))
               -(ui(1,j)*hum(1,j))*(dcos(theta0)-1.d0)/dcos(theta0))   !Jaap replaced hh with hu

             !vert = velocity of the reflected wave = total-specified
             vert = vu(1,j)-ui(1,j)*tan(theta0)
             alphanew = datan(vert*hh(1,j)/(qxr+1.d-16))                   ! wwvv can  use atan2 here
             if (alphanew .gt. (par%px*0.5d0)) alphanew=alphanew-par%px
             if (alphanew .le. (-par%px*0.5d0)) alphanew=alphanew+par%px
             if(dabs(alphanew-alpha2(j)).lt.0.001d0) goto 1000     ! wwvv can use exit here
             alpha2(j) = alphanew 
          end do
    1000  continue
          if (par%ARC==0) then
             uu(1,j) = (par%order-1)*ui(1,j)
             zs(1,j) = zs(2,j)
          else
		  	 ! jaap: to fix surge in 2DH case
	         ! s%umean(1,j) = (par%epsi*2.d0*qxr/(ht(1,j)+ht(2,j))+(1-par%epsi)*s%umean(1,j))  
			 !
             uu(1,j) = (par%order-1.d0)*ui(1,j)+2.d0*qxr/(ht(1,j)+ht(2,j)) + s%umean(1,j)
             !with a taylor expansion to get to the zs point at index 1 from uu(1) and uu(2)
             zs(1,j) = 1.5d0*((betanp1(1,j)-uu(1,j))**2/4.d0/par%g+.5d0*(zb(1,j)+zb(2,j)))- &
                       0.5d0*((beta(2,j)-uu(2,j))**2/4.d0/par%g+.5d0*(zb(2,j)+zb(3,j)))
					 ! Ad + Jaap: zs does indeed influence hydrodynamics at boundary --> do higher order taylor expansions to check influence
		     ! zs(1,j) = 13.d0/8.d0*((betanp1(1,j)-uu(1,j))**2.d0/4.d0/par%g+.5d0*(zb(1,j)+zb(2,j))) - &
		     !           0.75d0*((beta(2,j)-uu(2,j))**2.d0/4.d0/par%g+.5d0*(zb(2,j)+zb(3,j)))        + &
		     !		     0.125d0*0.5d0*(zs(3,j)+zs(4,j))
    
          end if
       end do
       vv(1,:)=vv(2,:)
    endif ! par%front
    ! uu, zs and umean shift horizontally in two directions (loop was 2..ny)
#ifdef USEMPI
    call xmpi_shift(s%umean,':1')
    call xmpi_shift(s%umean,':n')
    call xmpi_shift(uu,':1')
    call xmpi_shift(uu,':n')
    call xmpi_shift(zs,':1')
    call xmpi_shift(zs,':n')
#endif
  endif  ! xmpi_istop
  ! wwvv uu and zs and vv need to be communicated vertically 
#ifdef USEMPI
  call xmpi_shift(uu,'1:')
  call xmpi_shift(zs,'1:')
  call xmpi_shift(vv,'1:')
#endif
  !
  ! Radiating boundary at x=nx*dx
  !
  if (xmpi_isbot) then
     if (par%back==0) then ! set uu(nx+1,:)=0 
        uu(nx,:) = 0.d0   
        zs(nx+1,:) = zs(nx,:)
		! zs(nx+1,2:ny) = zs(nx+1,2:ny) + par%dt*hh(nx,2:ny)*uu(nx,2:ny)/(xu(nx+1)-xu(nx)) -par%dt*(hv(nx+1,2:ny+1)*vv(nx+1,2:ny+1)-hv(nx+1,1:ny)*vv(nx+1,1:ny))/(yv(2:ny+1)-yv(1:ny))
     elseif (par%back==1) then
        ! uu(nx+1,:)=sqrt(par%g/hh(nx+1,:))*(zs(nx+1,:)-max(zb(nx+1,:),s%zs0(nx+1,:))) ! cjaap: make sure if the last cell is dry no radiating flow is computed... 
        ! uu(nx+1,:)=sqrt(par%g/(s%zs0(nx+1,:)-zb(nx+1,:)))*(zs(nx+1,:)-max(zb(nx+1,:),s%zs0(nx+1,:)))
        s%umean(2,:) = factime*uu(nx,:)+(1-factime)*s%umean(2,:)    !Ap
        zs(nx+1,:)=max(s%zs0(nx+1,:),s%zb(nx+1,:))+(uu(nx,:)-s%umean(2,:))*sqrt(max((s%zs0(nx+1,:)-zb(nx+1,:)),par%eps)/par%g)    !Ap
     elseif (par%back==2) then
	   ht(1:2,:)=max(s%zs0(s%nx:s%nx+1,:)-zb(s%nx:s%nx+1,:),par%eps) !cjaap; make sure ht is always larger than zero

       beta=uu(s%nx-1:s%nx,:)+2.*dsqrt(par%g*hum(s%nx-1:s%nx,:)) !cjaap : replace hh with hum

       do j=2,ny
          if (wetu(s%nx,j)==1) then   ! Robert: dry back boundary points
            ! Compute gradients in u-points....
            dvdy(2,j)=(vu(s%nx,min(j,ny)+1)-vu(s%nx,max(j,2)-1))/(yz(min(j,ny)+1)-yz(max(j,2)-1))
            dhdx(2,j)=(ht(2,j)-ht(1,j))/(xz(s%nx+1)-xz(s%nx))
            dbetadx(2,j)=(beta(2,j)-beta(1,j))/(xu(s%nx)-xu(s%nx-1))
            dbetady(2,j)=(beta(1,min(j,ny)+1)-beta(1,max(j,2)-1))/(yz(min(j,ny)+1)-yz(max(j,2)-1))

            inv_ht(j) = 1.d0/hum(s%nx,j)                                           !Jaap replaced hh with hum

            bn(j)=-(uu(s%nx,j)+dsqrt(par%g*hum(s%nx,j)))*dbetadx(2,j) &            !Ap says plus  !Jaap replaced hh with hum 
                  -vu(s%nx,j)*dbetady(2,j)& !Ap vu
                  -dsqrt(par%g*hum(s%nx,j))*dvdy(2,j)&                             !Jaap replaced hh with hum 
                  +Fx(s%nx,j)*inv_ht(j)/par%rho-par%g/par%C**2.d0&                 !Ap
                  *sqrt(uu(s%nx,j)**2+vu(s%nx,j)**2)*uu(s%nx,j)/hum(s%nx,j)&       !Jaap replaced hh with hum
                  +par%g*dhdx(2,j)
          endif   ! Robert: dry back boundary points
       enddo

       do j=2,ny
         if (wetu(s%nx,j)==1) then                                                   ! Robert: dry back boundary points
          betanp1(1,j) = beta(2,j)+ bn(j)*par%dt                                   !Ap toch?
          alpha2(j)= theta0
          alphanew = 0.d0
          s%umean(2,j) = (factime*uu(s%nx,j)+(1-factime)*s%umean(2,j))           !Ap 
          do jj=1,50
             !---------- Lower order bound. cond. ---
             qxr = dcos(alpha2(j))/(dcos(alpha2(j))+1.d0)&  
               *(0.5*(ht(1,j)+ht(2,j))*(betanp1(1,j)-s%umean(2,j)-2.d0*DSQRT(par%g*0.5*(ht(1,j)+ht(2,j))))) 

             !vert = velocity of the reflected wave = total-specified
             vert = vu(s%nx,j)
             alphanew = datan(vert*hh(s%nx+1,j)/(qxr+1.d-16))                      !Ap  ! wwvv maybe better atan2
             if (alphanew .gt. (par%px*0.5d0)) alphanew=alphanew-par%px
             if (alphanew .le. (-par%px*0.5d0)) alphanew=alphanew+par%px
             if(dabs(alphanew-alpha2(j)).lt.0.001) goto 2000    ! wwvv can use exit here
             alpha2(j) = alphanew 
          end do
    2000  continue
          uu(s%nx,j) = 2.*qxr/(ht(1,j)+ht(2,j)) + s%umean(2,j)                       !Jaap: replaced ht(1,j) with 0.5*(ht(1,j)+ht(2,j))
              ! Ap replaced zs with extrapolation.
          zs(s%nx+1,j) = 1.5*((betanp1(1,j)-uu(s%nx,j))**2.d0/4.d0/par%g+.5*(zb(s%nx,j)+zb(s%nx+1,j)))-&
                                         0.5*((beta(1,j)-uu(s%nx-1,j))**2.d0/4.d0/par%g+.5*(zb(s%nx-1,j)+zb(s%nx,j)))
      
         endif   ! Robert: dry back boundary points
       enddo
     endif  !par%back
    ! fix first and last columns of s%umean and uu and zs
#ifdef USEMPI
    call xmpi_shift(s%umean,':1')
    call xmpi_shift(s%umean,':n')
    call xmpi_shift(uu,':1')
    call xmpi_shift(uu,':n')
    call xmpi_shift(zs,':1')
    call xmpi_shift(zs,':n')
#endif
endif !xmpi_istop
#ifdef USEMPI
  call xmpi_shift(uu,'m:')
  call xmpi_shift(zs,'m:')
  call xmpi_shift(vv,'m:')
#endif

endif   ! par%instat


!!! Wind boundary conditions

call LINEAR_INTERP(s%windinpt,s%windvel,par%windlen,par%t,s%windvnow,indt)
call LINEAR_INTERP(s%windinpt,s%winddir,par%windlen,par%t,s%winddirnow,indt)

end subroutine flow_bc

end module boundaryconditions
