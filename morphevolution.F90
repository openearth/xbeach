module morphevolution
contains
subroutine transus(s,par)
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
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars),target                   :: s
type(parameters)                         :: par

integer                                  :: i
integer                                  :: j,jg

real*8,dimension(:,:),allocatable,save   :: vmag2
real*8 , dimension(s%nx+1,s%ny+1)   :: dzbx,dzby
real*8,dimension(:,:),allocatable,save   :: cc,cu,cv,Su,Sv,Dc,termh

! wwvv removed s.ind and s.inc and put s% in appropriate locations
!      to be able to better analyse the code

if (.not. allocated(vmag2)) then
   allocate(vmag2 (s%nx+1,s%ny+1))
   allocate(cu (s%nx+1,s%ny+1))
   allocate(cv (s%nx+1,s%ny+1))
   allocate(cc (s%nx+1,s%ny+1))
   allocate(Su (s%nx+1,s%ny+1))
   allocate(Sv (s%nx+1,s%ny+1))
   allocate(Dc (s%nx+1,s%ny+1))
   allocate(termh (s%nx+1,s%ny+1))
endif
! use eulerian velocities
vmag2    = s%ue**2+s%ve**2
s%dcdx     = 0.0d0
s%dcdy     = 0.0d0
! calculate equilibrium concentration
if (par%form==1) then           ! soulsby van Rijn
        call sb_vr(s,par)
elseif (par%form==2) then       ! Van Rijn 2008
    call vr2008(s,par)
elseif (par%form==3) then       ! Van Thiel de Vries & Reniers 2008
    call sednew(s,par)
end if

! compute diffusion coefficient

if (par%nuhfac==1) then
   termh = s%hh/max(s%H,.01d0)
   ! termh = max(s%hh(i,j)/2.d0,0.01d0)
   termh = min(termh,10.d0);
   ! Dc = 5*(DR/par%rho)**(1.d0/3.d0)*s%hh/(exp(termh)-1.d0)
   Dc = par%nuh+par%nuhfac*s%hh*(s%DR/par%rho)**(1.d0/3.d0)
   ! Dc = par%dico
else
   Dc = par%dico
end if

do jg = 1,par%ngd
   cc = s%ccg(:,:,jg)
!   cc = ceqg(:,:,jg) ! Can be used to test total transport mode
   !
   ! X-direction
   do j=1,s%ny+1
      do i=1,s%nx
         if(s%ueu(i,j)>0.d0) then
   !test          cu(i,j)=cc(i,j)
           cu(i,j)=par%thetanum*cc(i,j)+(1.d0-par%thetanum)*cc(min(i+1,s%nx),j)
         elseif(s%ueu(i,j)<0.d0) then
            cu(i,j)=par%thetanum*cc(i+1,j)+(1.d0-par%thetanum)*cc(max(i,2),j)
         else
            cu(i,j)=0.5d0*(cc(i,j)+cc(i+1,j))
         endif
         s%dcdx(i,j)=(cc(i+1,j)-cc(i,j))/(s%xz(i+1)-s%xz(i)) !Jaap
      enddo 
   enddo
   ! wwvv dcdx(nx:1,:) is still untouched, correct this ofr the parallel case
#ifdef USEMPI
   call xmpi_shift(s%dcdx,'m:')
#endif
   cu(s%nx+1,:) = cc(s%nx+1,:) !Robert
   ! wwvv fix this in parallel case
#ifdef USEMPI
   call xmpi_shift(cu,'m:')
#endif
   !
   ! Bed slope terms
   !
   dzbx=0.d0
   do j=1,s%ny+1
       do i=1,s%nx
           dzbx(i,j)=(s%zb(i+1,j)-s%zb(i,j))/(s%xz(i+1)-s%xz(i))
       enddo
   enddo
   ! wwvv in parallel version, there will be a discrepancy between the values
   ! of dzbx(nx+1,:).
   !wwvv so fix that
#ifdef USEMPI
   call xmpi_shift(dzbx,'m:')
#endif
   Su=(cu*s%ueu*s%hu-Dc*s%hu*s%dcdx-par%facsl*cu*s%vmagu*s%hu*dzbx)*s%wetu   !
   !
   ! Y-direction
   !
   do j=1,s%ny
      do i=1,s%nx+1
         if(s%vev(i,j)>0) then
        !   cv(i,j)=cc(i,j)
            cv(i,j)=par%thetanum*cc(i,j)+(1.d0-par%thetanum)*cc(i,min(j+1,s%ny))
          else if(s%vev(i,j)<0) then
        !   cv(i,j)=cc(i,j+1)
            cv(i,j)=par%thetanum*cc(i,j+1)+(1.d0-par%thetanum)*cc(i,max(j,2))
        else
         cv(i,j)=0.5d0*(cv(i,j)+cv(i,j+1))
         end if
         s%dcdy(i,j)=(cc(i,j+1)-cc(i,j))/(s%yz(j+1)-s%yz(j)) !Jaap
      end do 
    end do
    ! wwvv s%dcdy(:,ny+1) is not filled in, so in parallel case:
#ifdef USEMPI
    call xmpi_shift(s%dcdy,':n')
#endif
    cv(:,s%ny+1) = cc(:,s%ny+1) !Robert
    !
    ! Bed slope terms
    !
    dzby=0.d0
    do j=1,s%ny
        do i=1,s%nx+1
            dzby(i,j)=(s%zb(i,j+1)-s%zb(i,j))/(s%yz(j+1)-s%yz(j))
        enddo
    enddo
   ! wwvv in parallel version, there will be a discrepancy between the values
   ! of dzby(:,ny+1).
    ! wwvv so fix that
#ifdef USEMPI
    call xmpi_shift(dzby,':n')
#endif
    Sv=(cv*s%vev*s%hv-Dc*s%hv*s%dcdy-par%facsl*cv*s%vmagv*s%hv*dzby)*s%wetv
    
    do j=2,s%ny+1
       do i=2,s%nx+1
          cc(i,j) = s%hold(i,j)*cc(i,j)-                               &
            par%dt*((Su(i,j)-Su(i-1,j))/(s%xu(i)-s%xu(i-1))+             &
            (Sv(i,j)-Sv(i,j-1))/(s%yv(j)-s%yv(j-1))-                     &
            s%hold(i,j)*(s%ceqg(i,j,jg)*s%graindistr(i,j,1,jg)-cc(i,j))/ &
            s%Tsg(i,j,jg)) ! Jaap: use Tsg instead of Ts

          cc(i,j)=max(cc(i,j),0.0d0)
       enddo
    enddo

    do j=1,s%ny+1
       do i=1,s%nx+1
          if(s%hh(i,j)>=par%hmin) then 
             cc(i,j)=cc(i,j)/s%hh(i,j)  
          else
             cc(i,j)=0.d0
          end if
       end do
    end do
    cc(1,:)=cc(2,:)
    cc(:,1)=cc(:,2)
    cc(s%nx+1,:)=cc(s%nx+1-1,:)
    cc(:,s%ny+1)=cc(:,s%ny+1-1)
    ! wwvv fix the first and last rows and columns of cc in parallel case
#ifdef USEMPI
    call xmpi_shift(cc,'1:')
    call xmpi_shift(cc,'m:')
    call xmpi_shift(cc,':1')
    call xmpi_shift(cc,':n')
#endif
    !Jaap
    cc=cc*s%wetz
    !
    ! wwvv border columns and rows of ccg Svg and Sug have to be communicated
    s%ccg(:,:,jg) = cc
    s%Svg(:,:,jg) = Sv
    s%Sug(:,:,jg) = Su
end do

s%vmag=sqrt(max(vmag2,par%umin))

end subroutine transus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bed_update(s,par)
use params
use spaceparams

IMPLICIT NONE

type(spacepars),target              :: s
type(parameters)                    :: par

integer                             :: i,ii
integer                             :: j,jg,jd
real*8                              :: dzb,dzmax
real*8 , dimension(s%nx+1,s%ny+1)   :: dzbx,dzby,Su,Sv 
real*8 , dimension(s%nx+1,s%ny+1)   :: dzbtot,fact0,fact1,fact2,fact3,fact4,fact5
real*8 , dimension(s%nx+1,s%ny+1,par%ngd) :: dzbdtg
real*8 , dimension(s%nx+1,s%ny+1,par%nd)  :: graindistrm
logical                                   :: aval

! wwvv removed s.ind and s.inc and put s% in appropriate locations
!      to be able to better analyse the code


dzbdtg = 0.d0
dzbtot = 0.d0

! Update bed level using continuity eq.
dzb=0.d0 !!!Ap
s%dzbdt=0.0d0

if (par%t>=par%morstart) then

do jg = 1,par%ngd
   ! Update bed level using continuity eq.
   Su = s%Sug(:,:,jg)
   Sv = s%Svg(:,:,jg)

   ! Ap and Jaap: Bottom update is morfac dependent....
   do j=2,s%ny !Jaap nx instead of s%ny+1
      do i=2,s%nx  !Jaap s%nx instead of s%nx+1
         s%dzbdt(i,j)=1/(1-par%por)*        &
             (-(Su(i,j)-Su(i-1,j))/(s%xu(i)-s%xu(i-1)) &
              -(Sv(i,j)-Sv(i,j-1))/(s%yv(j)-s%yv(j-1)))
         s%zb(i,j)=s%zb(i,j)+s%dzbdt(i,j)*par%morfac*par%dt
         s%sedero(i,j)=s%sedero(i,j)+s%dzbdt(i,j)*par%morfac*par%dt           ! total bed level changes for all sediment clases over all time steps so far
         dzbdtg(i,j,jg) = dzbdtg(i,j,jg) + s%dzbdt(i,j)*par%morfac*par%dt ! total bed level change for each sediment class at this time step
         dzbtot(i,j) = dzbtot(i,j) +  s%dzbdt(i,j)*par%morfac*par%dt      ! total bed level change for all sediment classes at this time step
       enddo
   enddo
enddo

fact0 = 0.d0
fact1 = 0.d0
fact2 = 0.d0
fact3 = 0.d0
fact4 = 0.d0
fact5 = 0.d0

! update sediment fractions in sediment layers
! cjaap
if (par%ngd>1) then

do j=2,s%ny 
   do i=2,s%nx  
      if (dzbtot(i,j)<0.d0) then
         fact0(i,j) = 0.d0      
         fact1(i,j) = 1.d0
         fact2(i,j) = -dzbtot(i,j)/par%dzg
         fact3(i,j) = 0.d0
         fact4(i,j) = (par%dzg+dzbtot(i,j))/par%dzg
         fact5(i,j) = -dzbtot(i,j)/par%dzg
      else
         fact0(i,j) = 1.d0  ! only in case of vegetated layer 
         fact1(i,j) = (par%dzg-dzbtot(i,j))/par%dzg
         fact2(i,j) = 0.d0
         fact3(i,j) = dzbtot(i,j)/par%dzg
         fact4(i,j) = (par%dzg-dzbtot(i,j))/par%dzg
         fact5(i,j) = 0.d0
      endif
   enddo
enddo

do jd = 1,par%nd-1
   do jg = 1,par%ngd
      if (jd == 1) then
         if (jg == 1) then
            s%graindistr(:,:,jd,jg) = dzbdtg(:,:,jg)/par%dzg + fact1*s%graindistr(:,:,jd,jg) + &
            fact2*s%graindistr(:,:,jd+1,jg) + max(dzbdtg(:,:,jg+1)/par%dzg,0.0d0)
         elseif (jg==2) then
            !s%graindistr(:,:,jd,jg) = dzbdtg(:,:,jg)/par%dzg + fact1*s%graindistr(:,:,jd,jg) + fact2*s%graindistr(:,:,jd+1,jg)
            s%graindistr(:,:,jd,jg) = min(dzbdtg(:,:,jg)/par%dzg,0.0d0) + fact1*s%graindistr(:,:,jd,jg) + &
            fact2*s%graindistr(:,:,jd+1,jg)
         else
            s%graindistr(:,:,jd,jg) = dzbdtg(:,:,jg)/par%dzg + fact1*s%graindistr(:,:,jd,jg) + fact2*s%graindistr(:,:,jd+1,jg) 
         endif
      else
         s%graindistr(:,:,jd,jg) = fact3*s%graindistr(:,:,jd-1,jg) + fact4*s%graindistr(:,:,jd,jg) + fact5*s%graindistr(:,:,jd+1,jg)
      endif
   enddo
enddo

! ensure total fractions equal 1
graindistrm = 0.d0
do jd = 1,par%nd
    do jg = 1,par%ngd
       graindistrm(:,:,jd) = graindistrm(:,:,jd) + s%graindistr(:,:,jd,jg)
    enddo
    do jg = 1,par%ngd
       s%graindistr(:,:,jd,jg) = s%graindistr(:,:,jd,jg)/max(graindistrm(:,:,jd),.0010d0)
    enddo
enddo

endif
!
! bed boundary conditions
! 
! Fix bed at back boundary, but allow sediment to pass though
!s%zb(s%nx+1,:) = s%zb(s%nx+1,:)+s%dzbdt(s%nx,:)*par%dt  !Ap
!s%zb(s%nx+1,:) = s%zb(s%nx+1,:)+s%dzbdt(s%nx,:)*par%dt  !Ap
!s%sedero(s%nx+1,:) = s%sedero(s%nx,:) !Ap
s%zb(:,1) = s%zb(:,2)
s%sedero(:,1) = s%sedero(:,2)
s%zb(:,s%ny+1) = s%zb(:,s%ny)
s%sedero(:,s%ny+1) = s%sedero(:,s%ny)
!
! Avalanching
!
do ii=1,nint(par%morfac)
   aval=.false.
   dzbx=0.d0
   dzby=0.d0
   do j=1,s%ny+1
      do i=1,s%nx
         dzbx(i,j)=(s%zb(i+1,j)-s%zb(i,j))/(s%xz(i+1)-s%xz(i))
      enddo
   enddo
   do j=1,s%ny
      do i=1,s%nx+1
         dzby(i,j)=(s%zb(i,j+1)-s%zb(i,j))/(s%yz(j+1)-s%yz(j))
      enddo
   enddo
   !

   do i=2,s%nx-1
      do j=1,s%ny+1
         if(max(s%hh(i,j),s%hh(i+1,j))>par%hswitch+par%eps) then
            dzmax=par%wetslp;
         else
            dzmax=par%dryslp;
         end if
         if(abs(dzbx(i,j))>dzmax ) then    
                    aval=.true.
            dzb=sign(1.0d0,dzbx(i,j))*(abs(dzbx(i,j))-dzmax)*(s%xz(i+1)-s%xz(i));;
            dzb=sign(1.0d0,dzb)*min(abs(dzb),0.05d0*par%dt); 
            s%zb(i,j)=s%zb(i,j)+dzb*(s%xu(i+1)-s%xu(i))/(s%xu(i)-s%xu(i-1)); !Jaap make sure there is continuity of sediment in non uniform grids;
            s%sedero(i,j) = s%sedero(i,j)+dzb*(s%xu(i+1)-s%xu(i))/(s%xu(i)-s%xu(i-1));
            s%zs(i,j)=s%zs(i,j)+dzb*(s%xu(i+1)-s%xu(i))/(s%xu(i)-s%xu(i-1));
            s%zb(i+1,j)=s%zb(i+1,j)-dzb
            s%sedero(i+1,j)=s%sedero(i+1,j)-dzb
            s%zs(i+1,j)=s%zs(i+1,j)-dzb
                        s%dzav(i+1,j)= s%dzav(i+1,j)-dzb  !Jaap compute total bed level change due to avalanching
                        s%dzav(i,j)= s%dzav(i,j)+dzb
         end if 
      end do
   end do

   do j=2,s%ny-1;
        do i=1,s%nx+1
            if(max(s%hh(i,j),s%hh(i,j+1))>par%hswitch+par%eps) then
                dzmax=par%wetslp
            else
                dzmax=par%dryslp
            end if
            if(abs(dzby(i,j))>dzmax ) then 
                aval=.true. 
                dzb=sign(1.0d0,dzby(i,j))*(abs(dzby(i,j))-dzmax)*(s%yz(j+1)-s%yz(j))       
                dzb=sign(1.0d0,dzb)*min(abs(dzb),0.05d0*par%dt) 
                s%zb(i,j)=s%zb(i,j)+dzb*(s%yv(j+1)-s%yv(j))/(s%yv(j)-s%yv(j-1))
                !Jaap make sure there is continuity of sediment in non uniform grids;;
                s%sedero(i,j)=s%sedero(i,j)+dzb*(s%yv(j+1)-s%yv(j))/(s%yv(j)-s%yv(j-1))
                s%zs(i,j)=s%zs(i,j)+dzb*(s%yv(j+1)-s%yv(j))/(s%yv(j)-s%yv(j-1))
                s%zb(i,j+1)=s%zb(i,j+1)-dzb
                s%sedero(i,j+1)=s%sedero(i,j+1)-dzb
                s%zs(i,j+1)=s%zs(i,j+1)-dzb
                s%dzav(i,j+1)= s%dzav(i,j+1)-dzb  !Jaap compute total bed level change due to avalanching
                s%dzav(i,j)= s%dzav(i,j)+dzb
            end if
        end do
   end do
   if (.not.aval) exit
end do
endif

end subroutine bed_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sb_vr(s,par)
use params
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars),target              :: s
type(parameters)                    :: par

integer                             :: i
integer                             :: j,jg
real*8                              :: dster,twothird
real*8                              :: z0,Ass,delta
real*8                              :: Te,kvis,Sster,c1,c2,wster


real*8 , dimension(:,:),allocatable,save   :: vmag2,Cd,Asb,dhdx,dhdy,Ts
real*8 , dimension(:,:),allocatable,save   :: urms2,Ucr,term1,term2
real*8 , dimension(:,:),allocatable,save   :: uandv,b,fslope,hloc,ceq
real*8 , dimension(:,:),allocatable,save   :: kb ! wwvv: added. After removing 
                                                 ! s.ind and s.inc, it appeared
! that with s.ind and s.inc kb was in spacepars

! wwvv removed s.ind and s.inc and put s% in appropriate locations
!      to be able to better analyse the code


if (.not. allocated(vmag2)) then
   allocate (vmag2 (s%nx+1,s%ny+1))
   allocate (Cd    (s%nx+1,s%ny+1))
   allocate (Asb   (s%nx+1,s%ny+1))
   allocate (dhdx  (s%nx+1,s%ny+1))   ! not used wwvv
   allocate (dhdy  (s%nx+1,s%ny+1))   ! not used wwvv
   allocate (urms2 (s%nx+1,s%ny+1))
   allocate (Ucr   (s%nx+1,s%ny+1))
   allocate (term1 (s%nx+1,s%ny+1))
   allocate (term2 (s%nx+1,s%ny+1))
   allocate (uandv (s%nx+1,s%ny+1))  ! not used wwvv
   allocate (b     (s%nx+1,s%ny+1))  ! not used wwvv
   allocate (fslope(s%nx+1,s%ny+1))  ! not used wwvv
   allocate (hloc  (s%nx+1,s%ny+1))
   allocate (kb    (s%nx+1,s%ny+1))
   allocate (Ts    (s%nx+1,s%ny+1))
   allocate (ceq   (s%nx+1,s%ny+1))
endif
! Soulsby van Rijn sediment transport formula
! Ad Reniers april 2006
!
! z is defined positive upward 
! x is defined positive toward the shore
!  Formal parameters:
!  ------------------
!
!   Var. I/O  Type Dimensions
!   -------------------------
!

hloc   = max(s%hh,par%hmin) !Jaap par%hmin instead of par%eps
twothird=2.d0/3.d0
delta = (par%rhos-par%rho)/par%rho
! use eulerian velocities
! cjaap: add turbulence near bottom
do j=1,s%ny+1 
    do i=1,s%nx+1
       kb(i,j) = (s%DR(i,j)/par%rho)**twothird/ &
         (exp( min( hloc(i,j)/max(s%H(i,j),0.1d0) ,100.d0)) -1.d0)
       vmag2(i,j) = s%ue(i,j)**2+s%ve(i,j)**2
    enddo
enddo

vmag2  = s%ue**2+s%ve**2  ! wwvv todo just to be sure ?
urms2  = s%urms**2+0.50d0*kb

do jg = 1,par%ngd

   ! cjaap: compute fall velocity with simple expression from Ahrens (2000)
   Te    = 20.d0
   kvis  = 4.d0/(20.d0+Te)*1d-5 ! Van rijn, 1993 
   Sster = s%D50(jg)/(4*kvis)*sqrt((par%rhos/par%rho-1)*par%g*s%D50(jg))
   c1    = 1.06d0*tanh(0.064d0*Sster*exp(-7.5d0/Sster**2))
   c2    = 0.22d0*tanh(2.34d0*Sster**(-1.18d0)*exp(-0.0064d0*Sster**2))
   wster = c1+c2*Sster
   par%w = wster*sqrt((par%rhos/par%rho-1.d0)*par%g*s%D50(jg))

   Ts       = par%tsfac*s%hh/par%w
   !Ts       = max(Ts,0.2d0)
   s%Tsg(:,:,jg) = max(Ts,0.2d0) 

   dster=25296*s%D50(jg)
   
   ! calculate treshold velocity Ucr
   if(s%D50(jg)<=0.0005d0) then
     Ucr=0.19d0*s%D50(jg)**0.1d0*log10(4*hloc/s%D90(jg))
   else if(s%D50(jg)<0.002d0) then
     Ucr=8.5d0*s%D50(jg)**0.6d0*log10(4*hloc/s%D90(jg))
   else
#ifdef USEMPI
     write(*,'(a,i4)') 'In process',xmpi_rank
#endif
     write(*,*) '  Remark from sb_vr: s%D50(jg) > 2mm, out of validity range'
   end if
   ! drag coefficient
   z0 = par%z0
   Cd=(0.40d0/(log(max(s%hh,par%hmin)/z0)-1.0d0))**2 !Jaap
   !Cd = par%g/par%C**2;   ! consistent with flow modelling 

   ! Diane Foster and Robert: limit Shields to par%smax -> vmag2 for transp. limited
  ! vmag2=min(vmag2,par%smax*par%C**2*s%D50(jg)*delta)
   vmag2=min(vmag2,par%smax*par%g/par%cf*s%D50(jg)*delta)       ! In terms of cf

   ! transport parameters
   Asb=0.005d0*hloc*(s%D50(jg)/hloc/(delta*par%g*s%D50(jg)))**1.2d0     ! bed load coefficent
   Ass=0.012d0*s%D50(jg)*dster**(-0.6d0)/(delta*par%g*s%D50(jg))**1.2d0 ! suspended load coeffient
   term1=sqrt(vmag2+0.018d0/Cd*urms2)     ! nearbed-velocity
      
   term2 = 0
   do j=1,s%ny+1
      do i=1,s%nx
         if(term1(i,j)>Ucr(i,j) .and. s%hh(i,j)>par%eps) then
            term2(i,j)=(term1(i,j)-Ucr(i,j))**2.4d0
!           term2(i,j)=(term1(i,j)-Ucr(i,j))**(1.7-0.7*tanh((s%ue(i,j)/sqrt(par%g*s%hh(i,j))-0.5)*10))
         end if
      end do
   end do
   ! wwvv in parallel version, there will be a discrepancy between the values
   ! of term2. term2(nx+1,:) is zero, while the corresponding row in the process
   ! below term2(2,:) has some value, different from zero.
   ! so we fix this:
#ifdef USEMPI
   call xmpi_shift(term2,'m:')
#endif
   ceq =(Asb+Ass)*term2  
   ceq = min(ceq,0.2d0)       ! equilibrium concentration
   ceq = ceq/hloc
   s%ceqg(:,:,jg) = ceq*s%sedcal(jg)
enddo  ! end og grain size classes

end subroutine sb_vr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine vr2008(s,par)
use params
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars),target             :: s
type(parameters)                   :: par

integer                            :: i,j,jg
real*8                             :: m1,m2,m3,m4,m5,m6
real*8                             :: bb1,bb2,a1,a2,a3,a4,p1,q1,p2,q2,c1,c2,c3,r3,r4,ksi,phi
real*8                             :: onethird,dster,Te,kvis,Sster,cc1,cc2,wster,Ass
real*8 , dimension(:,:),allocatable,save  :: Ur,Bm,B1,sigu2,amp2,amp1,As,Sk,hloc
real*8 , dimension(:,:),allocatable,save  :: tmin,tmax
real*8 , dimension(:,:),allocatable,save  :: vmg,uorb,Asb,Ts,Ucr,Ucrc,Ucrw,term1,B2,ceq

! wwvv removed s.ind and s.inc and put s% in appropriate locations
!      to be able to better analyse the code
! on the other hand, this code cannot run, because there are no allocations for 
! the allocatable arrays ...    todo

hloc   = max(s%hh,par%hmin)
onethird=1.d0/3.d0

!
! compute uon and uoff, analytical solution Ruessink and Van Rijn
!
m1 = 0;       ! a = 0
m2 = 0.7939;  ! b = 0.79 +/- 0.023
m3 = -0.6065; ! c = -0.61 +/- 0.041
m4 = 0.3539;  ! d = -0.35 +/- 0.032 
m5 = 0.6373;  ! e = 0.64 +/- 0.025
m6 = 0.5995;  ! f = 0.60 +/- 0.043

do j=1,s%ny+1     
   do i=1,s%nx+1
       if (s%k(i,j)*s%h(i,j)<par%px/2.d0 .and. s%H(i,j)>0.01d0) then
              Ur(i,j) = 3.d0/8.d0*sqrt(2.d0)*s%H(i,j)*s%k(i,j)/(s%k(i,j)*hloc(i,j))**3
                                         !Ursell number
          Bm(i,j) = m1+(m2-m1)/(1.d0+exp((m3-log10(Ur(i,j)))/m4))
                                         !Boltzmann sigmoid (eq 6) 
          B1(i,j) = -90.d0+90.d0*tanh(m5/Ur(i,j)**m6)    
                  B1(i,j) = B1(i,j)*par%px/180.d0

          Sk(i,j) = Bm(i,j)*cos(B1(i,j))
                                           !Skewness (eq 8)
          As(i,j) = Bm(i,j)*sin(B1(i,j))
                                           !Skewness (eq 9)

          sigu2(i,j) = (s%H(i,j)*s%sigm(i,j)/(2.d0*sqrt(2.d0)*sinh(s%k(i,j)*hloc(i,j))))**2
              !variance from linear theory (eq 13)

          p1 = -2.d0*sigu2(i,j)
          q1 = -4.d0/3.d0*Sk(i,j)*sigu2(i,j)**1.5d0/cos(B1(i,j))

          ksi = q1/2.d0/sqrt((-p1/3.d0)**3)
          phi = atan2(sqrt(1.d0-ksi**2.d0),ksi)

          amp2(i,j) = 2.d0*sqrt(-p1/3.d0)*cos(phi/3.d0+4.d0/3.d0*par%px)
          amp1(i,j) = sqrt(2.d0*sigu2(i,j)-amp2(i,j)**2)
          !
          ! compute uon and uoff
          !
          bb1 = -amp1(i,j)*s%sigm(i,j)
          bb2 = -2.d0*amp2(i,j)*s%sigm(i,j)
          a1 = bb1/bb2*cos(B1(i,j)/2); 
          a2 = (4.d0*bb2**2-bb1**2)/(-4.d0*bb2**2)
          a3 = -a1
          a4 = bb1**2*cos(B1(i,j)/2)**2/(-4*bb2**2)

          p1 = a2**2-3*a1*a3+12*a4
          q1 = 2*a2**3-9*a1*a2*a3+27*a3**2+27*a4*a1**2-72*a2*a4

          p2 = max(par%eps,q1 + dsqrt(max(0.d0,-4.d0*p1**3+q1**2)))
          q2 = p1*(2.d0/(27.d0*p2))**onethird+(p2/54.d0)**onethird

          c1 = sqrt(a1**2.d0/4.d0-2*a2/3.d0+q2)
          c2 = a1**2.d0/2.d0-4.d0*a2/3.d0-q2
          c3 = (-a1**3.d0+4*a1*a2-8.d0*a3)/(4.d0*c1)

          r3 = -(a1/4.d0)+c1/2.d0-dsqrt(max(c2+c3,0.d0))/2.d0    ! is maximum
          r4 = -(a1/4.d0)+c1/2.d0+dsqrt(max(c2+c3,0.d0))/2.d0    ! is minimum
          r3 = min(1.d0,max(r3,-1.d0))
                  r4 = min(1.d0,max(r4,-1.d0))
                   
          tmax(i,j) = (acos(r4)+B1(i,j)/2.d0)/s%sigm(i,j)
          tmin(i,j) = (-acos(r3)+B1(i,j)/2.d0)/s%sigm(i,j)
          s%uon(i,j)  = amp1(i,j)*cos(s%sigm(i,j)*tmax(i,j))+  &
                        amp2(i,j)*cos(2.d0*s%sigm(i,j)*tmax(i,j)-B1(i,j))
          s%uoff(i,j) = abs(amp1(i,j)*cos(s%sigm(i,j)*tmin(i,j))+  &
                        amp2(i,j)*cos(2.d0*s%sigm(i,j)*tmin(i,j)-B1(i,j)))

          s%ua(i,j) = (s%uon(i,j)**4-s%uoff(i,j)**4)/(s%uon(i,j)**3-s%uoff(i,j)**3)  
      endif
   enddo
enddo

vmg  = sqrt(s%ue**2+s%ve**2)
uorb  = (s%uon**3+s%uoff**3)**onethird

do jg = 1,par%ngd

   ! cjaap: compute fall velocity with simple expression from Ahrens (2000)
   Te    = 20.d0
   kvis  = 4.d0/(20.d0+Te)*1d-5 ! Van rijn, 1993 
   Sster = s%D50(jg)/(4*kvis)*sqrt((par%rhos/par%rho-1)*par%g*s%D50(jg))
   cc1   = 1.06d0*tanh(0.064d0*Sster*exp(-7.5d0/Sster**2))
   cc2   = 0.22d0*tanh(2.34d0*Sster**(-1.18d0)*exp(-0.0064d0*Sster**2))
   wster = cc1+cc2*Sster
   par%w = wster*sqrt((par%rhos/par%rho-1.d0)*par%g*s%D50(jg))
   
   Ts       = par%tsfac*s%hh/par%w
   s%Tsg(:,:,jg) = max(Ts,0.2d0) 

   dster=25296*s%D50(jg)
   
   ! calculate treshold velocity Ucr
   if(s%D50(jg)<=0.0005) then
     Ucrc=0.19d0*s%D50(jg)**0.1d0*log10(4.d0*hloc/s%D90(jg)) 
         Ucrw=0.24d0*(1.65d0*par%g)**0.66d0*s%D50(jg)**0.33d0*par%Trep**0.33d0 ! Van Rijn uses Tp 
   else if(s%D50(jg)<0.002) then
     Ucrc=8.5d0*s%D50(jg)**0.6d0*log10(4.d0*hloc/s%D90(jg))                ! Shields
         Ucrw=0.95d0*(1.65d0*par%g)**0.57d0*s%D50(jg)**0.43*par%Trep**0.14     ! Komar and Miller (1975)
   else
#ifdef USEMPI
     write(*,'(a,i4)') 'In process',xmpi_rank
#endif
     write(*,*) ' Remark from vr2008: s%D50(jg) > 2mm, out of validity range'
   end if
   B2 = vmg/max(vmg+uorb,par%eps)
   Ucr = B2*Ucrc + (1-B2)*Ucrw

   ! transport parameters
   Asb=0.015d0*hloc*(s%D50(jg)/hloc)**1.2d0/(1.65d0*par%g*s%D50(jg))**0.75d0         ! bed load coefficent
   Ass=0.012d0*s%D50(jg)*dster**(-0.6d0)/(1.65d0*par%g*s%D50(jg))**1.2d0             ! suspended load coeffient

   term1= vmg+0.4*uorb                                                     ! nearbed-velocity
      
   ceq = 0*term1
   do j=1,s%ny+1
      do i=1,s%nx
         if(term1(i,j)>Ucr(i,j) .and. s%hh(i,j)>par%eps) then
            ceq(i,j)=Asb(i,j)*(term1(i,j)-Ucr(i,j))**1.5 + Ass*(term1(i,j)-Ucr(i,j))**2.4
         end if
      end do
   end do
  
   ceq = min(ceq,0.20d0)                 ! equilibrium concentration
   ceq = ceq/hloc
   s%ceqg(:,:,jg) = ceq*s%sedcal(jg)
enddo          

end subroutine vr2008

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sednew(s,par)
use params
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars),target            :: s
type(parameters)                  :: par

integer                           :: i,nw
integer                           :: j,jg
real*8                            :: dster,onethird,twothird,Ass,dt,dcf
real*8                            :: m1,m2,m3,m4,m5,m6
real*8                            :: Te,kvis,Sster,cc1,cc2,wster
real*8                            :: bb1,bb2,a1,a2,a3,a4,p1,q1,c1,c2,c3,r3,r4,ksi,phi,q2
complex*8                         :: p2

real*8 , dimension(:),allocatable,save   :: sn,t
real*8 , dimension(:,:),allocatable,save :: vmag2,Asb,Ts
real*8 , dimension(:,:),allocatable,save :: uorb2,Ucr,Ucrc,Ucrw,term1,B2
real*8 , dimension(:,:),allocatable,save :: hloc,ceq
real*8 , dimension(:,:),allocatable,save :: Ur,Bm,B1,sigu2,amp2,amp1,As,Sk
real*8 , dimension(:,:),allocatable,save :: tmin,tmax,tmss
real*8 , dimension(:,:),allocatable,save :: kb ! wwvv: added. After removing 
                                               ! s.ind and s.inc, it appeared
! that with s.ind and s.inc kb was in spacepars

! wwvv removed s.ind and s.inc and put s% in appropriate locations
!      to be able to better analyse the code

nw = 50

if (.not. allocated(vmag2)) then
   allocate (vmag2 (s%nx+1,s%ny+1))
   allocate (term1 (s%nx+1,s%ny+1))
   allocate (B2    (s%nx+1,s%ny+1))
   allocate (Asb   (s%nx+1,s%ny+1))
   allocate (Ucr   (s%nx+1,s%ny+1))
   allocate (Ucrc  (s%nx+1,s%ny+1))
   allocate (Ucrw  (s%nx+1,s%ny+1))
   allocate (uorb2 (s%nx+1,s%ny+1))
   allocate (hloc  (s%nx+1,s%ny+1))
   allocate (kb    (s%nx+1,s%ny+1))
   allocate (Ts    (s%nx+1,s%ny+1))
   allocate (ceq   (s%nx+1,s%ny+1))
   allocate (Ur    (s%nx+1,s%ny+1))
   allocate (Bm    (s%nx+1,s%ny+1))
   allocate (B1    (s%nx+1,s%ny+1))
   allocate (sigu2 (s%nx+1,s%ny+1))
   allocate (As    (s%nx+1,s%ny+1))
   allocate (Sk    (s%nx+1,s%ny+1))
   allocate (tmax  (s%nx+1,s%ny+1))
   allocate (tmin  (s%nx+1,s%ny+1))
   allocate (tmss  (s%nx+1,s%ny+1))
   allocate (amp2  (s%nx+1,s%ny+1))
   allocate (amp1  (s%nx+1,s%ny+1))
   allocate (sn      (nw+1))            ! wwvv not used
   allocate (t       (50))
endif

hloc   = max(s%hh,par%hmin) !Jaap par%hmin instead of par%eps
onethird=1.d0/3.d0
twothird=2.d0/3.d0
!
! compute wave shape short waves (Ruessink and van Rijn 2007)
!
m1 = 0         ! a = 0
m2 = 0.7939d0  ! b = 0.79 +/- 0.023
m3 = -0.6065d0 ! c = -0.61 +/- 0.041
m4 = 0.3539d0  ! d = -0.35 +/- 0.032 
m5 = 0.6373d0  ! e = 0.64 +/- 0.025
m6 = 0.5995d0  ! f = 0.60 +/- 0.043

dt = par%Trep/50
do i = 1,50
   t(i) = (i-1)*dt
enddo


do j=1,s%ny+1     
    do i=1,s%nx+1
       if (s%k(i,j)*s%H(i,j)<par%px/2.d0 .and. s%H(i,j)>0.01d0) then
         Ur(i,j) = 3.d0/8.d0*sqrt(2.d0)*s%H(i,j)*s%k(i,j)/(s%k(i,j)*hloc(i,j))**3 !Ursell number
         Bm(i,j) = m1+(m2-m1)/(1.d0+exp((m3-log10(Ur(i,j)))/m4)) !Boltzmann sigmoid (eq 6) 
         B1(i,j) = -90.d0+90.d0*tanh(m5/Ur(i,j)**m6)    
         B1(i,j) = B1(i,j)*par%px/180.d0

         Sk(i,j) = Bm(i,j)*cos(B1(i,j)) !Skewness (eq 8)
         As(i,j) = Bm(i,j)*sin(B1(i,j)) !Skewness (eq 9)

         sigu2(i,j) = (s%H(i,j)*s%sigm(i,j)/(2.d0*sqrt(2.d0)*sinh(s%k(i,j)*hloc(i,j))))**2
                                        !variance from linear theory (eq 13)

         p1 = -2.d0*sigu2(i,j)
         q1 = -4.d0/3.d0*Sk(i,j)*sigu2(i,j)**1.5d0/cos(B1(i,j))

         ksi = q1/2.d0/sqrt((-p1/3.d0)**3)
         phi = atan2(sqrt(1.d0-ksi**2.d0),ksi)

         amp2(i,j) = 2.d0*sqrt(-p1/3.d0)*cos(phi/3.d0+4.d0/3.d0*par%px)
         amp1(i,j) = sqrt(2.d0*sigu2(i,j)-amp2(i,j)**2)
         !
         ! compute uon and uoff
         !
         bb1 = -amp1(i,j)*s%sigm(i,j)
         bb2 = -2.d0*amp2(i,j)*s%sigm(i,j)
         a1 = bb1/bb2*cos(B1(i,j)/2.d0); 
         a2 = (4.d0*bb2**2.d0-bb1**2)/(-4.d0*bb2**2)
         a3 = -a1
         a4 = bb1**2.d0*cos(B1(i,j)/2.d0)**2/(-4.d0*bb2**2)

         p1 = a2**2-3.d0*a1*a3+12.d0*a4
         q1 = 2.d0*a2**3-9.d0*a1*a2*a3+27.d0*a3**2+27.d0*a4*a1**2-72.d0*a2*a4

!        p2 = q1 + cdsqrt(cmplx(-4.d0*p1**3.d0+q1**2.d0))
         p2 = q1 + csqrt(cmplx(-4.d0*p1**3+q1**2))
         q2 = p1*(2.d0/(27.d0*p2))**onethird+(p2/54.d0)**onethird

         c1 = sqrt(a1**2/4.d0-2*a2/3.d0+q2)
         c2 = a1**2/2.d0-4.d0*a2/3.d0-q2
         c3 = (-a1**3+4*a1*a2-8.d0*a3)/(4.d0*c1)

         r3 = -(a1/4.d0)+c1/2.d0-dsqrt(max(c2+c3,0.d0))/2.d0    ! is maximum
         r4 = -(a1/4.d0)+c1/2.d0+dsqrt(max(c2+c3,0.d0))/2.d0    ! is minimum

         r3 = min(1.d0,max(r3,-1.d0))
         r4 = min(1.d0,max(r4,-1.d0))
                   
         tmax(i,j) = (acos(r4)+B1(i,j)/2.d0)/s%sigm(i,j)
         tmin(i,j) = (-acos(r3)+B1(i,j)/2.d0)/s%sigm(i,j)
         s%uon(i,j)  = amp1(i,j)*cos(s%sigm(i,j)*tmax(i,j))+ &
                         amp2(i,j)*cos(2.d0*s%sigm(i,j)*tmax(i,j)-B1(i,j))
         s%uoff(i,j) = abs(amp1(i,j)*cos(s%sigm(i,j)*tmin(i,j))+ &
                         amp2(i,j)*cos(2.d0*s%sigm(i,j)*tmin(i,j)-B1(i,j)))

         s%Tbore(i,j) = 0.5d0*min(abs(tmax(i,j)-tmin(i,j)),par%Trep-abs(tmax(i,j)-tmin(i,j)))        
         !
         ! compute maximum wave surface slope
         !
         bb1 = -amp1(i,j)*s%sigm(i,j)**2
         bb2 = -4.d0*amp2(i,j)*s%sigm(i,j)**2
         a1 = bb1/bb2*cos(B1(i,j)/2.d0)
         a2 = -1.d0+bb1**2/(4.d0*bb2**2)
         a3 = -0.5d0*a1
         a4 = (bb2**2-(bb1*sin(B1(i,j)/2.d0))**2)/(4*bb2**2)

         p1 = a2**2-3.d0*a1*a2+12.d0*a4
         q1 = 2.d0*a2**3-9.d0*a1*a2*a3+27.d0*a3**2+27.d0*a4*a1**2-72.d0*a2*a4

!        p2 = q1 + cdsqrt(cmplx(-4.d0*p1**3.d0+q1**2.d0))
         p2 = q1 + sqrt(cmplx(-4.d0*p1**3+q1**2))
         q2 = p1*(2.d0/(27.d0*p2))**onethird+(p2/54.d0)**onethird

         c1 = sqrt(a1**2/4.d0-2*a2/3.d0+q2)
         c2 = a1**2/2.d0-4.d0*a2/3.d0-q2
         c3 = (-a1**3+4*a1*a2-8.d0*a3)/(4.d0*c1)

         r4 = -(a1/4.d0)+c1/2.d0+dsqrt(max(c2+c3,0.d0))/2.d0    ! is minimum
         r4 = min(1.d0,max(r4,-1.d0))

         tmss(i,j) = (-acos(r4)+B1(i,j)/2)/s%sigm(i,j)
         s%BR(i,j) = -amp1(i,j)*s%sigm(i,j)*sin(s%sigm(i,j)*tmss(i,j))- &
                     2.d0*amp2(i,j)*s%sigm(i,j)*sin(2*s%sigm(i,j)*tmss(i,j)-B1(i,j))
         s%BR(i,j) = s%BR(i,j)*cosh(s%k(i,j)*hloc(i,j))/par%g
          !
                  ! compute velocity due to wave assymetry
                  !
          !u4 = sum(par%dt*(amp1(i,j)*cos(s%sigm(i,j)*t)+amp2(i,j)*cos(2.d0*s%sigm(i,j)*t-B1(i,j))*abs(amp1(i,j)*cos(s%sigm(i,j)*t)+ amp2(i,j)*cos(2.d0*s%sigm(i,j)*t-B1(i,j))))*2.4d0)
          !u3 = sum(par%dt*(abs(amp1(i,j)*cos(s%sigm(i,j)*t)+amp2(i,j)*cos(2.d0*s%sigm(i,j)*t-B1(i,j))))*2.4d0)
                  !ua(i,j) = u4/u3 
       else
         s%uon(i,j) = sqrt((s%H(i,j)*s%sigm(i,j)/(2.d0*sqrt(2.d0)*sinh(s%k(i,j)*hloc(i,j))))**2)
         s%uoff(i,j)  = s%uon(i,j)
         s%Tbore(i,j) = par%Trep/4
         s%BR(i,j)    = par%Beta
       endif
    enddo
enddo
s%Tbore = max(s%Tbore,par%Trep/50)

!dt = par%Trep/50
!do i = 1,50;
!   t = (i-1)*dt
!   u4 = par%dt*(amp1*cos(s%sigm*t)+amp2*cos(2.d0*s%sigm*t-B1)*abs(amp1*cos(s%sigm*t)+amp2*cos(2.d0*s%sigm*t-B1)))*2.4d0
!   u3 = par%dt*(abs(amp1*cos(s%sigm*t)+amp2*cos(2.d0*s%sigm*t-B1)))*2.4d0
!enddo
!ua = u4/u3

!
! compute near bed turbulence
!
do j=1,s%ny+1     
    do i=1,s%nx+1
           dcf = min(1.d0,1.d0/(exp(hloc(i,j)/max(s%H(i,j),0.1d0)) -1.d0))
           kb(i,j) = (s%DR(i,j)/par%rho)**twothird*dcf*par%Trep/s%Tbore(i,j)
           vmag2(i,j) = s%ue(i,j)**2+s%ve(i,j)**2
    enddo
enddo

vmag2  = s%ue**2+s%ve**2
!uorb2  = urms**2.d0+0.5*kb !Van Rijn use Peak orbital flow velocity
uorb2  = (s%uon**3+s%uoff**3)**twothird

do jg = 1,par%ngd

   ! cjaap: compute fall velocity with simple expression from Ahrens (2000)
   Te    = 20.d0
   kvis  = 4.d0/(20.d0+Te)*1d-5 ! Van rijn, 1993 
   Sster = s%D50(jg)/(4*kvis)*sqrt((par%rhos/par%rho-1)*par%g*s%D50(jg))
   cc1   = 1.06d0*tanh(0.064d0*Sster*exp(-7.5d0/Sster**2))
   cc2   = 0.22d0*tanh(2.34d0*Sster**(-1.18d0)*exp(-0.0064d0*Sster**2))
   wster = cc1+cc2*Sster
   par%w = wster*sqrt((par%rhos/par%rho-1.d0)*par%g*s%D50(jg))
   
   where (s%ceqg(:,:,jg)>s%ccg(:,:,jg))  Ts = s%hh/max(sqrt(kb),par%w)
   where (s%ceqg(:,:,jg)<=s%ccg(:,:,jg)) Ts = s%hh/par%w
   s%Tsg(:,:,jg) = max(Ts,0.2d0) 

   dster=25296*s%D50(jg)
   
   ! calculate treshold velocity Ucr
   if(s%D50(jg)<=0.0005) then
     Ucrc=0.19d0*s%D50(jg)**0.1d0*log10(4.d0*hloc/s%D90(jg)) 
         Ucrw=0.24d0*(1.65d0*par%g)**0.66d0*s%D50(jg)**0.33d0*par%Trep**0.33d0 ! Van Rijn uses Tp 
   else if(s%D50(jg)<0.002) then
     Ucrc=8.5d0*s%D50(jg)**0.6d0*log10(4.d0*hloc/s%D90(jg))                ! Shields
         Ucrw=0.95d0*(1.65d0*par%g)**0.57d0*s%D50(jg)**0.43*par%Trep**0.14     ! Komar and Miller (1975)
   else
#ifdef USEMPI
     write(*,'(a,i4)') 'In process',xmpi_rank
#endif
     write(*,*) ' Remark from sednew: s%D50(jg) > 2mm, out of validity range'
   end if
   B2 = sqrt(vmag2)/max(sqrt(vmag2)+sqrt(uorb2),par%eps)
   Ucr = B2*Ucrc + (1-B2)*Ucrw

   ! transport parameters
   Asb=0.015d0*hloc*(s%D50(jg)/hloc)**1.2d0/(1.65d0*par%g*s%D50(jg))**0.75d0         ! bed load coefficent
   Ass=0.012d0*s%D50(jg)*dster**(-0.6d0)/(1.65d0*par%g*s%D50(jg))**1.2d0             ! suspended load coeffient

   term1= sqrt(vmag2+0.5d0*(0.36d0*uorb2**2+2.d0*kb)) !+0.5*kb                    ! nearbed-velocity
      
   ceq = 0*term1
   do j=1,s%ny+1
      do i=1,s%nx
         if(term1(i,j)>Ucr(i,j) .and. s%hh(i,j)>par%eps) then
            ceq(i,j)=Asb(i,j)*(term1(i,j)-Ucr(i,j))**1.5 + Ass*(term1(i,j)-Ucr(i,j))**2.4
         end if
      end do
   end do
  
   ceq = min(ceq,0.20d0)                 ! equilibrium concentration
   ceq = ceq/hloc
   s%ceqg(:,:,jg) = ceq*s%sedcal(jg)
enddo                             ! end of grain size classes

end subroutine sednew

end module morphevolution
