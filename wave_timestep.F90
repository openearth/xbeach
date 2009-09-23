module wave_timestep_module
contains
subroutine wave_timestep(s,par)
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
    use roelvink_module
    use xmpi_module
    use mnemmodule

    IMPLICIT NONE

    type(spacepars), target     :: s
    type(parameters)            :: par

    integer                     :: i
    integer                     :: j
    integer                     :: itheta



    integer, dimension(:,:,:),allocatable,save  :: wete
    real*8 , dimension(:,:)  ,allocatable,save  :: dhdx,dhdy,dudx,dudy,dvdx,dvdy,ustw,Erfl
    real*8 , dimension(:,:)  ,allocatable,save  :: km,kmx,kmy,wm,xwadvec,ywadvec,sinh2kh
    real*8 , dimension(:,:,:),allocatable,save  :: xadvec,yadvec,thetaadvec,dd,drr
    real*8 , dimension(:,:,:),allocatable,save  :: xradvec,yradvec,thetaradvec
	real*8 , dimension(:,:)  ,allocatable,save  :: dkmxdx,dkmxdy,dkmydx,dkmydy,cgxm,cgym,arg,fac
	real*8 , dimension(:,:)  ,allocatable,save  :: wcifacu,wcifacv
	real*8                                      :: factime

    include 's.ind'
    include 's.inp'

    if (.not. allocated(wete)) then
       allocate(drr         (nx+1,ny+1,ntheta))
       allocate(wete        (nx+1,ny+1,ntheta))
       allocate(xadvec      (nx+1,ny+1,ntheta))
       allocate(yadvec      (nx+1,ny+1,ntheta))
       allocate(thetaadvec  (nx+1,ny+1,ntheta))
       allocate(xradvec     (nx+1,ny+1,ntheta))
       allocate(yradvec     (nx+1,ny+1,ntheta))
       allocate(thetaradvec (nx+1,ny+1,ntheta))
       allocate(dd          (nx+1,ny+1,ntheta))
       allocate(dhdx        (nx+1,ny+1))
       allocate(dhdy        (nx+1,ny+1))
       allocate(dudx        (nx+1,ny+1))
       allocate(dudy        (nx+1,ny+1))
       allocate(dvdx        (nx+1,ny+1))
       allocate(dvdy        (nx+1,ny+1))
       allocate(km          (nx+1,ny+1))
       allocate(kmx         (nx+1,ny+1))
       allocate(kmy         (nx+1,ny+1))
       allocate(wm          (nx+1,ny+1))
       allocate(ustw        (nx+1,ny+1))
       allocate(Erfl        (nx+1,ny+1)) ! wwvv not used
       allocate(xwadvec     (nx+1,ny+1))
       allocate(ywadvec     (nx+1,ny+1))
       allocate(sinh2kh     (nx+1,ny+1))
       allocate(dkmxdx      (nx+1,ny+1))
       allocate(dkmxdy      (nx+1,ny+1))
       allocate(dkmydx      (nx+1,ny+1))
       allocate(dkmydy      (nx+1,ny+1))
       allocate(cgxm        (nx+1,ny+1))
       allocate(cgym        (nx+1,ny+1))
	   allocate(arg         (nx+1,ny+1))
	   allocate(fac         (nx+1,ny+1))
	   allocate(wcifacu     (nx+1,ny+1))
	   allocate(wcifacv     (nx+1,ny+1))


! wwvv todo: I think these iniailization are superfluous
       drr         = 0.d0
       wete        = 0.d0
       xadvec      = 0.d0
       yadvec      = 0.d0
       thetaadvec  = 0.d0
       xradvec     = 0.d0
       yradvec     = 0.d0
       thetaradvec = 0.d0
       dd          = 0.d0
       dhdx        = 0.d0
       dhdy        = 0.d0
       dudx        = 0.d0
       dudy        = 0.d0
       dvdx        = 0.d0
       dvdy        = 0.d0
       km          = 0.d0
       kmx         = 0.d0
       kmy         = 0.d0
       wm          = 0.d0
       ustw        = 0.d0
       Erfl        = 0.d0
       xwadvec     = 0.d0
       ywadvec     = 0.d0
       sinh2kh     = 0.d0
       dkmxdx      = 0.d0
       dkmxdy      = 0.d0
       dkmydx      = 0.d0
       dkmydy      = 0.d0
       cgxm        = 0.d0
       cgym        = 0.d0
	   arg         = 0.d0
       fac         = 0.d0
       Fx          = 0.d0 ! in spacepars
       Fy          = 0.d0 ! in spacepars
    endif

    hh = max(hh,par%eps)

! Calculate once velocities used with and without wave current interaction
    wcifacu=u*par%wci*min(hh/par%hwci,1.d0)
    wcifacv=v*par%wci*min(hh/par%hwci,1.d0)

! Dispersion relation
    if (par%wci .ne. 0) then
       if (par%t==par%dt) then
          sigm = max((sum(sigt,3)/ntheta),0.01d0)
          call dispersion(par,s)
	      umwci = 0.d0
          vmwci = 0.d0
	      zswci = zs
	      km=k
          if (xmaster) write(*,*)'wave current interaction included'
       endif
       km(1,:) = k(1,:)   ! boundary condition *assuming zero flow at the boundary) 
#ifdef USEMPI
       call xmpi_shift(km,'1:')
#endif
       factime = 1.d0/par%cats/par%Trep*par%dt
       umwci   = factime*uu + (1-factime)*umwci
       vmwci   = factime*vv + (1-factime)*vmwci	
       zswci   = factime*zs + (1-factime)*zswci	
       arg     = min(100.0d0,km*(hh+par%delta*H))
       sigm(1,:) = sqrt( par%g*km(1,:)*tanh(arg(1,:))) ! *( 1.d0+ ((km(1,:)*H(1,:)/2.d0)**2)))
!  calculate change in intrinsic frequency
       tm  = (sum(ee*thet,3)/ntheta)/(max(sum(ee,3),0.00001d0)/ntheta)
       kmx = km*cos(tm)
       kmy = km*sin(tm)
       wm = sigm+kmx*umwci*par%wci*min((zswci-zb)/par%hwci,1.d0)+kmy*vmwci*par%wci*min((zswci-zb)/par%hwci,1.d0)

       cgym = cg*sin(tm) + vmwci*min((zswci-zb)/par%hwci,1.d0)
       cgxm = cg*cos(tm) + umwci*min((zswci-zb)/par%hwci,1.d0)

       call slope2D(kmx,nx,ny,xz,yz,dkmxdx,dkmxdy)
       call slope2D(kmy,nx,ny,xz,yz,dkmydx,dkmydy)
       call advecwx(wm,xwadvec,kmx,nx,ny,xz)   ! cjaap: xz or xu?
       kmx = kmx -par%dt*xwadvec  -1.0d0*par%dt*cgym*(dkmydx-dkmxdy)
       kmx(:,ny+1) = kmx(:,ny)  ! lateral bc
       kmx(:,1) = kmx(:,2)  ! lateral bc
   
! wwvv the following has consequences for the // version todo
#ifdef USEMPI
       call xmpi_shift(kmx,':n')  ! get column kml(:ny+1) from right neighbour
       call xmpi_shift(kmx,':1')
#endif
       call advecwy(wm,ywadvec,kmy,nx,ny,yz)   ! cjaap: yz or yv?
       kmy = kmy-par%dt*ywadvec  + 1.0*par%dt*cgxm*(dkmydx-dkmxdy)
       kmy(:,ny+1) = kmy(:,ny)   ! lateral bc
       kmy(:,1) = kmy(:,2)   ! lateral bc
! wwvv the following has consequences for the // version todo
#ifdef USEMPI
       call xmpi_shift(kmy,':n')
       call xmpi_shift(kmy,':1')
#endif
! update km
       km = sqrt(kmx**2+kmy**2)
! non-linear dispersion
       arg = min(100.0d0,km*((zswci-zb)+par%delta*H))
       arg = max(arg,0.0001)
!       fac = ( 1.d0 + ((km*H/2.d0)**2)*( (8.d0+(cosh(min(4.d0*arg,10.0d0)))**1.d0-2.d0*(tanh(arg))**2.d0 ) /(8.d0*(sinh(arg))**4.d0) ) )
       fac = ( 1.d0 + ((km*H/2.d0)**2))  ! use deep water correction instead of expression above (waves are short near blocking point anyway)
!       fac = 1.d0    ! Linear
       sigm = sqrt( par%g*km*tanh(arg)*fac)

!  update intrinsic frequency
       do itheta=1,ntheta
          sigt(:,:,itheta) = sigm
       enddo
       where(km>0.01d0)
          c  = sigm/km
!          cg = c*(0.5d0+arg/sinh(2.0d0*arg))    ! Linear
		  cg = c*(0.5d0+arg/sinh(2*arg))*sqrt(fac)  ! &  to include more
!		 	          + km*(H/2)**2*sqrt(max(par%g*km*tanh(arg),0.001d0))/sqrt(max(fac,0.001d0)) ! include wave steepness
          n=0.5d0+km*hh/sinh(2*max(km,0.00001d0)*hh)
       elsewhere
          c  = 0.01d0
	      cg = 0.01d0
		  n  = 1.d0
       endwhere
!  update k
       km = min(km,25.d0) ! limit to gravity waves
       k = km

    else  ! no wave current interaction
       sigm = max((sum(sigt,3)/ntheta),0.01d0)
       call dispersion(par,s)
    endif ! end wave current interaction

! Slopes of water depth
    call slope2D(hh+par%delta*H,nx,ny,xz,yz,dhdx,dhdy)
    call slope2D(wcifacu,nx,ny,xz,yz,dudx,dudy)
    call slope2D(wcifacv,nx,ny,xz,yz,dvdx,dvdy)
!
! Calculate once sinh(2kh)
    where(2*hh*k<=3000.d0)
       sinh2kh=sinh(min(2*k*(hh+par%delta*H),10.0d0))
    elsewhere
       sinh2kh = 3000.d0
    endwhere

    do itheta=1,ntheta
       cgx(:,:,itheta)= cg*cxsth(itheta)+wcifacu
       cgy(:,:,itheta)= cg*sxnth(itheta)+wcifacv
       cx(:,:,itheta) =  c*cxsth(itheta)+wcifacu
       cy(:,:,itheta) =  c*sxnth(itheta)+wcifacv
       ctheta(:,:,itheta)=  &
         sigm/sinh2kh*(dhdx*sxnth(itheta)-dhdy*cxsth(itheta)) + &
		 par%wci*(&
         cxsth(itheta)*(sxnth(itheta)*dudx - cxsth(itheta)*dudy) + &
         sxnth(itheta)*(sxnth(itheta)*dvdx - cxsth(itheta)*dvdy))
    enddo
!
! transform to wave action
!
    ee = ee/sigt
!
! Upwind Euler timestep propagation
!
    call advecxho(ee,cgx,xadvec,nx,ny,ntheta,xz,par%dt,par%scheme)
    call advecyho(ee,cgy,yadvec,nx,ny,ntheta,yz,par%dt,par%scheme)
    call advectheta(ee*ctheta,thetaadvec,nx,ny,ntheta,dtheta)
!
    ee=ee-par%dt*(xadvec+yadvec+thetaadvec)
!
! transform back to wave energy
!
    ee = ee*sigt
    ee=max(ee,0.0d0) !Jaap

!
! Energy integrated over wave directions,Hrms
!
    E=sum(ee,3)*dtheta
    H=sqrt(E/par%rhog8)
    do itheta=1,ntheta
       ee(:,:,itheta)=ee(:,:,itheta)/max(1.d0,(H/(par%gammax*hh))**2)
    enddo
    H=min(H,par%gammax*hh)
    E=par%rhog8*H**2

! Total dissipation
    if(par%break == 1 .or. par%break == 3)then
        call roelvink(par,s,km)
    else if(par%break == 2)then
        call baldock(par,s,km)
    else if (par%break == 4) then
        cgxm = cg*cos(tm) 
        cgym = cg*sin(tm)
        call advecqx(cgxm,Qb,xwadvec,nx,ny,xz)
        call advecqy(cgym,Qb,ywadvec,nx,ny,yz)
        Qb=Qb-par%dt*(xwadvec+ywadvec)
        call roelvink(par,s,km)        
    endif
!
! Distribution of dissipation over directions and frequencies
!
    do itheta=1,ntheta
! Only calculate for E>0 FB
        dd(:,:,itheta)=ee(:,:,itheta)*D/max(E,0.00001d0)
    enddo

    do j=1,ny+1
        do i=1,nx+1
! cjaap: replaced par%hmin by par%eps
            if(hh(i,j)+par%delta*H(i,j)>par%eps) then
                wete(i,j,1:ntheta)=1
            else
                wete(i,j,1:ntheta)=0
            end if
        end do
    end do
!
! Euler step dissipation
!
! calculate roller energy balance
!
    call advecxho(rr,cx,xradvec,nx,ny,ntheta,xz,par%dt,par%scheme)
    call advecyho(rr,cy,yradvec,nx,ny,ntheta,yz,par%dt,par%scheme)
    call advectheta(rr*ctheta,thetaradvec,nx,ny,ntheta,dtheta)

    rr=rr-par%dt*(xradvec+yradvec+thetaradvec)
    rr=max(rr,0.0d0)
!
! euler step roller energy dissipation (source and sink function)
!
    do itheta=1,ntheta
        do j=1,ny+1
            do i=1,nx+1
               if(wete(i,j,itheta)==1) then
                    ee(i,j,itheta)=ee(i,j,itheta)-par%dt*dd(i,j,itheta)
                    if(par%roller==1) then
                        rr(i,j,itheta)=rr(i,j,itheta)+par%dt*dd(i,j,itheta)           &
                                  -par%dt*2*par%g*BR(i,j)*rr(i,j,itheta)        &
                                  /sqrt(cx(i,j,itheta)**2+cy(i,j,itheta)**2)
                        drr(i,j,itheta) = 2*par%g*BR(i,j)*max(rr(i,j,itheta),0.0d0)/   &
                                  sqrt(cx(i,j,itheta)**2 +cy(i,j,itheta)**2)
                    else if (par%roller==0) then
                        rr(i,j,itheta)= 0.0d0
                        drr(i,j,itheta)= 0.0d0
                    endif
                    ee(i,j,itheta)=max(ee(i,j,itheta),0.0d0)
                    rr(i,j,itheta)=max(rr(i,j,itheta),0.0d0)
                elseif(wete(i,j,itheta)==0) then
                    ee(i,j,itheta)=0.0d0
                    rr(i,j,itheta)=0.0d0
                end if
            end do
        end do
    end do
!
! Bay boundary Robert + Jaap
! 
! wwvv 
! this has consequences for the parallel version, 
! but also, if we do  nothing, there are discrepancies
! between ee and rr in the different processes. We need to
! get valid values for ee(nx+1,:,:) and rr(nx+1,:,:) from
! the neighbour below. We cannot postpone this until this
! subroutine ends, because ee and rr are used in this subroutine
    ee(nx+1,:,:) =ee(nx,:,:)
    rr(nx+1,:,:) =rr(nx,:,:)
#ifdef USEMPI
    call xmpi_shift(ee,'m:')  ! fill in ee(nx+1,:,:)
    call xmpi_shift(rr,'m:')  ! fill in rr(nx+1,:,:)
#endif
!
!
! Compute mean wave direction
!
    thetamean=(sum(ee*thet,3)/size(ee,3))/(max(sum(ee,3),0.00001d0)/size(ee,3))

!
! Energy integrated over wave directions,Hrms
!
    E  = sum(ee,3)*dtheta
    R  = sum(rr,3)*dtheta
    DR = sum(drr,3)*dtheta
    H  = sqrt(E/par%rhog8)
!
! Radiation stresses and forcing terms
!
! n=cg/c   (Robert: calculated earlier in dispersion relation)
    Sxx=(n*sum((1.d0+(costhet)**2)*ee,3)-.5d0*sum(ee,3))*dtheta
    Syy=(n*sum((1.d0+(sinthet)**2)*ee,3)-.5d0*sum(ee,3))*dtheta
    Sxy=n*sum(sinthet*costhet*ee,3)*dtheta

! add roller contribution
    Sxx = Sxx + sum((costhet**2)*rr,3)*dtheta
    Syy = Syy + sum((sinthet**2)*rr,3)*dtheta
    Sxy = Sxy + sum(sinthet*costhet*rr,3)*dtheta

    do j=2,ny
        do i=1,nx
           Fx(i,j)=-(Sxx(i+1,j)-Sxx(i,j))/(xz(i+1)-xz(i))                                   &
                   -0.5*(Sxy(i,j+1)+Sxy(i+1,j+1)- Sxy(i,j-1)-Sxy(i+1,j-1))/(yz(j+1)-yz(j-1))
        enddo
    enddo

    do j=1,ny
        do i=2,nx
           Fy(i,j)=-(Syy(i,j+1)-Syy(i,j))/(yz(j+1)-yz(j))                                  &
                   -0.5d0*(Sxy(i+1,j)+Sxy(i+1,j+1)-Sxy(i-1,j)-Sxy(i-1,j+1))/(xz(i+1)-xz(i-1))
        enddo
    enddo
! wwvv in the previous, Fx and Fy are computed. The missing elements
!  elements are Fx(:,1), Fx(nx+1,:), Fx(:,ny+1)
!               Fy(1,:), Fy(nx+1,:), Fy(:,ny+1)

! wwvv todo the following has consequences for // version
    Fy(1,:)=Fy(2,:)
    Fx(nx+1,:) = 0.0d0
    Fy(nx+1,:) = 0.0d0
    Fx(:,1)=Fx(:,2)     
! Robert: Fix Neumann assumption for wave forcing on boudary, even if ee not Neumanned
    Fy(:,1)=Fy(:,2)
    Fy(:,ny+1)=Fy(:,ny)   

! wwvv so, Fx(:ny+1) and Fy(:ny+1) are left zero and Fx(nx+1,:) and Fy(nx+1,:)
! are made zero.  In the parallel case, Fx(:,1) and Fy(1,:) don't get a 
! value if the submatrices are not on suitable border. 
! I guess that it is necessary to communicate with neighbours the values of these elements

#ifdef USEMPI
    call xmpi_shift(Fx,':1')  ! shift in Fx(:,1)
    call xmpi_shift(Fx,'m:')  ! shift in Fx(nx+1,:)
    call xmpi_shift(Fx,':n')  ! shift in Fx(:,ny+1)
    call xmpi_shift(Fy,'1:')  ! shift in Fy(1,:)
    call xmpi_shift(Fy,'m:')  ! shift in Fy(nx+1,:)
    call xmpi_shift(Fy,':n')  ! shift in Fy(:,ny+1)
#endif

! Ad
    urms=par%px*H/par%Trep/(sqrt(2.d0)*sinh(min(max(k,0.01d0)*(hh+par%delta*H),10.0d0)))
!   urms=par%px*H/par%Trep/(sqrt(2.d0)*sinh(k*(hh+par%delta*H)))

    ustw= E/max(c,sqrt(par%hmin*par%g))/par%rho/max(hh,par%hmin)   ! Jaap
    uwf = ustw*cos(tm)
    vwf = ustw*sin(tm)
! roller contribution
    ustr=2.*R/max(c,sqrt(par%hmin*par%g))/par%rho/max(hh,par%hmin) ! Jaap
! introduce breaker delay
    call breakerdelay(par,s)
    ust=usd+ustw
!lateral boundaries
! wwvv todo the following has consequences for // version
    ust(1,:) = ust(2,:)
    ust(:,1) = ust(:,2)
    ust(:,ny+1) = ust(:,ny)

#ifdef USEMPI
    call xmpi_shift(ust,'1:')  ! get ust(1,:) from above
    call xmpi_shift(ust,':1')  ! get ust(:,1) from left
    call xmpi_shift(ust,':n')  ! get ust(:,ny+1) from right
#endif
end subroutine wave_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine slope2D(h,nx,ny,x,y,dhdx,dhdy)
use xmpi_module
IMPLICIT NONE

integer                           :: i,j,nx,ny
real*8, dimension(nx+1,ny+1)      :: h,dhdx,dhdy
real*8, dimension(nx+1)           :: x
real*8, dimension(ny+1)           :: y


! wwvv dhdx(2:nx,:) is computed, dhdx(1,:) and dhdx(nx+1,:) 
! get boundary values, so in the parallel case, we need
! to do something about that: get the boundaries from
! upper and lower neighbours

do j=1,ny+1
    if(nx+1>=2)then
      do i=2,nx
          dhdx(i,j)=(h(i+1,j)-h(i-1,j))/(x(i+1)-x(i-1))
      end do  
      dhdx(1,j)=(h(2,j)-h(1,j))/(x(2)-x(1))
      dhdx(nx+1,j)=(h(nx+1,j)-h(nx,j))/(x(nx+1)-x(nx))
    end if
end do
#ifdef USEMPI
call xmpi_shift(dhdx,'m:')  ! fill in dhdx(nx+1,:)
call xmpi_shift(dhdx,'1:')  ! fill in dhdx(1,:)
#endif

do i=1,nx+1
    if(ny+1>=2)then
      do j=2,ny
          dhdy(i,j)=(h(i,j+1)-h(i,j-1))/(y(j+1)-y(j-1))
      end do
      dhdy(i,1)=(h(i,2)-h(i,1))/(y(2)-y(1))
      dhdy(i,ny+1)=(h(i,ny+1)-h(i,ny))/(y(ny+1)-y(ny))
    end if
end do

#ifdef USEMPI
call xmpi_shift(dhdy,':n')  !  fill in dhdy(:,ny+1)
call xmpi_shift(dhdy,':1')  !  fill in dhdy(:,1)
#endif
end subroutine slope2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine advecx(arrin,xadvec,nx,ny,ntheta,xz)
use xmpi_module

IMPLICIT NONE

integer                                         :: i,j,nx,ny,ntheta
integer                                         :: itheta
real*8 , dimension(nx+1)                        :: xz
real*8 , dimension(nx+1,ny+1,ntheta)            :: xadvec,arrin
real*8                                          :: dxmin_i,dxplus_i,dxcent_i

xadvec = 0.d0

do itheta=1,ntheta
    do i=2,nx
        dxmin_i  = 1.d0/(xz(i)-xz(i-1))
        dxplus_i = 1.d0/(xz(i+1)-xz(i))
        dxcent_i = 1.d0/(xz(i+1)-xz(i-1))
        do j=1,ny+1
           if (arrin(i,j,itheta)>0) then
              xadvec(i,j,itheta)=(arrin(i,j,itheta)-arrin(i-1,j,itheta))*dxmin_i
           elseif (arrin(i,j,itheta)<0) then
              xadvec(i,j,itheta)=(arrin(i+1,j,itheta)-arrin(i,j,itheta))*dxplus_i
           else
              xadvec(i,j,itheta)=(arrin(i+1,j,itheta)-arrin(i-1,j,itheta))*dxcent_i
           endif
        end do
    end do
end do
! wwvv xadvec(1,:,:) and xadvec(nx+1,:,:) are not determined, in the parallel
! case this is not desirable. We copy these rows from top and bot neighbours
!
! However:
! wwvv todo: the parameters are somewhat misleading, not always is nx the
! nx from spaceparams, so this correction is not desirable here
! The same for the other advec subroutines. I commented the xmpi_shift
! calls out. ( it seems that nx is always 2)
!
!#ifdef USEMPI
!call xmpi_shift(xadvec,'m:')  ! fill in xadvec(nx+1,:,:)
!call xmpi_shift(xadvec,'1:')  ! fill in xadvec(1,:,:)
!#endif

end subroutine advecx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine advecxho(ee,cgx,xadvec,nx,ny,ntheta,xz,dt,scheme)

use xmpi_module
IMPLICIT NONE

integer                                         :: i,j,nx,ny,ntheta,scheme
integer                                         :: itheta
real*8 , dimension(nx+1)                        :: xz
real*8 , dimension(nx+1,ny+1,ntheta)            :: xadvec,ee,arrin,cgx
real*8                                          :: dxmin_i,dxplus_i,dxcent_i,dt

integer                                         :: istart, iend

xadvec = 0.d0
arrin=ee*cgx


select case (scheme)

case(1)
        do itheta=1,ntheta
                do i=2,nx
                        dxmin_i  = 1.d0/(xz(i)-xz(i-1))
                        dxplus_i = 1.d0/(xz(i+1)-xz(i))
                        dxcent_i = 1.d0/(xz(i+1)-xz(i-1))
                        do j=1,ny+1
                           if (arrin(i,j,itheta)>0) then
                                  xadvec(i,j,itheta)=(arrin(i,j,itheta)-arrin(i-1,j,itheta))*dxmin_i
                           elseif (arrin(i,j,itheta)<0) then
                                  xadvec(i,j,itheta)=(arrin(i+1,j,itheta)-arrin(i,j,itheta))*dxplus_i
                           else
                                  xadvec(i,j,itheta)=(arrin(i+1,j,itheta)-arrin(i-1,j,itheta))*dxcent_i
                           endif
                        end do
                end do
        end do
case(2)
        do itheta=1,ntheta

                if(xmpi_istop) then
                    do j=1,ny+1
                            i=2
                            ! this deserves special attention in the parallel case,
                            ! only meaningful at the top processes row
                            dxmin_i  = 1.d0/(xz(i)-xz(i-1))
                            dxplus_i = 1.d0/(xz(i+1)-xz(i))
                            dxcent_i = 1.d0/(xz(i+1)-xz(i-1))
                            if (arrin(i,j,itheta)>0) then
                                    xadvec(i,j,itheta)=(arrin(i,j,itheta)-arrin(i-1,j,itheta))*dxmin_i
                            elseif (arrin(i,j,itheta)<0) then
                                    xadvec(i,j,itheta)=(arrin(i+1,j,itheta)-arrin(i,j,itheta))*dxplus_i
                            else
                                    xadvec(i,j,itheta)=(arrin(i+1,j,itheta)-arrin(i-1,j,itheta))*dxcent_i
                            endif
                    end do
                endif

                if (xmpi_istop) then
                    istart = 3   ! wwvv take care for parallel case
                else
                    istart = 2
                endif
                if (xmpi_isbot) then
                  iend = nx-1
                else
                  iend = nx
                endif

               ! do i=3,nx-1    ! wwvv the bounds for the serial case
                do i=istart,iend     ! wwvv
                
                        do j=1,ny+1             ! Lax Wendroff
                           xadvec(i,j,itheta)=((arrin(i+1,j,itheta)-arrin(i-1,j,itheta))/(xz(i+1)-xz(i-1)))&
                                                                  -(0.5*dt/((xz(i+1)-xz(i))*(xz(i)-xz(i-1))))*&
                                                                   ((ee(i-1,j,itheta)*cgx(i-1,j,itheta)**2)+&
                                                                        (-2*ee(i,j,itheta)*cgx(i,j,itheta)**2)+&
                                                                        (ee(i+1,j,itheta)*cgx(i+1,j,itheta)**2))
                        end do
                end do

#ifdef USEMPI
                call xmpi_shift(xadvec(:,:,itheta),'1:')  ! wwvv
                call xmpi_shift(xadvec(:,:,itheta),'m:')  ! wwvv
#endif

                if (xmpi_isbot) then
                  do j=1,ny+1
                          i=nx
                          dxmin_i  = 1.d0/(xz(i)-xz(i-1))
                          dxplus_i = 1.d0/(xz(i+1)-xz(i))
                          dxcent_i = 1.d0/(xz(i+1)-xz(i-1))
                          if (arrin(i,j,itheta)>0) then
                                  xadvec(i,j,itheta)=(arrin(i,j,itheta)-arrin(i-1,j,itheta))*dxmin_i
                          elseif (arrin(i,j,itheta)<0) then
                                  xadvec(i,j,itheta)=(arrin(i+1,j,itheta)-arrin(i,j,itheta))*dxplus_i
                          else
                                  xadvec(i,j,itheta)=(arrin(i+1,j,itheta)-arrin(i-1,j,itheta))*dxcent_i
                          endif
                  end do
                endif
        end do
end select

end subroutine advecxho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine advecy(arrin,yadvec,nx,ny,ntheta,yz)

use xmpi_module
IMPLICIT NONE

integer                                         :: i,j,nx,ny,ntheta
integer                                         :: itheta
real*8 ,  dimension(ny+1)                       :: yz
real*8 ,  dimension(nx+1,ny+1,ntheta)           :: yadvec,arrin
real*8                                          :: dymin_i,dyplus_i,dycent_i

yadvec = 0.d0

do itheta=1,ntheta
    do j=2,ny
        dymin_i  = 1.d0/(yz(j)-yz(j-1))
        dyplus_i = 1.d0/(yz(j+1)-yz(j))
        dycent_i = 1.d0/(yz(j+1)-yz(j-1))
        do i=1,nx+1
           if (arrin(i,j,itheta)>0) then
              yadvec(i,j,itheta)=(arrin(i,j,itheta)-arrin(i,j-1,itheta))*dymin_i
           elseif (arrin(i,j,itheta)<0) then
              yadvec(i,j,itheta)=(arrin(i,j+1,itheta)-arrin(i,j,itheta))*dyplus_i
           else
              yadvec(i,j,itheta)=(arrin(i,j+1,itheta)-arrin(i,j-1,itheta))*dycent_i
           endif
        end do
    end do
end do
! lateral boundaries
do itheta=1,ntheta
   do i=1,nx+1
      if (arrin(i,ny+1,itheta)>0) then
             yadvec(i,ny+1,itheta)=(arrin(i,ny+1,itheta)-arrin(i,ny,itheta))*(1.d0/(yz(ny+1)-yz(ny)))
          elseif (arrin(i,ny+1,itheta)<=0) then
             yadvec(i,ny+1,itheta) = yadvec(i,ny,itheta)
      elseif (arrin(i,1,itheta)>=0) then
             yadvec(i,1,itheta)=yadvec(i,2,itheta) 
          elseif (arrin(i,1,itheta)<0) then
         yadvec(i,1,itheta)=(arrin(i,2,itheta)-arrin(i,1,itheta))*(1.d0/(yz(2)-yz(1)))
      endif
        enddo
enddo
                 
!yadvec(:,1,:)=yadvec(:,2,:)             !Ap
!yadvec(:,ny+1,:) = yadvec(:,ny,:)       !Ap

! wwvv: I guess that in the parallel case the first and last columns of yadvec
! have to be copied from neighbours

! It seems that this routine is always called with nx = 0

#ifdef USEMPI
call xmpi_shift(yadvec,':n')  ! fill in yadvec(:,nx+1,:)
call xmpi_shift(yadvec,':1')  ! fill in yadvec(:,1,:)
#endif

end subroutine advecy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine advecyho(ee,cgy,yadvec,nx,ny,ntheta,yz,dt,scheme)
use xmpi_module

IMPLICIT NONE

integer                                         :: i,j,nx,ny,ntheta,scheme
integer                                         :: itheta
real*8 ,  dimension(ny+1)                       :: yz
real*8 ,  dimension(nx+1,ny+1,ntheta)           :: yadvec,ee,arrin,cgy
real*8                                          :: dymin_i,dyplus_i,dycent_i,dt

integer                                         :: jstart, jend

yadvec = 0.d0
arrin=ee*cgy

select case(scheme)

case (1)   ! upwind
        do itheta=1,ntheta
                do j=2,ny
                        dymin_i  = 1.d0/(yz(j)-yz(j-1))
                        dyplus_i = 1.d0/(yz(j+1)-yz(j))
                        dycent_i = 1.d0/(yz(j+1)-yz(j-1))
                        do i=1,nx+1
                           if (arrin(i,j,itheta)>0) then
                                  yadvec(i,j,itheta)=(arrin(i,j,itheta)-arrin(i,j-1,itheta))*dymin_i
                           elseif (arrin(i,j,itheta)<0) then
                                  yadvec(i,j,itheta)=(arrin(i,j+1,itheta)-arrin(i,j,itheta))*dyplus_i
                           else
                                  yadvec(i,j,itheta)=(arrin(i,j+1,itheta)-arrin(i,j-1,itheta))*dycent_i
                           endif
                        end do
                end do
        end do
case(2)    ! Lax Wendroff
        do itheta=1,ntheta
                if(xmpi_isleft) then
                    do i=1,nx+1
                            j=2   ! wwvv only meaningful at left processorcolumn
                            dymin_i  = 1.d0/(yz(j)-yz(j-1))
                            dyplus_i = 1.d0/(yz(j+1)-yz(j))
                            dycent_i = 1.d0/(yz(j+1)-yz(j-1))
                            if (arrin(i,j,itheta)>0) then
                                    yadvec(i,j,itheta)=(arrin(i,j,itheta)-arrin(i,j-1,itheta))*dymin_i
                            elseif (arrin(i,j,itheta)<0) then
                                    yadvec(i,j,itheta)=(arrin(i,j+1,itheta)-arrin(i,j,itheta))*dyplus_i
                            else
                                    yadvec(i,j,itheta)=(arrin(i,j+1,itheta)-arrin(i,j-1,itheta))*dycent_i
                            endif
                    end do
                endif
                if(xmpi_isleft) then
                    jstart = 3
                else
                    jstart = 2
                endif
                if(xmpi_isright) then
                  jend = ny-1
                else
                  jend = ny
                endif
                
                do i=1,nx+1
                        !do j=3,ny-1  ! the serial bounds           ! Lax Wendroff
                        do j=jstart,jend             ! Lax Wendroff     wwvv
                           yadvec(i,j,itheta)=((arrin(i,j+1,itheta)-arrin(i,j-1,itheta))/(yz(j+1)-yz(j-1)))&
                                                                  -(0.5*dt/((yz(j+1)-yz(j))*(yz(j)-yz(j-1))))*&
                                                                   ((ee(i,j-1,itheta)*cgy(i,j-1,itheta)**2)+&
                                                                        (-2*ee(i,j,itheta)*cgy(i,j,itheta)**2)+&
                                                                        (ee(i,j+1,itheta)*cgy(i,j+1,itheta)**2))
                        end do
                end do
#ifdef USEMPI
                call xmpi_shift(yadvec(:,:,itheta),':1')  ! wwvv fill in yadvec(:,1,itheta)
                call xmpi_shift(yadvec(:,:,itheta),':n')  !      fill in yadvec(:,ny+1,itheta
#endif
                
                if (xmpi_isright) then
                  do i=1,nx+1
                          j=ny
                          dymin_i  = 1.d0/(yz(j)-yz(j-1))
                          dyplus_i = 1.d0/(yz(j+1)-yz(j))
                          dycent_i = 1.d0/(yz(j+1)-yz(j-1))
                          if (arrin(i,j,itheta)>0) then
                                  yadvec(i,j,itheta)=(arrin(i,j,itheta)-arrin(i,j-1,itheta))*dymin_i
                          elseif (arrin(i,j,itheta)<0) then
                                  yadvec(i,j,itheta)=(arrin(i,j+1,itheta)-arrin(i,j,itheta))*dyplus_i
                          else
                                  yadvec(i,j,itheta)=(arrin(i,j+1,itheta)-arrin(i,j-1,itheta))*dycent_i
                          endif
                  end do  
                endif
        end do
end select

! lateral boundaries
do itheta=1,ntheta
   do i=1,nx+1
        if (arrin(i,ny+1,itheta)>0) then
          if(xmpi_isright) then
               yadvec(i,ny+1,itheta)=(arrin(i,ny+1,itheta)-arrin(i,ny,itheta))*(1.d0/(yz(ny+1)-yz(ny)))
          endif
        elseif (arrin(i,ny+1,itheta)<=0) then
          if(xmpi_isright) then
               yadvec(i,ny+1,itheta) = yadvec(i,ny,itheta)
          endif
        elseif (arrin(i,1,itheta)>=0) then    ! wwvv is this ok?
                                              ! yadvec(i,1,itheta will not change if 
                                              ! yadvec(i,ny+1,itheta) did change
          if(xmpi_isleft) then
               yadvec(i,1,itheta)=yadvec(i,2,itheta) 
          endif
        elseif (arrin(i,1,itheta)<0) then
          if(xmpi_isleft) then
               yadvec(i,1,itheta)=(arrin(i,2,itheta)-arrin(i,1,itheta))*(1.d0/(yz(2)-yz(1)))
          endif
        endif
   enddo
enddo
                 
end subroutine advecyho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine advectheta(arrin,thetaadvec,nx,ny,ntheta,dtheta)

IMPLICIT NONE

integer                                         :: i,j,nx,ny,ntheta
integer                                         :: itheta
real*8                                          :: dtheta
real*8 ,  dimension(nx+1,ny+1,ntheta)           :: thetaadvec,arrin

thetaadvec = 0

! Ad: include all bins, but use min,max statements
do itheta=1,ntheta
    do j=1,ny+1
        do i=1,nx+1
           if (arrin(i,j,itheta)>0) then
              thetaadvec(i,j,itheta)=(arrin(i,j,itheta)-arrin(i,j,max(itheta-1,1)))/dtheta
           elseif (arrin(i,j,itheta)<0) then
              thetaadvec(i,j,itheta)=(arrin(i,j,min(itheta+1,ntheta))-arrin(i,j,itheta))/dtheta
           else
              thetaadvec(i,j,itheta)=(arrin(i,j,min(itheta+1,ntheta))-arrin(i,j,max(itheta-1,1)))/(2*dtheta)
           endif
        end do
    end do
end do

end subroutine advectheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine advecwx(arrin2d,xwadvec,kmx,nx,ny,xz)
use xmpi_module

IMPLICIT NONE

integer                                         :: i,nx,ny
integer                                         :: j
real*8 , dimension(nx+1)                        :: xz
real*8 , dimension(nx+1,ny+1)                   :: xwadvec,arrin2d,kmx

xwadvec = 0.d0

do j=2,ny
    do i=2,nx   
        if (kmx(i,j)>0) then
           xwadvec(i,j)=(arrin2d(i,j)-arrin2d(i-1,j))/(xz(i)-xz(i-1))
        elseif (kmx(i,j)<0) then
           xwadvec(i,j)=(arrin2d(i+1,j)-arrin2d(i,j))/(xz(i+1)-xz(i))
        else
           xwadvec(i,j)=(arrin2d(i+1,j)-arrin2d(i-1,j))/(xz(i+1)-xz(i-1))
        endif
    end do
end do

! wwvv here we miss the computations of the first and last columns and rows,
!  in the parallel case we shift these in form neighbours
#ifdef USEMPI
call xmpi_shift(xwadvec,'m:') ! fill in xwadvec(nx+1,:)
call xmpi_shift(xwadvec,'1:') ! fill in xwadvec(1,:)
call xmpi_shift(xwadvec,':n') ! fill in xwadvec(:,ny+1)
call xmpi_shift(xwadvec,':1') ! fill in xwadvec(:,1)
#endif

end subroutine advecwx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine advecwy(arrin2d,ywadvec,kmy,nx,ny,yz)
use xmpi_module
use xmpi_module
IMPLICIT NONE

integer                                         :: i,nx,ny
integer                                         :: j
real*8 , dimension(ny+1)                        :: yz
real*8 , dimension(nx+1,ny+1)                   :: ywadvec,arrin2d,kmy

ywadvec = 0.d0

do j=2,ny
    do i=2,nx
        if (kmy(i,j)>0) then
           ywadvec(i,j)=(arrin2d(i,j)-arrin2d(i,j-1))/(yz(j)-yz(j-1))
        elseif (kmy(i,j)<0) then
           ywadvec(i,j)=(arrin2d(i,j+1)-arrin2d(i,j))/(yz(j+1)-yz(j))
        else
           ywadvec(i,j)=(arrin2d(i,j+1)-arrin2d(i,j-1))/(yz(j+1)-yz(j-1))
        endif
    end do
end do

ywadvec(:,1)= ywadvec(:,2)          !Ap
ywadvec(:,ny+1) = ywadvec(:,ny)     !Ap

#ifdef USEMPI
call xmpi_shift(ywadvec,'m:') ! wwvv fill in yadvec(nx+1,:)
call xmpi_shift(ywadvec,'1:') !      fill in yadvec(1,:)
call xmpi_shift(ywadvec,':n') !      fill in yadvec(:,ny+1)
call xmpi_shift(ywadvec,':1') !      fill in yadvec(:,1)
#endif

end subroutine advecwy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine advecqx(c,arrin2d,xwadvec,nx,ny,xz)
use xmpi_module

IMPLICIT NONE

integer                                         :: i,nx,ny
integer                                         :: j
real*8 , dimension(nx+1)                        :: xz
real*8 , dimension(nx+1,ny+1)                   :: xwadvec,arrin2d,c

xwadvec = 0.d0

do j=2,ny
    do i=2,nx   
        if (c(i,j)>0) then
           xwadvec(i,j)=c(i,j)*(arrin2d(i,j)-arrin2d(i-1,j))/(xz(i)-xz(i-1))
        elseif (c(i,j)<0) then
           xwadvec(i,j)=c(i,j)*(arrin2d(i+1,j)-arrin2d(i,j))/(xz(i+1)-xz(i))
        else
           xwadvec(i,j)=c(i,j)*(arrin2d(i+1,j)-arrin2d(i-1,j))/(xz(i+1)-xz(i-1))
        endif
    end do
end do

xwadvec(:,1)= xwadvec(:,2)          !Ap
xwadvec(:,ny+1) = xwadvec(:,ny)     !Ap

#ifdef USEMPI
call xmpi_shift(xwadvec,'m:') ! fill in xwadvec(nx+1,:)
call xmpi_shift(xwadvec,'1:') ! fill in xwadvec(1,:)
call xmpi_shift(xwadvec,':n') ! fill in xwadvec(:,ny+1)
call xmpi_shift(xwadvec,':1') ! fill in xwadvec(:,1)
#endif

end subroutine advecqx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine advecqy(c,arrin2d,ywadvec,nx,ny,yz)
use xmpi_module
use xmpi_module
IMPLICIT NONE

integer                                         :: i,nx,ny
integer                                         :: j
real*8 , dimension(ny+1)                        :: yz
real*8 , dimension(nx+1,ny+1)                   :: ywadvec,arrin2d,c

ywadvec = 0.d0

do j=2,ny
    do i=2,nx
        if (c(i,j)>0) then
           ywadvec(i,j)=c(i,j)*(arrin2d(i,j)-arrin2d(i,j-1))/(yz(j)-yz(j-1))
        elseif (c(i,j)<0) then
           ywadvec(i,j)=c(i,j)*(arrin2d(i,j+1)-arrin2d(i,j))/(yz(j+1)-yz(j))
        else
           ywadvec(i,j)=c(i,j)*(arrin2d(i,j+1)-arrin2d(i,j-1))/(yz(j+1)-yz(j-1))
        endif
    end do
end do

#ifdef USEMPI
call xmpi_shift(ywadvec,'m:') ! wwvv fill in yadvec(nx+1,:)
call xmpi_shift(ywadvec,'1:') !      fill in yadvec(1,:)
call xmpi_shift(ywadvec,':n') !      fill in yadvec(:,ny+1)
call xmpi_shift(ywadvec,':1') !      fill in yadvec(:,1)
#endif

end subroutine advecqy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dispersion(par,s)
use params
use spaceparams

! Robert: iteration along L=L0tanh(2pih/L)

IMPLICIT NONE

type(spacepars)                     :: s
type(parameters)                    :: par

real*8, dimension(s%nx+1,s%ny+1)    :: h,L0,kh
real*8                              :: err,L2
!
! real*8,dimension(:,:),allocatable,save :: L1
! wwvv moved L1 to spaceparams, in the parallel version
! of the program L1 will be resized and distributed
!
integer                             :: i,j
real*8                              :: phi
real*8                              :: aphi
real*8                              :: bphi
real*8                              :: twopi

phi    = (1.0d0 + sqrt(5.0d0))/2
aphi   = 1/(phi+1)
bphi   = phi/(phi+1)
twopi  = 8*atan(1.0d0)
!
! In the original code, phi, aphi, bphi are saved
! variables, and calculated once. I think this is
! better. In the original code, these variables
! had to be broadcasted in the parallel version. 
! Also I rearranged some formulas, no need anymore
! for variables t and n.

! cjaap: replaced par%hmin by par%eps

h = max(s%hh + par%delta*s%H,par%eps)

L0 = twopi*par%g/(s%sigm**2)

if (par%t==0) then
  if (.not. associated(s%L1)) then
    allocate(s%L1(s%nx+1,s%ny+1))
    s%L1=L0
  end if
endif

do j = 1,s%ny+1
  do i = 1,s%nx+1
    err = huge(0.0d0)
    do while (err > 0.00001d0)
      L2        = L0(i,j)*tanh(2*par%px*h(i,j)/s%L1(i,j))
      err       = abs(L2 - s%L1(i,j))
      s%L1(i,j) = (s%L1(i,j)*aphi + L2*bphi)          ! Golden ratio
    end do
  end do
end do

s%k  = 2*par%px/s%L1
s%c  = s%sigm/s%k
!kh   = s%k*h
! Ad:
kh   = min(s%k*h,10.0d0)
s%n=0.5d0+kh/sinh(2*kh)
s%cg=s%c*s%n
!s%cg = s%c*(0.5d0+kh/sinh(2*kh))

end subroutine dispersion


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine breakerdelay(par,s)
use params
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars),target              :: s
type(parameters)                    :: par


real*8                              :: Lbr
real*8, dimension(s%nx+1)           :: utemp
integer                             :: jx,jy,i,nbr,tempxid
integer, dimension(s%nx+1)          :: ibr

include 's.ind'
include 's.inp'


do jy = 2,ny
  usd(1,jy) = ustr(1,jy)
  do jx = 2,nx+1
    nbr=0
    Lbr = sqrt(par%g*hh(jx,jy))*par%Trep
    i=jx-1
    do while (x(i,jy)>=(x(jx,jy)-Lbr).and. i>1)
        nbr = nbr+1
        i=i-1
    end do

    if(nbr.gt.1) then
        do i = 1,nbr+1
            ibr(i) = i
            tempxid = jx-nbr+i-1
            utemp(i) = ustr(tempxid,jy)
        enddo
        ! wwvv consequences for parallel version?
        ! todo, done: see later
        usd(jx,jy) = sum(ibr(1:nbr+1)*utemp(1:nbr+1))/sum(ibr(1:nbr+1)) 
    else
        usd(jx,jy) = ustr(jx,jy)
    end if
  end do
end do

!lateral boundaries
usd(1,:) = usd(2,:)
usd(:,1) = usd(:,2)
usd(:,ny+1) = usd(:,ny)
usd(nx+1,:) = usd(nx,:)

! wwvv for the parallel version, shift in the columns and rows
#ifdef USEMPI
call xmpi_shift(usd,'m:')  ! fill in usd(nx+1,:)
call xmpi_shift(usd,'1:')  ! fill in usd(1,:)
call xmpi_shift(usd,':n')  ! fill in usd(:,ny+1)
call xmpi_shift(usd,':1')  ! fill in usd(:,1)

#endif


end subroutine breakerdelay 
end module wave_timestep_module
