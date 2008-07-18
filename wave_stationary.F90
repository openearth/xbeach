module wave_stationary_module
contains
subroutine wave_stationary(s,par)
use params
use spaceparams
use roelvink_module
use wave_timestep_module
use xmpi_module

IMPLICIT NONE

type(spacepars), target     :: s
type(parameters)            :: par

integer                     :: i
integer                     :: j
integer                     :: itheta,iter
integer, dimension(:,:,:),allocatable,save  :: wete
real*8 , dimension(:,:)  ,allocatable,save  :: dhdx,dhdy,dudx,dudy,dvdx,dvdy,ustw
real*8 , dimension(:,:)  ,allocatable,save  :: km,kmx,kmy,wm,xwadvec,ywadvec,sinh2kh
real*8 , dimension(:,:,:),allocatable,save  :: xadvec,yadvec,thetaadvec,dd,drr
real*8 , dimension(:,:,:),allocatable,save  :: xradvec,yradvec,thetaradvec
real*8 , dimension(:),allocatable,save      :: Hprev
real*8                                      :: Herr,dtw
real*8 , dimension(:,:)  ,allocatable,save  :: dkmxdx,dkmxdy,dkmydx,dkmydy,cgxm,cgym

include 's.ind'
include 's.inp'

if (.not. allocated(wete)) then
   allocate(wete      (nx+1,ny+1,ntheta))
   allocate(xadvec    (nx+1,ny+1,ntheta))
   allocate(yadvec    (nx+1,ny+1,ntheta))
   allocate(thetaadvec(nx+1,ny+1,ntheta))
   allocate(xradvec    (nx+1,ny+1,ntheta))
   allocate(yradvec    (nx+1,ny+1,ntheta))
   allocate(thetaradvec(nx+1,ny+1,ntheta))
   allocate(dd        (nx+1,ny+1,ntheta))
   allocate(drr       (nx+1,ny+1,ntheta))
   allocate(dhdx(nx+1,ny+1))
   allocate(dhdy(nx+1,ny+1))
   allocate(dudx(nx+1,ny+1))
   allocate(dudy(nx+1,ny+1))
   allocate(dvdx(nx+1,ny+1))
   allocate(dvdy(nx+1,ny+1))
   allocate(km  (nx+1,ny+1))
   allocate(kmx (nx+1,ny+1))
   allocate(kmy (nx+1,ny+1))
   allocate(wm  (nx+1,ny+1))
   allocate(ustw(nx+1,ny+1))
   allocate(xwadvec(nx+1,ny+1))
   allocate(ywadvec(nx+1,ny+1))
   allocate(sinh2kh(nx+1,ny+1))
   allocate(Hprev(ny+1))

   allocate(dkmxdx  (nx+1,ny+1))
   allocate(dkmxdy  (nx+1,ny+1))
   allocate(dkmydx  (nx+1,ny+1))
   allocate(dkmydy  (nx+1,ny+1))
   allocate(cgxm    (nx+1,ny+1))
   allocate(cgym    (nx+1,ny+1))

endif

wete        = 0.0d0
xadvec      = 0.0d0
yadvec      = 0.0d0
thetaadvec  = 0.0d0
xradvec     = 0.0d0
yradvec     = 0.0d0
thetaradvec = 0.0d0
dd          = 0.0d0
drr         = 0.0d0
dhdx        = 0.0d0
dhdy        = 0.0d0
dudx        = 0.0d0
dudy        = 0.0d0
dvdx        = 0.0d0
dvdy        = 0.0d0
km          = 0.0d0
kmx         = 0.0d0
kmy         = 0.0d0
wm          = 0.0d0
ustw        = 0.0d0
xwadvec     = 0.0d0
ywadvec     = 0.0d0
sinh2kh     = 0.0d0

dkmxdx      = 0.0d0
dkmxdy      = 0.0d0
dkmydx      = 0.0d0
dkmydy      = 0.0d0
cgxm        = 0.0d0
cgym        = 0.0d0


! cjaap: replaced par%hmin by par%eps
hh = max(hh,par%eps)

!
! Dispersion relation
!
sigm = max((sum(sigt,3)/ntheta),0.01d0)
call dispersion(par,s)

! Slopes of water depth
call slope2D(hh+par%delta*H,nx,ny,xz,yz,dhdx,dhdy)
call slope2D(u*par%wci,nx,ny,xz,yz,dudx,dudy)
call slope2D(v*par%wci,nx,ny,xz,yz,dvdx,dvdy)

!
! Propagation speeds in x,y and theta space
do j=1,ny+1
    do i=1,nx+1
        if(2.d0*hh(i,j)*k(i,j)>=3000.d0) then
#ifdef USEMPI
            write(*,*) 'Process:', xmpi_rank ,i ,j, par%t, hh(i,j), k(i,j)
#else
            write(*,*) i ,j, par%t, hh(i,j), k(i,j)
#endif
        endif
        sinh2kh(i,j)=sinh(2*k(i,j)*(hh(i,j)+par%delta*H(i,j)))
    end do
end do
DO itheta=1,ntheta
    cgx(:,:,itheta)=cg*cxsth(itheta)+uu*par%wci*min(hh/par%hwci,1.d0)
    cgy(:,:,itheta)=cg*sxnth(itheta)+vv*par%wci*min(hh/par%hwci,1.d0)
    cx(:,:,itheta) =c*cxsth(itheta)+uu*par%wci*min(hh/par%hwci,1.d0)
    cy(:,:,itheta) =c*sxnth(itheta)+vv*par%wci*min(hh/par%hwci,1.d0)
    ctheta(:,:,itheta)= &
    sigm/sinh2kh*(dhdx*&
    sxnth(itheta)-dhdy*cxsth(itheta)) + &
    cxsth(itheta)*(sxnth(itheta)*dudx - cxsth(itheta)*dudy) + &
    sxnth(itheta)*(sxnth(itheta)*dvdx - cxsth(itheta)*dvdy)
END DO

    thetamean(1,:)=(sum(ee(1,:,:)*thet(1,:,:),2)/size(ee(1,:,:),2)) &
                  /(max(sum(ee(1,:,:),2),0.00001d0) /size(ee(1,:,:),2))
 if (par%wci==1) then
    tm = (sum(ee*thet,3)/ntheta)/(max(sum(ee,3),0.00001d0)/ntheta)
    km = k
    kmx = km*cos(tm)
    kmy = km*sin(tm)
    wm = sigm+kmx*uu*par%wci*min(hh/par%hwci,1.d0)+kmy*vv*par%wci*min(hh/par%hwci,1.d0)
        cgym = cg*sin(tm) + vv*min(hh/par%hwci,1.d0)
        cgxm = cg*cos(tm) + uu*min(hh/par%hwci,1.d0)

    call slope2D(kmx,nx,ny,xz,yz,dkmxdx,dkmxdy)
    call slope2D(kmy,nx,ny,xz,yz,dkmydx,dkmydy)
    call advecwx(wm,xwadvec,nx,ny,xz)   ! cjaap: xz or xu?

    kmx = kmx -par%dt*xwadvec -par%dt*cgym*(dkmydx-dkmxdy)
    kmx(:,ny+1) = kmx(:,ny)  ! lateral bc

    call advecwy(wm,ywadvec,nx,ny,yz)   ! cjaap: yz or yv?
    kmy = kmy-par%dt*ywadvec + par%dt*cgxm*(dkmydx-dkmxdy)
    kmy(:,ny+1) = kmy(:,ny)   ! lateral bc
    km = sqrt(kmx**2+kmy**2)
else
     km = k
endif

    E(1,:)=sum(ee(1,:,:),2)*dtheta
    R(1,:)=max(sum(rr(1,:,:),2)*dtheta,0.0d0)
    H(1,:)=sqrt(E(1,:)/par%rhog8)
                   
do i=2,nx
  dtw=.9*minval(xz(2:nx+1)-xz(1:nx))/sqrt(par%g*maxval(hh(i,:)))
  Herr=1.
  iter=0
  H(i,:)=H(i-1,:)
    ee(i,:,:)=ee(i-1,:,:)
    ee(i+1,:,:)=ee(i,:,:)

  do while(Herr>1.d-5.and. iter<10)
    iter=iter+1
    Hprev=H(i,:)
    !
    ! transform to wave action
    !
    ee(i-1:i+1,:,:) = ee(i-1:i+1,:,:)/sigt(i-1:i+1,:,:)
    !
    ! Upwind Euler timestep propagation
    !
    call advecx(ee(i-1:i+1,:,:)*cgx(i-1:i+1,:,:),xadvec(i-1:i+1,:,:),2,ny,ntheta,xz(i-1:i+1)) !Robert & Jaap
    call advecy(ee(i,:,:)*cgy(i,:,:),yadvec(i,:,:),0,ny,ntheta,yz)                   !Robert & Jaap, nx and ny increased with 1 in advecy
    call advectheta(ee(i,:,:)*ctheta(i,:,:),thetaadvec(i,:,:),0,ny,ntheta,dtheta)    !Robert & Jaap
    
    ee(i,:,:)=ee(i,:,:)-dtw*(xadvec(i,:,:) &
                            +yadvec(i,:,:) &
                            +thetaadvec(i,:,:))
    !
    ! transform back to wave energy
    !
    ee(i-1:i+1,:,:) = ee(i-1:i+1,:,:)*sigt(i-1:i+1,:,:)
    ee(i,:,:)=max(ee(i,:,:),0.0d0)

        
    !
    ! Energy integrated over wave directions,Hrms
    !
    E(i,:)=sum(ee(i,:,:),2)*dtheta
    H(i,:)=sqrt(E(i,:)/par%rhog8)
    do itheta=1,ntheta
       ee(i,:,itheta)=ee(i,:,itheta)/max(1.,(H(i,:)/(par%gammax*hh(i,:)))**2)
    enddo
    H(i,:)=min(H(i,:),par%gammax*hh(i,:))
    E(i,:)=par%rhog8*H(i,:)**2
    !write(*,*) 'sqrt in wave_timestep after H (line 97)'
    
    !
    ! calculate change in intrinsic frequency
    tm(i,:) = (sum(ee(i,:,:)*thet(i,:,:),2)/ntheta)/(max(sum(ee(i,:,:),2),0.00001)/ntheta)
    if (par%wci/=1) then
    km(i,:) = k(i,:)
    endif

    sigm(i,:) = sqrt(par%g*km(i,:)*tanh(km(i,:)*(hh(i,:)+par%delta*H(i,:))))

    DO itheta=1,ntheta
         sigt(i,:,itheta) = max(sigm(i,:),0.01)
    END DO
    !
    ! Total dissipation
    if(par%break==1 .or. par%break==3) THEN
        call roelvink1(E(i,:),hh(i,:), &
                       par%Trep,par%alpha,par%gamma,par%n, &
                       par%rho,par%g,par%delta,D(i,:),ny+1,par%break)
    else if(par%break==2) THEN
        call baldock1(E(i,:),hh(i,:),k(i,:), &
                       par%Trep,par%alpha,par%gamma, &
                       par%rho,par%g,par%delta,D(i,:),ny+1)
    end if
    !
    ! Distribution of dissipation over directions and frequencies
    !
    do itheta=1,ntheta
       dd(i,:,itheta)=ee(i,:,itheta)*D(i,:)/max(E(i,:),0.00001d0)
    end do
    do j=1,ny+1
        ! cjaap: replaced par%hmin by par%eps
        if(hh(i,j)+par%delta*H(i,j)>par%eps) then
            wete(i,j,1:ntheta)=1
        else
            wete(i,j,1:ntheta)=0
        end if
    end do
    !
    ! Euler step dissipation
    !
    ! calculate roller energy balance
    !
    call advecx(rr(i-1:i+1,:,:)*cx(i-1:i+1,:,:),xradvec(i-1:i+1,:,:),2,ny,ntheta,xz(i-1:i+1)) !Robert & Jaap
    call advecy(rr(i,:,:)*cy(i,:,:),yradvec(i,:,:),0,ny,ntheta,yz)                   !Robert & Jaap
    call advectheta(rr(i,:,:)*ctheta(i,:,:),thetaradvec(i,:,:),0,ny,ntheta,dtheta)   !Robert & Jaap
    
    rr(i,:,:)=rr(i,:,:)-dtw*(xradvec(i,:,:) &
                               +yradvec(i,:,:) &
                               +thetaradvec(i,:,:))
    rr(i,:,:)=max(rr(i,:,:),0.0d0)
    !
    ! euler step roller energy dissipation (source and sink function)
     do j=1,ny+1
        do itheta=1,ntheta
!            if(wete(i,j,itheta)==1) then
!                ee(i,j,itheta)=ee(i,j,itheta)-dtw*dd(i,j,itheta)
!                rr(i,j,itheta)=rr(i,j,itheta)+dtw*dd(i,j,itheta)&
!                              -dtw*2.*par%g*par%beta*rr(i,j,itheta)&
!                              /sqrt(cx(i,j,itheta)**2.+cy(i,j,itheta)**2.);
!                drr(i,j,itheta) = 2.*par%g*par%beta*max(rr(i,j,itheta),0.0d0)/           &
!                                  sqrt(cx(i,j,itheta)**2 +cy(i,j,itheta)**2.);
!            else if (par%roller==0) then
!                rr(i,j,itheta)= 0.0
!            end if
!            ee(i,j,itheta)=max(ee(i,j,itheta),0.0d0)
!            rr(i,j,itheta)=max(rr(i,j,itheta),0.0d0);
!            if(wete(i,j,itheta)==0) then
!                ee(i,j,itheta)=0.0d0
!                rr(i,j,itheta)=0.0d0;
!            end if
                        
              if(wete(i,j,itheta)==1) then
                ee(i,j,itheta)=ee(i,j,itheta)-dtw*dd(i,j,itheta)
                if (par%roller==1) then  !Christophe
                rr(i,j,itheta)=rr(i,j,itheta)+dtw*dd(i,j,itheta)&
                              -dtw*2.*par%g*par%beta*rr(i,j,itheta)&
                              /sqrt(cx(i,j,itheta)**2+cy(i,j,itheta)**2)
                drr(i,j,itheta) = 2*par%g*par%beta*max(rr(i,j,itheta),0.0d0)/           &
                                  sqrt(cx(i,j,itheta)**2 +cy(i,j,itheta)**2)
            else if (par%roller==0) then
                rr(i,j,itheta)= 0.0d0
                drr(i,j,itheta)= 0.0d0
            end if
            ee(i,j,itheta)=max(ee(i,j,itheta),0.0d0)
            rr(i,j,itheta)=max(rr(i,j,itheta),0.0d0)
            else if(wete(i,j,itheta)==0) then
                ee(i,j,itheta)=0.0d0
                rr(i,j,itheta)=0.0d0
                drr(i,j,itheta)=0.0d0
            end if
        end do
    end do
    ! Lateral boundary condition
    do itheta=1,ntheta 
       if (theta(itheta)>=0.) then
          ee(i,1,itheta)=ee(i,2,itheta)
          rr(i,1,itheta)=rr(i,2,itheta)
       endif
       if (theta(itheta)<=0.) then
          ee(i,ny+1,itheta)=ee(i,ny,itheta)
          rr(i,ny+1,itheta)=rr(i,ny,itheta)
       endif
    end do
    
    
    !
    ! Compute mean wave direction
    !
    thetamean(i,:)=(sum(ee(i,:,:)*thet(i,:,:),2)/size(ee(i,:,:),2)) &
                  /(max(sum(ee(i,:,:),2),0.00001) /size(ee(i,:,:),2))
    !
    ! Energy integrated over wave directions,Hrms
    !
    E(i,:)=sum(ee(i,:,:),2)*dtheta
    R(i,:)=sum(rr(i,:,:),2)*dtheta
    DR(i,:)=sum(drr(i,:,:),2)*dtheta
    H(i,:)=sqrt(E(i,:)/par%rhog8)
    Herr=maxval(abs(Hprev-H(i,:)))
  enddo
  if(xreader) then
    write(*,*)i,iter,Herr
  endif
enddo
!
! Radiation stresses and forcing terms
!
n=cg/c
Sxx=(n*sum((1.d0+(cos(thet))**2.d0)*ee,3)-.5d0*sum(ee,3))*dtheta
Syy=(n*sum((1.d0+(sin(thet))**2.d0)*ee,3)-.5d0*sum(ee,3))*dtheta
Sxy=n*sum(sin(thet)*cos(thet)*ee,3)*dtheta

! add roller contribution

Sxx = Sxx + sum((cos(thet)**2)*rr,3)*dtheta
Syy = Syy + sum((sin(thet)**2)*rr,3)*dtheta
Sxy = Sxy + sum(sin(thet)*cos(thet)*rr,3)*dtheta

do j=2,ny
    do i=1,nx
       Fx(i,j)=-(Sxx(i+1,j)-Sxx(i,j))/(xz(i+1)-xz(i))                                   &
               -0.5d0*(Sxy(i,j+1)+Sxy(i+1,j+1)- Sxy(i,j-1)-Sxy(i+1,j-1))/(yz(j+1)-yz(j-1))
    enddo
enddo

do j=1,ny
    do i=2,nx
       Fy(i,j)=-(Syy(i,j+1)-Syy(i,j))/(yz(j+1)-yz(j))                                  &
               -0.5d0*(Sxy(i+1,j)+Sxy(i+1,j+1)-Sxy(i-1,j)-Sxy(i-1,j+1))/(xz(i+1)-xz(i-1))
    enddo
enddo

Fx(:,1)=Fx(:,2)
Fy(:,1)=Fy(:,2)
Fx(:,ny+1)=Fx(:,ny+1-1)
Fy(:,ny+1)=Fy(:,ny+1-1)
Fx(1,:)=Fx(2,:)
Fy(1,:)=Fy(2,:)
urms=par%px*H/par%Trep/(sqrt(2.d0)*sinh(k*(hh+par%delta*H)))
!ust=E*k/sigm/par%rho/max(hh,0.001)

! wave induced mass flux
ustw=E*k/sigm/par%rho/max(hh,.001d0)
uwf = ustw*cos(tm)
vwf = ustw*sin(tm)
! roller contribution
ustr=2.*R*k/sigm/par%rho/max(hh,.001d0)
! introduce breaker delay
call breakerdelay(par,s)
!ust = usd
ust=usd+ustw
!lateral boundaries
ust(1,:) = ust(2,:)
ust(:,1) = ust(:,2)
ust(:,ny+1) = ust(:,ny)
D=2*par%g*par%beta*sum(rr/sqrt(cx**2+cy**2),3)*dtheta

end subroutine wave_stationary
end module wave_stationary_module
