module wave_functions_module

contains

subroutine slope2D(h,nx,ny,dsu,dnv,dhdx,dhdy)
use xmpi_module
IMPLICIT NONE

integer                           :: i,j,nx,ny
real*8, dimension(nx+1,ny+1)      :: h,dhdx,dhdy
real*8, dimension(nx+1,ny+1)           :: dsu
real*8, dimension(nx+1,ny+1)           :: dnv


! wwvv dhdx(2:nx,:) is computed, dhdx(1,:) and dhdx(nx+1,:) 
! get boundary values, so in the parallel case, we need
! to do something about that: get the boundaries from
! upper and lower neighbours

do j=1,ny+1
    if(nx+1>=2)then
      do i=2,nx
          dhdx(i,j)=(h(i+1,j)-h(i-1,j))/(dsu(i,j)+dsu(i-1,j))
      end do  
      dhdx(1,j)=(h(2,j)-h(1,j))/dsu(1,j)
      dhdx(nx+1,j)=(h(nx+1,j)-h(nx,j))/dsu(nx,j)
    end if
end do
#ifdef USEMPI
call xmpi_shift(dhdx,'m:')  ! fill in dhdx(nx+1,:)
call xmpi_shift(dhdx,'1:')  ! fill in dhdx(1,:)
#endif

do i=1,nx+1
    if(ny+1>=2)then
      do j=2,ny
          dhdy(i,j)=(h(i,j+1)-h(i,j-1))/(dnv(i,j)+dnv(i,j-1))
      end do
      dhdy(i,1)=(h(i,2)-h(i,1))/dnv(i,1)
      dhdy(i,ny+1)=(h(i,ny+1)-h(i,ny))/dnv(i,ny)
    end if
end do

#ifdef USEMPI
call xmpi_shift(dhdy,':n')  !  fill in dhdy(:,ny+1)
call xmpi_shift(dhdy,':1')  !  fill in dhdy(:,1)
#endif
end subroutine slope2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine advecxho(ee,cgx,xadvec,nx,ny,ntheta,dnu,dsu,dsdnzi,dt,scheme)

  IMPLICIT NONE

  integer                                         :: i,j,nx,ny,ntheta
  character(len=*), intent(in)                    :: scheme
  character(len=10)                               :: scheme_now
  integer                                         :: itheta
  real*8 , dimension(nx+1,ny+1)                   :: dnu,dsu,dsdnzi,fluxx
  real*8 , dimension(nx+1,ny+1,ntheta)            :: xadvec,ee,cgx
  real*8                                          :: dt,cgxu,eupw

  xadvec = 0.d0
  fluxx  = 0.d0


  ! split into schemes first, less split loops -> more efficiency 
  scheme_now=scheme
  select case(trim(scheme_now))
     case('upwind_1')
        do itheta=1,ntheta
           do j=1,ny+1
              do i=1,nx  ! Whole domain
                 cgxu=.5*(cgx(i+1,j,itheta)+cgx(i,j,itheta))
                 if (cgxu>0) then
                    fluxx(i,j)=ee(i,j,itheta)*cgxu*dnu(i,j)
                 else
                    fluxx(i,j)=ee(i+1,j,itheta)*cgxu*dnu(i,j)
                 endif
              enddo
           enddo
           do j=1,ny+1  ! 
              do i=2,nx 
                 xadvec(i,j,itheta)=(fluxx(i,j)-fluxx(i-1,j))*dsdnzi(i,j)
              enddo
           enddo
        enddo
     case('upwind_2')
        do itheta=1,ntheta
           do j=1,ny+1
              do i=2,nx-1
                 cgxu=.5*(cgx(i+1,j,itheta)+cgx(i,j,itheta))
                 if (cgxu>0) then
!                    eupw=((dsu(i,j)+.5*dsu(i-1,j))*ee(i,j,itheta)-.5*dsu(i-1,j)*ee(i-1,j,itheta))/dsu(i-1,j)
                     eupw=((dsu(i-1,j)+.5*dsu(i,j))*ee(i,j,itheta)-.5*dsu(i,j)*ee(i-1,j,itheta))/dsu(i-1,j)
                     if (eupw<0.d0) eupw=ee(i,j,itheta)
                     fluxx(i,j)=eupw*cgxu*dnu(i,j)
                 else
!                    eupw=((dsu(i+1,j)+.5*dsu(i+2,j))*ee(i+1,j,itheta)-.5*dsu(i+2,j)*ee(i+2,j,itheta))/dsu(i+1,j)
                     eupw=((dsu(i+1,j)+.5*dsu(i,j))*ee(i+1,j,itheta)-.5*dsu(i,j)*ee(i+2,j,itheta))/dsu(i+1,j)
                     if (eupw<0.d0) eupw=ee(i+1,j,itheta)
                     fluxx(i,j)=eupw*cgxu*dnu(i,j)
                 endif
              enddo
              i=1   ! only compute for i==1
              cgxu=.5*(cgx(i+1,j,itheta)+cgx(i,j,itheta))
              if (cgxu>0) then
                 fluxx(i,j)=ee(i,j,itheta)*cgxu*dnu(i,j)
              else
!                     eupw=((dsu(i+1,j)+.5*dsu(i+2,j))*ee(i+1,j,itheta)-.5*dsu(i+2,j)*ee(i+2,j,itheta))/dsu(i+1,j)
                      eupw=((dsu(i+1,j)+.5*dsu(i,j))*ee(i+1,j,itheta)-.5*dsu(i,j)*ee(i+2,j,itheta))/dsu(i+1,j)
                 if (eupw<0.d0) eupw=ee(i+1,j,itheta)
                 fluxx(i,j)=eupw*cgxu*dnu(i,j)
              endif
              i=nx  ! only compute for i==nx
              cgxu=.5*(cgx(i+1,j,itheta)+cgx(i,j,itheta))
              if (cgxu>0) then
!                    eupw=((dsu(i,j)+.5*dsu(i-1,j))*ee(i,j,itheta)-.5*dsu(i-1,j)*ee(i-1,j,itheta))/dsu(i-1,j)
                     eupw=((dsu(i-1,j)+.5*dsu(i,j))*ee(i,j,itheta)-.5*dsu(i,j)*ee(i-1,j,itheta))/dsu(i-1,j)
                 if (eupw<0.d0) eupw=ee(i,j,itheta)
                 fluxx(i,j)=eupw*cgxu*dnu(i,j)
              else
                 fluxx(i,j)=ee(i+1,j,itheta)*cgxu*dnu(i,j)
              endif
           enddo
           do j=1,ny+1  ! 
              do i=2,nx 
                 xadvec(i,j,itheta)=(fluxx(i,j)-fluxx(i-1,j))*dsdnzi(i,j)
              enddo
           enddo
        enddo
  end select

end subroutine advecxho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine advecthetaho(ee,ctheta,thetaadvec,nx,ny,ntheta,dtheta,scheme)

  IMPLICIT NONE

  integer                                         :: i,j,nx,ny,ntheta
  character(len=*), intent(in)                    :: scheme
  character(len=10)                               :: scheme_now
  integer                                         :: itheta
  real*8 , dimension(ntheta)                      :: fluxtheta
  real*8 , dimension(nx+1,ny+1,ntheta)            :: thetaadvec,ee,ctheta
  real*8                                          :: dtheta,ctheta_between,eupw

  thetaadvec = 0.d0
  fluxtheta  = 0.d0

  ! No refraction caan take place if ntheta==1 
  if (ntheta>1) then

  ! split into schemes first, less split loops -> more efficiency 
  scheme_now=scheme
  select case(trim(scheme_now))
     case('upwind_1')
        do j=1,ny+1
           do i=1,nx+1  
              do itheta=1,ntheta-1
                 ctheta_between=.5*(ctheta(i,j,itheta+1)+ctheta(i,j,itheta))
                 if (ctheta_between>0) then
                    fluxtheta(itheta)=ee(i,j,itheta)*ctheta_between
                 else
                    fluxtheta(itheta)=ee(i,j,itheta+1)*ctheta_between
                 endif
              enddo
              thetaadvec(i,j,1)=(fluxtheta(1)-0.d0)/dtheta ! No flux across lower boundary theta grid
              do itheta=2,ntheta-1
                 thetaadvec(i,j,itheta)=(fluxtheta(itheta)-fluxtheta(itheta-1))/dtheta
              enddo
              thetaadvec(i,j,ntheta)=(0.d0-fluxtheta(ntheta-1))/dtheta ! No flux across upper boundary theta grid
           enddo
        enddo
     case('upwind_2')
        do j=1,ny+1
           do i=1,nx+1  
              do itheta=2,ntheta-2
                 ctheta_between=.5*(ctheta(i,j,itheta+1)+ctheta(i,j,itheta))
                 if (ctheta_between>0) then
                     eupw=1.5d0*ee(i,j,itheta)-.5*ee(i,j,itheta-1)
                     if (eupw<0.d0) eupw=ee(i,j,itheta)
                     fluxtheta(itheta)=eupw*ctheta_between
                 else
                     eupw=1.5d0*ee(i,j,itheta+1)-.5*ee(i,j,itheta+2)
                     if (eupw<0.d0) eupw=ee(i,j,itheta+1)
                     fluxtheta(itheta)=eupw*ctheta_between
                 endif
              enddo
              itheta=1   ! only compute for itheta==1
                 ctheta_between=.5*(ctheta(i,j,itheta+1)+ctheta(i,j,itheta))
                 if (ctheta_between>0) then
                    fluxtheta(itheta)=ee(i,j,itheta)*ctheta_between
                 else
                     eupw=1.5d0*ee(i,j,itheta+1)-.5*ee(i,j,itheta+2)
                     if (eupw<0.d0) eupw=ee(i,j,itheta+1)
                     fluxtheta(itheta)=eupw*ctheta_between
                 endif
              itheta=ntheta-1  ! only compute for itheta==ntheta-1
                 ctheta_between=.5*(ctheta(i,j,itheta+1)+ctheta(i,j,itheta))
                 if (ctheta_between>0) then
                     eupw=1.5d0*ee(i,j,itheta)-.5*ee(i,j,itheta-1)
                     if (eupw<0.d0) eupw=ee(i,j,itheta)
                     fluxtheta(itheta)=eupw*ctheta_between
                 else
                     eupw=ee(i,j,itheta+1)
                     fluxtheta(itheta)=eupw*ctheta_between
                 endif
              thetaadvec(i,j,1)=(fluxtheta(1)-0.d0)/dtheta ! No flux across lower boundary theta grid
              do itheta=2,ntheta-1
                 thetaadvec(i,j,itheta)=(fluxtheta(itheta)-fluxtheta(itheta-1))/dtheta
              enddo
              thetaadvec(i,j,ntheta)=(0.d0-fluxtheta(ntheta-1))/dtheta ! No flux across upper boundary theta grid
           enddo
        enddo
  end select
  
  endif !ntheta>1
end subroutine advecthetaho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine advecyho(ee,cgy,yadvec,nx,ny,ntheta,dsv,dnv,dsdnzi,dt,scheme)

  IMPLICIT NONE

  integer                                         :: i,j,nx,ny,ntheta
  character(len=*), intent(in)                    :: scheme
  character(len=10)                               :: scheme_now
  integer                                         :: itheta
  real*8 ,  dimension(nx+1,ny+1)                  :: dsv,dnv,dsdnzi,fluxy
  real*8 ,  dimension(nx+1,ny+1,ntheta)           :: yadvec,ee,cgy
  real*8                                          :: dt,cgyv,eupw

  yadvec = 0.d0
  fluxy  = 0.d0

  ! split into schemes first, less split loops -> more efficiency 
  scheme_now=scheme
  select case(trim(scheme_now))
     case('upwind_1')
        do itheta=1,ntheta
           do j=1,ny  
              do i=1,nx+1  ! Whole domain
                 cgyv=.5*(cgy(i,j+1,itheta)+cgy(i,j,itheta))
                 if (cgyv>0) then
                    fluxy(i,j)=ee(i,j,itheta)*cgyv*dsv(i,j)
                 else
                    fluxy(i,j)=ee(i,j+1,itheta)*cgyv*dsv(i,j)
                 endif
              enddo
           enddo
           do j=2,ny
              do i=1,nx+1 
                 yadvec(i,j,itheta)=(fluxy(i,j)-fluxy(i,j-1))*dsdnzi(i,j)
              enddo
           enddo
        enddo
     case('upwind_2')
        do itheta=1,ntheta
           do j=2,ny-1  
              do i=1,nx+1
                 cgyv=.5*(cgy(i,j+1,itheta)+cgy(i,j,itheta))
                 if (cgyv>0) then
!                    eupw=((dnv(i,j)+.5*dnv(i,j-1))*ee(i,j,itheta)-.5*dnv(i,j-1)*ee(i,j-1,itheta))/dnv(i,j-1)
                     eupw=((dnv(i,j-1)+.5*dnv(i,j))*ee(i,j,itheta)-.5*dnv(i,j)*ee(i,j-1,itheta))/dnv(i,j-1)
                    if (eupw<0.d0) eupw=ee(i,j,itheta)
                    fluxy(i,j)=eupw*cgyv*dsv(i,j)
                 else
!                   eupw=((dnv(i,j+1)+.5*dnv(i,j+2))*ee(i,j+1,itheta)-.5*dnv(i,j+2)*ee(i,j+2,itheta))/dnv(i,j+1)
                    eupw=((dnv(i,j+1)+.5*dnv(i,j))*ee(i,j+1,itheta)-.5*dnv(i,j)*ee(i,j+2,itheta))/dnv(i,j+1)
                    if (eupw<0.d0) eupw=ee(i,j+1,itheta)
                    fluxy(i,j)=eupw*cgyv*dsv(i,j)
                 endif
              enddo
           enddo
           j=1   ! only compute for j==1
           do i=1,nx+1
              cgyv=.5*(cgy(i,j+1,itheta)+cgy(i,j,itheta))
              if (cgyv>0) then
                 fluxy(i,j)=ee(i,j,itheta)*cgyv*dsv(i,j)
              else
!                   eupw=((dnv(i,j+1)+.5*dnv(i,j+2))*ee(i,j+1,itheta)-.5*dnv(i,j+2)*ee(i,j+2,itheta))/dnv(i,j+1)
                    eupw=((dnv(i,j+1)+.5*dnv(i,j))*ee(i,j+1,itheta)-.5*dnv(i,j)*ee(i,j+2,itheta))/dnv(i,j+1)
                 if (eupw<0.d0) eupw=ee(i,j+1,itheta)
                 fluxy(i,j)=eupw*cgyv*dsv(i,j)
              endif
           enddo
           j=ny ! only compute for j==ny
           do i=1,nx+1
              cgyv=.5*(cgy(i,j+1,itheta)+cgy(i,j,itheta))
              if (cgyv>0) then
!                eupw=((dnv(i,j)+.5*dnv(i,j-1))*ee(i,j,itheta)-.5*dnv(i,j-1)*ee(i,j-1,itheta))/dnv(i,j-1)
                 eupw=((dnv(i,j-1)+.5*dnv(i,j))*ee(i,j,itheta)-.5*dnv(i,j)*ee(i,j-1,itheta))/dnv(i,j-1)
                 if (eupw<0.d0) eupw=ee(i,j,itheta)
                 fluxy(i,j)=eupw*cgyv*dsv(i,j)
              else
                 fluxy(i,j)=ee(i,j+1,itheta)*cgyv*dsv(i,j)
              endif
           enddo
           do j=2,ny
              do i=2,nx+1 
                 yadvec(i,j,itheta)=(fluxy(i,j)-fluxy(i,j-1))*dsdnzi(i,j)
              enddo
           enddo
        enddo
  end select

end subroutine advecyho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine advectheta(arrin,thetaadvec,nx,ny,ntheta,dtheta)

IMPLICIT NONE

integer                                         :: i,j,nx,ny,ntheta
integer                                         :: itheta
real*8                                          :: dtheta
real*8 ,  dimension(nx+1,ny+1,ntheta)           :: thetaadvec,arrin

thetaadvec = 0

! Dano: make sure no refraction is computed for ntheta==1
if (ntheta>1) then

! Ad: include all bins, but use min,max statements
do itheta=1,ntheta
    do j=1,ny+1
        do i=2,nx+1
           if (arrin(i,j,itheta)>0) then
              if (itheta==1) then
                 ! Dano: arrin (i,j,theta-1) is over the edge of the theta grid; energy there is zero
                 thetaadvec(i,j,itheta)=arrin(i,j,itheta)/dtheta               
              else
                 thetaadvec(i,j,itheta)=(arrin(i,j,itheta)-arrin(i,j,itheta-1))/dtheta
              endif
           elseif (arrin(i,j,itheta)<0) then
              if (itheta==ntheta) then
                 ! Dano: arrin (i,j,theta+1) is over the edge of the theta grid; energy there is zero
                 thetaadvec(i,j,itheta)=(-arrin(i,j,itheta))/dtheta
              else
                 thetaadvec(i,j,itheta)=(arrin(i,j,itheta+1)-arrin(i,j,itheta))/dtheta
              endif              
           else
              thetaadvec(i,j,itheta)=(arrin(i,j,min(itheta+1,ntheta))-arrin(i,j,max(itheta-1,1)))/(2*dtheta)
           endif
        end do
    end do
end do

endif

end subroutine advectheta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine advecwx(arrin2d,xwadvec,kmx,nx,ny,dsu)
use xmpi_module

IMPLICIT NONE

integer                                         :: i,nx,ny
integer                                         :: j
real*8 , dimension(nx+1,ny+1)                   :: dsu
real*8 , dimension(nx+1,ny+1)                   :: xwadvec,arrin2d,kmx

xwadvec = 0.d0

do j=2,ny
    do i=2,nx   
        if (kmx(i,j)>0) then
           xwadvec(i,j)=(arrin2d(i,j)-arrin2d(i-1,j))/dsu(i-1,j)
        elseif (kmx(i,j)<0) then
           xwadvec(i,j)=(arrin2d(i+1,j)-arrin2d(i,j))/dsu(i,j)
        else
           xwadvec(i,j)=(arrin2d(i+1,j)-arrin2d(i-1,j))/(dsu(i,j)+dsu(i-1,j))
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

subroutine advecwy(arrin2d,ywadvec,kmy,nx,ny,dnv)
use xmpi_module
use xmpi_module
IMPLICIT NONE

integer                                         :: i,nx,ny
integer                                         :: j
real*8 , dimension(nx+1,ny+1)                   :: dnv
real*8 , dimension(nx+1,ny+1)                   :: ywadvec,arrin2d,kmy

ywadvec = 0.d0

do j=2,ny
    do i=2,nx
        if (kmy(i,j)>0) then
           ywadvec(i,j)=(arrin2d(i,j)-arrin2d(i,j-1))/dnv(i,j-1)
        elseif (kmy(i,j)<0) then
           ywadvec(i,j)=(arrin2d(i,j+1)-arrin2d(i,j))/dnv(i,j)
        else
           ywadvec(i,j)=(arrin2d(i,j+1)-arrin2d(i,j-1))/(dnv(i,j)+dnv(i,j-1))
        endif
    end do
end do

if(ny>0) then
   ywadvec(:,1)= ywadvec(:,2)          !Ap
   ywadvec(:,ny+1) = ywadvec(:,ny)     !Ap
endif

#ifdef USEMPI
call xmpi_shift(ywadvec,'m:') ! wwvv fill in yadvec(nx+1,:)
call xmpi_shift(ywadvec,'1:') !      fill in yadvec(1,:)
call xmpi_shift(ywadvec,':n') !      fill in yadvec(:,ny+1)
call xmpi_shift(ywadvec,':1') !      fill in yadvec(:,1)
#endif

end subroutine advecwy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine advecqx(c,arrin2d,xwadvec,nx,ny,dsu)
use xmpi_module

IMPLICIT NONE

integer                                         :: i,nx,ny
integer                                         :: j
real*8 , dimension(nx+1,ny+1)                   :: dsu
real*8 , dimension(nx+1,ny+1)                   :: xwadvec,arrin2d,c

xwadvec = 0.d0

do j=2,ny
    do i=2,nx   
        if (c(i,j)>0) then
           xwadvec(i,j)=c(i,j)*(arrin2d(i,j)-arrin2d(i-1,j))/dsu(i-1,j)
        elseif (c(i,j)<0) then
           xwadvec(i,j)=c(i,j)*(arrin2d(i+1,j)-arrin2d(i,j))/dsu(i,j)
        else
           xwadvec(i,j)=c(i,j)*(arrin2d(i+1,j)-arrin2d(i-1,j))/(dsu(i,j)+dsu(i-1,j))
        endif
    end do
end do

if(ny>0) then
  xwadvec(:,1)= xwadvec(:,2)          !Ap
  xwadvec(:,ny+1) = xwadvec(:,ny)     !Ap
endif

#ifdef USEMPI
call xmpi_shift(xwadvec,'m:') ! fill in xwadvec(nx+1,:)
call xmpi_shift(xwadvec,'1:') ! fill in xwadvec(1,:)
call xmpi_shift(xwadvec,':n') ! fill in xwadvec(:,ny+1)
call xmpi_shift(xwadvec,':1') ! fill in xwadvec(:,1)
#endif

end subroutine advecqx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine advecqy(c,arrin2d,ywadvec,nx,ny,dnv)
use xmpi_module
use xmpi_module
IMPLICIT NONE

integer                                         :: i,nx,ny
integer                                         :: j
real*8 , dimension(nx+1,ny+1)                   :: dnv
real*8 , dimension(nx+1,ny+1)                   :: ywadvec,arrin2d,c

ywadvec = 0.d0

do j=2,ny
    do i=2,nx
        if (c(i,j)>0) then
           ywadvec(i,j)=c(i,j)*(arrin2d(i,j)-arrin2d(i,j-1))/dnv(i,j-1)
        elseif (c(i,j)<0) then
           ywadvec(i,j)=c(i,j)*(arrin2d(i,j+1)-arrin2d(i,j))/dnv(i,j)
        else
           ywadvec(i,j)=c(i,j)*(arrin2d(i,j+1)-arrin2d(i,j-1))/(dnv(i,j)+dnv(i,j-1))
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
  use logging_module

  ! Robert: iteration along L=L0tanh(2pih/L)

  IMPLICIT NONE

  type(spacepars)                     :: s
  type(parameters)                    :: par

  real*8, dimension(s%nx+1,s%ny+1)    :: h,L0,kh
  real*8, dimension(:,:),allocatable,save  :: Ltemp
  !
  ! real*8,dimension(:,:),allocatable,save :: L1
  ! wwvv moved L1 to spaceparams, in the parallel version
  ! of the program L1 will be resized and distributed
  !
  integer                             :: i,j,j1,j2
  real*8                              :: backdis,lback,disfac
  integer                             :: index
  
  if (s%ny==0) then
     j1=1
     j2=1
  else
     j1=2
     j2=s%ny
  endif
  !
  ! In the original code, phi, aphi, bphi are saved
  ! variables, and calculated once. I think this is
  ! better. In the original code, these variables
  ! had to be broadcasted in the parallel version. 
  ! Also I rearranged some formulas, no need anymore
  ! for variables t and n.

  ! cjaap: replaced par%hmin by par%eps

  h = max(s%hh + par%delta*s%H,par%eps)

  L0 = 2*par%px*par%g/(s%sigm**2)

     if (.not. associated(s%L1)) then
        allocate(s%L1(s%nx+1,s%ny+1))
        s%L1=L0
     endif
     if (.not. allocated(Ltemp)) then
        allocate(Ltemp(s%nx+1,s%ny+1))
        Ltemp = L0
     end if

  do j = j1,j2
     do i = 1,s%nx+1
        Ltemp(i,j) = iteratedispersion(L0(i,j),Ltemp(i,j),par%px,h(i,j))
        if (Ltemp(i,j)<0.d0) then   ! this is an error from iteratedispersion
           Ltemp(i,j) = -Ltemp(i,j)
           call writelog('lws','','Warning: no convergence in dispersion relation iteration at t = ', &
                par%t*max(par%morfac*par%morfacopt,1.d0))
        endif
     end do
  end do
  if (par%shoaldelay==1) then
     ! find Lmod looking back over distance par%facsd*L1
     ! presumes sigma direction is shore normal
     s%L1 = 0.d0           ! modified wave length, initially set to L1
     do j = j1,j2
        do i = 2,s%nx+1
           index = i       ! start index
           backdis = 0.d0  ! relative distance backward
           do while (backdis<1.d0)
              ! disfac = s%dsc(index,j)/(par%facsd*s%L1(index,j))
              ! use average wavelength over distance dsc
              disfac = s%dsc(index,j)/(par%facsd*0.5d0*(Ltemp(index,j)+Ltemp(max(index-1,1),j)))
              disfac = min(disfac,1.d0-backdis)
              
              !h(i,j) = h(i,j) + disfac*s%hh(index,j)
              s%L1(i,j) = s%L1(i,j)+disfac*0.5d0*(Ltemp(index,j)+Ltemp(max(index-1,1),j))
              backdis = backdis+disfac
              
              index = max(index-1,1)
           enddo
        enddo
     enddo
     s%L1(1,:) = Ltemp(1,:) 
  else
     s%L1 = Ltemp    
  endif

  ! boundary copies for non superfast 1D
  if (s%ny>0) then 
     s%L1(:,1)=s%L1(:,2)
     s%L1(:,s%ny+1)=s%L1(:,s%ny)
  endif
  s%k  = 2*par%px/s%L1
  s%c  = s%sigm/s%k
  !kh   = s%k*h
  ! Ad:
  kh   = min(s%k*h,10.0d0)
  s%n=0.5d0+kh/sinh(2*kh)
  s%cg=s%c*s%n
  !s%cg = s%c*(0.5d0+kh/sinh(2*kh))

end subroutine dispersion

elemental function iteratedispersion(L0,Lestimate,px,h) result(L)

  implicit none
  ! input
  real*8,intent(in)    :: L0
  real*8,intent(in)    :: Lestimate
  real*8,intent(in)    :: px
  real*8,intent(in)    :: h
  ! output
  real*8               :: L
  ! internal
  real*8               :: L1,L2
  integer              :: iter
  real*8               :: err
  real*8,parameter     :: aphi = 1.d0/(((1.0d0 + sqrt(5.0d0))/2)+1)
  real*8,parameter     :: bphi = ((1.0d0 + sqrt(5.0d0))/2)/(((1.0d0 + sqrt(5.0d0))/2)+1)
  integer,parameter    :: itermax = 150
  real*8,parameter     :: errmax = 0.00001d0
  
    
  err = huge(0.0d0)
  iter = 0
  L1 = Lestimate
  do while (err > errmax .and. iter < itermax)
     iter  = iter+1
     L2    = L0*tanh(2*px*h/L1)
     L1    = (L1*aphi + L2*bphi)          ! Golden ratio
     err   = abs(L2 - L1)
  end do
  
  if (iter<=itermax) then
     L = L1
  else
     ! signal this went wrong
     L = -L1
  endif

end function iteratedispersion

subroutine breakerdelay(par,s)

    use params
    use spaceparams
    use xmpi_module

    implicit none

    type(spacepars),target                          :: s
    type(parameters)                                :: par

    real*8                                          :: Lbr
    real*8, dimension(s%nx+1)                       :: utemp
    integer                                         :: jx,jy,i,j1,nbr,tempxid
    integer, dimension(s%nx+1)                      :: ibr

    include 's.ind'
    include 's.inp'

    ! Superfast 1D
    if (ny>0) then
        j1 = 2
    else
        j1 = 1
    endif

    do jy = j1,max(1,ny)
        usd(1,jy)   = ustr(1,jy)
        
        do jx = 2,nx+1
            nbr     = 0
            Lbr     = sqrt(par%g*hh(jx,jy))*par%Trep
            i       = jx-1
            do while (abs(xz(i,jy)-xz(jx,jy))<=Lbr .and. i>1)
                nbr = nbr+1
                i   = i-1
            end do

            if(nbr.gt.1) then
                do i = 1,nbr+1
                    ibr(i)      = i
                    tempxid     = jx-nbr+i-1
                    utemp(i)    = ustr(tempxid,jy)
                enddo
                
                usd(jx,jy)      = sum(ibr(1:nbr+1)*utemp(1:nbr+1))/sum(ibr(1:nbr+1)) 
            else
                usd(jx,jy)      = ustr(jx,jy)
            end if
        end do
    end do

    ! lateral boundaries
    if (xmpi_istop             ) then
        usd(1,:)    = usd(2,:)
    endif
    if (xmpi_isbot             ) then
        usd(nx+1,:) = usd(nx,:)
    endif
    if (xmpi_isleft  .and. ny>0) then
        usd(:,1)    = usd(:,2)
    endif
    if (xmpi_isright .and. ny>0) then
        usd(:,ny+1) = usd(:,ny)
    endif

    ! wwvv for the parallel version, shift in the columns and rows
#ifdef USEMPI
    call xmpi_shift(usd,'m:')  ! fill in usd(nx+1,:)
    call xmpi_shift(usd,'1:')  ! fill in usd(1,:)
    call xmpi_shift(usd,':n')  ! fill in usd(:,ny+1)
    call xmpi_shift(usd,':1')  ! fill in usd(:,1)
#endif

end subroutine breakerdelay 

end module wave_functions_module
