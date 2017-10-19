module wave_functions_module

   use paramsconst
   implicit none
   save

contains

   subroutine update_means_wave_flow(s,par)
      use params
      use spaceparams
      implicit none
      
      type(spacepars), target     :: s
      type(parameters)            :: par
      real*8                      :: factime
      
      
      ! For the "stationary" part:
      if (par%wavemodel == WAVEMODEL_SURFBEAT .and. par%single_dir == 1) then
         ! we need smoothed "stationary" values for the wave directions
         factime = 1.d0/(par%wavint*2)*par%dt
         s%hhws = max(factime*s%hhw + (1-factime)*s%hhws,par%hmin)
         if (par%wci==1) then 
            s%uws = factime*s%u + (1-factime)*s%uws
            s%vws = factime*s%v + (1-factime)*s%vws
         endif
      endif
      
      ! for the instationary part
      if (par%wavemodel == WAVEMODEL_SURFBEAT .and. par%wci==1) then
         ! we need smoothed water depth and velocities for wci in wave instationary
         if (par%single_dir == 1) then
            ! maintain consistency with stationary wave directions model
            factime = 1.d0/(par%wavint*2)*par%dt
         else
            ! maintain consistency with boundary conditions smoothing
            factime = 1.d0/par%cats/par%Trep*par%dt
         endif
         s%hhwcins = max(factime*s%hhw + (1-factime)*s%hhwcins,par%hmin)
         s%uwcins = factime*s%u + (1-factime)*s%uwcins
         s%vwcins = factime*s%v + (1-factime)*s%vwcins
      endif
   
   end subroutine update_means_wave_flow

   subroutine slope2D(h,nx,ny,dsu,dnv,dhdx,dhdy,wete)
      use xmpi_module
      implicit none

      integer                           :: i,j,nx,ny
      real*8, dimension(nx+1,ny+1)      :: h,dhdx,dhdy
      real*8, dimension(nx+1,ny+1)           :: dsu
      real*8, dimension(nx+1,ny+1)           :: dnv
      integer, dimension(nx+1,ny+1),intent(in)  :: wete


      ! wwvv dhdx(2:nx,:) is computed, dhdx(1,:) and dhdx(nx+1,:)
      ! get boundary values, so in the parallel case, we need
      ! to do something about that: get the boundaries from
      ! upper and lower neighbours
      
      ! u-gradients
      if(nx+1>=2) then
         forall(i=2:nx,j=1:ny+1,wete(i,j)==1)
            dhdx(i,j)=(h(i+1,j)-h(i-1,j))/(dsu(i,j)+dsu(i-1,j))
         endforall
         forall(j=1:ny+1,wete(1,j)==1)
            dhdx(1,j)=(h(2,j)-h(1,j))/dsu(1,j)
         endforall
         forall(j=1:ny+1,wete(nx+1,j)==1)
            dhdx(nx+1,j)=(h(nx+1,j)-h(nx,j))/dsu(nx,j)
         endforall
         where(wete==0)
            dhdx=0.d0
         endwhere
      else
         dhdx=0.d0
      endif
!#ifdef USEMPI 
!      call xmpi_shift(dhdx,SHIFT_X_U,3,4)
!      call xmpi_shift(dhdx,SHIFT_X_D,1,2)
!#endif 
      !
      ! v-gradients
      if(ny+1>=2) then
         forall(i=1:nx+1,j=2:ny,wete(i,j)==1)
            dhdy(i,j)=(h(i,j+1)-h(i,j-1))/(dnv(i,j)+dnv(i,j-1))
         endforall
         forall(i=1:nx+1,wete(i,1)==1)
            dhdy(i,1)=(h(i,2)-h(i,1))/dnv(i,1)
         endforall
         forall(i=1:nx+1,wete(i,ny+1)==1)
            dhdy(i,ny+1)=(h(i,ny+1)-h(i,ny))/dnv(i,ny)
         endforall
         where(wete==0)
            dhdy=0.d0
         endwhere
      else
         dhdy=0.d0
      endif
!#ifdef USEMPI
!      call xmpi_shift(dhdy,SHIFT_Y_R,1,2)
!      call xmpi_shift(dhdy,SHIFT_Y_L,3,4)
!#endif       
      
   end subroutine slope2D

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine advecxho(ee,cgx,xadvec,nx,ny,ntheta,dnu,dsu,dsdnzi,scheme,wete,dt,dsz)
      use spaceparams
      use xmpi_module

      implicit none

      integer                                         :: i,j,nx,ny,ntheta
      integer, intent(in)                             :: scheme
      integer, dimension(nx+1,ny+1),intent(in)        :: wete
      real*8,  intent(in)                             :: dt
      integer                                         :: itheta
      real*8 , dimension(nx+1,ny+1)                   :: dnu,dsu,dsz,dsdnzi,fluxx
      real*8 , dimension(nx+1,ny+1,ntheta)            :: xadvec,ee,cgx
      real*8                                          :: cgxu,eupw

      integer                                         :: scheme_now

      xadvec = 0.d0
      fluxx  = 0.d0


      ! split into schemes first, less split loops -> more efficiency
      scheme_now=scheme
      select case(scheme_now)
       case(SCHEME_UPWIND_1)
         do itheta=1,ntheta
            do j=1,ny+1
               do i=1,nx  ! Whole domain
                  if(wete(i,j)==1) then
                  cgxu=.5*(cgx(i+1,j,itheta)+cgx(i,j,itheta))
                  if (cgxu>0) then
                     fluxx(i,j)=ee(i,j,itheta)*cgxu*dnu(i,j)
                  else
                     fluxx(i,j)=ee(i+1,j,itheta)*cgxu*dnu(i,j)
                  endif
                  endif
               enddo
            enddo
            !do j=1,ny+1  !
            do j=jmin_ee,jmax_ee
               do i=2,nx
                  if(wete(i,j)==1) then
                  xadvec(i,j,itheta)=(fluxx(i,j)-fluxx(i-1,j))*dsdnzi(i,j)
                  endif
               enddo
            enddo
         enddo
       case(SCHEME_UPWIND_2,SCHEME_WARMBEAM)
         do itheta=1,ntheta
            do j=1,ny+1
               do i=2,nx-1
                  if(wete(i,j)==1) then
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
                  endif
               enddo
               if (xmpi_istop) then
                  i=1   ! only compute for i==1
                  if(wete(i,j)==1) then
                  cgxu=.5*(cgx(i+1,j,itheta)+cgx(i,j,itheta))
                  if (cgxu>0) then
                     fluxx(i,j)=ee(i,j,itheta)*cgxu*dnu(i,j)
                  else
                     !                     eupw=((dsu(i+1,j)+.5*dsu(i+2,j))*ee(i+1,j,itheta)-.5*dsu(i+2,j)*ee(i+2,j,itheta))/dsu(i+1,j)
                     eupw=((dsu(i+1,j)+.5*dsu(i,j))*ee(i+1,j,itheta)-.5*dsu(i,j)*ee(i+2,j,itheta))/dsu(i+1,j)
                     if (eupw<0.d0) eupw=ee(i+1,j,itheta)
                     fluxx(i,j)=eupw*cgxu*dnu(i,j)
                  endif
               endif
               endif
               if (xmpi_isbot) then
                  i=nx  ! only compute for i==nx0
                  if(wete(i,j)==1) then
                  cgxu=.5*(cgx(i+1,j,itheta)+cgx(i,j,itheta))
                  if (cgxu>0) then
                     !                    eupw=((dsu(i,j)+.5*dsu(i-1,j))*ee(i,j,itheta)-.5*dsu(i-1,j)*ee(i-1,j,itheta))/dsu(i-1,j)
                     eupw=((dsu(i-1,j)+.5*dsu(i,j))*ee(i,j,itheta)-.5*dsu(i,j)*ee(i-1,j,itheta))/dsu(i-1,j)
                     if (eupw<0.d0) eupw=ee(i,j,itheta)
                     fluxx(i,j)=eupw*cgxu*dnu(i,j)
                  else
                     fluxx(i,j)=ee(i+1,j,itheta)*cgxu*dnu(i,j)
                  endif
               endif
               endif
            enddo
            do j=jmin_ee,jmax_ee
               do i=2,nx
                  if(wete(i,j)==1) then
                  xadvec(i,j,itheta)=(fluxx(i,j)-fluxx(i-1,j))*dsdnzi(i,j)
                  endif
               enddo
            enddo
         enddo       
      end select
      if (scheme_now==SCHEME_WARMBEAM) then
         do itheta=1,ntheta
            do j=jmin_ee,jmax_ee
               do i=2,nx
                  if(wete(i,j)==1) then
                  xadvec(i,j,itheta)=   xadvec(i,j,itheta)             &
                                       -((ee(i+1,j,itheta)-ee(i  ,j,itheta))/dsu(i  ,j)   &
                                        -(ee(i  ,j,itheta)-ee(i-1,j,itheta))/dsu(i-1,j))/ &
                                          dsz(i,j)*dt/2*cgx(i,j,itheta)**2
                  endif
               enddo
            enddo
         enddo
      endif

   end subroutine advecxho

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine advecthetaho(ee,ctheta,thetaadvec,nx,ny,ntheta,dtheta,scheme,wete)
      use spaceparams
      use xmpi_module

      implicit none

      integer                                         :: i,j,nx,ny,ntheta
      integer, intent(in)                             :: scheme
      integer, dimension(nx+1,ny+1),intent(in)        :: wete
      integer                                         :: itheta
      real*8 , dimension(ntheta)                      :: fluxtheta
      real*8 , dimension(nx+1,ny+1,ntheta)            :: thetaadvec,ee,ctheta
      real*8                                          :: dtheta,ctheta_between,eupw

      integer                                         :: scheme_now

      thetaadvec = 0.d0
      fluxtheta  = 0.d0

      ! No refraction caan take place if ntheta==1
      if (ntheta>1) then

         ! split into schemes first, less split loops -> more efficiency
         scheme_now=scheme
         select case(scheme_now)
          case(SCHEME_UPWIND_1)
            do j=1,ny+1
               do i=1,nx+1
                  if(wete(i,j)==1) then
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
                  endif
               enddo
            enddo
          case(SCHEME_UPWIND_2,SCHEME_WARMBEAM)
            do j=1,ny+1
               do i=1,nx+1
                  if(wete(i,j)==1) then
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
                  endif
               enddo
            enddo
         end select

      endif !ntheta>1
   end subroutine advecthetaho

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine advecyho(ee,cgy,yadvec,nx,ny,ntheta,dsv,dnv,dsdnzi,scheme,wete,dt,dnz)

      implicit none

      integer                                         :: i,j,nx,ny,ntheta
      integer, intent(in)                             :: scheme
      integer, dimension(nx+1,ny+1),intent(in)        :: wete
      real*8,  intent(in)                             :: dt
      integer                                         :: itheta
      real*8 ,  dimension(nx+1,ny+1)                  :: dsv,dnv,dnz,dsdnzi,fluxy
      real*8 ,  dimension(nx+1,ny+1,ntheta)           :: yadvec,ee,cgy
      real*8                                          :: cgyv,eupw

      integer                                         :: scheme_now

      yadvec = 0.d0
      fluxy  = 0.d0

      ! split into schemes first, less split loops -> more efficiency
      scheme_now=scheme
      select case(scheme_now)
       case(SCHEME_UPWIND_1)
         do itheta=1,ntheta
            do j=1,ny
               do i=1,nx+1  ! Whole domain
                  if(wete(i,j)==1) then
                  cgyv=.5*(cgy(i,j+1,itheta)+cgy(i,j,itheta))
                  if (cgyv>0) then
                     fluxy(i,j)=ee(i,j,itheta)*cgyv*dsv(i,j)
                  else
                     fluxy(i,j)=ee(i,j+1,itheta)*cgyv*dsv(i,j)
                  endif
                  endif
               enddo
            enddo
            do j=2,ny
               do i=1,nx+1
                  if(wete(i,j)==1) then
                  yadvec(i,j,itheta)=(fluxy(i,j)-fluxy(i,j-1))*dsdnzi(i,j)
                  endif
               enddo
            enddo
         enddo
       case(SCHEME_UPWIND_2,SCHEME_WARMBEAM)
         do itheta=1,ntheta
            do j=2,ny-1
               do i=1,nx+1
                  if(wete(i,j)==1) then
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
                  endif
               enddo
            enddo
            j=1   ! only compute for j==1
            do i=1,nx+1
               if(wete(i,j)==1) then
               cgyv=.5*(cgy(i,j+1,itheta)+cgy(i,j,itheta))
               if (cgyv>0) then
                  fluxy(i,j)=ee(i,j,itheta)*cgyv*dsv(i,j)
               else
                  !                   eupw=((dnv(i,j+1)+.5*dnv(i,j+2))*ee(i,j+1,itheta)-.5*dnv(i,j+2)*ee(i,j+2,itheta))/dnv(i,j+1)
                  eupw=((dnv(i,j+1)+.5*dnv(i,j))*ee(i,j+1,itheta)-.5*dnv(i,j)*ee(i,j+2,itheta))/dnv(i,j+1)
                  if (eupw<0.d0) eupw=ee(i,j+1,itheta)
                  fluxy(i,j)=eupw*cgyv*dsv(i,j)
               endif
               endif
            enddo
            j=ny ! only compute for j==ny
            do i=1,nx+1
               if(wete(i,j)==1) then
               cgyv=.5*(cgy(i,j+1,itheta)+cgy(i,j,itheta))
               if (cgyv>0) then
                  !                eupw=((dnv(i,j)+.5*dnv(i,j-1))*ee(i,j,itheta)-.5*dnv(i,j-1)*ee(i,j-1,itheta))/dnv(i,j-1)
                  eupw=((dnv(i,j-1)+.5*dnv(i,j))*ee(i,j,itheta)-.5*dnv(i,j)*ee(i,j-1,itheta))/dnv(i,j-1)
                  if (eupw<0.d0) eupw=ee(i,j,itheta)
                  fluxy(i,j)=eupw*cgyv*dsv(i,j)
               else
                  fluxy(i,j)=ee(i,j+1,itheta)*cgyv*dsv(i,j)
               endif
               endif
            enddo
            do j=2,ny
               do i=2,nx+1
                  if(wete(i,j)==1) then
                  yadvec(i,j,itheta)=(fluxy(i,j)-fluxy(i,j-1))*dsdnzi(i,j)
                  endif
               enddo
            enddo
         enddo
      end select
      if (scheme_now==SCHEME_WARMBEAM) then
         do itheta=1,ntheta
            do j=2,ny
               do i=2,nx+1
                  if(wete(i,j)==1) then
                  yadvec(i,j,itheta) = yadvec(i,j,itheta)                                 &
                                       -((ee(i,j+1,itheta)-ee(i,j  ,itheta))/dnv(i,j  )   &
                                        -(ee(i,j  ,itheta)-ee(i,j-1,itheta))/dnv(i,j-1))/ &
                                          dnz(i,j)*dt/2*cgy(i,j,itheta)**2
                  endif
               enddo
            enddo
         enddo
      endif
   end subroutine advecyho


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine advecwx(arrin2d,xwadvec,kmx,nx,ny,dsu,wete)
      use xmpi_module

      implicit none

      integer                                         :: i,nx,ny
      integer                                         :: j
      real*8 , dimension(nx+1,ny+1)                   :: dsu
      real*8 , dimension(nx+1,ny+1)                   :: xwadvec,arrin2d,kmx
      integer, dimension(nx+1,ny+1),intent(in)        :: wete

      xwadvec = 0.d0

      do j=2,ny
         do i=2,nx
            if(wete(i,j)==1) then
            if (kmx(i,j)>0) then
               xwadvec(i,j)=(arrin2d(i,j)-arrin2d(i-1,j))/dsu(i-1,j)
            elseif (kmx(i,j)<0) then
               xwadvec(i,j)=(arrin2d(i+1,j)-arrin2d(i,j))/dsu(i,j)
            else
               xwadvec(i,j)=(arrin2d(i+1,j)-arrin2d(i-1,j))/(dsu(i,j)+dsu(i-1,j))
            endif
            endif
         end do
      end do

      ! wwvv here we miss the computations of the first and last columns and rows,
      !  in the parallel case we shift these in form neighbours

   end subroutine advecwx


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

   subroutine advecwy(arrin2d,ywadvec,kmy,nx,ny,dnv,wete)
      use xmpi_module
      use xmpi_module
      implicit none

      integer                                         :: i,nx,ny
      integer                                         :: j
      real*8 , dimension(nx+1,ny+1)                   :: dnv
      real*8 , dimension(nx+1,ny+1)                   :: ywadvec,arrin2d,kmy
      integer, dimension(nx+1,ny+1),intent(in)        :: wete

      ywadvec = 0.d0

      do j=2,ny
         do i=2,nx
            if(wete(i,j)==1) then
            if (kmy(i,j)>0) then
               ywadvec(i,j)=(arrin2d(i,j)-arrin2d(i,j-1))/dnv(i,j-1)
            elseif (kmy(i,j)<0) then
               ywadvec(i,j)=(arrin2d(i,j+1)-arrin2d(i,j))/dnv(i,j)
            else
               ywadvec(i,j)=(arrin2d(i,j+1)-arrin2d(i,j-1))/(dnv(i,j)+dnv(i,j-1))
            endif
            endif
         end do
      end do

      if(ny>0) then
         where(wete(:,1)==1)
         ywadvec(:,1)= ywadvec(:,2)          !Ap
         endwhere
         where(wete(:,ny+1)==1)
         ywadvec(:,ny+1) = ywadvec(:,ny)     !Ap
         endwhere
      endif


   end subroutine advecwy

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine advecqx(c,arrin2d,xwadvec,nx,ny,dsu,wete)
      use xmpi_module

      implicit none

      integer                                         :: i,nx,ny
      integer                                         :: j
      real*8 , dimension(nx+1,ny+1)                   :: dsu
      real*8 , dimension(nx+1,ny+1)                   :: xwadvec,arrin2d,c
      integer, dimension(nx+1,ny+1),intent(in)        :: wete

      xwadvec = 0.d0

      do j=2,ny
         do i=2,nx
            if(wete(i,j)==1) then
            if (c(i,j)>0) then
               xwadvec(i,j)=c(i,j)*(arrin2d(i,j)-arrin2d(i-1,j))/dsu(i-1,j)
            elseif (c(i,j)<0) then
               xwadvec(i,j)=c(i,j)*(arrin2d(i+1,j)-arrin2d(i,j))/dsu(i,j)
            else
               xwadvec(i,j)=c(i,j)*(arrin2d(i+1,j)-arrin2d(i-1,j))/(dsu(i,j)+dsu(i-1,j))
            endif
            endif
         end do
      end do

      if(ny>0) then
         where(wete(:,1)==1)
         xwadvec(:,1)= xwadvec(:,2)          !Ap
         endwhere
         where(wete(:,ny+1)==1)
         xwadvec(:,ny+1) = xwadvec(:,ny)     !Ap
         endwhere
      endif


   end subroutine advecqx


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

   subroutine advecqy(c,arrin2d,ywadvec,nx,ny,dnv,wete)
      use xmpi_module
      use xmpi_module
      implicit none

      integer                                         :: i,nx,ny
      integer                                         :: j
      real*8 , dimension(nx+1,ny+1)                   :: dnv
      real*8 , dimension(nx+1,ny+1)                   :: ywadvec,arrin2d,c
      integer, dimension(nx+1,ny+1),intent(in)        :: wete

      ywadvec = 0.d0

      do j=2,ny
         do i=2,nx
            if(wete(i,j)==1) then
            if (c(i,j)>0) then
               ywadvec(i,j)=c(i,j)*(arrin2d(i,j)-arrin2d(i,j-1))/dnv(i,j-1)
            elseif (c(i,j)<0) then
               ywadvec(i,j)=c(i,j)*(arrin2d(i,j+1)-arrin2d(i,j))/dnv(i,j)
            else
               ywadvec(i,j)=c(i,j)*(arrin2d(i,j+1)-arrin2d(i,j-1))/(dnv(i,j)+dnv(i,j-1))
            endif
            endif
         end do
      end do

   end subroutine advecqy

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine dispersion(par,s,h)
      use params
      use spaceparams
      use logging_module

      ! Robert: iteration along L=L0tanh(2pih/L)

      implicit none

      type(spacepars)                      :: s
      type(parameters)                     :: par
      real*8, dimension(1:s%nx+1,1:s%ny+1),intent(in) :: h  ! water depth for dispersion can vary for stationary, single_dir and 
                                                            ! instationary computations, with and without wci

      ! Robert: for some reason, MPI version does not like the "allocatable" version
      !         of this subroutine.
      real*8, dimension(1:s%nx+1,1:s%ny+1)  :: L0,kh,Ltemp
      integer                               :: i,j,j1,j2
      real*8                                :: backdis,disfac
      integer                               :: index

      if (s%ny==0) then
         j1=1
         j2=1
      else
         j1=1
         j2=s%ny+1
      endif
      !
      !
      where(s%wete==1)
         L0 = par%g*par%Trep**2/(2*par%px)
      elsewhere
         L0=par%eps
      endwhere

      if (.not. associated(s%L1)) then
         allocate(s%L1(s%nx+1,s%ny+1))
         s%L1=L0
      endif
      Ltemp = L0

      do j = j1,j2
         do i = 1,s%nx+1
            if(s%wete(i,j)==1) then
               Ltemp(i,j) = iteratedispersion(L0(i,j),Ltemp(i,j),par%px,h(i,j))
               if (Ltemp(i,j)<0.d0) then   ! this is an error from iteratedispersion
                  Ltemp(i,j) = -Ltemp(i,j)
                  call writelog('lws','','Warning: no convergence in dispersion relation iteration at t = ', &
                                       par%t*max(par%morfac*par%morfacopt,1.d0))
               endif
            endif
         end do
      end do
      if (par%shoaldelay==1) then
         ! find Lmod looking back over distance par%facsd*L1
         ! presumes sigma direction is shore normal
         s%L1 = 0.d0           ! modified wave length, initially set to L1
         do j = j1,j2
            do i = 2,s%nx+1
               if(s%wete(i,j)==1) then
                  index = i       ! start index
                  backdis = 0.d0  ! relative distance backward
                  do while (backdis<1.d0)
                     ! use average wavelength over distance dsc
                     disfac = s%dsc(index,j)/(par%facsd*0.5d0*(Ltemp(index,j)+Ltemp(max(index-1,1),j)))
                     disfac = min(disfac,1.d0-backdis)


                     s%L1(i,j) = s%L1(i,j)+disfac*0.5d0*(Ltemp(index,j)+Ltemp(max(index-1,1),j))
                     backdis = backdis+disfac

                     index = max(index-1,1)
                  enddo
               endif
            enddo
         enddo
         where(s%wete(1,:)==1)
            s%L1(1,:) = Ltemp(1,:)
         endwhere
      else
         where(s%wete==1)
            s%L1 = Ltemp
         endwhere
      endif

      ! boundary copies for non superfast 1D
      if (s%ny>0) then
         if (xmpi_isleft) then 
            where(s%wete(:,1)==1)
               s%L1(:,1)=s%L1(:,2)
            endwhere
         endif
         if (xmpi_isright) then
            where(s%wete(:,s%ny+1)==1)
               s%L1(:,s%ny+1)=s%L1(:,s%ny)
            endwhere
         endif
      endif
      where(s%wete==1)
         s%k  = 2*par%px/s%L1
         s%c  = s%sigm/s%k
         kh   = min(s%k*h,10.0d0)
         s%n=0.5d0+kh/sinh(2*kh)
         s%cg=s%c*s%n
      elsewhere
         s%k = 25.d0
         s%c = sqrt(par%g*par%eps)
         s%n = 1.d0
         s%cg= sqrt(par%g*par%eps)
      endwhere


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

      !include 's.ind'
      !include 's.inp'

      ! Superfast 1D
      if (s%ny>0) then
         j1 = 2
      else
         j1 = 1
      endif

      do jy = j1,max(1,s%ny)
         s%usd(1,jy)   = s%ustr(1,jy)

         do jx = 2,s%nx+1
            if(s%wete(jx,jy)==1) then
            nbr     = 0
            Lbr     = sqrt(par%g*s%hhw(jx,jy))*par%Trep*par%breakerdelay
            i       = jx-1
            do while (abs(s%xz(i,jy)-s%xz(jx,jy))<=Lbr .and. i>1)
               nbr = nbr+1
               i   = i-1
            end do

            if(nbr.gt.1) then
               do i = 1,nbr+1
                  ibr(i)      = i
                  tempxid     = jx-nbr+i-1
                  utemp(i)    = s%ustr(tempxid,jy)
               enddo

               s%usd(jx,jy)      = sum(ibr(1:nbr+1)*utemp(1:nbr+1))/sum(ibr(1:nbr+1))
            else
               s%usd(jx,jy)      = s%ustr(jx,jy)
            end if
            endif
         end do
      end do

      ! lateral boundaries
      if (xmpi_istop             ) then
         where(s%wete(1,:)==1)
         s%usd(1,:)    = s%usd(2,:)
         endwhere
      endif
      if (xmpi_isbot             ) then
         where(s%wete(s%nx+1,:)==1)
         s%usd(s%nx+1,:) = s%usd(s%nx,:)
         endwhere
      endif
      if (xmpi_isleft  .and. s%ny>0) then
         where(s%wete(:,1)==1)
         s%usd(:,1)    = s%usd(:,2)
         endwhere
      endif
      if (xmpi_isright .and. s%ny>0) then
         where(s%wete(:,s%ny+1)==1)
         s%usd(:,s%ny+1) = s%usd(:,s%ny)
         endwhere
      endif

      ! wwvv for the parallel version, shift in the columns and rows

   end subroutine breakerdelay
   
   subroutine wave_dispersion(s,par,useAverageDepthSwitch)
   
      use params
      use spaceparams
      
      implicit none
      
      type(spacepars), target     :: s
      type(parameters)            :: par
      integer,intent(in),optional :: useAverageDepthSwitch
   
      real*8,dimension(:,:),allocatable,save  :: km,kmx,kmy 
      real*8,dimension(:,:),allocatable,save  :: hhlocal,ulocal,vlocal,relangle 
      real*8,dimension(:,:),allocatable,save  :: arg,fac
      real*8,dimension(:,:),allocatable,save  :: cgym,cgxm
      real*8,dimension(:,:),allocatable,save  :: dkmxdx,dkmxdy,dkmydx,dkmydy
      real*8,dimension(:,:),allocatable,save  :: xwadvec,ywadvec
      real*8,dimension(:),allocatable,save    :: L0,L1
      real*8,save                             :: Trepold
      integer                                 :: itheta,j
      integer                                 :: luseAverageDepthSwitch
   
      if(.not.allocated(km)) then
         allocate(km(s%nx+1,s%ny+1))
         allocate(kmx(s%nx+1,s%ny+1))
         allocate(kmy(s%nx+1,s%ny+1))
         allocate(arg(s%nx+1,s%ny+1))
         allocate(fac(s%nx+1,s%ny+1))
         allocate(cgym(s%nx+1,s%ny+1))
         allocate(cgxm(s%nx+1,s%ny+1))
         allocate(dkmxdx(s%nx+1,s%ny+1))
         allocate(dkmxdy(s%nx+1,s%ny+1))
         allocate(dkmydx(s%nx+1,s%ny+1))
         allocate(dkmydy(s%nx+1,s%ny+1))
         allocate(xwadvec(s%nx+1,s%ny+1))
         allocate(ywadvec(s%nx+1,s%ny+1))
         allocate(hhlocal(s%nx+1,s%ny+1))
         if (par%wci==1) then
            allocate(ulocal(s%nx+1,s%ny+1))
            allocate(vlocal(s%nx+1,s%ny+1))
            allocate(relangle(s%nx+1,s%ny+1))
            allocate(L0(s%ny+1))
            allocate(L1(s%ny+1))
            L0 = 0.d0
            L1 = -huge(0.d0)
         endif
         Trepold = 0.d0
         call dispersion(par,s,s%hhw) ! at initialisation, water depth is always s%hhw
         km=s%k
      endif
   
      ! water depth and velocities to use depends on wave mode and presence of wci
      
      ! default to use instantaneous water depth and velocities
      if (present(useAverageDepthSwitch)) then
         luseAverageDepthSwitch = useAverageDepthSwitch
      else
         luseAverageDepthSwitch = 0
      endif
      
      if (luseAverageDepthSwitch==0) then  ! default: use instantaneous depth and velocities
         hhlocal = s%hhw
         if (par%wci==1) then
            ulocal = s%u
            vlocal = s%v
         endif
      elseif (luseAverageDepthSwitch==1) then ! used averaged depths for single_dir computation
         hhlocal = s%hhws
         if (par%wci==1) then
            ulocal = s%uws
            vlocal = s%vws
         endif
      elseif (luseAverageDepthSwitch==2) then ! used averaged depths for instationary computation with wci
         hhlocal = s%hhwcins
         ulocal = s%uwcins
         vlocal = s%vwcins
      else
         ! this should never occur
      endif
         
       
      if(par%wci==1) then
         ! Dano NEED TO CHECK THIS FOR CURVI
         !if (xmpi_istop) then
         !   ! requires boundary condition at the offshore boundary, where we assume zero flow
         !   L0 = par%g*par%Trep**2/(2*par%px)
         !   if (any(L1<=0.d0)) then
         !      L1 = L0
         !   endif
         !   do j=1,s%ny+1
         !      L1(j) = iteratedispersion(L0(j),L1(j),par%px,s%hhw(1,j))
         !   enddo
         !   km(1,:)  = 2*par%px/max(L1,0.00001d0)
         !endif
         arg     = min(100.0d0,km*hhlocal)
         fac = ( 1.d0 + ((km*s%H/2.d0)**2))  ! use deep water correction
         s%sigm(1,:) = sqrt( par%g*km(1,:)*tanh(arg(1,:))) ! *( 1.d0+ ((km(1,:)*s%H(1,:)/2.d0)**2)))
         !  calculate change in intrinsic frequency
         relangle = s%thetamean-s%alfaz
         kmx = km*dcos(relangle)
         kmy = km*dsin(relangle)
         s%wm = s%sigm+kmx*ulocal*min(&
                                  min(hhlocal/par%hwci,1.d0), &
                                  min(1.d0,(1.d0-hhlocal/par%hwcimax)) &
                                  )+ &
                       kmy*vlocal*min( &
                                  min(hhlocal/par%hwci,1.d0), &
                                  min(1.d0,(1.d0-hhlocal/par%hwcimax)) &
                                  )
         
         where(km>0.01d0)
            s%c  = s%sigm/km
            !          cg = c*(0.5d0+arg/sinh(2.0d0*arg))    ! Linear
            s%cg = s%c*(0.5d0+arg/sinh(2*arg))*sqrt(fac)  ! &  to include more
            !		 	          + km*(H/2)**2*sqrt(max(par%g*km*tanh(arg),0.001d0))/sqrt(max(fac,0.001d0)) ! include wave steepness
            s%n=0.5d0+km*hhlocal/sinh(2*max(km,0.00001d0)*hhlocal)
         elsewhere
            s%c  = 0.01d0
            s%cg = 0.01d0
            s%n  = 1.d0
         endwhere
         
         cgym = s%cg*dsin(relangle)+vlocal*min(min(hhlocal/par%hwci,1.d0),min(1.d0,(1.d0-hhlocal/par%hwcimax)))
         cgxm = s%cg*dcos(relangle)+ulocal*min(min(hhlocal/par%hwci,1.d0),min(1.d0,(1.d0-hhlocal/par%hwcimax)))

         call slope2D(kmx,s%nx,s%ny,s%dsu,s%dnv,dkmxdx,dkmxdy,s%wete)
         call slope2D(kmy,s%nx,s%ny,s%dsu,s%dnv,dkmydx,dkmydy,s%wete)
         call advecwx(s%wm,xwadvec,kmx,s%nx,s%ny,s%dsu,s%wete)   ! cjaap: s%xz or s%xu?         kmx = kmx -par%dt*xwadvec  -par%dt*cgym*(dkmydx-dkmxdy)
         kmx = kmx -par%dt*xwadvec  -par%dt*cgym*(dkmydx-dkmxdy)
         if (s%ny>0) then
            if (xmpi_isright) kmx(:,s%ny+1) = kmx(:,s%ny)  ! lateral bc
            if (xmpi_isleft)  kmx(:,1) = kmx(:,2)  ! lateral bc
         endif

         call advecwy(s%wm,ywadvec,kmy,s%nx,s%ny,s%dnv,s%wete)   ! cjaap: s%yz or s%yv?
         kmy = kmy-par%dt*ywadvec  + par%dt*cgxm*(dkmydx-dkmxdy)
         if (s%ny>0) then
            if (xmpi_isright) kmy(:,s%ny+1) = kmy(:,s%ny)   ! lateral bc
            if (xmpi_isleft)  kmy(:,1) = kmy(:,2)   ! lateral bc
         endif

#ifdef USEMPI
         call xmpi_shift_ee(kmx)
         call xmpi_shift_ee(kmy)
#endif

         ! update km
         km = sqrt(kmx**2+kmy**2)
         km = min(km,25.d0) ! limit to gravity waves
         ! non-linear dispersion
         arg = min(100.0d0,km*hhlocal)
         arg = max(arg,0.0001)
         !       fac = ( 1.d0 + ((km*H/2.d0)**2)*( (8.d0+(cosh(min(4.d0*arg,10.0d0)))**1.d0-2.d0*(tanh(arg))**2.d0 ) /(8.d0*(sinh(arg))**4.d0) ) )
         fac = ( 1.d0 + ((km*s%H/2.d0)**2))  ! use deep water correction instead of expression above (waves are short near blocking point anyway)
         !       fac = 1.d0    ! Linear
         s%sigm = sqrt( par%g*km*tanh(arg)*fac)
         s%sigm = max(s%sigm,0.010d0)
         !  update intrinsic frequency
         do itheta=1,s%ntheta
            s%sigt(:,:,itheta) = s%sigm
         enddo
         !
         !  update k
         s%k = km
      else
         ! check if we need to recompute sigm and sigt
         if (abs(par%Trep-Trepold)/par%Trep > 0.0001d0) then
            Trepold = par%Trep
            do itheta=1,s%ntheta
               s%sigt(:,:,itheta) = 2*par%px/par%Trep
            end do
            s%sigm = 2*par%px/par%Trep
         endif
         call dispersion(par,s,hhlocal)
      endif ! end wave current interaction
      
   
   end subroutine wave_dispersion
   
   subroutine compute_wave_direction_velocities(s,par,flavour,dhdx,dhdy,dudx,dudy,dvdx,dvdy,sinh2kh)
   
      use params
      use spaceparams
      
      implicit none
      
      type(spacepars), target     :: s
      type(parameters)            :: par
      integer,intent(in)          :: flavour ! 0 for stationary, 1 for instationary, 2 for directions part of single_dir
      real*8,dimension(s%nx+1,s%ny+1),intent(in) :: dhdx,dhdy,dudx,dudy,dvdx,dvdy,sinh2kh
      
      real*8,dimension(:,:,:),allocatable,save :: cgx,cgy,cx,cy,ctheta
      real*8,dimension(:,:),allocatable,save   :: uwci,vwci
      
      integer                     :: nthetalocal,nthetamax
      integer                     :: itheta,j,i
      real*8                      :: cs,sn
      
      if(.not. allocated(cgx)) then
         if (par%single_dir==1) then
            nthetamax = max(s%ntheta,s%ntheta_s)
         else
            nthetamax = s%ntheta
         endif
         allocate(cgx(s%nx+1,s%ny+1,nthetamax))
         allocate(cgy(s%nx+1,s%ny+1,nthetamax))
         allocate(cx(s%nx+1,s%ny+1,nthetamax))
         allocate(cy(s%nx+1,s%ny+1,nthetamax))
         allocate(ctheta(s%nx+1,s%ny+1,nthetamax))
         if (par%wci==1) then
            allocate(uwci(s%nx+1,s%ny+1))
            allocate(vwci(s%nx+1,s%ny+1))
         endif
      endif
      
      ! set local variables according to flavour
      select case (flavour)
      case (0)
         nthetalocal = s%ntheta
         if (par%wci==1) then
            uwci = s%u
            vwci = s%v
         endif
      case (1)
         nthetalocal = s%ntheta
         if (par%wci==1) then
            uwci = s%uwcins
            vwci = s%vwcins
         endif
      case (2)
         nthetalocal = s%ntheta_s
         if (par%wci==1) then
            uwci = s%uws
            vwci = s%vws
         endif
      end select
      !   
      ! split wave velocities in wave grid directions theta
      do itheta=1,nthetalocal
         do j=1,s%ny+1
            do i=1,s%nx+1
               if (s%wete(i,j) == 1) then
                  ! wave grid directions theta with respect to spatial grid s,n
                  if (flavour == 2) then
                     cs=s%costh_s(i,j,itheta)
                     sn=s%sinth_s(i,j,itheta)
                  else
                     cs=s%costh(i,j,itheta)
                     sn=s%sinth(i,j,itheta)
                  endif

                  ! split wave velocities over theta bins
                  ! note: cg is already correct for each flavour, as long as wave_dispersion is called before this subroutine
                  if (par%wci==1) then
                     cgx(i,j,itheta)= s%cg(i,j)*cs+uwci(i,j)
                     cgy(i,j,itheta)= s%cg(i,j)*sn+vwci(i,j)
                     cx(i,j,itheta) =  s%c(i,j)*cs+uwci(i,j)
                     cy(i,j,itheta) =  s%c(i,j)*sn+vwci(i,j)
                  else
                     cgx(i,j,itheta)= s%cg(i,j)*cs
                     cgy(i,j,itheta)= s%cg(i,j)*sn
                     cx(i,j,itheta) =  s%c(i,j)*cs
                     cy(i,j,itheta) =  s%c(i,j)*sn
                  endif

                  ! compute refraction velocity
                  ! note: sigm is already correct for each flavour, as long as wave_dispersion is called before this subroutine
                  if (par%wci==1) then
                     ctheta(i,j,itheta)= s%sigm(i,j)/sinh2kh(i,j)*(dhdx(i,j)*sn-dhdy(i,j)*cs) +  &
                                                               (cs*(sn*dudx(i,j)-cs*dudy(i,j)) + sn*(sn*dvdx(i,j)-cs*dvdy(i,j)))
                  else
                     ctheta(i,j,itheta)= s%sigm(i,j)/sinh2kh(i,j)*(dhdx(i,j)*sn-dhdy(i,j)*cs)
                  endif
               else
                  cgx(i,j,itheta) = 0.d0
                  cgy(i,j,itheta) = 0.d0
                  cx(i,j,itheta) = 0.d0
                  cy(i,j,itheta) = 0.d0
                  ctheta(i,j,itheta) = 0.d0
               endif
            enddo
         enddo
      enddo
      
      ! Dano Limit unrealistic refraction speed to 1/2 pi per wave period
      ctheta=sign(1.d0,ctheta)*min(abs(ctheta),.5*par%px/par%Trep)
      
      ! return to the right parts in s
      if (flavour == 2) then
         s%cgx_s(:,:,1:s%ntheta_s) = cgx(:,:,1:s%ntheta_s)
         s%cgy_s(:,:,1:s%ntheta_s) = cgy(:,:,1:s%ntheta_s)
         s%ctheta_s(:,:,1:s%ntheta_s) = ctheta(:,:,1:s%ntheta_s)
         ! no cx and cy needed for roller energy
      else
         s%cgx(:,:,1:s%ntheta) = cgx(:,:,1:s%ntheta)
         s%cgy(:,:,1:s%ntheta) = cgy(:,:,1:s%ntheta)
         s%ctheta(:,:,1:s%ntheta) = ctheta(:,:,1:s%ntheta)
         s%cx(:,:,1:s%ntheta) = cx(:,:,1:s%ntheta)
         s%cy(:,:,1:s%ntheta) = cy(:,:,1:s%ntheta)
      endif

   end subroutine compute_wave_direction_velocities


end module wave_functions_module
