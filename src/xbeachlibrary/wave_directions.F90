module wave_directions_module
   implicit none
   save
contains
   subroutine wave_directions(s,par)
      use params
      use spaceparams
      use roelvink_module
      use wave_functions_module
      use xmpi_module
      use logging_module
      use paramsconst

      ! wwvv in my testcase, this routine was not called, so it is not
      ! tested. Nevertheless, I put in code for the parallel version.

      IMPLICIT NONE

      type(spacepars), target     :: s
      type(parameters)            :: par

      integer                     :: i,imax,i1
      integer                     :: j
      integer                     :: itheta,iter
      real*8,dimension(:),allocatable,save        :: dist,factor,e01
      integer, dimension(:,:,:),allocatable,save  :: wete3d
      real*8 , dimension(:,:)  ,allocatable,save  :: dhdx,dhdy,dudx,dudy,dvdx,dvdy
      real*8 , dimension(:,:)  ,allocatable,save  :: km,uorb
      real*8 , dimension(:,:)  ,allocatable,save  :: kmx,kmy,sinh2kh ! ,wm
      real*8 , dimension(:,:,:),allocatable,save  :: xadvec,yadvec,thetaadvec,dd
      real*8 , dimension(:,:,:),allocatable,save  :: cgxu,cgyv
      real*8 , dimension(:),allocatable,save      :: Hprev
      real*8                                      :: Herr,dtw
      real*8 , dimension(:)  ,allocatable,save    :: dkmxdx,dkmxdy,dkmydx,dkmydy,cgxm,cgym,arg,fac,xwadvec,ywadvec
      real*8 , dimension(:,:),allocatable,save    :: wcifacu,wcifacv
      logical                                     :: stopiterate

      if (.not. allocated(wete3d)) then
         allocate(e01(1:s%ntheta_s))
         allocate(dist(1:s%ntheta_s))
         allocate(factor(1:s%ntheta_s))
         allocate(wete3d     (s%nx+1,s%ny+1,s%ntheta_s))
         allocate(xadvec    (s%nx+1,s%ny+1,s%ntheta_s))
         allocate(yadvec    (s%nx+1,s%ny+1,s%ntheta_s))
         allocate(thetaadvec(s%nx+1,s%ny+1,s%ntheta_s))
         allocate(dd        (s%nx+1,s%ny+1,s%ntheta_s))

         allocate(cgxu        (s%nx+1,s%ny+1,s%ntheta_s))
         allocate(cgyv        (s%nx+1,s%ny+1,s%ntheta_s))

         allocate(dhdx(s%nx+1,s%ny+1))
         allocate(dhdy(s%nx+1,s%ny+1))
         allocate(dudx(s%nx+1,s%ny+1))
         allocate(dudy(s%nx+1,s%ny+1))
         allocate(dvdx(s%nx+1,s%ny+1))
         allocate(dvdy(s%nx+1,s%ny+1))
         allocate(km  (s%nx+1,s%ny+1))
         allocate(uorb(s%nx+1,s%ny+1))
         allocate(kmx (3,s%ny+1))
         allocate(kmy (3,s%ny+1))
         allocate(xwadvec(s%ny+1))
         allocate(ywadvec(s%ny+1))
         allocate(sinh2kh(s%nx+1,s%ny+1))
         allocate(Hprev(s%ny+1))

         allocate(dkmxdx  (s%ny+1))
         allocate(dkmxdy  (s%ny+1))
         allocate(dkmydx  (s%ny+1))
         allocate(dkmydy  (s%ny+1))
         allocate(cgxm    (s%ny+1))
         allocate(cgym    (s%ny+1))
         allocate(arg     (s%ny+1))
         allocate(fac     (s%ny+1))
         allocate(wcifacu     (s%nx+1,s%ny+1))
         allocate(wcifacv     (s%nx+1,s%ny+1))

      endif

   
      xadvec      = 0.0d0
      yadvec      = 0.0d0
      thetaadvec  = 0.0d0
      dd          = 0.0d0
      dhdx        = 0.0d0
      dhdy        = 0.0d0
      dudx        = 0.0d0
      dudy        = 0.0d0
      dvdx        = 0.0d0
      dvdy        = 0.0d0
      km          = 0.0d0
      kmx         = 0.0d0
      kmy         = 0.0d0
      !wm          = 0.0d0
      xwadvec     = 0.0d0
      ywadvec     = 0.0d0
      sinh2kh     = 0.0d0

      dkmxdx      = 0.0d0
      dkmxdy      = 0.0d0
      dkmydx      = 0.0d0
      dkmydy      = 0.0d0
      cgxm        = 0.0d0
      cgym        = 0.0d0

      arg         = 0.0d0
      fac         = 0.0d0
      
      do i=1,s%ntheta_s
         where(s%wete==1)
            wete3d(:,:,i)=1
         elsewhere
            wete3d(:,:,i)=0
         endwhere
      enddo
          
      ! cjaap: replaced par%hmin by par%eps
      !s%hh = max(s%hh,par%eps)
      ! Robert: how is this mass conservative?
      !where(wete2d==1)
      !   s%hh =max(s%hh,par%hmin) ! note: replace with new limiter?
      !endwhere

      ! Slopes of water depth
      call slope2D(max(s%hh,par%delta*s%H),s%nx,s%ny,s%dsu,s%dnv,dhdx,dhdy,s%wete)
      ! Dano limit slopes used in refraction to avoid unrealistic refraction speeds
      dhdx=sign(1.d0,dhdx)*min(abs(dhdx),0.1d0)
      dhdy=sign(1.d0,dhdy)*min(abs(dhdy),0.1d0)
      call slope2D(s%u*par%wci,s%nx,s%ny,s%dsu,s%dnv,dudx,dudy,s%wete)
      call slope2D(s%v*par%wci,s%nx,s%ny,s%dsu,s%dnv,dvdx,dvdy,s%wete)

      ! wwvv these slope routines are in wave_timestep, and are
      !   MPI-aware
      !
      ! Calculate once sinh(2kh)
      
      where(s%wete==1 .and. 2*s%hh*s%k<=3000.d0)
         sinh2kh=sinh(min(2*s%k*max(s%hh,par%delta*s%H),10.0d0))
      elsewhere
         sinh2kh = 3000.d0
      endwhere

      !   Dano: This is already done in wave_bc
      !   call dispersion(par,s)

      !   if (par%wci==0) then
      !      sigm=2.d0*par%px/par%Trep
      !      sigm = max(sigm,0.010d0)
      !   endif

      ! Calculate once velocities used with and without wave current interaction
      where(s%wete==1)
         wcifacu=s%u*par%wci*min(s%hh/par%hwci,1.d0)
         wcifacv=s%v*par%wci*min(s%hh/par%hwci,1.d0)
      endwhere
      
      DO itheta=1,s%ntheta_s
         where(s%wete==1)
            s%cgx_s(:,:,itheta)= s%cg*s%costh_s(:,:,itheta)+wcifacu
            s%cgy_s(:,:,itheta)= s%cg*s%sinth_s(:,:,itheta)+wcifacv
            s%ctheta_s(:,:,itheta)=  &
               s%sigm/sinh2kh*(dhdx*s%sinth_s(:,:,itheta)-dhdy*s%costh_s(:,:,itheta)) + &
               par%wci*(&
                  s%costh_s(:,:,itheta)*(s%sinth_s(:,:,itheta)*dudx - s%costh_s(:,:,itheta)*dudy) + &
                  s%sinth_s(:,:,itheta)*(s%sinth_s(:,:,itheta)*dvdx - s%costh_s(:,:,itheta)*dvdy))
         endwhere
      END DO
      ! Dano Limit unrealistic refraction speed to 1/2 pi per wave period
      where(wete3d==1)
         s%ctheta_s=sign(1.d0,s%ctheta_s)*min(abs(s%ctheta_s),.5*par%px/par%Trep)
      endwhere
      
      where(s%wete==1)
         km = s%k
      endwhere
            
      forall (i=1:s%nx+1,j=1:s%ny+1,s%wete(i,j)==1)
         s%thetamean(i,j)=(sum(s%ee_s(i,j,:)*s%thet_s(i,j,:))/s%ntheta_s)/(max(sum(s%ee_s(i,j,:)),0.00001d0)/s%ntheta_s)
      endforall

      !dtw=.9*minval(xz(2:nx+1)-xz(1:nx))/maxval(cgx_s)

      s%E(1,:)=sum(s%ee_s(1,:,:),2)*s%dtheta_s
      s%H(1,:)=sqrt(s%E(1,:)*par%irhog8)

      imax=s%nx
      !Dano  This is ok, since we will set mpiboundary to y in stationary mode

      do i=2,imax
         dtw=.5*minval(s%dsu(i:i+1,jmin_ee:jmax_ee))/maxval(s%cgx_s(i-1:i+1,jmin_ee:jmax_ee,:))
         dtw=min(dtw,.5*minval(s%dnv(i,jmin_ee:jmax_ee))/maxval(abs(s%cgy_s(i,jmin_ee:jmax_ee,:))))
         dtw=min(dtw,.5*s%dtheta/maxval(abs(s%ctheta_s(i,jmin_ee:jmax_ee,:))))
         !Dano: need to make sure all processes use the same dtw, min of all processes
#ifdef USEMPI
         call xmpi_allreduce(dtw,MPI_MIN)
#endif
         Herr=1.
         iter=0
         where(s%wete(i,:)==1)
            arg = min(100.0d0,km(i,:)*(s%hh(i,:)+par%delta*s%H(i,:)))
            arg = max(arg,0.0001)
            fac = ( 1.d0 + ((km(i,:)*s%H(i,:)/2.d0)**2))  ! use deep water correction instead of expression above (waves are short near blocking point anyway)
         endwhere
         stopiterate=.false.
         do while (stopiterate .eqv. .false.)
            iter=iter+1
            Hprev=s%H(i,:)
            ! WCI
            if (par%wci==1) then
               ! Dano NEED TO CHECK THIS FOR CURVI
               kmx = km(i-1:i+1,:)*cos(s%thetamean(i-1:i+1,:)-s%alfaz(i-1:i+1,:))
               kmy = km(i-1:i+1,:)*sin(s%thetamean(i-1:i+1,:)-s%alfaz(i-1:i+1,:))
               s%wm(i-1:i+1,:) = s%sigm(i-1:i+1,:)+kmx*wcifacu(i-1:i+1,:)&
               +kmy*wcifacv(i-1:i+1,:)

               where(km(i,:)>0.01d0)
                  s%c(i,:)  = s%sigm(i,:)/km(i,:)
                  s%cg(i,:) = s%c(i,:)*(0.5d0+arg/sinh(2*arg))*sqrt(fac)
                  s%n(i,:)  = 0.5d0+km(i,:)*s%hh(i,:)/sinh(2*max(km(i,:),0.00001d0)*s%hh(i,:))
               elsewhere
                  s%c(i,:)  = 0.01d0
                  s%cg(i,:) = 0.01d0
                  s%n(i,:)  = 1.d0
               endwhere

               cgym = s%cg(i,:)*dsin(s%thetamean(i,:)-s%alfaz(i,:)) + wcifacv(i,:)
               cgxm = s%cg(i,:)*dcos(s%thetamean(i,:)-s%alfaz(i,:)) + wcifacu(i,:)

               dkmxdx       = (kmx(3,:)-kmx(1,:))/(s%dsu(i,:)+s%dsu(i+1,:))
               dkmxdy(2:s%ny) = (kmx(2,3:s%ny+1)-kmx(2,1:s%ny-1))/(s%dnv(i,2:s%ny)+s%dnv(i,3:s%ny+1))
               dkmxdy(1)    = dkmxdy(2)
               dkmxdy(s%ny+1) = dkmxdy(s%ny)
               dkmydx       = (kmy(3,:)-kmy(1,:))/(s%dsu(i,:)+s%dsu(i+1,:))
               dkmydy(2:s%ny) = (kmy(2,3:s%ny+1)-kmy(2,1:s%ny-1))/(s%dnv(i,2:s%ny)+s%dnv(i,3:s%ny+1))
               dkmydy(1)    = dkmydy(2)
               dkmydy(s%ny+1) = dkmydy(s%ny)

               xwadvec  = (s%wm(i,:)-s%wm(i-1,:))/s%dsu(i-1,:)
               kmx(2,:) = kmx(2,:) -dtw*xwadvec -dtw*cgym*(dkmydx-dkmxdy)

               ywadvec(2:s%ny) = (s%wm(i,3:s%ny+1)-s%wm(i,1:s%ny-1))/(s%dnv(i,2:s%ny)+s%dnv(i,3:s%ny+1))
               ywadvec(1)=ywadvec(2)
               ywadvec(s%ny+1)=ywadvec(s%ny)
               kmy(2,:) = kmy(2,:) -dtw*ywadvec + dtw*cgxm*(dkmydx-dkmxdy)

               ! Dano
#ifdef USEMPI
               call xmpi_shift(kmx(1:2,:),SHIFT_Y_R,1,2)
               call xmpi_shift(kmx(1:2,:),SHIFT_Y_L,3,4)
               call xmpi_shift(kmy(1:2,:),SHIFT_Y_R,1,2)
               call xmpi_shift(kmy(1:2,:),SHIFT_Y_L,3,4)
#endif
               km(i,:) = sqrt(kmx(2,:)**2+kmy(2,:)**2)
               km(i,:) = min(km(i,:),25.d0) ! limit to gravity waves

               !  non-linear dispersion
               arg = min(100.0d0,km(i,:)*(s%hh(i,:)+par%delta*s%H(i,:)))
               arg = max(arg,0.0001)
               fac = ( 1.d0 + ((km(i,:)*s%H(i,:)/2.d0)**2))
               s%sigm(i,:) = sqrt( par%g*km(i,:)*tanh(arg)*fac)
               s%sigm(i,:) = max(s%sigm(i,:),0.010d0)
            endif

            !
            ! transform to wave action
            !
            i1=max(i-2,1)
            do itheta=1,s%ntheta_s
               where(s%wete(i1:i+1,:)==1)
                  s%ee_s(i1:i+1,:,itheta) = s%ee_s(i1:i+1,:,itheta)/s%sigm(i1:i+1,:)
               endwhere
            enddo
            !
            ! Upwind Euler timestep propagation
            !
            if  (i>2.and. par%scheme==SCHEME_UPWIND_2) then
               call advecxho(s%ee_s(i-2:i+1,:,:),s%cgx_s(i-2:i+1,:,:),xadvec(i-2:i+1,:,:),    &
               3,s%ny,s%ntheta_s,s%dnu(i-2:i+1,:),s%dsu(i-2:i+1,:),s%dsdnzi(i-2:i+1,:),SCHEME_UPWIND_2, &
               s%wete(i-2:i+1,:))
            else
               call advecxho(s%ee_s(i-1:i+1,:,:),s%cgx_s(i-1:i+1,:,:),xadvec(i-1:i+1,:,:),    &
               2,s%ny,s%ntheta_s,s%dnu(i-1:i+1,:),s%dsu(i-1:i+1,:),s%dsdnzi(i-1:i+1,:),SCHEME_UPWIND_1,&
               s%wete(i-1:i+1,:))
            endif
            call advecyho(s%ee_s(i,:,:),s%cgy_s(i,:,:),yadvec(i,:,:),                                  &
            0,s%ny,s%ntheta_s,s%dsv(i,:),s%dnv(i,:),s%dsdnzi(i,:),SCHEME_UPWIND_1,s%wete(i,:))
            call advecthetaho(s%ee_s(i,:,:),s%ctheta_s(i,:,:),thetaadvec(i,:,:),0,s%ny,s%ntheta_s,s%dtheta,par%scheme,s%wete(i,:))

            where(wete3d(i,:,:)==1)
               s%ee_s(i,:,:)=s%ee_s(i,:,:)-dtw*(xadvec(i,:,:) + yadvec(i,:,:) + thetaadvec(i,:,:))
            elsewhere
               s%ee_s(i,:,:)=0.d0
            endwhere
#ifdef USEMPI
            call xmpi_shift(s%ee_s(i-1:i,:,:),SHIFT_Y_R,1,2)
            call xmpi_shift(s%ee_s(i-1:i,:,:),SHIFT_Y_L,3,4)
#endif
            !
            ! transform back to wave energy
            !
            do itheta=1,s%ntheta_s
               where(s%wete(i1:i+1,:)==1)
                  s%ee_s(i1:i+1,:,itheta) = s%ee_s(i1:i+1,:,itheta)*s%sigm(i1:i+1,:)
               endwhere
            enddo
            where(wete3d(i,:,:)==1)
               s%ee_s(i,:,:)=max(s%ee_s(i,:,:),0.0d0)
            elsewhere
               s%ee_s(i,:,:)=0.d0
            endwhere

            !
            ! Energy integrated over wave directions,Hrms
            !
            forall(j=1:s%ny+1,s%wete(i,j)==1)
               s%E(i,j)=sum(s%ee_s(i,j,:))*s%dtheta
            endforall
            where(s%wete(i,:)==1)
               s%H(i,:)=sqrt(s%E(i,:)*par%irhog8)
            endwhere
            do itheta=1,s%ntheta_s
               where(s%wete(i,:)==1)
                  s%ee_s(i,:,itheta)=s%ee_s(i,:,itheta)/max(1.0d0,(s%H(i,:)/(par%gammax*s%hh(i,:)))**2)
               endwhere
            enddo
            where(s%wete(i,:)==1)
               s%H(i,:)=min(s%H(i,:),par%gammax*s%hh(i,:))
               s%E(i,:)=par%rhog8*s%H(i,:)**2
            endwhere
            if (par%snells==0) then !Dano not for SNellius
               where(s%wete(i,:)==1)
                  s%thetamean(i,:) = (sum(s%ee_s(i,:,:)*s%thet_s(i,:,:),2)/s%ntheta_s)/ &
                                     (max(sum(s%ee_s(i,:,:),2),0.000010d0)/s%ntheta_s)
               endwhere
            endif
            !
            ! Total dissipation

            select case(par%break)
             case(BREAK_ROELVINK1,BREAK_ROELVINK2)
               call roelvink       (par,s,km(i,:),i)
             case(BREAK_BALDOCK)
               call baldock        (par,s,km(i,:),i)
             case(BREAK_JANSSEN)
               call janssen_battjes(par,s,km(i,:),i)
            end select


            ! Dissipation by bed friction
            if(par%fw>0.d0) then
               where(s%wete(i,:)==1 .and. s%hh(i,:)>par%fwcutoff)
                  uorb(i,:)=par%px*s%H(i,:)/par%Trep/sinh(min(max(s%k(i,:),0.01d0)*max(s%hh(i,:),par%delta*s%H(i,:)),10.0d0))
                  s%Df(i,:)=0.6666666d0/par%px*par%rho*par%fw*uorb(i,:)**3
               elsewhere
                  s%Df(i,:) = 0.d0
               end where
            else
               s%Df(i,:) = 0.d0
            endif
            !
            ! Distribution of dissipation over directions and frequencies
            !
            do itheta=1,s%ntheta_s
               where(s%wete(i,:)==1)
                  dd(i,:,itheta)=s%ee_s(i,:,itheta)*(s%D(i,:)+s%Df(i,:))/max(s%E(i,:),0.00001d0)
               elsewhere
                  dd(i,:,itheta)=0.d0
               endwhere
            end do
            !
            ! Euler step dissipation
            !
            do j=jmin_ee,jmax_ee
               do itheta=1,s%ntheta_s
                  if (dtw*dd(i,j,itheta)>s%ee_s(i,j,itheta) .and. s%wete(i,j)==1) then
                     dtw=min(dtw,.5*s%ee_s(i,j,itheta)/dd(i,j,itheta))
                  endif
               enddo
            enddo
            !Dano: need to make sure all processes use the same dtw, min of all processes
#ifdef USEMPI
            call xmpi_allreduce(dtw,MPI_MIN)
#endif
            do j=1,s%ny+1
               do itheta=1,s%ntheta_s
                  if(wete3d(i,j,itheta)==1) then
                     s%ee_s(i,j,itheta)=s%ee_s(i,j,itheta)-dtw*dd(i,j,itheta)
                     s%ee_s(i,j,itheta)=max(s%ee_s(i,j,itheta),0.0d0)
                  elseif(wete3d(i,j,itheta)==0) then
                     s%ee_s(i,j,itheta)=0.0d0
                  end if
               end do
            end do
            ! Lateral boundary condition
            if (xmpi_isleft .and. s%ny>0) then
               do itheta=1,s%ntheta_s
                  if (s%sinth_s(i,1,itheta)>=0.) then
                     s%ee_s(i,1,itheta)=s%ee_s(i,2,itheta)
                  endif
               enddo
               km(:,1)=km(:,2)
               s%sigm(:,1)=s%sigm(:,2)
            endif
            if (xmpi_isright .and. s%ny>0) then
               do itheta=1,s%ntheta_s
                  if (s%sinth_s(i,s%ny+1,itheta)<=0.) then
                     s%ee_s(i,s%ny+1,itheta)=s%ee_s(i,s%ny,itheta)
                  endif
               end do
               km(:,s%ny+1)=km(:,s%ny)
               s%sigm(:,s%ny+1)=s%sigm(:,s%ny)
            endif
            !
            ! Compute mean wave direction
            !
            if (par%snells==0) then
               where(s%wete(i,:)==1)
                  s%thetamean(i,:)=(sum(s%ee_s(i,:,:)*s%thet_s(i,:,:),2)/size(s%ee_s(i,:,:),2)) &
                  /(max(sum(s%ee_s(i,:,:),2),0.000010d0) /size(s%ee_s(i,:,:),2))
               endwhere
            endif
            !
            ! Energy integrated over wave directions,Hrms
            !
            where(s%wete(i,:)==1)
               s%E(i,:)=sum(s%ee_s(i,:,:),2)*s%dtheta
               s%H(i,:)=sqrt(s%E(i,:)*par%irhog8)
            elsewhere
               s%E(i,:)=0.d0
               s%H(i,:)=0.d0
            endwhere
            Herr=maxval(abs(Hprev(jmin_ee:jmax_ee)-s%H(i,jmin_ee:jmax_ee)),mask=s%wete(i,:)==1)

#ifdef USEMPI
            call xmpi_allreduce(Herr,MPI_MAX)
#endif
            ! Stopping criteria
            if (iter<par%maxiter) then
               if (Herr<par%maxerror) then
                  stopiterate=.true.
                  if(xmaster) call writelog('ls','(a,i4,a,i4)','Wave propagation row ',i,', iteration ',iter)
               endif
            else
               stopiterate=.true.
               if(xmaster) call writelog('ls','(a,i4,a,i4,a,f5.4)','Wave propagation row ',i,', iteration ',iter,', error: ',Herr)
            endif
         enddo ! End while loop
      enddo ! End do i=2:s%nx loop

      s%ee_s(s%nx+1,:,:) = s%ee_s(s%nx,:,:)
      s%E(s%nx+1,:)    = s%E(s%nx,:)
      s%H(s%nx+1,:)    = s%H(s%nx,:)
      km(s%nx+1,:)   = km(s%nx,:)
      s%sigm(s%nx+1,:) = s%sigm(s%nx,:)
      s%cg(s%nx+1,:)   = s%cg(s%nx,:)
      s%c(s%nx+1,:)    = s%c(s%nx,:)
      s%thet_s(s%nx+1,:,:) = s%thet_s(s%nx,:,:)
      s%k=km
   end subroutine wave_directions
end module wave_directions_module
