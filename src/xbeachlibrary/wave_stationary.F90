module wave_stationary_module
   implicit none
   save
contains
   subroutine wave_stationary(s,par)
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
      real*8 , dimension(:,:)  ,allocatable,save  :: dhdx,dhdy,dudx,dudy,dvdx,dvdy,ustw
      real*8 , dimension(:,:)  ,allocatable,save  :: sinh2kh ! ,wm
      real*8 , dimension(:,:,:),allocatable,save  :: xadvec,yadvec,thetaadvec,dd,drr,dder
      real*8 , dimension(:,:,:),allocatable,save  :: xradvec,yradvec,thetaradvec
      real*8 , dimension(:),allocatable,save      :: Hprev
      real*8                                      :: Herr,dtw
      real*8 , dimension(:,:),allocatable,save    :: uorb
      logical                                     :: stopiterate

      !include 's.ind'
      !include 's.inp'

      if (.not. allocated(xadvec)) then
         allocate(xadvec    (s%nx+1,s%ny+1,s%ntheta))
         allocate(yadvec    (s%nx+1,s%ny+1,s%ntheta))
         allocate(thetaadvec(s%nx+1,s%ny+1,s%ntheta))
         allocate(xradvec    (s%nx+1,s%ny+1,s%ntheta))
         allocate(yradvec    (s%nx+1,s%ny+1,s%ntheta))
         allocate(thetaradvec(s%nx+1,s%ny+1,s%ntheta))
         allocate(dd        (s%nx+1,s%ny+1,s%ntheta))
         allocate(drr       (s%nx+1,s%ny+1,s%ntheta))
         allocate(dder      (s%nx+1,s%ny+1,s%ntheta))

         allocate(dhdx(s%nx+1,s%ny+1))
         allocate(dhdy(s%nx+1,s%ny+1))
         allocate(dudx(s%nx+1,s%ny+1))
         allocate(dudy(s%nx+1,s%ny+1))
         allocate(dvdx(s%nx+1,s%ny+1))
         allocate(dvdy(s%nx+1,s%ny+1))
         allocate(ustw(s%nx+1,s%ny+1))
         allocate(sinh2kh(s%nx+1,s%ny+1))
         allocate(Hprev(s%ny+1))

         allocate(uorb        (s%nx+1,s%ny+1))

      endif

      xadvec      = 0.0d0
      yadvec      = 0.0d0
      thetaadvec  = 0.0d0
      xradvec     = 0.0d0
      yradvec     = 0.0d0
      thetaradvec = 0.0d0
      dd          = 0.0d0
      drr         = 0.0d0
      dder        = 0.0d0
      dhdx        = 0.0d0
      dhdy        = 0.0d0
      dudx        = 0.0d0
      dudy        = 0.0d0
      dvdx        = 0.0d0
      dvdy        = 0.0d0
      ustw        = 0.0d0
      sinh2kh     = 0.0d0

      uorb        = 0.0d0

      ! Slopes of water depth
      call slope2D(s%hhw,s%nx,s%ny,s%dsu,s%dnv,dhdx,dhdy,s%wete)
      ! Dano limit slopes used in refraction to avoid unrealistic refraction speeds
      dhdx=sign(1.d0,dhdx)*min(abs(dhdx),0.1d0)
      dhdy=sign(1.d0,dhdy)*min(abs(dhdy),0.1d0)
      if (par%wci==1) then
         call slope2D(s%u,s%nx,s%ny,s%dsu,s%dnv,dudx,dudy,s%wete)
         call slope2D(s%v,s%nx,s%ny,s%dsu,s%dnv,dvdx,dvdy,s%wete)
      else
         dudx = 0.d0
         dudy = 0.d0
         dvdx = 0.d0
         dvdy = 0.d0
      endif
      !
      ! wwvv these slope routines are in wave_timestep, and are
      !   MPI-aware
      !
      ! Calculate once sinh(2kh)
      where(2*s%hhw*s%k<=3000.d0)
         sinh2kh=sinh(min(2*s%k*s%hhw,10.0d0))
      elsewhere
         sinh2kh = 3000.d0
      endwhere
      !
      !
      if (par%snells==0) then
         s%thetamean=(sum(s%ee*s%thet,3)/s%ntheta)/(max(sum(s%ee,3),0.00001d0)/s%ntheta)
      else !Dano: Snellius
         s%thetamean=asin(sin(s%theta0-s%alfaz(1,1))*s%c/s%c(1,1))+s%alfaz(1,1)
         s%costh(:,:,1)=cos(s%thetamean-s%alfaz)
         s%sinth(:,:,1)=sin(s%thetamean-s%alfaz)
      endif
      !
      ! split wave velocities in wave grid directions theta
      call compute_wave_direction_velocities(s,par,0,dhdx,dhdy,dudx,dudy,dvdx,dvdy,sinh2kh)
      !
      !
      s%E(1,:)=max(sum(s%ee(1,:,:),2)*s%dtheta,0.0d0)
      s%R(1,:)=max(sum(s%rr(1,:,:),2)*s%dtheta,0.0d0)
      s%H(1,:)=sqrt(s%E(1,:)/par%rhog8)

      imax=s%nx
      !Dano  This is ok, since we will set mpiboundary to y in stationary mode
      !
      ! write to screen that waves are updated
      if(xmaster) call writelog('ls','(a,f0.2,a)','Computing wave transformation at t = ',par%t,' s')
      do i=2,imax
         dtw=.5*minval(s%dsu(i:i+1,jmin_ee:jmax_ee))/max(maxval(s%cgx(i-1:i+1,jmin_ee:jmax_ee,:)),1d-10)
         dtw=min(dtw,.5*minval(s%dnv(i,jmin_ee:jmax_ee))/max(maxval(abs(s%cgy(i,jmin_ee:jmax_ee,:))),1d-10))
         dtw=min(dtw,.5*s%dtheta/max(1.0d-6,maxval(abs(s%ctheta(i,jmin_ee:jmax_ee,:)))))
         !Dano: need to make sure all processes use the same dtw, min of all processes
#ifdef USEMPI
         call xmpi_allreduce(dtw,MPI_MIN)
#endif
         Herr=1.
         iter=0
         stopiterate=.false.
         do while (stopiterate .eqv. .false.)
            iter=iter+1
            Hprev=s%H(i,:)
            !
            ! transform to wave action
            !
            i1=max(i-2,1)
            s%ee(i1:i+1,:,:) = s%ee(i1:i+1,:,:)/s%sigt(i1:i+1,:,:)
            !
            ! Upwind Euler timestep propagation
            !
            if  (i>2 .and. (par%scheme==SCHEME_UPWIND_2 .or. par%scheme==SCHEME_WARMBEAM)) then
               call advecxho(s%ee(i-2:i+1,:,:),s%cgx(i-2:i+1,:,:),xadvec(i-2:i+1,:,:),    &
               3,s%ny,s%ntheta,s%dnu(i-2:i+1,:),s%dsu(i-2:i+1,:),s%dsdnzi(i-2:i+1,:),par%scheme, &
               s%wete(i-2:i+1,:),par%dt,s%dsz)
            else
               call advecxho(s%ee(i-1:i+1,:,:),s%cgx(i-1:i+1,:,:),xadvec(i-1:i+1,:,:),    &
               2,s%ny,s%ntheta,s%dnu(i-1:i+1,:),s%dsu(i-1:i+1,:),s%dsdnzi(i-1:i+1,:),SCHEME_UPWIND_1, &
               s%wete(i-1:i+1,:),par%dt,s%dsz)
            endif
            call advecyho(s%ee(i,:,:),s%cgy(i,:,:),yadvec(i,:,:),                                  &
            0,s%ny,s%ntheta,s%dsv(i,:),s%dnv(i,:),s%dsdnzi(i,:),SCHEME_UPWIND_1, &
            s%wete(i,:),par%dt,s%dnz)
            call advecthetaho(s%ee(i,:,:),s%ctheta(i,:,:),thetaadvec(i,:,:),0,s%ny,s%ntheta,s%dtheta,par%scheme,s%wete(i,:))

            s%ee(i,:,:)=s%ee(i,:,:)-dtw*(xadvec(i,:,:) + yadvec(i,:,:) &
            + thetaadvec(i,:,:))
#ifdef USEMPI
            call xmpi_shift(s%ee(i-1:i,:,:),SHIFT_Y_R,1,2)
            call xmpi_shift(s%ee(i-1:i,:,:),SHIFT_Y_L,3,4)
#endif
            !
            ! transform back to wave energy
            !
            s%ee(i1:i+1,:,:) = s%ee(i1:i+1,:,:)*s%sigt(i1:i+1,:,:)
            s%ee(i,:,:)=max(s%ee(i,:,:),0.0d0)


            !
            ! Energy integrated over wave directions,Hrms
            !
            s%E(i,:)=sum(s%ee(i,:,:),2)*s%dtheta
            s%H(i,:)=sqrt(s%E(i,:)/par%rhog8)
            do itheta=1,s%ntheta
               s%ee(i,:,itheta)=s%ee(i,:,itheta)/max(1.0d0,(s%H(i,:)/(par%gammax*s%hh(i,:)))**2)
            enddo
            s%H(i,:)=min(s%H(i,:),par%gammax*s%hh(i,:))
            s%E(i,:)=par%rhog8*s%H(i,:)**2

            if (par%snells==0) then !Dano not for SNellius
               s%thetamean(i,:) = (sum(s%ee(i,:,:)*s%thet(i,:,:),2)/s%ntheta)/(max(sum(s%ee(i,:,:),2),0.000010d0)/s%ntheta)
            endif
            !
            ! Total dissipation

            select case(par%break)
             case(BREAK_ROELVINK1,BREAK_ROELVINK2)
               call roelvink       (par,s,i)
             case(BREAK_BALDOCK)
               call baldock        (par,s,i)
             case(BREAK_JANSSEN)
               call janssen_battjes(par,s,i)
            end select

            ! Dissipation by bed friction
            uorb(i,:)=par%px*s%H(i,:)/par%Trep/sinh(min(max(s%k(i,:),0.01d0)*s%hhw(i,:),10.0d0))
            s%Df(i,:)=0.28d0*par%rho*s%fw(i,:)*uorb(i,:)**3
            where (s%hh>par%fwcutoff)
               s%Df = 0.d0
            end where
            !
            ! Distribution of dissipation over directions and frequencies
            !
            do itheta=1,s%ntheta
               dder(i,:,itheta)=s%ee(i,:,itheta)*s%D(i,:)/max(s%E(i,:),0.00001d0)
               ! Then all short wave energy dissipation, including bed friction and vegetation
               dd(i,:,itheta)=dder(i,:,itheta) + s%ee(i,:,itheta)*(s%Df(i,:)+s%Dveg(i,:))/max(s%E(i,:),0.00001d0)
            end do
            !
            ! Euler step dissipation
            !
            ! calculate roller energy balance
            !
            call advecxho(s%rr(i-1:i+1,:,:),s%cx(i-1:i+1,:,:),xradvec(i-1:i+1,:,:),   &
            2,s%ny,s%ntheta,s%dnu(i-1:i+1,:),s%dsu(i-1:i+1,:),s%dsdnzi(i-1:i+1,:),SCHEME_UPWIND_1, &
            s%wete(i-1:i+1,:),par%dt,s%dsz)
            call advecyho(s%rr(i,:,:),s%cy(i,:,:),yradvec(i,:,:),                                 &
            0,s%ny,s%ntheta,s%dsv(i,:),s%dnv(i,:),s%dsdnzi(i,:),SCHEME_UPWIND_1, &
            s%wete(i,:),par%dt,s%dnz)
            call advecthetaho(s%rr(i,:,:),s%ctheta(i,:,:),thetaradvec(i,:,:),0,s%ny,s%ntheta,s%dtheta,par%scheme,s%wete(i,:))

            s%rr(i,:,:)=s%rr(i,:,:)-dtw*(xradvec(i,:,:)+yradvec(i,:,:)+thetaradvec(i,:,:))
            s%rr(i,:,:)=max(s%rr(i,:,:),0.0d0)
#ifdef USEMPI
            call xmpi_shift(s%rr(i-1:i,:,:),SHIFT_Y_R,1,2)
            call xmpi_shift(s%rr(i-1:i,:,:),SHIFT_Y_L,3,4)
#endif
            !
            ! euler step roller energy dissipation (source and sink function)
            do j=jmin_ee,jmax_ee
               do itheta=1,s%ntheta
                  if (dtw*dd(i,j,itheta)>s%ee(i,j,itheta)) then
                     dtw=min(dtw,.5*s%ee(i,j,itheta)/dd(i,j,itheta))
                  endif
               enddo
            enddo
            !Dano: need to make sure all processes use the same dtw, min of all processes
#ifdef USEMPI
            call xmpi_allreduce(dtw,MPI_MIN)
#endif
            do j=1,s%ny+1
               do itheta=1,s%ntheta
                  if(s%wete(i,j)==1) then
                     s%ee(i,j,itheta)=s%ee(i,j,itheta)-dtw*dd(i,j,itheta)
                     if (par%roller==1) then  !Christophe
                        s%rr(i,j,itheta)=s%rr(i,j,itheta)+dtw*dder(i,j,itheta)&
                        -dtw*2.*par%g*par%beta*s%rr(i,j,itheta)&
                        /sqrt(s%cx(i,j,itheta)**2+s%cy(i,j,itheta)**2)
                        drr(i,j,itheta) = 2*par%g*par%beta*max(s%rr(i,j,itheta),0.0d0)/           &
                        sqrt(s%cx(i,j,itheta)**2 +s%cy(i,j,itheta)**2)
                     else
                        s%rr(i,j,itheta)= 0.0d0
                        drr(i,j,itheta)= 0.0d0
                     end if
                     s%ee(i,j,itheta)=max(s%ee(i,j,itheta),0.0d0)
                     s%rr(i,j,itheta)=max(s%rr(i,j,itheta),0.0d0)
                  else
                     s%ee(i,j,itheta)=0.0d0
                     s%rr(i,j,itheta)=0.0d0
                     drr(i,j,itheta)=0.0d0
                  end if
               end do
            end do
            ! Lateral boundary condition
            if (xmpi_isleft .and. s%ny>0) then
               do itheta=1,s%ntheta
                  if (s%sinth(i,1,itheta)>=0.) then
                     s%ee(i,1,itheta)=s%ee(i,2,itheta)
                     s%rr(i,1,itheta)=s%rr(i,2,itheta)
                  endif
               enddo
               s%k(:,1)=s%k(:,2)
               s%sigm(:,1)=s%sigm(:,2)
            endif
            if (xmpi_isright .and. s%ny>0) then
               do itheta=1,s%ntheta
                  if (s%sinth(i,s%ny+1,itheta)<=0.) then
                     s%ee(i,s%ny+1,itheta)=s%ee(i,s%ny,itheta)
                     s%rr(i,s%ny+1,itheta)=s%rr(i,s%ny,itheta)
                  endif
               end do
               s%k(:,s%ny+1)=s%k(:,s%ny)
               s%sigm(:,s%ny+1)=s%sigm(:,s%ny)
            endif
            !
            ! Compute mean wave direction
            !
            if (par%snells==0) then
               s%thetamean(i,:)=(sum(s%ee(i,:,:)*s%thet(i,:,:),2)/size(s%ee(i,:,:),2)) &
               /(max(sum(s%ee(i,:,:),2),0.000010d0) /size(s%ee(i,:,:),2))
            endif
            !
            ! Energy integrated over wave directions,Hrms
            !
            s%E(i,:)=sum(s%ee(i,:,:),2)*s%dtheta
            s%R(i,:)=sum(s%rr(i,:,:),2)*s%dtheta
            s%DR(i,:)=sum(drr(i,:,:),2)*s%dtheta
            s%H(i,:)=sqrt(s%E(i,:)/par%rhog8)
            Herr=maxval(abs(Hprev(jmin_ee:jmax_ee)-s%H(i,jmin_ee:jmax_ee)))
#ifdef USEMPI
            call xmpi_allreduce(Herr,MPI_MAX)
#endif
            ! Stopping criteria
            if (iter<par%maxiter) then
               if (Herr<par%maxerror) then
                  stopiterate=.true.
                  !if(xmaster) call writelog('ls','(a,i4,a,i4)','Wave propagation row ',i,', iteration ',iter)
               endif
            else
               stopiterate=.true.
               if(xmaster) call writelog('ls','(a,i4,a,i4,a,f5.4)','Wave propagation row ',i,', iteration ',iter,', error: ',Herr)
            endif
         enddo ! End while loop
      enddo ! End do i=2:s%nx loop

      s%ee(s%nx+1,:,:) = s%ee(s%nx,:,:)
      s%rr(s%nx+1,:,:) = s%rr(s%nx,:,:)
      s%E(s%nx+1,:)    = s%E(s%nx,:)
      s%R(s%nx+1,:)    = s%R(s%nx,:)
      s%DR(s%nx+1,:)   = s%DR(s%nx,:)
      s%H(s%nx+1,:)    = s%H(s%nx,:)
      s%k(s%nx+1,:)    = s%k(s%nx,:)
      s%sigm(s%nx+1,:) = s%sigm(s%nx,:)
      s%cg(s%nx+1,:)   = s%cg(s%nx,:)
      s%c(s%nx+1,:)    = s%c(s%nx,:)
      s%thet(s%nx+1,:,:) = s%thet(s%nx,:,:)

      !
      ! Radiation stresses and forcing terms
      !
      s%n=s%cg/max(s%c,1d-10)
      s%Sxx=(s%n*sum((1.d0+s%costh**2.d0)*s%ee,3)-.5d0*sum(s%ee,3))*s%dtheta
      s%Syy=(s%n*sum((1.d0+s%sinth**2.d0)*s%ee,3)-.5d0*sum(s%ee,3))*s%dtheta
      s%Sxy=s%n*sum(s%sinth*s%costh*s%ee,3)*s%dtheta

      ! add roller contribution

      s%Sxx = s%Sxx + sum((s%costh**2)*s%rr,3)*s%dtheta
      s%Syy = s%Syy + sum((s%sinth**2)*s%rr,3)*s%dtheta
      s%Sxy = s%Sxy + sum(s%sinth*s%costh*s%rr,3)*s%dtheta

      if (s%ny>0) then
         do j=2,s%ny
            do i=1,s%nx
               s%Fx(i,j)=-(s%Sxx(i+1,j)-s%Sxx(i,j))/s%dsu(i,j)                        &
               -(s%Sxy(i,j+1)+s%Sxy(i+1,j+1)-s%Sxy(i,j-1)-s%Sxy(i+1,j-1))/    &
               (s%dnv(i,j-1)+s%dnv(i,j)+s%dnv(i+1,j-1)+s%dnv(i+1,j))
            enddo
         enddo

         do j=1,s%ny
            do i=2,s%nx
               s%Fy(i,j)=-(s%Syy(i,j+1)-s%Syy(i,j))/s%dnv(i,j)            &
               -(s%Sxy(i+1,j)+s%Sxy(i+1,j+1)-s%Sxy(i-1,j)-s%Sxy(i-1,j+1))/                    &
               (s%dsu(i-1,j)+s%dsu(i,j)+s%dsu(i-1,j+1)+s%dsu(i,j+1))
            enddo
         enddo
      else
         j=1
         do i=1,s%nx
            s%Fx(i,j)=-(s%Sxx(i+1,j)-s%Sxx(i,j))/s%dsu(i,j)
         enddo
         do i=2,s%nx
            s%Fy(i,j)=-(s%Sxy(i+1,j)-s%Sxy(i-1,j))/ (s%dsu(i-1,j)+s%dsu(i,j))
         enddo
      endif

      if (s%ny>0) then
         s%Fx(:,1)=s%Fx(:,2)
         ! Fy(:,1)=Fy(:,2)
         s%Fx(:,s%ny+1)=s%Fx(:,s%ny)
         ! Fy(:,ny+1)=Fy(:,ny)
      endif
      s%Fx(1,:)=s%Fx(2,:)
      ! Fy(1,:)=Fy(2,:)

      s%urms=uorb/sqrt(2.d0)

      ! wave induced mass flux
      ustw=s%E*s%k/s%sigm/par%rho/max(s%hh,.001d0)   ! Robert: why different from wave_instationary ???
      s%uwf = ustw*dcos(s%thetamean-s%alfaz)
      s%vwf = ustw*dsin(s%thetamean-s%alfaz)
      
      ! roller contribution
      s%ustr=2.*s%R*s%k/s%sigm/par%rho/max(s%hh,.001d0) ! Robert: why different from wave_instationary ???
           
      !ust = usd
      s%ust=s%usd+ustw  ! Robert: why different from wave_instationary ???
      !lateral boundaries
      s%ust(1,:) = s%ust(2,:)
      if (s%ny>0) then
         s%ust(:,1) = s%ust(:,2)
         s%ust(:,s%ny+1) = s%ust(:,s%ny)
      endif

   end subroutine wave_stationary
end module wave_stationary_module
