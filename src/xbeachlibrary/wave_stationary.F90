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
    integer, dimension(:,:,:),allocatable,save  :: wete
    real*8 , dimension(:,:)  ,allocatable,save  :: dhdx,dhdy,dudx,dudy,dvdx,dvdy,ustw
    real*8 , dimension(:,:)  ,allocatable,save  :: km
      real*8 , dimension(:,:)  ,allocatable,save  :: kmx,kmy,sinh2kh ! ,wm
    real*8 , dimension(:,:,:),allocatable,save  :: xadvec,yadvec,thetaadvec,dd,drr
    real*8 , dimension(:,:,:),allocatable,save  :: xradvec,yradvec,thetaradvec
    real*8 , dimension(:,:,:),allocatable,save  :: cgxu,cgyv,cxu,cyv
    real*8 , dimension(:),allocatable,save      :: Hprev
    real*8                                      :: Herr,dtw
    real*8 , dimension(:)  ,allocatable,save    :: dkmxdx,dkmxdy,dkmydx,dkmydy,cgxm,cgym,arg,fac,xwadvec,ywadvec
    real*8 , dimension(:,:),allocatable,save    :: wcifacu,wcifacv,uorb
    logical                                     :: stopiterate

    !include 's.ind'
    !include 's.inp'

    if (.not. allocated(wete)) then
       allocate(wete      (s%nx+1,s%ny+1,s%ntheta))
       allocate(xadvec    (s%nx+1,s%ny+1,s%ntheta))
       allocate(yadvec    (s%nx+1,s%ny+1,s%ntheta))
       allocate(thetaadvec(s%nx+1,s%ny+1,s%ntheta))
       allocate(xradvec    (s%nx+1,s%ny+1,s%ntheta))
       allocate(yradvec    (s%nx+1,s%ny+1,s%ntheta))
       allocate(thetaradvec(s%nx+1,s%ny+1,s%ntheta))
       allocate(dd        (s%nx+1,s%ny+1,s%ntheta))
       allocate(drr       (s%nx+1,s%ny+1,s%ntheta))

       allocate(cgxu        (s%nx+1,s%ny+1,s%ntheta))
       allocate(cgyv        (s%nx+1,s%ny+1,s%ntheta))
       allocate(cxu         (s%nx+1,s%ny+1,s%ntheta))
       allocate(cyv         (s%nx+1,s%ny+1,s%ntheta))

       allocate(dhdx(s%nx+1,s%ny+1))
       allocate(dhdy(s%nx+1,s%ny+1))
       allocate(dudx(s%nx+1,s%ny+1))
       allocate(dudy(s%nx+1,s%ny+1))
       allocate(dvdx(s%nx+1,s%ny+1))
       allocate(dvdy(s%nx+1,s%ny+1))
       allocate(km  (s%nx+1,s%ny+1))
       allocate(kmx (3,s%ny+1))
       allocate(kmy (3,s%ny+1))
       !   allocate(wm  (3,ny+1))
       allocate(ustw(s%nx+1,s%ny+1))
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
       allocate(uorb        (s%nx+1,s%ny+1))

    endif

    wete        = 0
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
    !wm          = 0.0d0
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

    arg         = 0.0d0
    fac         = 0.0d0
    uorb        = 0.0d0

    ! cjaap: replaced par%hmin by par%eps
    s%hh = max(s%hh,par%eps)


    ! Slopes of water depth
    call slope2D(max(s%hh,par%delta*s%H),s%nx,s%ny,s%dsu,s%dnv,dhdx,dhdy)
    ! Dano limit slopes used in refraction to avoid unrealistic refraction speeds
    dhdx=sign(1.d0,dhdx)*min(abs(dhdx),0.1d0)
    dhdy=sign(1.d0,dhdy)*min(abs(dhdy),0.1d0)
    call slope2D(s%u*par%wci,s%nx,s%ny,s%dsu,s%dnv,dudx,dudy)
    call slope2D(s%v*par%wci,s%nx,s%ny,s%dsu,s%dnv,dvdx,dvdy)

    ! wwvv these slope routines are in wave_timestep, and are
    !   MPI-aware
    !
    ! Calculate once sinh(2kh)
    where(2*s%hh*s%k<=3000.d0)
       sinh2kh=sinh(min(2*s%k*max(s%hh,par%delta*s%H),10.0d0))
    elsewhere
       sinh2kh = 3000.d0
    endwhere

    call dispersion(par,s)   

    if (par%wci==0) then
       s%sigm=2.d0*par%px/par%Trep
         do itheta=1,s%ntheta
          s%sigt(:,:,itheta) = max(s%sigm,0.010d0)
         enddo
    endif

    if (par%snells==0) then
       s%thetamean=(sum(s%ee*s%thet,3)/s%ntheta)/(max(sum(s%ee,3),0.00001d0)/s%ntheta)
    else !Dano: Snellius
       s%thetamean=asin(sin(s%theta0-s%alfaz(1,1))*s%c/s%c(1,1))+s%alfaz(1,1)
       s%costh(:,:,1)=cos(s%thetamean-s%alfaz)
       s%sinth(:,:,1)=sin(s%thetamean-s%alfaz)
    endif

    ! Calculate once velocities used with and without wave current interaction
    wcifacu=s%u*par%wci*min(s%hh/par%hwci,1.d0)
    wcifacv=s%v*par%wci*min(s%hh/par%hwci,1.d0)

      do itheta=1,s%ntheta
       s%cgx(:,:,itheta)= s%cg*s%costh(:,:,itheta)+wcifacu
       s%cgy(:,:,itheta)= s%cg*s%sinth(:,:,itheta)+wcifacv
       s%cx(:,:,itheta) =  s%c*s%costh(:,:,itheta)+wcifacu
       s%cy(:,:,itheta) =  s%c*s%sinth(:,:,itheta)+wcifacv
       s%ctheta(:,:,itheta)=  &
            s%sigm/sinh2kh*(dhdx*s%sinth(:,:,itheta)-dhdy*s%costh(:,:,itheta)) + &
            par%wci*(&
            s%costh(:,:,itheta)*(s%sinth(:,:,itheta)*dudx - s%costh(:,:,itheta)*dudy) + &
            s%sinth(:,:,itheta)*(s%sinth(:,:,itheta)*dvdx - s%costh(:,:,itheta)*dvdy))
      enddo
    ! Dano Limit unrealistic refraction speed to 1/2 pi per wave period
    s%ctheta=sign(1.d0,s%ctheta)*min(abs(s%ctheta),.5*par%px/par%Trep)

    km = s%k


    !dtw=.9*minval(xz(2:nx+1)-xz(1:nx))/maxval(cgx)

       s%E(1,:)=sum(s%ee(1,:,:),2)*s%dtheta
       s%R(1,:)=max(sum(s%rr(1,:,:),2)*s%dtheta,0.0d0)
       s%H(1,:)=sqrt(s%E(1,:)/par%rhog8)

       imax=s%nx
!Dano  This is ok, since we will set mpiboundary to y in stationary mode

       do i=2,imax
          dtw=.5*minval(s%dsu(i:i+1,jmin_ee:jmax_ee))/maxval(s%cgx(i-1:i+1,jmin_ee:jmax_ee,:))
          dtw=min(dtw,.5*minval(s%dnv(i,jmin_ee:jmax_ee))/maxval(abs(s%cgy(i,jmin_ee:jmax_ee,:))))
          dtw=min(dtw,.5*s%dtheta/max(1.0d-6,maxval(abs(s%ctheta(i,jmin_ee:jmax_ee,:)))))
!Dano: need to make sure all processes use the same dtw, min of all processes
#ifdef USEMPI
          call xmpi_allreduce(dtw,MPI_MIN)
#endif 
          Herr=1.
          iter=0
          arg = min(100.0d0,km(i,:)*(s%hh(i,:)+par%delta*s%H(i,:)))
          arg = max(arg,0.0001)
          fac = ( 1.d0 + ((km(i,:)*s%H(i,:)/2.d0)**2))  ! use deep water correction instead of expression above (waves are short near blocking point anyway)
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
               do itheta=1,s%ntheta
                   s%sigt(i,:,itheta) = max(s%sigm(i,:),0.010d0)
               enddo
             endif

             !
             ! transform to wave action
             !
             i1=max(i-2,1)
             s%ee(i1:i+1,:,:) = s%ee(i1:i+1,:,:)/s%sigt(i1:i+1,:,:)
             !
             ! Upwind Euler timestep propagation
             !
             if  (i>2.and. par%scheme==SCHEME_UPWIND_2) then
                call advecxho(s%ee(i-2:i+1,:,:),s%cgx(i-2:i+1,:,:),xadvec(i-2:i+1,:,:),    &
                     3,s%ny,s%ntheta,s%dnu(i-2:i+1,:),s%dsu(i-2:i+1,:),s%dsdnzi(i-2:i+1,:),SCHEME_UPWIND_2)
             else
                call advecxho(s%ee(i-1:i+1,:,:),s%cgx(i-1:i+1,:,:),xadvec(i-1:i+1,:,:),    &
                     2,s%ny,s%ntheta,s%dnu(i-1:i+1,:),s%dsu(i-1:i+1,:),s%dsdnzi(i-1:i+1,:),SCHEME_UPWIND_1)
             endif
             call advecyho(s%ee(i,:,:),s%cgy(i,:,:),yadvec(i,:,:),                                  &
                  0,s%ny,s%ntheta,s%dsv(i,:),s%dnv(i,:),s%dsdnzi(i,:),SCHEME_UPWIND_1)
             !call advecx(ee(i-1:i+1,:,:)*cgx(i-1:i+1,:,:),xadvec(i-1:i+1,:,:),2,ny,ntheta,dnu(i-1:i+1,:),dsdnzi(i-1:i+1,:))
             !call advecy(ee(i,:,:)*cgy(i,:,:),yadvec(i,:,:),0,ny,ntheta,dsv(i,:),dsdnzi(i,:))
             !call advectheta(ee(i,:,:)*ctheta(i,:,:),thetaadvec(i,:,:),0,ny,ntheta,dtheta)
             call advecthetaho(s%ee(i,:,:),s%ctheta(i,:,:),thetaadvec(i,:,:),0,s%ny,s%ntheta,s%dtheta,par%scheme)

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
                call roelvink       (par,s,km(i,:),i)
               case(BREAK_BALDOCK)
                call baldock        (par,s,km(i,:),i)
               case(BREAK_JANSSEN)
                call janssen_battjes(par,s,km(i,:),i)
             end select


!             if(trim(par%break)=='roelvink1' .or. trim(par%break)=='roelvink2') then
!                call roelvink       (par,s,km(i,:),i)
!             else if(trim(par%break)=='baldock') then
!                call baldock        (par,s,km(i,:),i)
!             elseif (trim(par%break)=='janssen') then
!                call janssen_battjes(par,s,km(i,:),i)
!             end if
             ! Dissipation by bed friction
             uorb(i,:)=par%px*s%H(i,:)/par%Trep/sinh(min(max(s%k(i,:),0.01d0)*max(s%hh(i,:),par%delta*s%H(i,:)),10.0d0))
             s%Df(i,:)=0.28d0*par%rho*par%fw*uorb(i,:)**3
             where (s%hh>par%fwcutoff)
                s%Df = 0.d0
             end where
             !
             ! Distribution of dissipation over directions and frequencies
             !
             do itheta=1,s%ntheta
                dd(i,:,itheta)=s%ee(i,:,itheta)*(s%D(i,:)+s%Df(i,:))/max(s%E(i,:),0.00001d0)
             end do
             do j=1,s%ny+1
                ! cjaap: replaced par%hmin by par%eps
                if(s%hh(i,j)+par%delta*s%H(i,j)>par%eps) then
                   wete(i,j,1:s%ntheta)=1
                else
                   wete(i,j,1:s%ntheta)=0
                end if
             end do
             !
             ! Euler step dissipation
             !
             ! calculate roller energy balance
             !
             call advecxho(s%rr(i-1:i+1,:,:),s%cx(i-1:i+1,:,:),xradvec(i-1:i+1,:,:),   &
                  2,s%ny,s%ntheta,s%dnu(i-1:i+1,:),s%dsu(i-1:i+1,:),s%dsdnzi(i-1:i+1,:),SCHEME_UPWIND_1)
             call advecyho(s%rr(i,:,:),s%cy(i,:,:),yradvec(i,:,:),                                 &
                  0,s%ny,s%ntheta,s%dsv(i,:),s%dnv(i,:),s%dsdnzi(i,:),SCHEME_UPWIND_1)
             !call advecx(rr(i-1:i+1,:,:)*cx(i-1:i+1,:,:),xradvec(i-1:i+1,:,:),2,ny,ntheta,dnu(i-1:i+1,:),dsdnzi(i-1:i+1,:)) !Robert & Jaap
             !call advecy(rr(i,:,:)*cy(i,:,:),yradvec(i,:,:),0,ny,ntheta,dsv(i,:),dsdnzi(i,:))                   !Robert & Jaap
             !call advectheta(rr(i,:,:)*ctheta(i,:,:),thetaradvec(i,:,:),0,ny,ntheta,dtheta)   !Robert & Jaap
             call advecthetaho(s%rr(i,:,:),s%ctheta(i,:,:),thetaradvec(i,:,:),0,s%ny,s%ntheta,s%dtheta,par%scheme)

             s%rr(i,:,:)=s%rr(i,:,:)-dtw*(xradvec(i,:,:) &
                  +yradvec(i,:,:) &
                  +thetaradvec(i,:,:))
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
                   if(wete(i,j,itheta)==1) then
                      s%ee(i,j,itheta)=s%ee(i,j,itheta)-dtw*dd(i,j,itheta)
                      if (par%roller==1) then  !Christophe
                         s%rr(i,j,itheta)=s%rr(i,j,itheta)+dtw*dd(i,j,itheta)&
                              -dtw*2.*par%g*par%beta*s%rr(i,j,itheta)&
                              /sqrt(s%cx(i,j,itheta)**2+s%cy(i,j,itheta)**2)
                         drr(i,j,itheta) = 2*par%g*par%beta*max(s%rr(i,j,itheta),0.0d0)/           &
                              sqrt(s%cx(i,j,itheta)**2 +s%cy(i,j,itheta)**2)
                      else if (par%roller==0) then
                         s%rr(i,j,itheta)= 0.0d0
                         drr(i,j,itheta)= 0.0d0
                      end if
                      s%ee(i,j,itheta)=max(s%ee(i,j,itheta),0.0d0)
                      s%rr(i,j,itheta)=max(s%rr(i,j,itheta),0.0d0)
                   else if(wete(i,j,itheta)==0) then
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
                km(:,1)=km(:,2)
                s%sigm(:,1)=s%sigm(:,2)
             endif
             if (xmpi_isright .and. s%ny>0) then
                do itheta=1,s%ntheta 
                   if (s%sinth(i,s%ny+1,itheta)<=0.) then
                      s%ee(i,s%ny+1,itheta)=s%ee(i,s%ny,itheta)
                      s%rr(i,s%ny+1,itheta)=s%rr(i,s%ny,itheta)
                   endif
                end do
                km(:,s%ny+1)=km(:,s%ny)
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
                   if(xmaster) call writelog('ls','(a,i4,a,i4)','Wave propagation row ',i,', iteration ',iter)
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
    km(s%nx+1,:)   = km(s%nx,:)
    s%sigm(s%nx+1,:) = s%sigm(s%nx,:)
    s%cg(s%nx+1,:)   = s%cg(s%nx,:)
    s%c(s%nx+1,:)    = s%c(s%nx,:)
    s%thet(s%nx+1,:,:) = s%thet(s%nx,:,:)


    s%k=km
    !
    ! Radiation stresses and forcing terms
    !
    s%n=s%cg/s%c
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

    ! wwvv
    
    s%urms=par%px*s%H/par%Trep/(sqrt(2.d0)*sinh(s%k*max(s%hh,par%delta*s%H)))
    !ust=E*k/sigm/par%rho/max(hh,0.001)

    ! wave induced mass flux
    ustw=s%E*s%k/s%sigm/par%rho/max(s%hh,.001d0)
    s%uwf = ustw*dcos(s%thetamean-s%alfaz)
    s%vwf = ustw*dsin(s%thetamean-s%alfaz)
    ! roller contribution
    s%ustr=2.*s%R*s%k/s%sigm/par%rho/max(s%hh,.001d0)
    ! introduce breaker delay
    if (par%breakerdelay == 1) call breakerdelay(par,s)
    !ust = usd
    s%ust=s%usd+ustw
    !lateral boundaries
    s%ust(1,:) = s%ust(2,:)
    if (s%ny>0) then
       s%ust(:,1) = s%ust(:,2)
       s%ust(:,s%ny+1) = s%ust(:,s%ny)
    endif
    ! wwvv
    !D=2*par%g*par%beta*sum(rr/sqrt(cx**2+cy**2),3)*dtheta ! Arnold: commented this line, seems like a copy/paste error.

  end subroutine wave_stationary
end module wave_stationary_module
