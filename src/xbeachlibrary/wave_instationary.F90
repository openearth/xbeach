module wave_instationary_module

contains

  subroutine wave_instationary(s,par)

    use params
    use spaceparams
    use roelvink_module
    use wave_functions_module
    use xmpi_module
    use mnemmodule
    use interp
    use paramsconst

    IMPLICIT NONE

    type(spacepars), target     :: s
    type(parameters)            :: par

    integer                     :: i
    integer                     :: j
    integer                     :: itheta
    integer                     :: dummy



    integer, dimension(:,:,:),allocatable,save  :: wete
    real*8 , dimension(:,:)  ,allocatable,save  :: dhdx,dhdy,dudx,dudy,dvdx,dvdy,ustw,Erfl
    real*8 , dimension(:,:)  ,allocatable,save  :: km,kmx,kmy,xwadvec,ywadvec,sinh2kh !,wm
    real*8 , dimension(:,:,:),allocatable,save  :: xadvec,yadvec,thetaadvec,dd,drr,dder
    real*8 , dimension(:,:,:),allocatable,save  :: xradvec,yradvec,thetaradvec
    real*8 , dimension(:,:)  ,allocatable,save  :: dkmxdx,dkmxdy,dkmydx,dkmydy,cgxm,cgym,arg,fac
    real*8 , dimension(:,:)  ,allocatable,save  :: wcifacu,wcifacv,hrmsold,uorb
    real*8 , dimension(:)    ,allocatable,save  :: wcrestpos
    real*8                                      :: factime,cs,sn,coffshore
    real*8 , save                               :: waverr  


    !include 's.ind'
    !include 's.inp'

    if (.not. allocated(wete)) then
       allocate(drr         (s%nx+1,s%ny+1,s%ntheta))
       allocate(wete        (s%nx+1,s%ny+1,s%ntheta))
       allocate(xadvec      (s%nx+1,s%ny+1,s%ntheta))
       allocate(yadvec      (s%nx+1,s%ny+1,s%ntheta))
       allocate(thetaadvec  (s%nx+1,s%ny+1,s%ntheta))
       allocate(xradvec     (s%nx+1,s%ny+1,s%ntheta))
       allocate(yradvec     (s%nx+1,s%ny+1,s%ntheta))
       allocate(thetaradvec (s%nx+1,s%ny+1,s%ntheta))
       allocate(dd          (s%nx+1,s%ny+1,s%ntheta))
       allocate(dder        (s%nx+1,s%ny+1,s%ntheta))

       allocate(dhdx        (s%nx+1,s%ny+1))
       allocate(dhdy        (s%nx+1,s%ny+1))
       allocate(dudx        (s%nx+1,s%ny+1))
       allocate(dudy        (s%nx+1,s%ny+1))
       allocate(dvdx        (s%nx+1,s%ny+1))
       allocate(dvdy        (s%nx+1,s%ny+1))
       allocate(km          (s%nx+1,s%ny+1))
       allocate(kmx         (s%nx+1,s%ny+1))
       allocate(kmy         (s%nx+1,s%ny+1))
       !       allocate(wm          (nx+1,ny+1))
       allocate(ustw        (s%nx+1,s%ny+1))
       allocate(Erfl        (s%nx+1,s%ny+1)) ! wwvv not used
       allocate(xwadvec     (s%nx+1,s%ny+1))
       allocate(ywadvec     (s%nx+1,s%ny+1))
       allocate(sinh2kh     (s%nx+1,s%ny+1))
       allocate(dkmxdx      (s%nx+1,s%ny+1))
       allocate(dkmxdy      (s%nx+1,s%ny+1))
       allocate(dkmydx      (s%nx+1,s%ny+1))
       allocate(dkmydy      (s%nx+1,s%ny+1))
       allocate(cgxm        (s%nx+1,s%ny+1))
       allocate(cgym        (s%nx+1,s%ny+1))
       allocate(arg         (s%nx+1,s%ny+1))
       allocate(fac         (s%nx+1,s%ny+1))
       allocate(wcifacu     (s%nx+1,s%ny+1))
       allocate(wcifacv     (s%nx+1,s%ny+1))
       allocate(hrmsold     (s%nx+1,s%ny+1))
       allocate(uorb        (s%nx+1,s%ny+1))
       allocate(wcrestpos   (s%nx+1))


       ! wwvv todo: I think these iniailization are superfluous
       drr         = 0.d0
       wete        = 0
       xadvec      = 0.d0
       yadvec      = 0.d0
       thetaadvec  = 0.d0
       xradvec     = 0.d0
       yradvec     = 0.d0
       thetaradvec = 0.d0
       dd          = 0.d0
       dder        = 0.d0
       dhdx        = 0.d0
       dhdy        = 0.d0
       dudx        = 0.d0
       dudy        = 0.d0
       dvdx        = 0.d0
       dvdy        = 0.d0
       km          = 0.d0
       kmx         = 0.d0
       kmy         = 0.d0
       !     wm          = 0.d0
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
       uorb        = 0.d0
       s%Fx          = 0.d0 ! in spacepars
       s%Fy          = 0.d0 ! in spacepars
    endif

    s%hh = max(s%hh,par%eps)
    hrmsold=s%H
    ! Calculate once velocities used with and without wave current interaction
    wcifacu=s%u*par%wci*min(s%hh/par%hwci,1.d0)
    wcifacv=s%v*par%wci*min(s%hh/par%hwci,1.d0)

    if (par%single_dir==1) then
       s%costh(:,:,1)=cos(s%thetamean-s%alfaz)
       s%sinth(:,:,1)=sin(s%thetamean-s%alfaz)
        
    elseif (par%snells==0) then
       s%thetamean=(sum(s%ee*s%thet,3)/s%ntheta)/(max(sum(s%ee,3),0.00001d0)/s%ntheta)
    elseif(par%snells==1) then  !Dano: Snellius
       ! Check for borderline cases where critical c/c(1,1) is reached....
#ifdef USEMPI
       if (xmpi_istop .and. xmpi_isleft) then 
          coffshore = s%c(1,1)
       else
          coffshore = -huge(0.d0)
       endif
       call xmpi_allreduce(coffshore,MPI_MAX)
#else
       coffshore = s%c(1,1)
#endif       
       s%thetamean=asin(max(-1.0d0, min(1.0d0, sin(s%theta0-s%alfaz(1,1))*s%c/coffshore)))+s%alfaz(1,1)
       s%costh(:,:,1)=cos(s%thetamean-s%alfaz)
       s%sinth(:,:,1)=sin(s%thetamean-s%alfaz)
!       thetamean = modulo(thetamean,2*par%px)
!       costh(:,:,1) = modulo(costh(:,:,1),2*par%px)
!       sinth(:,:,1) = modulo(sinth(:,:,1),2*par%px)
    endif

    ! Dispersion relation
    if (par%wci .ne. 0) then
       if (par%t==par%dt) then
          s%sigm = max((sum(s%sigt,3)/s%ntheta),0.01d0)
          call dispersion(par,s)
          s%umwci = 0.d0
          s%vmwci = 0.d0
          s%zswci = s%zs
          km=s%k
       endif
       km(1,:) = s%k(1,:)   ! boundary condition *assuming zero flow at the boundary) 
       factime = 1.d0/par%cats/par%Trep*par%dt
       s%umwci   = factime*s%uu + (1-factime)*s%umwci
       s%vmwci   = factime*s%vv + (1-factime)*s%vmwci
       s%zswci   = factime*s%zs + (1-factime)*s%zswci
       arg     = min(100.0d0,km*max(s%hh,par%delta*s%H))
       s%sigm(1,:) = sqrt( par%g*km(1,:)*tanh(arg(1,:))) ! *( 1.d0+ ((km(1,:)*s%H(1,:)/2.d0)**2)))
       !  calculate change in intrinsic frequency
       kmx = km*dcos(s%thetamean)
       kmy = km*dsin(s%thetamean)
       s%wm = s%sigm+kmx*s%umwci*par%wci*min((s%zswci-s%zb)/par%hwci,1.d0)+kmy*s%vmwci*par%wci*min((s%zswci-s%zb)/par%hwci,1.d0)

       cgym = s%cg*dsin(s%thetamean) + s%vmwci*min((s%zswci-s%zb)/par%hwci,1.d0)
       cgxm = s%cg*dcos(s%thetamean) + s%umwci*min((s%zswci-s%zb)/par%hwci,1.d0)

       call slope2D(kmx,s%nx,s%ny,s%dsu,s%dnv,dkmxdx,dkmxdy)
       call slope2D(kmy,s%nx,s%ny,s%dsu,s%dnv,dkmydx,dkmydy)
       call advecwx(s%wm,xwadvec,kmx,s%nx,s%ny,s%dsu)   ! cjaap: s%xz or s%xu?
       kmx = kmx -par%dt*xwadvec  -par%dt*cgym*(dkmydx-dkmxdy)
       if (s%ny>0) then 
         if (xmpi_isright) kmx(:,s%ny+1) = kmx(:,s%ny)  ! lateral bc
         if (xmpi_isleft)  kmx(:,1) = kmx(:,2)  ! lateral bc
       endif

       call advecwy(s%wm,ywadvec,kmy,s%nx,s%ny,s%dnv)   ! cjaap: s%yz or s%yv?
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
       ! non-linear dispersion
       arg = min(100.0d0,km*((s%zswci-s%zb)+par%delta*s%H))
       arg = max(arg,0.0001)
       !       fac = ( 1.d0 + ((km*H/2.d0)**2)*( (8.d0+(cosh(min(4.d0*arg,10.0d0)))**1.d0-2.d0*(tanh(arg))**2.d0 ) /(8.d0*(sinh(arg))**4.d0) ) )
       fac = ( 1.d0 + ((km*s%H/2.d0)**2))  ! use deep water correction instead of expression above (waves are short near blocking point anyway)
       !       fac = 1.d0    ! Linear
       s%sigm = sqrt( par%g*km*tanh(arg)*fac)

       !  update intrinsic frequency
       do itheta=1,s%ntheta
          s%sigt(:,:,itheta) = s%sigm
       enddo
       where(km>0.01d0)
          s%c  = s%sigm/km
          !          cg = c*(0.5d0+arg/sinh(2.0d0*arg))    ! Linear
          s%cg = s%c*(0.5d0+arg/sinh(2*arg))*sqrt(fac)  ! &  to include more
          !		 	          + km*(H/2)**2*sqrt(max(par%g*km*tanh(arg),0.001d0))/sqrt(max(fac,0.001d0)) ! include wave steepness
          s%n=0.5d0+km*s%hh/sinh(2*max(km,0.00001d0)*s%hh)
       elsewhere
          s%c  = 0.01d0
          s%cg = 0.01d0
          s%n  = 1.d0
       endwhere
       !  update k
       km = min(km,25.d0) ! limit to gravity waves
       s%k = km

    else  ! no wave current interaction
       s%sigm = max((sum(s%sigt,3)/s%ntheta),0.01d0)
       call dispersion(par,s)
    endif ! end wave current interaction

    ! Slopes of water depth
    call slope2D(max(s%hh,par%delta*s%H),s%nx,s%ny,s%dsu,s%dnv,dhdx,dhdy)
    call slope2D(wcifacu,s%nx,s%ny,s%dsu,s%dnv,dudx,dudy)
    call slope2D(wcifacv,s%nx,s%ny,s%dsu,s%dnv,dvdx,dvdy)
    !
    ! Calculate once sinh(2kh)
    where(2*s%hh*s%k<=3000.d0)
       sinh2kh=sinh(min(2*s%k*max(s%hh,par%delta*s%H),10.0d0))
    elsewhere
       sinh2kh = 3000.d0
    endwhere

    ! split wave velocities in wave grid directions theta
    do j=1,s%ny+1
       do i=1,s%nx+1
          do itheta=1,s%ntheta
             ! wave grid directions theta with respect to spatial grid s,n
             cs=s%costh(i,j,itheta)
             sn=s%sinth(i,j,itheta)

             ! split wave velocities over theta bins
             s%cgx(i,j,itheta)= s%cg(i,j)*cs+wcifacu(i,j)
             s%cgy(i,j,itheta)= s%cg(i,j)*sn+wcifacv(i,j)
             s%cx(i,j,itheta) =  s%c(i,j)*cs+wcifacu(i,j)
             s%cy(i,j,itheta) =  s%c(i,j)*sn+wcifacv(i,j)

             ! compute refraction velocity
             s%ctheta(i,j,itheta)=                                             &
                  s%sigm(i,j)/sinh2kh(i,j)*(dhdx(i,j)*sn-dhdy(i,j)*cs) + par%wci*   &
                  ( cs*(sn*dudx(i,j) - cs*dudy(i,j)) +                            &
                  sn*(sn*dvdx(i,j) - cs*dvdy(i,j)) )
          enddo
       enddo
    enddo
    ! Dano Limit unrealistic refraction speed to 1/2 pi per wave period
    s%ctheta=sign(1.d0,s%ctheta)*min(abs(s%ctheta),.5*par%px/par%Trep)
    !
    ! transform to wave action
    !
    s%ee = s%ee/s%sigt
    !
    ! Upwind Euler timestep propagation
    !
    call advecxho(s%ee,s%cgx,xadvec,s%nx,s%ny,s%ntheta,s%dnu,s%dsu,s%dsdnzi,par%scheme)
    if (s%ny>0) then
       call advecyho(s%ee,s%cgy,yadvec,s%nx,s%ny,s%ntheta,s%dsv,s%dnv,s%dsdnzi,par%scheme)
    endif
    !call advectheta(ee*ctheta,thetaadvec,nx,ny,ntheta,dtheta)
    call advecthetaho(s%ee,s%ctheta,thetaadvec,s%nx,s%ny,s%ntheta,s%dtheta,par%scheme)!
    s%ee=s%ee-par%dt*(xadvec+yadvec+thetaadvec)
    !
    ! transform back to wave energy
    !
    s%ee = s%ee*s%sigt
    s%ee=max(s%ee,0.0d0) !Jaap
    !
    ! Energy integrated over wave directions,Hrms
    !
    s%E=sum(s%ee,3)*s%dtheta
    s%H=sqrt(s%E/par%rhog8)
    do itheta=1,s%ntheta
       s%ee(:,:,itheta)=s%ee(:,:,itheta)/max(1.d0,(s%H/(par%gammax*s%hh))**2)
    enddo
    s%H=min(s%H,par%gammax*s%hh)
    s%E=par%rhog8*s%H**2

    ! Total dissipation

    select case(par%break)
      case(BREAK_ROELVINK1,BREAK_ROELVINK2)
       call roelvink(par,s,km)
      case(BREAK_BALDOCK)
       call baldock(par,s,km)
      case(BREAK_JANSSEN)
       call janssen_battjes(par,s,km)
      case(BREAK_ROELVINK_DALY)
       cgxm = s%c*dcos(s%thetamean) 
       cgym = s%c*dsin(s%thetamean)
       call advecqx(cgxm,s%Qb,xwadvec,s%nx,s%ny,s%dsu)
       if (s%ny>0) then
          call advecqy(cgym,s%Qb,ywadvec,s%nx,s%ny,s%dnv)
          s%Qb  = s%Qb-par%dt*(xwadvec+ywadvec)
       else
          s%Qb  = s%Qb-par%dt*xwadvec
       endif
       call roelvink(par,s,km)        
    end select

    ! Dissipation by bed friction
    uorb=par%px*s%H/par%Trep/sinh(min(max(s%k,0.01d0)*max(s%hh,par%delta*s%H),10.0d0))
    s%Df=0.21d0*par%rho*par%fw*uorb**3
    where (s%hh>par%fwcutoff)
       s%Df = 0.d0
    end where
    !
    ! Distribution of dissipation over directions and frequencies
    !
    do itheta=1,s%ntheta
       ! Only calculate for E>0 FB
       ! First just the dissipation that is fed to the roller
       dder(:,:,itheta)=s%ee(:,:,itheta)*s%D/max(s%E,0.00001d0)  
       ! Then all short wave energy dissipation, including bed friction and vegetation
       dd(:,:,itheta)=dder(:,:,itheta) + s%ee(:,:,itheta)*(s%Df+s%Dveg)/max(s%E,0.00001d0)
    enddo

    do j=1,s%ny+1
       do i=1,s%nx+1
          ! cjaap: replaced par%hmin by par%eps
          if(s%hh(i,j)+par%delta*s%H(i,j)>par%eps) then
             wete(i,j,1:s%ntheta)=1
          else
             wete(i,j,1:s%ntheta)=0
          end if
       end do
    end do
    !
    ! Euler step dissipation
    !
    ! calculate roller energy balance
    !
    call advecxho(s%rr,s%cx,xradvec,s%nx,s%ny,s%ntheta,s%dnu,s%dsu,s%dsdnzi,par%scheme)
    if (s%ny>0) then
       call advecyho(s%rr,s%cy,yradvec,s%nx,s%ny,s%ntheta,s%dsv,s%dnv,s%dsdnzi,par%scheme)
    endif
    !call advectheta(rr*ctheta,thetaradvec,nx,ny,ntheta,dtheta)
    call advecthetaho(s%rr,s%ctheta,thetaradvec,s%nx,s%ny,s%ntheta,s%dtheta,par%scheme)
    s%rr=s%rr-par%dt*(xradvec+yradvec+thetaradvec)
    s%rr=max(s%rr,0.0d0)
    !
    ! euler step roller energy dissipation (source and sink function)
    !
    do itheta=1,s%ntheta
       do j=1,s%ny+1
          do i=1,s%nx+1
             if(wete(i,j,itheta)==1) then
                s%ee(i,j,itheta)=s%ee(i,j,itheta)-par%dt*dd(i,j,itheta)
                if(par%roller==1) then
                   drr(i,j,itheta) = 2*par%g*s%BR(i,j)*max(s%rr(i,j,itheta),0.0d0)/   &
                        sqrt(s%cx(i,j,itheta)**2 +s%cy(i,j,itheta)**2)
                   s%rr(i,j,itheta)=s%rr(i,j,itheta)+par%dt*(dder(i,j,itheta)         &  ! Robert: changed from dd to dder 
                        -drr(i,j,itheta))                                            ! (only from wave s%breaking,
                else if (par%roller==0) then                                         !  not vegetation or bed friction)
                   s%rr(i,j,itheta)= 0.0d0
                   drr(i,j,itheta)= 0.0d0
                endif
                s%ee(i,j,itheta)=max(s%ee(i,j,itheta),0.0d0)
                s%rr(i,j,itheta)=max(s%rr(i,j,itheta),0.0d0)
             elseif(wete(i,j,itheta)==0) then
                s%ee(i,j,itheta)=0.0d0
                s%rr(i,j,itheta)=0.0d0
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

    if (xmpi_isbot) then
       s%ee(s%nx+1,:,:) =s%ee(s%nx,:,:)
       s%rr(s%nx+1,:,:) =s%rr(s%nx,:,:)
    endif

    if (par%t>0.0d0) then
       if (s%ny>0) then
          if (xmpi_isleft)then ! Jaap
             if (par%lateralwave==LATERALWAVE_NEUMANN) then
                !
                ! Lateral boundary at y=0;
                !
                s%ee(1:s%nx+1,1,:)=s%ee(1:s%nx+1,2,:)
                s%rr(1:s%nx+1,1,:)=s%rr(1:s%nx+1,2,:)
             elseif (par%lateralwave==LATERALWAVE_WAVECREST) then
                !   wcrestpos=xz+tan(thetamean(:,2))*(yz(2)-yz(1))
                wcrestpos=s%sdist(:,1)+tan(s%thetamean(:,2)-s%alfaz(:,2))*s%dnv(:,1)
                do itheta=1,s%ntheta
                   do i=1,s%nx+1
                      call Linear_interp((/-huge(0.d0)   ,wcrestpos     ,huge(0.d0)/),&
                           (/s%ee(1,2,itheta),s%ee(:,2,itheta),s%ee(s%nx+1,2,itheta)/),&
                           s%nx+1,s%sdist(i,1),s%ee(i,1,itheta),dummy)
                      call Linear_interp((/-huge(0.d0)   ,wcrestpos     ,huge(0.d0)/),&
                           (/s%rr(1,2,itheta),s%rr(:,2,itheta),s%rr(s%nx+1,2,itheta)/),&
                           s%nx+1,s%sdist(i,1),s%rr(i,1,itheta),dummy)
                      s%ee(i,1,itheta)=max(s%ee(i,1,itheta),0.d0)
                      s%rr(i,1,itheta)=max(s%rr(i,1,itheta),0.d0)
                   enddo
                enddo
             endif
          endif
          if (xmpi_isright)then
             if (par%lateralwave==LATERALWAVE_NEUMANN) then
                !
                ! lateral; boundary at y=ny*dy
                !
                s%ee(1:s%nx+1,s%ny+1,:)=s%ee(1:s%nx+1,s%ny,:)
                s%rr(1:s%nx+1,s%ny+1,:)=s%rr(1:s%nx+1,s%ny,:)
             elseif (par%lateralwave==LATERALWAVE_WAVECREST) then
                !  wcrestpos=xz-tan(thetamean(:,ny))*(yz(ny+1)-yz(ny))
                wcrestpos=s%sdist(:,s%ny+1)-tan(s%thetamean(:,s%ny)-s%alfaz(:,s%ny))*s%dnv(:,s%ny)
                do itheta=1,s%ntheta
                   do i=1,s%nx+1
                      call Linear_interp((/-huge(0.d0)   ,wcrestpos     ,huge(0.d0)/),&
                           (/s%ee(1,s%ny,itheta),s%ee(:,s%ny,itheta),s%ee(s%nx+1,s%ny,itheta)/),&
                           s%nx+1,s%sdist(i,s%ny+1),s%ee(i,s%ny+1,itheta),dummy)
                      call Linear_interp((/-huge(0.d0)   ,wcrestpos     ,huge(0.d0)/),&
                           (/s%rr(1,s%ny,itheta),s%rr(:,s%ny,itheta),s%rr(s%nx+1,s%ny,itheta)/),&
                           s%nx+1,s%sdist(i,s%ny+1),s%rr(i,s%ny+1,itheta),dummy)
                      s%ee(i,s%ny+1,itheta)=max(s%ee(i,s%ny+1,itheta),0.d0)
                      s%rr(i,s%ny+1,itheta)=max(s%rr(i,s%ny+1,itheta),0.d0)
                   enddo
                enddo
             endif
          endif
       endif
    endif
    ! wwvv communicate ee(:,1,:)
#ifdef USEMPI
    call xmpi_shift_ee(s%ee)
    call xmpi_shift_ee(s%rr)
#endif

    !
    ! Energy integrated over wave directions,Hrms
    !
    s%E  = sum(s%ee,3)*s%dtheta
    s%R  = sum(s%rr,3)*s%dtheta
    s%DR = sum(drr,3)*s%dtheta
    s%H  = sqrt(s%E/par%rhog8)
    waverr=sum(abs(s%H-hrmsold))/((s%nx+1)*(s%ny+1))

    !
    ! Compute mean wave direction
    !
    if (par%snells==0 .and. par%single_dir==0) then
       s%thetamean=(sum(s%ee*s%thet,3)/s%ntheta)/(max(sum(s%ee,3),0.00001d0)/s%ntheta)
    endif
    !
    ! Radiation stresses and forcing terms
    !
    ! n=cg/c   (Robert: calculated earlier in dispersion relation)
    s%Sxx=(s%n*sum((1.d0+s%costh**2)*s%ee,3)-.5d0*sum(s%ee,3))*s%dtheta
    s%Syy=(s%n*sum((1.d0+s%sinth**2)*s%ee,3)-.5d0*sum(s%ee,3))*s%dtheta
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
    ! wwvv in the previous, Fx and Fy are computed. The missing elements
    !  elements are Fx(:,1), Fx(nx+1,:), Fx(:,ny+1)
    !               Fy(1,:), Fy(nx+1,:), Fy(:,ny+1)

    ! wwvv so, Fx(:ny+1) and Fy(:ny+1) are left zero and Fx(nx+1,:) and Fy(nx+1,:)
    ! are made zero.  In the parallel case, Fx(:,1) and Fy(1,:) don't get a 
    ! value if the submatrices are not on suitable border. 
    ! I guess that it is necessary to communicate with neighbours the values of these elements


    ! wwvv todo the following has consequences for // version
    if(xmpi_istop) then
       s%Fy(1,:)=s%Fy(2,:)
    endif
    if(xmpi_isbot) then 
       s%Fx(s%nx+1,:) = 0.0d0
       s%Fy(s%nx+1,:) = 0.0d0
    endif
    if(xmpi_isleft .and. s%ny>0) then 
       s%Fx(:,1)=s%Fx(:,2)     
       !   Fy(:,1)=Fy(:,2)
       !   ! where (Fy(:,1)>0.d0) Fy(:,1)=0.d0; !Jaap + Bas: Don't do this before calling Dano :-)
    endif
    if (xmpi_isright .and. s%ny>0) then
       !   Fy(:,ny+1)=Fy(:,ny)   ! only a dummy point in non mpi
       s%Fx(:,s%ny+1)=s%Fx(:,s%ny)   ! only a dummy point in non mpi
       !   ! where (Fy(:,ny+1)<0.d0) Fy(:,ny+1)=0.d0; !Jaap + Bas: Don't do this before calling Dano :-)
    endif

    ! Ad
    !    urms=par%px*H/par%Trep/(sqrt(2.d0)*sinh(min(max(k,0.01d0)*max(hh,par%delta*H),10.0d0)))
    s%urms=uorb/sqrt(2.d0)
    !   urms=par%px*H/par%Trep/(sqrt(2.d0)*sinh(k*max(hh,par%delta*H)))

    ustw= s%E/max(s%c,sqrt(par%hmin*par%g))/par%rho/max(s%hh,par%hmin)   ! Jaap
    s%uwf = ustw*cos(s%thetamean-s%alfaz)
    s%vwf = ustw*sin(s%thetamean-s%alfaz)
    ! roller contribution
    s%ustr=2.*s%R/max(s%c,sqrt(par%hmin*par%g))/par%rho/max(s%hh,par%hmin) ! Jaap
    ! introduce breaker delay
    if (par%breakerdelay == 1) then
       call breakerdelay(par,s)
       s%ust = ustw+s%usd
    else
       s%ust = ustw+s%ustr
    endif
    !lateral boundaries
    ! wwvv todo the following has consequences for // version
    if(xmpi_istop) then
       s%ust(1,:) = s%ust(2,:)
    endif
    if(xmpi_isleft.and.s%ny>0) then
       s%ust(:,1) = s%ust(:,2)
       
    endif
    if(xmpi_isright.and. s%ny>0) then
       s%ust(:,s%ny+1) = s%ust(:,s%ny)
    endif

  end subroutine wave_instationary

end module wave_instationary_module
