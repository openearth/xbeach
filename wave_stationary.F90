module wave_stationary_module
contains
  subroutine wave_stationary(s,par)
    use params
    use spaceparams
    use roelvink_module
    use wave_functions_module
    use xmpi_module
    use logging_module

    ! wwvv in my testcase, this routine was not called, so it is not
    ! tested. Nevertheless, I put in code for the parallel version.

    IMPLICIT NONE

    type(spacepars), target     :: s
    type(parameters)            :: par

    integer                     :: i,imax,it
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

#ifdef USEMPI
! some variables used for mpi only
    real*8                      :: imaxr
    integer                     :: iimpi
#endif

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
       
       allocate(cgxu        (nx+1,ny+1,ntheta))
       allocate(cgyv        (nx+1,ny+1,ntheta))
       allocate(cxu         (nx+1,ny+1,ntheta))
       allocate(cyv         (nx+1,ny+1,ntheta))
       
       allocate(dhdx(nx+1,ny+1))
       allocate(dhdy(nx+1,ny+1))
       allocate(dudx(nx+1,ny+1))
       allocate(dudy(nx+1,ny+1))
       allocate(dvdx(nx+1,ny+1))
       allocate(dvdy(nx+1,ny+1))
       allocate(km  (nx+1,ny+1))
       allocate(kmx (3,ny+1))
       allocate(kmy (3,ny+1))
       !   allocate(wm  (3,ny+1))
       allocate(ustw(nx+1,ny+1))
       allocate(xwadvec(ny+1))
       allocate(ywadvec(ny+1))
       allocate(sinh2kh(nx+1,ny+1))
       allocate(Hprev(ny+1))

       allocate(dkmxdx  (ny+1))
       allocate(dkmxdy  (ny+1))
       allocate(dkmydx  (ny+1))
       allocate(dkmydy  (ny+1))
       allocate(cgxm    (ny+1))
       allocate(cgym    (ny+1))
       allocate(arg     (ny+1))
       allocate(fac     (ny+1))
       allocate(wcifacu     (nx+1,ny+1))
       allocate(wcifacv     (nx+1,ny+1))
       allocate(uorb        (nx+1,ny+1))

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
    hh = max(hh,par%eps)


    ! Slopes of water depth
    call slope2D(max(hh,par%delta*H),nx,ny,dsu,dnv,dhdx,dhdy)
    call slope2D(u*par%wci,nx,ny,dsu,dnv,dudx,dudy)
    call slope2D(v*par%wci,nx,ny,dsu,dnv,dvdx,dvdy)

    ! wwvv these slope routines are in wave_timestep, and are
    !   MPI-aware
    !
    ! Calculate once sinh(2kh)
    where(2*hh*k<=3000.d0)
       sinh2kh=sinh(min(2*k*max(hh,par%delta*H),10.0d0))
    elsewhere
       sinh2kh = 3000.d0
    endwhere

    call dispersion(par,s)   

    if (par%wci==0) then
       sigm=2.d0*par%px/par%Trep
       DO itheta=1,ntheta
          sigt(:,:,itheta) = max(sigm,0.010d0)
       END DO
    endif

    thetamean=(sum(ee*thet,3)/size(ee,3))/(max(sum(ee,3),0.00001d0)/size(ee,3))
#ifdef USEMPI
    call xmpi_shift(thetamean,'1:')
#endif

    ! Calculate once velocities used with and without wave current interaction
    wcifacu=u*par%wci*min(hh/par%hwci,1.d0)
    wcifacv=v*par%wci*min(hh/par%hwci,1.d0)

    DO itheta=1,ntheta
       cgx(:,:,itheta)= cg*costh(:,:,itheta)+wcifacu
       cgy(:,:,itheta)= cg*sinth(:,:,itheta)+wcifacv
       cx(:,:,itheta) =  c*costh(:,:,itheta)+wcifacu
       cy(:,:,itheta) =  c*sinth(:,:,itheta)+wcifacv
       ctheta(:,:,itheta)=  &
            sigm/sinh2kh*(dhdx*sinth(:,:,itheta)-dhdy*costh(:,:,itheta)) + &
            par%wci*(&
            costh(:,:,itheta)*(sinth(:,:,itheta)*dudx - costh(:,:,itheta)*dudy) + &
            sinth(:,:,itheta)*(sinth(:,:,itheta)*dvdx - costh(:,:,itheta)*dvdy))
    END DO
    
    ! convert wave velocities from z to u points using the mean
    do j=1,ny+1
       do i=1,nx
          do itheta=1,ntheta
            cgxu(i,j,itheta) = 0.5*(cgx(i,j,itheta)+cgx(i+1,j,itheta))
            cxu(i,j,itheta) = 0.5*(cx(i,j,itheta)+cx(i+1,j,itheta))
          enddo
       enddo
    enddo
    cgxu(nx+1,:,:) = cgxu(nx,:,:)
    cxu(nx+1,:,:) = cxu(nx,:,:)
    
    do j=1,ny
       do i=1,nx+1
          do itheta=1,ntheta
            cgyv(i,j,itheta) = 0.5*(cgy(i,j,itheta)+cgy(i,j+1,itheta))
            cyv(i,j,itheta) = 0.5*(cy(i,j,itheta)+cy(i,j+1,itheta))
          enddo
       enddo
    enddo
    cgyv(:,ny+1,:) = cgyv(:,ny,:)
    cyv(:,ny+1,:) = cyv(:,ny,:)

    km = k

#ifdef USEMPI
    call xmpi_shift(km,'1:')
#endif

    !dtw=.9*minval(xz(2:nx+1)-xz(1:nx))/maxval(cgx)

#ifdef USEMPI
    do iimpi=1,xmpi_m    ! Start iteration for number of domains in x direction
       if (xmaster) call writelog('ls','(a,i0,a,i0,a)','Starting wave propagation over ',iimpi,' of ',xmpi_m,' MPI rows')
       call xmpi_shift(ee,'1:')
       call xmpi_shift(rr,'1:')
       call xmpi_shift(km,'1:')
       call xmpi_shift(sigm,'1:')
#endif 

       E(1,:)=sum(ee(1,:,:),2)*dtheta
       R(1,:)=max(sum(rr(1,:,:),2)*dtheta,0.0d0)
       H(1,:)=sqrt(E(1,:)/par%rhog8)

       imax=nx
#ifdef USEMPI
       imaxr=real(imax)
       call xmpi_allreduce(imaxr,MPI_MAX)
       imax=nint(imaxr)
#endif 

       do it=2,imax
          ! In the case of uneven grids, ensure that each process continues running
          ! until greatest grid is complete.
          ! Start shorter grid at offshore (+1) point again while larger grid finalizes last
          ! grid rows.
          if (it .le. nx) then
             i=it
          else
             i=it+1-nx
          endif
          !Dano
          dtw=.5*minval(dsu(i:i+1,:))/maxval(cgx(i-1:i+1,:,:))
          dtw=min(dtw,.5*minval(dnv(i,2:ny+1))/maxval(abs(cgy(i,:,:))))
          dtw=min(dtw,.5*dtheta/maxval(abs(ctheta(i,:,:))))
          Herr=1.
          iter=0
          arg = min(100.0d0,km(i,:)*(hh(i,:)+par%delta*H(i,:)))
          arg = max(arg,0.0001)
          fac = ( 1.d0 + ((km(i,:)*H(i,:)/2.d0)**2))  ! use deep water correction instead of expression above (waves are short near blocking point anyway)
          stopiterate=.false.
          do while (stopiterate .eqv. .false.)
             iter=iter+1
             Hprev=H(i,:)
             ! WCI
             if (par%wci==1) then
			    ! Dano NEED TO CHECK THIS FOR CURVI
                kmx = km(i-1:i+1,:)*cos(thetamean(i-1:i+1,:))
                kmy = km(i-1:i+1,:)*sin(thetamean(i-1:i+1,:))
                wm(i-1:i+1,:) = sigm(i-1:i+1,:)+kmx*wcifacu(i-1:i+1,:)&
                     +kmy*wcifacv(i-1:i+1,:)

                where(km(i,:)>0.01d0)
                   c(i,:)  = sigm(i,:)/km(i,:)
                   cg(i,:) = c(i,:)*(0.5d0+arg/sinh(2*arg))*sqrt(fac)  
                   n(i,:)  = 0.5d0+km(i,:)*hh(i,:)/sinh(2*max(km(i,:),0.00001d0)*hh(i,:))
                elsewhere
                   c(i,:)  = 0.01d0
                   cg(i,:) = 0.01d0
                   n(i,:)  = 1.d0
                endwhere

                cgym = cg(i,:)*dsin(thetamean(i,:)) + wcifacv(i,:)
                cgxm = cg(i,:)*dcos(thetamean(i,:)) + wcifacu(i,:)

                dkmxdx       = (kmx(3,:)-kmx(1,:))/(dsu(i,:)+dsu(i+1,:))
                dkmxdy(2:ny) = (kmx(2,3:ny+1)-kmx(2,1:ny-1))/(dnv(i,2:ny)+dnv(i,3:ny+1))
                dkmxdy(1)    = dkmxdy(2)
                dkmxdy(ny+1) = dkmxdy(ny)
                dkmydx       = (kmy(3,:)-kmy(1,:))/(dsu(i,:)+dsu(i+1,:))
                dkmydy(2:ny) = (kmy(2,3:ny+1)-kmy(2,1:ny-1))/(dnv(i,2:ny)+dnv(i,3:ny+1))
                dkmydy(1)    = dkmydy(2)
                dkmydy(ny+1) = dkmydy(ny)

                xwadvec  = (wm(i,:)-wm(i-1,:))/dsu(i-1,:)
                kmx(2,:) = kmx(2,:) -dtw*xwadvec -dtw*cgym*(dkmydx-dkmxdy)

                ywadvec(2:ny) = (wm(i,3:ny+1)-wm(i,1:ny-1))/(dnv(i,2:ny)+dnv(i,3:ny+1))
                ywadvec(1)=ywadvec(2)
                ywadvec(ny+1)=ywadvec(ny)
                kmy(2,:) = kmy(2,:) -dtw*ywadvec + dtw*cgxm*(dkmydx-dkmxdy)

                ! wwvv 
#ifdef USEMPI
                call xmpi_shift(kmx,':n')
                call xmpi_shift(kmy,':n')
                call xmpi_shift(kmx,':1')
                call xmpi_shift(kmy,':1')
#endif
                km(i,:) = sqrt(kmx(2,:)**2+kmy(2,:)**2)
                km(i,:) = min(km(i,:),25.d0) ! limit to gravity waves

                !  non-linear dispersion
                arg = min(100.0d0,km(i,:)*(hh(i,:)+par%delta*H(i,:)))
                arg = max(arg,0.0001)
                fac = ( 1.d0 + ((km(i,:)*H(i,:)/2.d0)**2)) 
                sigm(i,:) = sqrt( par%g*km(i,:)*tanh(arg)*fac)
                DO itheta=1,ntheta
                   sigt(i,:,itheta) = max(sigm(i,:),0.010d0)
                END DO
             endif

             !
             ! transform to wave action
             !
             ee(i-1:i+1,:,:) = ee(i-1:i+1,:,:)/sigt(i-1:i+1,:,:)
             !
             ! Upwind Euler timestep propagation
             !
             call advecxho(ee(i-1:i+1,:,:),cgx(i-1:i+1,:,:),xadvec(i-1:i+1,:,:),    &
                           2,ny,ntheta,dsu(i-1:i+1,:),dnu(i-1:i+1,:),dsdnzi(i-1:i+1,:),par%dt,'upwind_1')
             call advecyho(ee(i,:,:),cgy(i,:,:),yadvec(i,:,:),                                  &
                           0,ny,ntheta,dnv(i,:),dsv(i,:),dsdnzi(i,:),par%dt,'upwind_1')
             !call advecx(ee(i-1:i+1,:,:)*cgx(i-1:i+1,:,:),xadvec(i-1:i+1,:,:),2,ny,ntheta,dnu(i-1:i+1,:),dsdnzi(i-1:i+1,:))
             !call advecy(ee(i,:,:)*cgy(i,:,:),yadvec(i,:,:),0,ny,ntheta,dsv(i,:),dsdnzi(i,:))
             call advectheta(ee(i,:,:)*ctheta(i,:,:),thetaadvec(i,:,:),0,ny,ntheta,dtheta)

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
                ee(i,:,itheta)=ee(i,:,itheta)/max(1.0d0,(H(i,:)/(par%gammax*hh(i,:)))**2)
             enddo
             H(i,:)=min(H(i,:),par%gammax*hh(i,:))
             E(i,:)=par%rhog8*H(i,:)**2

             thetamean(i,:) = (sum(ee(i,:,:)*thet(i,:,:),2)/ntheta)/(max(sum(ee(i,:,:),2),0.000010d0)/ntheta)

             !
             ! Total dissipation
             if(trim(par%break)=='roelvink1' .or. trim(par%break)=='roelvink2') then
                call roelvink       (par,s,km(i,:),i)
             else if(trim(par%break)=='baldock') then
                call baldock        (par,s,km(i,:),i)
             elseif (trim(par%break)=='janssen') then
                call janssen_battjes(par,s,km(i,:),i)
             end if
             ! Dissipation by bed friction
             uorb(i,:)=par%px*H(i,:)/par%Trep/sinh(min(max(k(i,:),0.01d0)*max(hh(i,:),par%delta*H(i,:)),10.0d0))
             Df(i,:)=0.6666666d0/par%px*par%rho*par%fw*uorb(i,:)**3

             !
             ! Distribution of dissipation over directions and frequencies
             !
             do itheta=1,ntheta
                dd(i,:,itheta)=ee(i,:,itheta)*(D(i,:)+Df(i,:))/max(E(i,:),0.00001d0)
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
             call advecxho(rr(i-1:i+1,:,:),cx(i-1:i+1,:,:),xradvec(i-1:i+1,:,:),   &
                           2,ny,ntheta,dsu(i-1:i+1,:),dnu(i-1:i+1,:),dsdnzi(i-1:i+1,:),par%dt,'upwind_1')
             call advecyho(rr(i,:,:),cy(i,:,:),yradvec(i,:,:),                                 &
                           0,ny,ntheta,dnv(i,:),dsv(i,:),dsdnzi(i,:),par%dt,'upwind_1')
             !call advecx(rr(i-1:i+1,:,:)*cx(i-1:i+1,:,:),xradvec(i-1:i+1,:,:),2,ny,ntheta,dnu(i-1:i+1,:),dsdnzi(i-1:i+1,:)) !Robert & Jaap
             !call advecy(rr(i,:,:)*cy(i,:,:),yradvec(i,:,:),0,ny,ntheta,dsv(i,:),dsdnzi(i,:))                   !Robert & Jaap
             call advectheta(rr(i,:,:)*ctheta(i,:,:),thetaradvec(i,:,:),0,ny,ntheta,dtheta)   !Robert & Jaap

             rr(i,:,:)=rr(i,:,:)-dtw*(xradvec(i,:,:) &
                  +yradvec(i,:,:) &
                  +thetaradvec(i,:,:))
             rr(i,:,:)=max(rr(i,:,:),0.0d0)
             !
             ! euler step roller energy dissipation (source and sink function)
             do j=1,ny+1
                do itheta=1,ntheta        
                   if (dtw*dd(i,j,itheta)>ee(i,j,itheta)) then
                      dtw=min(dtw,.5*ee(i,j,itheta)/dd(i,j,itheta))                
                   endif
                enddo
             enddo
             do j=1,ny+1
                do itheta=1,ntheta                        
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
             if (xmpi_isleft) then
                do itheta=1,ntheta
                   if (sinth(i,1,itheta)>=0.) then
                      ee(i,1,itheta)=ee(i,2,itheta)
                      rr(i,1,itheta)=rr(i,2,itheta)
                   endif
                enddo
                km(:,1)=km(:,2)
                sigm(:,1)=sigm(:,2)
             endif
             if (xmpi_isright) then
                do itheta=1,ntheta 
                   if (sinth(i,1,itheta)<=0.) then
                      ee(i,ny+1,itheta)=ee(i,ny,itheta)
                      rr(i,ny+1,itheta)=rr(i,ny,itheta)
                   endif
                end do
                km(:,ny+1)=km(:,ny)
                sigm(:,ny+1)=sigm(:,ny)
             endif
             ! Pass internal boundaries
#ifdef USEMPI
             call xmpi_shift(ee,':1')
             call xmpi_shift(ee,':n')
             call xmpi_shift(rr,':1')
             call xmpi_shift(rr,':n')
             call xmpi_shift(km,':1')
             call xmpi_shift(km,':n')
             call xmpi_shift(sigm,':1')
             call xmpi_shift(sigm,':n')
#endif
             !
             ! Compute mean wave direction
             !
             thetamean(i,:)=(sum(ee(i,:,:)*thet(i,:,:),2)/size(ee(i,:,:),2)) &
                  /(max(sum(ee(i,:,:),2),0.000010d0) /size(ee(i,:,:),2))
             !
             ! Energy integrated over wave directions,Hrms
             !
             E(i,:)=sum(ee(i,:,:),2)*dtheta
             R(i,:)=sum(rr(i,:,:),2)*dtheta
             DR(i,:)=sum(drr(i,:,:),2)*dtheta
             H(i,:)=sqrt(E(i,:)/par%rhog8)
             Herr=maxval(abs(Hprev-H(i,:)))
#ifdef USEMPI
             call xmpi_allreduce(Herr,MPI_MAX)	 
#endif
             ! Stopping criteria
             if (iter<par%maxiter) then
                if (Herr<par%maxerror) then
                   stopiterate=.true.
                   if(xmaster) call writelog('ls','(a,i4,a,i4)','Wave propagation row ',it,', iteration ',iter)
                endif
             else
                stopiterate=.true.
                if(xmaster) call writelog('ls','(a,i4,a,i4,a,f5.4)','Wave propagation row ',it,', iteration ',iter,', error: ',Herr)
             endif
          enddo ! End while loop
       enddo ! End do i=2:nx loop  


#ifdef USEMPI
    enddo      ! End loop through mpi grids		 
#endif

    ee(nx+1,:,:) = ee(nx,:,:)
    rr(nx+1,:,:) = rr(nx,:,:)
    E(nx+1,:)    = E(nx,:)
    R(nx+1,:)    = R(nx,:)
    DR(nx+1,:)   = DR(nx,:)
    H(nx+1,:)    = H(nx,:)
    km(nx+1,:)   = km(nx,:)
    sigm(nx+1,:) = sigm(nx,:)
    cg(nx+1,:)   = cg(nx,:)
    c(nx+1,:)    = c(nx,:)
    thet(nx+1,:,:) = thet(nx,:,:)


    k=km
    !
    ! Radiation stresses and forcing terms
    !
n=cg/c
Sxx=(n*sum((1.d0+costh**2.d0)*ee,3)-.5d0*sum(ee,3))*dtheta
Syy=(n*sum((1.d0+sinth**2.d0)*ee,3)-.5d0*sum(ee,3))*dtheta
Sxy=n*sum(sinth*costh*ee,3)*dtheta

! add roller contribution

Sxx = Sxx + sum((costh**2)*rr,3)*dtheta
Syy = Syy + sum((costh**2)*rr,3)*dtheta
Sxy = Sxy + sum(sinth*costh*rr,3)*dtheta

    do j=2,ny 
       do i=1,nx
          Fx(i,j)=-(Sxx(i+1,j)-Sxx(i,j))/dsu(i,j)             &
                  -(Sxy(i,j+1)+Sxy(i+1,j+1)- Sxy(i,j-1)-Sxy(i+1,j-1))/ &
                   (dnv(i,j-1)+dnv(i,j)+dnv(i+1,j-1)+dnv(i+1,j))
       enddo
   enddo

   do j=1,ny 
      do i=2,nx
         Fy(i,j)=-(Syy(i,j+1)-Syy(i,j))/dnv(i,j)                                  &
                 -(Sxy(i+1,j)+Sxy(i+1,j+1)-Sxy(i-1,j)-Sxy(i-1,j+1))/ &
                  (dsu(i-1,j)+dsu(i-1,j+1)+dsu(i,j)+dsu(i,j+1))
      enddo
   enddo

    Fx(:,1)=Fx(:,2)
    Fy(:,1)=Fy(:,2)
    Fx(:,ny+1)=Fx(:,ny)
    Fy(:,ny+1)=Fy(:,ny)
    Fx(1,:)=Fx(2,:)
    Fy(1,:)=Fy(2,:)

    ! wwvv
#ifdef USEMPI
    call xmpi_shift(Fx,':1')
    call xmpi_shift(Fy,':1')
    call xmpi_shift(Fx,':n')
    call xmpi_shift(Fy,':n')
    call xmpi_shift(Fx,'1:')
    call xmpi_shift(Fy,'1:')
    call xmpi_shift(Fx,'m:')
    call xmpi_shift(Fy,'m:')
#endif
    urms=par%px*H/par%Trep/(sqrt(2.d0)*sinh(k*max(hh,par%delta*H)))
    !ust=E*k/sigm/par%rho/max(hh,0.001)

    ! wave induced mass flux
    ustw=E*k/sigm/par%rho/max(hh,.001d0)
    uwf = ustw*dcos(thetamean)
    vwf = ustw*dsin(thetamean)
    ! roller contribution
    ustr=2.*R*k/sigm/par%rho/max(hh,.001d0)
    ! introduce breaker delay
    if (par%breakerdelay == 1) call breakerdelay(par,s)
    !ust = usd
    ust=usd+ustw
    !lateral boundaries
    ust(1,:) = ust(2,:)
    ust(:,1) = ust(:,2)
    ust(:,ny+1) = ust(:,ny)
    ! wwvv
#ifdef USEMPI
    call xmpi_shift(ust,'1:')
    call xmpi_shift(ust,':1')
    call xmpi_shift(ust,':n')
#endif
    D=2*par%g*par%beta*sum(rr/sqrt(cx**2+cy**2),3)*dtheta

  end subroutine wave_stationary
end module wave_stationary_module
