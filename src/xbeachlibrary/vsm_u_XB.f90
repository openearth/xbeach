subroutine vsm_u_XB (  ue    ,ve    ,h     ,ks    ,               &
                       rho   ,tauwx ,tauwy ,theta ,               &
                       urms  ,omega ,k     ,diss  ,               &
                       e ,er ,fcvisc,facdel,ws    ,               &
                       fw    ,facdf ,swglm ,n     ,sigz  ,        &
                       uz    ,vz    ,ustz  ,nutz  ,               &
                       ue_sed,ve_sed )  
    
   implicit none

!  Input/output parameters      
   real*8,  intent(in)   :: ue                ! depth-averaged eulerian velocity, ksi comp
   real*8,  intent(in)   :: ve                ! depth-averaged eulerian velocity, eta comp
   real*8,  intent(in)   :: h                 ! water depth
   real*8,  intent(in)   :: ks                ! nikuradse roughness
   real*8,  intent(in)   :: rho               ! density
   real*8,  intent(in)   :: tauwx             ! wind shear stress, ksi comp
   real*8,  intent(in)   :: tauwy             ! wind shear stress, eta comp
   real*8,  intent(in)   :: theta             ! wave angle w.r.t. grid ksi-direction
   real*8,  intent(in)   :: urms              ! orbital velocity 
   real*8,  intent(in)   :: omega             ! angular peak frequency
   real*8,  intent(in)   :: k                 ! wave number
   real*8,  intent(in)   :: diss              ! breaking dissipation
   real*8,  intent(in)   :: e                 ! wave energy
   real*8,  intent(in)   :: er                ! roller energy
   real*8,  intent(in)   :: fcvisc            ! calibration factor breaking-induced viscosity (~0.05)
   real*8,  intent(in)   :: facdel            ! calibration factor bbl thickness (~20)
   real*8,  intent(in)   :: facdf             ! calibration factor bottom dissipation (~1)
   real*8,  intent(in)   :: ws                ! average fall velocity
   real*8,  intent(in)   :: fw                ! wave friction coefficient
   integer, intent(in)   :: swGLM             ! switch to use GLM formulation (1) or not (0)
   integer, intent(in)   :: n                 ! number of vertical layers
   real*8,  intent(out)  :: sigz(n)           ! vertical sigma grid (0 at bottom, 1 at surface)
   real*8,  intent(out)  :: uz(n)             ! vertically distributed velocity, ksi comp.
   real*8,  intent(out)  :: vz(n)             ! vertically distributed velocity, eta comp.
   real*8,  intent(out)  :: ustz(n)           ! vertically distributed stokes drift in wave direction
   real*8,  intent(out)  :: nutz(n)           ! vertically distributed viscosity
   real*8,  intent(out)  :: ue_sed            ! advection velocity sediment (ksi-dir)
   real*8,  intent(out)  :: ve_sed            ! advection velocity sediment (eta-dir)
!  Local variables      
   real*8     :: kappa,nut1,nut2,nut3,g,uorb
   real*8     :: nutda,nusur,nutb,numin
   real*8     :: hrms,z0,cw,aorb,delta,phiw
   real*8     :: Df,Dfx,Dfy,tauwav,Qw,Qwe,Qwx,Qwy
   real*8     :: ustokes,tausx,tausy,tswi,cf,sigma0,uold
   real*8     :: tbnw,dhdy,dhdx,sigmas,phis,phinub,sigs1
   real*8     :: facA,facA1,facG,facG1,facH,facH1
   real*8     :: facCx,facCy,facBx,facBy,facC1x,facC1y,facB1x,facB1y
   real*8     :: uxmean,uymean,uabs,err,sigmat,facln1,facln2
   real*8     :: vut,vvt,vud,vvd,Qstok,difstok,sigma,hd
   real*8     :: ustar,ast,rousepar
   real*8     :: crel(n),dsig(n)
   integer    :: order_st=2,itmax,iter,i
!-------------------------------------------------------------------
!- Modified version Dirk-Jan Walstra (08-06-2007 --> copied from V205)
!- It is assumed that the mass flux has a vertical distribution
!- as predicted by 2nd or 3rd order Stokes drift.
!- The bottom shear stress is assumed to react only to the Stokes
!- drift at the bottom. This shear stress is used to construct the 
!- undertow profile etc.
!- At the end the stokes drift is subtracted from the computed profile.
!- The roller contribution to the mass flux is included in a depth 
!- averaged manner.
!-------------------------------------------------------------------
   
   itmax  = 30 
   numin  = 1.d-6
   kappa  = 0.41
   g      = 9.81
   hrms   = sqrt(8.*e/rho/g)
   z0     = ks/33.
   uorb   =urms*sqrt(2.d0)
   cw=omega/k
   if (omega.gt.0.) then
      aorb   = uorb/omega
      delta  = 0.09*(aorb/ks)**0.82*ks/h
   !   fw     = min(1.39 * (uorb/(omega*z0))**(-.52),0.3)
   endif

   delta  = facdel*max(delta,exp(1.)*z0/h)
   delta  = min(delta,0.5)
   phiw   = 6./delta/delta
       
   Df    = facDf*.283*rho*fw*uorb**3
   Dfx   = Df*cos(theta)
   Dfy   = Df*sin(theta)

!  Separate terms in depth-averaged viscosity

   if (cw.gt.0.) then
      tauwav = diss/cw
      if (SWGLM.gt.0) then
         Qwe = (e+2.0*er)/cw/rho
!        Determine stokes velocity at the bottom
         call stokes(hrms,omega,k,h,0.d0,ustokes,swglm)
         Qw = ustokes*h 
      else
         Qw = (e+2.0*er)/cw/rho
      endif
   else
      tauwav=0.
      Qw =0.
      Qwe=0.
   endif

   tausx = tauwx + tauwav*cos(theta)
   tausy = tauwy + tauwav*sin(theta)
   tswi  = sqrt(tauwx**2+tauwy**2)
   Qwx   = Qw*cos(theta)
   Qwy   = Qw*sin(theta)

   nut2 = h*kappa/3.*sqrt(tswi/rho)
   nut3 = fcvisc*hrms*(Diss/rho)**(1./3.)

!  Viscosity at surface level

   nusur  = 1.5*sqrt(nut2*nut2+nut3*nut3)
   nusur  = max(nusur,numin)

!  Viscosity at bottom

   cf     = fw/2.
   nutb   = cf*cf*uorb*uorb/omega
 
   sigma0=z0/h
   uold=0.

! Iteration for u en dhdx
  
   do iter=1,itmax

      tbnw = abs(rho*g*h*sqrt(dhdy*dhdy+dhdx*dhdx))
     
      nut1 = h*kappa/6.*sqrt(tbnw/rho)
     
!     Total depth-averaged viscosity
     
      nutda = sqrt(nut1*nut1+nut2*nut2+nut3*nut3)
      nutda=max(nutda,numin)
     
      sigmas = (nutda-nusur/3.)/(nutda-nusur/2.)
      phis   = nusur/nutda/(sigmas-1.)
     
      facA   = h/(rho*phis*nutda)
      phinub =  phis*nutda        + phiw*nutb
      sigs1  = (phis*nutda*sigmas + phiw*nutb*delta)/ &
               (phis*nutda        + phiw*nutb      )
      facA1  = h/(rho*phinub)
     
      facG = facA/sigmas* (                                   &
              log(1./delta) +                                 &
              (sigmas-1.)*log((sigmas-1.)/(sigmas-delta))  )  
     
      facG1= facA1/sigs1* (                                   &
              log(delta/sigma0) +                             &
              (sigs1-1.)*log((sigs1-delta)/(sigs1-sigma0))  )
     
      facH = facA* (                                          &
              (sigmas-1.)*log((sigmas-1.)/(sigmas-delta)) +   &
              (1.-delta)   )
     
      facH1= facA1* (                                         &
              (sigs1-1.)*log((sigs1-delta)/(sigs1-sigma0)) +  &
              (delta-sigma0)   )
     
!     facCx = (-Qwx/h-(facG1+facG)*tausx-(facG1-facH1/delta)*Dfx/cw)/
      facCx = (  ue  -(facG1+facG)*tausx-(facG1-facH1/delta)*Dfx/cw)/ &
              (facH1+facH-facG1-facG)
      dhdx  = facCx/(rho*g*h)
!     facCy = rho*g*h*dhdy
      facCy = (  ve  -(facG1+facG)*tausy-(facG1-facH1/delta)*Dfy/cw)/ &
              (facH1+facH-facG1-facG)
      dhdy  = facCy/(rho*g*h)
     
      facBx = tausx-facCx
      facBy = tausy-facCy
     
      facC1x= facCx-Dfx/delta/cw
      facC1y= facCy-Dfy/delta/cw
     
      facB1x= facBx+Dfx/cw
      facB1y= facBy+Dfy/cw
     
      uxmean = facG*facBx+facH*facCx+facG1*facB1x+facH1*facC1x
      uymean = facG*facBy+facH*facCy+facG1*facB1y+facH1*facC1y
      uabs=sqrt(max(uxmean*uxmean+uymean*uymean,1.d-6))
     
      err  = abs((uold-uabs)/uabs)
     
      if (err.lt.1.e-3) exit
      uold=uabs
     
    enddo
    
    if (err.gt.1.e-3) then
       stop 'no convergence in vsm_u'
    endif

!   end iteration

!   generate vertical grid

    hd=1.d0
    do i=1,n
        sigz(i)=0.01d0*(hd/0.01d0)**(float(i-1)/float(n-1))/hd
    enddo

    
!   determine velocity at top boundary layer and at e*z0

      sigmat = 2.71828*sigma0

      facln1=log(sigmat/sigma0)
      facln2=log((sigs1-sigmat)/(sigs1-sigma0))

      vut = facA1/sigs1 * (                       &
            facB1x              * facln1 -        &
           (facB1x+facC1x*sigs1)* facln2 )
      vvt = facA1/sigs1 * (                       &
            facB1y              * facln1 -        &
           (facB1y+facC1y*sigs1)* facln2 )

      facln1=log(delta/sigma0)
      facln2=log((sigs1-delta)/(sigs1-sigma0))

      vud = facA1/sigs1 * (                       &
            facB1x              * facln1 -        &
           (facB1x+facC1x*sigs1)* facln2 )
      vvd = facA1/sigs1 * (                       &
            facB1y              * facln1 -        &
           (facB1y+facC1y*sigs1)* facln2 )

   if (SWGLM.gt.0) then
!- Construct vertical profile of 2nd or 3rd order stokes drift 
     do i=1,n
        call stokes(hrms,omega,k,h,sigz(i)*h,ustz(i),order_st)
     enddo
!- Determine mass flux of stokes drift
     Qstok=.5*sigz(1)*h*ustz(1)
     do i=2,n
       Qstok=Qstok+.5*h*((sigz(i)-sigz(i-1))*(ustz(i)+ustz(i-1)))
     enddo
!- Determine difference between stokes mass flux and total mass flux
!- including roller contribution
     difstok= Qwe - Qstok
!-   difstok=0.
   endif

!  Construct velocity profile for bottom layer (both sigma < sigmat and
!  sigma > sigmat) and middle layer

   do 10 i=1,n

      sigma=sigz(i)

      if (sigma .lt. sigmat) then

       uz(i) = (sigma/sigmat)*vut
       vz(i) = (sigma/sigmat)*vvt
       nutz(i) = nutb

      elseif (sigma .lt. delta) then

       facln1=log(sigma/sigma0)
       facln2=log((sigs1-sigma)/(sigs1-sigma0))

       uz(i)  = facA1/sigs1 * (                      &
                facB1x               * facln1 -      &
               (facB1x+facC1x*sigs1) * facln2 )
       vz(i)  = facA1/sigs1 * (                      &
                facB1y               * facln1 -      &
               (facB1y+facC1y*sigs1) * facln2 )      
       nutz(i) = phinub*sigma*(sigs1-sigma)

      else

       facln1=log(sigma/delta)
       facln2=log((sigmas-sigma)/(sigmas-delta))

       uz(i)  = vud + facA/sigmas * (                &
               facBx              * facln1 -         &
              (facBx+facCx*sigmas)* facln2 )
       vz(i)  = vvd + facA/sigmas * (                &
               facBy              * facln1 -         &
              (facBy+facCy*sigmas)* facln2 )
       nutz(i) = phis*nutda*sigma* ( sigmas - sigma )

      endif
!-   construct eulerian velocity profile (stokes corrected and roller 
!-   contribution included as a depth averaged flux)
!-   Ueuler=Uglm-Ustokes (i.e. Veloc(4,i)=uz(i)-Qroller/h+Ustokesbottom)
!-  NOTE that uz(i) now is the GLM velocities which is also used for 
!-  Suspended Transport Computations !!!
     if (SWGLM.gt.0) then
        uz(i) = uz(i) - (ustz(i)+difstok/h)*cos(theta) + Qwx/h
        vz(i) = vz(i) - (ustz(i)+difstok/h)*sin(theta) + Qwy/h
     endif
     

10 continue

     ! Construct representative velocity for sediment transport using Rouse profile
     ustar=sqrt(0.5d0*fw)*max(urms,1.e-6)
     rousepar=ws/kappa/ustar
     ast=ks/h
     do i=1,n
        if (sigz(i)>ast) then
           crel(i)=(sigz(i)/ast*(ast-1.d0)/(sigz(i)-1.d0))**(-rousepar)
        else
           crel(i)=1.d0
        endif
     enddo
     dsig(1)=0.5d0*sigz(2)
     do i=2,n-1
         dsig(i)=(sigz(i+1)-sigz(i-1))*0.5d0
     enddo
     dsig(n)=0.5d0*(sigz(n)-sigz(n-1))
     if (crel(1)>1.d-6) then
        ue_sed=sum(dsig*uz*crel)/sum(dsig*crel)
        ve_sed=sum(dsig*vz*crel)/sum(dsig*crel) 
     else
        ue_sed=uz(1)
        ve_sed=vz(1)
     endif
     
end subroutine vsm_u_XB
   
subroutine stokes (hrms,omega,k,h,z,ustokes,order_stk)
   real*8              :: hrms,omega,k,h,z,ustokes
   integer             :: order_stk
   
   real*8              :: a,f3,s
   a=hrms/2.
   if (omega.gt.0.) then
     if (order_stk.eq.2) then
!--------------------------------------------------------------------
!--Determine Second order Stokes drift
!--------------------------------------------------------------------
       ustokes=omega*k*a**2*(cosh(2*k*z)/(2.*sinh(k*h)*sinh(k*h)))
     else
!--------------------------------------------------------------------
!--Determine Third order Stokes drift
!--------------------------------------------------------------------
       f3 = 3./16.*(8.*(cosh(k*h))**6. + 1)/(sinh(k*h)**6.)
       s  = 108.*(hrms/(.5*k*k*f3)) + 12.*sqrt(12./((.25*k*k*f3)**3.) &
            + 81.*((hrms/(.5*k*k*f3))**2.))
       a  = s**(1./3.)/6. - 2.*1./(.25*k*k*f3*s**(1./3.))             
       ustokes=.5*(omega/k)*cosh(2*k*z)*(k*a/sinh(k*h))**2*           &
               (1.-(k*a/sinh(k*h))*cosh(k*z))
     endif
   else
     ustokes=0.
   endif
end subroutine stokes

