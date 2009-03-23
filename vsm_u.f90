subroutine vsmu(s,par,sig,fv,fd,ks,fwfac,um,vm)

use params
use spaceparams
use xmpi_module

IMPLICIT NONE

type(spacepars),target                   :: s
type(parameters)                         :: par

integer  :: i,j,ii,itmax,kmx
real*8   :: fv,fd,fwfac,ks
real*8   :: fac,omega,dhdx,dhdy,uxmean,uymean,err,htemp,tbnw,nut1,nutda,sigmas,phis,phnub1,phinub
real*8   :: sigs1,sigmat,sigma,vud,vvd
real*8   :: facA,facA1,facG,facG1,facH,facH1,facCx,facCy,facBx,facBy,facC1x,facC1y,facB1x,facB1y
real*8   :: fAS1,fBCS1x,fBCS1y,vut,vvt,facln1,facln2
real*8 , dimension(s%nx+1,s%ny+1)              :: Df,Dfx,Dfy,arms,uem,vem,um,vm 
real*8 , dimension(s%nx+1,s%ny+1)              :: uwindy,uwindx,uwind2,uwind1,taudw,tswi,tauwav 
real*8 , dimension(s%nx+1,s%ny+1)              :: z0,fw,delta,phiw,nutcst,tausx,tausy
real*8 , dimension(s%nx+1,s%ny+1)              :: nutb,nut2,nut3,nusur,sigma0,phnub2,rgh
real*8 , dimension(s%nx+1,s%ny+1)              :: logdelta,logdeltasig0 

include 's.ind'
include 's.inp'
    
itmax  = 30 
    
uwindx = 0.d0
uwindy = 0.d0

z0 = ks/33.d0

! uem and vem are the negative Stokes drift + the time averaged Langrangian flows
uem = ueu-uu+um
vem = vev-vv+vm

uwindx = cos(par%windth)*par%windv
uwindy = sin(par%windth)*par%windv
   
if (par%Trep > 0.d0) then
   omega  = 2.d0*par%px/(par%Trep+par%eps)
   arms   = urms/omega
   delta  = 0.09d0*(arms/ks)**0.82d0*ks/hh                      ! Eq 15
   fw     = min(1.39d0*(max(arms/z0,par%eps))**(-.52d0),0.3d0)  ! Eq 19
   nutb   = fw*fw*urms*urms/omega/4.d0                          ! Eq 21 (Viscosity at bottom)
else
   nutb   = 0.d0
   delta  = 0.d0
   fw     = 0.d0
endif
delta  = fd*max(delta,exp(1.d0)*z0/hh)                          ! Eq 15
delta  = min(delta,0.5d0)
phiw   = 6.d0/delta/delta                                       ! Eq A6

Df     = fwfac*0.283d0*par%rho*fw*urms**3                       ! Eq 18
Dfx    = Df*cos(thetamean)
Dfy    = Df*sin(thetamean)

! separate terms in depth-averaged viscosity

! shear stresses
uwind2 = uwindx**2 + uwindy**2
where (uwind2 > 0.0d0)
   uwind1 = dsqrt(uwind2)
   taudw =  1.9555e-05*uwind2 + 0.0024e+00*uwind1 - 0.0054e+00 + 0.0070e+00/uwind1
   tswi   = taudw*uwind1
elsewhere
   taudw  = 0.
   tswi   = 0.
endwhere

tauwav = DR/c                                 ! Eq 5

tausx = taudw*uwindx + tauwav*cos(thetamean)
tausy = taudw*uwindy + tauwav*sin(thetamean)

nut2 = hh*par%vonkar/3.d0*sqrt(tswi/par%rho)  ! Eq 10
nut3 = fv*H*(DR/par%rho)**(1.d0/3.d0)         ! Jaap: should we make this conssietent with computation of Dc and use hh instead of H?

! Viscosity at surface level

nusur  = 1.5*sqrt(nut2*nut2+nut3*nut3)        ! Eq 12 Jaap: nut1 is computed below...
nusur  = max(nusur,par%vicmol)
sigma0 = z0/hh
phnub2 = phiw*nutb                            ! part Eq A6
rgh    = par%rho*par%g*hh
logdelta     = log(1.d0/delta)
logdeltasig0 = log(delta/sigma0)
nutcst = hh*par%vonkar/6.d0
dhdy   = 0.d0
dhdx   = 0.d0

! Iteration for dhdx en dhdy on meanux and meanuy
do j=1,ny+1
   do i=1,nx+1
      do ii=1,itmax
         tbnw   = abs(rgh(i,j)*sqrt(dhdy*dhdy+dhdx*dhdx))
         nut1   = nutcst(i,j)*sqrt(tbnw/par%rho)                            ! Eq 11
         nutda  = max(sqrt(nut1**2+nut2(i,j)**2+nut3(i,j)**2),par%vicmol)   ! Eq 12
         sigmas = (nutda-nusur(i,j)/3.d0)/(nutda-nusur(i,j)/2.d0)           ! Eq A4
         phis   = nusur(i,j)/nutda/(sigmas-1.d0)                       
         phnub1 = phis*nutda
         phinub = phnub1+phnub2(i,j)                                        ! part of Eq A6
         sigs1  = (phnub1*sigmas + phnub2(i,j)*delta(i,j))/phinub
         facA   = hh(i,j)/(par%rho*phnub1)
         facA1  = hh(i,j)/(par%rho*phinub)
         facG   = facA /sigmas*(logdelta(i,j)+(sigmas-1.d0)*log((sigmas-1.d0)/(sigmas-delta(i,j))))
         facG1  = facA1/sigs1 * (logdeltasig0(i,j)+(sigs1-1.d0)*log((sigs1-delta(i,j))/(sigs1-sigma0(i,j))))
         facH   = facA* ((sigmas-1.d0)*log((sigmas-1.d0)/(sigmas-delta(i,j))) +(1.d0-delta(i,j)))
         facH1  = facA1*((sigs1-1.d0)*log((sigs1-delta(i,j))/(sigs1-sigma0(i,j))) +(delta(i,j)-sigma0(i,j)))
         facCx  = (uem(i,j)-(facG1+facG)*tausx(i,j)-(facG1-facH1/delta(i,j))*Dfx(i,j)/c(i,j))/(facH1+facH-facG1-facG)
         facCy  = (vem(i,j)-(facG1+facG)*tausy(i,j)-(facG1-facH1/delta(i,j))*Dfy(i,j)/c(i,j))/(facH1+facH-facG1-facG)
         dhdy   = facCy/rgh(i,j)
         dhdx   = facCx/rgh(i,j)
         facBx  = tausx(i,j)-facCx
         facBy  = tausy(i,j)-facCy
         facC1x = facCx-Dfx(i,j)/delta(i,j)/c(i,j)
         facC1y = facCy-Dfy(i,j)/delta(i,j)/c(i,j)
         facB1x = facBx+Dfx(i,j)/c(i,j)
         facB1y = facBy+Dfy(i,j)/c(i,j)
         uxmean = facG*facBx+facH*facCx+facG1*facB1x+facH1*facC1x
         uymean = facG*facBy+facH*facCy+facG1*facB1y+facH1*facC1y
         err = abs((uxmean-uem(i,j))/(uem(i,j)+par%eps)) + abs((uymean-vem(i,j))/(vem(i,j)+par%eps))
         if (err < 1.e-3) goto 101
    enddo ! end iteration
	stop 'no convergence'
	101 continue 	  

    ! determine velocity at top boundary layer and at e*z0

    sigmat = 2.71828*sigma0(i,j)

    facln2=log((sigs1-sigmat)/(sigs1-sigma0(i,j)))
    
    fas1   = facA1/sigs1
    fBCS1x = facB1x+facC1x*sigs1
    fBCS1y = facB1y+facC1y*sigs1

    vut = fas1 * (facB1x - fBCS1x * facln2 )
    vvt = fas1 * (facB1y - fBCS1y * facln2 )

    facln1=log(delta(i,j)/sigma0(i,j))
    facln2=log((sigs1-delta(i,j))/(sigs1-sigma0(i,j)))

    vud = fas1 * (facB1x * facln1 - fBCS1x * facln2 )
    vvd = fas1 * (facB1y * facln1 - fBCS1y * facln2 )

    ! Construct velocity profile for bottom layer (both sigma < sigmat and
    ! sigma > sigmat) and middle layer
    do ii=1,kmax
	 !JH 16-04-08 correction for velocity profile
     sigma= (1.0 + sig(i,j,ii))  !Jaap: why + sigma0? +sigma0(i,j)
	 
     if (sigma < sigmat) then
        veloc(i,j,ii,1) = (sigma/sigmat)*vut*wetu(i,j)
        veloc(i,j,ii,2) = (sigma/sigmat)*vvt*wetv(i,j)
        veloc(i,j,ii,3) = nutb(i,j)
     elseif (sigma < delta(i,j)) then
        facln1=log(sigma/sigma0(i,j))
        facln2=log((sigs1-sigma)/(sigs1-sigma0(i,j)))
        veloc(i,j,ii,1) = fAS1*(facB1x*facln1-fBCS1x*facln2)*wetu(i,j)
        veloc(i,j,ii,2) = fAS1*(facB1y*facln1-fBCS1y*facln2)*wetv(i,j)
        veloc(i,j,ii,3) = phinub*sigma*(sigs1-sigma)
     else
        facln1=log(sigma/delta(i,j))
        facln2=log((sigmas-sigma)/(sigmas-delta(i,j)))
        veloc(i,j,ii,1) = vud+facA/sigmas*(facBx*facln1-(facBx+facCx*sigmas)*facln2)*wetu(i,j)
        veloc(i,j,ii,2) = vvd+facA/sigmas*(facBy*facln1-(facBy+facCy*sigmas)*facln2)*wetv(i,j)
        veloc(i,j,ii,3) = phis*nutda*sigma*(sigmas-sigma)
     endif
   enddo

   ! define bottom layer (consistent with dwnvel.f90)
   !do ii = kmax, 1, -1
   !   htemp  = (1.0 + sig(i,j,ii))*hh(i,j)
   !   kmx = ii
   !   if (htemp>delta(i,j) .and. htemp<=delta(i-1,j)) then ! Jaap: waarom 0.05*hh?
   !      exit
   !   endif         
   !enddo
   kmx = 1

   ! uavg is depth averaged velocity
   ! uavg = sqrt(meanux*meanux+meanuy*meanuy)
   ubedm  = veloc(:,:,kmx,1)
   vbedm  = veloc(:,:,kmx,2)
   ! umod   = (ubedm*ubedm+vbedm*vbedm)**0.5d0
   ! zumod  = hh*(1.d0+sig(i,j,kmx))
  enddo
enddo


end
   
