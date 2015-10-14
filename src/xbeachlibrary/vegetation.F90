!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (C) 2011 UNESCO-IHE, WL|Delft Hydraulics and Delft University !
! Dano Roelvink, Ap van Dongeren, Ad Reniers, Jamie Lescinski,            !
! Jaap van Thiel de Vries, Robert McCall                                  !
!                                                                         !
! d.roelvink@unesco-ihe.org                                               !
! UNESCO-IHE Institute for Water Education                                !
! P.O. Box 3015                                                           !
! 2601 DA Delft                                                           !
! The Netherlands                                                         !
!                                                                         !
! This library is free software; you can redistribute it and/or           !
! modify it under the terms of the GNU Lesser General Public              !
! License as published by the Free Software Foundation; either            !
! version 2.1 of the License, or (at your option) any later version.      !
!                                                                         !
! This library is distributed in the hope that it will be useful,         !
! but WITHOUT ANY WARRANTY; without even the implied warranty of          !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        !
! Lesser General Public License for more details.                         !
!                                                                         !
! You should have received a copy of the GNU Lesser General Public        !
! License along with this library; if not, write to the Free Software     !
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307     !
! USA                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! VEGETATION MODULE XBEACH: ATTENUATION OF SHORT WAVES, IG WAVES, FLOW, AND WAVE SETUP ! ! !
! 
! Version 1.0: 
! Attenuation of short waves and IG waves 
! Jaap van Thiel de Vries, okt 2013, 
! see Linh K. Phan, Jaap S.M. van Thiel de Vries, and Marcel J.F. Stive (2015) Coastal Mangrove Squeeze in the Mekong Delta. Journal of Coastal Research: Volume 31, Issue 2: pp. 233 – 243.
! 
! Version 2.0: 
! Attenuation of short waves, IG waves and nonlinear wave effects 
! Arnold van Rooijen, okt 2015, 
! see Van Rooijen et al. submitted
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module vegetation_module
   use typesandkinds
   implicit none
   save
   private
   type vegie
   character(slen)                         :: name
   integer                                 :: npts        ! Number of points used in vertical schematization of vegetation [-]
   real*8  , dimension(:)    , allocatable :: zv          ! Vertical points used in vegetation schematization [m wrt zb_ini (zb0)]
   real*8  , dimension(:)    , allocatable :: Cd          ! Vertically integrated drag coefficient [-]
   real*8  , dimension(:)    , allocatable :: bv          ! Width of vegetation stands
   integer , dimension(:)    , allocatable :: N           ! Number of vegetation stands per unit horizontal area [m-2]
   real*8  , dimension(:,:,:), allocatable :: Dragterm1   ! Dragterm in wave action balance []
   real*8  , dimension(:,:,:), allocatable :: Dragterm2   ! Dragterm in momentum equations []
   integer , dimension(:,:)  , allocatable :: vegtype     ! spatial mapping of vegetation types [-]
   end type vegie

   type(vegie), dimension(:), allocatable    :: veg

   public vegie_init
   public vegatt

contains

subroutine vegie_init(s,par)
    use params
    use spaceparams
    use readkey_module
    use filefunctions
    use logging_module
    use interp

    IMPLICIT NONE

    type(parameters)                            :: par
    type(spacepars), target                     :: s

    character(1)                                :: ch
    integer                                     :: i,j,fid,ier
    !integer, dimension(s%nx+1,s%ny+1)           :: vegtype

    !include 's.ind'
    !include 's.inp'
    
    ! INITIALIZATION OF VEGETATION
    ! Read files with vegetation properties:
    ! file 1: list of species
    ! file 2: vegetation properties per specie (could be multiple files)
    ! file 3: distribution of species over space

    if(.not. xmaster) return ! wwvv todo: different from trunk

    if (par%vegetation == 0) then
       return
    endif

    !! 1) Check how many vegetation species are specified in vegiefile!!
    fid=create_new_fid() ! see filefunctions.F90
    call check_file_exist(par%vegiefile)
    open(fid,file=par%vegiefile)
        ier=0
        i=0
        do while (ier==0)
            read(fid,'(a)',iostat=ier)ch
            if (ier==0)i=i+1
        enddo
        par%nveg=i
        rewind(fid)

        allocate(veg(par%nveg))
        do i=1,par%nveg
            read(fid,'(a)')veg(i)%name      ! set vegetation name
        enddo
    close(fid)
    
    !! 2)  Allocate and read vegetation properties for every species!!! 
    do i=1,par%nveg  ! for each species
       call check_file_exist(veg(i)%name)
       veg(i)%npts = readkey_int(veg(i)%name,'npts',  2,        2,      10  )! Number of vertical points in vegetation schematization
       allocate (veg(i)%zv(veg(i)%npts))! Vertical coordinates (w.r.t. initial bathymetry)
       allocate (veg(i)%Cd(veg(i)%npts))
       allocate (veg(i)%bv(veg(i)%npts))
       allocate (veg(i)%N(veg(i)%npts))
       allocate (veg(i)%Dragterm1(veg(i)%npts,s%nx+1,s%ny+1))
       allocate (veg(i)%Dragterm2(veg(i)%npts,s%nx+1,s%ny+1))
       veg(i)%zv   = readkey_dblvec(veg(i)%name,'zv',veg(i)%npts,size(veg(i)%zv),0.1d0,0.05d0,20.d0)       ! Think about default values...
       veg(i)%Cd   = readkey_dblvec(veg(i)%name,'Cd',veg(i)%npts,size(veg(i)%Cd),0.1d0,0.05d0,20.d0)       ! Think about default values...
       veg(i)%bv   = readkey_dblvec(veg(i)%name,'bv',veg(i)%npts,size(veg(i)%bv),0.1d0,0.05d0,20.d0)       ! Think about default values...
       veg(i)%N    = nint(readkey_dblvec(veg(i)%name,'N',veg(i)%npts,size(veg(i)%bv),0.1d0,0.05d0,20.d0))  ! Jaap: quick&dirt transform real into integer
       
       ! If Cd is specified by user (constant), compute constant Dragterms 1 and 2 in initialization
       do j=1,veg(i)%npts ! for every vertical veg point
            if (veg(i)%Cd(j) > tiny(0.d0)) then
                veg(i)%Dragterm1(j,:,:) = 0.5d0/sqrt(par%px)*par%rho*veg(i)%Cd(j)*veg(i)%bv(j)*veg(i)%N(j) ! Drag coefficient based on first part equation 6.5 Suzuki, 2011
                veg(i)%Dragterm2(j,:,:) = 0.5d0*veg(i)%Cd(j)*veg(i)%bv(j)*veg(i)%N(j)
            endif
        enddo
    enddo
    
    !! 3)  Read spatial distribution of all vegetation species !!! 
    ! NB: vegtype = 1 corresponds to first vegetation specified in vegiefile etc.
    fid=create_new_fid() ! see filefunctions.F90
    call check_file_exist(par%vegiemapfile)
    open(fid,file=par%vegiemapfile)
    do j=1,s%ny+1    ! Is this the right way to do it in a module
       read(fid,*,iostat=ier)(s%vegtype,i=1,s%nx+1)
    enddo
    close(fid)

    ! TODO: interpolate vegetation to u-points(!)
    ! TODO: vertical profile of veg chars (linear interpolation)
    ! TODO: check 2d formulations
end subroutine vegie_init

subroutine vegatt(s,par)
    use params
    use spaceparams
    use readkey_module
    use filefunctions
    use interp

    type(parameters)                            :: par
    type(spacepars)                             :: s

    integer                                     :: i,j,ind,m
    real*8                                      :: Cdterm
    
    !include 's.ind'
    !include 's.inp'
 
    ! First compute drag coefficient (if not user-defined)
    s%Cdrag = 0d0 ! set drag coefficient to zero every timestep
    do ind=1,par%nveg  ! for each species
        do m=1,veg(ind)%npts ! for each vertical veg point
            if (veg(ind)%Cd(m)<=tiny(0.d0)) then ! If Cd is not user specified
                ! Call subroutine of M. Bendoni to compute Cd
                do j=1,s%ny+1
                    do i=1,s%nx+1
                        if (s%vegtype(i,j)>0) then ! only if vegetation is present
                            
                            ! Now fill every location with drag coef (dep on veg type)
                            call bulkdragcoeff(s,par,veg(ind)%zv(veg(ind)%npts)+s%zb0(i,j)-s%zb(i,j),ind,m,i,j,Cdterm)
                            s%Cdrag(i,j) = Cdterm
                
                            ! Compute new drag terms 1 and 2 with new Cd-values
                            veg(ind)%Dragterm1(m,i,j) = 0.5d0/sqrt(par%px)*par%rho*s%Cdrag(i,j)*veg(ind)%bv(m)*veg(ind)%N(m) ! Drag coefficient based on first part equation 6.5 Suzuki, 2011
                            veg(ind)%Dragterm2(m,i,j) = 0.5d0*s%Cdrag(i,j)*veg(ind)%bv(m)*veg(ind)%N(m)
                        else
                            veg(ind)%Dragterm1(m,i,j) = 0d0
                            veg(ind)%Dragterm2(m,i,j) = 0d0
                        endif
                    enddo
                enddo
                
            endif
        enddo
    enddo
    
    ! Attenuation by vegetation is computed in wave action balance (swvegatt) and the momentum balance (momeqveg); 
    !
    ! 1) Short wave dissipation by vegetation
    call swvegatt(s,par)
    ! 2) Mom.Eq.: Long wave dissipation, mean flow dissipation, nonlinear short wave effects, effect of emerged vegetation
    call momeqveg(s,par)
end subroutine vegatt

subroutine swvegatt(s,par)
    use params
    use spaceparams
    use readkey_module
    use filefunctions
    use interp

    type(parameters)                            :: par
    type(spacepars), target                     :: s

    integer                                     :: i,j,m,ind  ! indices of actual x,y point
    real*8                                      :: aht,hterm,htermold,Dvgt
    real*8, dimension(s%nx+1,s%ny+1)            :: Dvg,kmr

    !include 's.ind'
    !include 's.inp'

    kmr = min(max(s%k, 0.01d0), 100.d0)

    ! Set dissipation in vegetation to zero everywhere for a start
    Dvg = 0.d0
    do j=1,s%ny+1
        do i=1,s%nx+1
            ind = s%vegtype(i,j)
            htermold = 0.d0
            if (ind>0) then ! only if vegetation is present at (i,j)
                do m=1,veg(ind)%npts-1
             
                    ! Determine height of vegetation section (restricted to current bed level)
                    aht = veg(ind)%zv(m+1)+s%zb0(i,j)-s%zb(i,j)!(max(veg(ind)%zv(m)+s%zb0(i,j),s%zb(i,j)))
               
                    ! restrict vegetation height to local water depth
                    aht = min(aht,s%hh(i,j))
             
                    ! compute hterm based on ah
                    hterm = (sinh(kmr(i,j)*aht)**3+3*sinh(kmr(i,j)*aht))/(3.d0*kmr(i,j)*cosh(kmr(i,j)*s%hh(i,j))**3)
             
                    ! compute dissipation based on aht and correct for lower elevated dissipation layers (following Suzuki et al. 2012)
                    ! Take drag term = average of two vegetation points (not correct for emergent vegetation and roots below bed level)
                    Dvgt = 0.5d0*(veg(ind)%Dragterm1(m,i,j)+veg(ind)%Dragterm1(m+1,i,j))*(0.5d0*kmr(i,j)*par%g/s%sigm(i,j))**3*(hterm-htermold)*s%H(i,j)**3
             
                    ! save hterm to htermold to correct possibly in next vegetation section
                    htermold = hterm
             
                    ! add dissipation current layer
                    Dvg(i,j) = Dvg(i,j) + Dvgt
                enddo
            endif
        enddo
    enddo
    s%Dveg = Dvg

end subroutine swvegatt

subroutine momeqveg(s,par)
    use params
    use spaceparams
    use readkey_module
    use filefunctions
    use interp

    type(parameters)                            :: par
    type(spacepars)                             :: s

    integer                                     :: i,j,m,ind  ! indices of actual x,y point
    real*8                                      :: aht,ahtold,Fvgtu,Fvgtv,FvgStu,FvgStv,watr,wacr,uabsu,vabsv,Fvgnlu,Fvgnlv,FvgCau,ucan,uabsunl !uabsunl,vabsvnl,hterm,htermold,
    real*8, dimension(s%nx+1,s%ny+1)            :: Fvgu,Fvgv,kmr
    real*8, dimension(50)                       :: unl,etaw,hvegeff,Fvgnlu0

    !include 's.ind'
    !include 's.inp'

    ! Compute one force related to vegetation present in the water column:
    ! 1) Long wave velocity (ul)
    ! 2) Stokes velocity (us)
    ! 3) Non linear short wave velocity residual (ua)
    ! 4) return flow / undertow (ue)
    ! 5) wave-induced in-canopy flow (?)
    
    kmr = min(max(s%k, 0.01d0), 100.d0)

    Fvgu = 0.d0
    Fvgv = 0.d0
    Fvgnlu = 0.d0
    Fvgnlv = 0.d0
    FvgStu = 0.d0
    FvgStv = 0.d0
    ucan   = 0.d0
    FvgCau = 0.d0
    uabsunl = 0.d0

   do j=1,s%ny+1
       do i=1,s%nx+1
          ind = s%vegtype(i,j)
          ahtold = 0.d0
          if (ind>0) then ! Only if vegetation is present 
                  
            ! Compute uabsu for calculation of Fveg
            uabsu = 0.d0
            vabsv = 0.d0
            Fvgnlu0 = 0.d0
            if (par%vegnonlin == 1) then
                call swvegnonlin(s,par,i,j,unl,etaw)                
            endif
                  
            watr = 0d0
            wacr = 0d0
            do m=1,veg(ind)%npts-1   
                ! Determine height of vegetation section (restricted to current bed level)
                aht = veg(ind)%zv(m+1)+s%zb0(i,j)-s%zb(i,j)
                
                ! Determine which part of the vegetation is below the wave trough, and between trough and crest
                if (par%vegnonlin == 1) then
                    watr = minval(etaw)
                    watr = s%hh(i,j) + watr ! wave trough level
                    wacr = maxval(etaw)
                    wacr = s%hh(i,j) + wacr ! wave crest level
                else
                    watr = s%hh(i,j)
                    wacr = s%hh(i,j)
                endif
    
                if (ahtold > wacr) then ! if plant section is entirely above wave crest, then do nothing
                      
                    ! mean and long wave flow (ue)                      
                    Fvgtu = 0d0
                    Fvgtv = 0d0
                    
                    ! nonlinear waves 
                    Fvgnlu = 0.d0
                    Fvgnlv = 0.d0
                      
                else ! vegetation section is located (partly) in between wave trough and crest level                  
                    if (par%veguntow == 1) then
                        ! mean and long wave flow (ue, ve)
                        Fvgtu = max((min(aht,watr)-ahtold),0d0)/s%hh(i,j)*0.5d0*(veg(ind)%Dragterm2(m,i,j)+veg(ind)%Dragterm2(m+1,i,j))*(s%ueu(i,j)*s%vmageu(i,j))!*(s%hh(i,j)/(s%hh(i,j)-0.5d0*H(i,j)))**2
                        Fvgtv = max((min(aht,watr)-ahtold),0d0)/s%hh(i,j)*0.5d0*(veg(ind)%Dragterm2(m,i,j)+veg(ind)%Dragterm2(m+1,i,j))*(s%vev(i,j)*s%vmageu(i,j))!*(s%hh(i,j)/(s%hh(i,j)-0.5d0*H(i,j)))**2
                    else
                        ! Only long wave velocity (assume undertow is diverted over vegetation)
                        Fvgtu = max((min(aht,watr)-ahtold),0d0)/s%hh(i,j)*0.5d0*(veg(ind)%Dragterm2(m,i,j)+veg(ind)%Dragterm2(m+1,i,j))*(s%uu(i,j)*s%vmagu(i,j))!*(s%hh(i,j)/(s%hh(i,j)-0.5d0*H(i,j)))**2
                        Fvgtv = max((min(aht,watr)-ahtold),0d0)/s%hh(i,j)*0.5d0*(veg(ind)%Dragterm2(m,i,j)+veg(ind)%Dragterm2(m+1,i,j))*(s%vv(i,j)*s%vmagu(i,j))!*(s%hh(i,j)/(s%hh(i,j)-0.5d0*H(i,j)))**2
                    endif
                                                        
                    ! nonlinear waves (including emerged vegetation effect)
                    !etaw    = 0.d0
                    hvegeff = max(etaw + s%hh(i,j)-ahtold,0.d0) ! effective vegetation height over a wave cycle
                    Fvgnlu  = trapz((0.5d0*(veg(ind)%Dragterm2(m,i,j)+veg(ind)%Dragterm2(m+1,i,j))*min(hvegeff,aht)*unl*abs(unl)),par%Trep/50)/s%hh(i,j)
                    Fvgnlv  = 0.d0 !max((min(aht,watr)-ahtold),0d0)*0.5d0*(veg(ind)%Dragterm2(m,i,j)+veg(ind)%Dragterm2(m+1,i,j))*uabsunl
                    
                    ! wave induced incanopy flow (Luhar et al., 2010)
                    ucan   = sqrt(4.d0*kmr(i,j)*par%Trep*s%urms(i,j)**3/(6.d0*par%px**2))
                    FvgCau = max((min(aht,watr)-ahtold),0d0)/s%hh(i,j)*0.5d0*(veg(ind)%Dragterm2(m,i,j)+veg(ind)%Dragterm2(m+1,i,j))*ucan**2
                endif

                ! save aht to ahtold to correct possibly in next vegetation section
                ahtold = aht
                
                ! add Forcing current layer
                Fvgu(i,j) = Fvgu(i,j) + Fvgtu 

                if (par%vegnonlin == 1) then ! add nonlin wave effect
                    Fvgu(i,j) = Fvgu(i,j) + Fvgnlu
                endif
                if (par%vegcanflo == 1) then ! add in canopy flow (Luhar et al., 2010)
                    Fvgu(i,j) = Fvgu(i,j) + FvgCau
                endif
            enddo
          endif
       enddo
    enddo

    s%Fvegu = Fvgu
    s%Fvegv = Fvgv

end subroutine momeqveg

subroutine swvegnonlin(s,par,i,j,unl,etaw)
    use params
    use spaceparams
    
    IMPLICIT NONE
    
    type(parameters)                            :: par
    type(spacepars)                             :: s

    integer,intent(in)                          :: i,j
    integer                                     :: irf,jrf,ih0,it0,ih1,it1 !,m,ind,ih0,it0,ih1,it1,irf,jrf  ! indices of actual x,y point
    integer,  save                              :: nh,nt
    real*8                                      :: p,q,f0,f1,f2,f3 !,uabsunl,vabsvnl
    real*8,  save                               :: dh,dt,kmr,Urs,phi,w1,w2
    real*8, dimension(8),save                   :: urf0
    real*8, dimension(50),save                  :: urf2,urf !,urfueurfu
    real*8, dimension(50,8),save                :: cs,sn,urf1
    real*8, dimension(:,:),allocatable  ,save   :: h0,t0
    real*8, dimension(50),intent(out)           :: unl,etaw
   
    ! Subroutine to compute a net drag force due to wave skewness. Based on (matlab based) roller model with veggies by Ad.
    ! 
    ! Background:
    ! The drag force (Fveg) is a function of u*abs(u), which is zero for linear waves. For non-linear, skewed waves the 
    ! depth-averaged velocity integrated over the wave period is zero. However, due to the sharp peaks and flat troughs 
    ! the integral of u*abs(u) is non-zero, and can significantly reduce wave setup, or even lead to set-down (e.g. Dean & Bender,2006).
    !
    ! Here we use a method based on Rienecker & Fenton (1981), similar to the method used for onshore sediment transport due to wave asymmetry/
    ! skewness (see also morphevolution.F90 + Van Thiel de Vries Phd thesis par 6.2.3).
    !
        
    ! load Ad's RF-table (update for depth averaged velocities?)
    include 'RFveg.inc'
    
    !include 's.ind'
    !include 's.inp'
           
    ! Initialize/Prepare for interpolation of RF-value from RFveg-table
    if (.not. allocated(h0)) then
        allocate (h0         (s%nx+1,s%ny+1))
        allocate (t0         (s%nx+1,s%ny+1))

        dh = 0.03d0
        dt = 1.25d0
        nh = floor(0.54d0/dh);
        nt = floor(25.d0/dt);
        !construct velocity profile based on cosine/sine functions / Fourier components
        do irf=1,8
            do jrf=1,50                
                cs(jrf,irf) = cos((jrf*2*par%px/50)*irf)
                sn(jrf,irf) = sin((jrf*2*par%px/50)*irf)
            enddo
        enddo        
    endif

    h0 = min(nh*dh,max(dh,min(s%H,s%hh)/s%hh))
    t0 = min(nt*dt,max(dt,par%Trep*sqrt(par%g/s%hh)))
    
!    Initialize
    urf0     = 0.d0
    urf1     = 0.d0
    urf2     = 0.d0
    urf      = 0.d0
    w1       = 0.d0
    w2       = 0.d0
    phi      = 0.d0
    Urs      = 0.d0
    kmr      = 0.d0
    
    ! ! Now compute weight factors (w1,w2) for relative contribution of cosine and sine functions (for w1 = 1: only cosines -> 
    ! fully skewed Stokes wave, for w2 = 1: only sines -> fully asymmetric wave) based on Ruessink.
    kmr   = min(max(s%k(i,j), 0.01d0), 100.d0)
    Urs   = s%H(i,j)/kmr/kmr/(s%hh(i,j)**3)! Ursell number (check factor 3/8??, Ruessink et al 2012: Urs = 3/8*H*k/(kh)^3)
    
    ! What is the effect of vegetation on the Ursell number? Ruessink (emperical) relation only based on cases without vegetation...
    ! Need some kind of 'apparent depth'? !Ur = H(jh)/k/k/(hap^3)
       
    ! Compute phase and weight factors
    phi  = par%px/2*(1-tanh(0.815/(Urs**0.672)))! according to Ruessink et al 2012 (eq 10): p5 = 0.815 ipv 0.64; ip6 = 0.672 ipv 0.6, Dano&Ad book: 0.64 and 0.6
    w1   = 1-phi/(par%px/2)!w1 = 1.d0  if fully skewed waves
    w2   = 1.d0-w1
    ! or use relation between w1 and phi as in Phd thesis Jaap (eq 6.13)??
    
    ! Interpolate RieneckerFenton velocity from RFveg table from Ad
    ! in ftab-dimension, only read 4:11 and sum later
    
    ! interpolate RF table values....
    ih0=floor(h0(i,j)/dh)
    it0=floor(t0(i,j)/dt)
    ih1=min(ih0+1,nh)
    it1=min(it0+1,nt)
    p=(h0(i,j)-ih0*dh)/dh
    q=(t0(i,j)-it0*dt)/dt
    f0=(1-p)*(1-q)
    f1=p*(1-q)
    f2=q*(1-p)
    f3=p*q
           
    ! Compute velocity amplitude per component
    do irf=1,8
        urf0(irf) = f0*RFveg(irf+3,ih0,it0)+f1*RFveg(irf+3,ih1,it0)+ f2*RFveg(irf+3,ih0,it1)+f3*RFveg(irf+3,ih1,it1)
    enddo

    ! fill velocity amplitude matrix urf1([50 time points, 8 components])
    do irf=1,8
        urf1(:,irf) = urf0(irf)
    enddo
            
    ! Compute velocity profile matrix per component
    urf1 = urf1*(w1*cs+w2*sn)
    
    ! Add velocity components
    urf2 = sum(urf1,2)
    
    ! Scale the results to get velocity profile over wave period
    unl  = urf2*sqrt(par%g*s%hh(i,j))
    etaw = unl*sqrt(max(s%hh(i,j),0.d0)/par%g)
    
end subroutine swvegnonlin

function trapz(y,dx) result (value)
    implicit none
    real*8               :: integral,value,dx
    real*8, dimension(:) :: y
    integer              :: i,n

    integral = 0.d0
    n        = size(y)-1.d0
    do i=1,n
        integral = integral+dx*(y(i+1)+y(i))/2
    end do
    value = integral
    
end function trapz

subroutine bulkdragcoeff(s,par,ahh,ind,m,i,j,Cdterm)
!    Michele Bendoni: subroutine to calculate bulk drag coefficient for short wave
!    energy dissipation based on the Keulegan-Carpenter number
!    Ozeren et al. (2013) or Mendez and Losada (2004)

    use params
    use spaceparams

    implicit none

    type(parameters)     :: par
    type(spacepars)      :: s
    real*8,  intent(out) :: Cdterm
    real*8,  intent(in)  :: ahh    ! [m] plant (total) height
    integer, intent(in)  :: ind,m,i,j
    ! Local variables
    real*8               :: alfav  ! [-] ratio between plant height and water depth
    real*8               :: um     ! [m/s] typical velocity acting on the plant
    real*8               :: Tp     ! [s] reference wave period
    real*8               :: KC     ! [-] Keulegan-Carpenter number
    real*8               :: Q      ! [-] modified Keulegan-Carpenter number
    integer              :: myflag ! 1 => Ozeren et al. (2013); 2 => Mendez and Losada (2004)
    !
    !
    myflag = 2
    !
    ! Representative wave period
    Tp = 2*par%px/s%sigm(i,j)
    !
    ! Coefficient alfa
    if (ahh>=s%hh(i,j)) then
       alfav = 1.d0
    else
       alfav = ahh/s%hh(i,j)
    end if
    !
    ! Representative orbital velocity
    ! (Could we also use urms here?)
    um = 0.5d0*s%H(i,j)*s%sigm(i,j)*cosh(s%k(i,j)*alfav*s%hh(i,j))/sinh(s%k(i,j)*s%hh(i,j))
    !
    ! Keulegan-Carpenter number
    KC = um*Tp/veg(ind)%bv(m)
    if (KC > 0d0) then
        KC = KC
    endif
    !
    ! Bulk drag coefficient
    if (myflag == 1) then
       ! 
       ! Approach from Ozeren et al. (2013), eq?
       !
       if (KC>=10.d0) then
          Cdterm = 0.036d0+50.d0/(KC**0.926d0)
       else
          Cdterm = 0.036d0+50.d0/(10.d0**0.926d0)
       endif
    elseif (myflag == 2) then
       !
       ! Approach from Mendez and Losada (2004), eq. 40
       ! Only applicable for Laminaria Hyperborea (kelp)???
       !
       Q = KC/(alfav**0.76d0)
       if (Q>=7) then
          Cdterm = exp(-0.0138*Q)/(Q**0.3d0)
       else
          Cdterm = exp(-0.0138*7)/(7**0.3d0)
       endif
    endif
    !
    !Cdterm = 0.5d0/sqrt(par%px)*par%rho*Cdtemp*veg(ind)%bv(m)*veg(ind)%N(m)
    !
end subroutine bulkdragcoeff

end module vegetation_module
