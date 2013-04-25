!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Copyright (C) 2007 UNESCO-IHE, WL|Delft Hydraulics and Delft University !
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
module groundwaterflow
   private    ! set private default
   ! set these public to users of groundwater module
   public gw_init
   public gw_bc
   public gwflow
   
contains



subroutine gw_init(s,par)
   use params
   use xmpi_module
   use spaceparams
   use readkey_module

   IMPLICIT NONE

   type(parameters)                            :: par
   type(spacepars)                             :: s

   character(256)                              :: fname
   real*8                                      :: aquiferbot,temp
   integer                                     :: i,j


   allocate (s%gwhead(s%nx+1,s%ny+1))
   allocate (s%gwheadb(s%nx+1,s%ny+1))
   allocate (s%gwlevel(s%nx+1,s%ny+1))
   allocate (s%gwheight(s%nx+1,s%ny+1))
   allocate (s%gwu(s%nx+1,s%ny+1))
   allocate (s%gwv(s%nx+1,s%ny+1))
   allocate (s%gwqx(s%nx+1,s%ny+1))
   allocate (s%gwqy(s%nx+1,s%ny+1))
   allocate (s%gww(s%nx+1,s%ny+1))
   allocate (s%infil(s%nx+1,s%ny+1))
   allocate (s%gwbottom(s%nx+1,s%ny+1))
   allocate (s%dinfil(s%nx+1,s%ny+1))
   allocate (s%gw0back(2,s%ny+1))
   allocate (s%gwcurv(s%nx+1,s%ny+1))

   s%gww=0.d0
   s%infil = 0.d0

   if (par%gwflow==1) then
      if (par%aquiferbotfile==' ') then     ! Not a filename
         s%gwbottom=par%aquiferbot
      else
         open(31,file=trim(par%aquiferbotfile))
         do j=1,s%ny+1
            read(31,*)(s%gwbottom(i,j),i=1,s%nx+1)
         end do
         close(31)
      endif

      if (par%gw0file==' ') then     ! Not a filename
         s%gwhead=par%gw0
      else
         open(31,file=trim(par%gw0file))
         do j=1,s%ny+1
            read(31,*)(s%gwhead(i,j),i=1,s%nx+1)
         end do
         close(31)
      endif
      where(s%wetz==1)
         s%gwhead = s%zs0
      endwhere
      
      if (xmpi_istop) then
         s%gwbottom(1,:) = s%gwbottom(2,:)
      endif
      if (xmpi_isbot) then
         s%gwbottom(s%nx+1,:) = s%gwbottom(s%nx,:)
      endif
      if (xmpi_isleft .and. s%ny>0) then
         s%gwbottom(:,1) = s%gwbottom(:,2)
      endif
      if (xmpi_isright  .and. s%ny>0) then
         s%gwbottom(:,s%ny+1) = s%gwbottom(:,s%ny)
      endif
      
      s%gwbottom = min(s%gwbottom,s%zb)
      
      s%gwhead=max(s%gwhead,s%gwbottom) 
      s%gw0back=s%gwhead(s%nx:s%nx+1,:)
      
      s%gwlevel=min(s%zb,s%gwhead)
      s%gwlevel=max(s%gwlevel,s%gwbottom) !+par%eps)
            
      s%gwheight=s%gwlevel-s%gwbottom
      
      s%gwu=0.d0
      s%gwv=0.d0
      s%gww=0.d0
      s%gwqx = 0.d0
      s%gwqy = 0.d0
      s%gwcurv=0.d0
      s%infil = 0.d0
      s%dinfil=max(par%dwetlayer/3.d0,0.02)   ! Centroid of area influenced instantly by free surface level lies at dwetlayer/3
   endif
end subroutine gw_init

subroutine gw_bc(s,par)
   use params
   use xmpi_module
   use spaceparams

   IMPLICIT NONE

   type(parameters)                            :: par
   type(spacepars)                             :: s


   s%gwbottom=min(s%gwbottom,s%zb) !-par%eps)

   if(xmpi_istop) then
!      s%gwhead(1,:)=s%zs0(1,:)
      s%gwhead(1,:)=s%gwhead(2,:)
      s%gwbottom(1,:) = s%gwbottom(2,:)
   endif
   if (xmpi_isbot) then
      if (par%tideloc==4 .or. (par%tideloc==2 .and. trim(par%paulrevere)=='land')) then
!          s%gwhead(s%nx+1,:)=s%zs0(s%nx+1,:)
          s%gwhead(s%nx+1,:)=s%gwhead(s%nx,:)
      else
         s%gwhead(s%nx+1,:)=s%gw0back(2,:)
      endif
      s%gwbottom(s%nx+1,:) = s%gwbottom(s%nx,:)
   endif
   if (xmpi_isleft .and. s%ny>0) then
      s%gwbottom(:,1) = s%gwbottom(:,2)
      s%gwhead(:,1)=s%gwhead(:,2)
      s%gwlevel(:,1)=s%gwlevel(:,2)
   endif
   if (xmpi_isright .and. s%ny>0) then
      s%gwbottom(:,s%ny+1) = s%gwbottom(:,s%ny)
      s%gwhead(:,s%ny+1)=s%gwhead(:,s%ny)
      s%gwlevel(:,s%ny+1)=s%gwlevel(:,s%ny)
   endif

#ifdef USEMPI
   call xmpi_shift(s%gwhead,':1')
   call xmpi_shift(s%gwhead,':n')
   call xmpi_shift(s%gwhead,'1:')
   call xmpi_shift(s%gwhead,'m:')
   call xmpi_shift(s%gwlevel,':1')
   call xmpi_shift(s%gwlevel,':n')
   call xmpi_shift(s%gwlevel,'1:')
   call xmpi_shift(s%gwlevel,'m:')
#endif
   if(xmpi_istop) then
      s%gwlevel(1,:)=min(s%gwhead(1,:),s%zb(1,:))
   endif
   if (xmpi_isbot) then
      s%gwlevel(s%nx+1,:)=min(s%gwhead(s%nx+1,:),s%zb(s%nx+1,:))
   endif

end subroutine gw_bc

subroutine gwflow(s,par)

  use params
  use xmpi_module
  use spaceparams

  IMPLICIT NONE

  type(parameters)                            :: par
  type(spacepars)                             :: s
  
  ! internal variables
  integer                                     :: i,j
  integer                                     :: count
  integer,parameter                           :: maxiter = 50
  real*8                                      :: err,errx,erry,errz
  real*8                                      :: flux,dzr,wdt
  real*8,dimension(:,:),allocatable,save      :: gwhu,gwhv,gwheadtop
  real*8,dimension(:,:),allocatable,save      :: Kx,Ky,Kz
  real*8,dimension(:,:),allocatable,save      :: Kxupd,Kyupd,Kzupd
  logical,dimension(:,:),allocatable,save     :: connected
  real*8,dimension(:,:),allocatable,save      :: fracdt
  real*8,dimension(:,:),allocatable,save      :: zsupd,ratio
  real*8,dimension(:,:),allocatable,save      :: dynpresupd
  real*8,dimension(:,:),allocatable,save      :: infiluncon,infilcon
  real*8,dimension(:,:),allocatable,save      :: infilhorgw,infilhorsw
  real*8,dimension(:,:),allocatable,save      :: gwumean
  real*8                                      :: factime,w1,w2,dft
  logical                                     :: initial,turb
  real*8,save                                 :: connectcrit
    
  ! shortcut pointers
  include 's.ind'
  include 's.inp'
  
  initial = .false.

  ! allocate variables
  if (.not. allocated(gwhu)) then
     allocate(gwhu(1:nx+1,1:ny+1))
     allocate(gwhv(1:nx+1,1:ny+1))
     allocate(gwheadtop(1:nx+1,1:ny+1))
     
     allocate(Kx(1:nx+1,1:ny+1))
     allocate(Ky(1:nx+1,1:ny+1))
     allocate(Kz(1:nx+1,1:ny+1))
     allocate(Kxupd(1:nx+1,1:ny+1))
     allocate(Kyupd(1:nx+1,1:ny+1))
     allocate(Kzupd(1:nx+1,1:ny+1))
     
     allocate(connected(1:nx+1,1:ny+1))
     allocate(fracdt(1:nx+1,1:ny+1))
     allocate(zsupd(1:nx+1,1:ny+1))
     allocate(ratio(1:nx+1,1:ny+1))
     
     allocate(dynpresupd(1:nx+1,1:ny+1))
     
     allocate(infiluncon(1:nx+1,1:ny+1))
     allocate(infilcon(1:nx+1,1:ny+1))
     allocate(infilhorgw(1:nx+1,1:ny+1))
     allocate(infilhorsw(1:nx+1,1:ny+1))
     allocate(gwumean(1:nx+1,1:ny+1))
     
     ! initialise variables
     Kx = par%kx
     Ky = par%ky
     Kz = par%kz
     
     Kxupd = Kx
     Kyupd = Ky
     Kzupd = Kz
     
     infiluncon = 0.d0
     infilcon = 0.d0 
     gwumean = 0.d0    
     
     dynpresupd = 0.d0
     
     initial = .true.
     if (par%gwnonh==0) then
        connectcrit = par%dwetlayer
     else
        connectcrit = par%eps
     endif
    
     ratio = gw_calculate_smoothwetlayer(s,par%dwetlayer,par%px)
  endif
  !
  ! Definition of horizontal groundwater velocity:
  ! Because there is not simple way of measuring flow in pores, Darcy relates to the 
  ! flow velocity through the ground in terms of surface water (not pore water). 
  ! Therefore, after the flux has been calculated the groundwater level change is 
  ! divided by the porosity to get changes in terms of pore water.
  !
  ! Definition of vertical groundwater/surface water exchange:
  ! infil is defined positive from sea to groundwater, gww is defined positive up 
  ! (from groundwater to sea!). Both infil and gww are in volumes of surface water
  ! (not pore water). gwu and gwv are also in terms of surface water volume.
  !
  ! Initialize
  gww = 0.d0
  infil = 0.d0
  infiluncon = 0.d0
  infilcon = 0.d0
  infilhorgw = 0.d0
  infilhorsw = 0.d0
  fracdt = 0.d0
  !
  ! If using dynamic pressure (nonh==1), very short wave lengths mess things
  ! up because they are not resolved properly in the nonh pressure solver.
  ! To surpress this, carry out pressure averaging over time scale ~ 1/4 Trep
  if (par%nonh==1) then
     factime = min(4*par%dt/par%Trep,1.d0)
     dynpresupd = (1-factime) * dynpresupd + factime * pres
  endif
  !
  ! Infiltration and exfiltration are handled separately for places where the 
  ! groundwater and surface water are connected (Poisson equation), and where
  ! they are disconnected (swash infiltration, seepage). "Connected" is set
  ! by the difference between the groundwater level and the bed level, and a
  ! numerical constant for stability
  where(wetz==1 .and. zb-gwlevel<connectcrit)
     connected = .true.
  elsewhere
     connected = .false.
  endwhere
  !
  ! Handle infiltration in unconnected cells. Some cells may be
  ! unconnected for only part of the timestep, therefore a fraction
  ! of timestep is returned which was required to reach connection
  if (par%gwnonh==1) then
     call gw_unconnected_infil(s,par,Kz,connected,fracdt,infiluncon)
!     where ((.not. connected) .and. (gwlevel>zb))
!        infiluncon = infiluncon - (gwlevel-zb)/par%dt*par%por
!        fracdt = 1.d0
!     endwhere
     where (gwlevel>zb)
        infiluncon = infiluncon - (gwlevel-zb)/par%dt*par%por
     endwhere
  endif
  !
  !
  ! Update groundwater level and water surface level to account for 
  ! infiltration and exfiltration effects. Note infil is avarage
  ! infiltration rate over whole timestep
  gwlevel = gwlevel+par%dt*infiluncon/par%por
  zsupd = zs-par%dt*infiluncon
  ! 
  ! Recalculate connected property, this is needed further on in the code
  where(zsupd-zb>=par%eps .and. zb-gwlevel<connectcrit)
     connected = .true.
  elsewhere
     connected = .false.
  endwhere
  !
  ! Thickness of groundwater layer
  gwheight=gwlevel-gwbottom
  !
  ! Determine intermediate aquifer depths (upwind scheme)
  call gw_calculate_interfaceheight(s,gwhu,gwhv,initial)
  !
  ! Calculate head smoothing function over dwetlayer 
  ratio = gw_calculate_smoothwetlayer(s,par%dwetlayer,par%px)
  ! 
  ! Calculate top boundary condition for head (pressure on groundwater
  ! from surface water, or atmospheric pressure when dry). Takes into
  ! account a transition layer, specified by user to smooth transition
  ! from surface water head to atmospheric head
  gwheadtop = gwCalculateHeadTop(nx,ny,zsupd,gwlevel,wetz,ratio,par%nonh,par%g,dynpresupd)
  ! 
  ! Compute groundwater head in cell centres
  ! In case of hydrostatic this is the same as the surface water head,
  ! or groundwater level. In order to ensure small errors do not explode
  ! minimum time averaging is done.
  ! In case of non-hydrostatic computations, the head is solved by the
  ! Poisson solver. In the case of a turbulent/modflow groundwater scheme
  ! the solver is forced using the estimate of the hydraulic conductivity
  ! from the previous time step. The estimate from the previous time step
  ! is also used to compute velocities (else no convervation of mass). The
  ! computation of the velocity also returns a new update of the hydraulic
  ! conductivity to be used in the next time step
  if (par%gwnonh == 0) then
     ! hydrostatic pressure assumption, with time delay required for
     ! stability
!     factime = min(par%dt/0.2d0,1.d0)
!     gwhead = (1.d0-factime)*gwhead + factime*gwheadtop
      gwhead = (gwhead + gwheadtop)/2
  else
     ! non-hydrostatic pressure, solve using Poisson solver. Kxupd is the
     ! updated hydraulic conductivity from the previous time step
     call gw_solver(par,s,gwheadtop,Kxupd,Kyupd,Kzupd,gwhu,gwhv,fracdt)
  endif
  !
  !
  ! Calculate groundwater velocities
  ! We use Kxupd input which is the same value as sent into the Poisson solver in 
  ! the previous command. In the case of the non-hydrostatic solver, we must force
  ! the computation of the velocities to use the hydraulic conductivity estimate 
  ! put into the Poisson solver. In practice this is only important is the hydraulic
  ! conductivity varies due to the turbulent/modflow scheme.
  !
  ! Input Kx is the basic initial hydaulic conductivity
  ! Input Kxupd is the updated hydraulic conductivity from the previous time step
  !             used in the Poisson solver
  ! Output Kxupd is the new updated hydraulic conductivity for the next time step  
  call gw_calculate_velocities(s,par,fracdt,s%gwu,s%gwv,s%gww,Kx,Ky,Kz,Kxupd,Kyupd,Kzupd)
  !
  !
  ! Calculate horizontal fluxes
  gwqx=gwu*gwhu
  gwqy=gwv*gwhv
  !
  !
  ! Stop cells from drying up
  if (s%ny>0) then
     do i=2,nx
        do j=2,ny
           if (gwlevel(i,j)<=gwbottom(i,j)+par%eps) then
              gwqx(i,j) = min(gwqx(i,j),0.d0) ! let no water out, only in
              gwqx(i-1,j) = max(gwqx(i-1,j),0.d0) ! let no water out, only in
              gwqy(i,j) = min(gwqy(i,j),0.d0) ! let no water out, only in
              gwqy(i,j-1) = max(gwqy(i,j-1),0.d0) ! let no water out, only in
           endif
        enddo
     enddo 
  else
     do i=2,nx
        if (gwlevel(i,1)<=gwbottom(i,1)+par%eps) then
           gwqx(i,1) = min(gwqx(i,1),0.d0) ! let no water out, only in
           gwqx(i-1,1) = max(gwqx(i-1,1),0.d0) ! let no water out, only in
        endif
     enddo  
  endif
  !
  !
  ! Force continuity, where gww is computed for hydrostatic computation and 
  ! gww strictly enforced (no roundoff errors) for nonhydrostatic computation
  if (ny>2) then
     do j=2,s%ny
        do i=2,nx
           gww(i,j) = (-gwqx(i,j)*dnu(i,j)+gwqx(i-1,j)*dnu(i-1,j)   &
                       -gwqy(i,j)*dsv(i,j)+gwqy(i,j-1)*dsv(i,j-1))*dsdnzi(i,j)
        enddo
     enddo
  else
     if (ny>0) then
        do i=2,nx
           gww(i,2) = (-gwqx(i,2)+gwqx(i-1,2))/dsz(i,2)
        enddo
        gww(:,1) = gww(:,2)
        gww(:,3) = gww(:,2)
     else
        do i=2,nx
           gww(i,1) = (-gwqx(i,1)+gwqx(i-1,1))/dsz(i,1)
        enddo
     endif
  endif
  !
  !
  ! updating groundwater level is complicated: 
  ! - in sections that are unconnected : dz/dt = -w/por == grad(qx,qy)/por
  ! - in sections that are connected:
  !           if w > 0 (exfiltration) : dz/dt = 0 (rise until gwlevel == zb, remainder is 
  !                                                sent to zs) 
  !           if w < 0 (infiltration) : dz/dt = 0 (unless available zs is less than w*dt,
  !                                                then drop water level on remainder)
  !
  if (par%gwnonh==1) then 
     ! non-hydrostatic computation
     do j=min(2,ny+1),max(s%ny,1)
        do i=2,s%nx
           ! exfiltration terms, same for connected and unconnected cells
           if (gww(i,j)>0.d0) then
              if (connected(i,j)) then
                 dzr = max(zb(i,j)-gwlevel(i,j),0.d0)
                 wdt = gww(i,j)*par%dt/par%por
                 if (wdt>=dzr) then
                    infilcon(i,j) = -gww(i,j)*(1.d0-dzr/wdt)
                 else
                    ! nothing
                 endif
              endif
            ! connected infiltration terms, unconnected cells already computed
            ! infiltration at the start of the groundwater flow subroutine
            elseif (gww(i,j)<0.d0) then
               if (connected(i,j)) then
                  dzr = max(zs(i,j)-zb(i,j)-par%eps,0.d0)
                  wdt = -gww(i,j)*par%dt
                  if (wdt>=dzr) then
                     infilcon(i,j) = -gww(i,j)*(dzr/wdt)
                  else
                     infilcon(i,j) = -gww(i,j)
                  endif
               endif
            else
               infilcon(i,j) = 0.d0
            endif
         enddo
      enddo
  else
     ! hydrostatic computation
     where(s%wetz==0)
        s%dinfil=par%dwetlayer/3.d0        ! Centroid of area influenced instantly by free surface level lies at dwetlayer/3
     elsewhere
        s%dinfil=min(s%dinfil,s%zb-s%gwlevel)
        s%dinfil=max(s%dinfil,par%dwetlayer/3.d0)
     endwhere
     infilcon = gw_calculate_hydrostatic_w(s,par,ratio)
     where (infilcon*par%dt>hh-par%eps)
        infilcon=(hh-par%eps)/par%dt
     endwhere
     dinfil=dinfil+infilcon*par%dt/par%por
     if (par%gwhorinfil==1) then
        call gw_horizontal_infil_exfil(s,par,infilhorgw,infilhorsw,Kx,Ky,dynpresupd)
     endif
  endif
  ! Compute new water level
  gwlevel = gwlevel + (gww+infilcon+infilhorgw)*par%dt/par%por
  !
  ! In case of hydrostatic approach, try to fix the connected cells
 
  !
  ! Update total infiltration contribution
  infil = infiluncon+infilcon+infilhorsw
  !
  ! Model boundaries
  ! Robert: check if all these are actually needed
#ifdef USEMPI
  call xmpi_shift(gwlevel,'1:')   
  call xmpi_shift(gwlevel,'m:')   
  call xmpi_shift(gwlevel,':1')   
  call xmpi_shift(gwlevel,':n')
  call xmpi_shift(gwhead,'1:')   
  call xmpi_shift(gwhead,'m:')   
  call xmpi_shift(gwhead,':1')   
  call xmpi_shift(gwhead,':n')  
  call xmpi_shift(gwcurv,'1:')   
  call xmpi_shift(gwcurv,'m:')   
  call xmpi_shift(gwcurv,':1')   
  call xmpi_shift(gwcurv,':n') 
#endif
  if (xmpi_istop) then
     gwlevel(1,:) = gwlevel(2,:)
     gwhead(1,:) = gwhead(2,:)
     gwcurv(1,:) = gwcurv(2,:)
     gwbottom(1,:) = gwbottom(2,:)
     gww(1,:) = gww(2,:)
  endif
  if (xmpi_isbot) then
     gwlevel(nx+1,:) = gwlevel(nx,:)
     gwhead(nx+1,:) = gwhead(nx,:)
     gwcurv(nx+1,:) = gwcurv(nx,:)
     gwbottom(nx+1,:) = gwbottom(nx,:)
     gww(nx+1,:) = gww(nx,:)
  endif
  if (xmpi_isleft .and. ny>0) then
     gwlevel(:,1) = gwlevel(:,2)
     gwhead(:,1) = gwhead(:,2)
     gwcurv(:,1) = gwcurv(:,2)
     gwbottom(:,1) = gwbottom(:,2)
     gww(:,1) = gww(:,2)
  endif
  if (xmpi_isright .and. ny>0) then
     gwlevel(:,ny+1) = gwlevel(:,ny)
     gwhead(:,ny+1) = gwhead(:,ny)
     gwcurv(:,ny+1) = gwcurv(:,ny)
     gwbottom(:,ny+1) = gwbottom(:,ny)
     gww(:,ny+1) = gww(:,ny)
  endif
  !
  ! output pressure head at bottom
  !
  gwheadb = gwCalculateHeadBottom(gwhead,gwcurv,gwheadtop,gwheight,nx,ny,par%gwheadmodel)
  gwheadb(1,:) = gwheadb(2,:)
  gwheadb(nx+1,:) = gwheadb(nx,:)
end subroutine gwflow

subroutine gw_unconnected_infil(s,par,Kz,connected,fracdt,infil)
  use params
  use xmpi_module
  use spaceparams

  IMPLICIT NONE
  
  type(parameters),intent(in)                 :: par
  type(spacepars)                             :: s
  real*8,dimension(s%nx+1,s%ny+1),intent(in)  :: Kz
  logical,dimension(s%nx+1,s%ny+1),intent(in) :: connected
  real*8,dimension(s%nx+1,s%ny+1),intent(out) :: fracdt,infil
  ! local
  integer                                     :: i,j
  logical                                     :: turbapprox

  ! initialise
  fracdt = 0.d0
  infil = 0.d0
  !
  ! infiltration depends on turbulent or laminar method
  select case (par%gwscheme)
     case ('laminar')
        turbapprox = .false.
     case ('turbulent')
        turbapprox = .true.
  end select
  !   
  do j=1,s%ny+1
     do i=1,s%nx+1
        if (.not. connected(i,j)) then
           if (s%wetz(i,j)==1 .and. s%gwlevel(i,j)<s%zb(i,j)) then
              call gw_calc_local_infil(par,s%zb(i,j),s%hh(i,j),s%gwlevel(i,j),s%dinfil(i,j), &
                                       s%D50top(i,j),turbapprox, &
                                       fracdt(i,j),infil(i,j))
           endif
        endif
     enddo
  enddo
  !
  ! Update infiltration layer thickness
  where (infil>0.d0)
     s%dinfil = s%dinfil+par%dt*infil/par%por
  elsewhere 
!     s%dinfil = max(s%dinfil-par%dt*par%kz/par%por,0.d0)
     s%dinfil = 0.d0
  endwhere
  ! Ensure nothing went wrong here with fracdt
  fracdt = min(max(fracdt,0.d0),1.d0)
end subroutine gw_unconnected_infil

subroutine gw_calculate_interfaceheight(s,gwhu,gwhv,initial)
  use xmpi_module
  use spaceparams

  IMPLICIT NONE
  
  type(spacepars),intent(in)                  :: s
  logical,intent(in)                          :: initial
  real*8,dimension(s%nx+1,s%ny+1),intent(out) :: gwhu,gwhv
  ! local
  integer                                     :: i,j

  ! HU
  if (initial) then
     do j=1,s%ny+1
        do i=1,s%nx
           if (s%gwhead(i+1,j)>s%gwhead(i,j)) then
              gwhu(i,j) = s%gwheight(i+1,j)
           elseif (s%gwhead(i+1,j)<s%gwhead(i,j)) then
              gwhu(i,j) = s%gwheight(i,j)
           else              
              gwhu(i,j) = 0.5d0*(s%gwheight(i+1,j)+s%gwheight(i,j))
           endif
        enddo
     enddo
  else
     do j=1,s%ny+1
        do i=1,s%nx
           if (s%gwu(i,j)>0.d0) then
              gwhu(i,j) = min(s%gwheight(i,j),max(s%zb(i+1,j)-s%gwbottom(i+1,j),0.d0))
           elseif (s%gwu(i,j)<0.d0) then
              gwhu(i,j) = min(s%gwheight(i+1,j),max(s%zb(i,j)-s%gwbottom(i,j),0.d0))
           else
              gwhu(i,j) = 0.5d0*(s%gwheight(i+1,j)+s%gwheight(i,j))
           endif
        enddo
     enddo
  endif
  gwhu = max(gwhu,0.d0) ! in case of drying cells
  ! boundaries
#ifdef USEMPI
  call xmpi_shift(gwhu,'m:')
#endif 
  if (xmpi_isbot) then
     gwhu(s%nx+1,:) = gwhu(s%nx,:)
  endif
  
  ! HV
  if (s%ny>0) then
     if (initial) then
        do j=1,s%ny
           do i=1,s%nx+1
              if (s%gwhead(i,j+1)>s%gwhead(i,j)) then
                 gwhv(i,j) = s%gwheight(i,j+1)
              elseif (s%gwhead(i,j+1)<s%gwhead(i,j)) then
                 gwhv(i,j) = s%gwheight(i,j)
              else              
                 gwhv(i,j) = 0.5d0*(s%gwheight(i,j+1)+s%gwheight(i,j))
              endif
           enddo
        enddo
     else
        do j=1,s%ny
           do i=1,s%nx+1
              if (s%gwv(i,j)>0.d0) then
                 gwhv(i,j) = min(s%gwheight(i,j),max(s%zb(i,j+1)-s%gwbottom(i,j+1),0.d0))
              elseif (s%gwv(i,j)<0.d0) then
                 gwhv(i,j) = min(s%gwheight(i,j+1),max(s%zb(i,j)-s%gwbottom(i,j),0.d0))
              else
                 gwhv(i,j) = 0.5d0*(s%gwheight(i,j+1)+s%gwheight(i,j))
              endif
           enddo
        enddo
     endif
     gwhv = max(gwhv,0.d0) ! in case of drying cells
     ! Boundaries
#ifdef USEMPI
     call xmpi_shift(gwhv,':n')
#endif  
     if (xmpi_isright) then
        gwhv(:,s%ny+1) = gwhv(:,s%ny)
     endif   
  else  ! ny==0
     gwhv(:,1) = s%gwheight(:,1)
  endif
end subroutine gw_calculate_interfaceheight

pure function gw_calculate_smoothwetlayer(s,dwetlayer,px) result(r)
  use spaceparams

  IMPLICIT NONE

  type(spacepars),intent(in)                  :: s
  real*8,intent(in)                           :: dwetlayer,px
  real*8,dimension(s%nx+1,s%ny+1)             :: r

  ! Free surface head and ratio free surface to groundwater head to be used
  r = (s%zb-s%gwlevel)/dwetlayer
  ! This ratio should lie between 0 and 1
  r = min(max(r,0.d0),1.d0)
  ! Convert this to a smoother function at the edges
  r = 1.d0-(cos(r*px)+1)/2
end function gw_calculate_smoothwetlayer

pure function gwCalculateHeadTop(nx,ny,zs,gwlevel,wetz,r,nonh,g,pres) result(gwhead)
  IMPLICIT NONE
  ! input
  integer,intent(in)                         :: nx,ny,nonh
  real*8,intent(in)                          :: g
  real*8,dimension(nx+1,ny+1),intent(in)     :: zs,gwlevel,r,pres
  integer,dimension(nx+1,ny+1),intent(in)    :: wetz
  ! result
  real*8,dimension(nx+1,ny+1)                :: gwhead
  
  ! If this point is wet (wetz==1), and 0<r<=1, additional head is added to 
  ! the groundwater, described as ratio "r" times the difference in head between
  ! the surface water and the groundwater.
  if(nonh==1) then
     gwhead = gwlevel + wetz*(zs-gwlevel+pres/g)*(1.d0-r)
  else
     gwhead = gwlevel + wetz*(zs-gwlevel)*(1.d0-r)
  endif
end function gwCalculateHeadTop

pure subroutine gwCalculateHeadGradient(nx,ny,gwhead,dsu,dnv,dheaddx,dheaddy)
   IMPLICIT NONE
   ! input/output
   integer,intent(in)                       :: nx,ny
   real*8,dimension(nx+1,ny+1),intent(in)   :: gwhead,dsu,dnv
   real*8,dimension(nx+1,ny+1),intent(out)  :: dheaddx,dheaddy
   ! internal
   integer                                  :: i,j

   dheaddx=0.d0
   dheaddy=0.d0

   do j=1,ny+1
      do i=1,nx
         dheaddx(i,j)=(gwhead(i+1,j)-gwhead(i,j))/dsu(i,j)  
      end do
   end do

   do j=1,ny
      do i=1,nx+1
         dheaddy(i,j)=(gwhead(i,j+1)-gwhead(i,j))/dnv(i,j)   
      end do
   end do
end subroutine gwCalculateHeadGradient

subroutine gw_solver(par,s,hbc,Kx,Ky,Kz,hu,hv,fracdt)
  use params
  use xmpi_module
  use spaceparams
  use solver_module, only: solver_tridiag,solver_sip

  IMPLICIT NONE

  type(parameters),intent(in)                 :: par
  type(spacepars)                             :: s
  real*8,dimension(s%nx+1,s%ny+1),intent(in)  :: hbc
  real*8,dimension(s%nx+1,s%ny+1),intent(in)  :: Kx,Ky,Kz
  real*8,dimension(s%nx+1,s%ny+1),intent(in)  :: hu,hv
  real*8,dimension(s%nx+1,s%ny+1),intent(in)  :: fracdt
  ! internal
  real*8,dimension(:,:,:),allocatable,save    :: A,work
  real*8,dimension(:,:),allocatable,save      :: rhs,x,res
  integer                                     :: i,j,n,m,k,l,it
  integer                                     :: imin,imax,jmin,jmax
  real*8                                      :: dxum,dxup,dxpc
  real*8                                      :: dyvm,dyvp,dypc
  real*8                                      :: dyum,dyup,dxvm,dxvp,dA
  real*8                                      :: hcmx,hcc,hcpx
  real*8                                      :: hcmy,hcpy
  real*8                                      :: twothird
  real*8                                      :: hum,hup,hvm,hvp
  
  twothird = 2.d0/3
  
  ! use tridiagonal solver if ny<=2
  if (s%ny<3) then
     if (par%ny==0) then
        j=1
     else
        j=2
     endif
     n = s%nx-1
     ! allocate Matrix solver coefficients
     if(.not.allocated(A)) then
        allocate(A(5,n,1))
        allocate(work(1,n,1))
        work = 0.d0
        allocate(rhs(n,1))
        allocate(x(n,1))
     endif
     ! build Matrix solver coefficients
     select case (trim(par%gwheadmodel))
        case ('parabolic')
           do k=1,n
              i = k+1
              if (s%zb(i,1)>s%gwbottom(i,1)+par%eps) then
                  ! abbreviations
                  dxum = s%dsu(i-1,j)
                  dxup = s%dsu(i,j)
                  dxpc = s%dsz(i,j)
                  hcmx = s%gwheight(i-1,j)
                  hcc = s%gwheight(i,j)
                  hcpx = s%gwheight(i+1,j)
                  hum = hu(i-1,j)
                  hup = hu(i,j)
                  if (k==-1) then
    !                 A(2,k,1) = 0.d0  ! dummy
    !                 A(1,k,1) = +Kx(i,j)*hup/dxup*twothird*hcc**2 +Kz(i-1,j)*2*dxpc*hcc*(1.d0-fracdt(i,j))
    !                 A(3,k,1) = -Kx(i,j)*hup/dxup*twothird*hcpx**2
    !                 rhs(k,1) = - Kx(i,j)*hbc(i+1,j)*hup/dxup + Kx(i+1,j)*hbc(i,j)*hup/dxup
                  elseif (k==-n) then
    !                 A(2,k,1) = -Kx(i,j)*hum/dxum*twothird*hcmx**2
    !                 A(1,k,1) = Kx(i,j)*hum/dxum*twothird*hcc**2 +Kz(i,j)*2*dxpc*hcc*(1.d0-fracdt(i,j))
    !                 A(3,k,1) = 0.d0 ! dummy
    !                 rhs(k,1) = Kx(i,j)*hbc(i,j)*hum/dxum - Kx(i,j)*hbc(i-1,j)*hum/dxum
                  else
                     A(2,k,1) = -Kx(i-1,j)*hum/dxum*twothird*hcmx**2
                     A(1,k,1) = Kx(i-1,j)*hum/dxum*twothird*hcc**2 +Kx(i,j)*hup/dxup*twothird*hcc**2 + &
                                Kz(i,j)*2*dxpc*hcc*(1.d0-fracdt(i,j))
                     A(3,k,1) = -Kx(i,j)*hup/dxup*twothird*hcpx**2
                     rhs(k,1) = Kx(i-1,j)*hum*(hbc(i  ,j)-hbc(i-1,j))/dxum & 
                              - Kx(i  ,j)*hup*(hbc(i+1,j)-hbc(i  ,j))/dxup
                  endif
              else
                 A(:,k,1) = 0.d0
                 rhs(k,1) = 0.d0
              endif
           enddo
        case ('exponential')
           do k=1,n
              i = k+1
              ! abbreviations
              dxum = s%dsu(i-1,j)
              dxup = s%dsu(i,j)
              dxpc = s%dsz(i,j)
              hcmx = s%gwheight(i-1,j)
              hcc = s%gwheight(i,j)
              hcpx = s%gwheight(i+1,j)
              hum = hu(i-1,j)
              hup = hu(i,j)
              if (k==-1) then
                 A(2,k,1) = 0.d0  ! dummy
                 A(1,k,1) = +Kx(i+1,j)*hup/dxup*twothird*hcc**2 +Kz(i,j)*2*dxpc*hcc*(1.d0-fracdt(i,j))
                 A(3,k,1) = -Kx(i+1,j)*hup/dxup*twothird*hcpx**2
                 rhs(k,1) = - Kx(i+1,j)*hbc(i+1,j)*hup/dxup + Kx(i+1,j)*hbc(i,j)*hup/dxup
              elseif (k==-n) then
                 A(2,k,1) = -Kx(i,j)*hum/dxum*twothird*hcmx**2
                 A(1,k,1) = Kx(i,j)*hum/dxum*twothird*hcc**2 +Kz(i,j)*2*dxpc*hcc*(1.d0-fracdt(i,j))
                 A(3,k,1) = 0.d0 ! dummy
                 rhs(k,1) = Kx(i,j)*hbc(i,j)*hum/dxum - Kx(i,j)*hbc(i-1,j)*hum/dxum
              else
                 A(2,k,1) = -Kx(i-1,j)*hum/dxum*(cosh(hcmx)-1/hcmx*sinh(hcmx))
                 A(1,k,1) = Kx(i-1,j)*hum/dxum*(cosh(hcc)-1/hcc*sinh(hcc)) + &
                            Kx(i,j)*hup/dxup*(cosh(hcc)-1/hcc*sinh(hcc)) + &
                            Kz(i,j)*dxpc*sinh(hcc)*(1.d0-fracdt(i,j))
                 A(3,k,1) = -Kx(i,j)*hup/dxup*(cosh(hcpx)-1/hcpx*sinh(hcpx))
                 
                 rhs(k,1) = Kx(i-1,j)*hbc(i,j)*hum/dxum - Kx(i-1,j)*hbc(i-1,j)*hum/dxum &
                          - Kx(i,j)*hbc(i+1,j)*hup/dxup + Kx(i,j)*hbc(i,j)*hup/dxup
              endif
           enddo
     end select
     call solver_tridiag(A,rhs,x,work(1,:,1),n-1,0,fixshallow=.true.) 
     ! return to global variable
     s%gwcurv(2:s%nx,j) = x(:,1)
     select case (trim(par%gwheadmodel))
        case ('parabolic')
           s%gwhead(2:s%nx,j) = hbc(2:s%nx,j) - twothird*x(:,1)*s%gwheight(2:s%nx,j)**2
        case ('exponential')
           s%gwhead(2:s%nx,j) = hbc(2:s%nx,j)+ &
                                 x(:,1)/s%gwheight(2:s%nx,j)*sinh(s%gwheight(2:s%nx,j)) &
                                 -x(:,1)*cosh(s%gwheight(2:s%nx,j))
     end select
     s%gwhead(1,j) = s%gwhead(2,j)
     s%gwhead(s%nx+1,j) = s%gwhead(s%nx,j)
     ! spread left and right
     if (s%ny==2) then
        s%gwhead(:,1) = s%gwhead(:,2)
        s%gwhead(:,3) = s%gwhead(:,2)
     endif
  else
     if (par%gwfastsolve==1) then
        ! Use tridiagonal solver per row with explicit longshore velocities
        n = s%nx+1
        m = s%ny+1
        ! allocate Matrix solver coefficients
        if(.not.allocated(A)) then
           allocate(A(5,n,1))
           allocate(work(1,n,1))
           work = 0.d0
           allocate(rhs(n,1))
           allocate(x(n,1))
        endif
        do j=2,m-1
           ! build Matrix solver coefficients
           select case (trim(par%gwheadmodel))
              case ('parabolic')
                 do i=1,n
                    imin = max(i-1,1)
                    imax = min(i+1,s%nx+1)
                    ! abbreviations
                    dxum = s%dsu(imin,j)
                    dxup = s%dsu(i,j)
                    dxpc = s%dsz(i,j)
                    dA   = s%dsdnzi(i,j)
                    dyum = s%dnu(imin,j)
                    dyup = s%dnu(i,j)
                    dyvm = s%dnv(i,j-1)
                    dyvp = s%dnv(i,j)
                    dxvm = s%dsv(i,j-1)
                    dxvp = s%dsv(i,j)
                    hcmx = s%gwheight(imin,j)
                    hcc = s%gwheight(i,j)
                    hcpx = s%gwheight(imax,j)
                    hum = hu(imin,j)
                    hup = hu(i,j)
                    hvm = hv(i,j-1)
                    hvp = hv(i,j)
                    A(2,k,1) = -Kx(imin,j)*hum*dyum/dxum*twothird*hcmx**2
                    A(1,k,1) =  Kx(imin,j)*hum*dyum/dxum*twothird*hcc**2 + &
                                Kx(i,j)*hup*dyup/dxup*twothird*hcc**2 + &
                               Kz(i,j)*2*dA*hcc*(1.d0-fracdt(i,j))
                    A(3,k,1) = -Kx(i,j)*hup*dyup/dxup*twothird*hcpx**2
                    rhs(k,1) = Kx(imin,j)*hum*dyum*(hbc(i   ,j)-hbc(imin,j))/dxum & 
                             - Kx(i   ,j)*hup*dyup*(hbc(imax,j)-hbc(i  ,j))/dxup &
                             + hvm*dxvm*s%gwv(i,j-1) - hvp*dxvp*s%gwv(i,j)

                 enddo
              case ('exponential')
                 do k=1,n
                    imin = max(i-1,1)
                    imax = min(i+1,s%nx+1)
                    ! abbreviations
                    dxum = s%dsu(imin,j)
                    dxup = s%dsu(i,j)
                    dxpc = s%dsz(i,j)
                    dA   = s%dsdnzi(i,j)
                    dyum = s%dnu(imin,j)
                    dyup = s%dnu(i,j)
                    dyvm = s%dnv(i,j-1)
                    dyvp = s%dnv(i,j)
                    dxvm = s%dsv(i,j-1)
                    dxvp = s%dsv(i,j)
                    hcmx = s%gwheight(imin,j)
                    hcc = s%gwheight(i,j)
                    hcpx = s%gwheight(imax,j)
                    hum = hu(imin,j)
                    hup = hu(i,j)
                    hvm = hv(i,j-1)
                    hvp = hv(i,j)
                    A(2,k,1) = -Kx(imin,j)*hum*dyum/dxum*(cosh(hcmx)-1/hcmx*sinh(hcmx))
                    A(1,k,1) =  Kx(imin,j)*hum*dyum/dxum*(cosh(hcc)-1/hcc*sinh(hcc)) + &
                                Kx(i,j)*hup*dyup/dxup*(cosh(hcc)-1/hcc*sinh(hcc)) + &
                                Kz(i,j)*dA*sinh(hcc)*(1.d0-fracdt(i,j))
                    A(3,k,1) = -Kx(i,j)*hup*dyup/dxup*(cosh(hcpx)-1/hcpx*sinh(hcpx))
                    
                    rhs(k,1) = Kx(imin,j)*hbc(i,j)*hum*dyum/dxum - Kx(i,j)*hbc(i-1,j)*hum*dyum/dxum &
                             - Kx(i,j)*hbc(imax,j)*hup*dyup/dxup + Kx(i,j)*hbc(i,j)*hup*dyup/dxup
                 enddo
           end select
           call solver_tridiag(A,rhs,x,work(1,:,1),n-1,0,fixshallow=.true.) 
           s%gwcurv(:,j) = x(:,1)
           select case (trim(par%gwheadmodel))
              case ('parabolic')
                 s%gwhead(2:s%nx,j) = hbc(2:s%nx,j) - twothird*x(:,1)*s%gwheight(2:s%nx,j)**2
              case ('exponential')
                 s%gwhead(2:s%nx,j) = hbc(2:s%nx,j)+ &
                                      x(:,1)/s%gwheight(2:s%nx,j)*sinh(s%gwheight(2:s%nx,j)) &
                                     -x(:,1)*cosh(s%gwheight(2:s%nx,j))
           end select
        enddo
        ! spread left and right
        s%gwhead(:,1) = s%gwhead(:,2)
        s%gwhead(:,s%ny+1) = s%gwhead(:,s%ny)
     else
        n = s%nx+1
        m = s%ny+1
        ! allocate Matrix solver coefficients
        if(.not.allocated(A)) then
           allocate(A(5,n,m))
           allocate(work(5,n,m))
           work = 0.d0
           allocate(rhs(n,m))
           allocate(x(n,m))
           x = 0.d0
           allocate(res(n,m))
           res = 0.d0
        endif
        ! build Matrix solver coefficients
        select case (trim(par%gwheadmodel))
           case ('parabolic')
              do j=1,m
                 do i=1,n
                    imin = max(i-1,1)
                    imax = min(i+1,s%nx+1)
                    jmin = max(j-1,1)
                    jmax = min(j+1,s%ny+1)
                    ! abbreviations (need to include dx and dy for mass conservation
                    ! on curvilinear grid
                    dxum = s%dsu(imin,j)
                    dxup = s%dsu(i,j)
                    dxpc = s%dsz(i,j)
                    dA   = s%dsdnzi(i,j)
                    dyum = s%dnu(imin,j)
                    dyup = s%dnu(i,j)
                    dyvm = s%dnv(i,jmin)
                    dyvp = s%dnv(i,j)
                    dxvm = s%dsv(i,jmin)
                    dxvp = s%dsv(i,j)
                    hcmx = s%gwheight(imin,j)
                    hcpx = s%gwheight(imax,j)
                    hcc = s%gwheight(i,j)
                    hcmy = s%gwheight(i,jmin)
                    hcpy = s%gwheight(i,jmax)
                    hum = hu(imin,j)
                    hup = hu(i,j)
                    hvm = hv(i,jmin)
                    hvp = hv(i,j)
                    ! Matrix A
                    ! main diagonal
                    A(1,i,j) = Kx(imin,j)*hum*dyum/dxum*twothird*hcc**2+Kx(i,j)*hup*dyup/dxup*twothird*hcc**2 + &
                               Ky(i,jmin)*hvm*dxvm/dyvm*twothird*hcc**2+Ky(i,j)*hvp*dxvp/dyvp*twothird*hcc**2 + &
                               Kz(i,j)*2*dA*hcc*(1.d0-fracdt(i,j))
                    ! down x
                    A(2,i,j) = -Kx(imin,j)*hum*dyum/dxum*twothird*hcmx**2
                    ! up x
                    A(3,i,j) = -Kx(i,j)*hup*dyup/dxup*twothird*hcpx**2 
                    ! down y
                    A(4,i,j) = -Ky(i,jmin)*hvm*dxvm/dyvm*twothird*hcmy**2
                    ! up y
                    A(5,i,j) = -Ky(i,j)*hvp*dxvp/dyvp*twothird*hcpy**2
                    ! 
                    ! RHS
                    rhs(i,j) = Kx(imin,j)*dyum*hum*(hbc(i   ,j)-hbc(imin,j))/dxum & 
                             - Kx(i   ,j)*dyup*hup*(hbc(imax,j)-hbc(i   ,j))/dxup &
                             + Ky(i,jmin)*dxvm*hvm*(hbc(i,j   )-hbc(i,jmin))/dyvm &
                             - Ky(i,j   )*dxvp*hvp*(hbc(i,jmax)-hbc(i,j   ))/dyvp
                 enddo
              enddo
              res = 0.d0
              call solver_sip(A,rhs,x,res,work,it,s%nx,s%ny) 
              s%gwcurv(:,:) = x
              s%gwcurv(1,:) = s%gwcurv(2,:)
              s%gwcurv(s%nx+1,:) = s%gwcurv(s%nx,:)
              s%gwcurv(:,1) = s%gwcurv(:,2)
              s%gwcurv(:,s%ny+1) = s%gwcurv(:,s%ny)
        end select
     endif ! speedup 2D
    
     select case (trim(par%gwheadmodel))
        case ('parabolic')
           s%gwhead = hbc - twothird*s%gwcurv*s%gwheight**2
        case ('exponential')
           s%gwhead = hbc + s%gwcurv/s%gwheight*sinh(s%gwheight) &
                          - s%gwcurv*cosh(s%gwheight)
     end select
  endif
end subroutine gw_solver

subroutine gw_horizontal_infil_exfil(s,par,infilhorgw,infilhorsw,Kx,Ky,dynpres)
  ! compute horizontal part of infiltration/exfiltration
  use params
  use xmpi_module
  use spaceparams
  
  type(parameters),intent(in)                 :: par
  type(spacepars)                             :: s
  real*8,dimension(s%nx+1,s%ny+1),intent(in)  :: Kx,Ky,dynpres
  real*8,dimension(s%nx+1,s%ny+1),intent(out) :: infilhorgw,infilhorsw
  ! internal
  real*8,dimension(s%nx+1,s%ny+1)             :: dheaddx,dheaddy
  integer                                     :: i,j
  real*8                                      :: dz,dx,dy,vel,hsurf,hgw,dhead
  real*8,parameter                            :: visc = 1.0d-6
  real*8                                      :: vcr,dum,scaleareab,scaleareas
  
  infilhorgw = 0.d0
  infilhorsw = 0.d0
  vcr = par%gwReturb/par%D50(1)*par%por*visc
  !
  ! compute horizontal flow velocity between bed and surface water VERTICAL interfaces
  ! 
  ! loop over all u-interfaces
  do i=1,s%nx
     do j=1,s%ny+1
       dz = s%zb(i+1,j)-s%zb(i,j)
       if (dz>0.d0 .and. s%wetz(i,j)==1) then ! jump up
          hsurf = s%zs(i,j)+dynpres(i,j)/par%g
          hgw = s%gwhead(i+1,j)
          dx = s%dsz(i+1,j)
          dhead = hsurf-hgw         ! positive gradient if hsurf>hgw
                                    ! negative velocitity if infiltration
          call gw_calculate_velocities_local(vel,dum,dhead/dx,Kx(i+1,j),0.d0,par%gwscheme, &
                                             .false.,.false.,par%gwheadmodel,dum,dum,vcr)
          vel = -vel                ! convert negative down velocity into positive if infiltration,
                                    ! negative if exfiltration 
          ! horizontal infiltration occurs along small side section, so scale the horizontal
          ! infiltration area by the vertical area of the cell used in the mass balance equation
          ! and infiltration velocities
          ! 
          ! Scale area for bed cell and surface water cell
          scaleareab = dz/s%dsz(i+1,j)
          scaleareas = dz/s%dsz(i,j)
          ! Limit  flow to that available in space in ground and 
          ! in surface water volume
          if (vel>0.d0) then  ! infiltration
             vel = min(vel, &
                       s%hh(i,j)*s%dsz(i,j)/par%dt/dz, &
                       max(s%zb(i+1,j)-s%gwlevel(i+1,j),0.d0)*par%por*s%dsz(i+1,j)/par%dt/dz)
          else
             vel = max(vel, &
                       -(s%gwlevel(i+1,j)-s%gwbottom(i+1,j))*s%dsz(i+1,j)*par%por/par%dt/dz)
          endif                                           
          infilhorgw(i+1,j) = infilhorgw(i+1,j)+vel*scaleareab
          infilhorsw(i,j) = infilhorsw(i,j)+vel*scaleareas
       elseif (dz<0.d0 .and. s%wetz(i+1,j)==1) then  ! jump down
          dz = -dz
          hsurf = s%zs(i+1,j)+dynpres(i+1,j)/par%g
          hgw = s%gwhead(i,j)
          dx = s%dsz(i,j)
          dhead = hsurf-hgw         ! positive gradient if hsurf>hgw
                                    ! negative velocitity if infiltration
          call gw_calculate_velocities_local(vel,dum,dhead/dx,Kx(i,j),0.d0,par%gwscheme, &
                                             .false.,.false.,par%gwheadmodel,dum,dum,vcr)
          vel = -vel                ! convert negative down velocity into positive if infiltration,
                                    ! negative if exfiltration 
          ! horizontal infiltration occurs along small side section, so scale the horizontal
          ! infiltration area by the vertical area of the cell used in the mass balance equation
          ! and infiltration velocities
          ! 
          ! Scale area for bed cell and surface water cell
          scaleareab = dz/s%dsz(i,j)
          scaleareas = dz/s%dsz(i+1,j)
          ! Limit  flow to that available in space in ground and 
          ! in surface water volume
          if (vel>0.d0) then  ! infiltration
             vel = min(vel, &
                       s%hh(i+1,j)*s%dsz(i+1,j)/par%dt/dz, &
                       max(s%zb(i,j)-s%gwlevel(i,j),0.d0)*par%por*s%dsz(i,j)/par%dt/dz)
          else
             vel = max(vel, &
                       -(s%gwlevel(i,j)-s%gwbottom(i,j))*s%dsz(i,j)*par%por/par%dt/dz)
          endif                                           
          infilhorgw(i,j) = infilhorgw(i,j)+vel*scaleareab
          infilhorsw(i+1,j) = infilhorsw(i+1,j)+vel*scaleareas
       else                ! equal bed height
          ! nothing
       endif
    enddo
 enddo
 !
 !
 ! loop over all v-interfaces
 do i=1,s%nx+1
    do j=1,s%ny
       dz = s%zb(i,j+1)-s%zb(i,j)
       if (dz>0.d0 .and. s%wetz(i,j)==1) then ! jump up
          hsurf = s%zs(i,j)+dynpres(i,j)/par%g
          hgw = s%gwhead(i,j+1)
          dy = s%dnz(i,j+1)
          dhead = hsurf-hgw         ! positive gradient if hsurf>hgw
                                    ! negative velocitity if infiltration
          call gw_calculate_velocities_local(vel,dum,dhead/dy,Kx(i,j+1),0.d0,par%gwscheme, &
                                             .false.,.false.,par%gwheadmodel,dum,dum,vcr)
          vel = -vel                ! convert negative down velocity into positive if infiltration,
                                    ! negative if exfiltration 
          ! horizontal infiltration occurs along small side section, so scale the horizontal
          ! infiltration area by the vertical area of the cell used in the mass balance equation
          ! and infiltration velocities
          ! 
          ! Scale area for bed cell and surface water cell
          scaleareab = dz/s%dnz(i,j+1)
          scaleareas = dz/s%dnz(i,j)
          ! Limit  flow to that available in space in ground and 
          ! in surface water volume
          if (vel>0.d0) then  ! infiltration
             vel = min(vel, &
                       s%hh(i,j)*s%dsz(i,j)/par%dt/dz, &
                       max(s%zb(i,j+1)-s%gwlevel(i,j+1),0.d0)*par%por*s%dnz(i,j+1)/par%dt/dz)
          else
             vel = max(vel, &
                       -(s%gwlevel(i,j+1)-s%gwbottom(i,j+1))*s%dnz(i,j+1)*par%por/par%dt/dz)
          endif                                           
          infilhorgw(i,j+1) = infilhorgw(i,j+1)+vel*scaleareab
          infilhorsw(i,j) = infilhorsw(i,j)+vel*scaleareas
       elseif (dz<0.d0 .and. s%wetz(i,j+1)==1) then  ! jump down
          dz = -dz
          hsurf = s%zs(i,j+1)+dynpres(i,j+1)/par%g
          hgw = s%gwhead(i,j)
          dy = s%dnz(i,j)
          dhead = hsurf-hgw         ! positive gradient if hsurf>hgw
                                    ! negative velocitity if infiltration
          call gw_calculate_velocities_local(vel,dum,dhead/dy,Kx(i,j),0.d0,par%gwscheme, &
                                             .false.,.false.,par%gwheadmodel,dum,dum,vcr)
          vel = -vel                ! convert negative down velocity into positive if infiltration,
                                    ! negative if exfiltration 
          ! horizontal infiltration occurs along small side section, so scale the horizontal
          ! infiltration area by the vertical area of the cell used in the mass balance equation
          ! and infiltration velocities
          ! 
          ! Scale area for bed cell and surface water cell
          scaleareab = dz/s%dnz(i,j)
          scaleareas = dz/s%dnz(i,j+1)
          ! Limit  flow to that available in space in ground and 
          ! in surface water volume
          if (vel>0.d0) then  ! infiltration
             vel = min(vel, &
                       s%hh(i,j+1)*s%dnz(i,j+1)/par%dt/dz, &
                       max(s%zb(i,j)-s%gwlevel(i,j),0.d0)*par%por*s%dnz(i,j)/par%dt/dz)
          else
             vel = max(vel, &
                       -(s%gwlevel(i,j)-s%gwbottom(i,j))*s%dnz(i,j)*par%por/par%dt/dz)
          endif                                           
          infilhorgw(i,j) = infilhorgw(i,j)+vel*scaleareab
          infilhorsw(i,j+1) = infilhorsw(i,j+1)+vel*scaleareas
       else                ! equal bed height
          ! nothing
       endif
    enddo
 enddo  

end subroutine gw_horizontal_infil_exfil

subroutine gw_calculate_velocities(s,par,fracdt,gwu,gwv,gww,Kx,Ky,Kz,Kxupd,Kyupd,Kzupd)
  use params
  use xmpi_module
  use spaceparams

  IMPLICIT NONE

  type(parameters),intent(in)                 :: par
  type(spacepars)                             :: s
  real*8,dimension(s%nx+1,s%ny+1),intent(in)  :: fracdt
  real*8,dimension(s%nx+1,s%ny+1),intent(in)  :: Kx,Ky,Kz
  real*8,dimension(s%nx+1,s%ny+1),intent(out) :: gwu,gwv,gww
  real*8,dimension(s%nx+1,s%ny+1),intent(out) :: Kxupd,Kyupd,Kzupd
  ! internal
  real*8,dimension(s%nx+1,s%ny+1)             :: dheaddx,dheaddy,dheaddz
  real*8                                      :: vcr,vest,vupd,err,fac,kest  ! for turbulent/MODFLOW approximation
  real*8,parameter                            :: visc = 1.0d-6
  integer                                     :: i,j
    
  ! Calculate the head gradient
  call gwCalculateHeadGradient(s%nx,s%ny,s%gwhead,s%dsu,s%dnv,dheaddx,dheaddy)
#ifdef USEMPI
  call xmpi_shift(dheaddy,':1')
  call xmpi_shift(dheaddx,'1:')
#endif  
  if (xmpi_istop) then
     dheaddx(1,:) = dheaddx(2,:)
  endif
  if (xmpi_isleft .and. s%ny>0) then
     dheaddy(:,1) = dheaddy(:,2)
  endif
  ! initialise
  Kxupd = Kx
  Kyupd = Ky
  Kzupd = Kz
  vcr = par%gwReturb/par%D50(1)*par%por*visc
  ! Compute of all loops
  ! u-direction
  do j=1,s%ny+1
     do i=1,s%nx+1
        call gw_calculate_velocities_local(gwu(i,j),Kxupd(i,j),dheaddx(i,j),Kx(i,j),0.d0, &
                                           par%gwscheme,.false.,par%nonh==1,par%gwheadmodel, &
                                           s%gwcurv(i,j),s%gwheight(i,j),vcr)
     enddo
  enddo
  ! v-direction
  do j=1,s%ny+1
     do i=1,s%nx+1
        call gw_calculate_velocities_local(gwv(i,j),Kyupd(i,j),dheaddy(i,j),Ky(i,j),0.d0, &
                                           par%gwscheme,.false.,par%nonh==1,par%gwheadmodel, &
                                           s%gwcurv(i,j),s%gwheight(i,j),vcr)
     enddo
  enddo
  ! w-direction
  do j=1,s%ny+1
     do i=1,s%nx+1
        call gw_calculate_velocities_local(gww(i,j),Kzupd(i,j),0.d0,Kz(i,j),fracdt(i,j), &
                                           par%gwscheme,.true.,par%nonh==1,par%gwheadmodel, &
                                           s%gwcurv(i,j),s%gwheight(i,j),vcr)
     enddo
  enddo
end subroutine gw_calculate_velocities

pure subroutine gw_calculate_velocities_local(vel,Kupd,headgrad,Kin,fracdt,gwscheme,isvert,isnonh,headmodel,gwc,gwh,vcr)
! compute local groundwater velocity component
  IMPLICIT NONE
  
  ! input/output
  real*8,intent(out)          :: vel,Kupd
  real*8,intent(in)           :: headgrad,Kin,fracdt
  character(*),intent(in)     :: gwscheme,headmodel
  logical,intent(in)          :: isvert,isnonh
  real*8,intent(in)           :: gwc,gwh,vcr
  ! internal
  real*8                      :: Re,vest,err,fac,kest,vupd
    
  select case (trim(gwscheme))
     case ('laminar') 
        if (isvert .and. isnonh) then  ! vertical flow non-hydrostatic
           select case(trim(headmodel))
              case ('parabolic')
                 vel = -Kin*2*gwc*gwh*(1.d0-fracdt)
              case ('exponential')
                 vel = -Kin*gwc*sinh(gwh)*(1.d0-fracdt)
           end select
        elseif (isvert .and. .not. isnonh) then  ! vertical flow hydrostatic
           vel = 0.d0
        else       ! horizontal flow
           vel = -Kin*headgrad
        endif
        Kupd = Kin
     case ('turbulent')
        if (isvert .and. isnonh) then  ! vertical flow non-hydrostatic
           select case(trim(headmodel))
              case ('parabolic')
                 vest = -Kin*2*gwc*gwh
                 if (abs(vest)>vcr) then
                    err = 1.d0
                    do while (err>0.001d0)
                       fac = min(sqrt(vcr/max(abs(vest),vcr)),1.d0)
                       kest = Kin*fac
                       vupd = -kest*2*gwc*gwh
                       vest = (vest+vupd)/2
                       err = abs(vupd-vest)/abs(vest)               
                    enddo
                    Kupd = kest
                 else
                    Kupd = Kin
                 endif
                 if (isnonh) then
                    vel = -Kin*2*gwc*gwh ! result bound by pressure solver estimate, will
                                         ! update next time step with better K approximation
                 else
                    vel = -Kupd*2*gwc*gwh
                 endif
              case ('exponential')
                 vest = -Kin*gwc*sinh(gwh)
                 if (abs(vest)>vcr) then
                    err = 1.d0
                    do while (err>0.001d0)
                       fac = min(sqrt(vcr/max(abs(vest),vcr)),1.d0)
                       kest = Kin*fac
                       vupd = -kest*gwc*sinh(gwh)
                       vest = (vest+vupd)/2
                       err = abs(vupd-vest)/abs(vest)               
                    enddo
                    Kupd = kest
                 else
                    Kupd = Kin
                 endif
                 if (isnonh) then
                    vel = -Kin*gwc*sinh(gwh) ! result bound by pressure solver estimate, will
                                             ! update next time step with better K approximation
                 else
                    vel = -Kupd*gwc*sinh(gwh)
                 endif
           end select  ! head model
           vel = vel*(1.d0-fracdt)
        elseif (isvert .and. .not. isnonh) then  ! vertical flow hydrostatic
           vel = 0.d0
           Kupd = Kin
        else  ! horizontal flow
            vest = -Kin*headgrad
            if (abs(vest)>vcr) then
               err = 1.d0
               do while (err>0.001d0)
                  fac = min(sqrt(vcr/max(abs(vest),vcr)),1.d0)
                  kest = Kin*fac
                  vupd = -kest*headgrad
                  vest = (vest+vupd)/2
                  err = abs(vupd-vest)/abs(vest)               
               enddo
               Kupd = kest
            else
               Kupd = Kin
            endif
            if (isnonh) then
               vel = -Kin*headgrad ! result bound by pressure solver estimate, will
                                   ! update next time step with better K approximation
            else
               vel = -Kupd*headgrad
            endif  
        endif
  end select ! turbulence model
end subroutine gw_calculate_velocities_local

function gw_calculate_hydrostatic_w(s,par,ratio) result(infil)
  use params
  use xmpi_module
  use spaceparams

  IMPLICIT NONE

  type(parameters)                            :: par
  type(spacepars)                             :: s
  real*8,dimension(s%nx+1,s%ny+1)             :: ratio
  ! local
  integer                                     :: i,j
  real*8                                      :: flux,w1,w2,w3,dummy
  logical                                     :: turbapprox
 
  ! shortcut pointers
  include 's.ind'
  include 's.inp'
  
  ! Select turbulent approximation of groundwater flow
  if (trim(par%gwscheme)=='turbulent') then
     turbapprox = .true.
  else
     turbapprox = .false.
  endif
  !
  ! Velocity in connected cells is based on part instantaneous recovery 
  ! of volume (w1) and part darcy driven flow due to top pressure (w2). 
  ! These velocities are weighted according to 'ratio', the closeness
  ! of the groundwater surface to the bed.
  if (ny==0) then
     do i=2,nx
        if (s%gwlevel(i,1)>s%zb(i,1)) then
           infil(i,1) = par%por*(s%zb(i,1)-s%gwlevel(i,1))/par%dt
!           infil(i,1) = max(infil(i,1),-par%kz)
        elseif (s%wetz(i,1)==1) then
           ! subroutine returns infiltration rate                         
           call gw_calc_local_infil(par,s%zb(i,1),s%hh(i,1),s%gwlevel(i,1),s%dinfil(i,1), &
                                    s%D50top(i,1),turbapprox,dummy,w1)
           w1 = min(w1,(1.d0+s%hh(i,1)/(par%dwetlayer/3))*par%kz)
           ! part is instant, part is through infiltration
           w2 = par%por*(s%zb(i,1)-s%gwlevel(i,1))/par%dt
           w2 = min(w2,(s%hh(i,1)-par%eps)/par%dt)
           w3 = ratio(i,1)*w1+(1.d0-ratio(i,1))*w2
           ! this should not exceed the available space for water
           infil(i,1) = min(w3,(s%zb(i,1)-s%gwlevel(i,1))*par%por/par%dt)
        else
           infil(i,1) = 0.d0
        endif
     enddo
     ! boundaries
#ifdef USEMPI
     call xmpi_shift(infil,'1:')
     call xmpi_shift(infil,'m:')
#endif
     if (xmpi_istop) then
        infil(1,1) = infil(2,1)
     endif
     if (xmpi_isbot) then
        infil(nx+1,1) = infil(nx,1)
     endif
  else ! ny>0
     do j=2,ny
        do i=2,nx
           if (s%gwlevel(i,j)>s%zb(i,j)) then
              infil(i,j) = par%por*(s%zb(i,j)-s%gwlevel(i,j))/par%dt
           elseif (s%wetz(i,j)==1) then
              ! subroutine returns infiltration rate                         
              call gw_calc_local_infil(par,s%zb(i,j),s%hh(i,j),s%gwlevel(i,j),s%dinfil(i,j), &
                                       s%D50top(i,j),turbapprox,dummy,w1)
              w1 = min(w1,(1.d0+s%hh(i,j)/(par%dwetlayer/3))*par%kz)
              ! part is instant, part is through infiltration
              w2 = par%por*(s%zb(i,j)-s%gwlevel(i,j))/par%dt
              w2 = min(w2,(s%hh(i,j)-par%eps)/par%dt)
              w3 = ratio(i,j)*w1+(1.d0-ratio(i,j))*w2
              ! this should not exceed the available space for water
              infil(i,j) = min(w3,(s%zb(i,j)-s%gwlevel(i,j))*par%por/par%dt)
           else
              infil(i,j) = 0.d0
           endif
        enddo
     enddo
#ifdef USEMPI
     call xmpi_shift(infil,'1:')
     call xmpi_shift(infil,'m:')
     call xmpi_shift(infil,':1')
     call xmpi_shift(infil,':n')
#endif
     if (xmpi_istop) then
        infil(1,:) = infil(2,:)
     endif
     if (xmpi_isbot) then
        infil(nx+1,:) = infil(nx,:)
     endif
     if (xmpi_isleft) then
        infil(:,1) = infil(:,2)
     endif
     if (xmpi_isright) then
        infil(:,ny+1) = infil(:,ny)
     endif
  endif ! ny>0
end function gw_calculate_hydrostatic_w


pure subroutine gw_calc_local_infil(par,zb,hh,gwlevel,dinfil,D50top,turb,fracdt,infil)
  use params
  
  IMPLICIT NONE

  type(parameters),intent(in)                 :: par
  real*8,intent(in)                           :: zb,hh,gwlevel,dinfil,D50top
  logical,intent(in)                          :: turb
  real*8,intent(out)                          :: fracdt,infil
  ! local
  real*8,parameter                            :: visc = 1.0d-6
  real*8                                      :: dis,dtunsat,subdt
  real*8                                      :: vest,vupd,Re,err,kze,newinfil

  ! Determine infiltration velocity in single cell
  !
  ! The routine only solves Laminar/Turbulent unsaturated flow,
  ! saturated flow handled later by Poisson solver. The fraction of 
  ! time that unsaturated flow occurs is given by fracdt.
  !
  ! Unsaturated flow:
  ! w(n+1) = kz * (1+dp/dz) = kz * (1 + (hh)/(dz(n)+w(n+1)*dt/por) )
  ! w -> (-dzold + dt/por k + Sqrt[dzold^2 + 2 dt/por dzold k + 4 dt/por hold k + dt^2 k^2])/(2 dt/por)
  !
  ! The routine for turbulent flow is the same as for Darcy flow, 
  ! but iterate until correct value of vest found
  !
  ! Distance from bed to groundwater level
  dis = max(zb-gwlevel,0.d0)
  ! how long does it take to fill this distance under free vertical flow?
  dtunsat = (dis*par%por) / (par%kz*(1+hh/dis)) 
  ! what part of timestep is unsaturated?
  fracdt = min(dtunsat,par%dt)/par%dt
  subdt = fracdt*par%dt/par%por
  ! does this need any unsaturated flow?
  if (subdt>0.d0) then
     ! estimate vertical velocity due to unsaturated flow
     vest = ( (-dinfil + subdt*par%kz + &
                     sqrt(dinfil**2 + &
                          2*subdt*dinfil*par%kz + &
                          4*subdt*hh*par%kz + &
                          subdt**2*par%kz**2&
                          )&
                     )&
                     /(2*subdt) &
                   ) 
  else
     vest = 0.d0
  endif
  ! Only set error if using turbulent flow approximation AND 
  ! the Reynolds number is above the critical Re number
  if (turb) then
     ! Calculate Reynolds number
     Re = abs(vest)*D50top/par%por/visc
     if (Re>par%gwReturb) then
        err = 1.d0
     endif
  else
     err=-1.d0
  endif
  ! update estimate if needed
  kze = par%kz     ! note that this variable MUST remain separate from Kz(i,j)
                   ! which is used to update the Poisson part of the equation                   
  do while (err>0.001d0)
     Re = abs(vest)*D50top/par%por/visc
     kze = par%kz*sqrt(par%gwReturb/Re)
     ! how long does it take to fill distance with new kz
     dtunsat = (dis*par%por) / (kze*(1+hh/dis))
     ! what part of timestep unsaturated and saturated?
     fracdt = min(dtunsat,par%dt)/par%dt
     subdt = fracdt*par%dt/par%por
     ! new vertical velocity estimate
     vupd = ( (-dinfil + subdt*kze + &
                  sqrt(dinfil**2 + &
                       2*subdt*dinfil*kze + &
                       4*subdt*hh*kze + &
                       subdt**2*kze**2&
                       )&
                  )&
                  /(2*subdt) &
                ) 
     err = abs(vest-vupd)/abs(vest)
     vest = (vest+vupd)/2
  enddo  ! err>0.001
  infil = vest*fracdt  ! this is the average over the WHOLE TIMESTEP
  ! check to see that the amount of water available is larger than the amount
  ! infiltrating
  if (infil > (hh-par%eps)/par%dt) then
     newinfil = max((hh-par%eps)/par%dt,0.d0)
     fracdt = fracdt*newinfil/infil
     infil = newinfil
  endif
end subroutine gw_calc_local_infil

pure subroutine gw_calc_local_connected_infil(par,gwhead,gwlevel,zb,headtop,hh,D50top,turb,zs,infil)
  use params
  
  IMPLICIT NONE

  type(parameters),intent(in)                 :: par
  real*8,intent(in)                           :: gwhead,gwlevel,zb,headtop,hh,D50top,zs
  logical,intent(in)                          :: turb
  real*8,intent(out)                          :: infil
  ! local
  real*8,parameter                            :: visc = 1.0d-6
  real*8                                      :: headdif
  real*8                                      :: vest,vupd,Re,err,kze
  logical                                     :: infiltrate, exfiltrate
  
  ! Determine infiltration and exfiltration velocity in single cell in case of hydrostatic
  ! groundwater computation
  !
  ! Since in connected cells hydrostatic mode gwhead == headtop, we use the difference between
  ! headtop and gwhead when the cell needs filling, and the difference between gwlevel and zb 
  ! when the cell needs emptying 
  !
  ! Saturated flow:
  ! infil = kz * (dp/dz) = kz * headdif/dz = kx * headdif / (infil*dt/por)
  ! infil = sqrt(kz * headdif * por/dt)
  !
  if(gwlevel<zb) then
     infiltrate = .true.
     headdif = min(zs-gwhead,hh)
  else
     infiltrate = .false.
     if(gwlevel>zb) then
        exfiltrate = .true.
        headdif = gwlevel-zb
     else
        exfiltrate = .false.
     endif
  endif
  
  if (infiltrate .or. exfiltrate) then
     ! compute infiltration rate
     kze = par%kx
     vest = sqrt(kze*headdif*par%por/par%dt)
     if (turb) then
        ! Calculate Reynolds number
        Re = abs(vest)*D50top/par%por/visc
        if (Re>par%gwReturb) then
           err = 1.d0
        endif
     else
        err=-1.d0
     endif
     ! Compute turbulent flow velocity
     do while (err>0.001d0)
        Re = abs(vest)*D50top/par%por/visc
        kze = par%kz*sqrt(par%gwReturb/Re)
        ! new vertical velocity estimate
        vupd = sqrt(kze*headdif*par%por/par%dt)
        err = abs(vest-vupd)/abs(vest)
        vest = (vest+vupd)/2
     enddo  ! err>0.001
  else
     vest = 0.d0
  endif
  
  if (exfiltrate) then
     vest = min(vest,(gwlevel-zb)/par%dt*par%por)
     infil = -vest
  else
     infil = min(vest,(hh-par%eps)/par%dt)
     infil = min(infil,(zb-gwlevel)/par%dt*par%por)
  endif

end subroutine gw_calc_local_connected_infil

pure function gwCalculateHeadBottom(gwhead,gwcurv,gwheadtop,gwheight,nx,ny,model) result(gwheadb)
   
   IMPLICIT NONE
   integer,intent(in)                      :: nx,ny
   real*8,dimension(nx+1,ny+1),intent(in)  :: gwhead,gwcurv,gwheadtop,gwheight
   character(*),intent(in)                 :: model
   real*8,dimension(nx+1,ny+1)             :: gwheadb
   ! internal
   integer                                 :: i,j
   real*8                                  :: m,n,h,P,A,B,H1,H2
   
   select case (trim(model))
      case ('parabolic')
         gwheadb = gwheadtop-gwcurv*gwheight**2
      case ('exponential')
         gwheadb = gwcurv+gwheadtop-gwcurv*cosh(gwheight)
   end select
end function gwCalculateHeadBottom


end module
