module morphevolution

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
contains
  subroutine transus(s,par)
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
    use params
    use spaceparams
    use xmpi_module
    use interp
    ! use vsmumod

    IMPLICIT NONE

    type(spacepars),target                   :: s
    type(parameters)                         :: par

    integer                                  :: i,isig
    integer                                  :: j,jg
    real*8                                   :: exp_ero

    real*8,dimension(:),allocatable,save     :: chain,cumchain
    real*8,dimension(:,:),allocatable,save   :: vmag2,uau,uav,um,vm
    real*8,dimension(:,:),allocatable,save   :: ccvt,dcdz,dsigt,aref
    real*8,dimension(:,:),allocatable,save   :: cc,ccb,cu,cv,Sus,Svs
    real*8,dimension(:,:),allocatable,save   :: cub,cvb,Sub,Svb,pbbedu,pbbedv
    real*8,dimension(:,:),allocatable,save   :: suq3d,svq3d,eswmax,eswbed,sigs,deltas
    real*8,dimension(:,:,:),allocatable,save :: dsig,ccv,sdif,cuq3d,cvq3d,fac
    
    real*8,dimension(:,:),allocatable,save   :: sinthm,costhm

    include 's.ind'
    include 's.inp'

    if (.not. allocated(vmag2)) then
       allocate(vmag2 (nx+1,ny+1))
       allocate(uau (nx+1,ny+1))
       allocate(uav (nx+1,ny+1))
       allocate(cu  (nx+1,ny+1))
       allocate(cv  (nx+1,ny+1))
       allocate(cc  (nx+1,ny+1))
       allocate(ccb (nx+1,ny+1))
       allocate(fac (nx+1,ny+1,par%ngd))
       allocate(Sus (nx+1,ny+1))
       allocate(Svs (nx+1,ny+1))
       allocate(cub (nx+1,ny+1))
       allocate(cvb (nx+1,ny+1))
       allocate(Sub (nx+1,ny+1))
       allocate(Svb (nx+1,ny+1))
       allocate(pbbedu (nx+1,ny+1))
       allocate(pbbedv (nx+1,ny+1))
       allocate(ccvt (nx+1,ny+1))
       allocate(dcdz (nx+1,ny+1))
       allocate(dsigt (nx+1,ny+1))
       allocate(dsig(s%nx+1,s%ny+1,par%kmax))
       allocate(ccv(s%nx+1,s%ny+1,par%kmax))
       allocate(sdif(s%nx+1,s%ny+1,par%kmax))
       allocate(um (nx+1,ny+1))
       allocate(vm (nx+1,ny+1))
       allocate(deltas(nx+1,ny+1))
       allocate(sigs(nx+1,ny+1))
       allocate(eswmax(nx+1,ny+1))
       allocate(eswbed(nx+1,ny+1))
       allocate(suq3d(nx+1,ny+1))
       allocate(svq3d(nx+1,ny+1))
       allocate(cuq3d(nx+1,ny+1,par%kmax))
       allocate(cvq3d(nx+1,ny+1,par%kmax))
       allocate(aref(nx+1,ny+1))
       allocate(chain(par%kmax))
       allocate(cumchain(par%kmax))
       allocate (sinthm(nx+1,ny+1))
       allocate (costhm(nx+1,ny+1))
       uau       = 0.d0
       uav       = 0.d0
       um        = 0.d0
       vm        = 0.d0
       chain     = 0.0d0
       cumchain  = 0.0d0
       fac       = 1.d0
       exp_ero   = 0.d0
       ! generate sigma grid shape...
       do isig=2,par%kmax
          chain(isig) = par%sigfac**(isig-1)
          cumchain(isig) = cumchain(isig-1)+chain(isig)
       enddo
    endif
    
    ! use eulerian velocities
    vmag2     = ue**2+ve**2
    cu        = 0.0d0
    cv        = 0.0d0
    cub       = 0.0d0
    cvb       = 0.0d0
    dcsdx     = 0.0d0
    dcsdy     = 0.0d0

    sinthm = sin(thetamean-alfaz)
    costhm = cos(thetamean-alfaz)

    ! compute long wave turbulence due to breaking
    if (par%lwt==1) then
       call longwaveturb(s,par)
    endif

    if (par%swave==1) then
       ! include wave skewness and assymetry in sediment advection velocity
       if (trim(par%waveform)=='ruessink_vanrijn')then
          call RvR(s,par)
       elseif (trim(par%waveform)=='vanthiel') then
          call vT(s,par)
       endif

    endif

    ! calculate equilibrium concentration/sediment source
    if (trim(par%form)=='soulsby_vanrijn') then           ! Soulsby van Rijn
       call sb_vr(s,par)
    elseif (trim(par%form)=='vanthiel_vanrijn') then       ! Van Thiel de Vries & Reniers 2008
       call sednew(s,par)
    end if

    ! compute reduction factor for sediment sources due to presence of hard layers
    do jg = 1,par%ngd
       do j=1,ny+1
          do i=1,nx+1
             exp_ero = par%morfac*par%dt/(1.d0-par%por)*hh(i,j)*(ceqsg(i,j,jg)*pbbed(i,j,1,jg)/Tsg(i,j,jg) &
                                                               + ceqbg(i,j,jg)*pbbed(i,j,1,jg)/par%dt) 
             fac(i,j,jg) = min(1.d0,structdepth(i,j)*pbbed(i,j,1,jg)/max(tiny(0.d0),exp_ero) )         ! limit erosion to available sediment on top of hard layer
             !if (fac(i,j,jg)*exp_ero > dzbed(i,j,1)*pbbed(i,j,1,jg)) then
             !   fac(i,j,jg) = min(fac(i,j,jg),dzbed(i,j,1)*pbbed(i,j,1,jg)/max(tiny(0.d0),exp_ero) )  ! limit erosion to available sand in top layer
             !   write(*,*)'WARNING: expected erosion from top layer is larger than available sand in top layer'
             !endif
          enddo
       enddo
    enddo

    ! compute diffusion coefficient
    Dc = par%nuh+par%nuhfac*hh*(DR/par%rho)**(1.d0/3.d0)

    do jg = 1,par%ngd
       cc = ccg(:,:,jg)
       !ccb = ccbg(:,:,jg)
       if (D50(jg)>0.002d0) then
          ! RJ: set ceqsg to zero for gravel.
          ! Dano: try without this fix cc = 0.d0 ! Can be used to test total transport mode
       endif
       !
       ! X-direction
       !
       ! Get ua in u points and split out in u and v direction
       uau(1:nx,:) = 0.5*(ua(1:nx,:)*costhm(1:nx,:)+ua(2:nx+1,:)*costhm(1:nx,:))
       uav(1:nx,:) = 0.5*(ua(1:nx,:)*sinthm(1:nx,:)+ua(2:nx+1,:)*sinthm(1:nx,:))
       ! Compute vmagu including ua
       vmagu = sqrt((uu+uau)**2+(vu+uav)**2)
       ! sediment advection velocity for suspended load and bed load respectively
       ! REMARK: when vreps does not equal vv; no mass conservation 
       ureps = ueu+uau
       urepb = ueu+uau  ! RJ maybe reduce this velocity?
       !
       do j=1,ny+1
          do i=1,nx
             if(ureps(i,j)>0.d0) then
                ! test cu(i,j)=cc(i,j)
                cu(i,j)=par%thetanum*cc(i,j)+(1.d0-par%thetanum)*cc(min(i+1,nx),j)
                cub(i,j)=par%thetanum*pbbed(i,j,1,jg)*ceqbg(i,j,jg)+(1.d0-par%thetanum)&
                     *pbbed(min(i+1,nx),j,1,jg)*ceqbg(min(i+1,nx),j,jg)
                !cub(i,j)=par%thetanum*ccb(i,j)+(1.d0-par%thetanum)*ccb(min(i+1,nx),j)
             elseif(ureps(i,j)<0.d0) then
                cu(i,j)=par%thetanum*cc(i+1,j)+(1.d0-par%thetanum)*cc(max(i,2),j)
                cub(i,j)=par%thetanum*pbbed(i+1,j,1,jg)*ceqbg(i+1,j,jg)+(1.d0-par%thetanum)&
                     *pbbed(max(i,2),j,1,jg)*ceqbg(max(i,2),j,jg)
                !cub(i,j)=par%thetanum*ccb(i+1,j)+(1.d0-par%thetanum)*ccb(max(i,2),j)
             else
                cu(i,j)=0.5d0*(cc(i,j)+cc(i+1,j))
                cub(i,j)=0.5d0*(pbbed(i,j,1,jg)*ceqbg(i,j,jg)+pbbed(i+1,j,1,jg)*ceqbg(i+1,j,jg))
                !cub(i,j)=0.5d0*(ccb(i,j)+ccb(i+1,j))
             endif
             dcsdx(i,j)=(cc(i+1,j)-cc(i,j))/dsu(i,j)

          enddo
       enddo
       ! wwvv dcdx(nx:1,:) is still untouched, correct this ofr the parallel case
#ifdef USEMPI
       call xmpi_shift(dcsdx,'m:')

#endif
       cu(nx+1,:) = cc(nx+1,:) !Robert
       ! wwvv fix this in parallel case
#ifdef USEMPI
       call xmpi_shift(cu,'m:')
#endif
       ! wwvv in parallel version, there will be a discrepancy between the values
       ! of dzbdx(nx+1,:).
       !wwvv so fix that
#ifdef USEMPI
       call xmpi_shift(dzbdx,'m:')
#endif
       !
       Sus = 0.d0
       Sub = 0.d0

       ! suspended load
       Sus=par%sus*(cu*ureps*hu-Dc*hu*dcsdx-par%facsl*cu*vmagu*hu*dzbdx)*wetu   !No bed slope term in suspended transport?
       ! bed load
       Sub=par%bed*(cub*urepb*hu-par%facsl*cub*vmagu*hu*dzbdx)*wetu 
       ! 
       do j=1,ny+1
          do i=1,nx
             if(Sub(i,j)>0.d0) then
                pbbedu(i,j) = pbbed(i,j,1,jg)
             elseif(Sub(i,j)<0.d0) then
                pbbedu(i,j)= pbbed(i+1,j,1,jg)
             else
                pbbedu(i,j)=0.5d0*(pbbed(i,j,1,jg)+pbbed(i+1,j,1,jg))
             endif
          enddo
       enddo
       !
       Sub = pbbedu*Sub
       !
       ! Y-direction
       !
       ! Jaap: get ua in v points and split out in u and v direction
       if (ny>0) then
          uau(:,1:ny) = 0.5*(ua(:,1:ny)*costhm(:,1:ny)+ua(:,2:ny+1)*costhm(:,1:ny))
          uav(:,1:ny) = 0.5*(ua(:,1:ny)*sinthm(:,1:ny)+ua(:,2:ny+1)*sinthm(:,1:ny))
       else
          uau=ua*costhm
          uav=ua*sinthm
       endif
       ! Jaap: compute vmagv including ua
       vmagv = sqrt((uv+uau)**2+(vv+uav)**2)
       ! sediment advection velocity for suspended load and bed load respectively
       ! REMARK: when vreps does not equal vv; no mass conservation 
       vreps = vev+uav   
       vrepb = vev+uav   ! RJ maybe reduce this velocity?
       !
       if (ny>0) then
          do j=1,ny
             do i=1,nx+1
                if(vreps(i,j)>0) then
                   cv(i,j)=par%thetanum*cc(i,j)+(1.d0-par%thetanum)*cc(i,min(j+1,ny))
                   cvb(i,j)=par%thetanum*pbbed(i,j,1,jg)*ceqbg(i,j,jg)+(1.d0-par%thetanum)&
                        *pbbed(i,min(j+1,ny),1,jg)*ceqbg(i,min(j+1,ny),jg)
                   !cvb(i,j)=par%thetanum*ccb(i,j)+(1.d0-par%thetanum)*ccb(i,min(j+1,ny))
                elseif(vreps(i,j)<0) then
                   cv(i,j)=par%thetanum*cc(i,j+1)+(1.d0-par%thetanum)*cc(i,max(j,2))
                   cvb(i,j)=par%thetanum*pbbed(i,j+1,1,jg)*ceqbg(i,j+1,jg)+(1.d0-par%thetanum)&
                        *pbbed(i,max(j,2),1,jg)*ceqbg(i,max(j,2),jg)
                   !cvb(i,j)=par%thetanum*ccb(i,j+1)+(1.d0-par%thetanum)*ccb(i,max(j,2))
                else
                   cv(i,j)=0.5d0*(cc(i,j)+cc(i,j+1)) !Jaap: cc instead of cv
                   cvb(i,j)=0.5d0*(pbbed(i,j,1,jg)*ceqbg(i,j,jg)+pbbed(i,j+1,1,jg)*ceqbg(i,j+1,jg))
                   !cvb(i,j)=0.5d0*(ccb(i,j)+ccb(i,j+1))
                end if
                dcsdy(i,j)=(cc(i,j+1)-cc(i,j))/dnv(i,j) !Jaap

             end do
          end do
          ! wwvv dcdy(:,ny+1) is not filled in, so in parallel case:
#ifdef USEMPI
          call xmpi_shift(dcsdy,':n')

#endif
          cv(:,ny+1) = cc(:,ny+1) !Robert
          ! wwvv in parallel version, there will be a discrepancy between the values
          ! of dzbdy(:,ny+1).
          ! wwvv so fix that
#ifdef USEMPI
          call xmpi_shift(dzbdy,':n')
#endif
       else
          cv = cc
          cvb = ceqbg(:,:,jg)
       endif ! ny>0
       !
       ! Compute sedimnent transport in v-direction
       !
       Svs = 0.d0
       Svb = 0.d0
       ! Suspended load
       Svs=par%sus*(cv*vreps*hv-Dc*hv*dcsdy-par%facsl*cv*vmagv*hv*dzbdy)*wetv
       ! Bed load
       Svb=par%bed*(cvb*vrepb*hv-par%facsl*cvb*vmagv*hv*dzbdy)*wetv
       !
       do j=1,ny
          do i=1,nx+1
             if(Svb(i,j)>0) then
                pbbedv(i,j)=pbbed(i,j,1,jg)
             else if(Svb(i,j)<0) then
                pbbedv(i,j)=pbbed(i,j+1,1,jg)
             else
                pbbedv(i,j)=0.5d0*(pbbed(i,j,1,jg)+pbbed(i,j+1,1,jg))
             end if
          end do
       end do
       !
       Svb = pbbedv*Svb
       !
       ! BRJ: implicit concentration update (compute sources first, sink must be computed after updating actual sed.conc.)
       ! 
       if (ny>0) then
          do j=2,ny
             do i=2,nx
                ero(i,j,jg) = fac(i,j,jg)*hh(i,j)*ceqsg(i,j,jg)*pbbed(i,j,1,jg)/Tsg(i,j,jg)    ! Changed to hh from hold by RJ (13072009) !**2/max(hh(i,j),par%hmin)
                ! depo_ex(i,j,jg) = max(hold(i,j),0.01d0)*cc(i,j)/Tsg(i,j,jg)                    
                ! BRJ: the volume in the water column is updated and not the volume concentration.
                cc(i,j) = (par%dt*Tsg(i,j,jg))/(par%dt+Tsg(i,j,jg))* &
                     (hold(i,j)*cc(i,j)/par%dt -((Sus(i,j)*dnu(i,j)-Sus(i-1,j)*dnu(i-1,j)+&
                                                  Svs(i,j)*dsv(i,j)-Svs(i,j-1)*dsv(i,j-1))*dsdnzi(i,j)-&
                                                  ero(i,j,jg)))

                cc(i,j)=max(cc(i,j),0.0d0) ! Jaap: negative cc's are possible...   
                cc(i,j)=min(cc(i,j),par%cmax*hh(i,j))
                depo_ex(i,j,jg) = cc(i,j)/Tsg(i,j,jg)  
             enddo
          enddo
       else
          j=1
          do i=2,nx
             ero(i,j,jg) = fac(i,j,jg)*hh(i,j)*ceqsg(i,j,jg)*pbbed(i,j,1,jg)/Tsg(i,j,jg)    ! Changed to hh from hold by RJ (13072009) !**2/max(hh(i,j),par%hmin)
             ! depo_ex(i,j,jg) = max(hold(i,j),0.01d0)*cc(i,j)/Tsg(i,j,jg)                    
             ! BRJ: the volume in the water column is updated and not the volume concentration.
             cc(i,j) = (par%dt*Tsg(i,j,jg))/(par%dt+Tsg(i,j,jg))* &
                       (hold(i,j)*cc(i,j)/par%dt -((Sus(i,j)*dnu(i,j)-Sus(i-1,j)*dnu(i-1,j))*dsdnzi(i,j)-&
                        ero(i,j,jg)))

             cc(i,j)=max(cc(i,j),0.0d0) ! Jaap: negative cc's are possible...   
             cc(i,j)=min(cc(i,j),par%cmax*hh(i,j))

             depo_ex(i,j,jg) = cc(i,j)/Tsg(i,j,jg)  
          enddo
       endif

       ! Jaap: bed load --> Tsg = par%dt for bed load....
       !do j=2,ny
       !  do i=2,nx
       !     ero(i,j,jg) = ero(i,j,jg) + fac(i,j,jg)*hh(i,j)*ceqbg(i,j,jg)*pbbed(i,j,1,jg)/par%dt                 
       !     ccb(i,j) = (par%dt/2.d0)* &
       !                              (hold(i,j)*ccb(i,j)/par%dt -((Sub(i,j)-Sub(i-1,j))/(xu(i)-xu(i-1))+&
       !                                                           (Svb(i,j)-Svb(i,j-1))/(yv(j)-yv(j-1))-&
       !							                               fac(i,j,jg)*hh(i,j)*ceqbg(i,j,jg)*pbbed(i,j,1,jg)/par%dt) )	
       !
       !     ccb(i,j)=max(ccb(i,j),0.0d0) ! Jaap: negative cc's are possible...
       !	 depo_ex(i,j,jg) = depo_ex(i,j,jg) + ccb(i,j)/par%dt
       !  enddo
       !enddo

       cc = cc/hh
       !ccb = ccb/hh

       ! do lateral bounadries...
       if(xmpi_istop)then
          cc(1,:)=cc(2,:)
          !ccb(1,:)=ccb(2,:)
          ero(1,:,jg)=ero(2,:,jg)
          depo_ex(1,:,jg)=depo_ex(2,:,jg)
       endif
       if(xmpi_isleft .and. ny>0)then
          cc(:,1)=cc(:,2)
          !ccb(:,1)=ccb(:,2)
          ero(:,1,jg)=ero(:,2,jg)
          depo_ex(:,1,jg)=depo_ex(:,2,jg)
       endif
       if(xmpi_istop)then
          cc(nx+1,:)=cc(nx+1-1,:)
          !ccb(nx+1,:)=ccb(nx+1-1,:)
          ero(nx+1,:,jg)=ero(nx,:,jg)
          depo_ex(nx+1,:,jg)=depo_ex(nx,:,jg)
       endif
       if(xmpi_isright .and. ny>0)then
          cc(:,ny+1)=cc(:,ny+1-1)
          !ccb(:,ny+1)=ccb(:,ny+1-1)
          ero(:,ny+1,jg)=ero(:,ny,jg)
          depo_ex(:,ny+1,jg)=depo_ex(:,ny,jg)
       endif

       ! wwvv fix the first and last rows and columns of cc in parallel case
#ifdef USEMPI
       call xmpi_shift(cc,'1:')
       call xmpi_shift(cc,'m:')
       call xmpi_shift(cc,':1')   ! Maybe not necessary as zb calculated from 1:ny+1, but just in case...
       call xmpi_shift(cc,':n')   ! Dito
       !call xmpi_shift(ccb,'1:')
       !call xmpi_shift(ccb,'m:')
       !call xmpi_shift(ccb,':1')   ! Maybe not necessary as zb calculated from 1:ny+1, but just in case...
       !call xmpi_shift(ccb,':n')   ! Dito
#endif
       ! Jaap
       ! cc=cc*wetz
       !
       ! wwvv border columns and rows of ccg Svg and Sug have to be communicated
       ccg(:,:,jg) = cc
       !ccbg(:,:,jg) = ccb
       Svsg(:,:,jg) = Svs
       Susg(:,:,jg) = Sus
       Svbg(:,:,jg) = Svb
       Subg(:,:,jg) = Sub

    end do ! number of sediment fractions

    vmag=sqrt(max(vmag2,par%umin))

  end subroutine transus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine bed_update(s,par)
    use params
    use interp
    use spaceparams
    use xmpi_module

    IMPLICIT NONE

    type(spacepars),target              :: s
    type(parameters)                    :: par

    integer                                     :: i,j,j1,jg,ii,ie,id,je,jd,jdz,ndz,indx,di
    integer , dimension(s%nx+1)                 :: slopeind,bermind
    integer , dimension(:,:,:),allocatable,save :: indSus,indSub,indSvs,indSvb
    real*8                                      :: dzb,dzmax,dzt,dzleft,sdz,dzavt,fac,Savailable,dxfac,dyfac
    real*8                                      :: strucslope,bermwidth,rb,gamB,irrb,runup_old,runup_max,first
    real*8 , dimension(:,:),allocatable,save    :: dzbtot,Sout,hav
    real*8 , dimension(par%ngd)                 :: edg,edg1,edg2,dzg
    real*8 , dimension(:),pointer               :: dz
    real*8 , dimension(:,:),pointer             :: pb
    logical                                     :: aval

    include 's.ind'
    include 's.inp'

    if (.not. allocated(dzbtot)) then 
       allocate(dzbtot(s%nx+1,s%ny+1))
       allocate(Sout(s%nx+1,s%ny+1))
       allocate(hav(s%nx+1,s%ny+1))
       allocate(indSus(s%nx+1,s%ny+1,par%ngd))
       allocate(indSub(s%nx+1,s%ny+1,par%ngd))
       allocate(indSvs(s%nx+1,s%ny+1,par%ngd))
       allocate(indSvb(s%nx+1,s%ny+1,par%ngd))
    endif

    ! Super fast 1D
    if (ny==0) then
       j1 = 1
    else
       j1 = 2
    endif
    dzbtot = 0.d0
    dzbdt  = 0.d0
    dzb    = 0.d0

    if (par%t>=par%morstart .and. par%morfac > .999d0) then
       !
       ! bed_predict
       !
       ! reduce sediment transports when hard layer comes to surface
       ! this step is mainly necessary at the transition from hard layers to sand
       if (par%struct == 1) then

          do jg = 1,par%ngd
             indSus = 0
             indSub = 0
             indSvs = 0
             indSvb = 0
             Sout   = 0.d0 
             do j=j1,ny+1
                do i=2,nx+1
                   ! fluxes at i,j
                   if (Subg(i,j,jg) > 0.d0) then      ! bed load u-direction
                      indSub(i,j,jg) = 1
                      Sout(i,j) = Sout(i,j) + Subg(i,j,jg)*dnu(i,j)
                   endif
                   ! fluxes at i-1,j
                   if (Subg(i-1,j,jg) < 0.d0 ) then   ! bed load u-direction
                      Sout(i,j) = Sout(i,j) - Subg(i-1,j,jg)*dnu(i-1,j)
                   endif
                   if (par%sourcesink==0) then  
                      ! fluxes at i,j
                      if (Susg(i,j,jg) > 0.d0 ) then     ! suspended load u-direction
                         indSus(i,j,jg) = 1
                         Sout(i,j) = Sout(i,j) + Susg(i,j,jg)*dnu(i,j)
                      endif
                      ! fluxes at i-1,j
                      if (Susg(i-1,j,jg) < 0.d0 ) then   ! suspended load u-direction
                         Sout(i,j) = Sout(i,j) - Susg(i-1,j,jg)*dnu(i-1,j)
                      endif
                   endif
                enddo
             enddo
             if (ny>0) then
                do j=j1,ny+1
                   do i=2,nx+1
                      if (Svbg(i,j,jg) > 0.d0 ) then     ! bed load v-direction
                         indSvb(i,j,jg) = 1
                         Sout(i,j) = Sout(i,j) + Svbg(i,j,jg)*dsv(i,j)
                      endif
                      ! fluxes at i,j-1
                      if (Svbg(i,j-1,jg) < 0.d0 ) then   ! bed load v-direction
                         Sout(i,j) = Sout(i,j) - Svbg(i,j-1,jg)*dsv(i,j-1)
                      endif
                      if (par%sourcesink==0) then
                         if (Svsg(i,j,jg) > 0.d0 ) then     ! suspended load v-direction
                            indSvs(i,j,jg) = 1
                            Sout(i,j) = Sout(i,j) + Svsg(i,j,jg)*dsv(i,j)
                         endif
                         ! fluxes at i,j-1
                         if (Svsg(i,j-1,jg) < 0.d0 ) then   ! suspended load v-direction
                            Sout(i,j) = Sout(i,j) - Svsg(i,j-1,jg)*dsv(i,j-1)
                         endif
                      endif ! sourcesink = 0
                   enddo !nx+1
                enddo !ny+1
             endif !ny>0
             !
             do j=j1,ny+1
                do i=2,nx+1
                   Savailable = structdepth(i,j)*pbbed(i,j,1,jg)/par%morfac/par%dt*(1.d0-par%por)/dsdnzi(i,j)
                   ! reduction factor for cell outgoing sediment transports
                   fac  = min(1.d0,Savailable/max(Sout(i,j),tiny(0.d0)) )
                   ! fix sediment transports for the presence of a hard layer; remind indSus etc are 1 in cases of cell outgoing transports
                   ! updated S         oell outgoing transports                  cell incoming transports                   
                   if (fac < 1.d0)then
                      Subg(i,j,jg)   = fac*indSub(i,j,jg)*Subg(i,j,jg)         + (1-indSub(i,j,jg))*Subg(i,j,jg) 
                      Subg(i-1,j,jg) = fac*(1-indSub(i-1,j,jg))*Subg(i-1,j,jg) + indSub(i-1,j,jg)*Subg(i-1,j,jg)
                      if (ny>0) then
                         Svbg(i,j,jg)   = fac*indSvb(i,j,jg)*Svbg(i,j,jg)         + (1-indSvb(i,j,jg))*Svbg(i,j,jg)
                         Svbg(i,j-1,jg) = fac*(1-indSvb(i,j-1,jg))*Svbg(i,j-1,jg) + indSvb(i,j-1,jg)*Svbg(i,j-1,jg)
                      endif
                      if (par%sourcesink==0) then
                         Susg(i,j,jg)   = fac*indSus(i,j,jg)*Susg(i,j,jg)         + (1-indSus(i,j,jg))*Susg(i,j,jg) 
                         Susg(i-1,j,jg) = fac*(1-indSus(i-1,j,jg))*Susg(i-1,j,jg) + indSus(i-1,j,jg)*Susg(i-1,j,jg)
                         if (ny>0) then
                            Svsg(i,j,jg)   = fac*indSvs(i,j,jg)*Svsg(i,j,jg)         + (1-indSvs(i,j,jg))*Svsg(i,j,jg)
                            Svsg(i,j-1,jg) = fac*(1-indSvs(i,j-1,jg))*Svsg(i,j-1,jg) + indSvs(i,j-1,jg)*Svsg(i,j-1,jg)
                         endif !ny = 0
                      endif ! sourcesink = 0
                   endif !fac<1.d0 
                enddo ! nx+1
             enddo !ny + 1
          enddo !par%ngd
       endif !struct == 1

       if (ny>0) then
          do j=2,ny
             do i=2,nx

                ! bed level changes per fraction in this morphological time step in meters sand including pores
                ! positive in case of erosion
                if (par%sourcesink==0) then
                   dzg=par%morfac*par%dt/(1.d0-par%por)*( & ! Dano, dz from sus transport gradients      
                        ( Susg(i,j,:)*dnu(i,j)-Susg(i-1,j,:)*dnu(i-1,j) +&           
                          Svsg(i,j,:)*dsv(i,j)-Svsg(i,j-1,:)*dsv(i,j-1) +&
                                ! dz from bed load transport gradients
                          Subg(i,j,:)*dnu(i,j)-Subg(i-1,j,:)*dnu(i-1,j)+&           
                          Svbg(i,j,:)*dsv(i,j)-Svbg(i,j-1,:)*dsv(i,j-1) )*dsdnzi(i,j)    )       
                elseif (par%sourcesink==1) then
                   dzg=par%morfac*par%dt/(1.d0-par%por)*( &
                          ero(i,j,:)-depo_ex(i,j,:)   +&
                        ( Subg(i,j,:)*dnu(i,j)-Subg(i-1,j,:)*dnu(i-1,j)+&           
                          Svbg(i,j,:)*dsv(i,j)-Svbg(i,j-1,:)*dsv(i,j-1) )*dsdnzi(i,j)    )
                endif

                ! erosion/deposition rate of sand mass (m/s)
                ! positive in case of erosion
                edg = dzg*(1.d0-par%por)/par%dt

                dz=>dzbed(i,j,:)
                pb=>pbbed(i,j,:,:)           

                call update_fractions(par,s,i,j,dz,pb,edg,sum(dzg))

             enddo ! nx+1
          enddo ! ny+1
       else
          j=1
          do i=2,nx
             ! bed level changes per fraction in this morphological time step in meters sand including pores
             ! positive in case of erosion
             if (par%sourcesink==0) then
                       dzg=par%morfac*par%dt/(1.d0-par%por)*( & ! Dano, dz from sus transport gradients  
                     ( Susg(i,j,:)*dnu(i,j)-Susg(i-1,j,:)*dnu(i-1,j) +&           
                       Subg(i,j,:)*dnu(i,j)-Subg(i-1,j,:)*dnu(i-1,j) )*dsdnzi(i,j)    )           
             elseif (par%sourcesink==1) then
                       dzg=par%morfac*par%dt/(1.d0-par%por)*( &
                       ero(i,j,:)-depo_ex(i,j,:)       +&
                      ( Subg(i,j,:)*dnu(i,j)-Subg(i-1,j,:)*dnu(i-1,j) )*dsdnzi(i,j)    )
                      
             endif

             ! erosion/deposition rate of sand mass (m/s)
             ! positive in case of erosion
             edg = dzg*(1.d0-par%por)/par%dt

             dz=>dzbed(i,j,:)
             pb=>pbbed(i,j,:,:)           

             call update_fractions(par,s,i,j,dz,pb,edg,sum(dzg))

          enddo ! nx+1
       endif !ny>0      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       ! Avalanching
       !
       
       
       if (par%avalanching==1) then
           do ii=1,nint(par%morfac)

              aval=.false.
              dzbdx=0.d0
              dzbdy=0.d0
              do j=1,ny+1
                 do i=1,nx
                    dzbdx(i,j)=(zb(i+1,j)-zb(i,j))/dsu(i,j)  
                 enddo
              enddo
              
              ! Include short wave runup
              
              ! Jaap crude switch, need to do further testing
              
              hav = hh
              indx = nx+1
              
              if (par%struct==1 .and. par%shoaldelay==1) then
              
              do j=1,ny+1
                 first = 0;
                 do i=1,nx
                    if (wetz(i,j)-wetz(max(i-1,1),j)==-1 .and. first==0) then ! transition from wet to dry
                        ! only consider first dry point
                        first = 1
                        ! find wave height for runup at facsd*L1 meter from water line 
                        call linear_interp(xz(:,j),H(:,j),nx+1,xz(i-1,j)-par%facsd*L1(i-1,j),s%Hrunup(1,j),indx)
                        ! Find toe of runup slope if present (dzbdx > 0.15). 
                        ! If not present Hrunup will converge to H at the water line (where H = 0 per definition)
                        do j1=indx,i-1
                          if (dzbdx(j1,j)<0.15d0 .or. structdepth(j1,j)>0.1d0) then
                             indx = j1
                          endif
                        enddo 
                        ! update Hrunup and runup x-location
                        s%Hrunup(1,j) = H(indx,j)
                        s%xHrunup(1,j) = xz(indx,j);
                        ! now itteratively compute runup
                        hav(:,j) = hh(:,j)
                        runup_old = huge(0.d0)
                        s%runup(1,j) = 0;
                        do while (abs(s%runup(1,j)-runup_old)>0.01d0)
                          runup_old = s%runup(1,j)
                          slopeind = 0
                          where (hav(:,j)>par%eps .and. dzbdx(:,j)>0.15)
                            slopeind = 1
                          endwhere
                          !bermind = 0
                          !where (slopeind == 0 .and. wetz(:,j) == 0 .and. hav(:,j)>par%eps)
                          !  bermind = 1
                          !endwhere
                          strucslope = sum(dzbdx(indx:nx,j)*dsu(indx:nx,j)*slopeind(indx:nx))/ &
                          max(par%eps,sum(dsu(indx:nx,j)*slopeind(indx:nx)))
                          if (strucslope > 0.d0) then         
                             irrb = strucslope/sqrt(2*par%px*max(s%Hrunup(1,j),par%eps)/par%g/par%Trep**2)
                             !bermwidth  = sum(dsu(indx:nx,j)*bermind(indx:nx))
                             !rb = bermwidth/(bermwidth+sum(dsu(indx:nx,j)*slopeind(indx:nx)))
                             !gamB = max(0.6d0,1.d0-rb)
                             !runup_max = (4.3d0-1.6d0/sqrt(irrb))*s%Hrunup(1,j)
                             !s%runup(1,j) = min(runup_max,irrb*s%Hrunup(1,j))*cos(2*par%px/par%Trep*par%t)
                             s%runup(1,j) = min(irrb,2.3d0)*s%Hrunup(1,j)*cos(2*par%px/par%Trep*par%t)
                          else
                             s%runup(1,j) = 0.d0;
                          endif
        
                          hav(:,j) = hh(:,j) + wetz(:,j)*par%shoaldelay*s%runup(1,j) + &
                                              (1.d0-wetz(:,j))*max(par%eps, par%shoaldelay*(s%runup(1,j)-zb(:,j)));
                       enddo
                    endif
                 enddo
              enddo
              
              endif ! end crude switch
              !
              do i=2,nx-1
                 do j=1,ny+1
                    !if (max( max(hh(i,j),par%delta*H(i,j)), max(hh(i+1,j),par%delta*H(i+1,j)) )>par%hswitch+par%eps) then
                    if(max(hav(i,j),hav(i+1,j))>par%hswitch+par%eps) then ! Jaap instead of hh
                       dzmax=par%wetslp;
                       if (i>indx) then ! tricks: seaward of indx (transition from sand to structure) wetslope is set to 0.03;
                          dzmax = 0.03d0
                          !dzmax = max(0.01d0,abs(dzbdx(i,j))*0.9d0)
                       endif
                    else 
                       dzmax=par%dryslp;
                    end if

                    if(abs(dzbdx(i,j))>dzmax .and. structdepth(i+nint(max(0.d0,sign(1.d0,dzbdx(i,j)))),j)>par%eps) then 
                       aval=.true.     
                       dzb=sign(1.0d0,dzbdx(i,j))*(abs(dzbdx(i,j))-dzmax)*dsu(i,j);
                       !Dano: Need to make this mass-conserving for curved areas; good enough for now
                       if (dzb >= 0.d0) then
                          ie = i+1                                        ! index erosion point
                          id = i                                          ! index deposition point
                          dxfac = dsz(i+1,j)/dsz(i,j)                     ! take into account varying gridsize
                          dzb=min(dzb,par%dzmax*par%dt/dsu(i,j))          ! make sure dzb is not in conflict with maximum erosion rate par%dzmax
                          dzb=min(dzb,structdepth(i+1,j))                 ! make sure dzb is not larger than sediment layer thickness
                       else
                          ie = i                                          ! index erosion point
                          id = i+1                                        ! index deposition point
                          dxfac = dsz(i,j)/dsz(i+1,j)                     ! take into account varying gridsize
                          dzb=max(dzb,-par%dzmax*par%dt/dsu(i,j)) 
                          dzb=max(dzb,-structdepth(i,j))
                       endif
                       
                       ! now fix fractions....
                       dz => dzbed(ie,j,:) 
                       pb => pbbed(ie,j,:,:)

                       ! figure out how many depth layers (ndz) are eroded in point iii
                       sdz = 0
                       ndz = 0
                       do while (sdz<abs(dzb))
                          ndz = ndz+1
                          sdz = sdz+dz(ndz)
                       enddo

                       ! now update bed and fractions by stepping through each layer seperately
                       dzleft = abs(dzb)
                       dzavt  = 0.d0
                       
                       do jdz=1,ndz

                          dzt = min(dz(jdz),dzleft)
                          dzleft = dzleft-dzt;

                          ! erosion deposition per fraction (upwind or downwind); edg is positive in case of erosion   
                          do jg=1,par%ngd 
                             edg2(jg) =  sedcal(jg)*dzt*pb(jdz,jg)*(1.d0-par%por)/par%dt       ! erosion    (dzt always > 0 )
                             edg1(jg) = -sedcal(jg)*edg2(jg)*dxfac                             ! deposition (dzt always < 0 )
                          enddo

                          dzavt = dzavt + sum(edg2)*par%dt/(1.d0-par%por)

                          call update_fractions(par,s,ie,j,dzbed(ie,j,:),pbbed(ie,j,:,:),edg2,dzavt)           ! update bed in eroding point

                          call update_fractions(par,s,id,j,dzbed(id,j,:),pbbed(id,j,:,:),edg1,-dzavt*dxfac)    ! update bed in deposition point

                       enddo

                       ! update water levels and dzav
                       zs(ie,j)  = zs(ie,j)-dzavt
                       dzav(ie,j)= dzav(ie,j)-dzavt

                       zs(id,j)  = zs(id,j)+dzavt*dxfac
                       dzav(id,j)= dzav(id,j)+dzavt*dxfac

                    end if
                 end do
              end do
              !JJ: update y slopes after avalanching in X-direction seems more appropriate
              do j=1,ny
                 do i=1,nx+1
                    dzbdy(i,j)=(zb(i,j+1)-zb(i,j))/dnv(i,j)
                 enddo
              enddo

              do j=2,ny-1
                 do i=1,nx+1
                    if(max(hh(i,j),hh(i,j+1))>par%hswitch+par%eps) then
                       dzmax=par%wetslp
                    else
                       dzmax=par%dryslp
                    end if
                    if(abs(dzbdy(i,j))>dzmax .and. structdepth(i,j+nint(max(0.d0,sign(1.d0,dzbdy(i,j)))))>par%eps) then ! Jaap
                       aval=.true. 
                       dzb=sign(1.0d0,dzbdy(i,j))*(abs(dzbdy(i,j))-dzmax)*dnv(i,j)
                       !
                       if (dzb >= 0.d0) then
                          je = j+1                                        ! index erosion point
                          jd = j                                          ! index deposition point
                          dyfac = dnz(i,j+1)/dnz(i,j)                     ! take into account varying gridsize
                          dzb=min(dzb,par%dzmax*par%dt/dnv(i,j))
                          dzb=min(dzb,structdepth(i,j+1))
                       else
                          je = j                                          ! index erosion point
                          jd = j+1                                        ! index deposition point
                          dyfac = dnz(i,j)/dnz(i,j+1)                     ! take into account varying gridsize
                          dzb=max(dzb,-par%dzmax*par%dt/dnv(i,j))
                          dzb=max(dzb,-structdepth(i,j))
                       endif
                       
                       dz => dzbed(i,je,:) 
                       pb => pbbed(i,je,:,:)

                       ! figure out how many depth layers (ndz) are affected
                       sdz = 0
                       ndz = 0
                       do while (sdz<abs(dzb))
                          ndz = ndz+1
                          sdz = sdz+dz(ndz)
                       enddo

                       ! now update bed and fractions by stepping through each layer seperately
                       dzleft = abs(dzb)
                       dzavt  = 0.d0

                       do jdz=1,ndz
                          dzt = min(dz(jdz),dzleft)
                          dzleft = dzleft-dzt;

                          ! erosion deposition per fraction (upwind or downwind); edg is positive in case of erosion  
                          do jg=1,par%ngd 
                             edg2(jg) = sedcal(jg)*dzt*pb(jdz,jg)*(1.d0-par%por)/par%dt        ! erosion    (dzt always > 0 )
                             edg1(jg) = -sedcal(jg)*edg2(jg)*dyfac                             ! deposition (dzt always < 0 )
                          enddo

                          dzavt = dzavt + sum(edg2)*par%dt/(1.d0-par%por)

                          call update_fractions(par,s,i,je,dzbed(i,je,:),pbbed(i,je,:,:),edg2,dzavt)            ! upwind point

                          call update_fractions(par,s,i,jd,dzbed(i,jd,:),pbbed(i,jd,:,:),edg1,-dzavt*dyfac)    ! downwind point

                       enddo

                       ! update water levels and dzav
                       zs(i,je)  = zs(i,je)-dzavt
                       dzav(i,je)= dzav(i,je)-dzavt

                       zs(i,jd)  = zs(i,jd)+dzavt*dyfac
                       dzav(i,jd)= dzav(i,jd)+dzavt*dyfac
                       
                    end if
                 end do
              end do
              if (.not.aval) exit
           end do
       end if
       !
       ! bed boundary conditions
       ! 
       if(xmpi_isleft .and. ny>0) then
          zb(:,1) = zb(:,2)
          sedero(:,1) = sedero(:,2)
          structdepth(:,1) = structdepth(:,2)
          pbbed(:,1,:,:)=pbbed(:,2,:,:)
          z0bed(:,1)=z0bed(:,2)
          dzbed(:,1,:)=dzbed(:,2,:)
       endif

       if(xmpi_isright .and. ny>0) then
          zb(:,ny+1) = zb(:,ny)
          sedero(:,ny+1) = sedero(:,ny)
          structdepth(:,ny+1) = structdepth(:,ny)
          pbbed(:,ny+1,:,:)=pbbed(:,ny,:,:)
          z0bed(:,ny+1)=z0bed(:,ny)
          dzbed(:,ny+1,:)=dzbed(:,ny,:)
       endif

       ! Robert: in parallel version bed update must take place on internal boundaries:
#ifdef USEMPI
       call xmpi_shift(zb,'1:')
       call xmpi_shift(zb,'m:')
       call xmpi_shift(zb,':1')
       call xmpi_shift(zb,':n')
#endif

       ! Update representative sed.diameter at the bed for flow friction and output
       if (par%ngd>1) then
          do j=j1,max(ny,1)
             do i=2,nx
                s%D50top =  sum(s%pbbed(i,j,1,:)*s%D50)
                s%D90top =  sum(s%pbbed(i,j,1,:)*s%D90)
             enddo
          enddo
       endif

    endif ! if par%t>par%morstart

  end subroutine bed_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_fractions(par,s,i,j,dz,pb,edg,dzb)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     Copyright (C) 2009 Technische Universiteit Delft
    !        Bram van Prooijen
    !        b.c.vanprooijen@tudelft.nl
    !	   +31(0)15 2784070   
    !        Faculty of Civil Engineering and Geosciences
    !        department of Hydraulic Engineering
    !	   PO Box 5048
    !        2600 GA Delft
    !        The Netherlands
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use params
    use spaceparams
    use xmpi_module


    IMPLICIT NONE

    type(spacepars),target                          :: s
    type(parameters)                                :: par

    integer                                         :: i,j,jg,jd
    real*8                                          :: ED,zbold,dzbt,fac,dzb_loc
    real*8, intent(in)                              :: dzb
    real*8 , dimension(par%ngd),intent(in)          :: edg
    real*8 , dimension(s%nd(i,j)),intent(inout)     :: dz
    real*8 , dimension(s%nd(i,j),par%ngd),intent(inout) :: pb

    real*8 , dimension(:),allocatable,save                :: Ap,b
    real*8 , dimension(:,:),allocatable,save              :: Sm,A

    if (.not. allocated(Ap)) then
       allocate(Ap   (s%nd(1,1)))
       allocate(b    (s%nd(1,1)))
       allocate(Sm   (s%nd(1,1),par%ngd))
       allocate(A    (s%nd(1,1),3))
    endif
    !TODO, dzb_loc is not initialized can be Nan, leading to infinite loop
    dzb_loc = dzb

!!!initialize Sm      

    ED = sum(edg)

    ! do t_sub=1,nt_sub  !loop over subtimesteps
    do while (abs(dzb_loc) .gt. 0.d0)
       ! dzb can be nan, check...
   !    if (isnan(dzb_loc)) then
    !      write(*,*) dzbt, dzb_loc
    !      write(*,*) 'dzb is nan'
    !   end if
       dzbt     = min(dzb_loc,dz(par%nd_var))                 ! make sure erosion (dzg is positive) is limited to thickness of variable layer
       dzbt     = max(dzbt,-par%frac_dz*dz(par%nd_var+1))     ! make sure deposition (dzg is negative) is limited to thickness of first layer below variable layer

       fac = dzbt/dzb_loc                                     ! factor over mass change in cells to limit erosion and deposition

       dzb_loc = dzb_loc-dzbt                                 ! update dzb

       do jg=1,par%ngd

          A=0.
          Ap=0.
          b=0.

          Sm(:,jg)=pb(:,jg)*dz*(1.d0-par%por)

!!!build matrix A							
          select case(par%nd_var)
          case(1)
             !in this case: A=0
          case(2)
             A (1,1:3)= (/0.d0           ,  min(ED,0.d0) ,  max(ED,0.d0) /)   
             A (2,1:3)= (/-min(ED,0.d0)  , -max(ED,0.d0) ,  0.d0         /)
          case(3:1000)
             A (1,1:3)= (/0.d0           ,  min(ED,0.d0) ,  max(ED,0.d0) /)   

             A (2:par%nd_var-1,1)=-min(ED,0.d0)
             A (2:par%nd_var-1,2)=-abs(ED)
             A (2:par%nd_var-1,3)=max(ED,0.d0)

             A (par%nd_var,1:3)= (/-min(ED,0.d0) , -max(ED,0.d0) ,   0.d0    /)
          end select

!!!determine RHS 

          ! Ap makes sure that layer nd_var varies in thickness in time
          ! Ap = 0 with single fraction (in case(1))
          Ap(1) = sum(A(1,2:3)*pb(1:2,jg))                   
          do jd = 2,par%nd_var                               
             Ap(jd) = sum(A (jd,:)*pb(jd-1:jd+1,jg))
          enddo
          Ap(par%nd_var+1) = sum(A(par%nd_var+1,1:2)*pb(par%nd_var:par%nd_var+1,jg))

          ! b represents the actual erosion and deposition in the top layer. 
          ! However, the thickness of the top layer remains constant and instead layer nd_var will breath in thickness 
          b(1) = -edg(jg)

!!!update Sm

          ! Sm is the sediment mass per fraction per layer
          ! Sm(1:par%nd_var+1,jg) = Sm(1:par%nd_var+1,jg) + par%dt/nt_sub*(Ap(1:par%nd_var+1)+b(1:par%nd_var+1))
          Sm(1:par%nd_var+1,jg) = Sm(1:par%nd_var+1,jg) + par%dt*fac*(Ap(1:par%nd_var+1)+b(1:par%nd_var+1))

       enddo !fractions

       ! From Sm we can compute the new fraction ratios per layer and the layer thickness...
       do jd=1,par%nd_var+1
          if (sum(Sm(jd,:))>0) then
             pb(jd,:) = Sm(jd,:)/sum(Sm(jd,:))
          else
             pb(jd,:) = pb(jd,:)
          endif
          dz(jd) = sum(Sm(jd,:))/(1.d0-par%por)
       enddo

!!! modify grid

       !merge two upper layers in case of erosion
       if (dz(par%nd_var) .lt. par%merge*dz(par%nd_var+1)) then
          forall (jg=1:par%ngd)
             pb(par%nd_var,jg) = (dz(par%nd_var)*pb(par%nd_var,jg) + dz(par%nd_var+1)* &
                  pb(par%nd_var+1,jg))/(dz(par%nd_var)+dz(par%nd_var+1)) 
             pb(par%nd_var+1:s%nd(i,j)-1,jg) = pb(par%nd_var+2:s%nd(i,j),jg)
             pb(s%nd(i,j),jg) = pb(s%nd(i,j),jg)
          endforall
          s%z0bed(i,j) = s%z0bed(i,j)-dz(par%nd_var+1)
          dz(par%nd_var) = dz(par%nd_var+1)+dz(par%nd_var)   
       endif
       !split upper layer in case of sedimentation
       if (dz(par%nd_var)>par%split*dz(par%nd_var+1)) then
          pb(par%nd_var+1:s%nd(i,j),:) = pb(par%nd_var:s%nd(i,j)-1,:)
          s%z0bed(i,j) = s%z0bed(i,j)+dz(par%nd_var+1)
          dz(par%nd_var) = dz(par%nd_var)-dz(par%nd_var+1)  
       endif
    enddo ! nt_sub

    pb = max(0.d0,min(pb,1.d0))

    zbold = s%zb(i,j)
    s%zb(i,j) = s%z0bed(i,j)+sum(dz)
    s%dzbdt(i,j) = s%dzbdt(i,j)+(s%zb(i,j)-zbold)
    s%sedero(i,j) = s%sedero(i,j)+(s%zb(i,j)-zbold)
    s%structdepth(i,j) = max(0.d0,s%structdepth(i,j)+(s%zb(i,j)-zbold))

  end subroutine update_fractions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sb_vr(s,par)
    use params
    use spaceparams
    use xmpi_module

    IMPLICIT NONE

    type(spacepars),target              :: s
    type(parameters)                    :: par

    integer                             :: i
    integer                             :: j,jg
    real*8                              :: z0,Ass,dcf,dcfin,ML
    real*8                              :: Te,kvis,Sster,c1,c2,wster
    real*8,save                         :: delta,onethird,twothird

    real*8 , dimension(:)  ,allocatable,save   :: w,dster
    real*8 , dimension(:,:),allocatable,save   :: vmg,Cd,Asb,dhdx,dhdy,Ts,hfac
    real*8 , dimension(:,:),allocatable,save   :: urms2,Ucr,term1,term2
    real*8 , dimension(:,:),allocatable,save   :: uandv,b,fslope,hloc,ceqs,ceqb

    include 's.ind'
    include 's.inp'

    if (.not. allocated(vmg)) then
       allocate (vmg   (nx+1,ny+1))
       allocate (Cd    (nx+1,ny+1))
       allocate (Asb   (nx+1,ny+1))
       allocate (dhdx  (nx+1,ny+1))   ! not used wwvv
       allocate (dhdy  (nx+1,ny+1))   ! not used wwvv
       allocate (urms2 (nx+1,ny+1))
       allocate (Ucr   (nx+1,ny+1))
       allocate (term1 (nx+1,ny+1))
       allocate (term2 (nx+1,ny+1))
       allocate (uandv (nx+1,ny+1))  ! not used wwvv
       allocate (b     (nx+1,ny+1))  ! not used wwvv
       allocate (fslope(nx+1,ny+1))  ! not used wwvv
       allocate (hloc  (nx+1,ny+1))
       allocate (Ts    (nx+1,ny+1))
       allocate (hfac   (nx+1,ny+1))
       allocate (ceqs   (nx+1,ny+1))
       allocate (ceqb   (nx+1,ny+1))

       allocate (w     (par%ngd))
       allocate (dster (par%ngd))
       vmg = 0.d0
       onethird=1.d0/3.d0
       twothird=2.d0/3.d0
       ! Robert: do only once, not necessary every time step
       do jg=1,par%ngd
          ! cjaap: compute fall velocity with simple expression from Ahrens (2000)
          Te    = 20.d0
          kvis  = 4.d0/(20.d0+Te)*1d-5 ! Van rijn, 1993 
          Sster = D50(jg)/(4*kvis)*sqrt((par%rhos/par%rho-1)*par%g*D50(jg))
          c1    = 1.06d0*tanh(0.064d0*Sster*exp(-7.5d0/Sster**2))
          c2    = 0.22d0*tanh(2.34d0*Sster**(-1.18d0)*exp(-0.0064d0*Sster**2))
          wster = c1+c2*Sster
          w(jg) = wster*sqrt((par%rhos/par%rho-1.d0)*par%g*D50(jg))
          ! RJ: for modeling gravel
          delta = (par%rhos-par%rho)/par%rho
          dster(jg)=(delta*par%g/1.d-12)**onethird*s%D50(jg)      
       enddo
    endif
    !
    hloc = max(hh,0.01)
    !
    if (par%swave==1) then
       ! wave breaking induced turbulence due to short waves
       do j=1,ny+1
          do i=1,nx+1
             ! compute mixing length
             ! ML = 2*R(i,j)*par%Trep/(par%rho*c(i,j)*max(H(i,j),par%eps))
             ML = dsqrt(2*R(i,j)*par%Trep/(par%rho*c(i,j)))
             ! ML = 0.9d0*H(i,j)
             ML = min(ML,hloc(i,j));
             ! exponential decay turbulence over depth
             dcfin = exp(min(100.d0,hloc(i,j)/max(ML,0.01d0)))
             dcf = min(1.d0,1.d0/(dcfin-1.d0))
             if (trim(par%turb) == 'bore_averaged') then
                kb(i,j) = par%nuhfac*(DR(i,j)/par%rho)**twothird*dcf*par%Trep/Tbore(i,j)
             elseif (trim(par%turb) == 'wave_averaged') then
                kb(i,j) = par%nuhfac*(DR(i,j)/par%rho)**twothird*dcf
             elseif (trim(par%turb) == 'none') then
                kb(i,j) = 0.0d0
             endif
          enddo
       enddo
    endif

    ! switch to include long wave stirring
    if (par%lws==1) then
       vmg  = dsqrt(ue**2+ve**2)
    elseif (par%lws==0) then
       ! vmg lags on actual mean flow; but long wave contribution to mean flow is included... 
       vmg = (1.d0-1.d0/par%cats/par%Trep*par%dt)*vmg + (1.d0/par%cats/par%Trep*par%dt)*dsqrt(ue**2+ve**2) 
    endif

    urms2  = urms**2+1.45d0*(kb+kturb)

    do jg = 1,par%ngd

       Ts       = par%tsfac*hloc/w(jg)
       Tsg(:,:,jg) = max(Ts,par%Tsmin) 

       ! calculate treshold velocity Ucr
       if(D50(jg)<=0.0005d0) then
          Ucr=0.19d0*D50(jg)**0.1d0*log10(4.d0*hloc/D90(jg))
       else if(D50(jg)<0.05d0) then   !Dano see what happens with coarse material
          Ucr=8.5d0*D50(jg)**0.6d0*log10(4.d0*hloc/D90(jg))
       else

#ifdef USEMPI
          write(*,'(a,i4)') 'In process',xmpi_rank
#endif
       end if
       ! drag coefficient
       z0 = par%z0
       Cd=(0.40d0/(log(max(hloc,10.d0*z0)/z0)-1.0d0))**2

       ! transport parameters
       Asb=0.005d0*hloc*(D50(jg)/hloc/(delta*par%g*D50(jg)))**1.2d0         ! bed load coefficent
       Ass=0.012d0*D50(jg)*dster(jg)**(-0.6d0)/(delta*par%g*D50(jg))**1.2d0 ! suspended load coeffient

       term1=(vmg**2+0.018/Cd*par%sws*urms2) ! Make 0.018/Cd is always smaller than the flow friction coefficient 

       ! reduce sediment suspensions for (inundation) overwash conditions with critical flow velocitties
       ! vmag2=min(vmag2,par%smax*par%C**2*D50(jg)*delta)
       ! vmag2=min(vmag2,par%smax*par%g/par%cf*D50(jg)*delta)            ! In terms of cf
       ! term1=sqrt(vmag2+0.018d0/Cd*urms2)     ! nearbed-velocity
       ! the two above lines are comment out and replaced by a limit on total velocity u2+urms2, robert 1/9 and ap 28/11

       term1=min(term1,par%smax*par%g/par%cf*s%D50(jg)*delta)
       term1=sqrt(term1)      

       term2 = 0.d0
       do j=1,ny+1
          do i=1,nx
             if(term1(i,j)>Ucr(i,j) .and. hh(i,j)>par%eps) then
                term2(i,j)=(term1(i,j)-Ucr(i,j))**2.4d0
             end if
          end do
       end do
       ! wwvv in parallel version, there will be a discrepancy between the values
       ! of term2. term2(nx+1,:) is zero, while the corresponding row in the process
       ! below term2(2,:) has some value, different from zero.
       ! so we fix this:
#ifdef USEMPI
       call xmpi_shift(term2,'m:')
#endif
       ceqb = Asb*term2
       ceqb = min(ceqb/hloc,par%cmax/2)             ! maximum equilibrium bed concentration
       ceqbg(:,:,jg) = (1-par%bulk)*ceqb*sedcal(jg)*wetz
       ceqs = Ass*term2
       ceqs = min(ceqs/hloc,par%cmax/2)             ! maximum equilibrium suspended concentration
       ceqsg(:,:,jg) = (ceqs+par%bulk*ceqb)*sedcal(jg)*wetz

       ! Jaap: old brute method to prevent strong coastline erosion
       ! where (hloc<=par%hmin) ceqsg(:,:,jg) = 0.d0

    enddo  ! end og grain size classes

  end subroutine sb_vr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sednew(s,par)
    use params
    use spaceparams
    use readkey_module
    use xmpi_module

    IMPLICIT NONE

    type(spacepars),target                  :: s
    type(parameters)                        :: par

    integer                                 :: i,j,jg
    real*8                                  :: Ass,dcf,dcfin,ML
    real*8                                  :: Te,kvis,Sster,cc1,cc2,wster
    real*8 , save                           :: delta,onethird,twothird
    real*8 , dimension(:),allocatable    ,save     :: w,dster  
    real*8 , dimension(:,:),allocatable  ,save     :: vmg,Asb,Ts
    real*8 , dimension(:,:),allocatable  ,save     :: urms2,Ucr,Ucrc,Ucrw,term1,B2,Cd
    real*8 , dimension(:,:),allocatable  ,save     :: hloc,ceqs,ceqb

    include 's.ind'
    include 's.inp'

    if (.not. allocated(vmg)) then
       allocate (vmg   (nx+1,ny+1))
       allocate (term1 (nx+1,ny+1))
       allocate (B2    (nx+1,ny+1))
       allocate (Cd    (nx+1,ny+1))
       allocate (Asb   (nx+1,ny+1))
       allocate (Ucr   (nx+1,ny+1))
       allocate (Ucrc  (nx+1,ny+1))
       allocate (Ucrw  (nx+1,ny+1))
       allocate (urms2 (nx+1,ny+1))
       allocate (hloc  (nx+1,ny+1))
       allocate (Ts    (nx+1,ny+1))
       allocate (ceqs  (nx+1,ny+1))
       allocate (ceqb  (nx+1,ny+1))
       allocate (w     (par%ngd))
       allocate (dster (par%ngd))
       vmg = 0.d0
       onethird=1.d0/3.d0
       twothird=2.d0/3.d0
       ! Robert: do only once, not necessary every time step
       do jg=1,par%ngd
          ! cjaap: compute fall velocity with simple expression from Ahrens (2000)
          Te    = 20.d0
          kvis  = 4.d0/(20.d0+Te)*1d-5 ! Van rijn, 1993 
          Sster = s%D50(jg)/(4*kvis)*sqrt((par%rhos/par%rho-1)*par%g*s%D50(jg))
          cc1   = 1.06d0*tanh(0.064d0*Sster*exp(-7.5d0/Sster**2))
          cc2    = 0.22d0*tanh(2.34d0*Sster**(-1.18d0)*exp(-0.0064d0*Sster**2))
          wster = cc1+cc2*Sster
          w(jg) = wster*sqrt((par%rhos/par%rho-1.d0)*par%g*s%D50(jg))
          ! RJ: for modeling gravel
          delta = (par%rhos-par%rho)/par%rho
          dster(jg)=(delta*par%g/1.d-12)**onethird*s%D50(jg)
       enddo
    endif

    ! hloc   = max(hh,0.01d0) ! Jaap 
    hloc = max(hh,0.01)
    !
    ! compute near bed turbulence
    !
    ! due to short waves

    if (par%swave==1) then
       ! wave breaking induced turbulence due to short waves
       do j=1,ny+1
          do i=1,nx+1
             ! compute mixing length
             ! ML = 2*R(i,j)*par%Trep/(par%rho*c(i,j)*max(H(i,j),par%eps))
             ML = dsqrt(2*R(i,j)*par%Trep/(par%rho*c(i,j)))
             ! ML = 0.9d0*H(i,j)
             ML = min(ML,hloc(i,j));
             ! exponential decay turbulence over depth
             dcfin = exp(min(100.d0,hloc(i,j)/max(ML,0.01d0)))
             dcf = min(1.d0,1.d0/(dcfin-1.d0))
             if (trim(par%turb) == 'bore_averaged') then
                kb(i,j) = par%nuhfac*(DR(i,j)/par%rho)**twothird*dcf*par%Trep/Tbore(i,j)
             elseif (trim(par%turb) == 'wave_averaged') then
                kb(i,j) = par%nuhfac*(DR(i,j)/par%rho)**twothird*dcf
             elseif (trim(par%turb) == 'none') then
                kb(i,j) = 0.0d0
             endif
          enddo
       enddo
    endif !par%swave == 1

    ! switch to include long wave stirring
    if (par%lws==1) then
       vmg  = dsqrt(ue**2+ve**2)
    elseif (par%lws==0) then
       ! vmg lags on actual mean flow; but long wave contribution to mean flow is included... 
       vmg = (1.d0-1.d0/par%cats/par%Trep*par%dt)*vmg + (1.d0/par%cats/par%Trep*par%dt)*dsqrt(ue**2+ve**2) 
    endif

    urms2 = urms**2.d0+1.45d0*(kb+kturb)

    do jg = 1,par%ngd

       Ts       = par%tsfac*hloc/w(jg)
       Tsg(:,:,jg) = max(Ts,par%Tsmin) 
       !
       ! calculate treshold velocity Ucr
       !
       if(s%D50(jg)<=0.0005) then
          Ucrc=0.19d0*D50(jg)**0.1d0*log10(4.d0*hloc/D90(jg))                           !Shields
          Ucrw=0.24d0*(delta*par%g)**0.66d0*s%D50(jg)**0.33d0*par%Trep**0.33d0          !Komar and Miller (1975)
       else if(s%D50(jg)<=0.002) then
          Ucrc=8.5d0*D50(jg)**0.6d0*log10(4.d0*hloc/D90(jg))                            !Shields
          Ucrw=0.95d0*(delta*par%g)**0.57d0*D50(jg)**0.43*par%Trep**0.14                !Komar and Miller (1975)
       else if(s%D50(jg)>0.002) then
          Ucrc=1.3d0*sqrt(delta*par%g*D50(jg))*(hloc/D50(jg))**(0.5d0*onethird)         !Maynord (1978) --> also Neill (1968) where 1.3d0 = 1.4d0
          Ucrw=0.95d0*(delta*par%g)**0.57d0*D50(jg)**0.43*par%Trep**0.14                !Komar and Miller (1975)
       end if
       B2 = vmg/max(vmg+dsqrt(urms2),par%eps)
       Ucr = B2*Ucrc + (1-B2)*Ucrw                                                     !Van Rijn 2007 (Bed load transport paper)

       ! transport parameters
       Asb=0.015d0*hloc*(s%D50(jg)/hloc)**1.2d0/(delta*par%g*s%D50(jg))**0.75d0        !bed load coefficent
       Ass=0.012d0*s%D50(jg)*dster(jg)**(-0.6d0)/(delta*par%g*s%D50(jg))**1.2d0        !suspended load coeffient

       ! Jaap: par%sws to set short wave stirring to zero
       ! Jaap: Van Rijn use Peak orbital flow velocity --> 0.64 corresponds to 0.4 coefficient regular waves Van Rijn (2007)
       term1=vmg**2+0.64d0*par%sws*urms2
       ! reduce sediment suspensions for (inundation) overwash conditions with critical flow velocitties
       term1=min(term1,par%smax*par%g/par%cf*s%D50(jg)*delta)
       term1=sqrt(term1)                                        

       ceqb = 0.d0*term1                                                                     !initialize ceqb
       ceqs = 0.d0*term1                                                                     !initialize ceqs
       do j=1,ny+1
          do i=1,nx
             if(term1(i,j)>Ucr(i,j) .and. hloc(i,j)>par%eps) then
                ceqb(i,j)=Asb(i,j)*(term1(i,j)-Ucr(i,j))**1.5
                ceqs(i,j)=Ass*(term1(i,j)-Ucr(i,j))**2.4
             end if
          end do
       end do

       ceqb = min(ceqb/hloc,par%cmax/2) ! maximum equilibrium bed concentration
       ceqbg(:,:,jg) = (1-par%bulk)*ceqb*sedcal(jg)*wetz
       ceqs = min(ceqs/hloc,par%cmax/2) ! maximum equilibrium suspended concentration
       ceqsg(:,:,jg) = (ceqs+par%bulk*ceqb)*sedcal(jg)*wetz       

       ! Jaap: old brute method to prevent strong coastline erosion
       ! where (hloc<=par%hmin) ceqsg(:,:,jg) = 0.d0

    enddo                                 ! end of grain size classes
    ! end of grain size classes

  end subroutine sednew

  subroutine longwaveturb(s,par)
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
    use params
    use spaceparams
    use xmpi_module

    IMPLICIT NONE

    type(spacepars),target                   :: s
    type(parameters)                         :: par

    integer                                  :: i
    integer                                  :: j
    real                                     :: ksource
    real*8,dimension(:,:),allocatable,save   :: kturbu,kturbv,Sturbu,Sturbv,dzsdt_cr

    include 's.ind'
    include 's.inp'

    if (.not. allocated(kturbu)) then
       allocate(kturbu (nx+1,ny+1))
       allocate(kturbv (nx+1,ny+1))
       allocate(Sturbu (nx+1,ny+1))
       allocate(Sturbv (nx+1,ny+1))
       allocate(dzsdt_cr (nx+1,ny+1))
    endif
    ! use lagrangian velocities
    kturbu       = 0.0d0  !Jaap
    kturbv       = 0.0d0  !Jaap		
    dzsdt_cr=par%beta*c
    ! Update roller thickness
    rolthick=rolthick+par%dt*(abs(dzsdt)-dzsdt_cr)
    rolthick=max(rolthick,0.d0)
    !
    !  X-direction
    do j=1,ny+1
       do i=1,nx
          if(uu(i,j)>0.d0) then
             kturbu(i,j)=par%thetanum*kturb(i,j)+(1.d0-par%thetanum)*kturb(min(i+1,nx),j)
          elseif(uu(i,j)<0.d0) then
             kturbu(i,j)=par%thetanum*kturb(i+1,j)+(1.d0-par%thetanum)*kturb(max(i,2),j)
          else
             kturbu(i,j)=0.5d0*(kturb(i,j)+kturb(i+1,j))
          endif
       enddo
    enddo
    kturbu(nx+1,:) = kturb(nx+1,:) !Robert
    ! wwvv fix this in parallel case
#ifdef USEMPI
    call xmpi_shift(kturbu,'m:')
#endif
    !
    Sturbu=kturbu*uu*hu*wetu   !
    !
    ! Y-direction
    !
    do j=1,ny
       do i=1,nx+1
          if(vv(i,j)>0) then
             kturbv(i,j)=par%thetanum*kturb(i,j)+(1.d0-par%thetanum)*kturb(i,min(j+1,ny))
          else if(vv(i,j)<0) then
             kturbv(i,j)=par%thetanum*kturb(i,j+1)+(1.d0-par%thetanum)*kturb(i,max(j,2))
          else
             kturbv(i,j)=0.5d0*(kturbv(i,j)+kturbv(i,j+1))
          end if
       end do
    end do
    kturbv(:,ny+1) = kturb(:,ny+1) !Robert
    !
    Sturbv=kturbv*vv*hv*wetv
    do j=2,ny+1
       do i=2,nx+1
          ksource=par%g*rolthick(i,j)*par%beta*c(i,j)      !only important in shallow water, where c=sqrt(gh)  
          kturb(i,j) = hold(i,j)*kturb(i,j)-par%dt*(       &
	       (Sturbu(i,j)*dnu(i,j)-Sturbu(i-1,j)*dnu(i-1,j)+&
               Sturbv(i,j)*dsv(i,j)-Sturbv(i,j-1)*dsv(i,j-1))*dsdnzi(i,j)-&
               (ksource-par%betad*kturb(i,j)**1.5d0))
          kturb(i,j)=max(kturb(i,j),0.0d0)
       enddo
    enddo
    ! Jaap
    kturb = kturb/max(hh,0.01d0)
    kturb(1,:)=kturb(2,:)
    kturb(:,1)=kturb(:,2)
    kturb(nx+1,:)=kturb(nx+1-1,:)
    kturb(:,ny+1)=kturb(:,ny+1-1)
    ! wwvv fix the first and last rows and columns of kturb in parallel case
#ifdef USEMPI
    call xmpi_shift(kturb,'1:')
    call xmpi_shift(kturb,'m:')
    call xmpi_shift(kturb,':1')   ! Maybe not necessary as zb calculated from 1:ny+1, but just in case...
    call xmpi_shift(kturb,':n')   ! Dito
#endif
    kturb=kturb*wetz
    !

  end subroutine longwaveturb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine RvR(s,par)

    use params
    use spaceparams
    use xmpi_module

    IMPLICIT NONE

    type(spacepars),target                   :: s
    type(parameters)                         :: par

    real*8 , save                            :: m1,m2,m3,m4,m5,m6,alpha,beta

    real*8 , dimension(:,:),allocatable,save   :: Ur,Bm,B1

    include 's.ind'
    include 's.inp'

    ! only in first timestep..
    if (.not. allocated(Ur)) then

       allocate (Ur    (nx+1,ny+1))
       allocate (Bm    (nx+1,ny+1))
       allocate (B1    (nx+1,ny+1))

       m1 = 0;       ! a = 0
       m2 = 0.7939;  ! b = 0.79 +/- 0.023
       m3 = -0.6065; ! c = -0.61 +/- 0.041
       m4 = 0.3539;  ! d = -0.35 +/- 0.032 
       m5 = 0.6373;  ! e = 0.64 +/- 0.025
       m6 = 0.5995;  ! f = 0.60 +/- 0.043
       alpha = -log10(exp(1.d0))/m4
       beta  = exp(m3/m4)

    endif

    Ur = 3.d0/8.d0*sqrt(2.d0)*H*k/(k*hh)**3                    !Ursell number
    Ur = max(Ur,0.000000000001d0)
    Bm = m1 + (m2-m1)/(1.d0+beta*Ur**alpha)                    !Boltzmann sigmoid (eq 6)         
    B1 = (-90.d0+90.d0*tanh(m5/Ur**m6))*par%px/180.d0
    Sk = Bm*cos(B1)                                            !Skewness (eq 8)
    As = Bm*sin(B1)                                            !Asymmetry(eq 9)                                
    ua = par%sws*(par%facSk*Sk-par%facAs*As)*urms


  end subroutine RvR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

  subroutine vT(s,par)

    use params
    use spaceparams
    use readkey_module
    use xmpi_module

    IMPLICIT NONE

    type(spacepars),target                   :: s
    type(parameters)                         :: par

    integer                                  :: i,j
    integer , save                           :: nh,nt    
    integer                                  :: ih0,it0,ih1,it1
    real*8                                   :: p,q,f0,f1,f2,f3
    real*8                                   :: t0fac,siguref,duddtmax,dudtmax,duddtmean,dudtmean,detadxmean
    real*8 , save                            :: dh,dt

    real*8 , dimension(:,:),allocatable  ,save     :: h0,t0,detadxmax
    ! Robert: RF table now included in source code, rather than read from file
    ! Rienecker Fenton table with amongst others amplitudes non-linear components obtained with stream function theory
    include 'RF.inc'
    ! Robert: 'RF.inc' contains definition of RF as real*8(5,33,40) with "parameter" attribute
    !   so RF values may not be modified! To save memory, only rows 13,14,15,16 and 18 of the 
    !   original matrix are stored. So new row 1 corresponds with old row 13, etc.

    include 's.ind'
    include 's.inp'


    ! only in first timestep..
    if (.not. allocated(h0)) then
       allocate (h0    (nx+1,ny+1))
       allocate (t0    (nx+1,ny+1))
       allocate (detadxmax    (nx+1,ny+1))
       dh = 0.03d0
       dt = 1.25d0
       nh = floor(0.99d0/dh);
       nt = floor(50.d0/dt);
    endif

    ! non-linearity of short waves is listed in table as function of dimensionless wave height h0 and dimensionless wave period t0 

    ! compute dimensionless wave height and wave period in each grid point..
    h0 = min(nh*dh,max(dh,min(H,hh)/hh))
    t0 = min(nt*dt,max(dt,par%Trep*sqrt(par%g/hh)))

    ! estimate Sk, As and ua by interpolating table values
    do j=1,ny+1
       do i=1,nx+1
          ! interpolate table values....
          ih0=floor(h0(i,j)/dh);
          it0=floor(T0(i,j)/dt);
          ih1=min(ih0+1,nh);
          it1=min(it0+1,nt);
          p=(h0(i,j)-ih0*dh)/dh;
          q=(T0(i,j)-it0*dt)/dt;

          f0=(1-p)*(1-q);
          f1=p*(1-q);
          f2=q*(1-p);
          f3=p*q;

          ! Skewness and assymetry
          Sk(i,j) = f0*RF(1,ih0,it0)+f1*RF(1,ih1,it0)+ f2*RF(1,ih0,it1)+f3*RF(1,ih1,it1)
          As(i,j) = f0*RF(2,ih0,it0)+f1*RF(2,ih1,it0)+ f2*RF(2,ih0,it1)+f3*RF(2,ih1,it1)

          ! Sediment advection velocity from Skewness and Assymetry
          ! ua(i,j) = par%sws*par%facua*(Sk(i,j)-As(i,j))*urms(i,j)
          ua(i,j) = par%sws*(par%facSk*Sk(i,j)-par%facAs*As(i,j))*urms(i,j)

          ! Estimate bore period Tbore and mean slope bore front to feeded back in roller energy balance

          ! correct slope in case 1.25>T0>50
          if (t0(i,j)==50.d0) then
             t0fac = 50.d0/max((par%Trep*sqrt(par%g/hh(i,j))),50.d0) 
          elseif (t0(i,j)==1.25)then
             t0fac = 1.25d0/min((par%Trep*sqrt(par%g/hh(i,j))),1.25d0) 
          else
             t0fac = 1.d0
          endif

          ! detadxmax for Tbore...
          ! dimnesionless maximum acceleration under bore front 
          duddtmax = f0*RF(3,ih0,it0)+f1*RF(3,ih1,it0)+ f2*RF(3,ih0,it1)+f3*RF(3,ih1,it1)
          siguref = f0*RF(4,ih0,it0)+f1*RF(4,ih1,it0)+ f2*RF(4,ih0,it1)+f3*RF(4,ih1,it1)
          ! translate dimensionless duddtmax to real world dudtmax
          !         /scale with variance and go from [-] to [m/s^2]     /tableb./dimensionless dudtmax
          dudtmax = urms(i,j)/max(par%eps,siguref)*sqrt(par%g/hh(i,j))*t0fac*duddtmax
          detadxmax(i,j) = dudtmax*sinh(k(i,j)*hh(i,j))/max(c(i,j),sqrt(H(i,j)*par%g))/sigm(i,j)

          ! detadxmean for roller energy balance dissipation...
          if (par%rfb==1) then
             duddtmean = f0*RF(5,ih0,it0)+f1*RF(5,ih1,it0)+ f2*RF(5,ih0,it1)+f3*RF(5,ih1,it1)
             dudtmean = urms(i,j)/max(par%eps,siguref)*sqrt(par%g/hh(i,j))*t0fac*duddtmean
             detadxmean = dudtmean*sinh(k(i,j)*hh(i,j))/max(c(i,j),sqrt(H(i,j)*par%g))/sigm(i,j)
             BR(i,j) = par%BRfac*sin(atan(detadxmean))
          endif

       enddo
    enddo

    Tbore = max(par%Trep/25.d0,min(par%Trep/4.d0,H/(max(c,sqrt(H*par%g))*max(detadxmax,par%eps))))
    Tbore = par%Tbfac*Tbore

  end subroutine vT

#ifdef HAVE_CONFIG_H
#ifndef HAVE_FORTRAN_ISNAN
  ! define a isnan function based on the 2003 extension ieee_is_nan for portland group
  ! If your compiler doesn't support isnan and you are not using portland group, it's probably time to upgrade your compiler (see README/INSTALL for details)
  logical function isnan(a)
    use ieee_arithmetic 

    real*8 a
    if (ieee_is_nan(a)) then
       isnan = .true.
    else
       isnan = .false.
    end if
    return
  end function isnan
#endif
#endif

end module morphevolution
