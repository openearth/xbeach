module boundaryconditions
   use typesandkinds
   implicit none
   save
contains
   subroutine wave_bc(sg,sl,par)
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
      use waveparams
      use spaceparams
      use interp
      use wave_functions_module
      use xmpi_module
      use readkey_module
      use logging_module
      use spectral_wave_bc_module
      use nonh_module, only: nonh_init_wcoef
      use paramsconst

      implicit none

      type(spacepars), target                     :: sg, sl
      type(spacepars), pointer                    :: s
      type(parameters)                            :: par

      integer, save                               :: nt
      integer                                     :: i,new,reclen,wordsize
      integer, save                               :: old
      integer, save                               :: recpos, curline
      integer                                     :: j
      integer                                     :: itheta
      integer                                     :: E_idx
      integer                                     :: ier,ier2
      real*8                                      :: E1,ei,dum,Hm0, dum1, spreadpar, bcdur, dum2
      real*8, save                                :: dtbcfile,rt,bcendtime,bcstarttime
      real*8                                      :: em,tshifted,tnew
      real*8, save                                :: Emean,Llong
      real*8,dimension(:)     ,allocatable,save   :: e01,L0,L,Lest,kbw,wbw,tanhkhwb,kxmwt
      real*8,dimension(:)     ,allocatable,save   :: fac1,fac2
      real*8,dimension(:)     ,allocatable,save   :: tE,dataE,databi
      real*8,dimension(:,:)   ,allocatable,save   :: ht
      real*8,dimension(:,:)   ,allocatable,save   :: q1,q2,q
      real*8,dimension(:)     ,allocatable,save   :: wcrestpos
      real*8,dimension(:,:)   ,allocatable,save   :: ee1, ee2
      real*8,dimension(:,:)   ,allocatable,save   :: gq1,gq2,gq
      real*8,dimension(:,:)   ,allocatable,save   :: gee1, gee2
      character(len=1)                            :: bl
      character(slen),save                        :: ebcfname,qbcfname,nhbcfname,esbcfname
      real*8                                      :: E0,arms
      real*8,dimension(:),allocatable,save        :: dist,factor
      logical                                     :: startbcf
      logical                                     :: isSet_U,isSet_V,isSet_Z,isSet_W,isSet_dU,isSet_dV,isSet_Q
      real*8,dimension(:)     ,allocatable,save   :: uig,vig,zig,wig,duig,dvig
      real*8,dimension(:)     ,allocatable,save   :: tempu,tempv

      ! Initialize to false only once....
      logical, save                               :: bccreated = .false.

      s=>sl
#ifndef USEMPI
      s=>sg
#endif

      if (.not. allocated(fac1)) then
         allocate(ht      (2,s%ny+1))
         allocate(fac1    (s%nx))
         allocate(fac2    (s%nx))
         allocate(e01(1:s%ntheta))
         allocate(dist(1:s%ntheta))
         allocate(factor(1:s%ntheta))
         allocate(wcrestpos(s%nx+1))
         allocate(L(s%ny+1))
         allocate(L0(s%ny+1))
         allocate(Lest(s%ny+1))
         L0 = par%g*par%Trep**2/2/par%px
         L = L0
         allocate(kbw(s%ny+1))
         allocate(tanhkhwb(s%ny+1))
         allocate(kxmwt(s%ny+1))
         allocate(wbw(s%ny+1))
         wbw = 2*par%px/par%Trep
      endif
      !
      !  BLOCK I: GENERATE AND READ-IN WAVE BOUNDARY CONDITIONS
      !
      startbcf=.false.
      if(.not. bccreated ) then
         if (xmaster) then
            call writelog('ls','','Setting up boundary conditions')
         endif
         bccreated=.true.
         startbcf=.true.                     ! trigger read from bcf for instat 3,4,5,7
         bcendtime=huge(0.0d0)               ! initial assumption for instat 3,4,5,7
         s%newstatbc=1
         ! Part A) read-in / generate different wbctypes (and waveforms)
         ! 1) TS_1
         if (par%wbctype==WBCTYPE_TS_1) then
            if(xmaster) then
               open( unit=7, file='bc/gen.ezs')
            endif
            if (xmaster) then
5              continue
               read(7,'(a)',iostat=ier)bl
               if (ier .ne. 0) then
                  call report_file_read_error('bc/gen.ezs')
               endif
               if(bl.eq.'*') goto 5
               read(7,*,iostat=ier)nt
               if (ier .ne. 0) then
                  call report_file_read_error('bc/gen.ezs')
               endif
            endif
#ifdef USEMPI
            call xmpi_bcast(nt)
#endif
            allocate(dataE  (nt))
            allocate(tE     (nt))
            do i=1,nt
               if(xmaster) then
                  read(7,*,iostat=ier) tE(i),dum,dataE(i)
                  if (ier .ne. 0) then
                     call report_file_read_error('bc/gen.ezs')
                  endif
               endif
#ifdef USEMPI
               call xmpi_bcast(tE(i))
               call xmpi_bcast(dataE(i))
#endif
            end do
            if (xmaster) then
               close(7)
            endif
            Emean=sum(dataE)/nt
            call dispersion(par,s)
            ! 2) TS_2
         elseif (par%wbctype==WBCTYPE_TS_2) then
            if (xmaster) then
               open( unit=7, file='bc/gen.ezs')
6              continue
               read(7,'(a)',iostat=ier)bl
               if (ier .ne. 0) then
                  call report_file_read_error('bc/gen.ezs')
               endif
               if(bl.eq.'*') goto 6
               read(7,*,iostat=ier)nt
               if (ier .ne. 0) then
                  call report_file_read_error('bc/gen.ezs')
               endif
            endif
#ifdef USEMPI
            call xmpi_bcast(nt)
#endif
            allocate(dataE  (nt))
            allocate(databi (nt))
            allocate(tE     (nt))
            do i=1,nt
               if(xmaster) then
                  read(7,*,iostat=ier) tE(i),databi(i),dataE(i)
                  if (ier .ne. 0) then
                     call report_file_read_error('bc/gen.ezs')
                  endif
               endif
#ifdef USEMPI
               call xmpi_bcast(tE(i))
               call xmpi_bcast(databi(i))
               call xmpi_bcast(dataE(i))
#endif
            end do
            if (xmaster) then
               close(7)
            endif
            Emean=sum(dataE)/nt
            call dispersion(par,s)
            ! 3) jons_table in combination with stationary
         elseif (par%wbctype==WBCTYPE_JONS_TABLE .and. par%wavemodel==WAVEMODEL_STATIONARY) then
            if (xmaster) then
               open( unit=7, file=par%bcfile)
               read(7,*,iostat=ier) Hm0, par%Trep,par%dir0, dum1, spreadpar, bcendtime, dum2
               if (ier .ne. 0) then
                  call report_file_read_error(par%bcfile)
               endif
               par%Hrms = Hm0/sqrt(2.d0)
               par%m = int(2*spreadpar)
               if (par%morfacopt==1) bcendtime=bcendtime/max(par%morfac,1.d0)
               s%theta0=(1.5d0*par%px)-par%dir0*atan(1.d0)/45.d0
               do while(s%theta0<-par%px)
                  s%theta0=s%theta0+2.d0*par%px
               enddo
               do while(s%theta0>par%px)
                  s%theta0=s%theta0-2.d0*par%px
               enddo
            endif
            s%newstatbc=1
#ifdef USEMPI
            call xmpi_bcast(bcendtime)
            call xmpi_bcast(par%Hrms)
            call xmpi_bcast(par%Trep)
            call xmpi_bcast(par%m)
            call xmpi_bcast(s%theta0)
#endif
            do itheta=1,s%ntheta
               s%sigt(:,:,itheta) = 2.d0*par%px/par%Trep
            end do
            s%sigm = sum(s%sigt,3)/s%ntheta
            call dispersion(par,s)
            
            ! 4)  jons_table in combination with surfbeat or nonh
         elseif ((par%wbctype==WBCTYPE_PARAMETRIC .or. par%wbctype==WBCTYPE_JONS_TABLE) .and. &
         (par%wavemodel/=WAVEMODEL_STATIONARY .and. xmaster)) then
            call spectral_wave_bc(sg,par,curline)
            ! 5) SWAN
         elseif (par%wbctype==WBCTYPE_SWAN.and.xmaster) then
            call spectral_wave_bc(sg,par,curline)
            ! 6) vardens
         elseif (par%wbctype==WBCTYPE_VARDENS.and.xmaster) then
            call spectral_wave_bc(sg,par,curline)
            ! 7) reuse
         elseif (par%wbctype==WBCTYPE_REUSE.and.xmaster) then
            curline = 1
            ! 8) ts_nonh
         elseif (par%wbctype==WBCTYPE_TS_NONH) then
            if (.not. allocated(uig)) then
               if (xmaster) then
                  allocate(uig(sg%ny+1))
                  allocate(vig(sg%ny+1))
                  allocate(zig(sg%ny+1))
                  allocate(wig(sg%ny+1))
                  allocate(duig(sg%ny+1))
               else  ! only need valid address for MPI-distribution
                  allocate(uig(1))
                  allocate(vig(1))
                  allocate(zig(1))
                  allocate(wig(1))
                  allocate(duig(1))
               endif
            endif
            if (xmaster) then
               call velocity_Boundary('boun_U.bcf',uig,vig,zig,wig,duig,dvig,sg%ny,par%t, &
               isSet_U,isSet_V,isSet_Z,isSet_W,isSet_dU,isSet_dV,isSet_Q)
            endif
         endif  ! done with the different wbctypes
         !
         ! Part B) Directional distribution for non-spectral wbctypes
         if (par%wbctype==WBCTYPE_PARAMS .or. par%wbctype==WBCTYPE_TS_1 .or. par%wbctype==WBCTYPE_TS_2 .or. &
         (par%wbctype==WBCTYPE_JONS_TABLE .and. par%wavemodel==WAVEMODEL_STATIONARY)  ) then
            dist=(cos(s%theta-s%theta0))**par%m
            do i=1,s%ntheta
               if(cos(s%theta(i)-s%theta0)<0.d0) then
                  dist(i)=0.0d0
               end if
            end do
            if (par%wbctype==WBCTYPE_TS_1 .or. par%wbctype==WBCTYPE_TS_2) then
               par%Hrms=sqrt(8*Emean/(par%rho*par%g))
            endif
            E0=par%rhog8*par%Hrms**2
            ! Energy density distribution
            if (sum(dist)>0.d0) then
               factor = (dist/sum(dist))/s%dtheta
            else
               factor=0.d0
            endif
            e01    = factor*E0;
            e01    = max(e01,0.0d0);
            Llong  = par%Tlong*s%cg(1,1)
         endif
         if (xmaster) then
            call writelog('sl','','Boundary conditions complete, starting computation')
         endif
      end if ! if bc was not created
      !
      ! BLOCK II: Recalculate bcf-file
      !
      if (par%t .ge. bcendtime) then
         ! Go over the different wbctypes (and wavemodel)
         ! 1) jons_table in combination with stationary
         if (par%wbctype==WBCTYPE_JONS_TABLE .and. par%wavemodel==WAVEMODEL_STATIONARY) then
            if (xmaster) then
               call writelog('ls','','Reading new wave conditions')
               ! read new wave conditions until bcendtime exceeds current time step
               ! necessary in case of external time management
               do while (par%t .ge. bcendtime)
                  read(7,*,iostat=ier) Hm0, par%Trep,par%dir0, dum1, spreadpar, bcdur, dum2
                  if (ier .ne. 0) then
                     call report_file_read_error(par%bcfile)
                  endif
                  par%Hrms = Hm0/sqrt(2.d0)
                  ! set taper time to 1 second for new conditions (likewise in waveparams)
                  par%taper = 1.d0
                  par%m = int(2*spreadpar)
                  if (par%morfacopt==1) then
                     bcendtime=bcendtime+bcdur/max(par%morfac,1.d0)
                  else
                     bcendtime=bcendtime+bcdur
                  endif
                  s%theta0=(1.5d0*par%px)-par%dir0*atan(1.d0)/45.d0
                  ! Make shore theta0 is also set between -px and px for new wave conditions
                  do while(s%theta0<-par%px)
                     s%theta0=s%theta0+2.d0*par%px
                  enddo
                  do while(s%theta0>par%px)
                     s%theta0=s%theta0-2.d0*par%px
                  enddo
               enddo
            endif
            s%newstatbc=1
#ifdef USEMPI
            call xmpi_bcast(bcendtime)
            call xmpi_bcast(par%Hrms)
            call xmpi_bcast(par%Trep)
            call xmpi_bcast(par%m)
            call xmpi_bcast(s%theta0)
#endif
            do itheta=1,s%ntheta
               s%sigt(:,:,itheta) = 2*par%px/par%Trep
            end do
            s%sigm = sum(s%sigt,3)/s%ntheta
            call dispersion(par,s)
            ! Directional distribution
            dist=(cos(s%theta-s%theta0))**par%m
            do i=1,s%ntheta
               if(cos(s%theta(i)-s%theta0)<0.d0) then
                  dist(i)=0
               end if
            end do
            E0=par%rhog8*par%Hrms**2
            ! Energy density distribution
            if (sum(dist)>0.d0) then
               factor = (dist/sum(dist))/s%dtheta
            else
               factor=0.d0
            endif
            e01    = factor*E0;
            e01    = max(e01,0.0d0);
            ! 2) JONSWAP (table)
         elseif ((par%wbctype==WBCTYPE_PARAMETRIC .or. par%wbctype==WBCTYPE_JONS_TABLE) .and. &
         (par%wavemodel/=WAVEMODEL_STATIONARY .and. xmaster)) then
            close(71)
            close(72)
            call spectral_wave_bc(sg,par,curline)
            startbcf=.true.
            ! 3) SWAN
         elseif (par%wbctype==WBCTYPE_SWAN.and.xmaster) then
            close(71)
            close(72)
            call spectral_wave_bc(sg,par,curline)
            startbcf=.true.
            ! 4) vardens 
         elseif (par%wbctype==WBCTYPE_VARDENS.and.xmaster) then
            close(71)
            close(72)
            call spectral_wave_bc(sg,par,curline)
            startbcf=.true.
            ! 5) reuse 
         elseif (par%wbctype==WBCTYPE_REUSE.and.xmaster) then
            close(71)
            close(72)
            startbcf=.true.
            if (par%t <= (par%tstop-par%dt)) then
               ! this is not compatible with external time management
               curline = curline + 1
            end if
         end if ! go over the different wbctypes
#ifdef USEMPI
         call xmpi_bcast(startbcf)
#endif
      end if ! done recalculating bcf-file
      !
      ! BLOCK III: COMPUTE WAVE BOUNDARY CONDITIONS CURRENT TIMESTEP
      !  
      ! Part A) Stationary: go over the different waveforms
      if (par%wavemodel==WAVEMODEL_STATIONARY) then
         ! Go over all excepted wbctypes (all the same)
         if (par%wbctype == WBCTYPE_PARAMS .or. par%wbctype == WBCTYPE_JONS_TABLE) then
            if(par%nonhspectrum .ne. 1) then  ! Robert: needed because nonhspectrum can be -123 if not initialised
               do j=1,s%ny+1
                  s%ee(1,j,:)=e01*min(par%t/par%taper,1.0d0)
                  s%bi(1) = 0.0d0
                  s%ui(1,j) = 0.0d0
               end do
            else
               wbw = 2*par%px/par%Trep
               do j=1,s%ny+1
                  Lest(j) = L(j)
                  L(j) = iteratedispersion(L0(j),Lest(j),par%px,s%hh(1,j))
               enddo
               kbw = 2*par%px/L
               arms = par%Hrms/2
               if(par%order==1) then
                  s%zi(1,:) = arms*dsin(wbw*par%t-kbw*( dsin(s%theta0)*(s%yz(1,:)-s%yz(1,1)) &
                                                       +dcos(s%theta0)*(s%xz(1,:)-s%xz(1,1)) &
                                                       ) &
                                        )
                  s%ui(1,:) = 1.d0/s%hh(1,:)*wbw*arms/sinh(kbw*s%hh(1,:)) * &
                                             dsin(wbw*par%t &
                                                  -kbw*( dsin(s%theta0)*(s%yz(1,:)-s%yz(1,1)) &
                                                        +dcos(s%theta0)*(s%xz(1,:)-s%xz(1,1)) &
                                                        ) &
                                                  ) * &
                              1.d0/kbw*sinh(kbw*s%hh(1,:)) *cosd(s%theta0-s%alfaz(1,:))
               elseif(par%order==2) then
                  ! 2nd order Stokes wave: To do more extensive testing
                  tanhkhwb = tanh(kbw*s%hh(1,:))
                  kxmwt = kbw*(dsin(s%theta0)*(s%yz(1,:)-s%yz(1,1))+dcos(s%theta0)*(s%xz(1,:)-s%xz(1,1))) - wbw*par%t
                  s%zi(1,:) = arms*(dcos(kxmwt) + kbw*arms*(3-tanhkhwb**2)/(4*tanhkhwb**3)*dcos(2*kxmwt))
                  s%ui(1,:) = (wbw/kbw)/s%hh(1,:)*s%zi(1,:) *cosd(s%theta0-s%alfaz(1,:))
               elseif(par%order==3) then
                  ! 3rd order Stokes wave: only deep water
                  !s%zi(1,:) = arms*(dcos(wbw *par%t) + 0.5d0*kbw*arms*dcos(2*wbw*par%t) + 0.375d0*(kbw*arms)**2*dcos(3*wbw*par%t))
                  !s%ui(1,:) = (wbw/kbw)/s%hh(1,:)*s%zi(1,:)
               endif
            endif
         else
         endif   ! Go over the different wbctypes for wavemodel = stationary
         !
         ! Part B) Surfbeat: go over the different waveforms
      elseif (par%wavemodel == WAVEMODEL_SURFBEAT) then
         ! 1) Bichromatic waves
         if (    (par%wbctype == WBCTYPE_PARAMS) .and. (par%Tlong /= -123)   ) then
            do j=1,s%ny+1
               s%ee(1,j,:)=e01*0.5d0 * &
               (1.d0+cos(2*par%px*(par%t/par%Tlong-( sin(s%theta0)*(s%yz(1,j)-s%yz(1,1)) & !jaap : -sum(s%alfaz(1,1:j))/j
               +cos(s%theta0)*(s%xz(1,j)-s%xz(1,1)) )/Llong))) * & !jaap: -sum(s%alfaz(1,1:j))/j
               min(par%t/par%taper,1.d0)
               em = (sum(0.5d0*e01))*s%dtheta *min(par%t/par%taper,1.d0)
               ei =  sum(s%ee(1,j,1:s%ntheta))*s%dtheta
               s%bi(1) = -(2*s%cg(1,j)/s%c(1,j)-0.5d0)*(em-ei)/(s%cg(1,j)**2-par%g*s%hh(1,j))/par%rho
               ht=s%zs0(1:2,:)-s%zb(1:2,:)
               s%ui(1,j) = s%cg(1,j)*s%bi(1)/ht(1,j)*cos(s%theta0-s%alfaz(1,j))
            end do
            ! 2) TS_1: bound long wave is calculate by the bound long wave theory of LG & Steward (1964)
         elseif (par%wbctype==WBCTYPE_TS_1) then
            do j=1,s%ny+1
               if (abs(s%theta0-sum(s%alfaz(1,1:j))/j)<1e-3) then
                  call linear_interp(tE,dataE,nt,par%t,E1,E_idx)
               else
                  tshifted = max(par%t-(s%yz(1,j)-s%yz(1,1))*sin(s%theta0)/s%cg(1,1) & !-sum(s%alfaz(1,1:j))/j
                  +(s%xz(1,j)-s%xz(1,1))*cos(s%theta0)/s%cg(1,1),0.d0) !-sum(s%alfaz(1,1:j))/j
                  call linear_interp(tE,dataE,nt,tshifted,E1,E_idx)
               endif
               s%ee(1,j,:)=e01*E1/max(Emean,0.000001d0)*min(par%t/par%taper,1.d0)
               em = Emean *min(par%t/par%taper,1.d0)
               ei = sum(s%ee(1,j,:))*s%dtheta
               s%bi(1) = -(2*s%cg(1,j)/s%c(1,j)-0.5d0)*(em-ei)/(s%cg(1,j)**2-par%g*s%hh(1,j))/par%rho
               ht=s%zs0(1:2,:)-s%zb(1:2,:)
               s%ui(1,j) = s%cg(1,j)*s%bi(1)/ht(1,j)*cos(s%theta0-s%alfaz(1,j))
            end do
            ! 3) TS_2: bound long wave is specificed by the user
         elseif (par%wbctype==WBCTYPE_TS_2) then
            ht=s%zs0(1:2,:)-s%zb(1:2,:)
            do j=1,s%ny+1
               if (abs(s%theta0-sum(s%alfaz(1,1:j))/j)<1e-3) then
                  call linear_interp(tE,dataE,nt,par%t,E1,E_idx)
                  call linear_interp(tE,databi,nt,par%t,s%bi(1),E_idx)
               else
                  tshifted = max(par%t-(s%yz(1,j)-s%yz(1,1))*sin(s%theta0)/s%cg(1,1) & !-sum(s%alfaz(1,1:j))/j
                  +(s%xz(1,j)-s%xz(1,1))*cos(s%theta0)/s%cg(1,1),0.d0) !-sum(s%alfaz(1,1:j))/j
                  call linear_interp(tE,dataE,nt,tshifted,E1,E_idx)
                  call linear_interp(tE,databi,nt,tshifted,s%bi(1),E_idx)
               endif
               s%ee(1,j,:)=e01*E1/max(Emean,0.000001d0)*min(par%t/par%taper,1.d0)
               if (par%freewave==1) then
                  s%ui(1,j) = sqrt(par%g/ht(1,j))*s%bi(1)
               else
                  s%ui(1,j) = s%cg(1,j)*s%bi(1)/ht(1,j)*cos(s%theta0-s%alfaz(1,j))*min(par%t/par%taper,1.d0)
               endif
            end do
            ! 4) go over the different wbctypes (spectral)
         else
            if ((par%wbctype==WBCTYPE_PARAMETRIC) .or. (par%wbctype==WBCTYPE_JONS_TABLE) .or. &
            (par%wbctype==WBCTYPE_SWAN) .or. (par%wbctype==WBCTYPE_VARDENS) .or. (par%wbctype==WBCTYPE_REUSE)) then
               ! open file if first time
               if (startbcf) then
                  if(xmaster) then
                     open(53,file='ebcflist.bcf',form='formatted',position='rewind')
                     open(54,file='qbcflist.bcf',form='formatted',position='rewind')
                     if(par%single_dir==1 .and. par%wbctype==WBCTYPE_REUSE) then
                        open(55,file='esbcflist.bcf',form='formatted',position='rewind')
                     endif
                  endif
                  if (xmaster) then
                     do i=1,curline
                        read(53,*,iostat=ier)bcendtime,rt,dtbcfile,par%Trep,s%theta0,ebcfname
                        if (ier .ne. 0) then
                           call report_file_read_error('ebcflist.bcf')
                        endif
                        read(54,*,iostat=ier)bcendtime,rt,dtbcfile,par%Trep,s%theta0,qbcfname
                        if (ier .ne. 0) then
                           call report_file_read_error('qbcflist.bcf')
                        endif
                        if(par%single_dir==1 .and. par%wbctype==WBCTYPE_REUSE) then
                           read(55,*,iostat=ier)bcendtime,rt,dtbcfile,par%Trep,s%theta0,esbcfname
                           if (ier .ne. 0) then
                              call report_file_read_error('esbcflist.bcf')
                           endif
                        endif
                     enddo  ! wwvv strange
                     if(par%single_dir==1 .and. par%wbctype==WBCTYPE_REUSE) then
                        inquire(iolength=wordsize) 1.d0
                        reclen=wordsize*(sg%ny+1)*(sg%ntheta_s)
                        open(73,file=esbcfname,status='old',form='unformatted',access='direct',recl=reclen)
                        read(73,rec=1,iostat=ier)sg%ee_s(1,:,:)
                        if (ier .ne. 0) then
                           call report_file_read_error(esbcfname)
                        endif
                        close(73)
                     endif
                  endif
#ifdef USEMPI
                  call xmpi_bcast(bcendtime)
                  call xmpi_bcast(rt)
                  call xmpi_bcast(dtbcfile)
                  call xmpi_bcast(par%Trep)
                  call xmpi_bcast(s%theta0)
                  call xmpi_bcast(ebcfname)
                  if (par%single_dir==1) then
                     call space_distribute("y",sl,sg%ee_s(1,:,:),s%ee_s(1,:,:))
                  endif
#else
                  if (par%single_dir==1) then
                     s%ee_s=sg%ee_s
                  endif
#endif
                  if (xmaster) then
                     close(53)
                     close(54)
                     if(par%single_dir==1 .and. par%wbctype==WBCTYPE_REUSE) then
                        close(55)
                     endif
                  endif
                  ! Robert and Jaap : Initialize for new wave conditions
                  do itheta=1,s%ntheta
                     s%sigt(:,:,itheta) = 2*par%px/par%Trep
                  end do
                  s%sigm = sum(s%sigt,3)/s%ntheta
                  call dispersion(par,s)
                  ! End initialize
                  if (xmaster) then
                     inquire(iolength=wordsize) 1.d0
                     reclen=wordsize*(sg%ny+1)*(sg%ntheta)
                     open(71,file=ebcfname,status='old',form='unformatted',access='direct',recl=reclen)
                     reclen=wordsize*((sg%ny+1)*4)
                     open(72,file=qbcfname,status='old',form='unformatted',access='direct',recl=reclen)
                  endif
                  ! wwvv note that we need the global value of ny here
                  ! masterprocess reads and distributes
                  if (xmaster) then
                     if (.not. allocated(gq1) ) then
                        allocate(gq1(sg%ny+1,4),gq2(sg%ny+1,4),gq(sg%ny+1,4))
                        allocate(gee1(sg%ny+1,s%ntheta),gee2(sg%ny+1,s%ntheta))
                     endif
                  else
                     if (.not. allocated(gq1) ) then ! to get valid addresses for
                        ! gq1, gq2, gq, gee1, gee2
                        allocate(gq1(1,4),gq2(1,4),gq(1,4))
                        allocate(gee1(1,s%ntheta),gee2(1,s%ntheta))
                     endif
                  endif
                  if (.not. allocated(q1) ) then
                     allocate(q1(s%ny+1,4),q2(s%ny+1,4),q(s%ny+1,4))
                     allocate(ee1(s%ny+1,s%ntheta),ee2(s%ny+1,s%ntheta))
                  end if
                  if (xmaster) then
                     read(71,rec=1,iostat=ier)gee1       ! Earlier in time
                     read(71,rec=2,iostat=ier2)gee2       ! Later in time
                     if (ier+ier2 .ne. 0) then
                        call report_file_read_error(ebcfname)
                     endif
                     read(72,rec=1,iostat=ier)gq1        ! Earlier in time
                     read(72,rec=2,iostat=ier2)gq2        ! Later in time
                     if (ier+ier2 .ne. 0) then
                        call report_file_read_error(qbcfname)
                     endif
                  endif
#ifdef USEMPI
                  call space_distribute("y",sl,gee1,ee1)
                  call space_distribute("y",sl,gee2,ee2)
                  call space_distribute("y",sl,gq1,q1)
                  call space_distribute("y",sl,gq2,q2)
#else
                  ee1=gee1
                  ee2=gee2
                  q1=gq1
                  q2=gq2
#endif
                  old=floor((par%t/dtbcfile)+1)
                  recpos=1
               end if
               new=floor((par%t/dtbcfile)+1)
               ! Check for next level in boundary condition file
               if (new/=old) then
                  recpos=recpos+(new-old)
                  ! Check for how many bcfile steps are jumped
                  if (new-old>1) then  ! Many steps further in the bc file
                     if(xmaster) then
                        read(72,rec=recpos+1,iostat=ier)gq2
                        if (ier .ne. 0) then
                           call report_file_read_error(qbcfname)
                        endif
                        read(71,rec=recpos+1,iostat=ier)gee2
                        if (ier .ne. 0) then
                           call report_file_read_error(ebcfname)
                        endif
                        read(72,rec=recpos,iostat=ier)gq1
                        if (ier .ne. 0) then
                           call report_file_read_error(qbcfname)
                        endif
                        read(71,rec=recpos,iostat=ier)gee1
                        if (ier .ne. 0) then
                           call report_file_read_error(ebcfname)
                        endif
                     endif
#ifdef USEMPI
                     call space_distribute("y",sl,gee2,ee2)
                     call space_distribute("y",sl,gq2,q2)
                     call space_distribute("y",sl,gee1,ee1)
                     call space_distribute("y",sl,gq1,q1)
#else
                     ee1=gee1
                     ee2=gee2
                     q1=gq2
                     q2=gq2
#endif
                  else  ! Only one step further in the bc file
                     ee1=ee2
                     q1=q2
                     if(xmaster) then
                        read(72,rec=recpos+1,iostat=ier)gq2
                        if (ier .ne. 0) then
                           call report_file_read_error(qbcfname)
                        endif
                        read(71,rec=recpos+1,iostat=ier)gee2
                        if (ier .ne. 0) then
                           call report_file_read_error(ebcfname)
                        endif
                     endif
#ifdef USEMPI
                     call space_distribute("y",sl,gee2,ee2)
                     call space_distribute("y",sl,gq2,q2)
#else
                     ee2=gee2
                     q2=gq2
#endif
                  endif
                  old=new
               end if
               ht=s%zs0(1:2,:)-s%zb(1:2,:)
               tnew = dble(new)*dtbcfile
               if (xmpi_istop) s%ee(1,:,:) = (dtbcfile-(tnew-par%t))/dtbcfile*ee2 + & !Jaap
               (tnew-par%t)/dtbcfile*ee1
               q = (dtbcfile-(tnew-par%t))/dtbcfile*q2 + &          !Jaap
               (tnew-par%t)/dtbcfile*q1
               ! be aware: ui and vi are defined w.r.t. the grid, not w.r.t. the coordinate system!
               s%ui(1,:) = (q(:,1)*dcos(s%alfau(1,:)) + q(:,2)*dsin(s%alfau(1,:)))/ht(1,:)*min(par%t/par%taper,1.0d0)
               s%vi(1,:) = (-q(:,1)*dsin(s%alfau(1,:)) + q(:,2)*dcos(s%alfau(1,:)))/ht(1,:)*min(par%t/par%taper,1.0d0)
               s%ee(1,:,:)=s%ee(1,:,:)*min(par%t/par%taper,1.0d0)
            else
            endif
         endif ! Go over the different wbctypes for surfbeat
         ! Part C) Nonh: go over the different waveforms
      elseif (par%wavemodel == WAVEMODEL_NONH) then
         ! 1) spectral
         if ((par%wbctype==WBCTYPE_PARAMETRIC) .or. (par%wbctype==WBCTYPE_JONS_TABLE) .or. &
         (par%wbctype==WBCTYPE_SWAN) .or. (par%wbctype==WBCTYPE_VARDENS) .or. (par%wbctype==WBCTYPE_REUSE)) then
            if (startbcf) then
               if(xmaster) then
                  open(53,file='nhbcflist.bcf',form='formatted',position='rewind')
                  do i=1,curline
                     read(53,*,iostat=ier)bcstarttime,bcendtime,par%Trep,nhbcfname
                     if (ier .ne. 0) then
                        call report_file_read_error('nhbcflist.bcf')
                     endif
                  enddo
                  close(53)
               endif
               if (.not. allocated(uig)) then
                  if (xmaster) then
                     allocate(uig(sg%ny+1))
                     allocate(vig(sg%ny+1))
                     allocate(zig(sg%ny+1))
                     allocate(wig(sg%ny+1))
                     allocate(duig(sg%ny+1))
                     allocate(dvig(sg%ny+1))
                  else  ! only need valid address for MPI-distribution
                     allocate(uig(1))
                     allocate(vig(1))
                     allocate(zig(1))
                     allocate(wig(1))
                     allocate(duig(1))
                     allocate(dvig(1))
                  endif
               endif
               if (xmaster) then
                  ! Robert: optimise?
                  if (trim(nhbcfname)=='nh_reuse.bcf') then
                     ! Reuse time series will start at t=0, so we need to add the offset time to
                     ! the time vector in the boundary condition file using the optional bcst variable
                     call velocity_Boundary(nhbcfname,uig,vig,zig,wig,duig,dvig,sg%ny,par%t, &
                     isSet_U,isSet_V,isSet_Z,isSet_W,isSet_dU,isSet_dV,isSet_Q, &
                     force_init=.true.,bcst=bcstarttime)
                  else
                     ! If individual spectra are read, then the time vector in the boundary condition
                     ! file will already be corrected for the start of the spectrum condition. No need
                     ! to add the optional bcst variable
                     call velocity_Boundary(nhbcfname,uig,vig,zig,wig,duig,dvig,sg%ny,par%t, &
                     isSet_U,isSet_V,isSet_Z,isSet_W,isSet_dU,isSet_dV,isSet_Q, &
                     force_init=.true.)
                  endif
               endif
#ifdef USEMPI
               call space_distribute("y",sl,uig,s%ui(1,:))
               call space_distribute("y",sl,vig,s%ui(1,:))
               call space_distribute("y",sl,zig,s%zi(1,:))
               call space_distribute("y",sl,wig,s%wi(1,:))
               call space_distribute("y",sl,duig,s%dui(1,:))
               call xmpi_bcast(isSet_U)
               call xmpi_bcast(isSet_V)
               call xmpi_bcast(isSet_Z)
               call xmpi_bcast(isSet_W)
               call xmpi_bcast(isSet_dU)
               call xmpi_bcast(bcendtime)
#else
               s%ui(1,:) = uig
               s%vi(1,:) = vig
               s%zi(1,:) = zig
               s%wi(1,:) = wig
               s%dui(1,:) = duig
#endif
               if (.not. isSet_U) s%ui(1,:) = 0.d0
               if (.not. isSet_V) s%vi(1,:) = 0.d0
               if (.not. isSet_Z) s%zi(1,:) = s%zs(2,:)
               if (.not. isSet_W) s%wi(1,:) = s%ws(2,:)
               if (par%nonhq3d == 1 .and. par%nhlay > 0.0d0) then
                  ! Only if the reduced 2-layer model is initiated (P.B. Smit)
                  if (.not. isSet_dU) s%dui(1,:) = 0.
               endif
               call nonh_init_wcoef(s,par)
            else
               if (xmaster) then
                  if (trim(nhbcfname)=='nh_reuse.bcf') then
                     ! Reuse time series will start at t=0, so we need to add the offset time to
                     ! the time vector in the boundary condition file using the optional bcst variable
                     call velocity_Boundary(nhbcfname,uig,vig,zig,wig,duig,dvig,sg%ny,par%t,&
                     isSet_U,isSet_V,isSet_Z,isSet_W,isSet_dU,isSet_dV,isSet_Q, &
                     bcst=bcstarttime)
                  else
                     ! If individual spectra are read, then the time vector in the boundary condition
                     ! file will already be corrected for the start of the spectrum condition. No need
                     ! to add the optional bcst variable
                     call velocity_Boundary(nhbcfname,uig,vig,zig,wig,duig,dvig,sg%ny,par%t, &
                     isSet_U,isSet_V,isSet_Z,isSet_W,isSet_dU,isSet_dV,isSet_Q)
                  endif
               endif
#ifdef USEMPI
               call space_distribute("y",sl,uig,s%ui(1,:))
               call space_distribute("y",sl,vig,s%vi(1,:))
               call space_distribute("y",sl,zig,s%zi(1,:))
               call space_distribute("y",sl,wig,s%wi(1,:))
               call space_distribute("y",sl,duig,s%dui(1,:))
               call xmpi_bcast(isSet_U)
               call xmpi_bcast(isSet_V)
               call xmpi_bcast(isSet_Z)
               call xmpi_bcast(isSet_W)
               call xmpi_bcast(isSet_dU)
#else
               s%dui(1,:) = duig
               s%ui(1,:) = uig
               s%vi(1,:) = vig
               s%zi(1,:) = zig
               s%wi(1,:) = wig
#endif
               if (.not. isSet_U) s%ui(1,:) = 0.d0
               if (.not. isSet_V) s%vi(1,:) = 0.d0
               if (.not. isSet_Z) s%zi(1,:) = s%zs(2,:)
               if (.not. isSet_W) s%wi(1,:) = s%ws(2,:)
               if (par%nonhq3d == 1 .and. par%nhlay > 0.0d0) then
                  ! Only if the reduced 2-layer model is initiated (P.B. Smit)
                  if (.not. isSet_dU) s%dUi(1,:) = 0.
                  !
               endif
            endif
            ! be aware: ui and vi are defined w.r.t. the grid, not w.r.t. the coordinate system!
            if (.not. allocated(tempu)) then
               ! allocate to local size grid
               allocate(tempu(s%ny+1))
               allocate(tempv(s%ny+1))
            endif
            tempu = s%ui(1,:)
            tempv = s%vi(1,:)
            ! rotate to grid directions
            s%ui(1,:) = tempu*dcos(s%alfau(1,:)) + tempv*dsin(s%alfau(1,:))
            s%vi(1,:) = -tempu*dsin(s%alfau(1,:)) + tempv*dcos(s%alfau(1,:))
            s%ui(1,:) = s%ui(1,:)*min(par%t/par%taper,1.0d0)
            s%zi(1,:) = s%zi(1,:)*min(par%t/par%taper,1.0d0)
            s%wi(1,:) = s%wi(1,:)*min(par%t/par%taper,1.0d0)
            s%dUi(1,:) = s%dUi(1,:)*min(par%t/par%taper,1.0d0)
         elseif (par%wbctype==WBCTYPE_TS_NONH) then
            if (xmaster) then
               call velocity_Boundary('boun_U.bcf',uig,vig,zig,wig,duig,dvig,sg%ny,&
               par%t,isSet_U,isSet_V,isSet_Z,isSet_W,isSet_dU,isSet_dV,isSet_Q)
            endif
#ifdef USEMPI
            call space_distribute("y",sl,uig,s%ui(1,:))
            call space_distribute("y",sl,vig,s%vi(1,:))
            call space_distribute("y",sl,zig,s%zi(1,:))
            call space_distribute("y",sl,wig,s%wi(1,:))
            call space_distribute("y",sl,duig,s%dui(1,:))
            call xmpi_bcast(isSet_U)
            call xmpi_bcast(isSet_Z)
            call xmpi_bcast(isSet_W)
            call xmpi_bcast(isSet_dU)
#else
            if ( isSet_Q ) then
               if (par%nonhq3d == 1 .and. par%nhlay > 0.0d0) then
                  ! Calculate the difference velocity from the difference discharges (here both duig and uig are still discharges)
                  duig = ( (1.0d0 - 2.0d0*par%nhlay) * uig + duig ) / (  2.0d0 * par%nhlay * s%hu(1,:) *( 1.0d0 - par%nhlay )  )
               endif
               uig = uig/s%hu(1,:)
            endif
            s%ui(1,:) = uig
            s%vi(1,:) = vig
            s%zi(1,:) = zig
            s%wi(1,:) = wig
            s%dui(1,:) = duig
#endif
            if (.not. isSet_U) s%ui(1,:) = 0.d0
            if (.not. isSet_V) s%vi(1,:) = 0.d0
            if (.not. isSet_Z) s%zi(1,:) = s%zs(2,:)
            if (.not. isSet_W) s%wi(1,:) = s%ws(2,:)
            if (par%nonhq3d == 1 .and. par%nhlay > 0.0d0) then
               ! Only if the reduced 2-layer model is initiated (P.B. Smit)
               if (.not. isSet_dU) s%dui(1,:) = 0.
            endif
         else
         endif ! Go over the different wbctypes for nonh
      endif ! Go over the different wavemodels
      !
      ! Finalize: swave and lwave = 0 has influence on boundary
      s%ee(1,:,:) = par%swave*s%ee(1,:,:)
      if (par%nonhspectrum==0) then
         s%ui = par%lwave*(par%order-1.d0)*s%ui
      endif
   end subroutine wave_bc
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! FLOW BOUNDARY CONDITIONS
   !
   subroutine flow_bc(s,par)
      use params
      use spaceparams
      use interp
      use xmpi_module
      use paramsconst
      use compute_tide_module

      implicit none

      type(spacepars), target                     :: s
      type(parameters)                            :: par

      integer                                     :: i,j,jj,j1,indt
      real*8                                      :: alphanew,vert,windxnow,windynow,factime
      real*8                                      :: udnusum,dnusum,vdnvsum,dnvsum
      !real*8 , dimension(2)                       :: xzs0,yzs0,szs0
      real*8 , dimension(:,:)  ,allocatable,save  :: zs0old,dzs0
      real*8 , dimension(:,:)  ,allocatable,save  :: zsmean
      real*8 , dimension(:,:)  ,allocatable,save  :: ht,beta,betanp1
      real*8 , dimension(:)    ,allocatable,save  :: bn,alpha2,thetai
      real*8 , dimension(:,:)  ,allocatable,save  :: dhdx,dhdy,dvdx,dvdy
      real*8 , dimension(:,:)  ,allocatable,save  :: dbetadx,dbetady,dvudy
      real*8 , dimension(:)    ,allocatable,save  :: inv_ht


      if (.not. allocated(ht)) then
         allocate(ht      (2,s%ny+1))
         allocate(dhdx    (2,s%ny+1))
         allocate(dhdy    (2,s%ny+1))   ! wwvv not used
         allocate(dvdx    (2,s%ny+1))   ! wwvv not used
         allocate(dvdy    (2,s%ny+1))
         allocate(dvudy   (1,s%ny+1))   ! wwvv not used
         allocate(inv_ht  (s%ny+1))
         allocate(dbetadx (2,s%ny+1))
         allocate(dbetady (2,s%ny+1))
         allocate(beta    (2,s%ny+1))
         allocate(zsmean  (2,s%ny+1))
         allocate(bn      (s%ny+1))
         allocate(alpha2  (s%ny+1))
         allocate(thetai  (s%ny+1))
         allocate(betanp1 (1,s%ny+1))
         if (par%tidetype==TIDETYPE_INSTANT) then
            allocate(zs0old(s%nx+1,s%ny+1))
            allocate(dzs0  (s%nx+1,s%ny+1))
            zs0old = s%zs0
            dzs0 = 0.d0
         endif
         ! initialize zsmean
         zsmean(1,:) = s%zs(1,:)
         zsmean(2,:) = s%zs(s%nx,:)
         s%umean = s%uu
         s%vmean = s%vv
         thetai = 0.d0
      endif

      ! Super fast 1D
      if (s%ny==0) then
         j1 = 1
      else
         j1 = 2
      endif

      ! factime=1.d0/par%cats/par%Trep*par%dt !Jaap: not used anymore
      ! jaap: compute epsi based on offshore boundary conditions if epsi = -1
      if (par%epsi==-1.d0) then
         if (par%swave==1) then
            factime = 1.d0/par%cats/par%Trep*par%dt
         else
            ! need to do something here
            factime = 1.d0/60*par%dt
         endif
      else
         factime = par%epsi
      endif

      !allocate(szs0(1:2))  ! wwvv changed this, now defined as an automatic array
      !allocate(yzs0(1:2))
      !
      ! Sea boundary at x=0;
      !
      !
      ! UPDATE TIDE AND SURGE
      if (par%tidetype==TIDETYPE_INSTANT) then
         zs0old = s%zs0  ! store zs0 of old time step
         call tide_boundary_timestep(s,par) ! also updates zs0 across domain
         where (s%wetz==1)
            dzs0 = s%zs0-zs0old
            s%zs = s%zs+dzs0
         elsewhere
            dzs0 = 0.d0
         endwhere
      else
         call tide_boundary_timestep(s,par) ! only updates zs0 on boundaries
      endif
      !      !
      !      ! Need to interpolate input tidal signal to xbeach par%t to
      !      ! compute proper tide contribution
      !
      !      if (par%tideloc>0) then
      !
      !         ! read in first water surface time series
      !         call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,1), s%tidelen, par%t, s%zs01, indt)
      !
      !         if(par%tideloc.eq.1) then
      !            s%zs02=s%zs01
      !         end if
      !
      !         ! tideloc = 2, paulrevere = land
      !         if(par%tideloc.eq.2 .and. par%paulrevere==PAULREVERE_LAND) then
      !            ! read in second water surface time series
      !            call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,2), s%tidelen, par%t, s%zs03, indt)
      !            s%zs02=s%zs01 ! second offshore corner is equal to first offshore corner
      !            s%zs04=s%zs03 ! second bay corner is equal to first bay corner
      !         endif
      !
      !         ! tideloc = 2, paulrevere = sea
      !         if(par%tideloc.eq.2 .and. par%paulrevere==PAULREVERE_SEA) then
      !            call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,2), s%tidelen, par%t, s%zs02, indt)
      !            ! no timeseries at bay side, (and two different timeseries at offshore corners)
      !            s%zs03=0.d0
      !            s%zs04=0.d0
      !         endif
      !
      !         ! tideloc = 4: for each corner individual timeseries
      !         if(par%tideloc.eq.4) then
      !            call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,2), s%tidelen, par%t, s%zs02, indt)
      !            call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,3), s%tidelen, par%t, s%zs03, indt)
      !            call LINEAR_INTERP(s%tideinpt, s%tideinpz(:,4), s%tidelen, par%t, s%zs04, indt)
      !         endif
      !         !
      !         ! from here on set global variable zs0
      !         !
      !         if(par%tideloc.eq.1) s%zs0 = s%zs01
      !
      !         if(par%tideloc.eq.2 .and. par%paulrevere==PAULREVERE_SEA) then
      !            yzs0(1)=s%ndist(1,1)
      !            yzs0(2)=s%ndist(1,s%ny+1)
      !
      !            ! for MPI look at water level gradient over each subdomain
      !
      !            ! lonsghore water level difference over whole model domain at offshore boundary
      !            dzs0dy = (s%zs02-s%zs01)/(s%xyzs02(2)-s%xyzs01(2));
      !
      !            ! estimate water level at corners sub-domain
      !            szs0(1)=s%zs01+dzs0dy*(yzs0(1)-s%xyzs01(2))
      !            szs0(2)=s%zs01+dzs0dy*(yzs0(2)-s%xyzs01(2))
      !
      !            do j = 1,s%ny+1
      !               call LINEAR_INTERP(yzs0, szs0, 2, s%ndist(1,j), s%zs0(1,j), indt)
      !            enddo
      !
      !            s%zs0(2,:)=s%zs0(1,:)
      !            ! Dano: fix of Neumann boundaries with tide
      !            s%zs0(3:s%nx+1,1)=s%zs0(1,1)
      !            s%zs0(3:s%nx+1,s%ny+1)=s%zs0(1,s%ny+1)
      !            !do j = 1,ny+1
      !            !   do i = 1,nx+1
      !            !      zs0(i,j) = zs0(1,j)
      !            !   enddo
      !            !enddo
      !         elseif (par%tideloc.eq.2 .and. par%paulrevere==PAULREVERE_LAND) then
      !            if (xmpi_istop) then
      !               s%zs0(1:2,:)=s%zs01
      !            endif
      !            if (xmpi_isbot) then
      !               s%zs0(s%nx:s%nx+1,:)=s%zs04
      !            endif
      !         endif
      !
      !         if(par%tideloc.eq.4) then
      !            if (xmpi_istop) then
      !               yzs0(1)=s%xyzs01(2)
      !               yzs0(2)=s%xyzs02(2)
      !               szs0(1)=s%zs01
      !               szs0(2)=s%zs02
      !               do j = 1,s%ny+1
      !                  call LINEAR_INTERP(yzs0, szs0, 2, s%ndist(1,j), s%zs0(1,j), indt)
      !               enddo
      !               s%zs0(2,:)=s%zs0(1,:)
      !            endif
      !            if (xmpi_isbot) then
      !               yzs0(1)=s%xyzs04(2)
      !               yzs0(2)=s%xyzs03(2)
      !               szs0(1)=s%zs04
      !               szs0(2)=s%zs03
      !               do j = 1,s%ny+1
      !                  call LINEAR_INTERP(yzs0, szs0, 2, s%ndist(s%nx+1,j), s%zs0(s%nx+1,j), indt)
      !               enddo
      !               s%zs0(s%nx,:)=s%zs0(s%nx+1,:)
      !            endif
      !            if (xmpi_isleft) then
      !               xzs0(1)=s%xyzs01(1)
      !               xzs0(2)=s%xyzs04(1)
      !               szs0(1)=s%zs01
      !               szs0(2)=s%zs04
      !               do i = 1,s%nx+1
      !                  call LINEAR_INTERP(xzs0, szs0, 2, s%sdist(i,1), s%zs0(i,1), indt)
      !               enddo
      !            endif
      !            if (xmpi_isright) then
      !               xzs0(1)=s%xyzs02(1)
      !               xzs0(2)=s%xyzs03(1)
      !               szs0(1)=s%zs02
      !               szs0(2)=s%zs03
      !               do i = 1,s%nx+1
      !                  call LINEAR_INTERP(xzs0, szs0, 2, s%sdist(i,s%ny+1), s%zs0(i,s%ny+1), indt)
      !               enddo
      !            endif
      !         endif
      !
      !      else ! ie if tideloc=0
      !         s%zs0 = s%zs01
      !      endif
      !
      !      if (par%tidetype==TIDETYPE_INSTANT) then
      !         !
      !         ! RJ: 22-09-2010 Correct water level for surge and tide:
      !         !
      !         zsmean(1,:) = factime*s%zs(1,:)+(1-factime)*zsmean(1,:)
      !         zsmean(2,:) = factime*s%zs(s%nx,:)+(1-factime)*zsmean(2,:)
      !         ! compute difference between offshore/bay mean water level and imposed on tide
      !         dzs0(1,:) = s%zs0(1,:)-zsmean(1,:)
      !         dzs0(2,:) = s%zs0(s%nx,:)-zsmean(2,:)
      !#ifdef USEMPI
      !         call xmpi_getrow(dzs0,s%ny+1,'1',1,dzs0(1,:))
      !         call xmpi_getrow(dzs0,s%ny+1,'m',xmpi_m,dzs0(2,:))
      !#endif
      !
      !         do j = 1,s%ny+1
      !            do i = 1,s%nx+1
      !               s%zs(i,j) = s%zs(i,j) + (s%zs0fac(i,j,1)*dzs0(1,j) + s%zs0fac(i,j,2)*dzs0(2,j))*s%wetz(i,j)
      !            enddo
      !         enddo
      !
      !         ! RJ: 22-09-2010 end update tide and surge
      !
      !      endif ! tidetype = instant water level boundary

      !Dano: compute umean, vmean only at this location, same for all options; MPI compliant
      if (par%tidetype==TIDETYPE_VELOCITY) then
         if (xmpi_istop) then
            ! make sure we have same umean along whole offshore boundary
            udnusum   = sum(s%uu (1,jmin_uu:jmax_uu)*s%dnu(1,jmin_uu:jmax_uu))
            dnusum    = sum(s%dnu(1,jmin_uu:jmax_uu))
            vdnvsum   = sum(s%vv (1,jmin_vv:jmax_vv)*s%dnv(1,jmin_vv:jmax_vv))
            dnvsum    = sum(s%dnv(1,jmin_vv:jmax_vv))
         else
            udnusum   = 0.d0
            dnusum    = 0.d0
            vdnvsum   = 0.d0
            dnvsum    = 0.d0
         endif
#ifdef USEMPI
         call xmpi_allreduce(udnusum  ,MPI_SUM)
         call xmpi_allreduce(dnusum   ,MPI_SUM)
         call xmpi_allreduce(vdnvsum  ,MPI_SUM)
         call xmpi_allreduce(dnvsum   ,MPI_SUM)
#endif
         s%umean(1,:)=factime*udnusum/dnusum+(1.0d0-factime)*s%umean(1,:)
         s%vmean(1,:)=factime*vdnvsum/dnvsum+(1.0d0-factime)*s%vmean(1,:)

      else
         s%umean(1,:) = 0.d0
         s%vmean(1,:) = 0.d0
      endif
      !
      ! UPDATE (LONG) WAVES
      !
      if (par%wbctype/=WBCTYPE_OFF)then
         ! wwvv the following is probably only to do in the top processes, but take care for
         ! the mpi_shift calls in horizontal directions
         if(xmpi_istop) then
            if (par%front==FRONT_ABS_1D) then ! Ad's radiating boundary
               !ht(1:2,:)=max(zs0(1:2,:)-zb(1:2,:),par%eps)
               !uu(1,:)=2.d0*ui(1,:)-(sqrt(par%g/hh(1,:))*(zs(2,:)-zs0(2,:)))
               if (par%freewave==1) then
                  s%uu(1,:)=2.0d0*s%ui(1,:)-(sqrt(par%g/s%hh(1,:))*(s%zs(2,:)-s%zs0(2,:))) + s%umean(1,:)
               else
                  s%uu(1,:)=(1.0d0+sqrt(par%g*s%hh(1,:))/s%cg(1,:))*s%ui(1,:) - &
                  (sqrt(par%g/s%hh(1,:))*(s%zs(2,:)-s%zs0(2,:))) + s%umean(1,:)
               endif
               s%vv(1,:)=s%vv(2,:)
               s%zs(1,:)=s%zs(2,:)
               if (par%wavemodel==WAVEMODEL_NONH) s%ws(1,:) = s%ws(2,:)
            elseif (par%front==FRONT_WAVEFLUME) then ! Ad's radiating boundary
               if (par%tidetype==TIDETYPE_VELOCITY) then
                  zsmean(1,:) = factime*s%zs(1,:)+(1-factime)*zsmean(1,:)
               endif
               s%uu(1,:)=(1.0d0+sqrt(par%g*s%hh(1,:))/s%cg(1,:))*s%ui(1,:)-(sqrt(par%g/s%hh(1,:))*(s%zs(2,:)-zsmean(1,:)))
               s%vv(1,:)=s%vv(2,:)
               s%zs(1,:)=s%zs(2,:)
               if (par%wavemodel==WAVEMODEL_NONH) s%ws(1,:) = s%ws(2,:)
            elseif (par%front==FRONT_ABS_2D) then ! Van Dongeren (1997), weakly reflective boundary condition
               ht(1:2,:)=max(s%zs0(1:2,:)-s%zb(1:2,:),par%eps)
               beta=s%uu(1:2,:)-2.*dsqrt(par%g*s%hum(1:2,:))

               do j=j1,max(s%ny,1)
                  ! compute gradients in u-points
                  dhdx   (1,j) = ( ht  (2,j) - ht  (1,j) ) / s%dsu(1,j)
                  dbetadx(1,j) = ( beta(2,j) - beta(1,j) ) / s%dsz(2,j)

                  if (s%ny>0) then
                     dvdy   (1,j) = ( s%vu  (1,j+1) - s%vu  (1,j-1) ) / ( s%dnu(1,j)+(s%dnu(1,j-1)+s%dnu(1,j+1))/2 )
                     dbetady(1,j) = ( beta(1,j+1) - beta(1,j-1) ) / ( s%dnu(1,j)+(s%dnu(1,j-1)+s%dnu(1,j+1))/2 )
                  else
                     dvdy   (1,j) = 0.d0
                     dbetady(1,j) = 0.d0
                  endif

                  inv_ht(j) = 1.d0/s%hum(1,j)

                  bn    (j) = -( s%uu(1,j) - dsqrt(par%g*s%hum(1,j)) ) * dbetadx(1,j) &
                  -  s%vu(1,j)                           * dbetady(1,j) &
                  +            dsqrt(par%g*s%hum(1,j))   * dvdy   (1,j) &
                  +  par%g                             * dhdx   (1,j) &
                  +  s%Fx(1,j) * inv_ht(j) / par%rho                    &
                  - s%cfu(1,j) * sqrt(s%uu(1,j)**2 + s%vu(1,j)**2) * s%uu(1,j) / s%hum(1,j)
               end do

               ! Jaap: Compute angle of incominge wave
               thetai = datan(s%vi(1,:)/(s%ui(1,:)+1.d-16))

               do j=j1,max(s%ny,1)

                  betanp1(1,j) = beta(1,j)+ bn(j)*par%dt
                  alpha2(j)=-(s%theta0-s%alfaz(1,j)) ! Jaap: this is first estimate
                  alphanew = 0.d0

                  do jj=1,50
                     !---------- Lower order bound. cond. ---

                     if (par%freewave==1) then ! assuming incoming long wave propagates at sqrt(g*h)
                        s%ur(1,j) = dcos(alpha2(j))/(dcos(alpha2(j))+1.d0)&
                        *(betanp1(1,j)-s%umean(1,j)+2.d0*DSQRT(par%g*0.5d0*(ht(1,j)+ht(2,j))) &
                        -s%ui(1,j)*(dcos(thetai(j))-1.d0)/dcos(thetai(j)))
                     else                     ! assuming incoming long wave propagates at s%cg
                        s%ur(1,j) = dcos(alpha2(j))/(dcos(alpha2(j))+1.d0)&
                        *(betanp1(1,j)-s%umean(1,j)+2.d0*DSQRT(par%g*0.5d0*(ht(1,j)+ht(2,j)))&
                        -s%ui(1,j)*(s%cg(1,j)*dcos(thetai(j))- DSQRT(par%g*0.5d0*(ht(1,j)+ht(2,j))) )/&
                        (s%cg(1,j)*dcos(thetai(j))))
                     endif

                     !vert = velocity of the reflected wave = total-specified
                     vert = s%vu(1,j)-s%vmean(1,j)-s%vi(1,j)                    ! Jaap use s%vi instead of s%ui(1,j)*tan(s%theta0)
                     alphanew = datan(vert/(s%ur(1,j)+1.d-16))                   ! wwvv can  use atan2 here
                     if (alphanew .gt. (par%px*0.5d0)) alphanew=alphanew-par%px
                     if (alphanew .le. (-par%px*0.5d0)) alphanew=alphanew+par%px
                     if(dabs(alphanew-alpha2(j)).lt.0.001d0) goto 1000    ! wwvv can use exit here
                     alpha2(j) = alphanew
                  end do
1000              continue
                  if (par%ARC==0) then
                     s%uu(1,j) = (par%order-1)*s%ui(1,j) + s%umean(1,j)
                     s%zs(1,j) = s%zs(2,j)
                  else
                     !
                     s%uu(1,j) = (par%order-1.d0)*s%ui(1,j) + s%ur(1,j) + s%umean(1,j)
                     ! zs(1,:) is dummy variable
                     ! with a taylor expansion to get to the zs point at index 1 from uu(1) and uu(2)
                     s%zs(1,j) = 1.5d0*((betanp1(1,j)-s%uu(1,j))**2/4.d0/par%g+.5d0*(s%zb(1,j)+s%zb(2,j)))- &
                     0.5d0*((beta(2,j)   -s%uu(2,j))**2/4.d0/par%g+.5d0*(s%zb(2,j)+s%zb(3,j)))
                     ! Ad + Jaap: zs does indeed influence hydrodynamics at boundary --> do higher order taylor expansions to check influence
                     ! zs(1,j) = 13.d0/8.d0*((betanp1(1,j)-uu(1,j))**2.d0/4.d0/par%g+.5d0*(zb(1,j)+zb(2,j))) - &
                     !           0.75d0*((beta(2,j)-uu(2,j))**2.d0/4.d0/par%g+.5d0*(zb(2,j)+zb(3,j)))        + &
                     !           0.125d0*0.5d0*(zs(3,j)+zs(4,j))
                  end if
               end do
               if (xmpi_istop) s%vv(1,:)=s%vv(2,:)
               if (xmpi_isleft) s%uu(1,1)=s%uu(1,2)
               if (xmpi_isright) s%uu(1,s%ny+1)=s%uu(1,s%ny)
               if (par%wavemodel==WAVEMODEL_NONH.and.xmpi_istop) s%ws(1,:) = s%ws(2,:)
            else if (par%front==FRONT_WALL) then
               !       uu(1,:)=0.d0
               !      zs(1,:)=max(zs(2,:),zb(1,:))
            else if (par%front==FRONT_WLEVEL) then
               s%zs(1,:)=s%zs0(1,:)
            else if (par%front==FRONT_NONH_1D) then
               !Timeseries Boundary for nonh, only use in combination with instat=8
               if (par%ARC==0) then
                  s%uu(1,:) = s%ui(1,:)
                  s%zs(1,:) = s%zi(1,:)+s%zs0(1,:)
                  s%ws(1,:) = s%wi(1,:)
                  if (par%nonhq3d == 1 .and. par%nhlay > 0.0d0) then
                     ! Only if the reduced 2-layer model is initiated (P.B. Smit)
                     s%dU(1,:) = s%dui(1,:)
                     !
                  endif
                  !
               else
                  !Radiating boundary for short waves:
                  s%uu(1,:) = s%ui(1,:)-sqrt(par%g/s%hh(1,:))*(s%zs(2,:)-s%zi(1,:)-s%zs0(2,:)) + s%umean(1,:)
                  !
                  !                uu(1,:)=  2*ui(1,:)-sqrt(par%g/hh(1,:))*(zs(2,:)        -zs0(2,:))
                  if (par%nonhq3d == 1 .and. par%nhlay > 0.0d0) then
                     !

                     s%dU(1,:) = s%dui(1,:)
                     ! write(*,*) s%dU(1,%)
                     !
                  endif

                  s%zs(1,:) = s%zs(2,:)
                  s%ws(1,:) = s%ws(2,:)
               endif
            endif ! par%front

         endif  ! xmpi_istop
         !
         ! Radiating boundary at x=nx*dx
         !
         !Dano: compute umean, vmean only at this location, same for all options; MPI compliant
         if (par%tidetype==TIDETYPE_VELOCITY) then
            ! make sure we have same umean along whole offshore boundary
            if (xmpi_isbot) then
               udnusum   = sum(s%uu (s%nx,jmin_uu:jmax_uu)*s%dnu(s%nx,jmin_uu:jmax_uu))
               dnusum    = sum(s%dnu(s%nx,jmin_uu:jmax_uu))
               vdnvsum   = sum(s%vv (s%nx,jmin_vv:jmax_vv)*s%dnv(s%nx,jmin_vv:jmax_vv))
               dnvsum    = sum(s%dnv(s%nx,jmin_vv:jmax_vv))
            else
               udnusum   = 0.d0
               dnusum    = 0.d0
               vdnvsum   = 0.d0
               dnvsum    = 0.d0
            endif
#ifdef USEMPI
            call xmpi_allreduce(udnusum  ,MPI_SUM)
            call xmpi_allreduce(dnusum   ,MPI_SUM)
            call xmpi_allreduce(vdnvsum  ,MPI_SUM)
            call xmpi_allreduce(dnvsum   ,MPI_SUM)
#endif
            s%umean(s%nx,:)=factime*udnusum/dnusum+(1.0d0-factime)*s%umean(s%nx,:)
            s%vmean(s%nx,:)=factime*vdnvsum/dnvsum+(1.0d0-factime)*s%vmean(s%nx,:)

         else
            s%umean(s%nx,:) = 0.d0
            s%vmean(s%nx,:) = 0.d0
         endif

         if (par%back==BACK_WALL) then
            if (xmpi_isbot) then
               s%uu(s%nx,:) = 0.d0
               s%zs(s%nx+1,:) = s%zs(s%nx,:)
            endif
            ! zs(nx+1,2:ny) = zs(nx+1,2:ny) + par%dt*hh(nx,2:ny)*uu(nx,2:ny)/(xu(nx+1)-xu(nx)) -par%dt*(hv(nx+1,2:ny+1)*vv(nx+1,2:ny+1)-hv(nx+1,1:ny)*vv(nx+1,1:ny))/(yv(2:ny+1)-yv(1:ny))
         elseif (par%back==BACK_ABS_1D) then
            !        umean(nx,:) = factime*uu(nx,:)+(1.d0-factime)*umean(nx,:)
            ! After hack 3/6/2010, return to par%epsi :
            s%uu(s%nx,:)=sqrt(par%g/s%hh(s%nx,:))*(s%zs(s%nx,:)-max(s%zb(s%nx,:),s%zs0(s%nx,:)))+s%umean(s%nx,:) ! cjaap: make sure if the last cell is dry no radiating flow is computed...
            !uu(nx,:)=cg(nx,:)/hh(nx,:)*(zs(nx,:)-max(zb(nx,:),zs0(nx,:)))+umean(nx,:)   ! cjaap: make sure if the last cell is dry no radiating flow is computed...
            !uu(nx,:)=sqrt(par%g/(zs0(nx,:)-zb(nx,:)))*(zs(nx,:)-max(zb(nx,:),zs0(nx,:)))
            !umean(nx,:) = factime*uu(nx,:)+(1-factime)*umean(nx,:)    !Ap
            !zs(nx+1,:)=max(zs0(nx+1,:),zb(nx+1,:))+(uu(nx,:)-umean(nx,:))*sqrt(max((zs0(nx+1,:)-zb(nx+1,:)),par%eps)/par%g)    !Ap
            s%zs(s%nx+1,:)=s%zs(s%nx,:)
         elseif (par%back==BACK_ABS_2D) then
            ht(1:2,:)=max(s%zs0(s%nx:s%nx+1,:)-s%zb(s%nx:s%nx+1,:),par%eps) !cjaap; make sure ht is always larger than zero

            beta=s%uu(s%nx-1:s%nx,:)+2.*dsqrt(par%g*s%hum(s%nx-1:s%nx,:)) !cjaap : replace s%hh with s%hum

            do j=j1,max(s%ny,1)
               if (s%wetu(s%nx,j)==1) then   ! Robert: dry back boundary points
                  ! Compute gradients in u-points....
                  dvdy(2,j)=(s%vu(s%nx,j+1)-s%vu(s%nx,j-1))/(2.d0*s%dnu(s%nx,j))
                  dhdx(2,j)=(ht(2,j)-ht(1,j))/s%dsu(s%nx,j)
                  dbetadx(2,j)=(beta(2,j)-beta(1,j))/s%dsz(s%nx,j)
                  dbetady(2,j)=(beta(1,j+1)-beta(1,j-1))/(2.0*s%dnu(s%nx,j))

                  inv_ht(j) = 1.d0/s%hum(s%nx,j)

                  bn(j)=-(s%uu(s%nx,j)+dsqrt(par%g*s%hum(s%nx,j)))*dbetadx(2,j) &
                  -s%vu(s%nx,j)*dbetady(2,j)& !Ap s%vu
                  -dsqrt(par%g*s%hum(s%nx,j))*dvdy(2,j)&
                  +s%Fx(s%nx,j)*inv_ht(j)/par%rho &
                  -s%cfu(s%nx,j)*sqrt(s%uu(s%nx,j)**2+s%vu(s%nx,j)**2)*s%uu(s%nx,j)/s%hum(s%nx,j)&
                  +par%g*dhdx(2,j)
               endif    ! Robert: dry back boundary points
            enddo


            do j=j1,max(s%ny,1)
               if (s%wetu(s%nx,j)==1) then                                                     ! Robert: dry back boundary points
                  betanp1(1,j) = beta(2,j)+ bn(j)*par%dt                                   !Ap toch?
                  alpha2(j)= s%theta0
                  alphanew = 0.d0

                  do jj=1,50
                     !
                     !---------- Lower order bound. cond. ---
                     !
                     s%ur(s%nx,j) = dcos(alpha2(j))/(dcos(alpha2(j))+1.d0)&
                     *((betanp1(1,j)-s%umean(s%nx,j)-2.d0*DSQRT(par%g*0.5*(ht(1,j)+ht(2,j)))))

                     !vert = velocity of the reflected wave = total-specified
                     vert = s%vu(s%nx,j)-s%vmean(s%nx,j)                              !Jaap included s%vmean
                     alphanew = datan(vert/(s%ur(s%nx,j)+1.d-16))                      ! wwvv maybe better atan2
                     if (alphanew .gt. (par%px*0.5d0)) alphanew=alphanew-par%px
                     if (alphanew .le. (-par%px*0.5d0)) alphanew=alphanew+par%px
                     if(dabs(alphanew-alpha2(j)).lt.0.001) exit    ! wwvv can use exit here
                     alpha2(j) = alphanew
                  end do

                  s%uu(s%nx,j) = s%ur(s%nx,j) + s%umean(s%nx,j)
                  ! Ap replaced zs with extrapolation.
                  s%zs(s%nx+1,j) = 1.5*((betanp1(1,j)-s%uu(s%nx,j))**2.d0/4.d0/par%g+.5*(s%zb(s%nx,j)+s%zb(s%nx+1,j)))-&
                  0.5*((beta(1,j)-s%uu(s%nx-1,j))**2.d0/4.d0/par%g+.5*(s%zb(s%nx-1,j)+s%zb(s%nx,j)))

               endif   ! Robert: dry back boundary points
            enddo
         elseif (par%back==BACK_WLEVEL) then
            s%zs(s%nx+1,:)=s%zs0(s%nx+1,:)
         endif  !par%back

      endif !xmpi_isbot


      !
      ! Lateral tide boundary conditions
      !
      if (s%ny>0) then
         ! RJ: Neumann water levels in case of right = 1 or right = 0
         do i=1,s%nx+1
            !! Jaap multiply with wetz(i,ny+1)*wetz(i,1) here to prevent prsssure grdaient over land
            !dzsdnavg=s%wetz(i,s%ny+1)*s%wetz(i,1)*(s%zs0(i,s%ny+1)-s%zs0(i,1))/(s%ndist(i,s%ny+1)-s%ndist(i,1))
            !! Lateral boundary at y=0
            !if (xmpi_isleft) then
            !   s%zs(i,1)=max(s%zs(i,2) - dzsdnavg*s%dnv(i,1),s%zb(i,1))
            !endif
            !! Lateral boundary at y=ny+1
            !if (xmpi_isright) then
            !   s%zs(i,s%ny+1)=max(s%zs(i,s%ny) + dzsdnavg*s%dnv(i,s%ny),s%zb(i,s%ny+1))
            !endif
            if (xmpi_isleft) then
               s%zs(i,1)=max(s%zs(i,2) - s%dzs0dn(i,1)*s%dnv(i,1),s%zb(i,1))
            endif
            if (xmpi_isright) then
               s%zs(i,s%ny+1)=max(s%zs(i,s%ny) + s%dzs0dn(i,s%ny+1)*s%dnv(i,s%ny),s%zb(i,s%ny+1))
            endif
         enddo
      endif !s%ny>0



      !!! Wind boundary conditions
      if (s%windlen>1) then  ! only if non-stationary wind, otherwise waste of computational time
         call LINEAR_INTERP(s%windinpt,s%windxts,s%windlen,par%t,windxnow,indt)
         call LINEAR_INTERP(s%windinpt,s%windyts,s%windlen,par%t,windynow,indt)
         s%windsu = windxnow*dcos(s%alfau) + windynow*dsin(s%alfau)
         ! Can we do this more efficiently?
         s%windnv = windynow*dcos(s%alfav-0.5d0*par%px) - windxnow*dsin(s%alfav-0.5d0*par%px)
      endif
   end subroutine flow_bc

   function flow_lat_bc(s,par,bctype,jbc,jn,udvdxb,vdvdyb,viscvb) result(vbc)
      use params
      use spaceparams
      use xmpi_module
      use paramsconst

      type(spacepars),target                  :: s
      type(parameters)                        :: par
      integer                                 :: bctype
      integer                                 :: jbc,jn ! j-index of boundary and j-index of neighbouring cell
      ! varies for xmpi_isleft and xmpi_isright
      real*8,dimension(s%nx+1)                :: udvdxb,vdvdyb,viscvb   ! advection terms from flow timestep
      real*8,dimension(s%nx+1)                :: vbc    ! result vector for boundary flow
      integer                                 :: i      ! internal variables
      real*8                                  :: advterm
      real*8,save                             :: fc


      if (par%t<=par%dt) fc =2.d0*par%wearth*sin(par%lat)

      if (s%ny==0) then      ! 1D models have special case
         if (bctype==LR_WALL) then
            vbc(:) = 0.d0 ! Only this bc scheme can be applied to s%ny==0 models
         else
            vbc(:) = s%vv(:,jbc)  ! return whatever was already calculated in the flow routine for the cross shore row
         endif
      else
         select case (bctype)
          case (LR_WALL)
            ! No flow at boundary
            vbc(:) = 0.d0
          case (LR_NEUMANN_V)
            ! Calculate vv at boundary using identical forcing as (:,2), but with dzsdy determined from (:,1), which includes
            ! the tide-driven water level gradient as zs(:,1) = zs(:,2) + dzsdytide
            !
            ! in flow timestep we calculate:
            ! vv(jn) = vvold(jn)-dt(gDZS(jn) + R(jn))
            !   -> vvold(jn) = vv(jn) + dt(gDZS(jn) + R(jn))
            !
            ! we state in bc that:
            ! vv(jbc) = vvold(jn)-dt(gDZS(jbc) + R(jn))
            !   -> vv(jbc) = ( vv(jn) + dt(gDZS(jn) + R(jn)) )-dt(gDZS(jbc) + R(jn))
            !   -> vv(jbc) = vv(jn) + dt(gDZS(jn)) - dt(gDZS(jbc)) + dt(R(jn)) - dt(R(jn))
            !   -> vv(jbc) = vv(jn) + dt(gDZS(jn)-gDZS(jbc))
            !   -> vv(jbc) = vv(jn) + dtg(DZS(jn)-DZS(jbc))
            !                  vbc(:) = (vv(:,jn)+par%dt*par%g*(dzsdy(i,jn)-dzsdy(i,jbc))) * wetv(:,jbc)
            ! AANPASSING: Neumann zonder getij
            ! Dano/Jaap: Je kunt niet druk gradient a.g.v getij meenemen zonder ook andere forcerings termen (bijv bodem wrijving)
            ! uit te rekenen. In geval van getij gradient heb je een meer complexe neumann rvw die jij "free" noemt
            vbc(:) = s%vv(:,jn)
          case (LR_NO_ADVEC)
            ! We allow vv at the boundary to be calulated from NLSWE, but only include the advective terms if
            ! they decrease the flow velocity
            do i=2,s%nx
               if (s%wetv(i,jbc)==1) then
                  advterm = udvdxb(i)+vdvdyb(i)-viscvb(i)+ fc*s%uv(i,jbc)
                  if (s%vv(i,jbc)>0.d0) then
                     advterm = max(advterm,0.d0)
                  elseif (s%vv(i,jbc)<0.d0) then
                     advterm = min(advterm,0.d0)
                  else
                     advterm = 0.d0
                  endif
                  vbc(i) = s%vv(i,jbc)- par%dt*(advterm &
                  + par%g*s%dzsdy(i,jbc)&
                  + s%tauby(i,jbc)/(par%rho*s%hvm(i,jbc)) &
                  - par%lwave*s%Fy(i,jbc)/(par%rho*max(s%hvm(i,jbc),par%hmin)) &
                  - par%rhoa*par%Cd*s%windnv(i,jbc)*sqrt(s%windsu(i,jbc)**2+s%windnv(i,jbc)**2)/(par%rho*s%hvm(i,jbc)))  ! Kees: wind correction
                  ! Robert: no Coriolis because dependecy u in this combination of wall / Neumann
               else
                  vbc(i) = 0.d0
               endif
            enddo
            ! now fill i=1 and i=nx+1
            vbc(1) = vbc(2) * s%wetv(1,jbc)
            vbc(s%nx+1) = vbc(s%nx) * s%wetv(s%nx+1,jbc)
          case (LR_NEUMANN)
            ! Dano/Jaap: En dit is dan dus complexe vorm van Neumann (zoals in Delft3D) met getij
            ! We allow vv at the boundary to be calulated from NLSWE without any limitations
            do i=2,s%nx
               if (s%wetv(i,jbc)==1) then
                  vbc(i) = s%vv(i,jbc)- par%dt*(udvdxb(i)+vdvdyb(i)-viscvb(i)&
                  + par%g*s%dzsdy(i,jbc)&
                  + s%tauby(i,jbc)/(par%rho*s%hvm(i,jbc)) &
                  - par%lwave*s%Fy(i,jbc)/(par%rho*max(s%hvm(i,jbc),par%hmin)) &
                  + fc*s%uv(i,jbc) &
                  - par%rhoa*par%Cd*s%windnv(i,jbc)*sqrt(s%windsu(i,jbc)**2+s%windnv(i,jbc)**2)/(par%rho*s%hvm(i,jbc)))  ! Kees: wind correction
               else
                  vbc(i) = 0.d0
               endif
            enddo
            ! now fill i=1 and i=nx+1
            vbc(1) = vbc(2) * s%wetv(1,jbc)
            vbc(s%nx+1) = vbc(s%nx) * s%wetv(s%nx+1,jbc)
          case (LR_ABS_1D)
            ! Dano: weakly reflective boundary
            do i=2,s%nx
               if (s%wetv(i,jbc)==1) then
                  vbc(i) = sqrt(par%g/s%hh(i,jn))*(s%zs(i,jn)-max(s%zb(i,jn),s%zs0(i,jn)))*(jbc-jn)
                  s%zs(i,jbc)=s%zs(i,jn)
               else
                  vbc(i) = 0.d0
               endif
            enddo
            ! now fill i=1 and i=nx+1
            vbc(1) = vbc(2) * s%wetv(1,jbc)
            vbc(s%nx+1) = vbc(s%nx) * s%wetv(s%nx+1,jbc)
         end select
      endif
   end function flow_lat_bc


   subroutine discharge_boundary_h(s,par)
      use params
      use spaceparams
      use xmpi_module
      use readkey_module
      use interp
      use logging_module

      implicit none

      type(spacepars),target                  :: s
      type(parameters)                        :: par

      integer                                 :: i,indx
      integer                                 :: m1,m2,n1,n2
      integer                                 :: m1l,m2l,n1l,n2l
      real*8                                  :: qnow,A,CONST
      logical                                 :: indomain,isborder
      integer                                 :: ishift,jshift
      integer                                 :: iminl,imaxl,jminl,jmaxl

#ifdef USEMPI
      ishift  = s%is(xmpi_rank+1)-1
      jshift  = s%js(xmpi_rank+1)-1
#else
      ishift  = 0
      jshift  = 0
#endif
      iminl = imin_zs
      imaxl = imax_zs
      jminl = jmin_zs
      jmaxl = jmax_zs
      if (xmpi_istop) then
         iminl = imin_zs-1
      endif
      if (xmpi_isbot) then
         imaxl = imax_zs+1
      endif
      if (xmpi_isleft) then
         jminl = jmin_zs-1
      endif
      if (xmpi_isright) then
         jmaxl = jmax_zs+1
      endif

      ! loop through discharge location
      do i = 1,par%ndischarge
         if (s%pntdisch(i).eq.0) then

            ! interpolate discharge timeseries to current time
            call linear_interp(s%tdisch,s%qdisch(:,i),par%ntdischarge,par%t,qnow,indx)

            m1 = s%pdisch(i,1)
            n1 = s%pdisch(i,2)
            m2 = s%pdisch(i,3)
            n2 = s%pdisch(i,4)

            ! convert to local coordinates
            m1l = m1-ishift
            m2l = m2-ishift
            n1l = n1-jshift
            n2l = n2-jshift

            indomain =  (   m1l >= iminl .and. m1l <= imaxl .and.         &
            n1l >=jminl .and. n1l <= jmaxl       ) .or.  &
            (   m2l >=iminl .and. m2l <= imaxl .and.         &
            n2l >=jminl  .and. n2l <= jmaxl       )

            m1l = min(max(m1l,1),s%nx)
            m2l = min(max(m2l,1),s%nx)
            n1l = min(max(n1l,1),s%ny)
            n2l = min(max(n2l,1),s%ny)

            A = 0.d0

            ! add momentum in either u- or v-direction
            if (n1.eq.n2) then

               if (indomain) then
                  isborder = (n1l.eq.1 .and. xmpi_isleft) .or. (n1l.eq.s%ny .and. xmpi_isright)
                  m1l=max(m1l,iminl)
                  m2l=min(m2l,imaxl)

                  ! make sure discharge is an inflow at border and water levels are defined
                  if (isborder) then
                     if (n1l.gt.1) then
                        qnow = -1*qnow
                        s%hv(m1l:m2l,n1l) = s%hh(m1l:m2l,n1l)
                     elseif (n1l.lt.s%ny) then
                        s%hv(m1l:m2l,n1l) = s%hh(m1l:m2l,n1l+1)
                     endif
                  endif

                  A = sum(s%hv(m1l:m2l,n1l)**1.5*s%dsv(m1l:m2l,n1l))
               endif

#ifdef USEMPI
               call xmpi_allreduce(A,MPI_SUM)
#endif

               ! compute distribution of discharge over discharge cells according to:
               !     vv = CONST * hv^0.5
               !     qy = CONST * hv^1.5
               !     Q  = CONST * hv^1.5 * dsv

               if (indomain) then
                  CONST = qnow/max(A,par%eps)
                  s%vv(m1l:m2l,n1l) = CONST*sqrt(s%hv(m1l:m2l,n1l))
                  s%qy(m1l:m2l,n1l) = CONST*sqrt(s%hv(m1l:m2l,n1l))*s%hv(m1l:m2l,n1l)
               endif

            elseif (m1.eq.m2) then

               if (indomain) then
                  isborder = (m1l.eq.1 .and. xmpi_istop) .or. (m1l.eq.s%nx .and. xmpi_isbot)
                  n1l=max(n1l,jminl)
                  n2l=min(n2l,jmaxl)

                  ! make sure discharge is an inflow at border and water levels are defined
                  if (isborder) then
                     if (m1l.gt.1) then
                        qnow = -1*qnow
                        s%hu(m1l,n1l:n2l) = s%hh(m1l,n1l:n2l)
                     elseif (m1l.lt.s%nx) then
                        s%hu(m1l,n1l:n2l) = s%hh(m1l+1,n1l:n2l)
                     endif
                  endif

                  A = sum(s%hu(m1l,n1l:n2l)**1.5*s%dnu(m1l,n1l:n2l))
               endif

#ifdef USEMPI
               call xmpi_allreduce(A,MPI_SUM)
#endif

               ! compute distribution of discharge over discharge cells according to:
               !     uu = CONST * hu^0.5
               !     qx = CONST * hu^1.5
               !     Q  = CONST * hu^1.5 * dnu

               if (indomain) then
                  CONST = qnow/max(A,par%eps)
                  s%uu(m1l,n1l:n2l) = CONST*sqrt(s%hu(m1l,n1l:n2l))
                  s%qx(m1l,n1l:n2l) = CONST*sqrt(s%hu(m1l,n1l:n2l))*s%hu(m1l,n1l:n2l)
               endif
            endif
         endif
      enddo
   end subroutine discharge_boundary_h

   subroutine discharge_boundary_v(s,par)
      use params
      use spaceparams
      use xmpi_module
      use readkey_module
      use interp
      use logging_module

      implicit none

      type(spacepars),target                  :: s
      type(parameters)                        :: par

      integer                                 :: i,indx
      integer                                 :: m1,n1
      integer                                 :: m1l,n1l
      real*8                                  :: qnow
      integer                                 :: ishift,jshift


#ifdef USEMPI
      ishift  = s%is(xmpi_rank+1)-1
      jshift  = s%js(xmpi_rank+1)-1
#else
      ishift  = 0
      jshift  = 0
#endif

      do i = 1,par%ndischarge
         if (s%pntdisch(i).eq.1) then
            call linear_interp(s%tdisch,s%qdisch(:,i),par%ntdischarge,par%t,qnow,indx)

            m1 = s%pdisch(i,1)
            n1 = s%pdisch(i,2)

            ! convert to local coordinates
            m1l = min(max(m1-ishift,1),s%nx)
            n1l = min(max(n1-jshift,1),s%ny)

            ! add mass
            s%zs(m1l,n1l) = s%zs(m1l,n1l)+qnow*par%dt*s%dsdnzi(m1l,n1l)
         endif
      enddo
   end subroutine discharge_boundary_v

   !
   !==============================================================================
   subroutine velocity_Boundary(bcfile,ug,vg,zg,wg,dug,dvg,nyg,t,isSet_U,isSet_V,isSet_Z,isSet_W,isSet_dU, &
   isSet_dV,isSet_Q,force_init,bcst)
      !==============================================================================
      !

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !
      !--------------------------        PURPOSE         ----------------------------
      !
      !   Reads boundary values for the velocity from a file. Used as forcing in the
      !   nonh module

      !--------------------------        METHOD         ----------------------------


      !--------------------------     DEPENDENCIES       ----------------------------

      use params
      use xmpi_module, only: Halt_Program
      use spaceparams
      use logging_module
      use readkey_module, only: lowercase

      implicit none

      include 'nh_pars.inc'               !Default precision etc.

      !--------------------------     ARGUMENTS          ----------------------------
      !
      character(len=*)              ,intent(in)    :: bcfile
      integer(kind=iKind)           ,intent(in)    :: nyg
      real(kind=rKind),dimension(nyg+1)            :: ug,vg,zg,wg,dug,dvg
      real(kind=rKind)                ,intent(in)  :: t
      logical,intent(out)                          :: isSet_U,isSet_V,isSet_Z,isSet_W,isSet_dU,isSet_dV,isSet_Q
      logical,intent(in),optional                  :: force_init
      real(kind=rKind),intent(in),optional         :: bcst
      !

      !--------------------------     LOCAL VARIABLES    ----------------------------
      character(len=iFileNameLen),save          :: filename_U = ''
      integer                    ,save          :: unit_U = iUnit_U

      logical,save                              :: lVarU = .false.
      logical,save                              :: lIsEof = .false.
      logical,save                              :: lNH_boun_U  = .false.
      logical,save                              :: initialize  = .true.
      logical                                   :: lExists,lOpened
      character(slen)                          :: string
      character(len=2),allocatable,dimension(:),save :: header
      integer(kind=ikind)                       :: iAllocErr,ier
      integer(kind=ikind),save                  :: nvar

      integer(kind=ikind),save                  :: iZ = 0
      integer(kind=ikind),save                  :: iU = 0
      integer(kind=ikind),save                  :: iV = 0
      integer(kind=ikind),save                  :: idU = 0
      integer(kind=ikind),save                  :: idV = 0
      integer(kind=ikind),save                  :: iQ  = 0
      integer(kind=ikind),save                  :: idQ = 0
      integer(kind=ikind),save                  :: iW = 0
      integer(kind=ikind),save                  :: iT = 0
      integer(kind=ikind)                       :: i

      logical                                   :: fi_local
      real(kind=rKind)                          :: bcstarttime


      real(kind=rKind),allocatable,dimension(:,:)    :: tmp

      real(kind=rKind),allocatable,dimension(:),save :: u0 !
      real(kind=rKind),allocatable,dimension(:),save :: u1 !

      real(kind=rKind),allocatable,dimension(:),save :: v0 !
      real(kind=rKind),allocatable,dimension(:),save :: v1 !

      real(kind=rKind),allocatable,dimension(:),save :: z0 !
      real(kind=rKind),allocatable,dimension(:),save :: z1 !

      real(kind=rKind),allocatable,dimension(:),save :: w0 !
      real(kind=rKind),allocatable,dimension(:),save :: w1 !

      real(kind=rKind),allocatable,dimension(:),save :: du0 !
      real(kind=rKind),allocatable,dimension(:),save :: du1 !

      real(kind=rKind),allocatable,dimension(:),save :: dv0 !
      real(kind=rKind),allocatable,dimension(:),save :: dv1 !

      real(kind=rKind),save                          :: t0 = 0.
      real(kind=rKind),save                          :: t1 = 0.

      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------

      if (xmaster) then
         if (present(bcst)) then
            bcstarttime = bcst
         else
            bcstarttime = 0.d0
         endif

         if (present(force_init)) then
            fi_local = force_init
         else
            fi_local = .false.
         endif

         if (bcfile/=filename_U .or. fi_local) then
            filename_U=trim(bcfile)
            initialize = .true.
         endif
         !
         ! INITIALIZE NEW FILE
         !
         if (initialize) then
            initialize = .false.

            !Does the file exist and is it already opened (reuse files?)
            inquire(file=trim(filename_U),EXIST=lExists,OPENED=lOpened)

            ! Close previous file if needed
            if(lOpened) then
               close(unit_U)
            endif

            ! Deallocate previous arrays if needed
            if (allocated(header)) deallocate(header)
            if (allocated(tmp)) deallocate(tmp)
            if (allocated(u0)) deallocate(u0)
            if (allocated(u1)) deallocate(u1)
            if (allocated(v0)) deallocate(v0)
            if (allocated(v1)) deallocate(v1)
            if (allocated(z0)) deallocate(z0)
            if (allocated(z1)) deallocate(z1)
            if (allocated(w0)) deallocate(w0)
            if (allocated(w1)) deallocate(w1)
            if (allocated(du0)) deallocate(du0)
            if (allocated(dv0)) deallocate(dv0)
            if (allocated(du1)) deallocate(du1)
            if (allocated(dv1)) deallocate(dv1)

            if (lExists) then
               !If exists, open
               open(unit_U,file=filename_U,access='SEQUENTIAL',action='READ',form='FORMATTED')
               lNH_boun_U = .true.
            else
               lNH_boun_U = .false.
               call writelog('lswe','','Error bc file does not exist: ',trim(filename_U))
               call halt_program
               return
            endif


            !SCALAR OR VECTOR INPUT?
            read(unit_U,fmt=*,iostat=ier) string
            if (ier .ne. 0) then
               call report_file_read_error(filename_U)
            endif
            string = lowercase(string)
            if   (string == 'scalar') then
               lVaru = .false.
            elseif (string == 'vector') then
               lVaru = .true.
            else

            endif

            !READ NUMBER OF VARIABLES
            read(unit_U,fmt=*,iostat=ier) nvar
            if (ier .ne. 0) then
               call report_file_read_error(filename_U)
            endif

            !READ HEADER LINES
            allocate(header(nvar))
            read(unit_U,fmt=*,iostat=ier) header
            if (ier .ne. 0) then
               call report_file_read_error(filename_U)
            endif

            do i=1,nvar
               if (header(i) == 'z' .or. header(i) == 'Z' .or. &
               header(i) == 'zs' .or. header(i) == 'ZS' .or. &
               header(i) == 'Zs') then
                  iZ = i;
               elseif (header(i) == 't' .or. header(i) == 'T') then
                  iT = i;
               elseif (header(i) == 'u' .or. header(i) == 'U') then
                  iU = i;
               elseif (header(i) == 'v' .or. header(i) == 'V') then
                  iV = i;
               elseif (header(i) == 'w' .or. header(i) == 'W') then
                  iW = i;
               elseif (header(i) == 'dU' .or. header(i) == 'DU') then
                  idU = i;
               elseif (header(i) == 'dV' .or. header(i) == 'DV') then
                  idV = i;
               elseif (header(i) == 'q' .or. header(i) == 'Q') then
                  iQ = i;
                  iU = i;
               elseif (header(i) == 'dq' .or. header(i) == 'dQ') then
                  idQ = i;
                  idU = i
               endif
            enddo

            !Allocate arrays at old and new timelevel
            if (iU > 0) then
               allocate(u0(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               allocate(u1(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               u0 = 0.0d0
               u1 = 0.0d0
            endif

            if (iV > 0) then
               allocate(v0(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               allocate(v1(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               v0 = 0.0d0
               v1 = 0.0d0
            endif

            if (iZ > 0) then
               allocate(z0(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               allocate(z1(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               z0 = 0.0d0
               z1 = 0.0d0
            endif

            if (iW > 0) then
               allocate(w0(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               allocate(w1(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               w0 = 0.0d0
               w1 = 0.0d0
            endif

            if (idU > 0) then
               allocate(du0(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               allocate(du1(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               du0 = 0.0d0
               du1 = 0.0d0
            endif
            if (idV > 0) then
               allocate(dv0(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               allocate(dv1(nyg+1),stat=iAllocErr); if (iAllocErr /= 0) call Halt_program
               dv0 = 0.0d0
               dv1 = 0.0d0
            endif

            allocate(tmp(nyg+1,nvar-1))
            !Read two first timelevels
            call velocity_Boundary_read(t0,tmp,Unit_U,lVaru,lIsEof,nvar)
            t0 = t0+bcstarttime
            if (iU >0) u0 = tmp(:,iU-1)
            if (iV >0) v0 = tmp(:,iV-1)
            if (iZ >0) z0 = tmp(:,iZ-1)
            if (iW >0) w0 = tmp(:,iW-1)
            if (idU >0) dU0 = tmp(:,idU-1)
            if (idV >0) dV0 = tmp(:,idV-1)
            !
            if (lIsEof) then
               t1=t0
               if (iU >0) u1=u0
               if (iV >0) v1=v0
               if (iZ >0) z1=z0
               if (iW >0) w1=w0
               if (idU >0) du1=du0
               if (idV >0) dv1=dv0
            else
               call velocity_Boundary_read(t1,tmp,Unit_U,lVaru,lIsEof,nvar)
               t1 = t1+bcstarttime
               if (iU >0) u1 = tmp(:,iU-1)
               if (iV >0) v1 = tmp(:,iV-1)
               if (iZ >0) z1 = tmp(:,iZ-1)
               if (iW >0) w1 = tmp(:,iW-1)
               if (idu >0) du1 = tmp(:,idU-1)
               if (idv >0) dv1 = tmp(:,idV-1)
            endif
            deallocate(tmp)
            if (iU > 0) then
               ug  = u0 + (u1-u0)*(t-t0)/(t1-t0)
               isSet_U = .true.
            else
               ug = 0
               isSet_U = .false.
            endif

            if (iV > 0) then
               vg  = v0 + (v1-v0)*(t-t0)/(t1-t0)
               isSet_V = .true.
            else
               vg = 0
               isSet_V = .false.
            endif

            if (iZ > 0) then
               zg = z0 + (z1-z0)*(t-t0)/(t1-t0)
               isSet_Z = .true.
            else
               zg = 0
               isSet_Z = .false.
            endif
            if (iW > 0) then
               wg = w0 + (w1-w0)*(t-t0)/(t1-t0)
               isSet_W = .true.
            else
               wg = 0
               isSet_W = .false.
            endif

            if (idu > 0) then
               dug = dU0 + (dU1-dU0)*(t-t0)/(t1-t0)
               isSet_dU = .true.
            else
               dug = 0
               isSet_dU = .false.
            endif

            if (idv > 0) then
               dvg = dV0 + (dV1-dV0)*(t-t0)/(t1-t0)
               isSet_dV = .true.
            else
               dvg = 0
               isSet_dV = .false.
            endif

            if (iQ > 0) then
               !
               isSet_Q = .true.
            else
               isSet_Q = .false.
               !
            endif


            return
         endif   ! initialize
         !
         ! REGULAR TIMESTEP READ BOUNDARY CONDITIONS
         !
         if (lNH_boun_U) then
            if (.not. lIsEof) then
               allocate(tmp(nyg+1,nvar-1))
               !If current time not located within interval read next line
               do while (.not. (t>=t0 .and. t<t1))
                  t0 = t1
                  if (iU >0) u0 = u1
                  if (iV >0) v0 = v1
                  if (iZ >0) z0 = z1
                  if (iW >0) w0 = w1
                  if (idu >0) du0 = du1
                  if (idv >0) dv0 = dv1
                  call velocity_Boundary_read(t1,tmp,Unit_U,lVaru,lIsEof,nvar)
                  t1 = t1+bcstarttime
                  if (iU >0) u1 = tmp(:,iU-1)
                  if (iV >0) v1 = tmp(:,iV-1)
                  if (iZ >0) z1 = tmp(:,iZ-1)
                  if (iW >0) w1 = tmp(:,iW-1)
                  if (idu >0) du1 = tmp(:,idu-1)
                  if (idv >0) dv1 = tmp(:,idv-1)
                  if (lIsEof) exit !Exit on end of file condition
               end do
               deallocate(tmp)

               if(lIsEof) then
                  if (iU >0) ug = u1
                  if (iV >0) vg = v1
                  if (iZ >0) zg = z1
                  if (iW >0) wg = w1
                  if (idu >0) dug = du1
                  if (idv >0) dvg = dv1
                  return
               else
                  !Linear interpolation of u in time
                  if (iU > 0) then
                     ug  = u0 + (u1-u0)*(t-t0)/(t1-t0)
                     isSet_U = .true.
                  else
                     ug = 0
                     isSet_U = .false.
                  endif

                  if (iV > 0) then
                     vg  = v0 + (v1-v0)*(t-t0)/(t1-t0)
                     isSet_V = .true.
                  else
                     vg = 0
                     isSet_V = .false.
                  endif

                  if (iZ > 0) then
                     zg = z0 + (z1-z0)*(t-t0)/(t1-t0)
                     isSet_Z = .true.
                  else
                     zg = 0
                     isSet_Z = .false.
                  endif
                  if (iW > 0) then
                     wg = w0 + (w1-w0)*(t-t0)/(t1-t0)
                     isSet_W = .true.
                  else
                     wg = 0
                     isSet_W = .false.
                  endif
                  if (idu > 0) then
                     dug = dU0 + (dU1-dU0)*(t-t0)/(t1-t0)
                     isSet_dU = .true.
                  else
                     dug = 0
                     isSet_dU = .false.
                  endif
                  if (idv > 0) then
                     dvg = dV0 + (dV1-dV0)*(t-t0)/(t1-t0)
                     isSet_dV = .true.
                  else
                     dvg = 0
                     isSet_dV = .false.
                  endif
                  !
                  if (iQ > 0) then
                     !
                     isSet_Q = .true.
                     !
                  else
                     isSet_Q = .false.
                  endif
                  !
               endif
            else
               !If end of file the last value which was available is used until the end of the computation
               ug = u1
               vg = v1
               zg = z1
            endif ! .not. lIsEof
         endif ! lNH_boun_U
      endif ! xmaster
   end subroutine velocity_Boundary

   !==============================================================================
   subroutine velocity_Boundary_read(t,vector,iUnit,isvec,iseof,nvar)
      !==============================================================================

      !--------------------------        PURPOSE         ----------------------------

      ! Reads scalar/vector for timeseries

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------



      !                                - NONE -

      !--------------------------     DEPENDENCIES       ----------------------------

      use xmpi_module, only: Halt_Program
      use logging_module, only: report_file_read_error

      implicit none

      include 'nh_pars.inc'               !Default precision etc.
      !--------------------------     ARGUMENTS          ----------------------------

      real(kind=rKind),intent(inout)                   :: t
      real(kind=rKind),intent(inout),dimension(:,:)    :: vector
      integer(kind=iKind),intent(in)                   :: iUnit
      integer(kind=iKind),intent(in)                   :: nvar
      logical            ,intent(in)                   :: isvec
      logical            ,intent(out)                  :: iseof

      !--------------------------     LOCAL VARIABLES    ----------------------------

      real(kind=rKind)                            :: scalar(nvar-1)
      integer(kind=iKind)                         :: i
      integer(kind=iKind)                         :: ioStat
      character(len=slen)                         :: ername

      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------


      iseof = .false.

      if (isvec) then
         read(iUnit,ioStat=ioStat,fmt=*,end=9000) t,vector
      else
         read(iUnit,ioStat=ioStat,fmt=*,end=9000) t,scalar
         do i=1,nvar-1
            vector(:,i) = scalar(i)
         enddo
      endif

      if (iostat /= 0) then
         inquire(iUnit,name=ername)
         call report_file_read_error(ername)
      endif

      return
9000  iseof = .true.

   end subroutine velocity_Boundary_read

end module boundaryconditions
