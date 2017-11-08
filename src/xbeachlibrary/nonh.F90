!==============================================================================
!                               MODULE NH_MAIN
!==============================================================================

! DATE               AUTHOR               CHANGES
!
! october 2009       Pieter Bart Smit     New module
! october 2014       Pieter Bart Smit     Included an updated pressure correction scheme/general maintenance

! -- REMARKS --
! There are "two" mutually exclusive implementations of the non-hydrostatic model
! in this file. The first is the original non-hydrostatic code (see Draft report)
! which has been left untouched (except for fixing some consistency errors). The
! second implementation is the reduced two-layer formulation. The new implementation
! uses two-layers in the vertical (controlled by a thickness parameter) and should
! give better results for shorter wave components. Moreover, it contains the previous
! code as a special case (thickness of the lower layer is set to 0).
!
! INCLUDED SUBROUTINES
! + nonh_cor
!   This is the entry routine for the nonhydrostatic code when called from
!   flowtimestep. Depending on the parameter nh_red it calls the old nh-routines
!   or the new reduced routines.
!
! + nonh_init
!   Initializes the various resources required.
!
! + nonh_1lay_cor
!   The "old" non-hydrostatic correction routine, calculates the new dynamic
!   pressure.
!
! + nonh_1lay_pred
!
!
module nonh_module

   implicit none
   save


   ! If mpi is defined, the non-hydrostatic module is NOT included in the compilation
   ! to avoid unwanted side effects.

   !******************************************************************************
   !                                 INTERFACE
   !******************************************************************************

   private

   !----------------------------- PARAMETERS -----------------------------------

   include 'nh_pars.inc'               !Default precision etc.

   !----------------------------- VARIABLES  -----------------------------------

   logical                      :: initialized   = .false.    !Indicates whether or not the nh_module has been initialized

   !
   ! -- Numerical parameter --
   !
   real (kind=rKind), allocatable, dimension(:,:)   :: wcoef  !The relative thickness of the lower layer in the reduced
   !two-layer model at z-points (optimisation parameter in "old model")
   !
   ! -- Additional physical variables --
   !
   real (kind=rKind), allocatable, dimension(:,:)   :: dp     !Change in pressure between the old and
   !new timestep
   real (kind=rKind), allocatable, dimension(:,:)   :: Wm     !The mean vertical velocity in the top
   !layer (nonhq3d)
   real (kind=rKind), allocatable, dimension(:,:)   :: Wm_old !The mean vertical velocity in the top
   !layer on the previous timestep (nonhq3d)
   real (kind=rKind), allocatable, dimension(:,:)   :: Wm0    !The mean vertical velocity in the bottom
   !layer (nonhq3d)
   real (kind=rKind), allocatable, dimension(:,:)   :: dU0    !The velocity (s-dir) difference between the top
   !and bottom layers at the previous timestep
   real (kind=rKind), allocatable, dimension(:,:)   :: dV0    !The velocity (n-dir) difference between the top
   !and bottom layers at the previous timestep
   real (kind=rKind), allocatable, dimension(:,:)   :: omega  !The velocity vertical transport velocity

   !
   ! -- Arrays used for building the Poisson equation --
   !
   real (kind=rKind), allocatable, dimension(:,:,:) :: au     !Pressure coeficients for u(i,j), au(0,i,j) refers to
   !s%nhpres(i,j) and au(1,i,j) to s%nhpres(i+1,j)
   real (kind=rKind), allocatable, dimension(:,:,:) :: av     !Pressure coeficients for v(i,j), av(0,i,j) refers to
   !s%nhpres(i,j) and av(1,i,j) to s%nhpres(i,j+1)
   real (kind=rKind), allocatable, dimension(:,:,:) :: adu    !Pressure coeficients for du(i,j), adu(0,i,j) refers to
   !s%nhpres(i,j) and au(1,i,j) to s%nhpres(i+1,j)
   real (kind=rKind), allocatable, dimension(:,:,:) :: adv    !Pressure coeficients for dv(i,j), adv(0,i,j) refers to
   !s%nhpres(i,j) and adv(1,i,j) to s%nhpres(i,j+1)
   real (kind=rKind), allocatable, dimension(:,:,:) :: mat    !5-diagonal Pressure matrix.
   ! Diagonal Pressure Point
   ! 1        i  ,j
   ! 2        i-1,j
   ! 3        i+1,j
   ! 4        i  ,j-1
   ! 5        i  ,j+1
   real (kind=rKind), allocatable, dimension(:,:)   :: rhs    !RHS of the pressure matrix

   real (kind=rKind), allocatable, dimension(:,:,:) :: aws    !Pressure coeficients for ws(k,i,j) !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:,:,:) :: awb    !Pressure coeficients for wb(k,i,j) !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:,:)   :: aur    !RHS for du(i,j)                    !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:,:)   :: avr    !RHS for dv(i,j)                    !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:,:)   :: awbr   !Pressure coeficients for wb(k,i,j) !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:,:)   :: awsr   !Pressure coeficients for ws(k,i,j) !**OBSOLETE IN nonhq3d**

   !
   ! -- Additional grid parameters --
   !
   real (kind=rKind), allocatable, dimension(:,:)   :: zsu    !free surface-level in u-velocity points
   real (kind=rKind), allocatable, dimension(:,:)   :: zsv    !free surface-level in v-velocity points
   real (kind=rKind), allocatable, dimension(:,:)   :: zbu    !free bed-level in u-velocity points
   real (kind=rKind), allocatable, dimension(:,:)   :: zbv    !free bed-level in v-velocity points

   real (kind=rKind), allocatable, dimension(:)     ::  dxz   !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:)     ::  dyz   !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:)     ::  dxu   !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:)     ::  dyv   !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:)     ::  ddxz  !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:)     ::  ddyz  !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:)     ::  ddxu  !**OBSOLETE IN nonhq3d**
   real (kind=rKind), allocatable, dimension(:)     ::  ddyv  !**OBSOLETE IN nonhq3d**

   integer(kind=iKind),allocatable,dimension(:,:)   :: nonhU
   integer(kind=iKind),allocatable,dimension(:,:)   :: nonhV
   integer(kind=iKind),allocatable,dimension(:,:)   :: nonhZ
   real (kind=rKind), allocatable, dimension(:,:)   :: lbreakcond

   !--- PUBLIC SUBROUTINES ---
   public nonh_init
   public nonh_cor
   public nonh_init_wcoef

   !--- PRIVATE SUBROUTINES

   !        NONE
contains
   !
   !******************************************************************************
   !                             SUBROUTINES/FUNCTIONS
   !******************************************************************************
   !

   !
   !==============================================================================
   subroutine nonh_cor(s,par,ipredcor,uu0,vv0)
      !==============================================================================
      !
      !
      ! DATE               AUTHOR               CHANGES
      !
      ! October  2014     Pieter Bart Smit      new subroutine
      !
      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------

      !  This is the new entry routine for the non-hydrostatic module. It is called twice during each timestep by the
      !  flow module, (I) (withe path=0) before the second-order corrections are applied (if applicable) to include the known pressure
      !  explicitly in the momentum equations,  and (II) (path=1) just after the second order advection correction to calculate
      !  the new pressures implicitly. The codepath taken is specified with the path parameter (path=0 for step I, and path=1 for
      !  step II).
      !
      !  Moreover, depending on the parameter nonhq3d, which activates the new reduced 2 layer model, again different codepaths
      !  are taken:
      !
      !  IF nonhq3d=0
      !         call the "old" versions of the non-hydrostatic routines
      !             (path=0)  nonh_1lay_pred  (renamed from nonh_explicit)
      !             (path=1)  nonh_1lay_cor   (renamed from nonh_cor)
      !  IF nonhq3d=1
      !         call the "old" versions of the non-hydrostatic routine
      !            IF 1D
      !               (path=0)  nonh_2lay_pred_2DV  (renamed from nonh_explicit)
      !               (path=1)  nonh_2lay_cor_2DV   (renamed from nonh_cor)
      !            ELSE
      !               (path=0)  nonh_2lay_pred_3D  (renamed from nonh_explicit)
      !               (path=1)  nonh_2lay_cor_3D   (renamed from nonh_cor)
      !
      !--------------------------     DEPENDENCIES       ----------------------------

      use spaceparams
      use params
      use xmpi_module

      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars) ,intent(inout)                       :: s
      type(parameters),intent(in)                          :: par
      integer         ,intent(in)                          :: ipredcor !Select if we want to apply the explicit predicition / or the implicit correction

      real(kind=rkind),dimension(s%nx+1,s%ny+1),intent(in) :: uu0
      real(kind=rkind),dimension(s%nx+1,s%ny+1),intent(in) :: vv0

      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------

      select case ( par%nonhq3d )
         !
       case (1)
         !
         ! The "new" reduced two-layer formulations are used
         !
         if ( ipredcor == 0 ) then
            ! -- Predictor
            if ( s%ny == 0 ) then
               ! -- 2dv code path
               call nonh_2lay_pred_2dV(s,par,uu0)
               !
            else
               ! -- "3d" code path
               call nonh_2lay_pred_3d(s,par,uu0,vv0)
               !
            endif
            !
         else
            !
            if ( s%ny == 0 ) then
               !
               call nonh_2lay_cor_2dV(s,par)
               !
            else
               !
               call nonh_2lay_cor_3d(s,par)
               !
            endif
            !
         endif
         !
       case default
         !
         ! The "old" code is used
         !
         if ( ipredcor == 0 ) then
            ! -- Predictor
            call nonh_1lay_pred(s,par)
            !
         else
            ! -- Corrector
            call nonh_1lay_cor(s,par)
            !
         endif
         !
      end select
      !
   end subroutine nonh_cor


   !
   !==============================================================================
   subroutine nonh_init(s,par)
      !==============================================================================
      !

      ! DATE               AUTHOR               CHANGES
      !
      ! October 2010       Pieter Bart Smit     New Subroutine

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !
      !--------------------------        PURPOSE         ----------------------------
      !
      !   Initializes nh subroutines
      !
      !--------------------------     DEPENDENCIES       ----------------------------

      use spaceparams
      use params

      !--------------------------     ARGUMENTS          ----------------------------
      !
      type(spacepars) ,intent(inout) :: s
      type(parameters),intent(in)    :: par
      !

      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------

      !    open(unit=PrintFileUnit,file=trim(PrintFileName))

      if(.not. xcompute) return

      allocate(au (0:1,s%nx+1,s%ny+1)); au = 0.0_rKind
      allocate(av (0:1,s%nx+1,s%ny+1)); av = 0.0_rKind
      allocate(mat(5,s%nx+1,s%ny+1))  ;  mat = 0.0_rKind
      allocate(rhs(  s%nx+1,s%ny+1))  ;  rhs = 0.0_rKind
      allocate(zbu(  s%nx+1,s%ny+1))  ;  zbu = 0.0_rKind
      allocate(zbv(  s%nx+1,s%ny+1))  ;  zbv = 0.0_rKind
      allocate(zsu(  s%nx+1,s%ny+1))  ;  zsu = 0.0_rKind
      allocate(zsv(  s%nx+1,s%ny+1))  ;  zsv = 0.0_rKind
      allocate(nonhU(  s%nx+1,s%ny+1)); nonhU = 1
      allocate(nonhV(  s%nx+1,s%ny+1)); nonhV = 1
      allocate(nonhZ(  s%nx+1,s%ny+1)); nonhZ = 1
      allocate(lbreakcond(s%nx+1,s%ny+1)); lbreakcond = par%maxbrsteep
      allocate(dp(s%nx+1,s%ny+1))     ; dp = 0.0_rKind
      allocate(Wm(s%nx+1,s%ny+1))     ;Wm     = 0.0_rKind
      allocate(Wm_old(s%nx+1,s%ny+1)) ;Wm_old = 0.0_rKind
      allocate(Wcoef(s%nx+1,s%ny+1))  ;Wcoef = 1.0_rKind
      allocate(aws (5,s%nx+1,s%ny+1)) ; aws = 0.0_rKind
      !
      ! The reduced two-layer model
      !
      if ( par%nonhq3d == 1 ) then
         !
         allocate(omega(s%nx+1,s%ny+1));  omega = 0.0_rKind
         allocate(adu (0:1,s%nx+1,s%ny+1)); adu = 0.0_rKind
         allocate(adv (0:1,s%nx+1,s%ny+1)); adv = 0.0_rKind
         allocate(Wm0(s%nx+1,s%ny+1))     ;Wm0     = 0.0_rKind
         allocate(dU0(s%nx+1,s%ny+1))     ;dU0 = 0.0_rKind
         allocate(dV0(s%nx+1,s%ny+1))     ;dV0 = 0.0_rKind
         !
      else
         !
         ! These grid variables are obsolete in the new reduced model
         !
         allocate(dxz (  s%nx+1))        ;  dxz  = 0.0_rKind
         allocate(dyz (  s%ny+1))        ;  dyz  = 0.0_rKind
         allocate(dxu (  s%nx+1))        ;  dxu  = 0.0_rKind
         allocate(dyv (  s%ny+1))        ;  dyv  = 0.0_rKind

         allocate(ddxz (  s%nx+1))       ;  ddxz  = 0.0_rKind
         allocate(ddyz (  s%ny+1))       ;  ddyz  = 0.0_rKind
         allocate(ddxu (  s%nx+1))       ;  ddxu  = 0.0_rKind
         allocate(ddyv (  s%ny+1))       ;  ddyv  = 0.0_rKind
         allocate(awb (5,s%nx+1,s%ny+1)) ; awb = 0.0_rKind
         allocate(aur(   s%nx+1,s%ny+1)) ; aur = 0.0_rKind
         allocate(avr(   s%nx+1,s%ny+1)) ; avr = 0.0_rKind
         allocate(awbr(s%nx+1,s%ny+1))   ; awbr = 0.0_rKind
         allocate(awsr(s%nx+1,s%ny+1))   ; awsr = 0.0_rKind

         !Initialize grid variables
         dyz  = s%dnz(1,:)
         dyv  = s%dnv(1,:)
         ddyz = 1.0_rKind/dyz
         ddyv = 1.0_rKind/dyv

         dxu  = s%dsu(:,1)
         dxz  = s%dsz(:,1)
         ddxu = 1.0_rKind/dxu
         ddxz = 1.0_rKind/dxz
         !
      endif

      mat(1,1,:)      = 1.0_rKind
      mat(1,s%nx+1,:) = 1.0_rKind
      mat(1,:,1)      = 1.0_rKind
      mat(1,:,s%ny+1) = 1.0_rKind

      initialized   = .true.

      if ( par%nonhq3d == 1 ) then
         !
         wcoef = par%nhlay
         !
      else
         !
         call nonh_init_wcoef(s,par)
         !
      endif
      !

      !Initialize levels at u/v points (zu,zbu etc.)
      call zuzv(s)

   end subroutine nonh_init




   !==============================================================================
   !                        "ORIGINAL" NON-HYDROSTATIC ROUTINES
   !==============================================================================
   !
   !
   !
   !==============================================================================
   subroutine nonh_1lay_cor(s,par)
      !==============================================================================
      !

      ! DATE               AUTHOR               CHANGES
      !
      ! October 2010       Pieter Bart Smit     New Subroutine
      ! November  2010     Pieter Bart Smit     Added explicit prediction for 2nd order scheme

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------

      !  Releases resources

      !--------------------------     DEPENDENCIES       ----------------------------
      use spaceparams
      use params
      use solver_module
      use xmpi_module
      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars) ,intent(inout)                       :: s
      type(parameters),intent(in)                          :: par

      !--------------------------     LOCAL VARIABLES    ----------------------------

      !Indices
      integer(kind=iKind)                     :: i               !Index variables
      integer(kind=iKind)                     :: j               !Index variables
      integer(kind=iKind)                     :: jmin,jmax       !Index variables for superfast1D

      ! logical
      logical                                 :: sf1d            !Switch for superfast1D

      !Work

      real(kind=rKind)                        :: dzs_e,dzs_w
      real(kind=rKind)                        :: dzs_s,dzs_n
      real(kind=rKind)                        :: dzb_e,dzb_w
      real(kind=rKind)                        :: dzb_n,dzb_s


      !
      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------
      !
      dzs_e = 0
      dzs_w = 0
      dzs_s = 0
      dzs_n = 0
      dzb_e = 0
      dzb_w = 0
      dzb_n = 0
      dzb_s = 0
      if (s%ny>0) then
         sf1d = .false.
         jmin = 2
         jmax = s%ny
      else
         sf1d = .true.
         jmin = 1
         jmax = 1
      endif

      !call timer_start(timer_flow_nonh)
      if (.not. initialized) then
         call nonh_init(s,par)
      endif

      call zuzv(s)

      !Built pressure coefficients U
      aur  = s%uu
      avr  = s%vv


      aur(1,:)       = s%uu(1,:)
      aur(s%nx,:)    = s%uu(s%nx,:)

      avr(:,1)      = s%vv(:,1)
      if (s%ny>0) avr(:,s%ny)   = s%vv(:,s%ny)

      !Built pressure coefficients for W

      !AW Bottom
      if (sf1d) then
         do i=2,s%nx
            if (nonhZ(i,1) == 1) then
               if     (nonhU(i,1)*nonhU(i-1,1) == 1) then
                  dzb_e = .5_rKind*(zbu(i,1)-zbu(i-1,1))*ddxz(i)
                  dzb_w = dzb_e
               elseif (nonhU(i  ,1) == 0) then
                  dzb_e = .0_rKind
                  dzb_w = .5_rKind*(s%zb(i,1)-zbu(i-1,1))*ddxz(i)
               elseif (nonhU(i-1,1) == 0) then
                  dzb_e = .5_rKind*(zbu(i,1)-s%zb(i,1))  *ddxz(i)
                  dzb_w = .0_rKind
               endif

               if  (nonhV(i,1) == 1) then
                  dzb_s = .0_rKind
                  dzb_n = dzb_s
               elseif (nonhV(i  ,1  ) == 0) then
                  dzb_s = .0_rKind
                  dzb_n = .5_rKind*(s%zb(i,1)-zbv(i,1))*ddyz(1)
               endif

               awb(1,i,1) =  dzb_e * au(0,i  ,1  )+dzb_w*au(1,i-1,1  )  &           !main diagonal
               +  dzb_s * av(0,i  ,1  )+dzb_n*av(1,i  ,1)
               awb(2,i,1) =  dzb_w *  au(0,i-1,1  )                             !west
               awb(3,i,1) =  dzb_e *  au(1,i  ,1  )                             !east
               awb(4,i,1) =  dzb_n *  av(0,i  ,1  )                             !south
               awb(5,i,1) =  dzb_s *  av(1,i  ,1  )                             !north

               awbr(i,1)  =  dzb_e*s%uu(i,1)+dzb_w*s%uu(i-1,1  ) &
               +  dzb_s*s%vv(i,1)+dzb_n*s%vv(i  ,1  )
            else
               awb(:,i,1) = 0.0_rKind
               awbr(i,1)  = 0.0_rKind
            endif
         enddo
      else
         do j=2,s%ny
            do i=2,s%nx
               if (nonhZ(i,j) == 1) then
                  if     (nonhU(i,j)*nonhU(i-1,j) == 1) then
                     dzb_e = .5_rKind*(zbu(i,j)-zbu(i-1,j))*ddxz(i)
                     ! FB: this was dzs_e, which was not defined. Assumed it was a typo. TODO: check.
                     dzb_w = dzb_e
                  elseif (nonhU(i  ,j) == 0) then
                     dzb_e = .0_rKind
                     dzb_w = .5_rKind*(s%zb(i,j)-zbu(i-1,j))*ddxz(i)
                  elseif (nonhU(i-1,j) == 0) then
                     dzb_e = .5_rKind*(zbu(i,j)-s%zb(i,j))  *ddxz(i)
                     dzb_w = .0_rKind
                  endif

                  if     (nonhV(i,j)*nonhV(i,j-1) == 1) then
                     dzb_s = .5_rKind*(zbv(i,j)-zbv(i,j-1))*ddyz(j)
                     dzb_n = dzb_s
                  elseif (nonhV(i  ,j  ) == 0) then
                     dzb_s = .0_rKind
                     dzb_n = .5_rKind*(s%zb(i,j)-zbv(i,j-1))*ddyz(j)
                  elseif (nonhV(i  ,j-1) == 0) then
                     dzb_s = .5_rKind*(zbv(i,j)-s%zb(i,j))  *ddyz(j)
                     dzb_n = .0_rKind
                  endif

                  awb(1,i,j) =  dzb_e * au(0,i  ,j  )+dzb_w*au(1,i-1,j  )  &           !main diagonal
                  +  dzb_s * av(0,i  ,j  )+dzb_n*av(1,i  ,j-1)
                  awb(2,i,j) =  dzb_w *  au(0,i-1,j  )                             !west
                  awb(3,i,j) =  dzb_e *  au(1,i  ,j  )                             !east
                  awb(4,i,j) =  dzb_n *  av(0,i  ,j-1)                             !south
                  awb(5,i,j) =  dzb_s *  av(1,i  ,j  )                             !north

                  awbr(i,j)  =  dzb_e*s%uu(i,j)+dzb_w*s%uu(i-1,j  ) &
                  +  dzb_s*s%vv(i,j)+dzb_n*s%vv(i  ,j-1)
               else
                  awb(:,i,j) = 0.0_rKind
                  awbr(i,j)  = 0.0_rKind
               endif
            enddo
         enddo
      endif

      do j=jmin,jmax
         do i=2,s%nx
            if (nonhZ(i,j) == 1) then
               aws(1,i,j)   =  wcoef(i,j)*par%dt*2.0_rKind/s%hh(i,j)-awb(1,i,j)
               aws(2:5,i,j) =  -awb(2:5,i,j)
               awsr(i,j)    =  2.0_rKind * Wm(i,j)-awbr(i,j)
            else
               aws(:,i,j)   =  0.0_rKind
               awsr(i,j)    =  0.0_rKind
            endif
         enddo
      enddo

      !Substitute in the continuity equation
      if (sf1d) then
         do i=2,s%nx
            if (nonhZ(i,1)==1) then
               if     (nonhU(i,1)*nonhU(i-1,1) == 1) then
                  dzs_e = .5_rKind*(zsu(i,1)-zsu(i-1,1))*dyz(1)
                  dzs_w = dzs_e
               elseif (nonhU(i  ,1) == 0) then
                  dzs_e = .0_rKind
                  dzs_w = .5_rKind*(s%zs(i,1)-zsu(i-1,1))*dyz(1)
               elseif (nonhU(i-1,1) == 0) then
                  dzs_e = .5_rKind*(zsu(i,1)-s%zs(i,1))  *dyz(1)
                  dzs_w = .0_rKind
               endif

               ! wwvv problem the following construction make it impossible
               ! wwvv to determine if dzs_s gets a value
               ! wwvv (What other values than 0 and 1 can nonhV have)
               ! wwvv I put this kind of varibles to zero at the start
               ! wwvv of the subroutine
               if     (nonhV(i,1) == 1) then
                  dzs_s = .0_rKind
                  dzs_n = dzs_s
               elseif (nonhV(i  ,1  ) == 0) then
                  dzs_s = .0_rKind
                  dzs_n = .5_rKind*(s%zs(i,1)-zsv(i,1))*dxz(i)
               endif


               mat(1,i,1) =    dyz(1)*( s%hu(i,1)*au(0,i,1) - s%hu(i-1,1  )*au(1,i-1,  1) )  & !subs U left/right face
               +    dxz(i)*( s%hv(i,1)*av(0,i,1) - s%hv(i  ,1  )*av(1,i  ,  1) )  & !subs V left/rigth face
               -    dzs_e*au(0,i,1) - dzs_w*au(1,i-1,1)                           & !kin. boun. top
               -    dzs_s*av(0,i,1) - dzs_n*av(1,i-1,1)                           & !kin. boun. top
               +    dxz(i)*dyz(1)*aws(1,i,1)

               mat(2,i,1) = -  dyz(1)*s%hu(i-1,1)*au(0,i-1,  1)                              &
               -   dzs_w*au(0,i-1,1) + dxz(i)*dyz(1)*aws(2,i,1)

               mat(3,i,1) =    dyz(1)*s%hu(i,1)  * au(1,i,  1)                               &
               -    dzs_e*au(1,i  ,1) + dxz(i)*dyz(1)*aws(3,i,1)

               mat(4,i,1) = -  dxz(i)*s%hv(i,1)  *av(0,i,1)                              &
               -    dzs_n*av(0,i,1) + dxz(i)*dyz(1)*aws(4,i,1)

               mat(5,i,1) =    dxz(i)*s%hv(i,1)  * av(1,i,1)                                 &
               -    dzs_s*av(1,i,1)   + dxz(i)*dyz(1)*aws(5,i,1)

               rhs(i,1)   = -  dyz(1)*( s%hu(i,1)*aur(i,1) - s%hu(i-1,1  )*aur(i-1,  1) )    & !subs U left/right face
               -    dxz(i)*( s%hv(i,1)*avr(i,1) - s%hv(i  ,1  )*avr(i  ,  1) )    & !subs V left/rigth face
               +    dzs_e*aur(i,1) + dzs_w*aur(i-1,1)                             & !kin. boun. top
               +    dzs_s*avr(i,1) + dzs_n*avr(i,  1)                             & !kin. boun. top
               -    dxz(i)*dyz(1)*awsr(i,1)
            else
               mat(1  ,i,1) = 1.0_rKind
               mat(2:5,i,1) = 0.0_rKind
               rhs(i,1)     = 0.0_rKind
            endif
         enddo
      else
         do j=2,s%ny
            do i=2,s%nx
               if (nonhZ(i,j)==1) then
                  if     (nonhU(i,j)*nonhU(i-1,j) == 1) then
                     dzs_e = .5_rKind*(zsu(i,j)-zsu(i-1,j))*dyz(j)
                     dzs_w = dzs_e
                  elseif (nonhU(i  ,j) == 0) then
                     dzs_e = .0_rKind
                     dzs_w = .5_rKind*(s%zs(i,j)-zsu(i-1,j))*dyz(j)
                  elseif (nonhU(i-1,j) == 0) then
                     dzs_e = .5_rKind*(zsu(i,j)-s%zs(i,j))  *dyz(j)
                     dzs_w = .0_rKind
                  endif

                  if     (nonhV(i,j)*nonhV(i,j-1) == 1) then
                     dzs_s = .5_rKind*(zsv(i,j)-zsv(i,j-1))*dxz(i)
                     dzs_n = dzs_s
                  elseif (nonhV(i  ,j  ) == 0) then
                     dzs_s = .0_rKind
                     dzs_n = .5_rKind*(s%zs(i,j)-zsv(i,j-1))*dxz(i)
                  elseif (nonhV(i  ,j-1) == 0) then
                     dzs_s = .5_rKind*(zsv(i,j)-s%zs(i,j))  *dxz(i)
                     dzs_n = .0_rKind
                  endif


                  mat(1,i,j) =    dyz(j)*( s%hu(i,j)*au(0,i,j) - s%hu(i-1,j  )*au(1,i-1,  j) )  & !subs U left/right face
                  +    dxz(i)*( s%hv(i,j)*av(0,i,j) - s%hv(i  ,j-1)*av(1,i  ,j-1) )  & !subs V left/rigth face
                  -    dzs_e*au(0,i,j) - dzs_w*au(1,i-1,j)                           & !kin. boun. top
                  -    dzs_s*av(0,i,j) - dzs_n*av(1,i,j-1)                           & !kin. boun. top
                  +    dxz(i)*dyz(j)*aws(1,i,j)

                  mat(2,i,j) = -  dyz(j)*s%hu(i-1,j)*au(0,i-1,  j)                              &
                  -   dzs_w*au(0,i-1,j) + dxz(i)*dyz(j)*aws(2,i,j)

                  mat(3,i,j) =    dyz(j)*s%hu(i,j)  * au(1,i,  j)                               &
                  -    dzs_e*au(1,i  ,j) + dxz(i)*dyz(j)*aws(3,i,j)

                  mat(4,i,j) = -  dxz(i)*s%hv(i,j-1)  *av(0,i,j-1)                              &
                  -    dzs_n*av(0,i,j-1) + dxz(i)*dyz(j)*aws(4,i,j)

                  mat(5,i,j) =    dxz(i)*s%hv(i,j)  * av(1,i,j)                                 &
                  -    dzs_s*av(1,i,j)   + dxz(i)*dyz(j)*aws(5,i,j)

                  rhs(i,j)   = -  dyz(j)*( s%hu(i,j)*aur(i,j) - s%hu(i-1,j  )*aur(i-1,  j) )    & !subs U left/right face
                  -    dxz(i)*( s%hv(i,j)*avr(i,j) - s%hv(i  ,j-1)*avr(i  ,j-1) )    & !subs V left/rigth face
                  +    dzs_e*aur(i,j) + dzs_w*aur(i-1,j)                             & !kin. boun. top
                  +    dzs_s*avr(i,j) + dzs_n*avr(i,j-1)                             & !kin. boun. top
                  -    dxz(i)*dyz(j)*awsr(i,j)
               else
                  mat(1  ,i,j) = 1.0_rKind
                  mat(2:5,i,j) = 0.0_rKind
                  rhs(i,j)     = 0.0_rKind
               endif
            enddo
         enddo
      endif

      !Solve matrix
      !call timer_start(timer_flow_nonh_solv)
      dp = 0.0_rKind
      call solver_solvemat( mat  , rhs   , dp , s%nx, s%ny,par)
      !call timer_stop(timer_flow_nonh_solv)


      s%pres = s%pres + dp

      !Correct u/v/w

      !U
      do j=jmin,jmax
         do i=2,s%nx-1
            s%uu(i,j) = aur(i,j) + au(1,i,j)*dp(i+1,j)+au(0,i,j)*dp(i,j)
         enddo
      enddo

      !v
      if (sf1d) then
         do i=2,s%nx
            s%vv(i,1) = avr(i,1) + av(1,i,1)*dp(i,1)+av(0,i,1)*dp(i,1)
         enddo
      else
         do j=2,s%ny-1
            do i=2,s%nx
               s%vv(i,j) = avr(i,j) + av(1,i,j)*dp(i,j+1)+av(0,i,j)*dp(i,j)
            enddo
         enddo
      endif

      !W
      if (sf1d) then
         do i=2,s%nx
            if    (nonhZ(i,1) == 1) then
               s%ws(i,1) = awsr(i,1) + dp(i , 1)  * aws(1,i,1)                            &
               + dp(i-1,1)  * aws(2,i,1) + dp(i+1,1  ) * aws(3,i,1) &
               + dp(i  ,1)  * aws(4,i,1) + dp(i  ,  1) * aws(5,i,1)

            else
               !If the point is excluded in nh then get ws from the kinematic boundary condition
               if     (s%wetU(i,1)*s%wetU(i-1,1) == 1) then
                  dzs_e = .5_rKind*(zsu(i,1)-zsu(i-1,1))*ddxz(i)
                  dzs_w = dzs_e
               elseif (s%wetU(i  ,1) == 0) then
                  dzs_e = .0_rKind
                  dzs_w = .5_rKind*(s%zs(i,1)-zsu(i-1,1))*ddxz(i)
               elseif (s%wetU(i-1,1) == 0) then
                  dzs_e = .5_rKind*(zsu(i,1)-s%zs(i,1))  *ddxz(i)
                  dzs_w = .0_rKind
               else
                  dzs_e = .0_rKind
                  dzs_w = .0_rKind
               endif

               if     (s%wetV(i,1) == 1) then
                  dzs_s = .0_rKind
                  dzs_n = dzb_s
               elseif (s%wetV(i  ,1  ) == 0) then
                  dzs_s = .0_rKind
                  dzs_n = .5_rKind*(s%zs(i,1)-zsv(i,1))*ddyz(1)
               endif

               s%ws(i,1) = - (s%hu(i,1)*s%uu(i,1)-s%hu(i-1,1  )*s%uu(i-1,1  ))*ddxz(i) &
               - (s%hv(i,1)*s%vv(i,1)-s%hv(i  ,1  )*s%vv(i  ,1  ))*ddyz(1) &
               + dzs_e*s%uu(i,1)+dzs_w*s%uu(i-1,1)                         &
               + dzs_s*s%vv(i,1)+dzs_n*s%vv(i,  1)
            endif
         enddo
         do i=2,s%nx
            if (nonhZ(i,1) == 1) then
               s%wb(i,1) = awbr(i,1) + dp(i , 1)  * awb(1,i,1)                            &
               + dp(i-1,1)  * awb(2,i,1) + dp(i+1,1  ) * awb(3,i,1) &
               + dp(i  ,1)  * awb(4,i,1) + dp(i  ,1  ) * awb(5,i,1)
            else
               if     (s%wetU(i,1)*s%wetU(i-1,1) == 1) then
                  dzb_e = .5_rKind*(zbu(i,1)-zbu(i-1,1))*ddxz(i)
                  dzb_w = dzs_e
               elseif (s%wetU(i  ,1) == 0) then
                  dzb_e = .0_rKind
                  dzb_w = .5_rKind*(s%zb(i,1)-zbu(i-1,1))*ddxz(i)
               elseif (s%wetU(i-1,1) == 0) then
                  dzb_e = .5_rKind*(zbu(i,1)-s%zb(i,1))  *ddxz(i)
                  dzb_w = .0_rKind
               else
                  dzb_e = .0_rKind
                  dzb_w = .0_rKind
               endif

               if     (s%wetV(i,1) == 1) then
                  dzb_s = .0_rKind
                  dzb_n = dzb_s
               elseif (s%wetV(i  ,1  ) == 0) then
                  dzb_s = .0_rKind
                  dzb_n = .5_rKind*(s%zb(i,1)-zbv(i,1))*ddyz(1)
               endif

               s%wb(i,1) = + dzb_e*s%uu(i,1)+dzb_w*s%uu(i-1,1) &
               + dzb_s*s%vv(i,1)+dzb_n*s%vv(i,1)
            endif
         enddo
         !Assign boundaries
#ifdef USEMPI
         call xmpi_shift_ee(s%wb)
         call xmpi_shift_ee(s%ws)
         if (xmpi_istop) then
            s%wb(1,:) = s%wb(2,:)
         endif
         if (xmpi_isbot) then
            s%ws(s%nx+1,:) = s%ws(s%nx,:)
            s%wb(s%nx+1,:) = s%wb(s%nx,:)
         endif
#else
         s%ws(s%nx+1,:) = s%ws(s%nx,:)
         s%wb(1,:)      = s%wb(2,:)
         s%wb(s%nx+1,:) = s%wb(s%nx,:)
#endif
      else
         do j=2,s%ny
            do i=2,s%nx
               if    (nonhZ(i,j) == 1) then
                  s%ws(i,j) = awsr(i,j) + dp(i , j)  * aws(1,i,j)                            &
                  + dp(i-1,j)  * aws(2,i,j) + dp(i+1,j  ) * aws(3,i,j) &
                  + dp(i  ,j-1)* aws(4,i,j) + dp(i  ,j+1) * aws(5,i,j)

               else
                  !If the point is excluded in nh then get ws from the kinematic boundary condition
                  if     (s%wetU(i,j)*s%wetU(i-1,j) == 1) then
                     dzs_e = .5_rKind*(zsu(i,j)-zsu(i-1,j))*ddxz(i)
                     dzs_w = dzs_e
                  elseif (s%wetU(i  ,j) == 0) then
                     dzs_e = .0_rKind
                     dzs_w = .5_rKind*(s%zs(i,j)-zsu(i-1,j))*ddxz(i)
                  elseif (s%wetU(i-1,j) == 0) then
                     dzs_e = .5_rKind*(zsu(i,j)-s%zs(i,j))  *ddxz(i)
                     dzs_w = .0_rKind
                  else
                     dzs_e = .0_rKind
                     dzs_w = .0_rKind
                  endif

                  if     (s%wetV(i,j)*s%wetV(i,j-1) == 1) then
                     dzs_s = .5_rKind*(zsv(i,j)-zsv(i,j-1))*ddyz(j)
                     dzs_n = dzb_s
                  elseif (s%wetV(i  ,j  ) == 0) then
                     dzs_s = .0_rKind
                     dzs_n = .5_rKind*(s%zs(i,j)-zsv(i,j-1))*ddyz(j)
                  elseif (s%wetV(i  ,j-1) == 0) then
                     dzs_s = .5_rKind*(zsv(i,j)-s%zs(i,j))  *ddyz(j)
                     dzs_n = .0_rKind
                  else
                     dzs_e = .0_rKind
                     dzs_w = .0_rKind
                  endif

                  s%ws(i,j) = - (s%hu(i,j)*s%uu(i,j)-s%hu(i-1,j  )*s%uu(i-1,j  ))*ddxz(i) &
                  - (s%hv(i,j)*s%vv(i,j)-s%hv(i  ,j-1)*s%vv(i  ,j-1))*ddyz(j) &
                  + dzs_e*s%uu(i,j)+dzs_w*s%uu(i-1,j)                         &
                  + dzs_s*s%vv(i,j)+dzs_n*s%vv(i,j-1)
               endif
            enddo
         enddo

         do j=2,s%ny
            do i=2,s%nx
               if (nonhZ(i,j) == 1) then
                  s%wb(i,j) = awbr(i,j) + dp(i , j)  * awb(1,i,j)                            &
                  + dp(i-1,j)  * awb(2,i,j) + dp(i+1,j  ) * awb(3,i,j) &
                  + dp(i  ,j-1)* awb(4,i,j) + dp(i  ,j+1) * awb(5,i,j)
               else
                  if     (s%wetU(i,j)*s%wetU(i-1,j) == 1) then
                     dzb_e = .5_rKind*(zbu(i,j)-zbu(i-1,j))*ddxz(i)
                     dzb_w = dzs_e
                  elseif (s%wetU(i  ,j) == 0) then
                     dzb_e = .0_rKind
                     dzb_w = .5_rKind*(s%zb(i,j)-zbu(i-1,j))*ddxz(i)
                  elseif (s%wetU(i-1,j) == 0) then
                     dzb_e = .5_rKind*(zbu(i,j)-s%zb(i,j))  *ddxz(i)
                     dzb_w = .0_rKind
                  else
                     dzb_e = .0_rKind
                     dzb_w = .0_rKind
                  endif

                  if     (s%wetV(i,j)*s%wetV(i,j-1) == 1) then
                     dzb_s = .5_rKind*(zbv(i,j)-zbv(i,j-1))*ddyz(j)
                     dzb_n = dzb_s
                  elseif (s%wetV(i  ,j  ) == 0) then
                     dzb_s = .0_rKind
                     dzb_n = .5_rKind*(s%zb(i,j)-zbv(i,j-1))*ddyz(j)
                  elseif (s%wetV(i  ,j-1) == 0) then
                     dzb_s = .5_rKind*(zbv(i,j)-s%zb(i,j))  *ddyz(j)
                     dzb_n = .0_rKind
                  else
                     dzb_e = .0_rKind
                     dzb_w = .0_rKind
                  endif

                  s%wb(i,j) = + dzb_e*s%uu(i,j)+dzb_w*s%uu(i-1,j) &
                  + dzb_s*s%vv(i,j)+dzb_n*s%vv(i,j-1)
               endif
            enddo
         enddo
#ifdef USEMPI
         call xmpi_shift_ee(s%wb)
         call xmpi_shift_ee(s%ws)
         if (xmpi_istop) then
            s%wb(1,:) = s%wb(2,:)
         endif
         if (xmpi_isbot) then
            s%ws(s%nx+1,:) = s%ws(s%nx,:)
            s%wb(s%nx+1,:) = s%wb(s%nx,:)
         endif
         if (xmpi_isleft) then
            s%ws(:,1)      = s%ws(:,2)
            s%wb(:,1)      = s%wb(:,2)
         endif
         if (xmpi_isright) then
            s%wb(:,s%ny+1) = s%wb(:,s%ny)
            s%ws(:,s%ny+1) = s%ws(:,s%ny)
         endif
#else
         !Assign boundaries
         s%ws(:,1)      = s%ws(:,2)
         s%ws(:,s%ny+1) = s%ws(:,s%ny)
         !s%ws(1,:)      = s%ws(2,:)
         s%ws(s%nx+1,:) = s%ws(s%nx,:)

         s%wb(:,1)      = s%wb(:,2)
         s%wb(:,s%ny+1) = s%wb(:,s%ny)
         s%wb(1,:)      = s%wb(2,:)
         s%wb(s%nx+1,:) = s%wb(s%nx,:)
#endif
      endif

      Wm_old = .5_rKind*(s%ws+s%wb)
   end subroutine nonh_1lay_cor

   subroutine nonh_1lay_pred(s,par)
      !==============================================================================
      !

      ! DATE                AUTHOR               CHANGES
      !
      ! November  2010       Pieter Bart Smit     New Subroutine


      !------------------------------------------------------------------------------
      !                             DECLARATIONS
      !------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------
      !
      ! Include the pressure explicitly in the predictor step
      !
      !--------------------------     DEPENDENCIES       ----------------------------
      use spaceparams
      use params
      use flow_secondorder_module
      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars) ,intent(inout)                       :: s
      type(parameters),intent(in)                          :: par

      !--------------------------     LOCAL VARIABLES    ----------------------------

      !Indices
      integer(kind=iKind)                     :: i,ie,iee,iw               !Index variables
      integer(kind=iKind)                     :: j,js
      integer(kind=iKind)                     :: jmin,jmax

      real(kind=rKind)                        :: dwdx1    !Gradient of vertical velocity in x-dir at i+1/2,j
      real(kind=rKind)                        :: dwdx2    !Gradient of vertical velocity in x-dir at i-1/2,j
      real(kind=rKind)                        :: dwdy1    !Gradient of vertical velocity in x-dir at i    ,j+1/2
      real(kind=rKind)                        :: dwdy2    !Gradient of vertical velocity in x-dir at i    ,j-1/2
      real(kind=rKind)                        :: Vol
      real(kind=rKind)                        :: wmax,reformfac

      if (.not. initialized) then
         call nonh_init(s,par)
      endif

      if (s%ny>0) then
         jmin = 2
         jmax = s%ny
      else
         jmin = 1
         jmax = 1
      endif

      !
      ! Determine if a velocity point will be included in the nonh pressure matrix, do not include if:
      !
      ! (1)  The point is dry
      ! (2)  The relative wave length kd of the smallest possible wave (L=2dx) is smaller than kdmin
      ! (3)  The interpolated waterlevel zs is below the bottom (steep cliffs with overwash situations)
      ! (4)  Where Miche breaker criterium applies -> bores are hydrostatic
      !            max steepness = H/L = maxbrsteep
      !            dz/dx = maxbrsteep
      !            dz/dx = dz/dt/c = w/c = w/sqrt(gh)
      !            wmax = maxbrsteep*sqrt(gh)

      do j=1,s%ny+1
         do i=1,s%nx+1
            iw = max(i,i-1)
            ie = min(s%nx,i+1)
            iee = min(s%nx,i+2)

            if (  (s%wetU(i,j)==1                                      )  &
            .and. (0.5_rKind*(s%zs(i,j) + s%zs(ie,j))    > zbu(i,j)    )  &
            .and. ( s%dsu(i,1)*par%kdmin/par%px  < s%hum(i,j)  )  ) then
               nonhU(i,j) = 1
            else
               nonhU(i,j) = 0
            endif
         enddo
      enddo

      if (s%ny>2) then
         do j=1,s%ny+1
            js = min(s%ny,j+1)
            do i=1,s%nx+1
               if (  (s%wetV(i,j)==1                                      )  &
               .and. (0.5_rKind*(s%zs(i,j) + s%zs(i,js))    > zbv(i,j)    )  &
               .and. ( s%dnv(1,j)*par%kdmin/par%px  < s%hvm(i,j)  )  ) then
                  nonhV(i,j) = 1
               else
                  nonhV(i,j) = 0
               endif
            enddo
         enddo
      else
         nonhV = 0
      endif
      !
      ! Determine if a velocity point will be included in the nonh pressure matrix, include if
      ! any of the surrounding velocity points is included.
      !
      if (s%ny>0) then
         do j=2,s%ny
            do i=2,s%nx
               if (max(nonhV(i,j),nonhV(i,j-1),nonhU(i,j),nonhU(i-1,j)) > 0) then
                  nonhZ(i,j) = 1
               else
                  nonhZ(i,j) = 0
               endif
            enddo
         enddo
      else
         do i=2,s%nx
            if (max(nonhU(i,1),nonhU(i-1,1)) > 0) then
               nonhZ(i,1) = 1
            else
               nonhZ(i,1) = 0
            endif
         enddo
      endif

      if (par%nhbreaker == 1) then
         reformfac = par%reformsteep/par%maxbrsteep
         if (s%ny>0) then
            do j=2,s%ny
               do i=2,s%nx
                  wmax = par%maxbrsteep*sqrt(par%g*s%hh(i,j))
                  if (s%breaking(i,j) == 0) then
                     if (s%ws(i,j)>=wmax) then
                        s%breaking(i,j) = 1
                     elseif (s%ws(i,j)<=-wmax) then
                        s%breaking(i,j) = -1
                     endif
                  elseif (s%breaking(i,j)==1) then
                     if (s%ws(i,j)<reformfac*wmax) then
                        s%breaking(i,j) = 0
                     endif
                  elseif (s%breaking(i,j)==-1) then
                     if (s%ws(i,j)>reformfac*(-wmax)) then
                        s%breaking(i,j) = 0
                     endif
                  endif
               enddo
            enddo
         else
            do i=2,s%nx
               wmax = par%maxbrsteep*sqrt(par%g*s%hh(i,1))
               if (s%breaking(i,1) == 0) then
                  if (s%ws(i,1)>=wmax) then
                     s%breaking(i,1) = 1
                  elseif (s%ws(i,1)<=-wmax) then
                     s%breaking(i,1) = -1
                  endif
               elseif (s%breaking(i,1)==1) then
                  if (s%ws(i,1)<reformfac*wmax) then
                     s%breaking(i,1) = 0
                  endif
               elseif (s%breaking(i,1)==-1) then
                  if (s%ws(i,1)>reformfac*(-wmax)) then
                     s%breaking(i,1) = 0
                  endif
               endif
            enddo
         endif
         ! turn off non-hydrostatic pressure correction in areas with breaking and increase viscosity
         where (s%breaking/=0)
            nonhZ = 0
            s%pres = 0
         endwhere
      elseif (par%nhbreaker == 2) then
         ! First determine local breaker criterion
         lbreakcond = par%maxbrsteep
         if (s%ny==0) then
            do i=2,s%nx
               if (s%breaking(i,1)==1) then
                  lbreakcond(i-1:i+1,1) = par%secbrsteep
               endif
            enddo
         else
            do j=jmin,jmax
               do i=2,s%nx
                  if (s%breaking(i,j)==1) then
                     lbreakcond(i-1:i+1,j-1:j+1) = par%secbrsteep
                  endif
               enddo
            enddo
         endif
#ifdef USEMPI
         call xmpi_shift_ee(lbreakcond)
#endif
         ! Now find areas where main breaking criterion is exceeded
         do j=jmin,jmax
            do i=2,s%nx
               if (s%breaking(i,j)==1) then
                  if(s%ws(i,j)<=0.d0) then
                     s%breaking(i,j) = 0
                  endif
               else
                  wmax = lbreakcond(i,j)*sqrt(par%g*s%hh(i,j))  ! add current term in here too
                  if (s%ws(i,j)>=wmax) then
                     s%breaking(i,j) = 1
                  endif
               endif
            enddo
         enddo
         ! turn off non-hydrostatic pressure correction in areas with breaking and increase viscosity
         do j=jmin,jmax
            do i=2,s%nx
               if (s%breaking(i,j)==1) then
                  nonhZ(i,j) = 0
                  s%pres(i,j) = 0.d0
               endif
            enddo
         enddo
      elseif (par%nhbreaker == 3) then
         ! First determine local breaker criterion
         lbreakcond = par%maxbrsteep
         if (s%ny==0) then
            do i=2,s%nx
               if (s%breaking(i,1)==1) then
                  lbreakcond(i-1:i+1,1) = par%secbrsteep
               endif
            enddo
         else
            do j=jmin,jmax
               do i=2,s%nx
                  if (s%breaking(i,j)==1) then
                     lbreakcond(i-1:i+1,j-1:j+1) = par%secbrsteep
                  endif
               enddo
            enddo
         endif
         ! Now find areas where main breaking criterion is exceeded
         do j=jmin,jmax
            do i=2,s%nx
               s%wscrit(i,j) = lbreakcond(i,j)*sqrt(par%g*s%hh(i,j))  ! add current term in here too
               if (s%breaking(i,j)==1) then
                  if(s%ws(i,j)<=s%wscrit(i,j)) then
                     s%breaking(i,j) = 0
                  endif
               else
                  if (s%ws(i,j)>=s%wscrit(i,j)) then
                     s%breaking(i,j) = 1
                  endif
               endif
            enddo
         enddo
         ! turn off non-hydrostatic pressure correction in areas with breaking and increase viscosity
         do j=jmin,jmax
            do i=2,s%nx
               if (s%breaking(i,j)==1) then
                  nonhZ(i,j) = 0
                  s%pres(i,j) = 0.d0
               endif
            enddo
         enddo
      endif

      !Calculate explicit part average vertical momentum (advection)
      if (s%ny>0) then
         do j=2,s%ny
            do i=2,s%nx
               if (nonhZ(i,j) == 1) then
                  Wm(i,j) = Wm_old(i,j) - par%dt*( ddxu(i-1)*max(s%qx(i-1,j  ),0.0_rKind)*(Wm_old(i  ,j  ) &
                  - Wm_old(i-1,j  ))/s%hh(i,j)   &
                  + ddxu(i)  *min(s%qx(i  ,j  ),0.0_rKind)*(Wm_old(i+1,j  )-Wm_old(i  ,j  ))/s%hh(i,j)   &
                  + ddyv(j-1)*max(s%qy(i  ,j-1),0.0_rKind)*(Wm_old(i  ,j  )-Wm_old(i  ,j-1))/s%hh(i,j)   &
                  + ddyv(j  )*min(s%qy(i  ,j  ),0.0_rKind)*(Wm_old(i  ,j+1)-Wm_old(i  ,j  ))/s%hh(i,j) )
               else
                  Wm(i,j) = 0.0_rKind
                  Wm_old(i,j) = 0.0_rKind
                  s%ws(i,j)   = 0.0_rKind
                  s%wb(i,j)   = 0.0_rKind
                  s%pres(i,j) = 0.0_rKind
               endif
            enddo
         enddo
      else
         do i=2,s%nx
            if (nonhZ(i,1) == 1) then
               Wm(i,1) = Wm_old(i,1) - par%dt*( ddxu(i-1)*max(s%qx(i-1,1  ),0.0_rKind)*(Wm_old(i  ,1  ) &
               -Wm_old(i-1,1  ))/s%hh(i,1)   &
               + ddxu(i)  *min(s%qx(i  ,1  ),0.0_rKind)*(Wm_old(i+1,1  )-Wm_old(i  ,1  ))/s%hh(i,1) )
            else
               Wm(i,1) = 0.0_rKind
               Wm_old(i,1) = 0.0_rKind
               s%ws(i,1)   = 0.0_rKind
               s%wb(i,1)   = 0.0_rKind
               s%pres(i,1) = 0.0_rKind
            endif
         enddo
      endif

      !Calculate explicit part vertical viscocity
      if (s%ny>0) then
         do j=2,s%ny
            do i=2,s%nx
               dwdx1 = .5d0*(s%nuh(i-1,j  )+s%nuh(i,j  ))*s%hu(i-1,j  )*(Wm_old(i  ,j  )-Wm_old(i-1,j  ))*ddxu(i-1)
               dwdx2 = .5d0*(s%nuh(i+1,j  )+s%nuh(i,j  ))*s%hu(i  ,j  )*(Wm_old(i+1,j  )-Wm_old(i  ,j  ))*ddxu(i)
               dwdy1 = s%nuh(i  ,j-1)*s%hu(i  ,j-1)*(Wm_old(i  ,j  )-Wm_old(i  ,j-1))*ddyv(j-1)
               dwdy2 = s%nuh(i  ,j  )*s%hu(i  ,j  )*(Wm_old(i  ,j+1)-Wm_old(i  ,j  ))*ddyv(j)
               Wm(i,j) = Wm(i,j)   + (1.0_rKind/s%hh(i,j))*par%dt*(dwdx2-dwdx1)*ddxz(i)*real(s%wetu(i,j)*s%wetu(i-1,j),rKind) &
               + (1.0_rKind/s%hh(i,j))*par%dt*(dwdy2-dwdy1)*ddyz(j)*real(s%wetv(i,j)*s%wetv(i,j-1),rKind)
            enddo
         enddo
      else
         do i=2,s%nx
            dwdx1 = .5d0*(s%nuh(i-1,1  )+s%nuh(i,1  ))*s%hu(i-1,1  )*(Wm_old(i  ,1  )-Wm_old(i-1,1  ))*ddxu(i-1)
            dwdx2 = .5d0*(s%nuh(i+1,1  )+s%nuh(i,1  ))*s%hu(i  ,1  )*(Wm_old(i+1,1  )-Wm_old(i  ,1  ))*ddxu(i)
            dwdy1 = 0.d0
            dwdy2 = 0.d0
            Wm(i,1) = Wm(i,1)   + (1.0_rKind/s%hh(i,1))*par%dt*(dwdx2-dwdx1)*ddxz(i)*real(s%wetu(i,1)*s%wetu(i-1,1),rKind)
         enddo
      endif

      do j=jmin,jmax
         do i=2,s%nx
            if (nonhU(i,j)==1) then
               vol       = 0.5_rKind*par%dt/(s%hum(i,j)*dxu(i))
               au(1,i,j) = - (s%zs(i+1,j) - s%zb(i  ,j))*vol
               au(0,i,j) = + (s%zs(i  ,j) - s%zb(i+1,j))*vol
            else
               au(1,i,j) =  0.0_rKind
               au(0,i,j) =  0.0_rKind
            endif
         enddo
      enddo
      if (xmpi_istop) then               !wwvv: added test on istop
         au(:,1,:)      = 0.0_rKind
      endif
      if (xmpi_isbot) then               !wwvv: added test on isbot
         au(:,s%nx+1,:)   = 0.0_rKind
      endif


      !Built pressure coefficients V
      !call timer_start(timer_flow_nonh_av)
      if (s%ny>2) then
         do j=2,s%ny
            do i=2,s%nx
               if (nonhV(i,j)==1)then
                  vol       = 0.5_rKind*par%dt/(s%hvm(i,j)*dyv(j))
                  av(1,i,j)  = -(s%zs(i  ,j+1) - s%zb(i  ,j  ))*vol
                  av(0,i,j)  = +(s%zs(i  ,j  ) - s%zb(i  ,j+1))*vol
                  avr(i,j)   = s%vv(i,j)
               else
                  av(1,i,j) =  0.0_rKind
                  av(0,i,j) =  0.0_rKind
               endif
            enddo
         enddo
         if (xmpi_isleft) then               !wwvv: added test on isleft
            av(:,:,1)     = 0.0_rKind
         endif
         if (xmpi_isright) then              !wwvv: added test on isright
            av(:,:,s%ny+1)  = 0.0_rKind
         endif
      endif

      !Include explicit approximation for pressure in s%uu and s%vv   and Wm
      if (par%secorder == 1) then
         do j=jmin,jmax
            do i=2,s%nx-1
               s%uu(i,j) = s%uu(i,j) + au(1,i,j) * s%pres(i+1,j) + au(0,i,j) * s%pres(i  ,j)
            enddo
         enddo

         if (s%ny>0) then
            do j=2,s%ny-1
               do i=2,s%nx
                  s%vv(i,j) = s%vv(i,j) + av(1,i,j) * s%pres(i,j+1) + av(0,i,j) * s%pres(i,j  )
               enddo
            enddo
         else
            do i=2,s%nx
               s%vv(i,1) = s%vv(i,1) + av(1,i,1) * s%pres(i,1) + av(0,i,1) * s%pres(i,1  )
            enddo
         endif

         do j=jmin,jmax
            do i=2,s%nx
               if (nonhZ(i,j) == 1) then
                  Wm(i,j) = Wm(i,j)   + Wcoef(i,j)*par%dt * s%pres(i,j)/s%hh(i,j)
               endif
            enddo
         enddo
         call flow_secondorder_advW(s,par,Wm,Wm_old, nonhZ)
      endif

   end subroutine nonh_1lay_pred




   !==============================================================================
   !                     REDUCED NON-HYDROSTATIC ROUTINES 2DV
   !==============================================================================
   !
   !
   !==============================================================================
   subroutine nonh_2lay_cor_2dV(s,par)
      !==============================================================================
      !

      ! DATE               AUTHOR               CHANGES
      !
      ! October   2014     Pieter Bart Smit     New subroutine
      !
      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------

      !  In this subroutine we calculate the nonh-hydrostatic pressures by solving
      !  a discrete poisson equation, and subsequently correct the velocity field
      !  so that the velocities at the new timelevel are divergence free.
      !
      !  This subroutine uses the pressure coeficients calculated previously in
      !  the routine nonh_2lay_pred_2dV, and for each of the velocity components all
      !  processes have been accounted for (advection, bed-friction etc), including
      !  an explict estimation of the nonh press known from the previous timestep
      !  (unless we are calculating with first order accuracy).
      !
      !  Subsequently we will take the following steps:
      !
      !  1) Calculate the water and bed levels at u-points (subroutine call to zuzv).
      !  2) Built the poisson matrix by substituting the pressure dependencies of the
      !     momentum equations into the continuity equation.
      !  3) Solve the resulting system (call to solver_solvemat).
      !  4) Update the non-hydrostatic pressure. If second-order is active, dp is the
      !     difference between the old and the new pressure, so in that case we have
      !
      !     s%press = s%press + dp
      !
      !     else, it represents the new pressure directly, s%press = dp.
      !  5) Correct the velocity components with the nonhydrostatic pressures. The
      !     new velocity field will be divergence free.
      !
      !  After these steps, control is returned to flow_timestep.
      !
      !--------------------------     DEPENDENCIES       ----------------------------
      !
      use spaceparams
      use params
      use solver_module
      use xmpi_module
      !
      !--------------------------     ARGUMENTS          ----------------------------
      !
      type(spacepars) ,intent(inout)                       :: s
      type(parameters),intent(in)                          :: par
      !
      !--------------------------     LOCAL PARAMETERS   ----------------------------
      !
      integer(kind=iKind),parameter           :: j = 1           !Index variables
      !
      !--------------------------     LOCAL VARIABLES    ----------------------------
      !

      !Indices
      integer(kind=iKind)                     :: i,k !Index variables

      !Work
      real(kind=rKind)                        :: dzsdx        ! Free surface gradient (s-dir)
      real(kind=rKind)                        :: dzbdx        ! Bed level gradient    (s-dir)
      real(kind=rKind)                        :: alpha        ! Parameter controlling how the velocity at the
      ! layer interface is determined ( =-wcoef for energy
      ! conservative behaviour).
      real(kind=rKind)                        :: alpha_u(0:1) ! The layer distribution parameter interpolated at
      ! the east (=1) and west (=0) faces.
      real(kind=rKind)                        :: hfU(0:1)     ! Coefficient of the west U(i-1,j) (=0) or east U(i,j) (=1) velocity
      ! in the reduced continuity equation.
      real(kind=rKind)                        :: hfdU(0:1)    ! Coefficient of the west dU(i-1,j) (=0) or east dU(i,j) (=1) velocity
      !
      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------
      !

      !
      ! ------ Substitute in the continuity equation ------
      !
      do i=2,s%nx
         !
         if (nonhZ(i,j)==1) then
            !
            do k = 0,1
               !
               alpha_u(k) =  .5d0 * (wcoef(i+k,j) + wcoef(i-1+k,j))

               hfU(k)  = s%hu(i-1+k  ,j) * ( 1.0d0 + alpha_u(k) )
               hfdU(k) = s%hu(i-1+k  ,j) * alpha_u(k) * ( 1.0d0 - alpha_u(k) )
               !
            enddo
            !
            dzsdx = .5*( zsu(i,j) - zsu(i-1,j) )
            dzbdx = .5*( zbu(i,j) + alpha_u(1) * s%hu(i,j) - zbu(i-1,j) - alpha_u(0) * s%hu(i-1,j) )
            !


            !
            alpha = -wcoef(i,j)
            !
            mat(1,i,j) = 2.0d0* s%dsz(i,j) * aws(1,i,j) &
            + hfU( 1)*au( 0,i,j) - hfU( 0)*au(1, i-1,j)  &
            + hfdU(1)*adU(0,i,j) - hfdU(0)*adU(1,i-1,j)  &
            - dzsdx * ( au(0,i,j) + au(1,i-1,j) - ( adU(0,i,j) + adU(1,i-1,j) )* wcoef(i,j) )    &
            - dzbdx * ( au(0,i,j) + au(1,i-1,j) + ( adU(0,i,j) + adU(1,i-1,j) )* alpha )
            !
            mat(2,i,j) = - hfU( 0)*au (0, i-1,j) - hfdU( 0)*adU(0, i-1,j)      &
            - dzsdx * ( au(0,i-1,j) - adU(0,i-1,j) * wcoef(i,j) ) &
            - dzbdx * ( au(0,i-1,j) + adU(0,i-1,j) * alpha )
            !
            mat(3,i,j) =   hfU(1)*au (1,i  ,j) + hfdU(1)*adU(1,i,j)      &
            -    dzsdx * ( au(1,i,j)  - adU(1,i,j) * wcoef(i,j) ) &
            -    dzbdx * ( au(1,i,j)  + adU(1,i,j) * alpha )
            !
            rhs(i,j)   = - 2.0d0 * s%dsz(i,j) * Wm(i,j) &
            - hfU( 1)*s%uu( i  ,j) + hfU( 0)*s%uu( i-1 ,j) &
            - hfdU(1)*s%dU( i  ,j) + hfdU(0)*s%dU( i-1 ,j) &
            + dzsdx * ( s%uu(i,j) + s%uu(i-1,j) - ( s%dU(i,j) + s%dU(i-1,j) )* wcoef(i,j) ) &
            + dzbdx * ( s%uu(i,j) + s%uu(i-1,j) + ( s%dU(i,j) + s%dU(i-1,j) )* alpha )
         else
            !
            mat(1  ,i,j) = 1.0_rKind
            mat(2:5,i,j) = 0.0_rKind
            rhs(i,j)     = 0.0_rKind
            !
         endif
         !
      enddo

      !
      !---------------- Solve Linear System -----------------
      !
      dp = 0.0_rKind
      call solver_solvemat( mat  , rhs   , dp , s%nx, s%ny,par)
      !
      if (par%secorder == 1) then
         !
         s%pres = s%pres + dp
         !
      else
         !
         s%pres = dp
         !
      endif

      !
      !------------------ Correct u/v/w -------------------
      !
      !U
      do i=2,s%nx-1
         !
         s%uu(i,j) = s%uu(i,j) + au(1,i,j)*dp(i+1,j)+au(0,i,j)*dp(i,j)
         !
      enddo

      !dU
      if (par%nhlay > 0. ) then
         !
         do i=2,s%nx-1
            !
            s%dU(i,j) = s%dU(i,j) + adU(1,i,j)*dp(i+1,j)+adU(0,i,j)*dp(i,j)
            !
         enddo
         !
      endif
      !

#ifdef USEMPI
      ! Communicate velocities
      call xmpi_shift_ee(s%uu)
      call xmpi_shift_ee(s%dU)
#endif

      !
      ! Calculate the surface vertical velocity from the kinematic condition
      !
      !
      do i=2,s%nx
         !
         if ( s%wetz(i,j)==1 ) then
            !
            dzsdx = .5*( zsu(i,j) - zsu(i-1,j  ) )  / s%dsz(i,j)

            s%ws(i,j) = - (s%hu(i,j)*s%uu(i,j) -s%hu(i-1,j  )*s%uu(i-1,j  )) / s%dsz(i,j) &
            +  dzsdx * (  s%uu(i,j) + s%uu(i-1,j  ) - ( s%dU(i,j) + s%dU(i-1,j  ) ) * wcoef(i,j) )
            !
         else
            !
            s%ws(i,j) = 0.0d0
            !
         endif
         !
      enddo
      !

      !
      ! Calculate the bottom velocity from the kinematic condition
      !
      !
      do i=2,s%nx
         !
         if ( s%wetz(i,j)==1 ) then
            !
            dzbdx = .5*( zbu(i,j)  - zbu(i-1,j  ) ) / s%dsz(i,j)
            !
            s%wb(i,j) = +  dzbdx * (  s%uu(i,j) + s%uu(i-1,j  ) + ( s%dU(i,j) + s%dU(i-1,j  ) ) * (1-wcoef(i,j)) )
            !
         else
            !
            s%wb(i,j) = 0.0d0
            !
         endif
         !
      enddo
      !

      !
      do i=2,s%nx
         !
         !
         if ( nonhZ(i,j) == 1) then
            !
            Wm(i,j) = Wm(i,j) + aws(1,i,j)*dp(i,j)
            !
         else
            !
            Wm(i,j) = 0.0d0
            !
         endif
         !
      enddo
      !

      if (par%nhlay > 0. ) then
         !

         do i=2,s%nx
            !
            !
            if ( nonhZ(i,j) == 1 ) then
               !
               Wm0(i,j) = Wm(i,j) + 0.5d0* ( s%wb(i,j) - s%ws(i,j) )
               !
            else
               !
               Wm0(i,j) = 0.
               !
            endif
            !
         enddo
         !
      endif



#ifdef USEMPI
      call xmpi_shift_ee(Wm)
      call xmpi_shift_ee(s%dU)

      if (xmpi_istop) then
         Wm(1,:) = Wm(2,:)
      endif
      if (xmpi_isbot) then
         Wm(s%nx+1,:) = Wm(s%nx,:)
      endif
#else
      !Assign boundaries
      Wm(1,:) = Wm(2,:)
      Wm(s%nx+1,:) = Wm(s%nx,:)

#endif

      !
      ! Calculate the relative vertical advection
      !
      do i=2,s%nx
         !
         if (s%wetz(i,j) == 1) then
            !
            omega(i,j) = - ( s%hu(i  ,j) * dU0(i  ,j) - s%hu(i-1,j) * dU0(i-1,j) ) / s%dsz( i , j )
            !
         else
            !
            omega(i,j) = 0.0d0
            !
         endif
         !
      enddo

      !
      Wm_old = Wm
      !s%ws   = s%dU
      dU0    = s%dU
      !1111111
   end subroutine nonh_2lay_cor_2dv

   !
   subroutine nonh_2lay_pred_2dV(s,par,uu0)
      !==============================================================================
      !

      ! DATE               AUTHOR               CHANGES
      !
      ! October  2014       Pieter Bart Smit     New Subroutine


      !------------------------------------------------------------------------------
      !                             DECLARATIONS
      !------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------
      !
      ! Include the pressure explicitly in the predictor step. In this subroutine
      ! the following steps are taken:
      !
      ! 1) Calculate whether or not breaking is active (handled in subroutine
      !    nonh_masks ).
      ! 2) Calculate the coeficients of the nonh. pres. as they occur in the momentum
      !    equations. (these are use to calculate the explicit pressure contributions
      !    here, and the implicit contributions in the correction routine later).
      ! 3) Calculate the first order additional contributions to the advection due to
      !    the reduced two-layer formulation (if applicable).
      ! 4) Include an explicit estimate for the nonh press. in the momentum equations
      !    using the earlier calculated coef. (if second-order approximations for the
      !    advection are used, else this step is not necessary).
      ! 5) Include the contributions due to internal stress terms and bottom friction
      !    in the equation for the velocity difference (implicitly for stability).
      !
      ! After these steps, control flows back to flow_timestep, and - if active -
      ! second-order approximations are calculated. Consequently, flow_timestep calls
      ! nonh_cor again and the nonhydrostatic computations are continued in
      ! nonh_2lay_cor.
      !
      !--------------------------     DEPENDENCIES       ----------------------------
      use spaceparams
      use params
      use flow_secondorder_module
      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars) ,intent(inout)                       :: s
      type(parameters),intent(in)                          :: par
      real(kind=rkind),dimension(s%nx+1,s%ny+1),intent(in) :: uu0 !The velocity at the previous timestep

      !--------------------------     LOCAL PARAMETERS   ----------------------------

      integer(kind=iKind), parameter  :: j = 1        !Y-index, constant in 1d

      !--------------------------     LOCAL VARIABLES    ----------------------------
      !General
      real(kind=rKind)                :: fac          !A factor

      !Indices
      integer(kind=iKind)             :: i,k          !Loop indices

      !Used in calculation of pressure gradients
      real(kind=rKind)                :: Vol          !The volume of a momentum cell

      !Used in calculation of advection
      real(kind=rKind)                :: qx(0:1)      !Discharge through the faces of a momentum cell (0: left face
      !1: right face)
      real(kind=rKind)                :: faceval(0:1) !Value of dU at the faces of a momentum cell (0: left face
      !1: right face)

      !Used in calculation of vertical exchange of momentum
      real(kind=rKind)                :: fric         !Bottom friction coeficient (= dt * cf * |U| / h)
      real(kind=rKind)                :: stress       !Coeficient used in calculation of stress
      real(kind=rKind)                :: ub           !The (Eulerian) velocity at the bottom
      real(kind=rKind)                :: umag         !The (Eulerian) magnitude of the velocity vector at the bottom
      real(kind=rKind)                :: ve           !The eddy viscosity used for vertical mixing
      real(kind=rKind)                :: advec        !Coeficient used in calculation of vertical advection
      real(kind=rKind)                :: Amat,rhs     !The 1*1 "matrix" and RHS of the vertical "system".

      !------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !------------------------------------------------------------------------------

      !
      ! Determine whether or not to use the nh correction method (dry/wet/breaking etc.)
      !
      call nonh_break( s , par )
      call zuzv(s)
      !
      ! -- Calculate the pressure dependencies for U --
      !
      do i=2,s%nx
         !
         if (nonhU(i,j)==1) then
            !
            vol       = 0.5d0 * par%dt/(s%hum(i,j)*s%dsu(i,j))
            fac       = (s%zb(i+1,j) - s%zb(i  ,j)) * vol
            !
            au(1,i,j) = - vol * (1.d0 + wcoef(i+1,j) ) * s%hh(i+1,j) - fac
            au(0,i,j) = + vol * (1.d0 + wcoef(i  ,j) ) * s%hh(i,j)   - fac
            !
         else
            !
            au(1,i,j) =  0.0_rKind
            au(0,i,j) =  0.0_rKind
            !
         endif
         !
      enddo

      !
      ! -- Calculate the pressure dependencies for dU --
      !

      if ( par%nhlay > 0. ) then
         !
         do i=2,s%nx
            !
            if (nonhU(i,j)*(1-s%breaking(i,j))*(1-s%breaking(i-1,j))==1 ) then
               !
               vol       = 0.5_rKind*par%dt/(s%hum(i,j)*s%dsu(i,j))
               !
               fac = (  s%zs( i+1,j ) -  s%zs(i,j) ) * vol / ( 1.0d0 - .5d0*wcoef(i,j) - .5d0*wcoef(i+1,j) ) &
               - ( wcoef( i+1,j ) - wcoef(i,j) ) * s%hum(i,j) * vol / (2.0d0 - wcoef(i,j) - wcoef(i+1,j) )
               !
               adU(1,i,j) = - s%hh(i+1,j)*vol + fac
               adU(0,i,j) = + s%hh(i,j)*vol   + fac
               !
            else
               !
               adU(1,i,j) =  0.0_rKind
               adU(0,i,j) =  0.0_rKind
               !
            endif
            !
         enddo
         !
      endif

      !
      ! -- Calculate the pressure dependencies for Wm --
      !

      do i=2,s%nx
         !

         if (nonhZ(i,j) == 1) then
            !
            aws(1,i,j)   = par%dt/( s%hh(i,j) * ( 1.0d0 - wcoef( i , j ) ) )
            !
         else
            !
            aws(1,i,j)   =  0.0_rKind
            !
         endif
      enddo
      !
      !
      !----------------- Advection Wm (FOU) -------------------
      !
      do i=2,s%nx
         !
         if ( nonhZ(i,j) == 1) then
            !

            do k = 0,1
               !
               qx(k) = ( s%qx(i - 1 + k , j ) )

               if (qx(k) > 0.) then
                  !
                  faceval(k) = Wm_old( i - 1 + k ,j  ) * real( nonhZ(i - 1 + k,j) , kind=rKind )
                  !
               else
                  !
                  faceval(k) = Wm_old( i + k ,j  ) * real( nonhZ(i + k,j) , kind=rKind )
                  !
               endif
               !
            enddo

            fac = par%dt / s%hh(i,j) / s%dsz(i,j)

            Wm(i,j) = Wm_old(i,j) + fac * Wm_old(i,j)*( s%qx(i,j) - s%qx(i-1,j) ) - fac*( qx(1) * faceval(1) - qx(0) * faceval(0) )
            !
         else
            !
            Wm(i,j) = 0.0_rKind
            Wm_old(i,j) = 0.0_rKind
            s%pres(i,j) = 0.0_rKind
            !
         endif
         !
      enddo
      !
      !----------------- Advection dU (FOU) -------------------
      !
      if ( par%nhlay > 0. ) then
         !
         do i=2,s%nx-1
            !
            if ( s%wetu( i,j) == 1) then
               !
               do k = 0,1
                  !
                  qx(k) = .5 * ( s%qx(i - 1 + k , j ) + s%qx(i     + k , j ) )

                  if (qx(k) > 0.) then
                     !
                     faceval(k) = dU0( i - 1 + k , j )
                     !
                  else
                     !
                     faceval(k) = dU0( i+ k , j )
                     !
                  endif
                  !
               enddo

               fac = par%dt / s%hum(i,j) / s%dsu(i,j)
                       
               s%dU(i,j) = dU0(i,j) + .5 * fac *dU0(i,j) * ( s%qx(i+1,j) - s%qx(i-1,j) ) &
                    - fac * ( qx(1)*faceval(1) - qx(0)*faceval(0) )                  
               !
            else
               !
               s%dU(i,j) = 0.0
               !
            endif
            !
         enddo
         !
      endif
      !

      !
      if (xmpi_istop) then
         !
         au(:,1,:)      = 0.0_rKind
         adU(:,1,:)     = 0.0_rKind
         !
      endif
      !
      if (xmpi_isbot) then
         !
         au(:,s%nx,:)   = 0.0_rKind
         adU(:,s%nx,:)  = 0.0_rKind
         !
      endif
      !
      if (par%nhlay > 0.) then
         !
         do i=2,s%nx
            !
            if ( nonhZ(i,j) == 1 ) then
               !
               ! Vertical advection (First order upwind)
               !

               if ( omega(i,j) > 0.) then
                  !
                  s%wm(i,j) = s%wm(i,j) + par%dt * omega(i,j) * wm0(i,j)/s%hh(i,j) *  wcoef(i,j)
                  !
               else
                  !
                  s%wm(i,j) = s%wm(i,j) + par%dt * omega(i,j) * wm_old(i,j)/s%hh(i,j) * wcoef(i,j)
               endif
            endif

         enddo
      endif

      !Include explicit approximation for pressure in s%uu du and Wm
      if (par%secorder == 1) then
         !
         do i=2,s%nx-1
            !
            s%uu(i,j) = s%uu(i,j) + au(1,i,j) * s%pres(i+1,j) + au(0,i,j) * s%pres(i  ,j)
            !
         enddo
         !

         if ( par%nhlay > 0. ) then
            !
            do i=2,s%nx-1
               !
               s%dU(i,j) = s%dU(i,j) + adU(1,i,j) * s%pres(i+1,j) + adU(0,i,j) * s%pres(i  ,j)
               !
            enddo
            !
         endif

         !
         do i=2,s%nx
            !
            Wm(i,j) = Wm(i,j) + aws(1,i,j) * s%pres(i,j)
            !
         enddo
         !

#ifdef USEMPI
         !
         call xmpi_shift_ee(Wm)
#endif
         !
         call flow_secondorder_advW(s,par,Wm,Wm_old, nonhZ)
         !
#ifdef USEMPI
         !
         call xmpi_shift_ee(Wm)
         !
#endif
      endif
      !

      if (par%nhlay > 0.) then
         !
         !
         !
         !
         !----------------- Bottom stress, internal stresses and vertical advection -------------------
         !
         ! Integration of these terms occurs
         ! implicity to ensure that they do not introduce any additional stability
         ! concerns.
         !
         ! Moreover, the vertical exchange of momentum due to advection between the
         ! layers also influences dU (but not U), so we have
         !
         ! the equation reads dU - dU*
         !                    -------- =  dU ( fric + stress )  + s%uu * fric - omega * U_face
         !                       dt
         !
         do i=2,s%nx-1
            !
            if ( s%wetu(i,j) == 1 ) then
               !
               !
               ! Bottom friction
               !
               ub = uu0(i,j) + (1.-par%nhlay) * dU0(i,j)
               umag = abs(ub)
               fric   = par%dt * s%cfu(i,j) * umag / s%hu(i,j)

               !
               ! Internal stresses
               !
               ve     = 1.0e-4 !Vertical eddy viscocity
               stress = 2.* par%dt * ve / ( s%hum(i,j)**2 * (par%nhlay - par%nhlay**2) )

               !
               ! Vertical advection
               !
               advec = 0.5d0* ( omega(i,j) + omega(i+1,j) )* par%dt / s%hum(i,j)

               !
               ! Built rhs/diagonal
               !
               fac   = 0.5d0*( wcoef(i,j) + wcoef(i+1,j) )
               rhs  = s%dU(i,j) - uu0(i,j) * ( fric + advec )
               Amat = 1.0d0 +  (1.0d0-fac) * fric + stress &
               + (1.0d0 - fac)     * max( advec , 0. ) &
               - fac               * min( advec , 0. )

               !
               ! "Solve" system
               !
               s%dU(i,j) = rhs / Amat
               !
            endif
            !
         enddo
         !
      endif

#ifdef USEMPI
      !
      call xmpi_shift_ee(s%dU)
      call xmpi_shift_ee(s%dV)
      !
#endif
      !
   end subroutine nonh_2lay_pred_2dv


   !==============================================================================
   !                     REDUCED NON-HYDROSTATIC ROUTINES 3D
   !==============================================================================
   !


   !
   !==============================================================================
   subroutine nonh_2lay_cor_3d(s,par)
      !==============================================================================
      !

      ! DATE               AUTHOR               CHANGES
      !
      ! October 2010       Pieter Bart Smit     New Subroutine
      ! November  2010     Pieter Bart Smit     Added explicit prediction for 2nd order scheme

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------

      !  In this subroutine we calculate the nonh-hydrostatic pressures by solving
      !  a discrete poisson equation, and subsequently correct the velocity field
      !  so that the velocities at the new timelevel are divergence free.
      !
      !  This subroutine uses the pressure coeficients calculated previously in
      !  the routine nonh_2lay_pred_3d, and for each of the velocity components all
      !  processes have been accounted for (advection, bed-friction etc), including
      !  an explict estimation of the nonh press known from the previous timestep
      !  (unless we are calculating with first order accuracy).
      !
      !  Subsequently we will take the following steps:
      !
      !  1) Calculate the water and bed levels at u-points (subroutine call to zuzv).
      !  2) Built the poisson matrix by substituting the pressure dependencies of the
      !     momentum equations into the continuity equation.
      !  3) Solve the resulting system (call to solver_solvemat).
      !  4) Update the non-hydrostatic pressure. If second-order is active, dp is the
      !     difference between the old and the new pressure, so in that case we have
      !
      !     s%press = s%press + dp
      !
      !     else, it represents the new pressure directly, s%press = dp.
      !  5) Correct the velocity components with the nonhydrostatic pressures. The
      !     new velocity field will be divergence free.
      !
      !  After these steps, control is returned to flow_timestep.
      !
      !--------------------------     DEPENDENCIES       ----------------------------
      use spaceparams
      use params
      use solver_module
      use xmpi_module
      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars) ,intent(inout)                       :: s
      type(parameters),intent(in)                          :: par

      !--------------------------     LOCAL VARIABLES    ----------------------------

      !Indices
      integer(kind=iKind)                     :: i,j,k !Index variables

      !Work
      real(kind=rKind)                        :: dzsdx        ! Free surface gradient (s-dir)
      real(kind=rKind)                        :: dzbdx        ! Bed level gradient    (s-dir)
      real(kind=rKind)                        :: dzsdy        ! Free surface gradient (n-dir)
      real(kind=rKind)                        :: dzbdy        ! Bed level gradient    (n-dir)
      real(kind=rKind)                        :: alpha        ! Parameter controlling how the velocity at the
      ! layer interface is determined ( =-wcoef for energy
      ! conservative behaviour).
      real(kind=rKind)                        :: alpha_u(0:1) ! The layer distribution parameter interpolated at
      ! the east (=1) and west (=0) faces.
      real(kind=rKind)                        :: alpha_v(0:1) ! The layer distribution parameter interpolated at
      ! the north (=1) and south (=0) faces.
      real(kind=rKind)                        :: hfU(0:1)     ! Coefficient of the west U(i-1,j) (=0) or east U(i,j) (=1) velocity
      ! in the reduced continuity equation.
      real(kind=rKind)                        :: hfdU(0:1)    ! Coefficient of the west dU(i-1,j) (=0) or east dU(i,j) (=1) velocity
      ! difference in the reduced continuity equation.
      real(kind=rKind)                        :: hfV(0:1)     ! Coefficient of the south V(i,j-1) (=0) or North V(i,j) (=1) velocity
      ! in the reduced continuity equation.
      real(kind=rKind)                        :: hfdV(0:1)    ! Coefficient of the south dV(i,j-1) (=0) or North dV(i,j) (=1) velocity
      ! difference in the reduced continuity equation.

      !
      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------
      !

      !
      !Substitute in the continuity equation
      !
      do j=2,s%ny
         !
         do i=2,s%nx
            !
            if (nonhZ(i,j)==1) then
               !
               do k = 0,1
                  !
                  alpha_u(k) =  .5d0 * (wcoef(i+k,j) + wcoef(i-1+k,j))

                  hfU(k)  = s%dsdnzi(i,j) * s%dnu(i-1+k,j) * s%hu(i-1+k  ,j) * ( 1.0d0 + alpha_u(k) )
                  hfdU(k) = s%dsdnzi(i,j) * s%dnu(i-1+k,j) * s%hu(i-1+k  ,j) * alpha_u(k) * ( 1.0d0 - alpha_u(k) )
                  !
               enddo

               do k = 0,1
                  !
                  alpha_v(k) =  .5d0 * (wcoef(i,j+k) + wcoef(i,j-1+k))

                  hfV(k)  = s%dsdnzi(i,j) * s%dsv(i,j-1+k) * s%hv(i,j-1+k  ) * ( 1.0d0 + alpha_v(k) )
                  hfdV(k) = s%dsdnzi(i,j) * s%dsv(i,j-1+k) * s%hv(i,j-1+k  ) * alpha_v(k) * ( 1.0d0 - alpha_v(k) )
                  !
               enddo
               !

               if (nonhu(i,j)*nonhu(i-1,j) == 1  ) then
                  !
                  dzsdx = .5*( zsu(i,j) - zsu(i-1,j  ) ) / s%dsz(i,j)
                  dzbdx = .5*( zbu(i,j) + alpha_u(1) * s%hu(i,j) - zbu(i-1,j  ) - alpha_u(0) * s%hu(i-1,j  )  ) / s%dsz(i,j)
                  !
               else
                  dzsdx  = 0.0d0
                  dzbdx  = 0.0
               endif

               if (nonhv(i,j)*nonhv(i,j-1) == 1  ) then
                  !
                  dzsdy = .5*( zsv(i,j) - zsv(i  ,j-1) ) / s%dnz(i,j)
                  dzbdy = .5*( zbv(i,j) + alpha_v(1) * s%hv(i,j) - zbv(i  ,j-1) - alpha_v(0) * s%hv(i  ,j-1)  ) / s%dnz(i,j)
               else
                  !
                  dzsdy = 0.0d0
                  dzbdy = 0.0d0
                  !
               endif

               alpha = -wcoef(i,j)
               !
               ! Main diagonal ( i,j )
               !
               mat(1,i,j) = 2.0d0 * aws(1,i,j) &
               + hfU(1)   * au(0, i  ,j) - hfU(0)  * au(1, i-1,j  )  &
               + hfV(1)   * av(0, i  ,j) - hfV(0)  * av(1, i  ,j-1)  &
               + hfdU(1)  * adU(0,i  ,j) - hfdU(0) * adU(1,i-1 ,j )  &
               + hfdV(1)  * adV(0,i  ,j) - hfdV(0) * adV(1,i   ,j-1) &
               - dzsdx * ( au( 0,i,j) + au( 1,i-1,j  ) - ( adU(0,i,j) + adU(1,i-1,j  ) ) * wcoef(i,j)  ) &
               - dzsdy * ( av( 0,i,j) + av( 1,i  ,j-1) - ( adV(0,i,j) + adV(1,i  ,j-1) ) * wcoef(i,j)  ) &
               - dzbdx * ( au( 0,i,j) + au( 1,i-1,j  ) + ( adU(0,i,j) + adU(1,i-1,j  ) ) * alpha   ) &
               - dzbdy * ( av( 0,i,j) + av( 1,i  ,j-1) + ( adV(0,i,j) + adV(1,i  ,j-1) ) * alpha   )
               !
               ! West diagonal ( i-1,j )
               !
               mat(2,i,j) = - hfU(0)  * au (0, i-1,j)    &
               - hfdU(0) * adU(0, i-1,j)    &
               - dzsdx * (  au(0,i-1,j) - adU(0,i-1,j) * wcoef(i,j)  ) &
               - dzbdx * (  au(0,i-1,j) + adU(0,i-1,j) * alpha   )
               !
               ! South diagonal (i,j-1 )
               !
               mat(4,i,j) = - hfV(0)  * av (0, i,j-1)    &
               - hfdV(0) * adV(0, i,j-1)    &
               - dzsdy * (  av(0,i,j-1) - adV(0,i,j-1) * wcoef(i,j) ) &
               - dzbdy * (  av(0,i,j-1) + adV(0,i,j-1) * alpha  )
               !
               ! East diagonal ( i+1,j )
               !
               mat(3,i,j) = hfU(1)  * au (1,i  ,j)    &
               + hfdU(1) * adU(1,i  ,j)    &
               -  dzsdx * (  au(1,i,j) - adU(1,i,j) * wcoef(i,j)  ) &
               -  dzbdx * (  au(1,i,j) + adU(1,i,j) * alpha   )
               !
               ! North diagonal ( i,j+1 )
               !
               mat(5,i,j) = hfV(1)   * av (1,i  ,j)    &
               + hfdV(1)  * adV(1,i  ,j)    &
               - dzsdy * (  av(1,i,j) - adV(1,i,j) * wcoef(i,j)  ) &
               - dzbdy * (  av(1,i,j) + adV(1,i,j) * alpha   )
               !
               ! RHS
               !
               rhs(i,j)   = - 2.0d0  * Wm(i,j) &
               - hfU(  1 ) * s%uu( i  ,j)  + hfU(  0 ) * s%uu( i-1 ,j  )   &
               - hfdU( 1 ) * s%dU( i  ,j)  + hfdU( 0 ) * s%dU( i-1 ,j  )   &
               - hfV(  1 ) * s%vv( i  ,j)  + hfV(  0 ) * s%vv( i   ,j-1)   &
               - hfdV( 1 ) * s%dV( i  ,j)  + hfdV( 0 ) * s%dV( i   ,j-1)   &
               + dzsdx * ( s%uu(i,j) + s%uu(i-1,j  ) - ( s%dU(i,j) + s%dU(i-1,j  ) ) * wcoef(i,j)  ) &
               + dzsdy * ( s%vv(i,j) + s%vv(i  ,j-1) - ( s%dV(i,j) + s%dV(i  ,j-1) ) * wcoef(i,j)  ) &
               + dzbdx * ( s%uu(i,j) + s%uu(i-1,j  ) + ( s%dU(i,j) + s%dU(i-1,j  ) ) * alpha   ) &
               + dzbdy * ( s%vv(i,j) + s%vv(i  ,j-1) + ( s%dV(i,j) + s%dV(i  ,j-1) ) * alpha   )
               !
            else
               !
               mat(1  ,i,j) = 1.0_rKind
               mat(2:5,i,j) = 0.0_rKind
               rhs(i,j)     = 0.0_rKind
               !
            endif
         enddo
      enddo

      !------------------ Solve matrix -----------------
      dp = 0.0_rKind
      call solver_solvemat( mat  , rhs   , dp , s%nx, s%ny,par)

      if ( par%secorder == 1 ) then
         !
         ! If second-order is active, we are calculating with the pressure difference
         !
         s%pres = s%pres + dp
         !
      else
         !
         s%pres = dp
         !
      endif


      !Correct u/v/w

      !U
      do j=2,s%ny
         !
         do i=2,s%nx-1
            !
            s%uu(i,j) = s%uu(i,j) + au(1,i,j)*dp(i+1,j)+au(0,i,j)*dp(i,j)
            !
         enddo
         !
      enddo

      !V
      do j=2,s%ny-1
         !
         do i=2,s%nx
            !
            s%vv(i,j) = s%vv(i,j) + av(1,i,j)*dp(i,j+1)+av(0,i,j)*dp(i,j)
            !
         enddo
         !
      enddo
      !

      !
      if (par%nhlay > 0.) then
         !
         !dU
         !
         do j=2,s%ny
            !
            do i=2,s%nx-1
               !
               s%dU(i,j) = s%dU(i,j) + adU(1,i,j)*dp(i+1,j)+adU(0,i,j)*dp(i,j)
               !
            enddo
            !
         enddo
         !
         !dV
         !
         do j=2,s%ny-1
            !
            do i=2,s%nx
               !
               s%dV(i,j) = s%dV(i,j) + adV(1,i,j)*dp(i,j+1)+adV(0,i,j)*dp(i,j)
               !
            enddo
            !
         enddo
         !
      endif
      !

      !
      ! Communicate results...
      !
#ifdef USEMPI
      call xmpi_shift_ee(s%uu)
      call xmpi_shift_ee(s%vv)
      call xmpi_shift_ee(s%dU)
      call xmpi_shift_ee(s%dV)
#endif

      !
      ! Calculate the surface vertical velocity from the kinematic condition
      !
      do j=2,s%ny
         !
         do i=2,s%nx
            !
            if ( s%wetz(i,j)==1 ) then
               !
               dzsdx = .5*( zsu(i,j) - zsu(i-1,j  ) )  / s%dsz(i,j)
               dzsdy = .5*( zsv(i,j) - zsv(i  ,j-1) )  / s%dnz(i,j)

               s%ws(i,j) = - (s%dnu(i,j)*s%hu(i,j)*s%uu(i,j) -s%dnu(i-1,j)*s%hu(i-1,j  )*s%uu(i-1,j  )) * s%dsdnzi(i,j) &
               - (s%dsv(i,j)*s%hv(i,j)*s%vv(i,j) -s%dnv(i,j-1)*s%hv(  i,j-1)*s%vv(i  ,j-1)) * s%dsdnzi(i,j) &
               +  dzsdx * (  s%uu(i,j) + s%uu(i-1,j  ) - ( s%dU(i,j) + s%dU(i-1,j  ) ) * wcoef(i,j) ) &
               +  dzsdy * (  s%vv(i,j) + s%vv(i  ,j-1) - ( s%dV(i,j) + s%dV(i  ,j-1) ) * wcoef(i,j) )
               !
            else
               !
               s%ws(i,j) = 0.0d0
               !
            endif
            !
         enddo
         !
      enddo

      !
      ! Calculate the bottom velocity from the kinematic condition
      !
      do j=2,s%ny
         !
         do i=2,s%nx
            !
            if ( s%wetz(i,j)==1 ) then
               !
               dzbdx = .5*( zbu(i,j)  - zbu(i-1,j  ) )
               dzbdy = .5*( zbv(i,j)  - zbv(i  ,j-1) )

               s%wb(i,j) = + dzbdx * ( s%uu(i,j) + s%uu(i-1,j  ) + ( s%dU(i,j) + s%dU(i-1,j  ) ) * (1- wcoef(i,j)) ) / s%dsz(i,j) &
               +  dzbdy * (  s%vv(i,j) + s%vv(i  ,j-1) + ( s%dV(i,j) + s%dV(i  ,j-1) ) * (1- wcoef(i,j)) ) / s%dnz(i,j)
               !
            else
               !
               s%wb(i,j) = 0.0d0
               !
            endif
            !
         enddo
         !
      enddo

      !
      ! Calculate the mean vertical velocity in the top layer
      !
      do j=2,s%ny
         !
         do i=2,s%nx
            !
            if ( nonhZ(i,j) == 1) then
               !
               Wm(i,j) = Wm(i,j) + aws(1,i,j)*dp(i,j)
               !
            else
               !
               Wm(i,j) = 0.0d0
               !
            endif
            !
         enddo
         !
      enddo

      !
      ! Calculate the relative vertical advection
      !
      do j=2,s%ny
         !
         do i=2,s%nx
            !
            if (s%wetz(i,j)==1) then
               !
               omega(i,j) = - ( s%dnu(i,j) * s%hu(i  ,j) * dU0(i  ,j) - s%dnu(i-1,j) * s%hu(i-1,j) * dU0(i-1,j) &
               +   s%dsv(i,j) * s%hv(i  ,j) * dV0(i  ,j) - s%dsv(i,j-1) * s%hv(i,j-1) * dV0(i,j-1)  ) * s%dsdnzi( i , j )
               !
            else
               !
               omega(i,j) = 0.0d0
               !
            endif
            !
         enddo
         !
      enddo

      !
      ! Calculate the mean vertical velocity in the lower layer
      !
      do j=2,s%ny
         !
         do i=2,s%nx
            !
            !
            if ( s%wetz(i,j)==1 ) then
               !
               Wm0(i,j) =  Wm(i,j) + .5d0 * ( s%wb(i,j) - s%ws(i,j) )
               !
            else
               !
               Wm0(i,j) = 0.
               !
            endif
            !
         enddo
         !
      enddo


#ifdef USEMPI
      if (xmpi_istop) then
         Wm(1,:) = Wm(2,:)
      endif
      if (xmpi_isbot) then
         Wm(s%nx+1,:) = Wm(s%nx,:)
      endif
      if (xmpi_isleft) then
         Wm(1,:) = Wm(2,:)
      endif
      if (xmpi_isright) then
         Wm(:,s%ny+1) = Wm(:,s%ny)
      endif
      !
      call xmpi_shift_ee(Wm)
      !
#else
      !Assign boundaries
      Wm(1,:) = Wm(2,:)
      Wm(s%nx+1,:) = Wm(s%nx,:)
      Wm(:,1) = Wm(:,2)
      Wm(:,s%ny+1) = Wm(:,s%ny)

#endif
      Wm_old = Wm
      dU0  = s%dU
      dV0  = s%dV
      !
   end subroutine nonh_2lay_cor_3d


   subroutine nonh_2lay_pred_3d(s,par,uu0,vv0)
      !==============================================================================
      !

      ! DATE                AUTHOR               CHANGES
      !
      ! October  2014       Pieter Bart Smit     New Subroutine


      !------------------------------------------------------------------------------
      !                             DECLARATIONS
      !------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------
      !
      ! Include the pressure explicitly in the predictor step. In this subroutine
      ! the following steps are taken:
      !
      ! 1) Calculate whether or not breaking is active (handled in subroutine
      !    nonh_masks ).
      ! 2) Calculate the coeficients of the nonh. pres. as they occur in the momentum
      !    equations. (these are use to calculate the explicit pressure contributions
      !    here, and the implicit contributions in the correction routine later).
      ! 3) Calculate the first order additional contributions to the advection due to
      !    the reduced two-layer formulation (if applicable).
      ! 4) Include an explicit estimate for the nonh press. in the momentum equations
      !    using the earlier calculated coef. (if second-order approximations for the
      !    advection are used, else this step is not necessary).
      ! 5) Include the contributions due to internal stress terms and bottom friction
      !    in the equation for the velocity difference (implicitly for stability).
      !
      ! After these steps, control flows back to flow_timestep, and - if active -
      ! second-order approximations are calculated. Consequently, flow_timestep calls
      ! nonh_cor again and the nonhydrostatic computations are continued in
      ! nonh_2lay_cor.
      !
      !--------------------------     DEPENDENCIES       ----------------------------
      use spaceparams
      use params
      use flow_secondorder_module
      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars) ,intent(inout)                       :: s
      type(parameters),intent(in)                          :: par
      real(kind=rkind),dimension(s%nx+1,s%ny+1),intent(in) :: uu0
      real(kind=rkind),dimension(s%nx+1,s%ny+1),intent(in) :: vv0

      !--------------------------     VARIABLES          ----------------------------

      !General
      real(kind=rKind)                :: Fac          !A factor

      !Indices
      integer(kind=iKind)             :: i,j,k          !Loop indices

      !Used in calculation of pressure gradients
      real(kind=rKind)                :: Vol          !The volume of a momentum cell

      !Used in calculation of advection dU/dV/Wm
      real(kind=rKind)                :: qx(0:1)        !Discharge through the faces of a momentum cell (0: left face
      !1: right face)
      real(kind=rKind)                :: qy(0:1)        !Discharge through the faces of a momentum cell (0: bottom face
      !1: top face)
      real(kind=rKind)                :: faceval_x(0:1) !Value of dU at the faces of a momentum cell (0: left face
      !1: right face)
      real(kind=rKind)                :: faceval_y(0:1) !Value of dU at the faces of a momentum cell (0: left face
      !1: right face)
      real(kind=rKind)                :: DU_face
      real(kind=rKind)                :: DV_face
      real(kind=rKind) :: cosd,sind                     !cosine or sine of the difference angle
      real(kind=rKind) :: sig(0:1)                      !sign function

      !Used in calculation of vertical exchange of momentum
      real(kind=rKind)                :: fric         !Bottom friction coeficient (= dt * cf * |U| / h)
      real(kind=rKind)                :: stress       !Coeficient used in calculation of stress
      real(kind=rKind)                :: ub           !The (Eulerian) u-velocity at the bottom
      real(kind=rKind)                :: vb           !The (Eulerian) v-velocity at the bottom
      real(kind=rKind)                :: umag         !The (Eulerian) magnitude of the velocity vector at the bottom
      real(kind=rKind)                :: ve           !The eddy viscosity used for vertical mixing
      real(kind=rKind)                :: advec        !Coeficient used in calculation of vertical advection
      real(kind=rKind)                :: Amat,rhs     !The 1*1 "matrix" and RHS of the vertical "system".

      !------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !------------------------------------------------------------------------------

      !
      ! Determine whether or not to use the nh correction method (dry/wet/breaking etc.)
      !
      call nonh_break( s , par )
      call zuzv(s)

      !
      ! Calculate the pressure dependencies for U
      !
      do j=2,s%ny
         !
         do i=1,s%nx
            !
            if (nonhU(i,j)==1) then
               !
               vol       = 0.5d0 * par%dt/(s%hum(i,j)*s%dsu(i,j))
               Fac       = (s%zb(i+1,j) - s%zb(i  ,j)) * vol
               au(1,i,j) = - vol * (1.d0 + wcoef(i,j) ) * s%hh(i+1,j) - Fac
               au(0,i,j) = + vol * (1.d0 + wcoef(i,j) ) * s%hh(i,j)   - Fac
               !
            else
               !
               au(1,i,j) =  0.0_rKind
               au(0,i,j) =  0.0_rKind
               !
            endif
            !
         enddo
         !
      enddo
      !
      ! Calculate the pressure dependencies for V
      !
      do j=1,s%ny
         !
         do i=2,s%nx
            !
            if (nonhV(i,j)==1) then
               !
               vol       = 0.5d0 * par%dt/(s%hvm(i,j)*s%dnv(i,j))
               Fac       = (s%zb(i,j+1) - s%zb(i  ,j)) * vol
               av(1,i,j) = - vol * (1.d0 + wcoef(i,j+1) ) * s%hh(i,j+1) - Fac
               av(0,i,j) = + vol * (1.d0 + wcoef(i,j  ) ) * s%hh(i,j)   - Fac
               !
            else
               !
               av(1,i,j) =  0.0_rKind
               av(0,i,j) =  0.0_rKind
               !
            endif
            !
         enddo
         !
      enddo

      if ( par%nhlay > 0.0d0 ) then
         !
         ! Calculate the pressure dependencies for dU
         !
         do j=2,s%ny
            !
            do i=1,s%nx
               !
               if (nonhU(i,j)*(1-s%breaking(i,j))*(1-s%breaking(max(i-1,1),j))==1 ) then
                  !
                  vol       = 0.5_rKind*par%dt/(s%hum(i,j)*s%dsu(i,j))

                  fac = (  s%zs( i+1,j ) -  s%zs(i,j) ) * vol / ( 1.0d0 - .5d0*wcoef(i,j) - .5d0*wcoef(i+1,j) ) &
                  - ( wcoef( i+1,j ) - wcoef(i,j) ) * s%hum(i,j) * vol / (2.0d0 - wcoef(i,j) - wcoef(i+1,j) )

                  adU(1,i,j) = - s%hh(i+1,j)*vol + Fac
                  adU(0,i,j) = + s%hh(i,j)*vol   + Fac
                  !
               else
                  !
                  adU(1,i,j) =  0.0_rKind
                  adU(0,i,j) =  0.0_rKind
                  !
               endif
               !
            enddo
            !
         enddo
         !
         ! Calculate the pressure dependencies for dV
         !
         do j=1,s%ny
            !
            do i=2,s%nx
               !
               if (nonhV(i,j)*(1-s%breaking(i,j))*(1-s%breaking(i,max(j-1,1)))==1 ) then
                  !
                  vol       = 0.5_rKind*par%dt/(s%hvm(i,j)*s%dnv(i,j))
                  fac = (  s%zs( i,j ) -  s%zs(i,j+1) ) * vol / ( 1.0d0 - .5d0*wcoef(i,j) - .5d0*wcoef(i,j+1) ) &
                  - ( wcoef( i,j ) - wcoef(i,j+1) ) * s%hvm(i,j) * vol / (2.0d0 - wcoef(i,j) - wcoef(i,j+1) )
                  adV(1,i,j) = - s%hh(i,j+1)*vol + Fac
                  adV(0,i,j) = + s%hh(i,j)*vol   + Fac
                  !
               else
                  !
                  adV(1,i,j) =  0.0_rKind
                  adV(0,i,j) =  0.0_rKind
                  !
               endif
               !
            enddo
            !
         enddo
         !
      endif

      !
      !Built pressure coefficients W
      !
      do j=2,s%ny
         !
         do i=2,s%nx
            !
            if (nonhZ(i,j) == 1) then
               !
               aws(1,i,j)   = par%dt/( s%hh(i,j) * ( 1 - wcoef(i,j) ) )
               !
            else
               !
               aws(1,i,j)   =  0.0_rKind
               !
            endif
            !
         enddo
         !
      enddo
      !

      !
      if (xmpi_istop) then
         !
         au(:,1,:)       = 0.0_rKind
         adU(:,1,:)      = 0.0_rKind
         !
      endif
      !
      if (xmpi_isbot) then
         !
         au(:,s%nx,:)    = 0.0_rKind
         adU(:,s%nx,:)   = 0.0_rKind
         !
      endif
      !

      !
      if (xmpi_isleft) then
         !
         av(:,:,1)      = 0.0_rKind
         adV(:,:,1)     = 0.0_rKind
         !
      endif
      !

      !
      if (xmpi_isright) then
         !
         av(:,:,s%ny)   = 0.0_rKind
         adV(:,:,s%ny)  = 0.0_rKind
         !
      endif

      !
      !----------------- Advection Wm (FOU) -------------------
      !
      !  d ( h U Wm )    d ( h V Wm )
      !  ------------ +  ------------
      !      ds              dn
      !
      ! Advection according to Stelling & Duinmeijer 2003, momentum conservative
      !
      do j=2,s%ny
         !
         do i=2,s%nx
            !
            if ( nonhZ(i,j)==1 ) then
               !
               do k = 0,1
                  !
                  qx(k) = s%qx(i - 1 + k , j ) * s%dnu(i - 1 + k , j )
                  qy(k) = s%qy(i , j - 1 + k)  * s%dsv(i , j - 1 + k )

                  if (qx(k) > 0.) then
                     !
                     faceval_x(k) = Wm_old( i - 1 + k ,j  ) * real( nonhZ(i - 1 + k,j) , kind=rKind )
                     !
                  else
                     !
                     faceval_x(k) = Wm_old( i + k ,j  ) * real( nonhZ(i + k,j) , kind=rKind )
                     !
                  endif

                  if (qy(k) > 0.) then
                     !
                     faceval_y(k) = Wm_old( i ,j - 1 + k ) * real( nonhZ(i,j - 1 + k) , kind=rKind )
                     !
                  else
                     !
                     faceval_y(k) = Wm_old( i,j + k ) * real( nonhZ(i,j + k) , kind=rKind )
                     !
                  endif
                  !
               enddo

               fac = par%dt / s%hh(i,j) * s%dsdnzi(i,j)

               Wm(i,j) = Wm_old(i,j) +  fac * Wm_old(i,j) * ( qx(1) + qy(1) - qx(0) - qy(0) ) &
                    -  fac * ( qx(1) * faceval_x(1) +  qy(1) * faceval_y(1) &
                    - qx(0) * faceval_x(0) - qy(0) * faceval_y(0) )
               !
            else
               !
               Wm(i,j) = 0.0_rKind
               Wm_old(i,j) = 0.0_rKind
               s%pres(i,j) = 0.0_rKind
               !
            endif
            !
         enddo
         !
      enddo

      if ( par%nhlay > 0.0d0 ) then
         !
         !
         !----------------- Advection dU (FOU) -------------------
         !
         !  d ( h U dU )    d ( h V dU )
         !  ----------   +  ----------
         !      dx              dy
         !
         ! Advection according to Stelling & Duinmeijer 2003, momentum conservative (second term in their equation 24)
         !
         sig(0) = -1
         sig(1) = 1
         !
         do j=2,s%ny
            !
            do i=2,s%nx-1
               !
               if ( s%wetu( i,j) == 1) then
                  !
                  ! Calculate transport velocities / face values
                  !
                  do k = 0,1
                     !
                     qx(k) = .5 * ( s%qx(i - 1 + k , j   ) + s%qx(i     + k , j ) ) * s%dnz(i+k,j)
                     qy(k) = .5 * ( s%qy(i+1, j - 1 + k  ) + s%qy(i, j  - 1 + k ) ) * s%dsc(i,j+k)
                     !

                     dv_face = 0.5d0* ( dV0( i + k , j - 1 ) + dV0( i + k , j ) )

                     cosd = cos( s%alfau(i+k,j) - s%alfau(i+k-1,j) )
                     sind = sig(k)*sin( s%alfau(i+k,j) - s%alfaz(i+k-1,j) )
                     !
                     if (qx(k) > 0.) then
                        !
                        faceval_x(k) = dU0( i - 1 + k , j ) * cosd + dV_face * sind
                        !
                     else
                        !
                        faceval_x(k) = dU0( i+ k , j ) * cosd + dV_face * sind
                        !
                     endif
                     !

                     !
                     cosd = cos( s%alfau(i,j+k-1) - s%alfau(i,j+k) )
                     sind = sig(k)*sin( s%alfau(i,j+k-1) - s%alfau(i,j+k) )

                     dV_face = 0.5d0* ( dV0( i  , j + k - 1) + dV0( i + 1 , j + k - 1 ) )
                     if (qy(k) > 0.) then
                        !
                        faceval_y(k) = dU0( i , j -1 + k ) * cosd + dV_face*sind
                        !
                     else
                        !
                        faceval_y(k) = dU0( i , j + k) * cosd + dV_face*sind
                        !
                     endif
                     !
                  enddo

                  fac = par%dt / s%hum(i,j) * s%dsdnui(i,j)

                  s%dU(i,j) = dU0(i,j) + fac *dU0(i,j) * ( qx(1) + qy(1) - qx(0) - qy(0) ) &
                  - fac *( qx(1)*faceval_x(1) + qy(1)*faceval_y(1)    &
                  - qx(0)*faceval_x(0) - qy(0)*faceval_y(0) )
                  !
               else
                  !
                  s%dU(i,j) = 0.0
                  !
               endif
               !
            enddo
            !
         enddo

         !
         !----------------- Advection dV (FOU) -------------------
         !
         !  d ( h V dV )   d ( h U dV )
         !  ----------   +  ----------
         !      dy             dx
         !
         !
         ! Advection according to Stelling & Duinmeijer 2003, momentum conservative (second term in their equation 25)
         !
         do j=2,s%ny-1
            !
            do i=2,s%nx
               !
               if ( s%wetv( i,j) == 1) then
                  !
                  ! Calculate transport velocities / face values
                  !
                  do k = 0,1
                     !
                     qy(k) = .5 * ( s%qy(i , j- 1 + k    ) + s%qy(i , j     + k  ) ) * s%dsz(i,j+k)
                     qx(k) = .5 * ( s%qx(i- 1 + k, j +1  ) + s%qx(i  - 1 + k , j) )  * s%dnc(i+k,j)
                     !

                     cosd = cos( s%alfav(i+k-1,j) - s%alfav(i+k,j) )
                     sind = sig(k)*sin( s%alfav(i+k-1,j) - s%alfav(i+k,j) )

                     dU_face = 0.5d0* ( dU0( i - 1 + k , j ) + dU0( i - 1 + k , j +1 ) )
                     if (qx(k) > 0.) then
                        !
                        faceval_x(k) = dV0( i- 1 + k, j   ) * cosd + dU_face * sind
                        !
                     else
                        !
                        faceval_x(k) = dV0( i + k, j  ) *cosd + dU_face * sind
                        !
                     endif
                     !

                     cosd = cos( s%alfav(i,j+k) - s%alfav(i,j+k-1) )
                     sind = sig(k)*sin( s%alfav(i,j+k) - s%alfav(i,j+k-1) )
                     !
                     dU_face = 0.5d0* ( dU0( i - 1 , j + k ) + dU0( i  , j + k ) )
                     if (qy(k) > 0.) then
                        !
                        faceval_y(k) = dV0( i , j -1 + k ) * cosd+ dU_face * sind
                        !
                     else
                        !
                        faceval_y(k) = dV0( i , j + k) * cosd + dU_face * sind
                        !
                     endif
                     !
                  enddo

                  fac = par%dt / s%hvm(i,j) * s%dsdnvi(i,j)

                  s%dV(i,j) = dV0(i,j) + fac *dV0(i,j) * ( qx(1) + qy(1) - qx(0) - qy(0) ) &
                  - fac *( qx(1)*faceval_x(1) + qy(1)*faceval_y(1)    &
                  - qx(0)*faceval_x(0) - qy(0)*faceval_y(0) )
                  !
               else
                  !
                  s%dV(i,j) = 0.0
                  !
               endif
               !
            enddo
            !
         enddo
         !
      endif

      !
      ! ------------- Include explicit approximations --------------------
      !
      if (par%secorder == 1) then
         !
         ! Use explicit estimate for the nonhydrostatic pressure in U-momentum
         !
         do j=2,s%ny
            !
            do i=2,s%nx-1
               !
               s%uu(i,j) = s%uu(i,j) + au(1,i,j) * s%pres(i+1,j) + au(0,i,j) * s%pres(i  ,j)
               !
            enddo
            !
         enddo
         !
         ! Use explicit estimate for the nonhydrostatic pressure in V-momentum
         !
         do j=2,s%ny-1
            !
            do i=2,s%nx
               !
               s%vv(i,j) = s%vv(i,j) + av(1,i,j) * s%pres(i,j+1) + av(0,i,j) * s%pres(i  ,j)
               !
            enddo
            !
         enddo

         !
         ! Use explicit estimate for the nonhydrostatic pressure in dU-momentum
         !
         if ( par%nhlay > 0.0d0 ) then
            !
            do j=2,s%ny
               !
               do i=2,s%nx-1
                  !
                  s%dU(i,j) = s%dU(i,j) + adU(1,i,j) * s%pres(i+1,j) + adU(0,i,j) * s%pres(i  ,j)
                  !
               enddo
               !
            enddo
            !
            ! Use explicit estimate for the nonhydrostatic pressure in dV-momentum
            !
            do j=2,s%ny-1
               !
               do i=2,s%nx
                  !
                  s%dV(i,j) = s%dV(i,j) + adV(1,i,j) * s%pres(i,j+1) + adV(0,i,j) * s%pres(i  ,j)
                  !
               enddo
               !
            enddo
            !
         endif
         !
         !
         ! Use explicit estimate for the nonhydrostatic pressure in Wm-momentum
         !
         do j=2,s%ny
            !
            do i=2,s%nx
               !
               Wm(i,j) = Wm(i,j)   + aws(1,i,j) *  s%pres(i,j)
               !
            enddo
            !
         enddo
         !

#ifdef USEMPI
         !
         call xmpi_shift_ee(Wm)
#endif
         !
         call flow_secondorder_advW(s,par,Wm,Wm_old, nonhZ)

#ifdef USEMPI
         call xmpi_shift_ee(Wm)
#endif
         !
      endif

      !
      if (par%nhlay > 0.) then
         !
         !
         !
         !
         !----------------- Bottom stress, internal stresses and vertical advection -------------------
         !
         ! Integration of bottom stress and internal stress terms occurs implicity to ensure that they
         ! do not introduce any additional stability concerns.
         !
         ! Moreover, the vertical exchange of momentum due to advection between the
         ! layers also influences dU (but not U), so we have
         !
         ! the equation reads dU - dU*
         !                    -------- =  dU ( fric + stress )  + s%uu * fric - omega * U_face
         !                       dt
         ! To ensure that the system is solvable, we approximate U_face by it's upwind value.
         !
         do j=2,s%ny-1
            !
            do i=2,s%nx-1
               !
               if ( s%wetu(i,j) == 1 .and. wcoef(i,j) > 0.0d0 ) then
                  !
                  !
                  ! Bottom friction
                  !
                  ub = uu0(i,j) + (1.- wcoef(i,j)) * s%dU(i,j)
                  vb = s%vu(i,j) + .25 * (1.- wcoef(i,j)) * (  dV0(i,j) + dV0(i+1,j) + dV0(i,j-1) + dV0(i+1,j-1) )
                  umag = sqrt( ub**2 + vb**2 )
                  fric   = par%dt * s%cfu(i,j) * umag / s%hu(i,j)

                  !
                  ! Internal stresses
                  !
                  ve     = 0.0001 !Vertical eddy viscocity
                  stress = 2.* par%dt * ve / ( s%hum(i,j)**2 * (par%nhlay - par%nhlay**2) )

                  !
                  ! Vertical advection
                  !
                  advec = .5d0 * ( omega(i,j) + omega(i+1,j) ) * par%dt / s%hum(i,j)

                  !
                  ! Built rhs/diagonal
                  !
                  fac   = 0.5d0*( wcoef(i,j) + wcoef(i+1,j) )
                  rhs  = s%dU(i,j) - uu0(i,j) * ( fric + advec )

                  Amat = 1.0d0 +  (1.0d0-fac) * fric + stress &
                  + (1.0d0 - fac) * max( advec , 0. ) &
                  - fac           * min( advec , 0. )

                  !
                  ! "Solve" system
                  !
                  s%dU(i,j) = rhs / Amat
                  !
               endif
               !
            enddo
            !
         enddo
         !

         !
         do j=2,s%ny-1
            !
            do i=2,s%nx
               !
               ! Velocity in the lower layer at u-point
               !
               if (s%wetv(i,j) == 1 .and. wcoef(i,j) > 0.0d0 ) then
                  !
                  ! Bottom friction
                  !
                  vb = vv0(i,j) + (1.- wcoef(i,j)) * s%dV(i,j)
                  ub = s%uv(i,j) + .25 * (1.- wcoef(i,j)) * (  dU0(i,j) + dU0(i,j+1) + dU0(i-1,j) + dU0(i-1,j+1) )
                  umag = sqrt( ub**2 + vb**2 )
                  fric   = par%dt * s%cfu(i,j) * umag / s%hv(i,j)

                  !
                  ! Internal stresses
                  !
                  ve     = 0.0001 !Vertical eddy viscocity
                  stress = 2.* par%dt * ve / ( s%hvm(i,j)**2 * (par%nhlay - par%nhlay**2) )

                  !
                  ! Vertical advection
                  !
                  advec = .5d0 * ( omega(i,j) + omega(i,j+1) ) * par%dt / s%hvm(i,j)

                  !
                  ! Built rhs/diagonal
                  !
                  fac   = 0.5d0*( wcoef(i,j) + wcoef(i+1,j) )
                  rhs  = s%dV(i,j) - vv0(i,j) * ( fric + advec )

                  Amat = 1.0d0 +  (1.0d0-fac) * fric + stress &
                  + (1.0d0 - fac) * max( advec , 0. ) &
                  - fac           * min( advec , 0. )

                  !
                  ! "Solve" system
                  !
                  s%dV(i,j) = rhs / Amat
                  !
               endif
               !
            enddo
            !
         enddo
         !
      endif

#ifdef USEMPI
      !ARRAYS DU and DV need to be communicated
      call xmpi_shift_ee(s%dU)
      call xmpi_shift_ee(s%dV)
#endif
      !
   end subroutine nonh_2lay_pred_3d

   subroutine nonh_break( s , par )
      !--------------------------        PURPOSE         ----------------------------
      !
      ! Include the pressure explicitly in the predictor step
      !
      !--------------------------     DEPENDENCIES       ----------------------------
      use spaceparams
      use params
      use flow_secondorder_module
      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars) ,intent(inout)                       :: s
      type(parameters),intent(in)                          :: par

      !--------------------------     LOCAL VARIABLES    ----------------------------

      !Indices
      integer(kind=iKind)                     :: i,ie
      integer(kind=iKind)                     :: j,js
      integer(kind=iKind)                     :: jmin,jmax,ineighbour
      real(kind=rKind)                        :: wmax,dzdt


      if (s%ny>0) then
         jmin = 2
         jmax = s%ny
      else
         jmin = 1
         jmax = 1
      endif

      !
      ! Determine if a velocity point will be included in the nonh pressure matrix, do not include if:
      !
      ! (1)  The point is dry
      ! (2)  The relative wave length kd of the smallest possible wave (L=2dx) is smaller than kdmin
      ! (3)  The interpolated waterlevel zs is below the bottom (steep cliffs with overwash situations)
      ! (4)  Where Miche breaker criterium applies -> bores are hydrostatic
      !            max steepness = H/L = maxbrsteep
      !            dz/dx = maxbrsteep
      !            dz/dx = dz/dt/c = w/c = w/sqrt(gh)
      !            wmax = maxbrsteep*sqrt(gh)

      !
      do j=1,s%ny+1
         !
         do i=1,s%nx+1
            !
            ie = min(s%nx,i+1)
            !
            if (  (s%wetU(i,j)==1 ) .and. ( 0.5_rKind*( s%zs(i,j) + s%zs(ie,j) ) > zbu(i,j) ) ) then
               !
               nonhU(i,j) = 1
               !
            else
               !
               nonhU(i,j) = 0
               !
            endif
            !
         enddo
         !
      enddo

      !
      if (s%ny>2) then
         !
         do j=1,s%ny+1
            !
            js = min(s%ny,j+1)
            !
            do i=1,s%nx+1
               !
               if (  (s%wetV(i,j)==1 ) .and. (0.5_rKind*(s%zs(i,j) + s%zs(i,js)) > zbv(i,j) ) ) then
                  !
                  nonhV(i,j) = 1
                  !
               else
                  !
                  nonhV(i,j) = 0
                  !
               endif
               !
            enddo
            !
         enddo
         !
      else
         !
         nonhV = 0
         !
      endif
      !
      ! Determine if a velocity point will be included in the nonh pressure matrix, include if
      ! any of the surrounding velocity points is included.
      !

      !
      do j=jmin,jmax
         !
         js = max( j-1 , 1 )
         !
         do i=2,s%nx
            !
            if (max( nonhV(i,j) , nonhV(i,js) ,nonhU(i,j), nonhU(i-1,j)) > 0) then
               !
               nonhZ(i,j) = 1
               !
            else
               !
               nonhZ(i,j) = 0
               !
            endif
            !
         enddo
         !
      enddo
      !

      if (par%nhbreaker == 1) then
         !
         ! Breaking is active
         !
         do j = jmin , jmax
            !
            js = max( j-1 , 1 )
            !
            do i = 2 , s%nx
               !
               dzdt = - ( s%uu(i,j) * s%hu(i,j) - s%uu(i-1,j) * s%hu(i-1,j) ) / s%dsz(i,j)
               !
               ineighbour = 0
               if ( s%breaking(i-1,j)==1 .or. s%breaking(i+1,j) == 1 ) then
                  !
                  ineighbour = 1
                  !
               endif

               if ( s%ny > 1 ) then
                  !
                  dzdt = dzdt - ( s%vv(i,j) * s%hv(i,j) - s%vv(i,j-1) * s%hv(i,j-1) ) / s%dnz(i,j)

                  if ( s%breaking(i,j-1) == 1 .or. s%breaking(i,j+1) == 1 ) then
                     !
                     ineighbour = 1
                     !
                  endif
                  !
               endif
               !



               wmax =  sqrt(par%g*s%hh(i,j))

               if ( s%breaking(i,j) == 0 ) then
                  !
                  if ( dzdt > par%maxbrsteep*wmax ) then
                     !
                     s%breaking(i,j) = 2
                     !
                  elseif (dzdt > par%reformsteep*wmax .and. ineighbour == 1) then
                     !
                     s%breaking(i,j) = 2
                     !
                  endif
                  !
               elseif ( s%breaking(i,j) == 1 ) then
                  !
                  if ( dzdt < 0.0d0 ) then
                     !
                     s%breaking(i,j) = -1
                     !
                  endif
                  !
               endif
               !
            enddo
            !
         enddo

         where (s%breaking==2)
            s%breaking = 1
         endwhere

         where (s%breaking==-1)
            s%breaking = 0
         endwhere

         ! turn off non-hydrostatic pressure correction in areas with breaking and increase viscosity
         where (s%breaking/=0)
            nonhZ = 0
            s%pres = 0
         endwhere
         !
      endif
      !
   end subroutine nonh_break






   subroutine nonh_masks( s , par )
      !--------------------------        PURPOSE         ----------------------------
      !
      ! Include the pressure explicitly in the predictor step
      !
      !--------------------------     DEPENDENCIES       ----------------------------
      use spaceparams
      use params
      use flow_secondorder_module
      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars) ,intent(inout)                       :: s
      type(parameters),intent(in)                          :: par

      !--------------------------     LOCAL VARIABLES    ----------------------------

      !Indices
      integer(kind=iKind)                     :: i,ie,iw               !Index variables
      integer(kind=iKind)                     :: j,js
      integer(kind=iKind)                     :: jmin,jmax
      real(kind=rKind)                        :: wmax,reformfac


      if (s%ny>0) then
         jmin = 2
         jmax = s%ny
      else
         jmin = 1
         jmax = 1
      endif

      !
      ! Determine if a velocity point will be included in the nonh pressure matrix, do not include if:
      !
      ! (1)  The point is dry
      ! (2)  The relative wave length kd of the smallest possible wave (L=2dx) is smaller than kdmin
      ! (3)  The interpolated waterlevel zs is below the bottom (steep cliffs with overwash situations)
      ! (4)  Where Miche breaker criterium applies -> bores are hydrostatic
      !            max steepness = H/L = maxbrsteep
      !            dz/dx = maxbrsteep
      !            dz/dx = dz/dt/c = w/c = w/sqrt(gh)
      !            wmax = maxbrsteep*sqrt(gh)

      do j=1,s%ny+1
         do i=1,s%nx+1
            iw = max(i,i-1)
            ie = min(s%nx,i+1)

            if (  (s%wetU(i,j)==1                                      )  &
            .and. (0.5_rKind*(s%zs(i,j) + s%zs(ie,j))    > zbu(i,j)    )  &
            .and. ( s%dsu(i,1)*par%kdmin/par%px  < s%hum(i,j)  )  ) then
               nonhU(i,j) = 1
            else
               nonhU(i,j) = 0
            endif
         enddo
      enddo

      if (s%ny>2) then
         do j=1,s%ny+1
            js = min(s%ny,j+1)
            do i=1,s%nx+1
               if (  (s%wetV(i,j)==1                                      )  &
               .and. (0.5_rKind*(s%zs(i,j) + s%zs(i,js))    > zbv(i,j)    )  &
               .and. ( s%dnv(1,j)*par%kdmin/par%px  < s%hvm(i,j)  )  ) then
                  nonhV(i,j) = 1
               else
                  nonhV(i,j) = 0
               endif
            enddo
         enddo
      else
         nonhV = 0
      endif
      !
      ! Determine if a velocity point will be included in the nonh pressure matrix, include if
      ! any of the surrounding velocity points is included.
      !
      if (s%ny>0) then
         do j=2,s%ny
            do i=2,s%nx
               if (max(nonhV(i,j),nonhV(i,j-1),nonhU(i,j),nonhU(i-1,j)) > 0) then
                  nonhZ(i,j) = 1
               else
                  nonhZ(i,j) = 0
               endif
            enddo
         enddo
      else
         do i=2,s%nx
            if (max(nonhU(i,1),nonhU(i-1,1)) > 0) then
               nonhZ(i,1) = 1
            else
               nonhZ(i,1) = 0
            endif
         enddo
      endif


      if (par%nhbreaker == 1) then
         reformfac = par%reformsteep/par%maxbrsteep
         if (s%ny>0) then
            do j=2,s%ny
               do i=2,s%nx
                  wmax = par%maxbrsteep*sqrt(par%g*s%hh(i,j))
                  if (s%breaking(i,j) == 0) then
                     if (s%ws(i,j)>=wmax) then
                        s%breaking(i,j) = 1
                     elseif (s%ws(i,j)<=-wmax) then
                        s%breaking(i,j) = -1
                     endif
                  elseif (s%breaking(i,j)==1) then
                     if (s%ws(i,j)<reformfac*wmax) then
                        s%breaking(i,j) = 0
                     endif
                  elseif (s%breaking(i,j)==-1) then
                     if (s%ws(i,j)>reformfac*(-wmax)) then
                        s%breaking(i,j) = 0
                     endif
                  endif
               enddo
            enddo
         else
            do i=2,s%nx
               wmax = par%maxbrsteep*sqrt(par%g*s%hh(i,1))
               if (s%breaking(i,1) == 0) then
                  if (s%ws(i,1)>=wmax) then
                     s%breaking(i,1) = 1
                  elseif (s%ws(i,1)<=-wmax) then
                     s%breaking(i,1) = -1
                  endif
               elseif (s%breaking(i,1)==1) then
                  if (s%ws(i,1)<reformfac*wmax) then
                     s%breaking(i,1) = 0
                  endif
               elseif (s%breaking(i,1)==-1) then
                  if (s%ws(i,1)>reformfac*(-wmax)) then
                     s%breaking(i,1) = 0
                  endif
               endif
            enddo
         endif
         ! turn off non-hydrostatic pressure correction in areas with breaking and increase viscosity
         where (s%breaking/=0)
            nonhZ = 0
            s%pres = 0
         endwhere
      elseif (par%nhbreaker == 2) then
         ! First determine local breaker criterion
         lbreakcond = par%maxbrsteep
         if (s%ny==0) then
            do i=2,s%nx
               if (s%breaking(i,1)==1) then
                  lbreakcond(i-1:i+1,1) = par%secbrsteep
               endif
            enddo
         else
            do j=jmin,jmax
               do i=2,s%nx
                  if (s%breaking(i,j)==1) then
                     lbreakcond(i-1:i+1,j-1:j+1) = par%secbrsteep
                  endif
               enddo
            enddo
         endif
#ifdef USEMPI
         call xmpi_shift_ee(lbreakcond)
#endif
         ! Now find areas where main breaking criterion is exceeded
         do j=jmin,jmax
            do i=2,s%nx
               if (s%breaking(i,j)==1) then
                  if(s%ws(i,j)<=0.d0) then
                     s%breaking(i,j) = 0
                  endif
               else
                  wmax = lbreakcond(i,j)*sqrt(par%g*s%hh(i,j))  ! add current term in here too
                  if (s%ws(i,j)>=wmax) then
                     s%breaking(i,j) = 1
                  endif
               endif
            enddo
         enddo
         ! turn off non-hydrostatic pressure correction in areas with breaking and increase viscosity
         do j=jmin,jmax
            do i=2,s%nx
               if (s%breaking(i,j)==1) then
                  nonhZ(i,j) = 0
                  s%pres(i,j) = 0.d0
               endif
            enddo
         enddo
      elseif (par%nhbreaker == 3) then
         ! First determine local breaker criterion
         lbreakcond = par%maxbrsteep
         if (s%ny==0) then
            do i=2,s%nx
               if (s%breaking(i,1)==1) then
                  lbreakcond(i-1:i+1,1) = par%secbrsteep
               endif
            enddo
         else
            do j=jmin,jmax
               do i=2,s%nx
                  if (s%breaking(i,j)==1) then
                     lbreakcond(i-1:i+1,j-1:j+1) = par%secbrsteep
                  endif
               enddo
            enddo
         endif
         ! Now find areas where main breaking criterion is exceeded
         do j=jmin,jmax
            do i=2,s%nx
               s%wscrit(i,j) = lbreakcond(i,j)*sqrt(par%g*s%hh(i,j))  ! add current term in here too
               if (s%breaking(i,j)==1) then
                  if(s%ws(i,j)<=0.d0) then
                     s%breaking(i,j) = 0
                  endif
               else
                  if (s%ws(i,j)>=s%wscrit(i,j)) then
                     s%breaking(i,j) = 1
                  endif
               endif
            enddo
         enddo
         ! turn off non-hydrostatic pressure correction in areas with breaking and increase viscosity
         do j=jmin,jmax
            do i=2,s%nx
               if (s%breaking(i,j)==1) then
                  nonhZ(i,j) = 0
                  s%pres(i,j) = 0.d0
               endif
            enddo
         enddo
      endif
      !
   end subroutine nonh_masks

   !
   !==============================================================================
   subroutine nonh_init_wcoef(s,par)
      !==============================================================================
      !

      ! DATE               AUTHOR               CHANGES
      !
      ! March 2013         Robert McCall        Move computation of wcoef to own function so callable separately

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------

      !  Recompute wcoef optimisation

      !--------------------------     DEPENDENCIES       ----------------------------
      use spaceparams
      use params
      use wave_functions_module, only: dispersion

      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars) ,intent(in)                          :: s
      type(parameters),intent(in)                          :: par

      !--------------------------     LOCAL VARIABLES    ----------------------------

      integer(kind=iKind)        :: i,j
      real   (kind=rKind)        :: d,sigma,k


      !
      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------
      !

      if (allocated(wcoef)) then
         if ( par%nonhq3d == 1 ) then
            wcoef = par%nhlay
         else
            sigma = 2.0_rkind*par%px/par%Topt
            if (par%dispc <= 0.0_rKind) then
               !Calculate the optimum value of the dispersion coeficiant.
               do j=1,s%ny+1
                  do i=1,s%nx+1
                     d = max(s%zs0(i,j)-s%zb(i,j),par%eps)
                     k = disper(sigma,d,par%g,2.0_rKind*par%px,0.0001_rKind)
                     if (d>0.0_rKind) then
                        wcoef(i,j) = 1.0_rKind / (  4.0_rKind*par%g / (d*sigma**2) - 4.0_rKind / ((k*d)**2)  )
                     else
                        wcoef(i,j) = 1.d0
                     endif
                  enddo
               enddo
            else
               wcoef = par%dispc
            endif
         endif
      endif


   end subroutine nonh_init_wcoef
   !




   subroutine zuzv(s)
      !==============================================================================
      !

      ! DATE               AUTHOR               CHANGES
      !
      ! November 2010       Pieter Bart Smit     New Subroutine


      !------------------------------------------------------------------------------
      !                             DECLARATIONS
      !------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------
      !
      !   Interpolate bottom and free surface location to u/v points
      !
      !--------------------------     DEPENDENCIES       ----------------------------
      use spaceparams
      use xmpi_module
      !--------------------------     ARGUMENTS          ----------------------------

      type(spacepars) ,intent(inout)                       :: s
      !--------------------------     LOCAL VARIABLES    ----------------------------

      integer(kind=iKind)                                  :: i
      integer(kind=iKind)                                  :: j

      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------
      !Bottom location in U point
      do j=1,s%ny+1
         do i=1,s%nx
            zbu(i,j) = max(s%zb(i,j),s%zb(min(s%nx,i)+1,j))
         enddo
      enddo
#ifdef USEMPI
      if (xmpi_isbot) then
         zbu(s%nx+1,:) = s%zb(s%nx+1,:)
      endif
#else
      zbu(s%nx+1,:) = s%zb(s%nx+1,:)
#endif
      !Bottom location in V point
      if (s%ny>2) then
         do j=1,s%ny
            do i=1,s%nx+1
               zbv(i,j) = max(s%zb(i,j),s%zb(i,min(s%ny,j)+1))
            enddo
         enddo
         zbv(:,s%ny+1) = s%zb(:,s%ny+1)
      else
         zbv = s%zb
      endif

      !Free surface location in u-point
      do j=1,s%ny+1
         do i=1,s%nx
            !
            zsu(i,j) = s%hu(i,j) + zbu(i,j)
            !
         enddo
      enddo

      !Free surface location in v-point
      if (s%ny>2) then
         do j=1,s%ny
            do i=1,s%nx+1
               zsv(i,j) = s%hv(i,j) + zbv(i,j)
            enddo
         enddo
      else
         do j=1,s%ny+1
            zsv(:,j) = s%zs(:,min(2,s%ny+1))
         enddo
      endif
      !
   end subroutine zuzv

   !
   !==============================================================================
   real(kind=rKind) function disper(w,d,g,pi2,accuracy)
      !==============================================================================
      !


      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !--------------------------        PURPOSE         ----------------------------

      !  Calculate k for a given intrinsic frequency w and depth d. First use the fenton
      !  approximation and then iterate for better accuracy (when necessary)


      !--------------------------     ARGUMENTS          ----------------------------

      real(kind=rKind),intent(in) :: w
      real(kind=rKind),intent(in) :: d
      real(kind=rKind),intent(in) :: g
      real(kind=rKind),intent(in) :: pi2
      real(kind=rKind),intent(in) :: accuracy

      !--------------------------     LOCAL VARIABLES    ----------------------------

      real(kind=rKind) :: k
      real(kind=rKind) :: alpha
      real(kind=rKind) :: L,Ldeep,pi2d,w2
      real(kind=rKind) :: error

      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------
      w2    = w**2
      alpha = d*w2/g
      k     = alpha*(1.0_rKind/sqrt(tanh(alpha)))/d
      L    = pi2/k
      Ldeep = g*pi2/(w2)
      pi2d     = pi2*d

      error = abs(g*k*tanh(k*d)-w2)/w2
      do while (error > accuracy)
         L    = Ldeep*tanh(pi2d/L)
         k     = pi2/L
         error = abs(g*k*tanh(k*d)-w2)/w2
      enddo

      disper = k
   end function disper


end module nonh_module
