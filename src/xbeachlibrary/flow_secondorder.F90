!==============================================================================
!                               FLOW_SECONDORDER_MODULE
!==============================================================================
!
! DATE               AUTHOR               CHANGES
!
! october 2009       Pieter Bart Smit     New module
!
!******************************************************************************
!                                 INTERFACE
!******************************************************************************

module flow_secondorder_module

   implicit none
   save
   private

   !----------------------------- PARAMETERS -----------------------------------

   include 'nh_pars.inc'               !Default precision etc.

   !----------------------------- VARIABLES  -----------------------------------


   real(kind=rKind),dimension(:,:),allocatable  :: wrk1 !Temporary storage for:
   !    (i )  du-velocity at i    ,j+1/2
   !    (ii)  dv-velocity at i    ,j+1/2
   !    (iii) dqx         at i+1/2,j
   real(kind=rKind),dimension(:,:),allocatable  :: wrk2 !Temporary storage for:
   !    (i )  du-velocity at i+1/2,j
   !    (ii)  dv-velocity at i+1/2,j
   !    (iii) dqx         at i    ,j+1/2

   logical                                      :: initialized = .false.

   public flow_secondorder_advUV
   public flow_secondorder_advW
   public flow_secondorder_con
   public minmod
   public  flow_secondorder_huhv
   !
   !******************************************************************************
   !                             SUBROUTINES/FUNCTIONS
   !******************************************************************************

contains

   subroutine flow_secondorder_init(s)
      !

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !
      !--------------------------        PURPOSE         ----------------------------
      !
      !   Initializes the resources needed for the second order mcCormack scheme
      !
      !
      !--------------------------     DEPENDENCIES       ----------------------------
      !
      !
      ! -- MODULES --
      use spaceparams
      use params

      !
      !--------------------------     ARGUMENTS          ----------------------------
      !
      type(spacepars)    ,intent(inout) :: s

      !Allocate resources
      allocate (  wrk1(s%nx+1,s%ny+1))
      allocate (  wrk2(s%nx+1,s%ny+1))

      wrk1   = 0.0_rKind
      wrk2   = 0.0_rKind

      initialized = .true.
   end subroutine flow_secondorder_init

   !
   !==============================================================================
   subroutine flow_secondorder_advUV(s,par,uu_old,vv_old)
      !==============================================================================
      !

      ! DATE               AUTHOR               CHANGES
      !
      ! November 2010       Pieter Bart Smit     New Subroutine
      ! November 2014       Pieter Bart Smit     Updated for curvilinear and mpi

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !
      !--------------------------        PURPOSE         ----------------------------
      !
      !   Calculates second order correction to the advection terms for U and V
      !

      !
      !--------------------------     DEPENDENCIES       ----------------------------
      !
      use spaceparams
      use params

      !
      !--------------------------      METHOD            ----------------------------
      !
      ! First a prediction for the velocity on the new timelevel is done in flow_timestep
      ! using an explicit first order euler step. We use this prediction (s%uu,s%vv), and the velocity
      ! at the beginning of the step (uu_old, vv_old), to do a slope limited correction.
      !

      !
      !--------------------------     ARGUMENTS          ----------------------------
      !
      type(spacepars)  ,intent(inout)  :: s
      type(parameters) ,intent(in)     :: par

      real(kind=rkind),dimension(s%nx+1,s%ny+1),intent(in) :: vv_old  !"Old" v-velocity
      real(kind=rkind),dimension(s%nx+1,s%ny+1),intent(in) :: uu_old  !"Old" u-velocity
      !

      !--------------------------     LOCAL VARIABLES    ----------------------------
      !
      integer(kind=iKind)        :: iww  !i-2 :Two points 'west' of current point i
      integer(kind=iKind)        :: iw   !i-1 :One ...
      integer(kind=iKind)        :: i    !Current point
      integer(kind=iKind)        :: ie   !i+1
      integer(kind=iKind)        :: iee  !i+2
      integer(kind=iKind)        :: jnn  !j-2
      integer(kind=iKind)        :: jn   !j-1
      integer(kind=iKind)        :: j    !Current point
      integer(kind=iKind)        :: js   !j+1
      integer(kind=iKind)        :: jss  !j+2
      integer(kind=iKind)        :: jmin1d,jmax1d,jmin,jmax,imin,imax  ! index for superfast1D

      real(kind=rKind)           :: mindepth ! Near the dry/wet interface the higher order interpolations
      ! can cause unwanted effects. To avoid this any extrapolation/interpolation
      ! is disabled when the lowest surface elevation within the molecule is lower then
      ! the highers bottom elevation.

      real(kind=rKind)           :: delta1   ! "Central" difference
      real(kind=rKind)           :: delta2   ! "Upwind"  difference

      real(kind=rKind)           :: qx,qy    !discharge south
      real(kind=rKind)           :: fac      !A factor

      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------

      !Initialize/allocate arrays on first entry
      if (.not. initialized) then
         call flow_secondorder_init(s)
      endif

      !
      imin = 2
      if (xmpi_istop) then
         !
         imin = 3
         !
      endif

      imax = s%nx
      if (xmpi_isbot) then
         !
         imax = s%nx-1
         !
      endif

      jmin = 2
      if (xmpi_isleft) then
         !
         jmin = 3
         !
      endif

      jmax = s%ny
      if (xmpi_isright) then
         !
         jmax = s%ny-1
         !
      endif

      !
      if (s%ny>0) then
         !
         jmin1d = 2
         jmax1d = s%ny
         !
      else
         !
         jmin1d = 1
         jmax1d = 1
         !
      endif

      !
      !-- calculate qx * du in z-points --
      !
      do j=jmin1d,jmax1d
         !
         do i=2,s%nx
            !
            wrk1(i,j) = 0.0_rKind
            ie   = min(i+1,imax+1)
            iee  = min(i+2,imax+1)
            iw   = max(i-1,1)
            iww  = max(i-2,1)
            mindepth = minval(s%zs(iww:iee,j))-maxval(s%zb(iww:iee,j))

            qx = .5 * ( s%qx(i,j) + s%qx(iw,j) ) * cos( s%alfau( i ,j ) - s%alfau( iw,j) )
            !
            if (mindepth > par%eps) then
               !
               if   ( qx > 0.d0 .and. i>imin )    then
                  !
                  delta1    = (s%uu(i,j) -uu_old(iw ,j)) / s%dsz(i,j)
                  delta2    = (s%uu(iw,j)-uu_old(iww,j)) / s%dsz(iw,j)
                  wrk1(i,j) = 0.5d0*s%dsu(iw,j) * minmod(delta1,delta2) * qx
                  !
               elseif ( qx < 0.d0 .and. i < imax + 1  ) then
                  !
                  delta1    = (uu_old(i,j) -s%uu(iw ,j)) / s%dsz(i,j)
                  delta2    = (uu_old(ie,j)-s%uu(i,j))   / s%dsz(ie,j)
                  wrk1(i,j) = -0.5d0*s%dsu(i,j) * minmod(delta1,delta2) *qx
                  !
               endif
               !
            endif
            !
         enddo
         !
      enddo

      !
      !-- calculate qy * du in c-points --
      !
      if ( s%ny > 0 ) then
         !
         do j=2,s%ny
            !
            do i=2,s%nx-1
               !
               wrk2(i,j) = 0.0_rkind
               js   = min(j+1,jmax+1)
               jss  = min(j+2,jmax+1)
               jn   = max(j-1,1)
               mindepth = minval(s%zs(i:i+1,jn:jss))-maxval(s%zb(i:i+1,jn:jss))
               !
               if ((mindepth > par%eps)) then
                  !
                  qy = ( s%qy(i+1,j) + s%qy(i,j) ) * .5d0 * cos( s%alfau( i , j ) - s%alfau( i , j + 1 ) )
                  !
                  if ( qy  > 0.d0 .and. j > jmin - 1 ) then
                     !
                     delta1    = (s%uu(i,js) -uu_old(i ,j )) / s%dnv(i,j)
                     delta2    = (s%uu(i,j)  -uu_old(i ,jn)) / s%dnv(i,jn)
                     wrk2(i,j) = 0.5d0*s%dnv(i,j)*minmod(delta1,delta2) * qy
                     !
                  elseif ( qy < 0.d0 .and. j < jmax ) then
                     !
                     delta1    = (uu_old(i,js) -s%uu(i ,j)) / s%dnv(i,j)
                     delta2    = (uu_old(i,jss)-s%uu(i,js)) / s%dnv(i,js)
                     wrk2(i,j)   = -0.5d0*s%dnv(i,j)*minmod(delta1,delta2) * qy
                     !
                  endif
               endif
               !
            enddo
            !
         enddo
         !
         wrk2(:,1   ) = 0.0_rKind
         !
      endif

      !CORRECTION TO U
      if (s%ny>0) then
         !
         do j=2,s%ny
            !
            do i=2,s%nx-1
               !
               fac = par%dt / s%hum(i,j) * s%dsdnui(i,j)
               s%uu(i,j) = s%uu(i,j) - fac * ( s%dnz(i+1,j) *  wrk1(i+1,j) - s%dnz(i,  j) *  wrk1(i  ,j  ) &
               +   s%dsc(i,j  ) *  wrk2(i  ,j) - s%dsc(i,j-1) *  wrk2(i  ,j-1) )
               !
            enddo
            !
         enddo
         !
      else
         !
         do i=2,s%nx-1
            !
            s%uu(i,1) = s%uu(i,1)-par%dt/s%hum(i,1)*(  ( wrk1(i+1,1) - wrk1(i  ,1  ) ) / s%dsu(i,1) )
            !
         enddo
         !
      endif

      !== SECOND ORDER EXPLICIT CORRECTION TO V ==
      if (s%ny>2) then
         !
         do j=2,s%ny
            !
            do i=2,s%nx
               !
               wrk1(i,j) = 0.
               js   = min(j+1,jmax+1)
               jss  = min(j+2,jmax+1)
               jn   = max(j-1,1)
               jnn  = max(j-2,1)
               mindepth = minval(s%zs(i,jnn:jss))-maxval(s%zb(i,jnn:jss))
               !
               qy = 0.5d0 * ( s%qy(i,j)+s%qy(i,jn) ) * cos( s%alfav(i,j) - s%alfav(i,jn) ) * s%dsz(i,j)
               !
               if (mindepth > par%eps) then
                  !
                  if   (( qy > 0.0) .and. (j>jmin)) then
                     !
                     delta1    = (s%vv(i,j ) - vv_old(i ,jn )) / s%dnz(i,j)
                     delta2    = (s%vv(i,jn) - vv_old(i ,jnn)) / s%dnz(i,jn)
                     wrk1(i,j)   = 0.5d0*s%dnv(i,jn)*minmod(delta1,delta2)*qy
                     !
                  elseif ( qy < -0.0d0 .and. j<jmax+1) then
                     !
                     delta1    = (vv_old(i,j ) - s%vv(i,jn)) / s%dnz(i,j)
                     delta2    = (vv_old(i,js) - s%vv(i,j )) / s%dnz(i,js)
                     wrk1(i,j)   = -0.5d0*s%dnv(i,j)*minmod(delta1,delta2)*qy
                     !
                  endif
                  !
               endif
            enddo
         enddo

         do j=2,s%ny-1
            !
            do i=2,s%nx
               !
               wrk2(i,j) = 0.0d0
               ie   = min(i+1,s%nx)
               iee  = min(i+2,s%nx)
               iw   = max(i-1,1)
               mindepth = minval(s%zs(iw:iee,j:j+1))-maxval(s%zb(iw:iee,j:j+1))
               !
               if (mindepth > par%eps) then
                  !
                  qx = .5d0 * ( s%qx(i,j+1) + s%qx(i,j) ) * cos( s%alfav( i , j ) - s%alfav( i + 1 , j ) ) * s%dsc(i,j)
                  !
                  if     ( qx > par%Umin .and. i > imin - 1) then
                     delta1    = (s%vv(ie,j) - vv_old(i ,j )) / s%dsu(i,1)
                     delta2    = (s%vv(i ,j) - vv_old(iw,j))  / s%dsu(iw,1)
                     wrk2(i,j) = 0.5d0*s%dsu(i,1)*minmod(delta1,delta2)*qx
                  elseif ( qx < -par%Umin .and. i < imax ) then
                     delta1    = (vv_old(ie,j) -s%vv(i ,j)) / s%dsu(ie,1)
                     delta2    = (vv_old(iee,j)-s%vv(ie,j)) / s%dsu(ie,1)
                     wrk2(i,j) = -0.5d0*s%dsu(i,1)*minmod(delta1,delta2)*qx
                  endif
               endif
            enddo
         enddo
         wrk2(1,:)    = 0.0_rKind

         !CORRECTION TO V
         do j=2,s%ny-1
            !
            do i=2,s%nx
               !
               fac = par%dt / s%hvm(i,j) * s%dsdnvi(i,j)
               s%vv(i,j) = s%vv(i,j) - fac*( wrk2(i,j) - wrk2(i-1,j) + wrk1(i,j+1) - wrk1(i  ,j) )
               !
            enddo
            !
         enddo
      endif
      !
   end subroutine flow_secondorder_advUV

   !
   !==============================================================================
   subroutine flow_secondorder_advW(s,par,w,w_old,mskZ)
      !==============================================================================
      !

      ! DATE               AUTHOR               CHANGES
      !
      ! November 2010       Pieter Bart Smit     New Subroutine
      ! November 2014       Pieter Bart Smit     Modified for nonhq3d

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !
      !--------------------------        PURPOSE         ----------------------------
      !
      !   Calculates second order correction to the advection terms for W. Only used
      !   in combination WITH the non-hydrostatic module

      !--------------------------     DEPENDENCIES       ----------------------------
      !

      ! -- MODULES --
      use spaceparams
      use params
      use xmpi_module


      !--------------------------     ARGUMENTS          ----------------------------
      !

      type(spacepars)  ,intent(inout)  :: s
      type(parameters) ,intent(in)     :: par

      real(kind=rKind),dimension(s%nx+1,s%ny+1),intent(inout) :: w      !The velocity at the star level
      real(kind=rKind),dimension(s%nx+1,s%ny+1),intent(in)    :: w_old  !The vertical velocity at the n-level
      integer(kind=iKind),dimension(s%nx+1,s%ny+1),intent(in) :: mskZ   !A mask that determines whether or not the w-equation is active

      !
      !--------------------------     LOCAL VARIABLES    ----------------------------
      !
      integer(kind=iKind)        :: iw   !i-1 :One ...
      integer(kind=iKind)        :: i    !Current point
      integer(kind=iKind)        :: ie   !i+1
      integer(kind=iKind)        :: iee  !i+2
      integer(kind=iKind)        :: jn   !j+1
      integer(kind=iKind)        :: jnn   !j+2
      integer(kind=iKind)        :: j    !Current point
      integer(kind=iKind)        :: js   !j-1
      integer(kind=iKind)        :: jmin,jmax,imin,imax

      real(kind=rKind)           :: rfac
      real(kind=rKind)           :: mindepth
      real(kind=rKind)           :: delta1
      real(kind=rKind)           :: delta2,qx,qy

      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------


      !
      imin = 1
      if (xmpi_istop) then
         !
         imin = 2
         !
      endif

      imax = s%nx
      if (xmpi_isbot) then
         !
         imax = s%nx-1
         !
      endif

      jmin = 1
      if (xmpi_isleft) then
         !
         jmin = 2
         !
      endif

      jmax = s%ny
      if (xmpi_isright) then
         !
         jmax = s%ny-1
         !
      endif
      !

      !Initialize/allocate arrays on first entry
      if (.not. initialized) then
         !
         call flow_secondorder_init(s)
         !
      endif

      if (s%ny > 0) then
         !
         ! 2D Codepath
         !

         !
         ! Calculate the updated fluxes of w-mom. in the s-dir
         !
         wrk1 = 0.0d0
         do j=2,s%ny
            !
            do i=1,s%nx
               !
               if ( mskZ(i,j) == 0 ) cycle

               wrk1(i,j) = 0.0_rKind
               !
               ie   = min(i+1,imax+1)
               iee  = min(i+2,imax+1)
               iw   = max(i-1,1)

               rfac = real( mskZ( ie , j ) * mskZ( iee , j ) * mskZ( iw , j ), kind=rKind )

               mindepth = minval(s%zs(iw:iee,j))-maxval(s%zb(iw:iee,j))
               qx = s%qx(i,j) * s%dnu( i,j ) * rfac
               !
               if (mindepth > par%eps) then
                  !
                  if   ( qx > 0.0_rKind .and. i>imin ) then
                     !
                     delta1    = (w(ie,j ) - w_old(i  ,j )) / s%dsu(i,j)
                     delta2    = (w(i,j )  - w_old(iw ,j )) / s%dsu(iw,j)
                     wrk1(i,j)   = 0.5d0*s%dsu(i,j)*minmod(delta1,delta2) * qx
                     !
                  elseif (qx < 0.0_rKind .and. i<imax) then
                     !
                     delta1    = (w_old(ie ,j)  - w(i ,j )) / s%dsu(i,j)
                     delta2    = (w_old(iee,j ) - w(ie,j )) / s%dsu(ie,j)
                     wrk1(i,j) = -0.5d0*s%dsu(i,j)*minmod(delta1,delta2)*qx
                     !
                  endif
                  !
               endif
               !
            enddo
            !
         enddo

         !
         ! Calculate the updated fluxes of w-mom. in the n-dir
         !
         wrk2 = 0.
         do j=1,s%ny
            !
            do i=2,s%nx
               !
               if ( mskZ(i,j) == 0 ) cycle

               wrk2(i,j) = 0.0_rKind
               !
               jn   = min(j+1,jmax+1)
               jnn  = min(j+2,jmax+1)
               js   = max(j-1,1)

               rfac = real( mskZ( i , jn ) * mskZ( i , jnn ) * mskZ( i , js ), kind=rKind )

               mindepth = minval(s%zs(i,js:jnn))-maxval(s%zb(i,js:jnn))
               qy = s%qy(i,j) * s%dsv( i,j ) * rfac
               !
               if (mindepth > par%eps) then
                  !
                  if   ( qy > 0.0_rKind .and. j>jmin ) then
                     !
                     delta1    = (w(i,jn ) - w_old(i  ,j )) / s%dnv(i,j)
                     delta2    = (w(i,j )  - w_old(i ,js )) / s%dnv(i,js)
                     wrk2(i,j)   = 0.5d0*s%dnv(i,j)*minmod(delta1,delta2) * qy
                     !
                  elseif (qy < 0.0_rKind .and. j<jmax) then
                     !
                     delta1    = (w_old(i ,jn)  - w(i ,j )) / s%dnv(i,j)
                     delta2    = (w_old(i,jnn ) - w(i,jn )) / s%dnv(i,jn)
                     wrk2(i,j) = -0.5d0*s%dnv(i,j)*minmod(delta1,delta2)*qy
                     !
                  endif
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
            wrk1(1,:) = 0.0d0
            !
         endif
         if (xmpi_isbot) then
            !
            wrk1(s%nx,:) = 0.0d0
            !
         endif
         if (xmpi_isleft) then
            !
            wrk2(:,1) = 0.0d0
            !
         endif
         if (xmpi_isright) then
            !
            wrk2(:,s%ny) =0.
            !
         endif
         !

         !
         ! Update w-mom
         !
         do j=2,s%ny
            !
            do i=2,s%nx
               !
               if ( mskZ(i,j) == 0 ) cycle
               !
               w(i,j) = w(i,j) - par%dt * s%dsdnzi(i,j) *( wrk1(i,j) + wrk2(i,j) - wrk1(i-1,j) - wrk2(i,j-1) ) /s%hh(i,j)
               !
            enddo
            !
         enddo
         !
      else
         !
         ! 1D Codepath
         !
         j=1
         do i=2,s%nx
            !
            if ( mskZ(i,j) == 0 ) cycle

            wrk1(i,j) = 0.0_rKind
            !
            ie   = min(i+1,s%nx)
            iee  = min(i+2,s%nx)
            iw   = max(i-1,1)

            rfac = real( mskZ( ie , j ) * mskZ( iee , j ) * mskZ( iw , j ), kind=rKind )

            mindepth = minval(s%zs(iw:iee,j))-maxval(s%zb(iw:iee,j))
            qx = s%qx(i,j) * rfac
            !
            if (mindepth > par%eps) then
               !
               if   ( qx > 0.0_rKind .and. i>2 ) then
                  !
                  delta1    = (w(ie,j ) - w_old(i  ,j )) / s%dsu(i,j)
                  delta2    = (w(i,j )  - w_old(iw ,j )) / s%dsu(iw,j)
                  wrk1(i,j)   = 0.5d0*s%dsu(i,j)*minmod(delta1,delta2) * qx
                  !
               elseif (qx < 0.0_rKind .and. i<s%nx-1) then
                  !
                  delta1    = (w_old(ie ,j)  - w(i ,j )) / s%dsu(i,j)
                  delta2    = (w_old(iee,j ) - w(ie,j )) / s%dsu(ie,j)
                  wrk1(i,j)   = -0.5d0*s%dsu(i,j)*minmod(delta1,delta2)*qx
                  !
               endif
               !
            endif
            !
         enddo

         !
         if (xmpi_istop) then
            !
            wrk1(1,:) = 0.0d0
            !
         endif
         if (xmpi_isbot) then
            !
            wrk1(s%nx,:) = 0.0d0
            !
         endif
         !
         ! Update w-mom to ** level
         !
         do i=2,s%nx
            !
            if ( mskZ(i,1) == 0 ) cycle
            !
            w(i,1) = w(i,1)-par%dt/s%hh(i,1)* ( wrk1(i,1)- wrk1(i-1,1) )/ s%dsz(i,1)
            !
         enddo
         !
      endif
      !
   end subroutine flow_secondorder_advW

   !
   !==============================================================================
   subroutine flow_secondorder_con(s,par,zs_old)
      !==============================================================================
      !
      ! DATE               AUTHOR               CHANGES
      !
      ! November 2010       Pieter Bart Smit     New Subroutine

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !
      !--------------------------        PURPOSE         ----------------------------
      !
      !   Calculates second order correction to the continuity terms
      !

      !--------------------------     DEPENDENCIES       ----------------------------
      !
      use spaceparams
      use params

      !--------------------------     ARGUMENTS          ----------------------------
      !
      type(spacepars)  ,intent(inout)  :: s
      type(parameters) ,intent(in)     :: par

      real(kind=rKind),dimension(s%nx+1,s%ny+1),intent(in) :: zs_old
      !

      !--------------------------     LOCAL VARIABLES    ----------------------------
      !

      integer(kind=iKind)        :: iw   !i-1 :One ...
      integer(kind=iKind)        :: i    !Current point
      integer(kind=iKind)        :: ie   !i+1
      integer(kind=iKind)        :: iee  !i+2
      integer(kind=iKind)        :: jn   !j-1
      integer(kind=iKind)        :: j    !Current point
      integer(kind=iKind)        :: js   !j+1
      integer(kind=iKind)        :: jss  !j+2
      integer(kind=iKind)        :: jmin,jmax  !index for superfast1D

      real(kind=rKind)           :: mindepth
      real(kind=rKind)           :: delta1
      real(kind=rKind)           :: delta2
      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------

      if (s%ny>0) then
         jmin = 2
         jmax = s%ny
      else
         jmin = 1
         jmax = 1
      endif

      !== SECOND ORDER EXPLICIT CORRECTION TO CONTINUITY ==
      !Initialize/allocate arrays on first entry
      if (.not. initialized) then
         call flow_secondorder_init(s)
      endif
      !return
      !correction to mass flux qx
      do j=jmin,jmax
         do i=2,s%nx-1
            wrk1(i,j) = 0.0_rkind
            ie  = min(i+1,s%nx)
            iee = min(i+2,s%nx)
            iw  = max(i-1,2)
            mindepth = minval(s%zs(iw:iee,j))-maxval(s%zb(iw:iee,j))
            if (mindepth>par%eps) then
               if     (s%uu(i,j) >  par%umin .and. i>2     ) then
                  delta1    =  (s%zs(ie,j)- zs_old(i,j ))/ s%dsu(i,1)
                  delta2    =  (s%zs(i,j) - zs_old(iw,j))/ s%dsu(i-1,1)
                  wrk1(i,j)  =   s%uu(i,j)*0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
               elseif (s%uu(i,j) < -par%umin .and. i<s%nx-1) then
                  delta1    =  (zs_old(iee,j) - s%zs(ie,j)) / s%dsu(i+1,1)
                  delta2    =  (zs_old(ie ,j) - s%zs(i ,j)) / s%dsu(i,1)
                  wrk1(i,j)  = - s%uu(i,j)*0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
               endif
            endif
         enddo
      enddo
      wrk1(1   ,:) = 0.0_rkind
      wrk1(s%nx,:) = 0.0_rkind

      !Correction to mass flux qy
      if (s%ny > 2) then
         do j=2,s%ny-1
            do i=2,s%nx
               wrk2(i,j) = 0.0_rKind
               js  = min(j+1,s%ny)
               jss = min(j+2,s%ny)
               jn  = max(j-1,2)
               mindepth = minval(s%zs(i,jn:jss))-maxval(s%zb(i,jn:jss))
               if (mindepth> par%eps) then
                  if     (s%vv(i,j) >  par%Umin .and. j>2) then
                     delta1    = (s%zs(i,js) - zs_old(i,j ))/ s%dnv(1,j)
                     delta2    = (s%zs(i,j)  - zs_old(i,jn))/ s%dnv(1,j-1)
                     wrk2(i,j) =  s%vv(i,j)*0.5d0*s%dnv(1,j)*minmod(delta1,delta2)
                  elseif (s%vv(i,j) < -par%Umin .and. j<s%ny-1) then
                     delta1    = (zs_old(i,jss) - s%zs(i,js)) / s%dnv(1,j+1)
                     delta2    = (zs_old(i ,js) - s%zs(i ,j)) / s%dnv(1,j)
                     wrk2(i,j) =  s%vv(i,j)*0.5d0*s%dnv(1,j)*minmod(delta1,delta2)
                  endif
               endif
            enddo
         enddo
         wrk2(:,1   ) = 0.0_rKind
         wrk2(:,s%ny) = 0.0_rKind
      else
         wrk2 = 0.0_rKind
      endif

      !Update waterlevels
      if (s%ny>0) then
         do j=2,s%ny
            do i=2,s%nx
               s%zs(i,j) = s%zs(i,j)-par%dt*(  (wrk1(i,j)-wrk1(i-1,j))/ s%dsz(i,1)  &
               +  (wrk2(i,j)-wrk2(i,j-1))/ s%dnz(1,j)  )
            enddo
         enddo
      else
         do i=2,s%nx
            s%zs(i,1) = s%zs(i,1)-par%dt*(  (wrk1(i,1)-wrk1(i-1,1))/ s%dsz(i,1) )
         enddo
      endif

      !Update fluxes
      if (s%ny>0) then
         s%qx(2:s%nx,2:s%ny) = s%qx(2:s%nx,2:s%ny) +wrk1(2:s%nx,2:s%ny)
         s%qy(2:s%nx,2:s%ny) = s%qy(2:s%nx,2:s%ny) +wrk2(2:s%nx,2:s%ny)
      else
         s%qx(2:s%nx,1) = s%qx(2:s%nx,1) +wrk1(2:s%nx,1)
      endif

   end subroutine flow_secondorder_con

   !
   !==============================================================================
   subroutine flow_secondorder_huhv(s,par)
      !==============================================================================
      !
      ! DATE               AUTHOR               CHANGES
      !
      ! November 2014       Pieter Bart Smit     New Subroutine

      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !
      !--------------------------        PURPOSE         ----------------------------
      !
      !   Calculates second order correction to the waterlevels
      !

      !--------------------------     DEPENDENCIES       ----------------------------
      !
      use spaceparams
      use params
      use xmpi_module

      !--------------------------     ARGUMENTS          ----------------------------
      !
      type(spacepars)  ,intent(inout)  :: s
      type(parameters) ,intent(in)     :: par

      !    real(kind=rKind),dimension(s%nx+1,s%ny+1),intent(in) :: zs_old
      !

      !--------------------------     LOCAL VARIABLES    ----------------------------
      !

      integer(kind=iKind)        :: iw   !i-1 :One ...
      integer(kind=iKind)        :: i    !Current point
      integer(kind=iKind)        :: ie   !i+1
      integer(kind=iKind)        :: iee  !i+2
      integer(kind=iKind)        :: j    !Current point
      integer(kind=iKind)        :: jmin,jmax,imin,imax
      real(kind=rKind)           :: mindepth
      real(kind=rKind)           :: delta1
      real(kind=rKind)           :: delta2
      !-------------------------------------------------------------------------------
      !                             IMPLEMENTATION
      !-------------------------------------------------------------------------------

      !
      imin = 1
      if (xmpi_istop) then
         !
         imin = 2
         !
      endif

      imax = s%nx
      if (xmpi_isbot) then
         !
         imax = s%nx-1
         !
      endif

      jmin = 1
      if (xmpi_isleft) then
         !
         jmin = 2
         !
      endif

      jmax = s%ny
      if (xmpi_isright) then
         !
         jmax = s%ny-1
         !
      endif

      if (s%ny > 0) then
         !
         ! 2DH Code
         !
         do j= 2 , s%ny
            !
            do i = 1 , s%nx
               !
               ie  = min(i+1,imax+1)
               iee = min(i+2,imax+1)
               iw  = max(i-1,1)
               !
               mindepth = minval(s%zs(iw:iee,j))-maxval(s%zb(iw:iee,j))
               !
               if (mindepth>par%eps) then
                  !
                  if     (s%uu(i,j) >  par%umin .and. i > imin     ) then
                     !
                     delta1    =  (s%zs(ie,j)- s%zs(i,j ))/ s%dsu(i,1)
                     delta2    =  (s%zs(i,j) - s%zs(iw,j))/ s%dsu(i-1,1)
                     s%hu(i,j)  = s%hu(i,j) +  0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
                     !
                  elseif (s%uu(i,j) < -par%umin .and. i < imax ) then
                     !
                     delta1    =  (s%zs(iee,j) - s%zs(ie,j)) / s%dsu(i+1,1)
                     delta2    =  (s%zs(ie ,j) - s%zs(i ,j)) / s%dsu(i,1)
                     s%hu(i,j)  = s%hu(i,j) - 0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
                     !
                  endif
                  !
               endif
               !
            enddo
            !
         enddo

         do j= 1,s%ny
            !
            ie  = min(j+1,jmax+1)
            iee = min(j+2,jmax+1)
            iw  = max(j-1,1)
            !
            do i=2,s%nx
               !
               mindepth = minval(s%zs(i,iw:iee))-maxval(s%zb(i,iw:iee))
               !
               if (mindepth>par%eps) then
                  !
                  if     (s%vv(i,j) >  par%umin .and. j>jmin     ) then
                     !
                     delta1    =  (s%zs(i,ie)- s%zs(i,j ))/ s%dnv(i,j)
                     delta2    =  (s%zs(i,j) - s%zs(i,iw))/ s%dnv(i,j-1)
                     s%hv(i,j)  = s%hv(i,j) +  0.5d0*s%dnv(i,j)*minmod(delta1,delta2)
                     !
                  elseif (s%vv(i,j) < -par%umin .and. j<jmax) then
                     !
                     delta1    =  (s%zs(i,iee) - s%zs(i,ie)) / s%dnv(i,j+1)
                     delta2    =  (s%zs(i ,ie) - s%zs(i ,j)) / s%dnv(i,j)
                     s%hv(i,j)  = s%hv(i,j) - 0.5d0*s%dnv(i,j)*minmod(delta1,delta2)
                     !
                  endif
                  !
               endif
               !
            enddo
         enddo
      else
         !
         ! 1DH Code
         !
         j= 1
         do i=1,s%nx-1
            !
            ie  = min(i+1,imax+1)
            iee = min(i+2,imax+1)
            iw  = max(i-1,1)
            mindepth = minval(s%zs(iw:iee,j))-maxval(s%zb(iw:iee,j))
            !
            if (mindepth>par%eps) then
               !
               if     (s%uu(i,j) >  par%umin .and. i>imin    ) then
                  !
                  delta1    =  (s%zs(ie,j)- s%zs(i,j ))/ s%dsu(i,1)
                  delta2    =  (s%zs(i,j) - s%zs(iw,j))/ s%dsu(i-1,1)
                  s%hu(i,j)  = s%hu(i,j) +  0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
                  !
               elseif (s%uu(i,j) < -par%umin .and. i<imax) then
                  !
                  delta1    =  (s%zs(iee,j) - s%zs(ie,j)) / s%dsu(i+1,1)
                  delta2    =  (s%zs(ie ,j) - s%zs(i ,j)) / s%dsu(i,1)
                  s%hu(i,j)  = s%hu(i,j) - 0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
               endif
               !
            endif
            !
         enddo
         !
      endif

   end subroutine flow_secondorder_huhv
   !
   !==============================================================================
   real(kind=rKind) pure function minmod(delta1,delta2)
      !==============================================================================
      !
      ! DATE               AUTHOR               CHANGES
      !
      ! November 2010       Pieter Bart Smit     New Subroutine
      !
      !-------------------------------------------------------------------------------
      !                             DECLARATIONS
      !-------------------------------------------------------------------------------

      !
      !--------------------------        PURPOSE         ----------------------------
      !
      !   Limits succesive gradients using the TVD minmod limiter.
      !
      !--------------------------        METHOD          ----------------------------
      !
      ! - MINMOD LIMITER -
      !
      ! The minmod slope limiter essentiallycorrects the first order upwind
      ! estimate for hu/u/v with a second order upwind or central approximation.
      ! If the succesive gradients are opposite in sign the correction is set to zero to
      ! avoid unwanted oscillations.
      !
      ! Note: Declared as pure to explicitly denote that the function has no side effects.
      !       If this is not supported in the compiler the pure attribute can be safely removed.
      !
      ! ( see for instance Hirch [2001] for a better introduction to limiters. )

      !
      !--------------------------     DEPENDENCIES       ----------------------------
      !
      !                                 - NONE -
      !
      !--------------------------     ARGUMENTS          ----------------------------
      !
      real(kind=rKind),intent(in) :: delta1  ! first gradient
      real(kind=rKind),intent(in) :: delta2  ! secpond gradient
      !
      !--------------------------     LOCAL VARIABLES    ----------------------------
      !
      if (delta1*delta2 <= 0.0_rKind) then
         minmod = 0.0_rKind
         return
      endif

      if     (delta1 > 0.0_rKind) then
         minmod = min(delta1,delta2)
      elseif (delta1  < 0.0_rKind) then
         minmod = max(delta1,delta2)
      else
         minmod = 0.0_rKind
      endif
   end function minmod


end module flow_secondorder_module
