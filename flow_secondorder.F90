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
    
    real(kind=rKind)           :: mindepth ! Near the dry/wet interface the higher order interpolations
                                           ! can cause unwanted effects. To avoid this any extrapolation/interpolation
                                           ! is disabled when the lowest surface elevation within the molecule is lower then 
                                           ! the highers bottom elevation.

    real(kind=rKind)           :: delta1   ! "Central" difference
    real(kind=rKind)           :: delta2   ! "Upwind"  difference   
    
    real(kind=rKind)           :: qe       !discharge east
    real(kind=rKind)           :: qw       !discharge west
    real(kind=rKind)           :: qn       !discharge north
    real(kind=rKind)           :: qs       !discharge south

!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------

    !== SECOND ORDER EXPLICIT CORRECTION TO U ==
    
    !Initialize/allocate arrays on first entry
    if (.not. initialized) then
      call flow_secondorder_init(s)
    endif  
    
    !-- calculate du in z-points --
    do j=2,s%ny
      do i=2,s%nx
        wrk1(i,j) = 0.0_rKind
        ie   = min(i+1,s%nx)
        iee  = min(i+2,s%nx)
        iw   = max(i-1,1)
        iww  = max(i-2,1)
        mindepth = minval(s%zs(iww:iee,j))-maxval(s%zb(iww:iee,j))
        if (mindepth > par%eps) then
          if   ((s%qx(i,j)+s%qx(iw,j) > 0.d0) .and. (i>2))    then
            delta1    = (s%uu(i,j) -uu_old(iw ,j)) / s%dsz(i,1)
            delta2    = (s%uu(iw,j)-uu_old(iww,j)) / s%dsz(iw,1)
            wrk1(i,j) = 0.5d0*s%dsu(iw,1)*minmod(delta1,delta2)
          elseif (s%qx(i,j)+s%qx(iw,j) < 0.d0 .and. (i<s%nx-2)) then
            delta1    = (uu_old(i,j) -s%uu(iw ,j)) / s%dsz(i,1)
            delta2    = (uu_old(ie,j)-s%uu(i,j))   / s%dsz(ie,1)
            wrk1(i,j)   = -0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
          endif
        endif
      enddo
    enddo
    
    do j=2,s%ny-1
      do i=2,s%nx-1
        wrk2(i,j) = 0.0_rkind
        js   = min(j+1,s%ny)
        jss  = min(j+2,s%ny)
        jn   = max(j-1,1)
        mindepth = minval(s%zs(i:i+1,jn:jss))-maxval(s%zb(i:i+1,jn:jss))
        if ((mindepth > par%eps)) then
          if     ((s%qy(i+1,j) + s%qy(i,j)   > 0.d0) .and. j>2) then
            delta1    = (s%uu(i,js) -uu_old(i ,j )) / s%dnv(1,j)
            delta2    = (s%uu(i,j)  -uu_old(i ,jn)) / s%dnv(1,jn)
            wrk2(i,j) = 0.5d0*s%dnv(1,j)*minmod(delta1,delta2)
          elseif ( s%qy(i+1,j) + s%qy(i,j)   < 0.d0   .and. j<s%ny-2) then
            delta1    = (uu_old(i,js) -s%uu(i ,j)) / s%dnv(1,j)
            delta2    = (uu_old(i,jss)-s%uu(i,js)) / s%dnv(1,js)
            wrk2(i,j)   = -0.5d0*s%dnv(1,j)*minmod(delta1,delta2)
          endif
        endif          
      enddo
    enddo
    wrk2(:,1   ) = 0.0_rKind    
    wrk2(:,s%ny) = 0.0_rKind    
    
    
    
    !CORRECTION TO U
    do j=2,s%ny 
      do i=2,s%nx-1
        qe = .5*(s%qx(i+1,j  ) + s%qx(i  ,j  ))
        qw = .5*(s%qx(i  ,j  ) + s%qx(i-1,j  ))
        qn = .5*(s%qy(i  ,j-1) + s%qy(i+1,j-1))
        qs = .5*(s%qy(i  ,j  ) + s%qy(i+1,j  ))          
        s%uu(i,j) = s%uu(i,j)-par%dt/s%hum(i,j)*(  (qe*wrk1(i+1,j)-qw*wrk1(i  ,j  ))/s%dsu(i,1) &
                                                +  (qs*wrk2(i  ,j)-qn*wrk2(i  ,j-1))/s%dnz(1,j) )
      enddo
    enddo
     
    !== SECOND ORDER EXPLICIT CORRECTION TO V ==
    if (s%ny>2) then
      do j=2,s%ny
        do i=2,s%nx
          wrk1(i,j) = 0.  
          js   = min(j+1,s%ny)
          jss  = min(j+2,s%ny)
          jn   = max(j-1,1)
          jnn  = max(j-2,1)
          mindepth = minval(s%zs(i,jnn:jss))-maxval(s%zb(i,jnn:jss))
          if (mindepth > par%eps) then
            if   ((s%qy(i,j)+s%qy(i,jn) > par%Umin) .and. (j>2)) then
              delta1    = (s%vv(i,j ) - vv_old(i ,jn )) / s%dnz(1,j)
              delta2    = (s%vv(i,jn) - vv_old(i ,jnn)) / s%dnz(1,jn)
              wrk1(i,j)   = 0.5d0*s%dnv(1,jn)*minmod(delta1,delta2)
            elseif (s%qy(i,j)+s%qy(i,jn) < -par%Umin .and. j<s%ny-2) then
              delta1    = (vv_old(i,j ) - s%vv(i,jn)) / s%dnz(1,j)
              delta2    = (vv_old(i,js) - s%vv(i,j )) / s%dnz(1,js)
              wrk1(i,j)   = -0.5d0*s%dnv(1,j)*minmod(delta1,delta2)
            endif
          endif
        enddo
      enddo
      
       
      do j=2,s%ny-1
        do i=2,s%nx-1
          wrk2(i,j) = 0.0d0
          ie   = min(i+1,s%nx)
          iee  = min(i+2,s%nx)
          iw   = max(i-1,1)
          mindepth = minval(s%zs(iw:iee,j:j+1))-maxval(s%zb(iw:iee,j:j+1))
          if (mindepth > par%eps) then
            if     (s%qx(i,j+1) + s%qx(i,j) > par%Umin .and. i>2) then
              delta1    = (s%vv(ie,j) - vv_old(i ,j )) / s%dsu(i,1)
              delta2    = (s%vv(i ,j) - vv_old(iw,j))  / s%dsu(iw,1)
              wrk2(i,j) = 0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
            elseif (s%qx(i,j+1) + s%qx(i,j) < -par%Umin .and. i<s%nx-2) then
              delta1    = (vv_old(ie,j) -s%vv(i ,j)) / s%dsu(ie,1)
              delta2    = (vv_old(iee,j)-s%vv(ie,j)) / s%dsu(ie,1)
              wrk2(i,j) = -0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
            endif
          endif          
        enddo
      enddo    
      wrk2(1,:)    = 0.0_rKind
      wrk2(s%nx,:) = 0.0_rKind
      
      !CORRECTION TO V
      do j=2,s%ny-1
        do i=2,s%nx
          qe = .5_rKind*(s%qx(i  ,j+1) + s%qx(i  ,j  ))
          qw = .5_rKind*(s%qx(i-1,j  ) + s%qx(i-1,j+1))
          qn = .5_rKind*(s%qy(i  ,j  ) + s%qy(i  ,j-1))
          qs = .5_rKind*(s%qy(i  ,j+1) + s%qy(i  ,j  ))
          s%vv(i,j) = s%vv(i,j)-par%dt/s%hvm(i,j)*(  (qe*wrk2(i,j  )- qw*wrk2(i-1,j))/ s%dsz(i,1)  &
                                                  +  (qs*wrk1(i,j+1)- qn*wrk1(i  ,j))/ s%dnv(1,j)  )
        enddo
      enddo
    endif

end subroutine flow_secondorder_advUV

!
!==============================================================================
subroutine flow_secondorder_advW(s,par,w,w_old)
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
!   Calculates second order correction to the advection terms for W. Only used
!   in combination WITH the non-hydrostatic module

!--------------------------     DEPENDENCIES       ----------------------------
!

! -- MODULES --
  use spaceparams
  use params 
  use xmpi_module, only: halt_program 


!--------------------------     ARGUMENTS          ----------------------------
! 

    type(spacepars)  ,intent(inout)  :: s
    type(parameters) ,intent(in)     :: par
        
    real(kind=rKind),dimension(s%nx+1,s%ny+1),intent(inout) :: w
    real(kind=rKind),dimension(s%nx+1,s%ny+1),intent(in)    :: w_old    

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
    
    real(kind=rKind)           :: mindepth
    real(kind=rKind)           :: delta1
    real(kind=rKind)           :: delta2        

!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------
    !Initialize/allocate arrays on first entry
    if (.not. initialized) then
      call flow_secondorder_init(s)
    endif  

   !== SECOND ORDER EXPLICIT CORRECTION TO W ==
    do j=2,s%ny
      do i=2,s%nx-1
        wrk1(i,j) = 0.0_rKind
        ie   = min(i+1,s%nx)
        iee  = min(i+2,s%nx)
        iw   = max(i-1,1)
        mindepth = minval(s%zs(iw:iee,j))-maxval(s%zb(iw:iee,j))
        if (mindepth > par%eps) then
          if   (s%qx(i,j) > 0.0_rKind  .and. i>2) then
            delta1    = (w(ie,j ) - w_old(i  ,j )) / s%dsu(i,1)
            delta2    = (w(i,j )  - w_old(iw ,j )) / s%dsu(iw,1)
            wrk1(i,j)   = 0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
          elseif (s%qx(i,j) < 0.0_rKind .and. i<s%nx-1) then
            delta1    = (w_old(ie ,j)  - w(i ,j )) / s%dsu(i,1)
            delta2    = (w_old(iee,j ) - w(ie,j )) / s%dsu(ie,1)
            wrk1(i,j)   = -0.5d0*s%dsu(i,1)*minmod(delta1,delta2)
          endif
        endif
      enddo
    enddo
    wrk1(: ,1) = 0.0_rKind
    wrk1(:,s%ny+1) = 0.0_rKind
    wrk1(1,:)  = 0.0_rKind
    wrk1(s%nx,:) = 0.0_rKind
     
    do j=2,s%ny-1
      do i=2,s%nx-1      
        wrk2(i,j) = 0.0_rKind
        js   = min(j+1,s%ny)
        jss  = min(j+2,s%ny)
        jn   = max(j-1,1)
        mindepth = minval(s%zs(i,jn:jss))-maxval(s%zb(i,jn:jss))
        if (mindepth > par%eps) then
          if   (s%qy(i,j) > 0.0_rKind .and. j>2) then
            delta1    = (w(i,js ) - w_old(i  ,j )) / s%dnv(1,j)
            delta2    = (w(i,j  ) - w_old(i  ,jn)) / s%dnv(1,jn)
            wrk2(i,j) = 0.5d0*s%dnv(1,j)*minmod(delta1,delta2)
          elseif (s%qy(i,j) < 0.0_rKind .and. j<s%ny-1) then
            delta1    = (w_old(i,js)   - w(i ,j )) / s%dnv(1,j)
            delta2    = (w_old(i,jss ) - w(i ,js)) / s%dnv(1,js)
            wrk2(i,j) = -0.5d0*s%dnv(1,j)*minmod(delta1,delta2)
          endif
        endif        
      enddo
    enddo    
    wrk2(:,1)      = 0.0_rKind
    wrk2(:,s%ny)   = 0.0_rKind

    !CORRECTION TO W
    do j=2,s%ny
      do i=2,s%nx
        w(i,j) = w(i,j)-par%dt/s%hh(i,j)*(  (s%qx(i,j)*wrk1(i,j)- s%qx(i-1,j)*wrk1(i-1,j))/ s%dsz(i,1)  &
                                         +  (s%qy(i,j)*wrk2(i,j)- s%qy(i,j-1)*wrk2(i,j-1))/ s%dnz(1,j)  )
      enddo
    enddo

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

    real(kind=rKind)           :: mindepth    
    real(kind=rKind)           :: delta1
    real(kind=rKind)           :: delta2    
!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------

    !== SECOND ORDER EXPLICIT CORRECTION TO CONTINUITY ==
    !Initialize/allocate arrays on first entry
      if (.not. initialized) then
        call flow_secondorder_init(s)
      endif  
      !return
      !correction to mass flux qx
      do j=2,s%ny
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
      do j=2,s%ny
        do i=2,s%nx
          s%zs(i,j) = s%zs(i,j)-par%dt*(  (wrk1(i,j)-wrk1(i-1,j))/ s%dsz(i,1)  &
                                       +  (wrk2(i,j)-wrk2(i,j-1))/ s%dnz(1,j)  )
        enddo
      enddo
      
      !Update fluxes
      s%qx(2:s%nx,2:s%ny) = s%qx(2:s%nx,2:s%ny) +wrk1(2:s%nx,2:s%ny)
      s%qy(2:s%nx,2:s%ny) = s%qy(2:s%nx,2:s%ny) +wrk2(2:s%nx,2:s%ny)

end subroutine flow_secondorder_con


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
