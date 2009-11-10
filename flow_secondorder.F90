!==============================================================================
!                               FLOW_SECONDORDER_MODULE                        
!==============================================================================



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

subroutine flow_secondorder_init(s,iUnit)
!

!-------------------------------------------------------------------------------
!                             DECLARATIONS
!-------------------------------------------------------------------------------

!
!--------------------------        PURPOSE         ----------------------------
!
!   Initializes the resources needed for the second order mcCormack scheme
!

!--------------------------     DEPENDENCIES       ----------------------------
!

! -- MODULES --
  use spaceparams
  use params 
  use xmpi_module, only: Halt_Program 

!
!--------------------------     ARGUMENTS          ----------------------------
!
  type(spacepars)    ,intent(inout) :: s
  integer(kind=iKind),intent(in)    :: iUnit    !Unit number for error messages
                                                ! iUnit>0  : write to indicated unit
                                                ! iUnit<=0 : write to screen
!         
!--------------------------     LOCAL VARIABLES    ----------------------------
! 
  integer(kind=iKind)              :: iAllocErr   ! Error code returned (if /=0)
  integer(kind=iKind)              :: iWriteUnit  ! Write unit (see also iUnit)
!
!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------
!
  

  
 
  
   if (iAllocErr /= 0) call Halt_program
  
  !Allocate resources
  allocate (  wrk1(s%nx+1,s%ny+1),stat=iAllocErr);   if (iAllocErr /= 0) goto 10
  allocate (  wrk2(s%nx+1,s%ny+1),stat=iAllocErr);   if (iAllocErr /= 0) goto 10
    
  wrk1   = 0.0_rKind
  wrk2   = 0.0_rKind
  
10 if (iAllocErr /=0) then
    if   (iUnit > 0 ) then
      iWriteUnit = iUnit
    else 
      iWriteUnit = iScreen
    endif
    write(iWriteUnit,'(a,i8)') 'Could not allocate array in flow_secondorder_init, error code returned: ', iAllocErr
    call Halt_program
  endif
  initialized = .true.
end subroutine flow_secondorder_init  

!
!==============================================================================    
subroutine flow_secondorder_advUV(s,par,uu_old,vv_old)
!==============================================================================    
!

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
      call flow_secondorder_init(s,0)
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
            delta1    = (s%uu(i,j) -uu_old(iw ,j)) / (s%xu(i) -s%xu(iw))
            delta2    = (s%uu(iw,j)-uu_old(iww,j)) / (s%xu(iw)-s%xu(iww))
            wrk1(i,j) = (s%xz(i)-s%xu(iw))  *minmod(delta1,delta2)
          elseif (s%qx(i,j)+s%qx(iw,j) < 0.d0 .and. (i<s%nx)) then
            delta1    = (uu_old(i,j) -s%uu(iw ,j)) / (s%xu(i) -s%xu(iw))
            delta2    = (uu_old(ie,j)-s%uu(i,j))   / (s%xu(ie)-s%xu(i))
            wrk1(i,j)   = -(s%xu(i)-s%xz(i))*minmod(delta1,delta2)
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
            delta1    = (s%uu(i,js) -uu_old(i ,j )) / (s%yz(js)-s%yz(j ))
            delta2    = (s%uu(i,j)  -uu_old(i ,jn)) / (s%yz(j )-s%yz(jn))
            wrk2(i,j) = (s%yv(j)-s%yz(j))*minmod(delta1,delta2)
          elseif ( s%qy(i+1,j) + s%qy(i,j)   < 0.d0   .and. j<s%ny) then
            delta1    = (uu_old(i,js) -s%uu(i ,j)) / (s%yz(js )-s%yz(j ))
            delta2    = (uu_old(i,jss)-s%uu(i,js)) / (s%yz(jss)-s%yz(js))
            wrk2(i,j)   = -(s%yz(js)-s%yv(j))*minmod(delta1,delta2)
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
        s%uu(i,j) = s%uu(i,j)-par%dt/s%hum(i,j)*(  (qe*wrk1(i+1,j)-qw*wrk1(i  ,j  ))/(s%xz(i+1)-s%xz(i  )) &
                                                +  (qs*wrk2(i  ,j)-qn*wrk2(i  ,j-1))/(s%yv(j  )-s%yv(j-1)) )
      enddo
    enddo
     
    !== SECOND ORDER EXPLICIT CORRECTION TO V ==
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
            delta1    = (s%vv(i,j ) - vv_old(i ,jn )) / (s%yv(j )-s%yv(jn ))
            delta2    = (s%vv(i,jn) - vv_old(i ,jnn)) / (s%yv(jn)-s%yv(jnn))
            wrk1(i,j)   = (s%yz(j)-s%yv(jn))*minmod(delta1,delta2)
          elseif (s%qy(i,j)+s%qy(i,jn) < -par%Umin .and. j<s%ny) then
            delta1    = (vv_old(i,j ) - s%vv(i,jn)) / (s%yv(j) -s%yv(jn))
            delta2    = (vv_old(i,js) - s%vv(i,j )) / (s%yv(js)-s%yv(j ))
            wrk1(i,j)   = -(s%yv(j)-s%yz(j))*minmod(delta1,delta2)
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
            delta1    = (s%vv(ie,j) - vv_old(i ,j )) / (s%xz(ie)-s%xz(i ))
            delta2    = (s%vv(i ,j) - vv_old(iw,j))  / (s%xz(i )-s%xz(iw))
            wrk2(i,j) = (s%xu(i)-s%xz(i))*minmod(delta1,delta2)
          elseif (s%qx(i,j+1) + s%qx(i,j) < -par%Umin .and. i<s%nx) then
            delta1    = (vv_old(ie,j) -s%vv(i ,j)) / (s%xz(ie )-s%xz(i ))
            delta2    = (vv_old(iee,j)-s%vv(ie,j)) / (s%xz(iee)-s%xz(ie))
            wrk2(i,j) = -(s%xz(ie)-s%xu(i))*minmod(delta1,delta2)
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
        s%vv(i,j) = s%vv(i,j)-par%dt/s%hvm(i,j)*(  (qe*wrk2(i,j  )- qw*wrk2(i-1,j))/(s%xu(i  )-s%xu(i-1))  &
                                                +  (qs*wrk1(i,j+1)- qn*wrk1(i  ,j))/(s%yz(j+1)-s%yz(j  ))  )
      enddo
    enddo
    
!        if (par%nonh == 1) then
!      !Include explicit approximation for pressure in s%uu and s%vv
!    
!      do j=2,s%ny
!        do i=2,s%nx-1
!          if (s%wetu(i,j) == 1) then
!            s%uu(i,j) = s%uu(i,j) - 0.5_rKind*par%dt/s%hum(i,j) * ( (s%zs(i+1,j)-s%zb(i  ,j)) * s%pres(i+1,j)   &
!                                                                  - (s%zs(i  ,j)-s%zb(i+1,j)) * s%pres(i  ,j)  ) &
!                                                                  / (s%xz(i+1)-s%xz(i))
!          endif
!        enddo
!      enddo
!
!      do j=2,s%ny-1
!        do i=2,s%nx
!          if (s%wetv(i,j) == 1) then        
!            s%vv(i,j) = s%vv(i,j) - 0.5_rKind*par%dt/s%hvm(i,j) * ( (s%zs(i,j+1)-s%zb(i,j  )) * s%pres(i,j+1)   &
!                                                                  - (s%zs(i,j  )-s%zb(i,j+1)) * s%pres(i,j  ) ) &
!                                                                  / (s%yz(j+1)-s%yz(j))
!          endif
!        enddo
!      enddo
!    endif
end subroutine flow_secondorder_advUV

!
!==============================================================================
subroutine flow_secondorder_advW(s,par,w,w_old)
!==============================================================================
!

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
  use xmpi_module, only: Halt_Program 


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
      call flow_secondorder_init(s,0)
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
            delta1    = (w(ie,j ) - w_old(i  ,j )) / (s%xz(ie )-s%xz(i ))
            delta2    = (w(i,j )  - w_old(iw ,j )) / (s%xz(i)  -s%xz(iw))
            wrk1(i,j)   = (s%xu(i)-s%xz(i))*minmod(delta1,delta2)
          elseif (s%qx(i,j) < 0.0_rKind .and. i<s%nx-1) then
            delta1    = (w_old(ie ,j)  - w(i ,j )) / (s%xz(ie ) -s%xz(i ))
            delta2    = (w_old(iee,j ) - w(ie,j )) / (s%xz(iee) -s%xz(ie))
            wrk1(i,j)   = -(s%xz(ie)-s%xu(i))*minmod(delta1,delta2)
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
            delta1    = (w(i,js ) - w_old(i  ,j )) / (s%yz(js)-s%yz(j ))
            delta2    = (w(i,j  ) - w_old(i  ,jn)) / (s%yz(j )-s%yz(jn))
            wrk2(i,j) = (s%yv(j)-s%yz(j))*minmod(delta1,delta2)
          elseif (s%qy(i,j) < 0.0_rKind .and. j<s%ny-1) then
            delta1    = (w_old(i,js)   - w(i ,j )) / (s%yz(js ) -s%yz(j ))
            delta2    = (w_old(i,jss ) - w(i ,js)) / (s%yz(jss) -s%yz(js))
            wrk2(i,j) = -(s%yz(js)-s%yv(j))*minmod(delta1,delta2)
          endif
        endif        
      enddo
    enddo    
    wrk2(:,1)      = 0.0_rKind
    wrk2(:,s%ny)   = 0.0_rKind
    
   ! write(*,*) maxval(abs(wrk1))
 
    !CORRECTION TO W
    do j=2,s%ny
      do i=2,s%nx
        w(i,j) = w(i,j)-par%dt/s%hh(i,j)*(  (s%qx(i,j)*wrk1(i,j)- s%qx(i-1,j)*wrk1(i-1,j))/(s%xu(i  )-s%xu(i-1))  &
                                         +  (s%qy(i,j)*wrk2(i,j)- s%qy(i,j-1)*wrk2(i,j-1))/(s%yv(j)  -s%yv(j-1))  )
      enddo
    enddo

end subroutine flow_secondorder_advW

!
!==============================================================================    
subroutine flow_secondorder_con(s,par,zs_old)
!==============================================================================    
!

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
        call flow_secondorder_init(s,0)
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
              delta1    =  (s%zs(ie,j)- zs_old(i,j ))/(s%xz(i+1)-s%xz(i))
              delta2    =  (s%zs(i,j) - zs_old(iw,j))/(s%xz(i)  -s%xz(i-1))
              wrk1(i,j)  =   s%uu(i,j)*(s%xu(i)-s%xz(i))*minmod(delta1,delta2)
            elseif (s%uu(i,j) < -par%umin .and. i<s%nx-1) then
              delta1    =  (zs_old(iee,j) - s%zs(ie,j)) / (s%xz(i+2)-s%xz(i+1))
              delta2    =  (zs_old(ie ,j) - s%zs(i ,j)) / (s%xz(i+1)-s%xz(i))
              wrk1(i,j)  = - s%uu(i,j)*(s%xz(ie)-s%xu(i))*minmod(delta1,delta2)
            endif          
          endif
        enddo
      enddo
      wrk1(1   ,:) = 0.0_rkind
      wrk1(s%nx,:) = 0.0_rkind

      
      
!      do j=2,s%ny
!        do i=2,s%nx-1
!          wrk1(i,j) = 0.0_rKind
!          ie  = min(i+1,s%nx)
!          iee = min(i+2,s%nx)
!          iw  = max(i-1,2)
!          mindepth = minval(s%zs(iw:iee,j))-maxval(s%zb(iw:iee,j))
!          if (mindepth>par%eps) then
!            if     (s%uu(i,j) >  par%Umin .and. i>2     ) then
!              delta1    =  (s%zs(ie,j)- zs_old(i,j ))/(s%xz(i+1)-s%xz(i))
!              delta2    =  (s%zs(i,j) - zs_old(iw,j))/(s%xz(i)  -s%xz(i-1))
!              wrk1(i,j)  = s%uu(i,j)*(.5d0*(s%zs(i,j)+(s%xu(i)-s%xz(i))*minmod(delta1,delta2)-s%zb(i,j))-.5d0*s%hu(i,j))
!            elseif (s%uu(i,j) < -par%Umin .and. i<s%nx-1) then
!              delta1    =  (zs_old(iee,j) - s%zs(ie,j)) / (s%xz(i+2)-s%xz(i+1))
!              delta2    =  (zs_old(ie ,j) - s%zs(i ,j)) / (s%xz(i+1)-s%xz(i))
!              wrk1(i,j)  = s%uu(i,j)* (.5d0*(s%zs(ie,j)-(s%xz(ie)-s%xu(i))*minmod(delta1,delta2)-s%zb(ie,j))-.5d0*s%hu(i,j))
!            endif          
!          endif
!        enddo
!      enddo
!      wrk1(1   ,:) = 0.0_rKind
!      wrk1(s%nx,:) = 0.0_rKind


      !Correction to mass flux qy
      do j=2,s%ny-1
        do i=2,s%nx
          wrk2(i,j) = 0.0_rKind
          js  = min(j+1,s%ny)
          jss = min(j+2,s%ny)
          jn  = max(j-1,2)
          mindepth = minval(s%zs(i,jn:jss))-maxval(s%zb(i,jn:jss))
          if (mindepth> par%eps) then
            if     (s%vv(i,j) >  par%Umin .and. j>2) then
              delta1    = (s%zs(i,js) - zs_old(i,j ))/(s%yz(j+1)-s%yz(j))
              delta2    = (s%zs(i,j)  - zs_old(i,jn))/(s%yz(j)  -s%yz(j-1))
              wrk2(i,j) =  s%vv(i,j)*(s%yv(j)-s%yz(j))*minmod(delta1,delta2)
            elseif (s%vv(i,j) < -par%Umin .and. j<s%ny-1) then
              delta1    = (zs_old(i,jss) - s%zs(i,js)) / (s%yz(j+2)-s%yz(j+1))
              delta2    = (zs_old(i ,js) - s%zs(i ,j)) / (s%yz(j+1)-s%yz(j))
              wrk2(i,j) =  s%vv(i,j)*(s%yv(j)-s%yz(j))*minmod(delta1,delta2)
            endif          
          endif
        enddo
      enddo
      wrk2(:,1   ) = 0.0_rKind
      wrk2(:,s%ny) = 0.0_rKind

      !Update waterlevels
      do j=2,s%ny
        do i=2,s%nx
          s%zs(i,j) = s%zs(i,j)-par%dt*(  (wrk1(i,j)-wrk1(i-1,j))/(s%xu(i)-s%xu(i-1))  &
                                       +  (wrk2(i,j)-wrk2(i,j-1))/(s%yv(j)-s%yv(j-1))  )
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
! The minmod slope limiter essentially essentially corrects the first order upwind 
! estimate for hu/uz/yz with either a second order upwind or central approximation.
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
