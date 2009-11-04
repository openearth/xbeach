!==============================================================================
!                               MODULE NH_MAIN                        
!==============================================================================

module nonh_module


  implicit none
  
#ifndef USEMPI

! If mpi is defined, the non-hydrostatic module is NOT included in the compilation
! to avoid unwanted side effects.  

!******************************************************************************
!                                 INTERFACE                                
!******************************************************************************
    
  save                            
  private
  
  !----------------------------- PARAMETERS -----------------------------------
  
  include 'nh_pars.inc'               !Default precision etc.

  !----------------------------- VARIABLES  -----------------------------------
 


  !--- PUBLIC VARIABLES ---
  character(len=80)            :: PrintFileName = 'Print'
  integer(kind=iKind),public   :: PrintFileUnit = iPrint
  
  logical                      :: initialized   = .false.
  
  !        NONE

  real (kind=rKind), allocatable, dimension(:,:,:)   :: au     !Pressure coeficients for u(i,j)
  real (kind=rKind), allocatable, dimension(:,:,:)   :: av     !Pressure coeficients for v(i,j)
  real (kind=rKind), allocatable, dimension(:,:,:,:) :: aw     !Pressure coeficients for w(k,i,j)
  
  real (kind=rKind), allocatable, dimension(:,:)     :: aur    !RHS for u(i,j)
  real (kind=rKind), allocatable, dimension(:,:)     :: avr    !RHS for v(i,j)
  real (kind=rKind), allocatable, dimension(:,:,:)   :: awr    !RHS for w(k,i,j)
  real (kind=rKind), allocatable, dimension(:,:,:,:) :: dpdz   !RHS for w(k,i,j)
  real (kind=rKind), allocatable, dimension(:,:,:)   :: dpdzr   !RHS for w(k,i,j)  
    
  real (kind=rKind), allocatable, dimension(:,:,:)   :: mat   !Pressure matrix
  real (kind=rKind), allocatable, dimension(:,:)     :: rhs   !RHS of the pressure matrix
  
  real (kind=rKind), allocatable, dimension(:,:)     :: zsu    !free surface/bottom in velocity points
  real (kind=rKind), allocatable, dimension(:,:)     :: zsv    !free surface/bottom in velocity points
  
  real (kind=rKind), allocatable, dimension(:,:)     :: zbu    !free surface/bottom in velocity points
  real (kind=rKind), allocatable, dimension(:,:)     :: zbv    !free surface/bottom in velocity points
  
  real (kind=rKind), allocatable, dimension(:,:)     :: dp    !free surface/bottom in velocity points  
    
  real (kind=rKind), allocatable, dimension(:)       ::  dxz    !free surface/bottom in velocity points
  real (kind=rKind), allocatable, dimension(:)       ::  dyz    !free surface/bottom in velocity points  
  real (kind=rKind), allocatable, dimension(:)       ::  dxu    !free surface/bottom in velocity points
  real (kind=rKind), allocatable, dimension(:)       ::  dyv    !free surface/bottom in velocity points  
  real (kind=rKind), allocatable, dimension(:)       ::  ddxz    !free surface/bottom in velocity points
  real (kind=rKind), allocatable, dimension(:)       ::  ddyz    !free surface/bottom in velocity points  
  real (kind=rKind), allocatable, dimension(:)       ::  ddxu    !free surface/bottom in velocity points
  real (kind=rKind), allocatable, dimension(:)       ::  ddyv    !free surface/bottom in velocity points  
  
  real (kind=rKind), allocatable, dimension(:,:)     :: ws_old    !Velocity at the old time level  
  real (kind=rKind), allocatable, dimension(:,:)     :: wb_old    !Velocity at the old time level    

  !--- PUBLIC SUBROUTINES ---
  public nonh_init
  public nonh_cor
  public nonh_free
  
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
  subroutine nonh_init(s)
!==============================================================================    
!

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

!--------------------------     ARGUMENTS          ----------------------------
!
    type(spacepars) ,intent(inout) :: s
!

!--------------------------     LOCAL VARIABLES    ----------------------------
 
    integer(kind=iKind)        :: iw,i,ie,jn,j,js,iww,iee
!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------
   
    open(unit=PrintFileUnit,file=trim(PrintFileName))
        
    allocate(au (0:1,s%nx+1,s%ny+1)); au = 0.0_rKind
    allocate(av (0:1,s%nx+1,s%ny+1)); av = 0.0_rKind
    
    allocate(aw (5,0:1,s%nx+1,s%ny+1)); aw = 0.0_rKind
    
    allocate(aur(   s%nx+1,s%ny+1)); aur = 0.0_rKind
    allocate(avr(   s%nx+1,s%ny+1)); avr = 0.0_rKind
    allocate(awr(0:1,s%nx+1,s%ny+1)); awr = 0.0_rKind
    
    allocate(mat(5,s%nx+1,s%ny+1));  mat = 0.0_rKind
    allocate(rhs(  s%nx+1,s%ny+1));  rhs = 0.0_rKind
    
    allocate(zbu(  s%nx+1,s%ny+1));  zbu = 0.0_rKind
    allocate(zbv(  s%nx+1,s%ny+1));  zbv = 0.0_rKind
    allocate(zsu(  s%nx+1,s%ny+1));  zsu = 0.0_rKind
    allocate(zsv(  s%nx+1,s%ny+1));  zsv = 0.0_rKind
        
    allocate(dxz (  s%nx+1));  dxz  = 0.0_rKind
    allocate(dyz (  s%ny+1));  dyz  = 0.0_rKind
    allocate(dxu (  s%nx+1));  dxu  = 0.0_rKind
    allocate(dyv (  s%ny+1));  dyv  = 0.0_rKind
    allocate(ddxz (  s%nx+1));  ddxz  = 0.0_rKind
    allocate(ddyz (  s%ny+1));  ddyz  = 0.0_rKind
    allocate(ddxu (  s%nx+1));  ddxu  = 0.0_rKind
    allocate(ddyv (  s%ny+1));  ddyv  = 0.0_rKind      
    
    allocate(ws_old(s%nx+1,s%ny+1)); ws_old = 0.0_rKind
    allocate(wb_old(s%nx+1,s%ny+1)); wb_old = 0.0_rKind      
    allocate(dp(s%nx+1,s%ny+1)); dp = 0.      
    allocate(dpdz(5,0:1,s%nx+1,s%ny+1)); dpdz = 0.
    allocate(dpdzr(0:1,s%nx+1,s%ny+1)); dpdzr = 0.        
    
    if (.not. associated(s%pres)) then   
      allocate(s%pres(s%nx+1,s%ny+1)); s%pres = 0.0_rKind   
      allocate(s%ws(s%nx+1,s%ny+1)); s%ws = 0.0_rKind
      allocate(s%wb(s%nx+1,s%ny+1)); s%wb = 0.0_rKind
    endif  
       
    do j=1,s%ny+1
      js = min(j+1,s%ny+1)
      jn = max(j-1,1     )      
      dyz(j) = s%yv(j) -s%yv(jn)
      dyv(j) = s%yz(js)-s%yz(j)      
    enddo    
    dyv(s%ny+1)   = dyv(s%ny)
    dyz(1)      = dyz(2)
    ddyz        = 1.0_rKind/dyz
    ddyv        = 1.0_rKind/dyv
    
    do i=1,s%nx+1
      ie = min(i+1,s%nx+1)
      iw = max(i-1,1     )      
      dxz(i) = s%xu(i) -s%xu(iw)
      dxu(i) = s%xz(ie)-s%xz(i)
    enddo    
    dxu(s%nx+1)   = dxu(s%nx)
    dxz(1)      = dxz(2)
    ddxu        = 1.0_rKind/dxu
    ddxz        = 1.0_rKind/dxz
    
    mat(1,1,:)      = 1.0_rKind
    mat(1,s%nx+1,:) = 1.0_rKind
    mat(1,:,1)      = 1.0_rKind
    mat(1,:,s%ny+1) = 1.0_rKind
    
    initialized   = .true. 
  end subroutine nonh_init

!
!==============================================================================
  subroutine nonh_cor(s,par)
!==============================================================================
!

!-------------------------------------------------------------------------------
!                             DECLARATIONS
!-------------------------------------------------------------------------------

!--------------------------        PURPOSE         ----------------------------

!  Releases resources

!--------------------------     DEPENDENCIES       ----------------------------  
    use spaceparams
    use params
    !use timer_module
    use solver_module
    use flow_secondorder_module
!--------------------------     ARGUMENTS          ----------------------------

    type(spacepars) ,intent(inout)                       :: s
    type(parameters),intent(in)                          :: par 

!--------------------------     LOCAL VARIABLES    ----------------------------
 
    !Indices
    integer(kind=iKind)                     :: i,ie,iee,iw,iww               !Index variables
    integer(kind=iKind)                     :: j               !Index variables
    
    !Work
    real(kind=rKind)                        :: dwdx1,dwdx2
    real(kind=rKind)                        :: dwdy1,dwdy2    
    real(kind=rKind)                        :: dzdx,dzdy,delta1,delta2
    real(kind=rKind)                        :: dzsu,dzsv    
    real(kind=rKind)                        :: dzbu,dzbv
    real(kind=rKind)                        :: vol    
    real(kind=rKind)                        :: mindepth
   

!
!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------   
!
  !call timer_start(timer_flow_nonh)
  
  if (.not. initialized) then
    call nonh_init(s)
  endif  
  
  
  !Calculate explicit part vertical momentum (advection)
  do j=2,s%ny
    do i=2,s%nx
      s%ws(i,j) = s%ws(i,j) -   par%dt*( ddxu(i-1)*max(s%qx(i-1,j  ),0.0_rKind)*(ws_old(i  ,j  )-ws_old(i-1,j  ))/s%hh(i,j)   &
                                       + ddxu(i)  *min(s%qx(i  ,j  ),0.0_rKind)*(ws_old(i+1,j  )-ws_old(i  ,j  ))/s%hh(i,j)   &
                                       + ddyv(j-1)*max(s%qy(i  ,j-1),0.0_rKind)*(ws_old(i  ,j  )-ws_old(i  ,j-1))/s%hh(i,j)   &
                                       + ddyv(j  )*min(s%qy(i  ,j  ),0.0_rKind)*(ws_old(i  ,j+1)-ws_old(i  ,j  ))/s%hh(i,j) )
    enddo
  enddo

  !Calculate explicit part vertical viscocity
  do j=2,s%ny
    do i=2,s%nx
      dwdx1 = (ws_old(i  ,j  )-ws_old(i-1,j  ))*ddxu(i-1)
      dwdx2 = (ws_old(i+1,j  )-ws_old(i  ,j  ))*ddxu(i)
      dwdy1 = (ws_old(i  ,j  )-ws_old(i  ,j-1))*ddyv(j-1)
      dwdy2 = (ws_old(i  ,j+1)-ws_old(i  ,j  ))*ddyv(j)
      s%ws(i,j) = s%ws(i,j)   + par%dt*par%nuh*(dwdx2-dwdx1)*ddxz(i)*real(s%wetu(i,j)*s%wetu(i-1,j),rKind) &
                              + par%dt*par%nuh*(dwdy2-dwdy1)*ddyz(j)*real(s%wetv(i,j)*s%wetv(i,j-1),rKind)
    enddo
  enddo



  !Calculate explicit part vertical momentum (advection)
  do j=2,s%ny
    do i=2,s%nx
      s%wb(i,j) = s%wb(i,j) -   par%dt*( ddxu(i-1)*max(s%qx(i-1,j  ),0.0_rKind)*(s%wb(i  ,j  )-s%wb(i-1,j  ))/s%hh(i,j) &
                                       + ddxu(i)  *min(s%qx(i  ,j  ),0.0_rKind)*(s%wb(i+1,j  )-s%wb(i  ,j  ))/s%hh(i,j) &
                                       + ddyv(j-1)*max(s%qy(i  ,j-1),0.0_rKind)*(s%wb(i  ,j  )-s%wb(i  ,j-1))/s%hh(i,j) &
                                       + ddyv(j  )*min(s%qy(i  ,j  ),0.0_rKind)*(s%wb(i  ,j+1)-s%wb(i  ,j  ))/s%hh(i,j) )                                      
    enddo
  enddo
    
 !Calculate explicit part vertical viscocity
  do j=2,s%ny
    do i=2,s%nx
      dwdx1 = (wb_old(i  ,j  )-wb_old(i-1,j  ))*ddxu(i-1)
      dwdx2 = (wb_old(i+1,j  )-wb_old(i  ,j  ))*ddxu(i)
      dwdy1 = (wb_old(i  ,j  )-wb_old(i  ,j-1))*ddyv(j-1)
      dwdy2 = (wb_old(i  ,j+1)-wb_old(i  ,j  ))*ddyv(j)
      s%wb(i,j) = s%wb(i,j)   + par%dt*par%nuh*(dwdx2-dwdx1)*ddxz(i)*real(s%wetu(i,j)*s%wetu(i-1,j),rKind) &
                              + par%dt*par%nuh*(dwdy2-dwdy1)*ddyz(j)*real(s%wetv(i,j)*s%wetv(i,j-1),rKind)
    enddo
  enddo
 
 
 
  !Free surface location in u-point
  do j=1,s%ny+1
    do i=1,s%nx
      if     (s%uu(i,j) > 0.0_rKind) then
        zsu(i,j) = s%zs(i  ,j)      
      elseif (s%uu(i,j) < 0.0_rKind) then
        zsu(i,j) = s%zs(i+1,j)
      else
        zsu(i,j) = max(s%zs(i  ,j),s%zs(i+1,j))
      endif
    enddo
  enddo  

  !Free surface location in v-point
  do j=1,s%ny
    do i=1,s%nx+1
      if     (s%vv(i,j) > 0.0_rKind) then
        zsv(i,j) = s%zs(i  ,j)      
      elseif (s%vv(i,j) < 0.0_rKind) then
        zsv(i,j) = s%zs(i,j+1)
      else
        zsv(i,j) = max(s%zs(i  ,j),s%zs(i,j+1))
      endif
    enddo
  enddo

  !Bottom location in U point
  do j=1,s%ny+1
    do i=1,s%nx
      zbu(i,j) = zsu(i,j)-s%hu(i,j)
    enddo
  enddo
  zbu(s%nx+1,:) = s%zb(s%nx+1,:)

  !Bottom location in V point  
  do j=1,s%ny
    do i=1,s%nx+1
      zbv(i,j) = zsv(i,j)-s%hv(i,j)
    enddo
  enddo  
  zbv(:,s%ny+1) = s%zb(:,s%ny+1)

  !call timer_start(timer_flow_nonh_mat)
  
  !Built pressure coefficients U  
  !call timer_start(timer_flow_nonh_au)  
  do j=2,s%ny
    do i=2,s%nx
      if (s%wetU(i,j)==1) then
        vol       = 0.5_rKind*par%dt/(s%hum(i,j)*dxu(i))      
        au(1,i,j) = - (s%zs(i+1,j) - s%zb(i  ,j))*vol
        au(0,i,j) = + (s%zs(i  ,j) - s%zb(i+1,j))*vol
        aur(i,j)  = s%uu(i,j)
      else
        au(1,i,j) =  0.0_rKind
        au(0,i,j) =  0.0_rKind
        aur(i,j)  =  0.0_rKind
      endif  
    enddo    
  enddo
  au(:,1,:)      = 0.0_rKind
  au(:,s%nx,:)   = 0.0_rKind
  aur(1,:)       = s%uu(1,:)
  aur(s%nx,:)    = s%uu(s%nx,:)
  !call timer_stop(timer_flow_nonh_au)
  
  !Built pressure coefficients V
  !call timer_start(timer_flow_nonh_av)    
  do j=2,s%ny
    do i=2,s%nx
      vol       = 0.5_rKind*par%dt/(s%hvm(i,j)*dyv(j))
      if (s%wetV(i,j)==1)then
        av(1,i,j)  = -(s%zs(i  ,j+1) - s%zb(i  ,j  ))*vol
        av(0,i,j)  = +(s%zs(i  ,j  ) - s%zb(i  ,j+1))*vol
        avr(i,j)   = s%vv(i,j)
      else
        av(1,i,j) =  0.0_rKind
        av(0,i,j) =  0.0_rKind
        avr(i,j)  =  0.0_rKind
      endif  
    enddo    
  enddo
  av(:,:,1)     = 0.0_rKind
  av(:,:,s%ny)  = 0.0_rKind
  avr(:,1)      = s%vv(:,1)
  avr(:,s%ny)   = s%vv(:,s%ny)
  !call timer_stop(timer_flow_nonh_av)

  !Built pressure coefficients for W
  !call timer_start(timer_flow_nonh_aw)    
  
  do j=2,s%ny
    do i=2,s%nx
        !Bottom gradients
        dzdx = .5_rKind*(  zbu(i,  j) - zbu(i-1,j) )  *ddxz(i)
        dzdy = .5_rKind*(  zbv(i,  j) - zbv(i  ,j-1) )*ddyz(j)
     
        dpdz(1,0,i,j) = - dzdx *( au(0,i  ,j  )+au(1,i-1,j  ) ) &           !main diagonal
                        - dzdy *( av(0,i  ,j  )+av(1,i  ,j-1) )
        dpdz(2,0,i,j) = - dzdx *  au(0,i-1,j  )                             !west
        dpdz(3,0,i,j) = - dzdx *  au(1,i  ,j  )                             !east
        dpdz(4,0,i,j) = - dzdy *  av(0,i  ,j-1)                             !south
        dpdz(5,0,i,j) = - dzdy *  av(1,i  ,j  )                             !north

        dpdzr(0,i,j)  = - dzdx*( s%uu(i,j)+s%uu(i-1,j  ) ) &
                        - dzdy*( s%vv(i,j)+s%vv(i  ,j-1) ) &
                        + s%wb(i,j)
                      
        dpdz(1,1,i,j)   = - 2.0_rKind*par%dt/s%hh(i,j) - dpdz(1,0,i,j)
        dpdz(2,1,i,j)   =                              - dpdz(2,0,i,j)
        dpdz(3,1,i,j)   =                              - dpdz(3,0,i,j)
        dpdz(4,1,i,j)   =                              - dpdz(4,0,i,j)
        dpdz(5,1,i,j)   =                              - dpdz(5,0,i,j)
        dpdzr(1,i,j)    =  - dpdzr(0,i,j)
    enddo
  enddo  
   

  if (par%secorder == 1) then
    !Include non-hydrostatic pressure explicitly
    do j=2,s%ny
      do i=2,s%nx
        s%ws(i,j) = s%ws(i,j)  - dpdzr(1,i,j)+ 2.0_rKind*par%dt/s%hh(i,j)*s%pres(i,j)
      enddo
    enddo
  else  
    do j=2,s%ny
      do i=2,s%nx
        s%ws(i,j) = s%ws(i,j)  - dpdzr(1,i,j)
      enddo
    enddo  
  endif  

  do j=2,s%ny
    do i=2,s%nx
      s%wb(i,j) = s%wb(i,j)  - dpdzr(0,i,j)      
    enddo
  enddo

    
  if (par%secorder == 1) then    
    call flow_secondorder_advW(s,par,s%ws,ws_old)
    call flow_secondorder_advW(s,par,s%wb,wb_old)    
  endif   
   
  do j=2,s%ny
    do i=2,s%nx
      if (s%wetZ(i,j)==1) then     
        aw(:,:,i,j) = - dpdz(:,:,i,j)
      else
        aw(:,:,i,j)  = 0.0_rKind
      endif        
    enddo
  enddo
  !call timer_stop(timer_flow_nonh_aw)  

  do j=2,s%ny
    do i=2,s%nx
      if (s%wetZ(i,j)==1) then
        awr(0,i,j) = s%wb(i,j) !- dpdzr(0,i,j)
        awr(1,i,j) = s%ws(i,j) !- dpdzr(1,i,j)
      else
        awr(:,i,j) = 0.0_rKind
      endif
    enddo
  enddo
 
  !Substitute in the continuity equation
  !call timer_start(timer_flow_nonh_subs)
  do j=2,s%ny
    do i=2,s%nx
      if (s%wetZ(i,j)==1) then
        dzsu = .5*dyz(j)*(zsu(i,j)-zsu(i-1,j))
        dzsv = .5*dxz(i)*(zsv(i,j)-zsv(i,j-1))
        dzbu = .5*dyz(j)*(zbu(i,j)-zbu(i-1,j))
        dzbv = .5*dxz(i)*(zbv(i,j)-zbv(i,j-1))
    
        mat(1,i,j) =    dyz(j)*( s%hu(i,j)*au(0,i,j) - s%hu(i-1,j  )*au(1,i-1,  j) )  & !subs U left/right face
                   +    dxz(i)*( s%hv(i,j)*av(0,i,j) - s%hv(i  ,j-1)*av(1,i  ,j-1) )  & !subs V left/rigth face
                   -    dzsu*(au(0,i,j) + au(1,i-1,j)) - dzsv*(av(0,i,j) + av(1,i-1,j)) & !kin. boun. top
                   +    dzbu*(au(0,i,j) + au(1,i-1,j)) + dzbv*(av(0,i,j) + av(1,i-1,j)) & !kin. boun. bot V
                   +    dxz(i)*dyz(j)*(aw(1,1,i,j)-aw(1,0,i,j))                       !Vert velocity top/bottom face

        mat(2,i,j) = -  dyz(j)*s%hu(i-1,j)*au(0,i-1,  j)                              &
                    -   dzsu*au(0,i-1,j) +    dzbu*au(0,i-1,j)                          & !kin. boun. top/bot U
                    +   dxz(i)*dyz(j)*(aw(2,1,i,j)-aw(2,0,i,j))                       !Vert velocity top/bottom face
                  
        mat(3,i,j) =    dyz(j)*s%hu(i,j)  *au(1,i,  j)                                &
                   -    dzsu*au(1,i  ,j) +  dzbu*au(1,i,j)                              & !kin. boun. top/bot U
                   +    dxz(i)*dyz(j)*(aw(3,1,i,j)-aw(3,0,i,j))                       !Vert velocity top/bottom face                  
                  
        mat(4,i,j) = -  dxz(i)*s%hv(i,j-1)  *av(0,i,j-1)                              &
                   -    dzsv*av(0,i,j-1) + dzbv*av(0,i,j-1)                             & !kin. boun. top/bot U
                   +    dxz(i)*dyz(j)*(aw(4,1,i,j)-aw(4,0,i,j))

        mat(5,i,j) =    dxz(i)*s%hv(i,j)  *av(1,i,j)                                  &
                   -    dzsv*av(1,i,j) + dzbv*av(1,i,j)                                 & !kin. boun. top/bot U
                   +    dxz(i)*dyz(j)*(aw(5,1,i,j)-aw(5,0,i,j))
       
        rhs(i,j)   = -  dyz(j)*( s%hu(i,j)*aur(i,j) - s%hu(i-1,j  )*aur(i-1,  j) )    & !subs U left/right face
                   -    dxz(i)*( s%hv(i,j)*avr(i,j) - s%hv(i  ,j-1)*avr(i  ,j-1) )    & !subs V left/rigth face
                   +    dzsu*(aur(i,j) + aur(i-1,j))  + dzsv*(avr(i,j) + avr(i,j-1))    & !kin. boun. top
                   -    dzbu*(aur(i,j) + aur(i-1,j))  - dzbv*(avr(i,j) + avr(i,j-1))    & !kin. boun. bot
                   -    dxz(i)*dyz(j)*(awr(1,i,j)-awr(0,i,j))                         !Vert velocity top/bottom face
      else
        mat(1  ,i,j) = 1.0_rKind
        mat(2:5,i,j) = 0.0_rKind
        rhs(i,j)     = 0.0_rKind
      endif     
    enddo
  enddo
  !call timer_stop(timer_flow_nonh_subs)
  
  !call timer_stop(timer_flow_nonh_mat)
  
  !Solve matrix
  !call timer_start(timer_flow_nonh_solv)
  dp = 0.0_rKind
  call solver_solvemat( mat  , rhs   , dp , s%nx, s%ny,par)
  !call timer_stop(timer_flow_nonh_solv)
  
  s%pres = s%pres + dp
  
  !Correct u/v/w

  !U
  !call timer_start(timer_flow_nonh_coru)
  do j=2,s%ny
    do i=2,s%nx-1
      s%uu(i,j) = aur(i,j) + au(1,i,j)*dp(i+1,j)+au(0,i,j)*dp(i,j)
    enddo
  enddo    
  !call timer_stop(timer_flow_nonh_coru)

  !v
  !call timer_start(timer_flow_nonh_corv)
  do j=2,s%ny-1
    do i=2,s%nx
     s%vv(i,j) = avr(i,j) + av(1,i,j)*dp(i,j+1)+av(0,i,j)*dp(i,j)  
    enddo
  enddo    
  !call timer_stop(timer_flow_nonh_corv)
  
  !W
  !call timer_start(timer_flow_nonh_corw)
  do j=2,s%ny
    do i=2,s%nx
      s%ws(i,j) = awr(1,i,j) + dp(i , j)  * aw(1,1,i,j)                             &
                             + dp(i-1,j)  * aw(2,1,i,j) + dp(i+1,j  ) * aw(3,1,i,j) &
                             + dp(i  ,j-1)* aw(4,1,i,j) + dp(i  ,j+1) * aw(5,1,i,j)
      s%wb(i,j) = awr(0,i,j) + dp(i , j)  * aw(1,0,i,j)                             &
                             + dp(i-1,j)  * aw(2,0,i,j) + dp(i+1,j  ) * aw(3,0,i,j) &
                             + dp(i  ,j-1)* aw(4,0,i,j) + dp(i  ,j+1) * aw(5,0,i,j)                                    
    enddo
  enddo    

  !call timer_stop(timer_flow_nonh_corw)    
  
  !Assign boundaries
  s%ws(:,1)      = s%ws(:,2)
  s%ws(:,s%ny+1) = s%ws(:,s%ny)
  s%ws(1,:)      = s%ws(2,:)
  s%ws(s%nx+1,:) = s%ws(s%nx,:)
  
  ws_old=s%ws
  wb_old=s%wb
  !call timer_stop(timer_flow_nonh)
end subroutine nonh_cor





!
!==============================================================================
  subroutine nonh_free()
!==============================================================================    
!

!-------------------------------------------------------------------------------
!                             DECLARATIONS
!-------------------------------------------------------------------------------

!--------------------------        PURPOSE         ----------------------------

!  Releases resources


!--------------------------     ARGUMENTS          ----------------------------

!                                - NONE -


!--------------------------     LOCAL VARIABLES    ----------------------------
 
!                                - NONE -

!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------   

    if (allocated(au  )) deallocate(au)
    if (allocated(av  )) deallocate(av)
    if (allocated(aw  )) deallocate(aw)
    if (allocated(aur )) deallocate(aur)
    if (allocated(avr )) deallocate(avr)
    if (allocated(awr )) deallocate(awr)
    if (allocated(mat))  deallocate(mat)
    if (allocated(rhs))  deallocate(rhs)
    if (allocated(dp))  deallocate(dp)           
  end subroutine nonh_free
  
#endif
end module nonh_module
