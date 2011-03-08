!==============================================================================
!                               MODULE NH_MAIN                        
!==============================================================================

! DATE               AUTHOR               CHANGES        
!
! october 2009       Pieter Bart Smit     New module


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
  real (kind=rKind), allocatable, dimension(:,:,:)   :: aws    !Pressure coeficients for ws(k,i,j)
  real (kind=rKind), allocatable, dimension(:,:,:)   :: awb    !Pressure coeficients for wb(k,i,j)

  
  real (kind=rKind), allocatable, dimension(:,:)     :: aur    !RHS for u(i,j)
  real (kind=rKind), allocatable, dimension(:,:)     :: avr    !RHS for v(i,j)
  real (kind=rKind), allocatable, dimension(:,:)     :: awbr   !Pressure coeficients for wb(k,i,j)  
  real (kind=rKind), allocatable, dimension(:,:)     :: awsr   !Pressure coeficients for ws(k,i,j)
    
  real (kind=rKind), allocatable, dimension(:,:,:)   :: mat   !Pressure matrix
  real (kind=rKind), allocatable, dimension(:,:)     :: rhs   !RHS of the pressure matrix
  
  real (kind=rKind), allocatable, dimension(:,:)     :: zsu    !free surface/bottom in velocity points
  real (kind=rKind), allocatable, dimension(:,:)     :: zsv    !free surface/bottom in velocity points
  
  real (kind=rKind), allocatable, dimension(:,:)     :: zbu    !free surface/bottom in velocity points
  real (kind=rKind), allocatable, dimension(:,:)     :: zbv    !free surface/bottom in velocity points
  
  real (kind=rKind), allocatable, dimension(:,:)     :: dp     !Change in pressure
    
  real (kind=rKind), allocatable, dimension(:)       ::  dxz   
  real (kind=rKind), allocatable, dimension(:)       ::  dyz   
  real (kind=rKind), allocatable, dimension(:)       ::  dxu   
  real (kind=rKind), allocatable, dimension(:)       ::  dyv   
  real (kind=rKind), allocatable, dimension(:)       ::  ddxz  
  real (kind=rKind), allocatable, dimension(:)       ::  ddyz   
  real (kind=rKind), allocatable, dimension(:)       ::  ddxu  
  real (kind=rKind), allocatable, dimension(:)       ::  ddyv   
  
  real (kind=rKind), allocatable, dimension(:,:)     :: wcoef
  
  real (kind=rKind), allocatable, dimension(:,:)     :: Wm
  real (kind=rKind), allocatable, dimension(:,:)     :: Wm_old
  
  integer(kind=iKind),allocatable,dimension(:,:)     :: nonhU
  integer(kind=iKind),allocatable,dimension(:,:)     :: nonhV
  integer(kind=iKind),allocatable,dimension(:,:)     :: nonhZ

  !--- PUBLIC SUBROUTINES ---
  public nonh_init
  public nonh_cor
  public nonh_free
  public nonh_explicit
  
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
    use wave_timestep_module, only: dispersion

!--------------------------     ARGUMENTS          ----------------------------
!
    type(spacepars) ,intent(inout) :: s
    type(parameters),intent(in)    :: par     
!

!--------------------------     LOCAL VARIABLES    ----------------------------
 
    integer(kind=iKind)        :: iw,i,ie,jn,j,js
    real   (kind=rKind)        :: d,sigma,k
!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------
   
    open(unit=PrintFileUnit,file=trim(PrintFileName))
        
    allocate(au (0:1,s%nx+1,s%ny+1)); au = 0.0_rKind
    allocate(av (0:1,s%nx+1,s%ny+1)); av = 0.0_rKind  
    allocate(aws (5,s%nx+1,s%ny+1)) ; aws = 0.0_rKind
    allocate(awb (5,s%nx+1,s%ny+1)) ; awb = 0.0_rKind        
    allocate(aur(   s%nx+1,s%ny+1)) ; aur = 0.0_rKind
    allocate(avr(   s%nx+1,s%ny+1)) ; avr = 0.0_rKind
    allocate(awsr(s%nx+1,s%ny+1))   ; awsr = 0.0_rKind
    allocate(awbr(s%nx+1,s%ny+1))   ; awbr = 0.0_rKind        
    allocate(mat(5,s%nx+1,s%ny+1))  ;  mat = 0.0_rKind
    allocate(rhs(  s%nx+1,s%ny+1))  ;  rhs = 0.0_rKind    
    allocate(zbu(  s%nx+1,s%ny+1))  ;  zbu = 0.0_rKind
    allocate(zbv(  s%nx+1,s%ny+1))  ;  zbv = 0.0_rKind
    allocate(zsu(  s%nx+1,s%ny+1))  ;  zsu = 0.0_rKind
    allocate(zsv(  s%nx+1,s%ny+1))  ;  zsv = 0.0_rKind        
    allocate(dxz (  s%nx+1))        ;  dxz  = 0.0_rKind
    allocate(dyz (  s%ny+1))        ;  dyz  = 0.0_rKind
    allocate(dxu (  s%nx+1))        ;  dxu  = 0.0_rKind
    allocate(dyv (  s%ny+1))        ;  dyv  = 0.0_rKind
    allocate(ddxz (  s%nx+1))       ;  ddxz  = 0.0_rKind
    allocate(ddyz (  s%ny+1))       ;  ddyz  = 0.0_rKind
    allocate(ddxu (  s%nx+1))       ;  ddxu  = 0.0_rKind
    allocate(ddyv (  s%ny+1))       ;  ddyv  = 0.0_rKind          
    allocate(nonhU(  s%nx+1,s%ny+1)); nonhU = 1
    allocate(nonhV(  s%nx+1,s%ny+1)); nonhV = 1
    allocate(nonhZ(  s%nx+1,s%ny+1)); nonhZ = 1    
    allocate(dp(s%nx+1,s%ny+1))     ; dp = 0.0_rKind
    allocate(Wm(s%nx+1,s%ny+1))     ;Wm     = 0.0_rKind
    allocate(Wm_old(s%nx+1,s%ny+1)) ;Wm_old = 0.0_rKind
    allocate(Wcoef(s%nx+1,s%ny+1))  ;Wcoef = 1.0d0
    
    if (.not. associated(s%pres)) then   
      allocate(s%pres(s%nx+1,s%ny+1)); s%pres = 0.0_rKind   
      allocate(s%ws(s%nx+1,s%ny+1)); s%ws = 0.0_rKind
      allocate(s%wb(s%nx+1,s%ny+1)); s%wb = 0.0_rKind
    endif  
    
    
    !Initialize grid variables
    dyz  = s%dnz(1,:)
    dyv  = s%dnv(1,:)
    ddyz = 1.0_rKind/dyz
    ddyv = 1.0_rKind/dyv

    dxu  = s%dsu(:,1)
    dxz  = s%dsz(:,1)
    ddxu = 1.0_rKind/dxu
    ddxz = 1.0_rKind/dxz
    
    mat(1,1,:)      = 1.0_rKind
    mat(1,s%nx+1,:) = 1.0_rKind
    mat(1,:,1)      = 1.0_rKind
    mat(1,:,s%ny+1) = 1.0_rKind
    
    initialized   = .true. 
    
    sigma = 2.0_rkind*par%px/par%Topt
    if (par%dispc <= 0.0_rKind) then
      !Calculate the optimum value of the dispersion coeficiant.
      do j=1,s%ny+1
        do i=1,s%nx+1      
          d          = par%zs0-s%zb(i,j)
          k = disper(sigma,d,par%g,2.0_rKind*par%px,0.0001_rKind)
          if (d>0.0_rKind) then          
            wcoef(i,j) = 1.0_rKind / (  4.0_rKind*par%g / (d*sigma**2) - 4.0_rKind / ((k*d)**2)  )
          else
            wcoef(i,j) = par%dispc
          endif  
        enddo
      enddo
    endif
   
    !Initialize levels at u/v points (zu,zbu etc.)
    call zuzv(s)

  end subroutine nonh_init

!
!==============================================================================
  subroutine nonh_cor(s,par)
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
!--------------------------     ARGUMENTS          ----------------------------

    type(spacepars) ,intent(inout)                       :: s
    type(parameters),intent(in)                          :: par 

!--------------------------     LOCAL VARIABLES    ----------------------------
 
    !Indices
    integer(kind=iKind)                     :: i               !Index variables
    integer(kind=iKind)                     :: j               !Index variables
    
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
  avr(:,s%ny)   = s%vv(:,s%ny)

  !Built pressure coefficients for W
  
  !AW Bottom
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
  
  do j=2,s%ny
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
                   -    dzs_s*av(0,i,j) - dzs_n*av(1,i-1,j)                           & !kin. boun. top
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
  
  !Solve matrix
  !call timer_start(timer_flow_nonh_solv)
  dp = 0.0_rKind
  call solver_solvemat( mat  , rhs   , dp , s%nx, s%ny,par)
  !call timer_stop(timer_flow_nonh_solv)

  
  s%pres = s%pres + dp
  
  !Correct u/v/w

  !U
  do j=2,s%ny
    do i=2,s%nx-1
      s%uu(i,j) = aur(i,j) + au(1,i,j)*dp(i+1,j)+au(0,i,j)*dp(i,j)
    enddo
  enddo    

  !v
  do j=2,s%ny-1
    do i=2,s%nx
     s%vv(i,j) = avr(i,j) + av(1,i,j)*dp(i,j+1)+av(0,i,j)*dp(i,j)  
    enddo
  enddo    
  
  !W
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
  
  !Assign boundaries
  s%ws(:,1)      = s%ws(:,2)
  s%ws(:,s%ny+1) = s%ws(:,s%ny)
  !s%ws(1,:)      = s%ws(2,:)
  s%ws(s%nx+1,:) = s%ws(s%nx,:)
  
  s%wb(:,1)      = s%wb(:,2)
  s%wb(:,s%ny+1) = s%wb(:,s%ny)
  s%wb(1,:)      = s%wb(2,:)
  s%wb(s%nx+1,:) = s%wb(s%nx,:)  
  
  Wm_old = .5_rKind*(s%ws+s%wb)
end subroutine nonh_cor

subroutine nonh_explicit(s,par,nuh)
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
    
    real(kind=rKind),dimension(s%nx+1,s%ny+1)            :: nuh

!--------------------------     LOCAL VARIABLES    ----------------------------
 
    !Indices
    integer(kind=iKind)                     :: i,ie,iee,iw               !Index variables
    integer(kind=iKind)                     :: j,js       
    
    real(kind=rKind)                        :: dwdx1    !Gradient of vertical velocity in x-dir at i+1/2,j
    real(kind=rKind)                        :: dwdx2    !Gradient of vertical velocity in x-dir at i-1/2,j
    real(kind=rKind)                        :: dwdy1    !Gradient of vertical velocity in x-dir at i    ,j+1/2
    real(kind=rKind)                        :: dwdy2    !Gradient of vertical velocity in x-dir at i    ,j-1/2   
    real(kind=rKind)                        :: Vol    

  if (.not. initialized) then
    call nonh_init(s,par)
  endif      
   
!
! Determine if a velocity point will be included in the nonh pressure matrix, do not include if:
!
! (1)  The point is dry
! (2)  The relative wave length kd of the smallest possible wave (L=2dx) is smaller than kdmin 
! (3)  The interpolated waterlevel zs is below the bottom (steep cliffs with overwash situations)


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
  do j=2,s%ny
    do i=2,s%nx
      if (max(nonhV(i,j),nonhV(i,j-1),nonhU(i,j),nonhU(i-1,j)) > 0) then
        nonhZ(i,j) = 1
      else
        nonhZ(i,j) = 0      
      endif
    enddo
  enddo


 
 !Calculate explicit part average vertical momentum (advection)
  do j=2,s%ny
    do i=2,s%nx
      if (nonhZ(i,j) == 1) then
       Wm(i,j) = Wm_old(i,j) - par%dt*( ddxu(i-1)*max(s%qx(i-1,j  ),0.0_rKind)*(Wm_old(i  ,j  )-Wm_old(i-1,j  ))/s%hh(i,j)   &
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

  !Calculate explicit part vertical viscocity
  do j=2,s%ny
    do i=2,s%nx
      dwdx1 = .5d0*(nuh(i-1,j  )+nuh(i,j  ))*s%hu(i-1,j  )*(Wm_old(i  ,j  )-Wm_old(i-1,j  ))*ddxu(i-1)
      dwdx2 = .5d0*(nuh(i+1,j  )+nuh(i,j  ))*s%hu(i  ,j  )*(Wm_old(i+1,j  )-Wm_old(i  ,j  ))*ddxu(i)
      dwdy1 = nuh(i  ,j-1)*s%hu(i  ,j-1)*(Wm_old(i  ,j  )-Wm_old(i  ,j-1))*ddyv(j-1)
      dwdy2 = nuh(i  ,j  )*s%hu(i  ,j  )*(Wm_old(i  ,j+1)-Wm_old(i  ,j  ))*ddyv(j)
      Wm(i,j) = Wm(i,j)   + (1.0_rKind/s%hh(i,j))*par%dt*(dwdx2-dwdx1)*ddxz(i)*real(s%wetu(i,j)*s%wetu(i-1,j),rKind) &
                          + (1.0_rKind/s%hh(i,j))*par%dt*(dwdy2-dwdy1)*ddyz(j)*real(s%wetv(i,j)*s%wetv(i,j-1),rKind)
    enddo
  enddo 

  do j=2,s%ny
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
  au(:,1,:)      = 0.0_rKind
  au(:,s%nx,:)   = 0.0_rKind

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
    av(:,:,1)     = 0.0_rKind
    av(:,:,s%ny)  = 0.0_rKind
  endif  
       

 !Include explicit approximation for pressure in s%uu and s%vv   and Wm
  if (par%secorder == 1) then 
    do j=2,s%ny
      do i=2,s%nx-1
        s%uu(i,j) = s%uu(i,j) + au(1,i,j) * s%pres(i+1,j) + au(0,i,j) * s%pres(i  ,j)  
      enddo
    enddo

    do j=2,s%ny-1
      do i=2,s%nx
        s%vv(i,j) = s%vv(i,j) + av(1,i,j) * s%pres(i,j+1) + av(0,i,j) * s%pres(i,j  )
      enddo
    enddo
    
    do j=2,s%ny
      do i=2,s%nx
        if (nonhZ(i,j) == 1) then
          Wm(i,j) = Wm(i,j)   + Wcoef(i,j)*par%dt * s%pres(i,j)/s%hh(i,j)
        endif  
      enddo
    enddo
    call flow_secondorder_advW(s,par,Wm,Wm_old)
  endif

end subroutine nonh_explicit

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
!--------------------------     ARGUMENTS          ----------------------------

    type(spacepars) ,intent(inout)                       :: s
!--------------------------     LOCAL VARIABLES    ----------------------------
 
    integer(kind=iKind)                                  :: i
    integer(kind=iKind)                                  :: j    

!-------------------------------------------------------------------------------
!                             IMPLEMENTATION
!-------------------------------------------------------------------------------      
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
    if (s%ny>2) then
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
    else
      zsv(:,1) = s%zs(:,2)
      zsv(:,2) = s%zs(:,2)
      zsv(:,3) = s%zs(:,2)
    endif

    !Bottom location in U point
    do j=1,s%ny+1
      do i=1,s%nx
        zbu(i,j) = zsu(i,j)-s%hu(i,j)
      enddo
    enddo
    zbu(s%nx+1,:) = s%zb(s%nx+1,:)

    !Bottom location in V point  
    if (s%ny>2) then
      do j=1,s%ny
        do i=1,s%nx+1
          zbv(i,j) = zsv(i,j)-s%hv(i,j)
        enddo
      enddo  
      zbv(:,s%ny+1) = s%zb(:,s%ny+1)    
    else
      zbv = s%zb
    endif  
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





!
!==============================================================================
  subroutine nonh_free()
!==============================================================================    
!

! DATE               AUTHOR               CHANGES        
!
! October 2010       Pieter Bart Smit     New Subroutine

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
    
if (allocated(au)    ) deallocate(au  )  
if (allocated(av)    ) deallocate(av  )  
if (allocated(aws)   ) deallocate(aws )  
if (allocated(awb)   ) deallocate(awb )  
if (allocated(aur)   ) deallocate(aur )  
if (allocated(avr)   ) deallocate(avr )  
if (allocated(awbr)  ) deallocate(awbr)  
if (allocated(awsr)  ) deallocate(awsr)  
if (allocated(mat)   ) deallocate(mat )  
if (allocated(rhs)   ) deallocate(rhs )  
if (allocated(zsu)   ) deallocate(zsu )  
if (allocated(zsv)   ) deallocate(zsv )  
if (allocated(zbu)   ) deallocate(zbu )  
if (allocated(zbv)   ) deallocate(zbv )  
if (allocated(dp)    ) deallocate(dp  )  
if (allocated(dxz)   ) deallocate(dxz )  
if (allocated(dyz)   ) deallocate(dyz )  
if (allocated(dxu)   ) deallocate(dxu )  
if (allocated(dyv)   ) deallocate(dyv )  
if (allocated(ddxz)  ) deallocate(ddxz) 
if (allocated(ddyz)  ) deallocate(ddyz)  
if (allocated(ddxu)  ) deallocate(ddxu) 
if (allocated(ddyv)  ) deallocate(ddyv)   
if (allocated(Wm)    ) deallocate(Wm)
if (allocated(Wm_old)) deallocate(Wm_old)
if (allocated(nonhU) ) deallocate(nonhU)
if (allocated(nonhV) ) deallocate(nonhV)
if (allocated(nonhZ) ) deallocate(nonhZ)
    
    
    
    
  end subroutine nonh_free
  
#endif
end module nonh_module
