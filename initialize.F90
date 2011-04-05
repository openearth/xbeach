module initialize 
contains
  subroutine grid_bathy(s,par)                                         

  use params   
  use spaceparams                       
  use xmpi_module
  use general_mpi_module
  use readkey_module
  use logging_module

  implicit none                                                            

  type(spacepars)                :: s               
  type(parameters)               :: par               

  integer                        :: i
  integer                        :: j
  integer                        :: itheta
  real*8                         :: degrad

  s%nx      = par%nx
  s%ny      = par%ny
  s%dx      = par%dx
  s%dy      = par%dy
  s%xori    = par%xori
  s%yori    = par%yori
  s%alfa    = par%alfa
  s%posdwn  = par%posdwn
  s%vardx   = par%vardx

  if (s%alfa.lt.0) then 
    s%alfa = 360.d0+s%alfa
  endif

  s%alfa  = s%alfa*atan(1.0d0)/45.d0
  ! Robert: huh?
  s%posdwn = s%posdwn*sign(s%posdwn,1.d0)
  ! end huh?

  if (xmaster) then
    allocate(s%x(1:s%nx+1,1:s%ny+1))
    allocate(s%y(1:s%nx+1,1:s%ny+1))
    allocate(s%xz(1:s%nx+1,1:s%ny+1))
    allocate(s%yz(1:s%nx+1,1:s%ny+1))
    allocate(s%xu(1:s%nx+1,1:s%ny+1))
    allocate(s%yu(1:s%nx+1,1:s%ny+1))
    allocate(s%xv(1:s%nx+1,1:s%ny+1))
    allocate(s%yv(1:s%nx+1,1:s%ny+1))
    allocate(s%zb(1:s%nx+1,1:s%ny+1))
    allocate(s%zb0(1:s%nx+1,1:s%ny+1))
    allocate(s%dzbdx(1:s%nx+1,1:s%ny+1))
    allocate(s%dzbdy(1:s%nx+1,1:s%ny+1))
    allocate(s%dsu(1:s%nx+1,1:s%ny+1))
    allocate(s%dsv(1:s%nx+1,1:s%ny+1))
    allocate(s%dsz(1:s%nx+1,1:s%ny+1))
    allocate(s%dsc(1:s%nx+1,1:s%ny+1))
    allocate(s%dnu(1:s%nx+1,1:s%ny+1))
    allocate(s%dnv(1:s%nx+1,1:s%ny+1))
    allocate(s%dnz(1:s%nx+1,1:s%ny+1))
    allocate(s%dnc(1:s%nx+1,1:s%ny+1))
    allocate(s%dsdnui(1:s%nx+1,1:s%ny+1))
    allocate(s%dsdnvi(1:s%nx+1,1:s%ny+1))
    allocate(s%dsdnzi(1:s%nx+1,1:s%ny+1))
    allocate(s%alfau(1:s%nx+1,1:s%ny+1))
    allocate(s%alfav(1:s%nx+1,1:s%ny+1))
    allocate(s%alfaz(1:s%nx+1,1:s%ny+1))
    allocate(s%sdist(1:s%nx+1,1:s%ny+1))
    allocate(s%ndist(1:s%nx+1,1:s%ny+1))
  endif
  !
  ! Create grid
  !
  if(s%vardx==0)then
    if (xmaster) then
      open(31,file=par%depfile)
      do j=1,s%ny+1
        read(31,*)(s%zb(i,j),i=1,s%nx+1)
      end do
      close(31)
      do j=1,s%ny+1
        do i=1,s%nx+1
          s%x(i,j)=(i-1)*s%dx
          s%y(i,j)=(j-1)*s%dy
        end do
      end do
    endif
  elseif(s%vardx==1)then   ! Robert keep vardx == 1 for backwards compatibility??
    if (xmaster) then
      open(31,file=par%depfile)
      open(32,file=par%xfile)
      open(33,file=par%yfile)
      read(31,*)((s%zb(i,j),i=1,s%nx+1),j=1,s%ny+1)
      read(32,*)((s%x(i,j),i=1,s%nx+1),j=1,s%ny+1)
      read(33,*)((s%y(i,j),i=1,s%nx+1),j=1,s%ny+1)
      close(31)
      close(32)
      close(33)
    endif
  else 
    call writelog('esl','','Invalid value for vardx: ',par%vardx)
    call halt_program
  endif

  if(xmaster) then
      
       s%zb=-s%zb*s%posdwn
       ! Make sure that at the lateral boundaries the bathymetry is alongshore uniform
       if (s%ny>0) then
         s%zb(:,1) = s%zb(:,2)
         s%zb(:,s%ny+1) = s%zb(:,s%ny)
       endif
       s%zb(1,:)=s%zb(2,:)
       s%zb(s%nx+1,:)=s%zb(s%nx,:)
      
       call gridprops (s)
	  
       s%zb0 = s%zb
       !
       ! Specify theta-grid
       !
       !
       ! from Nautical wave directions in degrees to Cartesian wave directions in radian !!!
       !
!      s%theta0=(1.5d0*par%px-s%alfa)-par%dir0*atan(1.d0)/45.d0 ! Updated in waveparams.f90 for instat 4,5,6,7
       s%theta0=(1.5d0*par%px)-par%dir0*atan(1.d0)/45.d0 ! Updated in waveparams.f90 for instat 4,5,6,7
       if (s%theta0<-par%px) s%theta0=s%theta0+2.d0*par%px
       if (s%theta0> par%px) s%theta0=s%theta0-2.d0*par%px
       
       degrad=par%px/180.d0
       
       if (par%thetanaut==1) then  
            s%thetamin=(270-par%thetamax)*degrad
            s%thetamax=(270-par%thetamin)*degrad
       else
            ! rotate theta grid to world coordinates for backwards compatibility
            s%thetamin=par%thetamin+s%alfa/degrad
            s%thetamax=par%thetamax+s%alfa/degrad
            
            s%thetamin=s%thetamin*degrad
            s%thetamax=s%thetamax*degrad
       endif
       
       if (s%thetamax>par%px) then
            s%thetamax=s%thetamax-2*par%px
            s%thetamin=s%thetamin-2*par%px
       endif
       if (s%thetamin<-par%px) then
            s%thetamax=s%thetamax+2*par%px
            s%thetamin=s%thetamin+2*par%px
       endif
       
       s%dtheta=par%dtheta*degrad
       s%ntheta=nint((s%thetamax-s%thetamin)/s%dtheta)

       allocate(s%theta(1:s%ntheta))
       allocate(s%thet(1:s%nx+1,1:s%ny+1,1:s%ntheta))
       allocate(s%costh(1:s%nx+1,1:s%ny+1,1:s%ntheta))
       allocate(s%sinth(1:s%nx+1,1:s%ny+1,1:s%ntheta))

       do itheta=1,s%ntheta
           s%theta(itheta)=s%thetamin+s%dtheta/2+s%dtheta*(itheta-1)
       end do

       do itheta=1,s%ntheta
          do j=1,s%ny+1
             do i=1,s%nx+1
                s%thet(i,j,itheta) = s%theta(itheta)
                s%costh(i,j,itheta)=cos(s%theta(itheta)-s%alfaz(i,j))
                s%sinth(i,j,itheta)=sin(s%theta(itheta)-s%alfaz(i,j))
             enddo
          enddo
       enddo

    endif

    if (xmaster) then
       ! Initialize dzbdx, dzbdy
       do j=1,s%ny+1
          do i=1,s%nx
             s%dzbdx(i,j)=(s%zb(i+1,j)-s%zb(i,j))/s%dsu(i,j)
          enddo
       enddo
       ! dummy, needed to keep compiler happy
       s%dzbdx(s%nx+1,:)=s%dzbdx(s%nx,:)

       do j=1,s%ny
          do i=1,s%nx+1
             s%dzbdy(i,j)=(s%zb(i,j+1)-s%zb(i,j))/s%dnv(i,j)
          enddo
       enddo
       if (s%ny>0) then
          s%dzbdy(:,s%ny+1)=s%dzbdy(:,s%ny)
       endif
    endif
  end subroutine grid_bathy


  subroutine wave_init (s,par)
    use params
    use spaceparams
    use readkey_module
    use xmpi_module
    use wave_timestep_module

    IMPLICIT NONE

    type(spacepars)                     :: s
    type(parameters)                    :: par

    integer                             :: i
    integer                             :: j
    integer                             :: itheta

    include 's.ind'
    include 's.inp'
    
    allocate(s%thetamean(1:nx+1,1:ny+1))
    allocate(s%Fx(1:nx+1,1:ny+1))
    allocate(s%Fy(1:nx+1,1:ny+1))
    allocate(s%Sxx(1:nx+1,1:ny+1))
    allocate(s%Sxy(1:nx+1,1:ny+1))
    allocate(s%Syy(1:nx+1,1:ny+1))
    allocate(s%n(1:nx+1,1:ny+1))
    allocate(s%H(1:nx+1,1:ny+1))
    allocate(s%cgx(1:nx+1,1:ny+1,1:ntheta))
    allocate(s%cgy(1:nx+1,1:ny+1,1:ntheta))
    allocate(s%cx(1:nx+1,1:ny+1,1:ntheta))
    allocate(s%cy(1:nx+1,1:ny+1,1:ntheta))
    allocate(s%ctheta(1:nx+1,1:ny+1,1:ntheta))
    allocate(s%sigt(1:nx+1,1:ny+1,1:ntheta))
    allocate(s%ee(1:nx+1,1:ny+1,1:ntheta))
    allocate(s%rr(1:nx+1,1:ny+1,1:ntheta))
    allocate(s%sigm(1:nx+1,1:ny+1))
    allocate(s%c(1:nx+1,1:ny+1))
    allocate(s%cg(1:nx+1,1:ny+1))
    allocate(s%k(1:nx+1,1:ny+1))
    allocate(s%ui(1:nx+1,1:ny+1))
    allocate(s%vi(1:nx+1,1:ny+1))
    allocate(s%E(1:nx+1,1:ny+1)) 
    allocate(s%R(1:nx+1,1:ny+1)) 
    allocate(s%urms(1:nx+1,1:ny+1)) 
    allocate(s%D(1:nx+1,1:ny+1)) 
    allocate(s%Df(1:nx+1,1:ny+1)) 
    allocate(s%Dp(1:nx+1,1:ny+1)) 
    allocate(s%Qb(1:nx+1,1:ny+1)) 
    allocate(s%ust(1:nx+1,1:ny+1)) 
    allocate(s%tm(1:nx+1,1:ny+1)) 
    allocate(s%uwf(1:nx+1,1:ny+1)) 
    allocate(s%vwf(1:nx+1,1:ny+1)) 
    allocate(s%ustr(1:nx+1,1:ny+1)) 
    allocate(s%usd(1:nx+1,1:ny+1))
    allocate(s%bi(1:ny+1))
    allocate(s%DR(1:nx+1,1:ny+1)) 
    allocate(s%umwci       (1:nx+1,1:ny+1))
    allocate(s%vmwci       (1:nx+1,1:ny+1))
    allocate(s%zswci       (1:nx+1,1:ny+1))
    allocate(s%BR(1:nx+1,1:ny+1))
    !
    ! Initial condition
    !
    do itheta=1,ntheta
       do j=1,ny+1
          do i=1,nx+1
             s%ee(i,j,itheta)=0.d0
          end do
       end do
    end do

    s%thetamean = 0.d0
    s%Fx        = 0.d0
    s%Fy        = 0.d0
    s%Sxx       = 0.d0
    s%Sxy       = 0.d0
    s%Syy       = 0.d0
    s%n         = 0.d0
    s%H         = 0.d0
    s%cgx       = 0.d0
    s%cgy       = 0.d0
    s%cx        = 0.d0
    s%cy        = 0.d0
    s%ctheta    = 0.d0
    s%sigt      = 0.d0
    s%rr        = 0.d0
    s%sigm      = 0.d0
    s%c         = 0.d0
    s%cg        = 0.d0
    s%k         = 0.d0
    s%ui        = 0.d0
    s%vi        = 0.d0
    s%E         = 0.d0
    s%R         = 0.d0
    s%urms      = 0.d0
    s%D         = 0.d0
    s%Qb        = 0.d0
    s%ust       = 0.d0
    s%tm        = 0.d0
    s%uwf       = 0.d0
    s%vwf       = 0.d0
    s%ustr      = 0.d0
    s%usd       = 0.d0
    s%bi        = 0.d0
    s%DR        = 0.d0
    s%BR        = par%Beta


    ! introduce intrinsic frequencies for wave action
    if (  trim(par%instat)=='jons' .or. &
         trim(par%instat)=='swan' .or. &
         trim(par%instat)=='vardens' .or. &
         trim(par%instat)=='reuse' .or. &
         trim(par%instat)=='nonh' &
         ) par%Trep=10.d0 
    !Robert
    ! incorrect values are computed below for instat = 4/5/6/7
    ! in this case right values are computed in wave params.f90
    do itheta=1,ntheta
       s%sigt(:,:,itheta) = 2*par%px/par%Trep
    end do
    s%sigm = sum(s%sigt,3)/ntheta
    call dispersion(par,s)


  end subroutine wave_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine flow_init (s,par)
    use params
    use spaceparams
    use readkey_module
    use logging_module
    use interp
    use xmpi_module

    IMPLICIT NONE

    type(spacepars)                         :: s
    type(parameters), intent(in)            :: par

    integer*4                               :: i,j,ig,indt
    logical                                 :: exists
    logical                                 :: offshoreregime
    integer                                 :: indoff,indbay

    real*8,dimension(:),allocatable         :: xzs0,yzs0,szs0
    
    include 's.ind'
    include 's.inp'
    
    allocate(s%zs(1:s%nx+1,1:s%ny+1))
    allocate(s%dzsdt(1:s%nx+1,1:s%ny+1))
    allocate(s%dzsdx(1:s%nx+1,1:s%ny+1))
    allocate(s%dzsdy(1:s%nx+1,1:s%ny+1))
    allocate(s%dzbdt(1:s%nx+1,1:s%ny+1))
    allocate(s%uu(1:s%nx+1,1:s%ny+1))
    allocate(s%vv(1:s%nx+1,1:s%ny+1))
    allocate(s%qx(1:s%nx+1,1:s%ny+1))
    allocate(s%qy(1:s%nx+1,1:s%ny+1))
    allocate(s%sedero(1:s%nx+1,1:s%ny+1))
    allocate(s%ueu(1:s%nx+1,1:s%ny+1)) 
    allocate(s%vev(1:s%nx+1,1:s%ny+1)) 
    allocate(s%vmagu(1:s%nx+1,1:s%ny+1)) 
    allocate(s%vmagv(1:s%nx+1,1:s%ny+1)) 
    allocate(s%vmageu(1:s%nx+1,1:s%ny+1)) 
    allocate(s%vmagev(1:s%nx+1,1:s%ny+1)) 
    allocate(s%u(1:s%nx+1,1:s%ny+1)) 
    allocate(s%v(1:s%nx+1,1:s%ny+1)) 
    allocate(s%ue(1:s%nx+1,1:s%ny+1)) 
    allocate(s%ve(1:s%nx+1,1:s%ny+1)) 
    allocate(s%hold(1:s%nx+1,1:s%ny+1)) 
    allocate(s%wetu(1:s%nx+1,1:s%ny+1)) 
    allocate(s%wetv(1:s%nx+1,1:s%ny+1)) 
    allocate(s%wetz(1:s%nx+1,1:s%ny+1))
    allocate(s%hh(1:s%nx+1,1:s%ny+1))
    allocate(s%hu(1:s%nx+1,1:s%ny+1)) 
    allocate(s%hv(1:s%nx+1,1:s%ny+1))
    allocate(s%hum(1:s%nx+1,1:s%ny+1)) 
    allocate(s%hvm(1:s%nx+1,1:s%ny+1))
    allocate(s%vu(1:s%nx+1,1:s%ny+1))
    allocate(s%uv(1:s%nx+1,1:s%ny+1))
    allocate(s%maxzs(1:s%nx+1,1:s%ny+1))
    allocate(s%minzs(1:s%nx+1,1:s%ny+1))
    allocate(s%taubx(1:s%nx+1,1:s%ny+1))
    allocate(s%tauby(1:s%nx+1,1:s%ny+1))
    allocate(s%ws(1:s%nx+1,1:s%ny+1))
    allocate(s%wb(1:s%nx+1,1:s%ny+1))
    allocate(s%nuh(1:s%nx+1,1:s%ny+1))  
    allocate(s%pres(1:s%nx+1,1:s%ny+1))
    allocate(s%wi(2,1:s%ny+1))
    allocate(s%zi(2,1:s%ny+1))
    allocate(s%cf(1:s%nx+1,1:s%ny+1))
    allocate(s%zs0(1:s%nx+1,1:s%ny+1)) 
    allocate(s%zs0fac(1:s%nx+1,1:s%ny+1,2))
    allocate(s%wm(1:s%nx+1,1:s%ny+1))
    allocate(s%umean(2,1:s%ny+1))
    allocate(s%vmean(2,1:s%ny+1))
    allocate(s%xyzs01(2))
    allocate(s%xyzs02(2))
    allocate(s%xyzs03(2))
    allocate(s%xyzs04(2))

    allocate(szs0(1:2)) 
    allocate(xzs0(1:2)) 
    allocate(yzs0(1:2)) 
    

    ! Just to be sure!
    s%zs = 0.0d0
    s%dzsdt = 0.0d0
    s%dzsdx = 0.0d0
    s%dzsdy = 0.0d0
    s%dzbdt = 0.0d0
    s%uu = 0.0d0
    s%vv = 0.0d0
    s%qx = 0.0d0
    s%qy = 0.0d0
    s%sedero = 0.0d0
    s%ueu = 0.0d0
    s%vev = 0.0d0
    s%vmagu = 0.0d0
    s%vmagv = 0.0d0
    s%vmageu = 0.0d0
    s%vmagev = 0.0d0
    s%u = 0.0d0
    s%v = 0.0d0
    s%ue = 0.0d0
    s%ve = 0.0d0
    s%hold = 0.0d0
    s%wetu = 0
    s%wetv = 0
    s%wetz = 0
    s%hh = 0.0d0
    s%hu = 0.0d0
    s%hv = 0.0d0
    s%hum = 0.0d0
    s%hvm = 0.0d0
    s%vu = 0.0d0
    s%uv = 0.0d0
    s%maxzs = 0.0d0
    s%minzs = 0.0d0
    s%taubx = 0.0d0
    s%tauby = 0.0d0
    s%ws = 0.0d0
    s%wb = 0.0d0 
    s%nuh = 0.0d0
    s%pres = 0.0d0 
    s%wi = 0.0d0 
    s%zi = 0.0d0 
    s%cf = 0.0d0
    s%zs0 = 0.0d0
    s%zs0fac = 0.0d0
    s%wm = 0.0d0
    s%umean = 0.0d0
    s%vmean = 0.0d0
    s%xyzs01 = 0.0d0
    s%xyzs02 = 0.0d0
    s%xyzs03 = 0.0d0
    s%xyzs04 = 0.0d0

    szs0 = 0.0d0
    xzs0 = 0.0d0
    yzs0 = 0.0d0

    ! TODO: do this properly....
    ! All variables above, should be initialized below (for all cells)
    s%zs0  = 0.0d0 
    s%ue   = 0.0d0
    s%ve   = 0.0d0
    s%ws   = 0.0d0
    s%wb   = 0.0d0
    s%pres = 0.0d0
    s%zi   = 0.0d0
    s%wi   = 0.0d0
    s%nuh  = 0.0d0
    s%cf   = par%cf
    s%wm   =0.d0
    s%umean=0.d0
    s%vmean=0.d0
    s%zs0fac = 0.0d0



    !
    ! set-up tide and surge waterlevels
    s%zs01=par%zs0
    ! I: read zs0 at model corners using zs0file
    if (par%tideloc>0) then

       s%zs01=s%tideinpz(1,1)

       if(par%tideloc.eq.1) s%zs02=s%zs01

       if(par%tideloc.eq.2 .and. trim(par%paulrevere)=='land') then
          s%zs03=s%tideinpz(1,2)
          s%zs02=s%zs01
          s%zs04=s%zs03
       endif

       if(par%tideloc.eq.2 .and. trim(par%paulrevere)=='sea') then
          s%zs02=s%tideinpz(1,2)
          s%zs03=0.d0
          s%zs04=0.d0
       endif

       if(par%tideloc.eq.4) then
          s%zs02=s%tideinpz(1,2)
          s%zs03=s%tideinpz(1,3)
          s%zs04=s%tideinpz(1,4)
       endif

       ! Set global domain corners for MPI simulations
       s%xyzs01(1) = s%sdist(1,1)          !x(1,1)
       s%xyzs01(2) = s%ndist(1,1)          !y(1,1)
       s%xyzs02(1) = s%sdist(1,s%ny+1)       !x(1,s%ny+1)
       s%xyzs02(2) = s%ndist(1,s%ny+1)       !y(1,s%ny+1)
       s%xyzs03(1) = s%sdist(s%nx+1,s%ny+1)    !x(s%nx+1,s%ny+1)
       s%xyzs03(2) = s%ndist(s%nx+1,s%ny+1)    !y(s%nx+1,s%ny+1)
       s%xyzs04(1) = s%sdist(s%nx+1,1)       !x(s%nx+1,1)
       s%xyzs04(2) = s%ndist(s%nx+1,1)               !y(s%nx+1,1) 
            
       !
       ! Fill in matrix zs0
       !
       if(par%tideloc.eq.1) s%zs0 = s%xz*0.0d0 + s%zs01

       if(par%tideloc.eq.2 .and. trim(par%paulrevere)=='sea') then
          yzs0(1)=s%ndist(1,1)
          yzs0(2)=s%ndist(1,s%ny+1)
          szs0(1)=s%zs01
          szs0(2)=s%zs02

          do j = 1,s%ny+1
             call LINEAR_INTERP(yzs0, szs0, 2, s%ndist(1,j), s%zs0(1,j), indt)
          enddo
          do j = 1,s%ny+1 
             do i = 1,s%nx+1
                s%zs0(i,j) = s%zs0(1,j)
             enddo
          enddo
       endif

       if(par%tideloc.eq.2 .and. trim(par%paulrevere)=='land') then
          yzs0(1)=s%sdist(1,1)
          yzs0(2)=s%sdist(s%nx+1,1)
          szs0(1)=s%zs01
          szs0(2)=s%zs04
          s%zs0(1,:)=s%zs01
          s%zs0(s%nx+1,:)=s%zs03
          do j = 1,s%ny+1 
             do i = 1,s%nx+1
                if (s%zb(i,j).gt.s%zs0(1,j)+par%eps) goto 302
                s%zs0(i,j) = s%zs0(1,j)
             enddo

302          if (i.lt.s%nx+1) then
                do ig=i+1,s%nx+1
                   s%zs0(ig,j) = s%zs0(s%nx+1,j) 
                enddo
             else
                do i = 1,s%nx+1
                  yzs0(1)=s%sdist(1,j)
                  yzs0(2)=s%sdist(s%nx+1,j)
                  call LINEAR_INTERP(yzs0, szs0, 2, s%sdist(i,j), s%zs0(i,j), indt)
                enddo
             endif
          enddo
       endif

       if(par%tideloc.eq.4) then
          szs0(1)=s%zs01
          szs0(2)=s%zs02
          do j = 1,s%ny+1
             yzs0(1)=s%ndist(1,1)
             yzs0(2)=s%ndist(1,s%ny+1)
             call LINEAR_INTERP(yzs0, szs0, 2, s%ndist(1,j), s%zs0(1,j), indt)
          enddo
          szs0(1)=s%zs04
          szs0(2)=s%zs03
          do j = 1,s%ny+1
             yzs0(1)=s%ndist(s%nx+1,1)
             yzs0(2)=s%ndist(s%nx+1,s%ny+1)
             call LINEAR_INTERP(yzs0, szs0, 2, s%ndist(s%nx+1,j), s%zs0(s%nx+1,j), indt)
          enddo

          do j = 1,s%ny+1 
             do i = 1,s%nx+1
                if (s%zb(i,j).gt.s%zs0(1,j)+par%eps) goto 303
                s%zs0(i,j) = s%zs0(1,j)
             enddo

303          if (i.lt.s%nx+1) then
                do ig=i+1,s%nx+1
                   s%zs0(ig,j) = s%zs0(s%nx+1,j) 
                enddo
             else
                do i = 1,s%nx+1
				   xzs0(1)=s%sdist(1,j)
				   xzs0(2)=s%sdist(s%nx+1,j)
				   szs0(1)=s%zs0(1,j)
				   szs0(2)=s%zs0(s%nx+1,j)
                   call LINEAR_INTERP(xzs0, szs0, 2, s%sdist(i,j), s%zs0(i,j), indt)
                enddo
             endif
          enddo
       endif
    else
       s%zs0 = s%zs01
    endif
     
    inquire(file=par%zsinitfile,exist=exists)
    if (exists) then
       open(723,file=par%zsinitfile)
       do j=1,s%ny+1
          read(723,*)(s%zs0(i,j),i=1,s%nx+1)
       enddo
       close(723)
    endif
    
    inquire(file=par%bedfricfile,exist=exists)
    if ((exists) .and. trim(par%bedfriction)=='chezy') then
       open(723,file=par%bedfricfile)
       do j=1,s%ny+1
          read(723,*)(s%cf(i,j),i=1,s%nx+1)
       enddo
       close(723)
       ! convert from C to cf
       s%cf=par%g/s%cf**2
    endif
    !
    ! set zs, hh, wetu, wetv, wetz
    !
    
    s%zs0 = max(s%zs0,s%zb)
    ! cjaap: replaced par%hmin by par%eps
    s%hh=max(s%zs0-s%zb,par%eps)
    s%zs=0.d0
    s%zs=max(s%zb,s%zs0)

    !Initialize hu correctly to prevent spurious initial flow (Pieter)
    do j=1,s%ny+1
       do i=1,s%nx
          s%hu(i,j) = max(s%zs(i,j),s%zs(i+1,j))-max(s%zb(i,j),s%zb(i+1,j))
       enddo
    enddo
    s%hu(s%nx+1,:)=s%hu(s%nx,:)

    if (s%ny>0) then
        do j=1,s%ny
           do i=1,s%nx+1
              s%hv(i,j) = max(s%zs(i,j),s%zs(i,j+1))-max(s%zb(i,j),s%zb(i,j+1))
           enddo
        enddo
        s%hv(:,s%ny+1)=s%hv(:,s%ny)
    else
        s%hv(:,1)=s%zs(:,1)-s%zb(:,1)
    endif
    
    s%hum(1:s%nx,:) = 0.5d0*(s%hh(1:s%nx,:)+s%hh(2:s%nx+1,:))
    s%hum(s%nx+1,:)=s%hh(s%nx+1,:)

    s%dzsdt=0.d0
    s%dzsdx=0.d0
    s%dzsdy=0.d0
    s%dzbdt=0.d0
    s%uu=0.d0
    s%u=0.d0
    s%vv=0.d0
    s%v=0.d0
    s%vu=0.d0
    s%uv=0.d0
    s%ueu=0.d0
    s%vev=0.d0
    s%qx=0.d0
    s%qy=0.d0
    s%sedero=0.d0
    s%vmagu=0.d0
    s%vmagv=0.d0
    s%vmageu=0.d0
    s%vmagev=0.d0
    s%taubx=0.d0
    s%tauby=0.d0
    s%maxzs=-999.d0
    s%minzs=999.d0
    where(s%zs>s%zb+par%eps)
       s%wetz=1
       s%wetu=1
       s%wetv=1
    elsewhere
       s%wetz=0
       s%wetu=0
       s%wetv=0
    endwhere
    !
    if (trim(par%tidetype)=='instant') then
        ! RJ: 22-09-2010
        ! Check for whole domain whether a grid cell should be associated with
        ! 1) offshore tide and surge
        ! 2) bay tide and surge
        ! 3) no tide and surge
        ! 4) weighted tide and surge (for completely wet arrays)
        ! relative weight of offshore boundary and bay boundary for each grid point is stored in zs0fac 
        !
        s%zs0fac = 0.d0
        do j = 1,s%ny+1 
           offshoreregime = .true.   
           indoff = s%nx+1 ! ind of last point (starting at offshore boundary) that should be associated with offshore boundary
           indbay = 1      ! ind of first point (starting at offshore boundary) that should be associated with bay boundary
           do i = 1,s%nx+1
              if (offshoreregime .and. s%wetz(i,j)==0) then
                 indoff = max(i-1,1)
                 offshoreregime = .false.
              endif
              if (s%wetz(i,j)==0 .and. s%wetz(min(i+1,s%nx+1),j)==1) then
                 indbay = min(i+1,s%nx+1)
              endif
           enddo
              
           if (indbay==1 .and. indoff==s%nx+1) then ! in case of completely wet arrays linear interpolation for zs0fac
! Dano: don't know how to fix this for curvilinear
!          zs0fac(:,j,2) = (xz-xz(1))/(xz(s%nx+1)-xz(1))
!          zs0fac(:,j,1) = 1-zs0fac(:,j,2)
           else                                    ! in all other cases we assume three regims offshore, dry and bay
              s%zs0fac(1:indoff,j,1) = 1.d0
              s%zs0fac(1:indoff,j,2) = 0.d0
              if (indbay > 1) then
                 s%zs0fac(indoff+1:indbay-1,j,1) = 0.d0
                 s%zs0fac(indbay:s%nx+1,j,1) = 0.d0
                 s%zs0fac(indoff+1:indbay-1,j,2) = 0.d0
                 s%zs0fac(indbay:s%nx+1,j,2) = 1.d0
              endif   
           endif       
        enddo
 
    endif ! tidetype = instant water level boundary
      
  end subroutine flow_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine sed_init (s,par)
    use params
    use spaceparams
    use readkey_module
    use xmpi_module
    use logging_module

    IMPLICIT NONE

    type(spacepars)                     :: s
    type(parameters)                    :: par

    integer                             :: i,j,m,jg,start
    character*80                        :: fnameg,line
    character*4                         :: tempc
    real*8                              :: tempr

    include 's.ind'
    include 's.inp'

    allocate(s%ccg(1:nx+1,1:ny+1,par%ngd))
    allocate(s%ccbg(1:nx+1,1:ny+1,par%ngd))
    allocate(s%dcbdy(1:nx+1,1:ny+1))
    allocate(s%dcbdx(1:nx+1,1:ny+1))
    allocate(s%dcsdy(1:nx+1,1:ny+1))
    allocate(s%dcsdx(1:nx+1,1:ny+1))
    allocate(s%Tsg(1:nx+1,1:ny+1,par%ngd)) 
    allocate(s%Susg(1:nx+1,1:ny+1,par%ngd)) 
    allocate(s%Svsg(1:nx+1,1:ny+1,par%ngd)) 
    allocate(s%Subg(1:nx+1,1:ny+1,par%ngd)) 
    allocate(s%Svbg(1:nx+1,1:ny+1,par%ngd)) 
    allocate(s%vmag(1:nx+1,1:ny+1)) 
    allocate(s%ceqsg(1:nx+1,1:ny+1,par%ngd))
    allocate(s%ceqbg(1:nx+1,1:ny+1,par%ngd))
    allocate(s%D50(1:par%ngd))
    allocate(s%D90(1:par%ngd))
    allocate(s%D50top(1:nx+1,1:ny+1))
    allocate(s%D90top(1:nx+1,1:ny+1))
    allocate(s%sedcal(1:par%ngd))
    allocate(s%ucrcal(1:par%ngd))
    allocate(s%nd(1:nx+1,1:ny+1))
    allocate(s%dzbed(1:nx+1,1:ny+1,1:max(par%nd,3))) 
    allocate(s%pbbed(1:nx+1,1:ny+1,1:max(par%nd,3),1:par%ngd)) 
    allocate(s%z0bed(1:nx+1,1:ny+1))
    allocate(s%ureps(1:nx+1,1:ny+1))
    allocate(s%urepb(1:nx+1,1:ny+1))
    allocate(s%vreps(1:nx+1,1:ny+1))
    allocate(s%vrepb(1:nx+1,1:ny+1))
    allocate(s%ero(1:nx+1,1:ny+1,1:par%ngd))
    allocate(s%depo_ex(1:nx+1,1:ny+1,1:par%ngd))
    allocate(s%depo_im(1:nx+1,1:ny+1,1:par%ngd))
    allocate(s%kb(1:nx+1,1:ny+1))
    allocate(s%Tbore(1:nx+1,1:ny+1))
    allocate(s%ua(1:nx+1,1:ny+1))  
    allocate(s%dzav(1:nx+1,1:ny+1))  
    allocate(s%Sk(1:nx+1,1:ny+1))
    allocate(s%As(1:nx+1,1:ny+1))
    allocate(s%kturb(1:nx+1,1:ny+1))
    allocate(s%rolthick(1:nx+1,1:ny+1))
    allocate(s%Sutot(1:nx+1,1:ny+1))     ! Only really for easy output 
    allocate(s%Svtot(1:nx+1,1:ny+1))     ! Only really for easy output
    allocate(s%cctot(1:nx+1,1:ny+1))     ! Only really for easy output

    ! Initialize so structures can be implemented more easily
    s%pbbed = 0.d0


    !
    ! Set grain size(s)
    !
    call readkey('params.txt','D50',line)
    if (line=='') then
       s%D50=0.0002d0   ! Default
    else
       read(line,*) s%D50(1:par%ngd)
    endif
    if (maxval(s%D50)>0.002) call writelog('ls','', 'D50 > 2mm, out of validity range')
    call readkey('params.txt','D90',line)
    if (line=='') then
       s%D90=0.0003d0   ! Default
    else
       read(line,*) s%D90(1:par%ngd)
    endif
    call readkey('params.txt','sedcal',line)
    if (line=='') then
       s%sedcal=1.d0    ! Default
    else
       read(line,*) s%sedcal(1:par%ngd)
    endif
    call readkey('params.txt','ucrcal',line)
    if (line=='') then
       s%ucrcal=1.d0    ! Default
    else
       read(line,*) s%ucrcal(1:par%ngd)
    endif

    if (par%ngd==1) then

       ! No multi sediment, but we do need some data to keep the script running

       s%pbbed(:,:,:,1)=1.d0   ! set sand fraction everywhere, not structure fraction (if exist) which is still 0.d0
       par%nd_var=2

       s%dzbed(:,:,1:par%nd_var-1)       = max(par%dzg1,10.d0)
       s%dzbed(:,:,par%nd_var)           = max(par%dzg2,10.d0)
       s%dzbed(:,:,par%nd_var+1:par%nd)  = max(par%dzg3,10.d0)

    else

       ! Fill pbed en dzbed
       ! 
       do jg=1,par%ngd
          write(tempc,'(i4)')jg
          start=4-floor(log10(real(jg)))
          write(fnameg,'(a,a,a)')'gdist',tempc(start:4),'.inp'
          open(31,file=fnameg)
          do m=1,par%nd
             do j=1,ny+1
                read(31,*)(s%pbbed(i,j,m,jg),i=1,nx+1)
             enddo
          enddo
          close(31)
       enddo
       ! Rework pbbed so that sum fractions = 1
       do m=1,par%nd
          do j=1,ny+1     !Jaap instead of 2:ny
             do i=1,nx+1 !Jaap instead of 2:nx

                tempr=sum(s%pbbed(i,j,m,1:par%ngd))
                if (abs(1.d0-tempr)>0.d0) then
                   ! Maybe fix this warning if in combination with structures
                   call writelog('ls','ai0ai0ai0a',' Warning: Resetting sum of sediment fractions in point (',&
                        i,',',j,') layer ,',m,&
                        ' to equal unity.')
                   if (tempr<=tiny(0.d0)) then    ! In case cell has zero sediment (i.e. only hard structure)
                      s%pbbed(i,j,m,:)=1.d0/dble(par%ngd) 
                   else
                      s%pbbed(i,j,m,:)=s%pbbed(i,j,m,:)/tempr
                   endif
                endif
             enddo
          enddo
       enddo
       ! boundary neumann --> Jaap not necessary already done in loop above

       ! sediment thickness				   
       s%dzbed(:,:,1:par%nd_var-1)       = par%dzg1
       s%dzbed(:,:,par%nd_var)           = par%dzg2
       s%dzbed(:,:,par%nd_var+1:par%nd)  = par%dzg3
    endif

    ! Initialize representative sed.diameter at the bed for flow friction and output
    do j=1,ny+1
       do i=1,nx+1
          s%D50top(i,j) =  sum(s%pbbed(i,j,1,:)*s%D50)
          s%D90top(i,j) =  sum(s%pbbed(i,j,1,:)*s%D90)
       enddo
    enddo
    ! 
    ! Set non-erodable layer
    !
    allocate(s%structdepth(nx+1,ny+1))

    s%structdepth = 100.d0

    if (par%struct==1) then
       !call readkey('params.txt','ne_layer',fnameh)
       !open(31,file=fnameh)
       open(31,file=par%ne_layer)

       do j=1,ny+1
          read(31,*)(s%structdepth(i,j),i=1,nx+1)
       end do

       close(31)

    endif

    ! bottom of sediment model
    s%z0bed = s%zb - sum(s%dzbed,DIM=3)

    s%nd = max(par%nd,2)

    s%ureps      = 0.d0
    s%vreps      = 0.d0
    s%ccg        = 0.d0
    s%ccbg       = 0.d0
    s%ceqbg      = 0.d0
    s%ceqsg      = 0.d0
    s%Susg       = 0.d0
    s%Svsg       = 0.d0
    s%Subg       = 0.d0
    s%Svbg       = 0.d0
    s%dcsdx      = 0.d0
    s%dcsdy      = 0.d0
    s%dcbdx      = 0.d0
    s%dcbdy      = 0.d0
    s%ero        = 0.d0
    s%depo_im    = 0.d0
    s%depo_ex    = 0.d0
    s%kb         = 0.d0
    s%Tbore      = 0.d0
    s%ua         = 0.d0
    s%dzav       = 0.d0
    s%Sk         = 0.d0
    s%As         = 0.d0
    s%kturb      = 0.d0
    s%rolthick   = 0.d0
    s%Sutot      = 0.d0
    s%Svtot      = 0.d0
    s%cctot      = 0.d0

    ! Initialize dzbdx, dzbdy
    do j=1,ny+1
       do i=1,nx
          s%dzbdx(i,j)=(s%zb(i+1,j)-s%zb(i,j))/dsu(i,j)
       enddo
    enddo
    ! dummy, needed to keep compiler happy
    s%dzbdx(nx+1,:)=s%dzbdx(nx,:)

    if (ny>0) then
       do j=1,ny
          do i=1,nx+1
             s%dzbdy(i,j)=(s%zb(i,j+1)-s%zb(i,j))/dnv(i,j)
          enddo
       enddo
       s%dzbdy(:,ny+1)=s%dzbdy(:,ny)
    else
	   s%dzbdy=0.d0
	endif

  end subroutine sed_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine discharge_init(s, par)
  
    use params
    use spaceparams
    use readkey_module
    use xmpi_module
    use logging_module

    implicit none

    type(spacepars)                         :: s
    type(parameters)                        :: par

    integer                                 :: i,j
    integer                                 :: io
    integer                                 :: m1,m2,n1,n2
    real*8                                  :: dxd,dyd
    real*8, dimension(:),allocatable        :: xdb,ydb,xde,yde
    integer,dimension(2)                    :: mnb,mne

    include 's.ind'
    include 's.inp'
    
    io              = 0

    allocate(xdb    (par%ndischarge)      )
    allocate(ydb    (par%ndischarge)      )
    allocate(xde    (par%ndischarge)      )
    allocate(yde    (par%ndischarge)      )
    
    allocate(s%pntdisch     (1:par%ndischarge)                      )
    allocate(s%pdisch       (1:par%ndischarge   , 1:4)              )
    allocate(s%tdisch       (1:par%ntdischarge)                     )
    allocate(s%qdisch       (1:par%ntdischarge  , 1:par%ndischarge) )
    
    s%pntdisch  = 0.d0
    s%tdisch    = 0.d0
    s%pdisch    = 0.d0
    s%qdisch    = 0.d0
    
    if (xmaster) then
        if (par%ndischarge>0) then
        
            ! read discharge locations
            open(10,file=par%disch_loc_file)
            do i=1,par%ndischarge
                read(10,*,IOSTAT=io) xdb(i),ydb(i),xde(i),yde(i)
                
                ! distinguish between horizontal and vertical discharge
                if (xdb(i).eq.xde(i) .and. ydb(i).eq.yde(i)) then
                    s%pntdisch(i) = 1
                else
                    s%pntdisch(i) = 0
                endif
                
            enddo
            close(10)
        
            if (par%ntdischarge>0) then
            
                ! read time series
                open(10,file=par%disch_timeseries_file)
                do i=1,par%ntdischarge
                    read(10,*,IOSTAT=io) s%tdisch(i),(s%qdisch(i,j),j=1,par%ndischarge)
                enddo
                close(10)
            endif
        endif
    endif
    
    if (xmaster) then
    
        ! initialise each discharge location
        do i=1,par%ndischarge
            
            dxd = abs(xde(i)-xdb(i))
            dyd = abs(yde(i)-ydb(i))
            
            ! convert discharge location to cell indices depending on type of discharge:
            !     point discharge, in v-direction or in u-direction
            if (s%pntdisch(i).eq.1) then
            
                ! point discharge (no orientation, no added momentum, just mass)
                
                mnb = minloc(sqrt((s%xz-xdb(i))**2+(s%yz-ydb(i))**2))
                mne = minloc(sqrt((s%xz-xde(i))**2+(s%yz-yde(i))**2))
                
                s%pdisch(i,:) = (/mnb(1),mnb(2),0,0/)
                
            elseif (dxd.gt.dyd) then
            
                ! discharge through v-points
                
                mnb = minloc(sqrt((s%xv-xdb(i))**2+(s%yv-ydb(i))**2))
                mne = minloc(sqrt((s%xv-xde(i))**2+(s%yv-yde(i))**2))
            
                m1 = minval((/mnb(1),mne(1)/))
                m2 = maxval((/mnb(1),mne(1)/))
                n1 = nint(0.5*(mnb(2)+mne(2)))
                
                if (n1.lt.1)    n1 = 1
                if (n1.gt.s%ny) n1 = s%ny
                
                s%pdisch(i,:) = (/m1,n1,m2,n1/)
            else
            
                ! discharge through u-points
                
                mnb = minloc(sqrt((s%xu-xdb(i))**2+(s%yu-ydb(i))**2))
                mne = minloc(sqrt((s%xu-xde(i))**2+(s%yu-yde(i))**2))
                
                m1 = nint(0.5*(mnb(1)+mne(1)))
                n1 = minval((/mnb(2),mne(2)/))
                n2 = maxval((/mnb(2),mne(2)/))
                
                if (m1.lt.1)    m1 = 1
                if (m1.gt.s%nx) m1 = s%nx
                
                s%pdisch(i,:) = (/m1,n1,m1,n2/)
            endif
        enddo
    endif
  end subroutine discharge_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine drifter_init(s, par)
  
    use params
    use spaceparams
    use readkey_module
    use xmpi_module
    use logging_module

    implicit none

    type(spacepars)                         :: s
    type(parameters)                        :: par
    
    character(256)                          :: drifterfile
    integer                                 :: i
    real*8                                  :: xdrift,ydrift
    real*8                                  :: ds,dn
    integer,dimension(2)                    :: mn

    include 's.ind'
    include 's.inp'
    
    if (par%ndrifter>0) then
   
        if (xmaster) then
            
            allocate(s%idrift   (par%ndrifter))
            allocate(s%jdrift   (par%ndrifter))
            allocate(s%tdriftb  (par%ndrifter))
            allocate(s%tdrifte  (par%ndrifter))
        
            ! read drifter file
            drifterfile = readkey_name('params.txt','drifterfile',bcast=.false.)
            open(10,file=drifterfile)
            do i=1,par%ndrifter
                read(10,*)xdrift,ydrift,s%tdriftb(i),s%tdrifte(i)
                
                mn          = minloc(sqrt((s%xz-xdrift)**2+(s%yz-ydrift)**2))
                
                ds          =  (xdrift - s%xz(mn(1),mn(2)))*cos(s%alfaz(mn(1),mn(2))) &
                              +(ydrift - s%yz(mn(1),mn(2)))*sin(s%alfaz(mn(1),mn(2)))
                dn          = -(xdrift - s%xz(mn(1),mn(2)))*sin(s%alfaz(mn(1),mn(2))) &
                              +(ydrift - s%yz(mn(1),mn(2)))*cos(s%alfaz(mn(1),mn(2)))
                       
                s%idrift(i) = mn(1) + ds/dsu(mn(1),mn(2))
                s%jdrift(i) = mn(2) + dn/dnv(mn(1),mn(2))
            enddo
            close(10)
        endif
    endif
  end subroutine drifter_init

end module initialize
