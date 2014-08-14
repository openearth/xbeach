   module wave_bc_nextgen
      ! New setup to handle wave boundary conditions in XBeach, that can also be ported to 
      ! Delft3D / FlowFM   
      ! The following modules from XBeach need to be included in other models to allow this
      ! module to work
      use paramsconst
      implicit none
      private 
      !public wave_bc_generate
      
      logical      :: initialized = .false.
      
   contains
!   
!   subroutine wave_bc_generate(t,wbctype,xref, yref, ntheta,dtheta,thetamin, theta, seed,  &
!                                nbc, xbc, ybc, hbc, isBoundary,                     &
!                                eebc, uibc, vibc, zsbc, isRecomputed )
!      ! This subroutine is the main (only) call from the external program (XBeach/Delft3D)
!      ! The subroutine generates one value of wave energy and/or wave flux per offshore boundary 
!      ! point, depending on the wbctype specified in the input of this call, for time = t. If 
!      ! necessary, this subroutine will initialise variables and/or generate new time series
!      ! files.
!      !
!      !
!      ! Input variables into subroutine are:
!      !
!      ! t             = current time (s)
!      ! wbctype       = switch for type of wave boundary condition (==par%instat in XBeach)
!      ! xref          = world x-coordinate of reference point from which waves are unrolled 
!      !                 (in case of XBeach, this is xw(1,1))
!      ! yref          = world y-coordinate of reference point from which waves are unrolled 
!      !                 (in case of XBeach, this is yw(1,1))
!      ! 
!      ! 
!      !   
!      ! Output variables of the subroutine are:
!      !
!      !
!      ! Input variables
!      real*8, intent(in)                        :: t,xref,yref,dtheta,thetamin
!      integer, intent(in)                       :: wbctype,ntheta,nbc
!      real*8,dimension(nbc), intent(in)         :: xbc,ybc,hbc
!      integer,dimension(40),intent(in)          :: seed
!      logical,intent(in)                        :: isBoundary
!      ! Output variables
!      real*8,dimension(nbc,ntheta),intent(out)  :: eebc
!      real*8,dimension(nbc),intent(out)         :: uibc,vibc,zsbc
!      logical,intent(out)                       :: isRecomputed
!      ! Internal variables
!      integer                                   :: i,j,itheta
!      
!   
!      ! First check if this is a subdomain with a wave boundary
!      if (isBoundary .eqv. false) then
!         ! This domain does not have to do any work
!         eebc = 0.d0
!         uibc = 0.d0
!         vibc = 0.d0
!         zsbc = 0.d0
!         isRecomputed = .false.
!      else
!         
!      endif
!         
!         
!         
!      
!   
!
!!!!   New setup for wave_bc, spectral wave conditions (surfbeat mode only)
!!!! 
!!!!   Independent of XBeach where possible
!!!!   Parallel implementation: each process has to be independent
!!!!   Information common to all domains needed:
!!!!      x,y reference point
!!!!      mean depth along boundary
!!!!      ntheta,dtheta,thetamin, theta
!!!!      seed for random phases (is this enough?)
!!!!   Local information needed:
!!!!      nbc              - number of boundary points
!!!!      xbc, ybc         - locations of boundary points (preferably ordered)
!!!!      distbc           - cum. distance along boundary points (or computed within?)
!!!!      hbc              - depth at boundary points
!!!
!!!
!!!    integer, intent(in)             :: instat
!!!    real*8,  intent(in)             :: xref
!!!    real*8,  intent(in)             :: yref
!!!
!!!
!!!    startbcf=.false.
!!!
!!!    if(.not. bccreated ) then
!!!       bccreated=.true.
!!!       startbcf=.true.                     ! trigger read from bcf for instat 3,4,5,7
!!!       bcendtime=huge(0.0d0)               ! initial assumption for instat 3,4,5,7
!!!       s%newstatbc=1
!!!       
!!!       
!!!       if ((instat==INSTAT_JONS.or.instat==INSTAT_JONS_TABLE &
!!!         & .or. instat==INSTAT_SWAN.or.instat==INSTAT_VARDENS).and.xmaster) then
!!!          call spectral_wave_bc(sg,par,curline)
!!!       elseif (instat==INSTAT_REUSE.and.xmaster) then
!!!          curline = 1
!!!       endif
!!!       
!!!    end if
!!!
!!!    if (t .ge. bcendtime) then  ! Recalculate bcf-file 
!!!       close(71)
!!!       close(72)
!!!       if ((instat==INSTAT_JONS.or.instat==INSTAT_JONS_TABLE &
!!!         & .or. instat==INSTAT_SWAN.or.instat==INSTAT_VARDENS).and.xmaster) then
!!!          call spectral_wave_bc(sg,par,curline)
!!!       elseif (instat==INSTAT_REUSE.and.xmaster) then
!!!          startbcf=.true. 
!!!          curline = curline + 1
!!!       endif
!!!    end if
!!!    !
!!!    ! COMPUTE WAVE BOUNDARY CONDITIONS CURRENT TIMESTEP
!!!    if (  (instat==INSTAT_JONS).or. &
!!!          (instat==INSTAT_JONS_TABLE).or. &
!!!          (instat==INSTAT_SWAN) .or. &
!!!          (instat==INSTAT_VARDENS) .or. &
!!!          (instat==INSTAT_REUSE) ) then
!!!       ! open file if first time
!!!       if (startbcf) then
!!!
!!!          bcendtime=wbcseries(curline)%bcendtime
!!!          rt       =wbcseries(curline)%rt
!!!          dtbcfile =wbcseries(curline)%dtbcfile
!!!          Trep     =wbcseries(curline)%Trep
!!!          theta0   =wbcseries(curline)%theta0
!!!          ebcfname =wbcseries(curline)%ebcfname
!!!          qbcfname =wbcseries(curline)%qbcfname       
!!!       
!!!       
!!!       endif
!!!
!!!!Dano DOES THIS  BELONG HERE?       
!!!       do itheta=1,ntheta
!!!          sigt(:,:,itheta) = 2*par%px/par%Trep
!!!       end do
!!!       sigm = sum(sigt,3)/ntheta
!!!       call dispersion(par,s)     
!!!       ! End initialize
!!!!Dano END QUESTION       
!!!  
!!!       inquire(iolength=wordsize) 1.d0
!!!       reclen=wordsize*(sg%ny+1)*(sg%ntheta)
!!!       open(71,file=ebcfname,status='old',form='unformatted',access='direct',recl=reclen)
!!!       reclen=wordsize*((sg%ny+1)*4)
!!!       open(72,file=qbcfname,status='old',form='unformatted',access='direct',recl=reclen)
!!!       if (.not. allocated(q1) ) then
!!!          allocate(q1(ny+1,4),q2(ny+1,4),q(ny+1,4))
!!!          allocate(ee1(ny+1,ntheta),ee2(ny+1,ntheta))
!!!       end if
!!!       read(71,rec=1,iostat=ier )ee1      ! Earlier in time
!!!       read(71,rec=2,iostat=ier2)ee2      ! Later in time
!!!       read(72,rec=1,iostat=ier )q1       ! Earlier in time
!!!       read(72,rec=2,iostat=ier2)q2       ! Later in time
!!!       old=floor((par%t/dtbcfile)+1)
!!!       recpos=1
!!!    end if
!!!
!!!    new=floor((par%t/dtbcfile)+1)
!!!
!!!    ! Check for next level in boundary condition file
!!!    if (new/=old) then
!!!       recpos=recpos+(new-old)
!!!       ! Check for how many bcfile steps are jumped
!!!       if (new-old>1) then  ! Many steps further in the bc file
!!!          read(72,rec=recpos+1,iostat=ier)q2
!!!          read(71,rec=recpos+1,iostat=ier)ee2
!!!          read(72,rec=recpos  ,iostat=ier)q1
!!!          read(71,rec=recpos  ,iostat=ier)ee1
!!!       else  ! Only one step further in the bc file
!!!          ee1=ee2
!!!          q1=q2
!!!          read(72,rec=recpos+1,iostat=ier)q2
!!!          read(71,rec=recpos+1,iostat=ier)ee2
!!!          endif
!!!          old=new
!!!       end if
!!!    endif
!!!
!!!    tnew = dble(new)*dtbcfile
!!!    facinterp=(tnew-par%t)/dtbcfile
!!!    eebc = (1.d0-facinterp)*ee2 + facinterp*ee1
!!!    qbc  = (1.d0-facinterp)*q2  + facinterp*q1 
!!!    uibc = qbc/hbc*min(par%t/par%taper,1.0d0)
!!!    vibc = qbc/hbc*min(par%t/par%taper,1.0d0)
!!!    eebc = eebc*min(par%t/par%taper,1.0d0)
!!!
!!!
!      end subroutine wave_bc_generate
      
      end module wave_bc_nextgen
     