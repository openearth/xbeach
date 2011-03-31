!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   MODULE POSTPROCESS    !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This module contains transformation routines from xbeach internal data structures to output data structures.

module postprocessmod
    implicit none
  interface gridrotate
     ! rotate grids, given s, t and outputs x
     module procedure gridrotate_r0
     module procedure gridrotate_r1
     module procedure gridrotate_r2
     module procedure gridrotate_r3
     module procedure gridrotate_r4
     module procedure gridrotate_i0
     module procedure gridrotate_i1
     module procedure gridrotate_i2
     module procedure gridrotate_i3
     module procedure gridrotate_i4
  end interface gridrotate

contains
  subroutine snappointstogrid(par, s, xpoints, ypoints,scrprintinp,dmin)
    ! Lookup the nearest neighbour grid coordinates of the output points specified in params file
    ! Convert world coordinates of points to nearest (lsm) grid point
    use spaceparams
    use params
    use logging_module

    type(spacepars), intent(in)      :: s
    type(parameters), intent(in)     :: par

    integer*4,dimension(:),intent(inout) :: xpoints     ! model x-coordinate of output points
    integer*4,dimension(:),intent(inout) :: ypoints     ! model y-coordinate of output points

    real*8,dimension(s%nx+1,s%ny+1)     :: mindist
    integer,dimension(2)                :: minlocation
    integer                             :: i
    logical,optional,intent(in)         :: scrprintinp
    logical                             :: scrprint
    real*8,dimension(:),allocatable     :: mindistr
    real*8,dimension(:),optional,intent(out)  :: dmin
    !
    if (present(scrprintinp)) then
       scrprint=scrprintinp
    else
       scrprint=.true.
    endif
    !
    allocate(mindistr(par%npoints+par%nrugauge))
    ! Let's hope that the s%xz and s%yz are already available....
    ! Compute the minimum distances for each point
    if (par%npoints + par%nrugauge > 0) then
       do i=1,(par%npoints+par%nrugauge)
          mindist=sqrt((par%xpointsw(i)-s%xz)**2+(par%ypointsw(i)-s%yz)**2)
          ! look up the location of the found minimum
          minlocation=minloc(mindist)
          ! minimum distance
          mindistr(i) = mindist(minlocation(1),minlocation(2))
          ! The y coordinate is always the same
          ypoints(i)=minlocation(2)
          ! For rugauges the xpoint is always 1
          if (par%pointtypes(i) == 1) then
             xpoints(i)=1
             if (scrprint) call writelog('ls','(a,i0)','Runup gauge at grid line iy=',ypoints(i))
          else
             xpoints(i)=minlocation(1)
             if (scrprint) call writelog('ls','(a,i0,a,i0,a,f0.2,a)',' Distance output point to nearest grid point ('&
                                                                    ,minlocation(1),',',minlocation(2),') is '&
                                                                    ,mindistr(i), ' meters')
          endif
       end do
    end if
    if (present(dmin)) then
       dmin = mindistr
    endif
    deallocate(mindistr)
  end subroutine snappointstogrid


  subroutine gridrotate_r0(s, t, x)
    use spaceparams
    use mnemmodule

    type(spacepars), intent(in) :: s
    type(arraytype), intent(in) :: t
    real*8                      :: x
    x = t%r0
  end subroutine gridrotate_r0

  subroutine gridrotate_r1(s, t, x)
    use spaceparams
    use mnemmodule

    type(spacepars), intent(in)   :: s
    type(arraytype), intent(in)   :: t
    real*8, dimension(:)          :: x
    real*8                        :: pi
    pi = 4*atan(1.0d0)
    select case(t%name)
    case(mnem_theta)
       x=270-(t%r1*(180/pi))
    case(mnem_theta0)
       x=270-(t%r1*(180/pi))
    case default
       x=t%r1
    end select
  end subroutine gridrotate_r1


  subroutine gridrotate_r2(s, t, x)
    use spaceparams
    use mnemmodule

    type(spacepars), intent(in)   :: s
    type(arraytype), intent(in)   :: t
    real*8, dimension(:,:), intent(out) :: x
    real*8                        :: pi
    pi = 4*atan(1.0d0)
    select case(t%name)
    case(mnem_thetamean)
       x=270-((t%r2+s%alfaz)*(180/pi))
    case(mnem_Fx)
       x=t%r2*cos(s%alfaz)-s%Fy*sin(s%alfaz)
    case(mnem_Fy)
       x=s%Fx*sin(s%alfaz)+t%r2*cos(s%alfaz)
    case(mnem_u)
       x=t%r2*cos(s%alfaz)-s%v*sin(s%alfaz)
    case(mnem_gwu)
       x=t%r2*cos(s%alfaz)-s%gwv*sin(s%alfaz)
    case(mnem_v)
       x=s%u*sin(s%alfaz)+t%r2*cos(s%alfaz)
    case(mnem_gwv)
       x=s%gwu*sin(s%alfaz)+t%r2*cos(s%alfaz)
    case(mnem_ue)
       x=t%r2*cos(s%alfaz)-s%ve*sin(s%alfaz)
    case(mnem_ve)
       x=s%ue*sin(s%alfaz)+t%r2*cos(s%alfaz)
    case(mnem_uwf)
       x=t%r2*cos(s%alfaz)-s%vwf*sin(s%alfaz)
    case(mnem_vwf)
       x=s%uwf*sin(s%alfaz)+t%r2*cos(s%alfaz)
    case(mnem_Sutot)
       x=(sum(s%Subg,dim=3)+sum(s%Susg,dim=3))*cos(s%alfaz) - (sum(s%Svbg,dim=3)+sum(s%Svsg,dim=3))*sin(s%alfaz)
    case(mnem_Svtot)
       x=(sum(s%Subg,dim=3)+sum(s%Susg,dim=3))*sin(s%alfaz) + (sum(s%Svbg,dim=3)+sum(s%Svsg,dim=3))*cos(s%alfaz)
    case(mnem_cctot)
       x=sum(s%ccg,dim=3)
    case default
       x=t%r2
    end select

  end subroutine gridrotate_r2

  subroutine gridrotate_r3(s, t, x)
    use spaceparams
    use mnemmodule
    use logging_module

    implicit none
    
    type(spacepars), intent(in)   :: s
    type(arraytype), intent(in)   :: t
    real*8, dimension(:,:,:)      :: x
    ! no need to allocate
    real*8, dimension(size(s%alfaz,1), size(s%alfaz,2), size(t%r3,3)) :: alfazr3
    real*8                        :: pi
    integer                       :: i
    pi = 4*atan(1.0d0)
    
    ! fill variable alfazr3. We presume first 2 dimension are related to nx+1, ny+1 respectively
    ! Better double check    
    if (size(s%alfaz,1) .ne. size(t%r3,1)) call writelog('ls', '', &
           & 'Assertion error, s%alfaz and t%r3 do not align on 1st dimension', size(s%alfaz,1), size(t%r3,1))
    if (size(s%alfaz,2) .ne. size(t%r3,2)) call writelog('ls', '', & 
           & 'Assertion error, s%alfaz and t%r3 do not align on 2nd dimension', size(s%alfaz,2), size(t%r3,2))
    ! This should be something like:
    ! alfazr3 = (/(s%alfaz, i=1,size(t%r3,3)) /)
    do i=1,size(t%r3,3)
        alfazr3(:,:,i) = s%alfaz
    end do  
    
    select case(t%name)
    case(mnem_cgx)
       x=t%r3*cos(alfazr3)-s%cgy*sin(alfazr3)
    case(mnem_cgy)
       x=s%cgx*sin(alfazr3)+t%r3*cos(alfazr3)
    case(mnem_cx)
       x=t%r3*cos(alfazr3)-s%cy*sin(alfazr3)
    case(mnem_cy)
       x=s%cx*sin(alfazr3)+t%r3*cos(alfazr3)
    case(mnem_thet)
       x=270-((s%thet+alfazr3)*(180/pi))
    case(mnem_Susg)
       x=t%r3*cos(alfazr3)-s%Svsg*sin(alfazr3)
    case(mnem_Svsg)
       x=s%Susg*sin(alfazr3)+t%r3*cos(alfazr3)
    case(mnem_Subg)
       x=t%r3*cos(alfazr3)-s%Svbg*sin(alfazr3)
    case(mnem_Svbg)
       x=s%Subg*sin(alfazr3)+t%r3*cos(alfazr3)
    case default
       x=t%r3
    end select
  end subroutine gridrotate_r3

  subroutine gridrotate_r4(s, t, x)
    use spaceparams
    use mnemmodule

    type(spacepars), intent(in)   :: s
    type(arraytype), intent(in)   :: t
    real*8, dimension(:,:,:,:)    :: x
    x = t%r4
  end subroutine gridrotate_r4


  subroutine gridrotate_i0(s, t, x)
    use spaceparams
    use mnemmodule

    type(spacepars), intent(in) :: s
    type(arraytype), intent(in) :: t
    integer                     :: x
    x = t%i0
  end subroutine gridrotate_i0

  subroutine gridrotate_i1(s, t, x)
    use spaceparams
    use mnemmodule

    type(spacepars), intent(in) :: s
    type(arraytype), intent(in) :: t
    integer, dimension(:),intent(out) :: x
    x = t%i1
  end subroutine gridrotate_i1

  subroutine gridrotate_i2(s, t, x)
    use spaceparams
    use mnemmodule

    type(spacepars), intent(in) :: s
    type(arraytype), intent(in) :: t
    integer, dimension(:,:),intent(out) :: x
    x = t%i2
  end subroutine gridrotate_i2
  subroutine gridrotate_i3(s, t, x)
    use spaceparams
    use mnemmodule

    type(spacepars), intent(in) :: s
    type(arraytype), intent(in) :: t
    integer, dimension(:,:,:),intent(out) :: x
    x = t%i3
  end subroutine gridrotate_i3
  subroutine gridrotate_i4(s, t, x)
    use spaceparams
    use mnemmodule

    type(spacepars), intent(in) :: s
    type(arraytype), intent(in) :: t
    integer, dimension(:,:,:,:),intent(out) :: x
    x = t%i4
  end subroutine gridrotate_i4


  function runup(par, s)
    use spaceparams
    use params
    type(spacepars) :: s
    type(parameters)    :: par
    integer, dimension(par%npoints + par%nrugauge) :: runup

    ! This is not correct....
    where (par%pointtypes == 1)
       runup = max(maxval(minloc(s%wetz(:,:)))-1,1)
    end where
  end function runup

end module postprocessmod
