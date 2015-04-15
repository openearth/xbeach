!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! wwvv-todo many changes wrt trunk
!!!!!!!!!!!!!!!!!!!!!!!   MODULE POSTPROCESS    !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This module contains transformation routines from xbeach internal data structures to output data structures.

module postprocessmod
   implicit none
   save
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
      character(1000)                           :: txt
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
               if(scrprint) then  ! wwvv-todo
                  write(txt,"('Distance output point',i3.3,'(',f0.2,',',f0.2,') to gridpoint(',i0,',',i0'): ',f0.2,' m')") &
                  i,par%xpointsw(i),par%ypointsw(i),xpoints(i),ypoints(i),mindistr(i)
                  call writelog('ls','(a)',txt)
               endif
            endif
         end do
      end if
      if (present(dmin)) then
         dmin = mindistr
      endif
      deallocate(mindistr)
   end subroutine snappointstogrid


   subroutine gridrotate_r0(t, x)
      use params
      use spaceparams
      use mnemmodule

      type(arraytype), intent(in)  :: t
      real*8                       :: x
      x = t%r0
   end subroutine gridrotate_r0

   subroutine gridrotate_r1(par, s, t, x)
      use params
      use spaceparams
      use mnemmodule

      type(parameters), intent(in) :: par
      type(spacepars), intent(in)  :: s
      type(arraytype), intent(in)  :: t
      real*8, dimension(:)         :: x
      real*8, parameter            :: pi = 4*atan(1.0d0)

      if (par%rotate .eq. 1 .and. s%nx .ne. -1) then
         select case(t%name)
          case(mnem_theta)
            x=270-(t%r1*(180/pi))
          case(mnem_theta0)
            x=270-(t%r1*(180/pi))
          case default
            x=t%r1
         end select
      else
         x=t%r1
      endif
   end subroutine gridrotate_r1


   subroutine gridrotate_r2(par, s, t, x &
#ifdef USEMPI
   , sl &
#endif
   )
      !
      ! sl: if present then this subroutine needs to be called by everyone
      !  in xmpi_ocomm, the extra needed s% variables needed will be collected
      !  x will not assigned to
      !
      use params
      use spaceparams
      use mnemmodule

      type(parameters), intent(in)        :: par
      type(spacepars), intent(inout)      :: s
      type(arraytype), intent(in)         :: t
      real*8, dimension(:,:), intent(out) :: x
#ifdef USEMPI
      type(spacepars), intent(in),optional :: sl
#else
      type(spacepars)                      :: sl
#endif

      logical                              :: getonly
      real*8, parameter                    :: pi = 4*atan(1.0d0)
      real*8, dimension(size(s%alfaz,1), size(s%alfaz,2)) :: Sutot, Svtot

#ifdef USEMPI
      getonly = present(sl)
#else
      getonly = .false.
#endif
      if(.not. getonly) then
         Sutot = 0;
         Svtot = 0;
      endif
      if (par%rotate .eq. 1) then
         select case(t%name)
          case(mnem_thetamean)
            x=270-((t%r2)*(180/pi))
          case(mnem_Fx)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_Fy)
            else
               x=t%r2*cos(s%alfaz)-s%Fy*sin(s%alfaz)
            endif
          case(mnem_Fy)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_Fx)
            else
               x=s%Fx*sin(s%alfaz)+t%r2*cos(s%alfaz)
            endif
          case(mnem_u)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_v)
            else
               x=t%r2*cos(s%alfaz)-s%v*sin(s%alfaz)
            endif
          case(mnem_gwu)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_gwv)
            else
               x=t%r2*cos(s%alfaz)-s%gwv*sin(s%alfaz)
            endif
          case(mnem_v)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_u)
            else
               x=s%u*sin(s%alfaz)+t%r2*cos(s%alfaz)
            endif
          case(mnem_gwv)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_gwu)
            else
               x=s%gwu*sin(s%alfaz)+t%r2*cos(s%alfaz)
            endif
          case(mnem_ue)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_ve)
            else
               x=t%r2*cos(s%alfaz)-s%ve*sin(s%alfaz)
            endif
          case(mnem_ve)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_ue)
            else
               x=s%ue*sin(s%alfaz)+t%r2*cos(s%alfaz)
            endif
          case(mnem_ui)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_vi)
            else
               x=t%r2*cos(s%alfaz)-s%vi*sin(s%alfaz)
            endif
          case(mnem_vi)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_ui)
            else
               x=s%ui*sin(s%alfaz)+t%r2*cos(s%alfaz)
            endif
          case(mnem_umean)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_vmean)
            else
               x=t%r2*cos(s%alfaz)-s%vmean*sin(s%alfaz)
            endif
          case(mnem_vmean)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_umean)
            else
               x=s%umean*sin(s%alfaz)+t%r2*cos(s%alfaz)
            endif
          case(mnem_uwf)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_vwf)
            else
               x=t%r2*cos(s%alfaz)-s%vwf*sin(s%alfaz)
            endif
          case(mnem_vwf)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_uwf)
            else
               x=s%uwf*sin(s%alfaz)+t%r2*cos(s%alfaz)
            endif
          case(mnem_Sutot)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_Subg)
               call space_collect_mnem(s,sl,par,mnem_Svbg)
               call space_collect_mnem(s,sl,par,mnem_Susg)
               call space_collect_mnem(s,sl,par,mnem_Svsg)
            else
               ! Jaap interpolate transports to water level points before rotating to real world
               Sutot = sum(s%Subg,dim=3)+sum(s%Susg,dim=3)
               Svtot = sum(s%Svbg,dim=3)+sum(s%Svsg,dim=3)
               Sutot(2:s%nx,:)=0.5d0*(Sutot(1:s%nx-1,:)+Sutot(2:s%nx,:))
               Svtot(:,2:s%ny)=0.5d0*(Svtot(:,1:s%ny-1)+Svtot(:,2:s%ny))
               x=Sutot*cos(s%alfaz) - Svtot*sin(s%alfaz)
            endif
          case(mnem_Svtot)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_Subg)
               call space_collect_mnem(s,sl,par,mnem_Svbg)
               call space_collect_mnem(s,sl,par,mnem_Susg)
               call space_collect_mnem(s,sl,par,mnem_Svsg)
            else
               Sutot = sum(s%Subg,dim=3)+sum(s%Susg,dim=3)
               Svtot = sum(s%Svbg,dim=3)+sum(s%Svsg,dim=3)
               Sutot(2:s%nx,:)=0.5d0*(Sutot(1:s%nx-1,:)+Sutot(2:s%nx,:))
               Svtot(:,2:s%ny)=0.5d0*(Svtot(:,1:s%ny-1)+Svtot(:,2:s%ny))
               x=Sutot*sin(s%alfaz) + Svtot*cos(s%alfaz)
            endif
          case(mnem_cctot)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_ccg)
            else
               x=sum(s%ccg,dim=3)
            endif
          case default
            if(getonly) then
               continue
            else
               x=t%r2
            endif
         end select
      else
         if(getonly) then
            continue
         else
            x=t%r2
         endif
      endif

   end subroutine gridrotate_r2
#ifndef USEMPI
   ! dummy subroutine for serial version to avoit many ifdefs in above and below
   ! subroutines
   subroutine space_collect_mnem(sg,sl,par,mnem)
      use params
      use spaceparams
      type(spacepars)   sg,sl
      type(parameters)  par
      character(*)      mnem
      print *,'This should not happen'
      print *,sg%nx,sl%nx,par%swave,mnem
   end subroutine space_collect_mnem
#endif

   subroutine gridrotate_r3(par, s, t, x &
#ifdef USEMPI
   , sl &
#endif
   )
      ! sl: see gridrotate_r2
      use params
      use spaceparams
      use mnemmodule
      use logging_module

      implicit none

      type(parameters), intent(in)  :: par
      type(spacepars), intent(in)   :: s
      type(arraytype), intent(in)   :: t
      real*8, dimension(:,:,:)      :: x
#ifdef USEMPI
      type(spacepars), intent(in),optional :: sl
#else
      type(spacepars)                      :: sl
#endif

      real*8, dimension(size(s%alfaz,1), size(s%alfaz,2), size(t%r3,3)) :: alfazr3,Susg,Svsg,Subg,Svbg
      real*8, parameter             :: pi = 4*atan(1.0d0)
      integer                       :: i
      logical                       :: getonly

#ifdef USEMPI
      getonly = present(sl)
#else
      getonly = .false.
#endif

      ! Jaap: initialize local transport variables (used for interpolation to water level points)
      Susg = 0;
      Svsg = 0;
      Subg = 0;
      Svbg = 0;

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

      if (par%rotate .eq. 1) then
         select case(t%name)
          case(mnem_cgx)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_cgy)
            else
               x=t%r3*cos(alfazr3)-s%cgy*sin(alfazr3)
            endif
          case(mnem_cgy)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_cgx)
            else
               x=s%cgx*sin(alfazr3)+t%r3*cos(alfazr3)
            endif
          case(mnem_cx)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_cy)
            else
               x=t%r3*cos(alfazr3)-s%cy*sin(alfazr3)
            endif
          case(mnem_cy)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_cx)
            else
               x=s%cx*sin(alfazr3)+t%r3*cos(alfazr3)
            endif
          case(mnem_thet)
            x=270-((s%thet+alfazr3)*(180/pi))
          case(mnem_Susg)
            ! Jaap: interpolate transports to water level points before rotating to world coordinates
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_Svsg)
            else
               Susg = t%r3
               Svsg = s%Svsg
               Susg(2:s%nx,:,:)=0.5d0*(t%r3(1:s%nx-1,:,:)+t%r3(2:s%nx,:,:))
               Svsg(:,2:s%ny,:)=0.5d0*(s%Svsg(:,1:s%ny-1,:)+s%Svsg(:,2:s%ny,:))
               x=Susg*cos(alfazr3)-Svsg*sin(alfazr3)
            endif
          case(mnem_Svsg)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_Susg)
            else
               Susg = s%Susg
               Svsg = t%r3
               Susg(2:s%nx,:,:)=0.5d0*(s%Susg(1:s%nx-1,:,:)+s%Susg(2:s%nx,:,:))
               Svsg(:,2:s%ny,:)=0.5d0*(t%r3(:,1:s%ny-1,:)+t%r3(:,2:s%ny,:))
               x=s%Susg*sin(alfazr3)+t%r3*cos(alfazr3)
            endif
          case(mnem_Subg)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_Svbg)
            else
               Subg = t%r3
               Svbg = s%Svbg
               Subg(2:s%nx,:,:)=0.5d0*(t%r3(1:s%nx-1,:,:)+t%r3(2:s%nx,:,:))
               Svbg(:,2:s%ny,:)=0.5d0*(s%Svbg(:,1:s%ny-1,:)+s%Svbg(:,2:s%ny,:))
               x=Subg*cos(alfazr3)-Svbg*sin(alfazr3)
            endif
          case(mnem_Svbg)
            if(getonly) then
               call space_collect_mnem(s,sl,par,mnem_Subg)
            else
               Subg = s%Subg
               Svbg = t%r3
               Subg(2:s%nx,:,:)=0.5d0*(s%Subg(1:s%nx-1,:,:)+s%Subg(2:s%nx,:,:))
               Svbg(:,2:s%ny,:)=0.5d0*(t%r3(:,1:s%ny-1,:)+t%r3(:,2:s%ny,:))
               x=Subg*sin(alfazr3)+Svbg*cos(alfazr3)
            endif
          case default
            if(getonly) then
               continue
            else
               x=t%r3
            endif
         end select
      else

         if(getonly) then
            continue
         else
            x = t%r3
         endif
      endif
   end subroutine gridrotate_r3

   subroutine gridrotate_r4(t, x)
      use params
      use spaceparams
      use mnemmodule

      type(arraytype), intent(in)   :: t
      real*8, dimension(:,:,:,:)    :: x
      x = t%r4
   end subroutine gridrotate_r4


   subroutine gridrotate_i0(t, x)
      use params
      use spaceparams
      use mnemmodule

      type(arraytype), intent(in)   :: t
      integer                       :: x
      x = t%i0
   end subroutine gridrotate_i0

   subroutine gridrotate_i1(t, x)
      use params
      use spaceparams
      use mnemmodule

      type(arraytype), intent(in)   :: t
      integer, dimension(:),intent(out) :: x
      x = t%i1
   end subroutine gridrotate_i1

   subroutine gridrotate_i2(t, x)
      use params
      use spaceparams
      use mnemmodule

      type(arraytype), intent(in)   :: t
      integer, dimension(:,:),intent(out) :: x
      x = t%i2
   end subroutine gridrotate_i2

   subroutine gridrotate_i3(t, x)
      use params
      use spaceparams
      use mnemmodule

      type(arraytype), intent(in)   :: t
      integer, dimension(:,:,:),intent(out) :: x
      x = t%i3
   end subroutine gridrotate_i3

   subroutine gridrotate_i4(t, x)
      use params
      use spaceparams
      use mnemmodule

      type(arraytype), intent(in)   :: t
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
