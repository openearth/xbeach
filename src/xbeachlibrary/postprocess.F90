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
                  call writelog('ls','(a)',txt,xomaster)
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


   subroutine gridrotate_r2(par, s, t, x , row, col )
      !
      ! sl: if present then this subroutine needs to be called by everyone
      !  in xmpi_ocomm, the extra needed s% variables needed will be collected
      !  x will not assigned to
      !
      ! if row is present, and row>0, the rotation will be done only for
      ! s%some_arry(row,col) and assigned to x(1,1)

      use params
      use spaceparams
      use mnemmodule

      type(parameters), intent(in)        :: par
      type(spacepars), intent(inout)      :: s
      type(arraytype), intent(in)         :: t
      real*8, dimension(:,:), intent(out) :: x
      integer, intent(in), optional       :: row, col

      real*8, parameter                    :: pi = 4*atan(1.0d0)
      real*8, dimension(size(s%alfaz,1), size(s%alfaz,2)) :: Sutot, Svtot

      integer is,  ie,  js,  je
      integer isx, iex, jsx, jex

         Sutot = 0;
         Svtot = 0;

      is  = 1
      ie  = size(x,1)
      js  = 1
      je  = size(x,2)
      isx = is
      iex = ie
      jsx = js
      jex = je

      if (present(row)) then
         if (row .gt. 0 ) then
            is  = row
            ie  = row
            js  = col
            je  = col
            isx = 1
            iex = 1
            jsx = 1
            jex = 1
      endif
      endif

      if (par%rotate .eq. 1) then
         select case(t%name)
          case(mnem_thetamean)
            x(isx:iex,jsx:jex) = 270-((t%r2(is:ie,js:je))*(180/pi))
          case(mnem_Fx)
            x(isx:iex,jsx:jex) = t%r2(is:ie,js:je)*cos(s%alfaz(is:ie,js:je))- &
            &                    s%Fy(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))
          case(mnem_Fy)
            x(isx:iex,jsx:jex) = s%Fx(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))+ &
            &                    t%r2(is:ie,js:je)*cos(s%alfaz(is:ie,js:je))
          case(mnem_u)
            x(isx:iex,jsx:jex) = t%r2(is:ie,js:je)*cos(s%alfaz(is:ie,js:je))- &
            &                    s%v (is:ie,js:je)*sin(s%alfaz(is:ie,js:je))
          case(mnem_gwu)
            x(isx:iex,jsx:jex) = t%r2 (is:ie,js:je)*cos(s%alfaz(is:ie,js:je))- &
            &                    s%gwv(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))
          case(mnem_v)
            x(isx:iex,jsx:jex) = s%u (is:ie,js:je)*sin(s%alfaz(is:ie,js:je))+ &
            &                    t%r2(is:ie,js:je)*cos(s%alfaz(is:ie,js:je))
          case(mnem_gwv)
            x(isx:iex,jsx:jex) = s%gwu(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))+ &
            &                    t%r2 (is:ie,js:je)*cos(s%alfaz(is:ie,js:je))
          case(mnem_ue)
            x(isx:iex,jsx:jex) = t%r2(is:ie,js:je)*cos(s%alfaz(is:ie,js:je))- &
            &                    s%ve(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))
          case(mnem_ve)
            x(isx:iex,jsx:jex) = s%ue(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))+ &
            &                    t%r2(is:ie,js:je)*cos(s%alfaz(is:ie,js:je))
          case(mnem_ui)
            x(isx:iex,jsx:jex) = t%r2(is:ie,js:je)*cos(s%alfaz(is:ie,js:je))- &
            &                    s%vi(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))
          case(mnem_vi)
            x(isx:iex,jsx:jex) = s%ui(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))+ &
            &                    t%r2(is:ie,js:je)*cos(s%alfaz(is:ie,js:je))
          case(mnem_umean)
            x(isx:iex,jsx:jex) = t%r2   (is:ie,js:je)*cos(s%alfaz(is:ie,js:je))- &
            &                    s%vmean(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))
          case(mnem_vmean)
            x(isx:iex,jsx:jex) = s%umean(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))+ &
            &                    t%r2   (is:ie,js:je)*cos(s%alfaz(is:ie,js:je))
          case(mnem_uwf)
            x(isx:iex,jsx:jex) = t%r2 (is:ie,js:je)*cos(s%alfaz(is:ie,js:je))- &
            &                    s%vwf(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))
          case(mnem_vwf)
            x(isx:iex,jsx:jex) = s%uwf(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))+ &
            &                    t%r2 (is:ie,js:je)*cos(s%alfaz(is:ie,js:je))
          case(mnem_Sutot)
               ! Jaap interpolate transports to water level points before rotating to real world
               Sutot = sum(s%Subg,dim=3)+sum(s%Susg,dim=3)
               Svtot = sum(s%Svbg,dim=3)+sum(s%Svsg,dim=3)
               Sutot(2:s%nx,:)=0.5d0*(Sutot(1:s%nx-1,:)+Sutot(2:s%nx,:))
               Svtot(:,2:s%ny)=0.5d0*(Svtot(:,1:s%ny-1)+Svtot(:,2:s%ny))
            x(isx:iex,jsx:jex) = Sutot(is:ie,js:je)*cos(s%alfaz(is:ie,js:je)) - &
            &                    Svtot(is:ie,js:je)*sin(s%alfaz(is:ie,js:je))
          case(mnem_Svtot)
               Sutot = sum(s%Subg,dim=3)+sum(s%Susg,dim=3)
               Svtot = sum(s%Svbg,dim=3)+sum(s%Svsg,dim=3)
               Sutot(2:s%nx,:)=0.5d0*(Sutot(1:s%nx-1,:)+Sutot(2:s%nx,:))
               Svtot(:,2:s%ny)=0.5d0*(Svtot(:,1:s%ny-1)+Svtot(:,2:s%ny))
            x(isx:iex,jsx:jex) = Sutot(is:ie,js:je)*sin(s%alfaz(is:ie,js:je)) + &
            &                    Svtot(is:ie,js:je)*cos(s%alfaz(is:ie,js:je))
          case(mnem_cctot)
            x(isx:iex,jsx:jex) = sum(s%ccg(is:ie,js:je,:),dim=3)
          case default
            x(isx:iex,jsx:jex) = t%r2(is:ie,js:je)
         end select
      else
         x(isx:iex,jsx:jex) = t%r2(is:ie,js:je)
      endif

   end subroutine gridrotate_r2

   subroutine gridrotate_r3(par, s, t, x, row, col )
      ! if row is present, and row>0, the rotation will be done only for
      ! s%some_arry(row,col,:) and assigned to x(1,1,:)

      use params
      use spaceparams
      use mnemmodule
      use logging_module

      implicit none

      type(parameters), intent(in)  :: par
      type(spacepars), intent(in)   :: s
      type(arraytype), intent(in)   :: t
      real*8, dimension(:,:,:)      :: x
      integer, intent(in), optional :: row, col

      real*8, dimension(size(s%alfaz,1), size(s%alfaz,2), size(t%r3,3)) :: alfazr3,Susg,Svsg,Subg,Svbg
      real*8, parameter             :: pi = 4*atan(1.0d0)
      integer                       :: i

      integer is,  ie,  js,  je
      integer isx, iex, jsx, jex

      if (present(row)) then
         if (row .gt. 0 ) then
            is  = row
            ie  = row
            js  = col
            je  = col
            isx = 1
            iex = 1
            jsx = 1
            jex = 1
         endif
      else
         is  = 1
         ie  = size(x,1)
         js  = 1
         je  = size(x,2)
         isx = is
         iex = ie
         jsx = js
         jex = je
      endif

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
            x(isx:iex,jsx:jex,:) = t%r3 (is:ie,js:je,:)*cos(alfazr3(is:ie,js:je,:))- &
            &                      s%cgy(is:ie,js:je,:)*sin(alfazr3(is:ie,js:je,:))
          case(mnem_cgy)
            x(isx:iex,jsx:jex,:) = s%cgx(is:ie,js:je,:)*sin(alfazr3(is:ie,js:je,:))+ &
            &                      t%r3 (is:ie,js:je,:)*cos(alfazr3(is:ie,js:je,:))
          case(mnem_cx)
            x(isx:iex,jsx:jex,:) = t%r3(is:ie,js:je,:)*cos(alfazr3(is:ie,js:je,:))- &
            &                      s%cy(is:ie,js:je,:)*sin(alfazr3(is:ie,js:je,:))
          case(mnem_cy)
            x(isx:iex,jsx:jex,:) = s%cx(is:ie,js:je,:)*sin(alfazr3(is:ie,js:je,:))+ &
            &                      t%r3(is:ie,js:je,:)*cos(alfazr3(is:ie,js:je,:))
          case(mnem_thet)
            x(isx:iex,jsx:jex,:) = 270-((s%thet(is:ie,js:je,:)+alfazr3(is:ie,js:je,:))*(180/pi))
          case(mnem_Susg)
            ! Jaap: interpolate transports to water level points before rotating to world coordinates
               Susg = t%r3
               Svsg = s%Svsg
               Susg(2:s%nx,:,:)=0.5d0*(t%r3(1:s%nx-1,:,:)+t%r3(2:s%nx,:,:))
               Svsg(:,2:s%ny,:)=0.5d0*(s%Svsg(:,1:s%ny-1,:)+s%Svsg(:,2:s%ny,:))
            x(isx:iex,jsx:jex,:) = Susg(is:ie,js:je,:)*cos(alfazr3(is:ie,js:je,:))- &
            &                      Svsg(is:ie,js:je,:)*sin(alfazr3(is:ie,js:je,:))
          case(mnem_Svsg)
               Susg = s%Susg
               Svsg = t%r3
               Susg(2:s%nx,:,:)=0.5d0*(s%Susg(1:s%nx-1,:,:)+s%Susg(2:s%nx,:,:))
               Svsg(:,2:s%ny,:)=0.5d0*(t%r3(:,1:s%ny-1,:)+t%r3(:,2:s%ny,:))
            x(isx:iex,jsx:jex,:) = Susg(is:ie,js:je,:)*sin(alfazr3(is:ie,js:je,:))+ &
            &                      Svsg(is:ie,js:je,:)*cos(alfazr3(is:ie,js:je,:))
          case(mnem_Subg)
               Subg = t%r3
               Svbg = s%Svbg
               Subg(2:s%nx,:,:)=0.5d0*(t%r3(1:s%nx-1,:,:)+t%r3(2:s%nx,:,:))
               Svbg(:,2:s%ny,:)=0.5d0*(s%Svbg(:,1:s%ny-1,:)+s%Svbg(:,2:s%ny,:))
            x(isx:iex,jsx:jex,:) = Subg(is:ie,js:je,:)*cos(alfazr3(is:ie,js:je,:))- &
            &                      Svbg(is:ie,js:je,:)*sin(alfazr3(is:ie,js:je,:))
          case(mnem_Svbg)
               Subg = s%Subg
               Svbg = t%r3
               Subg(2:s%nx,:,:)=0.5d0*(s%Subg(1:s%nx-1,:,:)+s%Subg(2:s%nx,:,:))
               Svbg(:,2:s%ny,:)=0.5d0*(t%r3(:,1:s%ny-1,:)+t%r3(:,2:s%ny,:))
            x(isx:iex,jsx:jex,:) = Subg(is:ie,js:je,:)*sin(alfazr3(is:ie,js:je,:))+ &
            &                      Svbg(is:ie,js:je,:)*cos(alfazr3(is:ie,js:je,:))
          case default
            x(isx:iex,jsx:jex,:) = t%r3(is:ie,js:je,:)
         end select
      else

         x(isx:iex,jsx:jex,:) = t%r3(is:ie,js:je,:)
      endif
   end subroutine gridrotate_r3

   subroutine gridrotate_r4(t, x, row, col)
      use params
      use spaceparams
      use mnemmodule

      type(arraytype), intent(in)   :: t
      real*8, dimension(:,:,:,:)    :: x
      integer, optional, intent(in) :: row,col
      integer is,  ie,  js,  je
      integer isx, iex, jsx, jex

      if (present(row)) then
         if (row .gt. 0 ) then
            is  = row
            ie  = row
            js  = col
            je  = col
            isx = 1
            iex = 1
            jsx = 1
            jex = 1
         endif
      else
         is  = 1
         ie  = size(x,1)
         js  = 1
         je  = size(x,2)
         isx = is
         iex = ie
         jsx = js
         jex = je
      endif

      x(isx:iex,jsx:jex,:,:) = t%r4(is:ie,js:je,:,:)

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
   
   function get_sister_mnem(mnem) result(sistermnem)
      use mnemmodule
      character(maxnamelen)     :: sistermnem
      character(maxnamelen)     :: mnem
      
      select case (mnem)
          case(mnem_Fx)
             sistermnem = mnem_Fy
          case(mnem_Fy)
             sistermnem = mnem_Fx
          case(mnem_u)
            sistermnem = mnem_v
          case(mnem_gwu)
            sistermnem = mnem_gwv
          case(mnem_v)
            sistermnem = mnem_u
          case(mnem_gwv)
            sistermnem = mnem_gwu
          case(mnem_ue)
           sistermnem = mnem_ve
          case(mnem_ve)
            sistermnem = mnem_ue
          case(mnem_ui)
            sistermnem = mnem_vi
          case(mnem_vi)
            sistermnem = mnem_ui
          case(mnem_umean)
            sistermnem = mnem_vmean
          case(mnem_vmean)
            sistermnem = mnem_umean
          case(mnem_uwf)
            sistermnem = mnem_vwf
          case(mnem_vwf)
            sistermnem = mnem_uwf
          case(mnem_Sutot)
              sistermnem = mnem_Svtot
          case(mnem_Svtot)
              sistermnem = mnem_Sutot
          case(mnem_cctot)
              sistermnem = mnem_ccg
          case(mnem_cgx)
            sistermnem = mnem_cgy
          case(mnem_cgy)
            sistermnem = mnem_cgx
          case(mnem_cx)
            sistermnem = mnem_cy
          case(mnem_cy)
            sistermnem = mnem_cx
          case(mnem_Susg)
            sistermnem = mnem_Svsg
          case(mnem_Svsg)
             sistermnem = mnem_Susg
          case(mnem_Subg)
               sistermnem = mnem_Svbg
          case(mnem_Svbg)
              sistermnem = mnem_Subg
      end select
   end function get_sister_mnem
end module postprocessmod
