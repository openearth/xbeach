!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   MODULE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Note: all code here has been made dead by the following #if 0
! the file is preserved for reference
#if 0
!
! if defined OUTSINGLE, varoutput will output all reals in single precision
!  otherwise: double precision
!
#define OUTSINGLE
#ifdef OUTSINGLE
#define CONVREAL sngl
#else
#define CONVREAL
#endif

module fortoutput_module
   use xmpi_module
   use mnemmodule
   use means_module

   implicit none
   private
   public var_output_init, var_output

   ! Robert: Add choice of output variables
   !logical,dimension(999)              :: outputindex ! [-]      tracks which  global variables are to be outputted.

   integer*4,dimension(:),allocatable,save  :: crosstype   ! 0 = cross shore (x), 1 = longshore (y)
   integer*4,dimension(:),allocatable,save  :: xpoints     ! model x-coordinate of output points
   integer*4,dimension(:),allocatable,save  :: ypoints     ! model y-coordinate of output points
   integer*4,dimension(:),allocatable,save  :: xcross      ! model x-coordinate of output cross sections
   integer*4,dimension(:),allocatable,save  :: ycross      ! model y-coordinate of output cross sections
   integer*4,dimension(:),allocatable,save  :: nvarcross   ! vector with number of output variable per output cross section
   integer*4,dimension(:,:),allocatable,save:: Avarpoint   ! Array with associated index of output variables per point
   integer*4,dimension(:,:),allocatable,save:: Avarcross   ! Array with associated index of output variables per cross section
   integer*4,dimension(:,:),allocatable,save:: rugmaskg,rugmaskl   ! Mask with 1 for row with runup gauge, 0 without. Dimensions: nx+1,ny+1
   ! One for global field, one for mpi subdomain if using MPI
   integer*4,dimension(:),allocatable,save  :: rugrowindex ! Array with row index where runup gauge can be found
   ! Only alive at xmaster
   integer*4,save                           :: stpm        ! size of tpm

   ! Store the global variables in numbers....
   integer,save                             :: noutnumbers = 0  ! the number of outnumbers
   integer, dimension(numvars),save         :: outnumbers  ! numbers, corrsponding to mnemonics, which are to be output


   integer                             :: itg,itp,itc,itm,itd,day,ot
   type(arraytype)                     :: At

   interface outarray
      module procedure outarray_r0
      module procedure outarray_r1
      module procedure outarray_r2
      module procedure outarray_r3
      module procedure outarray_r4
      module procedure outarray_i0
      module procedure outarray_i1
      module procedure outarray_i2
      module procedure outarray_i3
      module procedure outarray_i4
   end interface outarray

contains



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!   INITIALISE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   subroutine var_output_init(s,sl,par,tpar)
      use params
      use spaceparams
      use readkey_module
      use indextos_module
      use timestep_module
      use logging_module
      use postprocessmod
      use filefunctions
#ifdef USEMPI
      use general_mpi_module
#endif

      IMPLICIT NONE

      type(spacepars),intent(in)          :: s,sl
      type(parameters),intent(in)         :: par
      type(timepars),intent(in)           :: tpar

      integer                             :: i,j
      integer                             :: i1,i2,i3
      integer                             :: reclen,reclenp,wordsize,reclenm
      integer                             :: fid
      character(99)                       :: fname,fnamemean,fnamevar,fnamemin,fnamemax
      type(arraytype)                     :: t

#ifdef USEMPI
      logical                             :: toall = .true.
#endif

      reclenm = -123

      ! Initialize places in output files
      itg = 0
      itm = 0
      itp = 0
      itc = 0
      itd = 0
      stpm = size(tpar%tpm)

      ! Record size for global and mean output
      inquire(iolength=wordsize) CONVREAL(1.d0)
      reclen=wordsize*(s%nx+1)*(s%ny+1)

      !!!!! XY.DAT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      if (xomaster) then
         open(100,file='xy.dat',form='unformatted',access='direct',recl=reclen,status='REPLACE')
         write(100,rec=1)CONVREAL(s%xz)
         write(100,rec=2)CONVREAL(s%yz)
         write(100,rec=3)CONVREAL(s%x)
         write(100,rec=4)CONVREAL(s%y)
         close(100)
      endif

      !!!!! GLOBAL VARS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! This module identifies output variables by their index
      ! store this at the module level so we don't have to pass par around
      ! do i=1,size(par%globalvars)
      !   if (trim(par%globalvars(i))=='abc') then
      !      exit
      !      endif
      ! enddo
      noutnumbers = par%nglobalvar
      ! store all indices for the global variables
      do i= 1,noutnumbers
         outnumbers(i) = chartoindex(par%globalvars(i))
      enddo

#ifdef USEMPI
      call xmpi_bcast(noutnumbers,toall)
      call xmpi_bcast(outnumbers,toall)
#endif

      !!!!! OUTPUT POINTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! read from par to local
      if (par%npoints + par%nrugauge > 0) then
         allocate(xpoints(par%npoints+par%nrugauge))
         allocate(ypoints(par%npoints+par%nrugauge))
         ! 3 for rugauge, npointvar for points
         allocate(Avarpoint(par%npoints+par%nrugauge,max(par%npointvar,par%nrugdepth*3)))
         xpoints=0
         ypoints=0

         ! Convert world coordinates of points to nearest (lsm) grid point
         if (xmaster) then
            call snappointstogrid(par, s, xpoints, ypoints)
            do i=1,par%npoints
               do j=1,par%npointvar
                  Avarpoint(i,j) = chartoindex(trim(par%pointvars(j)))
               end do
            end do
            do i=par%npoints+1,par%npoints+par%nrugauge
               Avarpoint(i,1) = chartoindex('xz')
               Avarpoint(i,2) = chartoindex('yz')
               Avarpoint(i,3) = chartoindex('zs')
            enddo
         endif
#ifdef USEMPI
         call xmpi_bcast(xpoints,toall)
         call xmpi_bcast(ypoints,toall)
         call xmpi_bcast(Avarpoint,toall)
#endif
         ! make mask for grid rows which include
         if (par%nrugauge>0) then
            allocate(rugrowindex(par%nrugauge))
#ifdef USEMPI
            allocate(rugmaskg(s%nx+1,s%ny+1))
            allocate(rugmaskl(sl%nx+1,sl%ny+1))
#endif
            ! Make rugrowindex with row number (per subprocess) with runup gauge
            do i=1,par%nrugauge
#ifndef USEMPI
               ! very easy
               rugrowindex(i)=ypoints(par%npoints+i)
#else
               ! very complicated
               if (xmaster .or. xomaster) then
                  ! generate rugmask on global grid
                  do j=1,s%ny+1
                     if (j==ypoints(par%npoints+i)) then
                        rugmaskg(:,j)=1
                     else
                        rugmaskg(:,j)=0
                     endif
                  enddo
               endif
               ! now distribute rugmask global to rugmask local on all subgrids
               if(xcompute) then
                  call matrix_distr(rugmaskg,rugmaskl,sl%is,sl%lm,sl%js,sl%ln,xmpi_master,xmpi_comm)
               endif
               ! everybody has their own rugmaskl. Look to see if your domain has "1" in rugmaskl
               ! first assume that no runup gauge exists in this domain, so we set rowindex to zero
               rugrowindex(i)=0    ! rugrowindex is local on all subgrid
               if (xcompute) then  ! rugmaskl has no meaning on xomaster
                  do j=1,sl%ny+1
                     if (rugmaskl(1,j)==1) rugrowindex(i)=j   ! okay, there is a runup gauge this row
                  enddo
               endif
#endif
            enddo  ! i=1,par%nrugauge
#ifdef USEMPI
            deallocate(rugmaskg)   ! not needed anymore
            deallocate(rugmaskl)   ! not needed anymore
#endif
         endif  ! runup gauge > 0
         !
         ! First time file opening for point output
         !
         if (xomaster) then
            do i=1,par%npoints+par%nrugauge
               fname = ''
               if (par%pointtypes(i)==0) then
                  fname(1:5)='point'
                  i1=floor(real(i)/100.d0)
                  i2=floor(real(i-i1*100)/10.d0)
                  i3=i-i1*100-i2*10
               else
                  fname(1:5)='rugau'
                  i1=floor(real(i-par%npoints)/100.d0)
                  i2=floor(real((i-par%npoints)-i1*100)/10.d0)
                  i3=(i-par%npoints)-i1*100-i2*10
               endif
               fname(6:6)=char(48+i1)
               fname(7:7)=char(48+i2)
               fname(8:8)=char(48+i3)
               fname(9:12)='.dat'
               if (par%pointtypes(i)==0) then
                  reclenp=wordsize*(par%npointvar+1)*1
               else
                  reclenp=wordsize*(1+par%nrugdepth*3)*1
               endif
               open(indextopointsunit(i),file=fname,&
               form='unformatted',access='direct',recl=reclenp,status='REPLACE')
            enddo
            if (par%npoints>0) then
               ! write index file of point output variables
               fid=create_new_fid()
               open(fid,file='pointvars.idx',status='replace',action='write')
               do i=1,par%npointvar
                  write(fid,*)trim(par%pointvars(i))
               enddo
               close(fid)
            endif
         endif ! xomaster
      end if ! npoints+nrugauge>0


      !!!!! TIME-AVEARGE, VARIANCE and MIN-MAX ARRAYS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      if (par%nmeanvar>0) then
         !! First time file opening for time-average output
         if(xomaster) then
            do i=1,par%nmeanvar
               call makeaveragenames(chartoindex(trim(par%meanvars(i))),fnamemean,fnamevar,fnamemin,fnamemax)
               fnamemean=trim(fnamemean)
               fnamevar =trim(fnamevar)
               fnamemin =trim(fnamemin)
               fnamemax =trim(fnamemax)
               call indextos(s,chartoindex(trim(par%meanvars(i))),t)
               reclenm = wordsize
               select case(t%rank)
                case (2)
                  reclenm = wordsize*size(t%r2)
                case (3)
                  reclenm=wordsize*size(t%r3)
                case (4)
                  reclenm=wordsize*size(t%r4)
               end select
               open(indextomeanunit(i),file=fnamemean,form='unformatted',access='direct',recl=reclenm,status='REPLACE')
               open(indextovarunit(i) ,file=fnamevar ,form='unformatted',access='direct',recl=reclenm,status='REPLACE')
               open(indextominunit(i) ,file=fnamemin ,form='unformatted',access='direct',recl=reclenm,status='REPLACE')
               open(indextomaxunit(i),file=fnamemax,  form='unformatted',access='direct',recl=reclenm,status='REPLACE')
            enddo
         endif
      endif ! par%nmeanvar > 0


      !
      ! drifter output files
      !
      if (par%ndrifter>0) then
         if (xomaster) then
            inquire(iolength=wordsize) CONVREAL(1.d0)

            reclen=wordsize*3
            do i=1,par%ndrifter
               write(fname(7:10),'(i4)')i+1000
               fname(1:7)='drifter'
               fname(11:14)='.dat'
               open(indextodrifterunit(i),file=fname,form='unformatted',access='direct',recl=reclen,status='REPLACE')
            enddo
         endif
      endif ! par%ndrifter >0
      ! wwvv to avoid warning about unused sl:
      if (sl%nx .eq. -1) return

   end subroutine var_output_init

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!   OUTPUT AT EVERY TIMESTEP    !!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine var_output(s,sl,par,tpar)
      use params
      use spaceparams
      use timestep_module
      use indextos_module
      use logging_module
#ifdef USEMPI
      use xmpi_module
#endif

      IMPLICIT NONE

      type(spacepars)                         :: s,sl
      type(parameters)                        :: par
      type(timepars),intent(in)               :: tpar
      integer                                 :: i,ii,ird
      !  integer                                 :: i1,i2,i3
      integer                                 :: wordsize,idumhl
      integer                                 :: xmax
#ifdef USEMPI
      integer                                 :: xrank
#endif
      real*8,dimension(numvars)               :: intpvector
      real*8,dimension(numvars,s%nx+1)        :: crossvararray0
      real*8,dimension(numvars,s%ny+1)        :: crossvararray1
      integer,dimension(:),allocatable        :: tempvectori
      real*8,dimension(:),allocatable         :: tempvectorr
      real*8,dimension(size(tpar%tpg)+size(tpar%tpp)+size(tpar%tpc)+size(tpar%tpm)) :: outputtimes
      type(arraytype)                         :: t

      integer                                 :: iz,jz
      real*8                                  :: di,dj,dx,dy

#ifdef USEMPI
      logical                                 :: toall = .true.
#endif

#ifdef USEMPI
      xrank = huge(xrank)
#endif
      inquire(iolength=wordsize) CONVREAL(1.d0)
      !  reclen=wordsize*(s%nx+1)*(s%ny+1)
      !  reclen2=wordsize*(s%nx+1)*(s%ny+1)*(par%ngd)*(par%nd)

      ! Determine if this is an output timestep
      if (tpar%output) then
         !!! Write at every output timestep

         !!! Write point variables
         !!! Only write if it is output time for points
         if (par%npoints+par%nrugauge>0) then

            !! Runup gauge depth

            if (tpar%outputp) then
               itp=itp+1
               !
               ! Set up runup gauge output vector
               allocate(tempvectorr(1+par%nrugdepth*3))
               tempvectorr=huge(0.d0)
               do i=1,par%nrugauge
                  do ird=1,par%nrugdepth
                     ! in MPI we only want to cycle through sub matrix, else through whole matrix
#ifdef USEMPI
                     xmax = sl%nx+1
                     xrank  = huge(xrank) ! Set default, so processes not involved in
                     ! runup gauge do not affect all_reduce statement
                     idumhl = xmax        ! Set default
                     if (rugrowindex(i)>0) then  ! this (sub) domain contains this runup gauge
                        ! local index of minimum location where hh<rugdepth
                        do ii=2,xmax
                           if ((sl%hh(ii,rugrowindex(i))<=par%rugdepth(ird)) .and. &
                           (sl%hh(ii-1,rugrowindex(i))>par%rugdepth(ird)) ) then
                              idumhl= ii-1
                              xrank = xmpi_rank  ! the row number of this process in the MPI grid of subdomains
                              exit
                           endif
                        enddo
                     endif
#else
                     xmax = s%nx+1
                     idumhl = xmax        ! Set default
                     if (rugrowindex(i)>0) then  ! master domain always contains this runup gauge
                        ! local index of minimum location where hh<rugdepth
                        do ii=2,xmax
                           if ((s%hh(ii,rugrowindex(i))<=par%rugdepth(ird)) .and. &
                           (s%hh(ii-1,rugrowindex(i))>par%rugdepth(ird)) ) then
                              idumhl=ii-1
                              exit
                           endif
                        enddo
                     endif
#endif

#ifdef USEMPI
                     ! In MPI multiple domains may have a non-zero value for idumhl, so we choose the one with the
                     ! lowest MPI rank (closest to xmpi_top, or the offshore boundary)
                     !call xmpi_allreduce(xrank,MPI_MIN) ! wwvv-todo
                     ! now only look at this process
                     !if (xmpi_rank==xrank) then
                     if(xmaster) then
                        if (par%morfacopt==1) then
                           tempvectorr(1)=par%t*max(par%morfac,1.d0)
                        else
                           tempvectorr(1)=par%t
                        endif
                        tempvectorr((ird-1)*3+2)=sl%x(idumhl,rugrowindex(i))
                        tempvectorr((ird-1)*3+3)=sl%y(idumhl,rugrowindex(i))
                        tempvectorr((ird-1)*3+4)=sl%zs(idumhl,rugrowindex(i))
                     endif
                     ! Reduce the whole set to only the real numbers in tempvectori in xmpi_rank
                     ! wwvv this only works, because the lowest rank belongs to xmaster
                     if(xcompute) call xmpi_allreduce(tempvectorr,MPI_MIN)
                     ! wwvv xomaster needs this tempvectorr:
                     call xmpi_bcast(tempvectorr,toall)  ! wwvv should make a allreduce toall
#else
                     if (par%morfacopt==1) then
                        tempvectorr(1)=par%t*max(par%morfac,1.d0)
                     else
                        tempvectorr(1)=par%t
                     endif
                     tempvectorr((ird-1)*3+2)=s%xz(idumhl,rugrowindex(i))
                     tempvectorr((ird-1)*3+3)=s%yz(idumhl,rugrowindex(i))
                     tempvectorr((ird-1)*3+4)=s%zs(idumhl,rugrowindex(i))
#endif
                  enddo
                  if (xomaster) write(indextopointsunit(i+par%npoints),rec=itp)CONVREAL(tempvectorr)
               enddo  ! i=1,par%nrugauge
               deallocate(tempvectorr)

               !! point output

               do i=1,par%npoints
                  !!! Make vector of all s% values at n,m grid coordinate
                  call makeintpvector(par,sl,intpvector,xpoints(i),ypoints(i))
                  if (xomaster) then
                     allocate(tempvectori(par%npointvar))
                     tempvectori=Avarpoint(i,1:par%npointvar)
                     allocate(tempvectorr(par%npointvar+1))
                     tempvectorr=0.d0
                     do ii=1,par%npointvar
                        tempvectorr(ii+1)=intpvector(tempvectori(ii))
                     enddo
                     if (par%morfacopt==1) then
                        tempvectorr(1)=par%t*max(par%morfac,1.d0)
                     else
                        tempvectorr(1)=par%t
                     endif
                     write(indextopointsunit(i),rec=itp)CONVREAL(tempvectorr)
                     deallocate(tempvectori)
                     deallocate(tempvectorr)
                  endif
               enddo
            endif
         endif

         !!! Write cross section variables
         !!! Only write if it is output time for cross sections
         if (par%ncross>0) then
            if (tpar%outputc) then
               itc=itc+1
               do i=1,par%ncross
                  !!! Make array of all s% values at n or m grid line
                  if (crosstype(i)==0) then
                     !    makecrossvector(par, s, local s, output array, no of variables, index of variables in output, m or n coordinate, cross section type)
                     call makecrossvector(s,sl,par,crossvararray0,nvarcross(i),Avarcross(i,1:nvarcross(i)),ycross(i),0)
                     if (xomaster) write(indextocrossunit(i),rec=itc)CONVREAL(crossvararray0(1:nvarcross(i),:))
                  else
                     call makecrossvector(s,sl,par,crossvararray1,nvarcross(i),Avarcross(i,1:nvarcross(i)),xcross(i),1)
                     if (xomaster) write(indextocrossunit(i),rec=itc)CONVREAL(crossvararray1(1:nvarcross(i),:))
                  endif
               enddo
            endif
         endif

         !!! Collect mean variable to global grid
         if (par%nmeanvar>0) then
            ! Not at the first in tpm as this is the start of averaging. Only output after second in tpm
            if (tpar%outputm .and. tpar%itm>1) then
               do i=1,par%nmeanvar
#ifdef USEMPI
                  call means_collect(sl,meansparsglobal(i),meansparslocal(i))
#else
                  meansparsglobal(i)=meansparslocal(i)
#endif
               enddo
            endif
         endif
         !!! Write average variables

         if (par%nmeanvar>0) then
            ! Not at the first in tpm as this is the start of averaging. Only output after second in tpm
            if (tpar%outputm .and. tpar%itm>1) then
               itm=itm+1  ! Note, this is a local counter, used to position in output file
               do i=1,par%nmeanvar
                  select case (meansparsglobal(i)%rank)
                   case (2)
                     if(xomaster) then
                        if (trim(par%meanvars(i))=='H') then                ! Hrms changed to H
                           write(indextomeanunit(i),rec=itm)CONVREAL(sqrt(meansparsglobal(i)%variancesquareterm2d))
                        elseif (trim(par%meanvars(i))=='urms') then       ! urms
                           write(indextomeanunit(i),rec=itm)CONVREAL(sqrt(meansparsglobal(i)%variancesquareterm2d))
                        elseif (trim(par%meanvars(i))=='thetamean') then       ! thetamean
                           write(indextomeanunit(i),rec=itm) &
                           CONVREAL( &
                           mod(2.d0*par%px + atan2(nint(meansparsglobal(i)%mean2d)/1d7, &
                           mod(meansparsglobal(i)%mean2d,1.d0)*1d1), 2.d0*par%px) / par%px * 180 &
                           )
                        else                                                    ! non-rms variables
                           write(indextomeanunit(i),rec=itm)CONVREAL(meansparsglobal(i)%mean2d)
                        endif
                        write(indextovarunit(i),rec=itm)CONVREAL(meansparsglobal(i)%variance2d)
                        where(meansparsglobal(i)%min2d>0.99d0*huge(0.d0))
                           meansparsglobal(i)%min2d=-999.d0
                        endwhere
                        where(meansparsglobal(i)%max2d<-0.99d0*huge(0.d0))
                           meansparsglobal(i)%max2d=-999.d0
                        endwhere
                        write(indextominunit(i),rec=itm)CONVREAL(meansparsglobal(i)%min2d)
                        write(indextomaxunit(i),rec=itm)CONVREAL(meansparsglobal(i)%max2d)
                     endif
                   case (3)
                     if(xomaster) then
                        write(indextomeanunit(i),rec=itm)CONVREAL(meansparsglobal(i)%mean3d)
                        write(indextovarunit(i),rec=itm)CONVREAL(meansparsglobal(i)%variance3d)
                        where(meansparsglobal(i)%min3d>0.99d0*huge(0.d0))
                           meansparsglobal(i)%min3d=-999.d0
                        endwhere
                        where(meansparsglobal(i)%max3d<-0.99d0*huge(0.d0))
                           meansparsglobal(i)%max3d=-999.d0
                        endwhere
                        write(indextominunit(i),rec=itm)CONVREAL(meansparsglobal(i)%min3d)
                        write(indextomaxunit(i),rec=itm)CONVREAL(meansparsglobal(i)%max3d)
                     endif
                   case (4)
                     if(xomaster) then
                        write(indextomeanunit(i),rec=itm)CONVREAL(meansparsglobal(i)%mean4d)
                        write(indextovarunit(i),rec=itm)CONVREAL(meansparsglobal(i)%variance4d)
                        where(meansparsglobal(i)%min4d>0.99d0*huge(0.d0))
                           meansparsglobal(i)%min4d=-999.d0
                        endwhere
                        where(meansparsglobal(i)%max4d<-0.99d0*huge(0.d0))
                           meansparsglobal(i)%max4d=-999.d0
                        endwhere
                        write(indextominunit(i),rec=itm)CONVREAL(meansparsglobal(i)%min4d)
                        write(indextomaxunit(i),rec=itm)CONVREAL(meansparsglobal(i)%max4d)
                     endif
                  end select
               enddo
            endif  ! t output
            par%tintm=tpar%tpm(min(itm+2,stpm))-tpar%tpm(itm+1)  ! Next averaging period (min to stop array out of bounds)
            par%tintm=max(par%tintm,tiny(0.d0))        ! to prevent par%tintm=0 after last output
         endif  ! nmeanvar > 0


         !!! Write global variables
         if (par%nglobalvar/=0) then
            if (tpar%outputg) then
               itg=itg+1
               do i = 1,par%nglobalvar
#ifdef USEMPI
                  call space_collect_index(s,sl,par,outnumbers(i))
#endif
                  if(xomaster) then
                     call indextos(s,outnumbers(i),t)
                     select case(t%type)
                      case ('r')
                        select case(t%rank)
                         case(0)
                           call outarray(i,t%r0)
                         case(1)
                           call outarray(i,t%r1)
                         case(2)
                           call outarray(s,i,t%r2)
                         case(3)
                           call outarray(s,i,t%r3)
                         case(4)
                           call outarray(i,t%r4)
                        end select
                      case('i')
                        select case(t%rank)
                         case(0)
                           call outarray(i,t%i0)
                         case(1)
                           call outarray(i,t%i1)
                         case(2)
                           call outarray(i,t%i2)
                         case(3)
                           call outarray(i,t%i3)
                         case(4)
                           call outarray(i,t%i4)
                        end select
                     end select
                  endif ! xomaster
               enddo   ! outnumber loop
            endif
         endif  ! end global file writing

         !
         ! write drifter output
         !

#ifdef USEMPI
         ! to send sl%idrift from xmaster to xomaster
         ! wwvv todo use xmpi_send for this
         call xmpi_bcast(sl%idrift,toall)
         ! wwvv todo shouldn't we send also sl%jdrift ?
         call xmpi_bcast(sl%jdrift,toall)
#endif
         if (xomaster) then

#ifdef USEMPI
            s%idrift = sl%idrift
            s%jdrift = sl%jdrift
#endif

            if (abs(mod(par%t,par%tintp))<1.d-6) then
               itd = itd+1
               do i=1,par%ndrifter
                  if (    par%t>=s%tdriftb(i) .and. par%t<=s%tdrifte(i) .and. &
                  s%idrift(i)>1       .and. s%idrift(i)<=s%nx   .and. &
                  s%jdrift(i)>1       .and. s%jdrift(i)<=s%ny             ) then

                     iz = int(s%idrift(i))
                     jz = int(s%jdrift(i))

                     di = mod(s%idrift(i),1.d0)
                     dj = mod(s%jdrift(i),1.d0)

                     dx = di*s%dsu(iz,jz)*cos(s%alfaz(iz,jz)) - &
                     dj*s%dnv(iz,jz)*sin(s%alfaz(iz,jz))
                     dy = di*s%dsu(iz,jz)*sin(s%alfaz(iz,jz)) + &
                     dj*s%dnv(iz,jz)*cos(s%alfaz(iz,jz))

                     write(indextodrifterunit(i),rec=itd)    &
                     CONVREAL(s%xz(iz,jz)+dx),                     &
                     CONVREAL(s%yz(iz,jz)+dy),                     &
                     CONVREAL(par%t)
                  else
                     write(indextodrifterunit(i),rec=itd)    &
                     CONVREAL(-999d0),                    &
                     CONVREAL(-999d0),                    &
                     CONVREAL(par%t)
                  endif
               enddo
            endif
         endif

         if(xomaster) then
            outputtimes=-999.d0
            outputtimes(1:itg)=tpar%tpg(1:itg)
            outputtimes(itg+1:itg+itp)=tpar%tpp(1:itp)
            outputtimes(itg+itp+1:itg+itp+itc)=tpar%tpc(1:itc)
            outputtimes(itg+itp+itc+1:itg+itp+itc+itm)=tpar%tpm(2:itm+1)          ! mean output always shifted by 1
            if (par%morfacopt==1) outputtimes=outputtimes*max(par%morfac,1.d0)
            open(999,file='dims.dat',form='unformatted',access='direct',recl=wordsize*(10+size(outputtimes)))
            write(999,rec=1) CONVREAL(itg*1.d0),&
            CONVREAL(s%nx*1.d0),&
            CONVREAL(s%ny*1.d0),&
            CONVREAL(s%ntheta*1.d0),&
            CONVREAL(par%kmax*1.d0),&
            CONVREAL(par%ngd*1.d0),&
            CONVREAL(par%nd*1.d0), &
            CONVREAL(tpar%itp*1.d0),&
            CONVREAL(itc*1.d0),&
            CONVREAL(itm*1.d0),&
            CONVREAL(outputtimes)
            call flush(999)
            ! Just output for MICORE for backwards compat
            !      open(999,file='dims.dat',form='unformatted',access='direct',recl=wordsize*(7+size(outputtimes)))
            !      write(999,rec=1)itg*1.d0,&
            ! s%nx*1.d0,&
            ! s%ny*1.d0,&
            !                      par%ngd*1.d0,&
            ! par%nd*1.d0, &
            !                      itp*1.d0,&
            ! itm*1.d0,&
            ! outputtimes
            ! close(999)
         endif  !xomaster

      end if  !

      !!! Close files

      if(xomaster) then
         if(par%t>=par%tstop) then
            do i=1,par%npoints+par%nrugauge
               close(indextopointsunit(i))
            enddo

            do i=1,par%nmeanvar
               close(indextomeanunit(i))
               close(indextovarunit(i))
               close(indextominunit(i))
               close(indextomaxunit(i))
            enddo

            do i=1,par%ndrifter
               close (indextodrifterunit(i))
            enddo

            do i=1,noutnumbers
               close (indextoglobalunit(i))
            enddo
         end if
      endif
      ! wwvv to avoid warning about unused xrank:
#ifdef USEMPI
      if (xrank .eq. -1) return
#endif

   end subroutine var_output

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!   INTERNAL SUBROUTINE   !!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



   subroutine makeintpvector(par,s,intpvector,mg,ng)

      use params
      use spaceparams
      use indextos_module
      use logging_module

      IMPLICIT NONE

      type(parameters), intent(in)   :: par
      type(spacepars), intent(in)    :: s    ! need slocal here, because
      ! the is,js lm and ln arrays
      ! on the not-master nodes
      ! are only in slocal
      real*8,dimension(:)            :: intpvector
      integer,intent(in)             :: mg,ng

      type(arraytype)                :: t
      integer                        :: j,m,n
      real*8                         :: value

#ifdef USEMPI
      integer                        :: i,p,ierr
      logical                        :: toall = .true.
#endif

      m = mg
      n = ng

#ifdef USEMPI
      !
      ! wwvv argh! for performance, we will not collect each
      ! matrix or block to master. Each process has to
      ! determine if it has the values, and then send them.
      !
      ! locate the process who has index m

      p=-1
      do i=1,xmpi_size
         !    if (mg .ge. s%is(i) .and. mg .le. s%is(i) + s%lm(i) .and. &
         !         ng .ge. s%js(i) .and. ng .le. s%js(i) + s%ln(i)) then
         ! Change suggested by Jim Gunson:
         if (mg .ge. s%is(i) .and. mg .lt. s%is(i) + s%lm(i) .and. &
         ng .ge. s%js(i) .and. ng .lt. s%js(i) + s%ln(i)) then
            !   if (mg .ge. s%icgs(i) .and. mg .le. s%icge(i) .and.  ng .ge. s%jcgs(i) .and. ng .le. s%jcge(i)) then
            ! wwvv-todo ?
            m = mg - s%is(i) + 1  ! m and n now the indices for process p
            n = ng - s%js(i) + 1
            p = i-1
            exit
         endif
      enddo
      if (p .lt. 0) then
         call writelog('els','','Catastrophic error in makeintpvector')
         call writelog('els','(a,i0,",",i0)','Cannot find processor for indexes ',mg,ng)
         call halt_program
      endif


      if (xmpi_rank .eq. p) then
#endif
      do j=1,numvars
         call indextos(s,j,t)
         value = -999   ! default
         select case (t%name)
          case(mnem_thetamean)
            !   value = 270-((s%thetamean(m,n)+s%alfaz(m,n))*(180/par%px))  ! iwwvv thetamean: is that
            value = 270-((s%thetamean(m,n))*(180/par%px))  ! iwwvv thetamean: is that
            ! different on each
            ! process?
          case(mnem_Fx)
            value = s%Fx(m,n)*cos(s%alfaz(m,n))-s%Fy(m,n)*sin(s%alfaz(m,n))
          case(mnem_Fy)
            value = s%Fx(m,n)*sin(s%alfaz(m,n))+s%Fy(m,n)*cos(s%alfaz(m,n))
          case(mnem_u)
            value = s%u(m,n)*cos(s%alfaz(m,n))-s%v(m,n)*sin(s%alfaz(m,n))
          case(mnem_gwu)
            value = s%gwu(m,n)*cos(s%alfaz(m,n))-s%gwv(m,n)*sin(s%alfaz(m,n))
          case(mnem_v)
            value = s%u(m,n)*sin(s%alfaz(m,n))+s%v(m,n)*cos(s%alfaz(m,n))
          case(mnem_gwv)
            value = s%gwu(m,n)*sin(s%alfaz(m,n))+s%gwv(m,n)*cos(s%alfaz(m,n))
          case(mnem_ue)
            value = s%ue(m,n)*cos(s%alfaz(m,n))-s%ve(m,n)*sin(s%alfaz(m,n))
          case(mnem_ve)
            value = s%ue(m,n)*sin(s%alfaz(m,n))+s%ve(m,n)*cos(s%alfaz(m,n))
          case(mnem_uwf)
            value = s%uwf(m,n)*cos(s%alfaz(m,n))-s%vwf(m,n)*sin(s%alfaz(m,n))
          case(mnem_vwf)
            value = s%uwf(m,n)*sin(s%alfaz(m,n))+s%vwf(m,n)*cos(s%alfaz(m,n))
          case(mnem_ui)
            value = s%ui(m,n)*cos(s%alfaz(m,n))-s%vi(m,n)*sin(s%alfaz(m,n))
          case(mnem_vi)
            value = s%ui(m,n)*sin(s%alfaz(m,n))+s%vi(m,n)*cos(s%alfaz(m,n))
          case(mnem_umean)
            value = s%umean(m,n)*cos(s%alfaz(m,n))-s%vmean(m,n)*sin(s%alfaz(m,n))
          case(mnem_vmean)
            value = s%umean(m,n)*sin(s%alfaz(m,n))+s%vmean(m,n)*cos(s%alfaz(m,n))
          case default
            select case(t%rank)
             case (0)
               if (t%type .eq.'i') then
                  value = dble(t%i0)
               else
                  value = t%r0
               endif
             case(1)
               if (t%type .eq. 'r') then
                  select case(t%name)
                   case default
                     continue
                     !                case (mnem_yz,mnem_yv)
                     !                   value = t%r1(n)
                     !                case (mnem_xz,mnem_xu)
                     !                   value = t%r1(m)
                  end select
               endif
             case (2)
               if (t%type .eq. 'i') then
                  if (m .le. size(t%i2,1) .and. n .le. size(t%i2,2)) then
                     value = dble(t%i2(m,n))
                  endif
               else
                  if (m .le. size(t%r2,1) .and. n .le. size(t%r2,2)) then
                     value = t%r2(m,n)
                  endif
               endif
            end select
         end select
         intpvector(j) = value
      enddo
#ifdef USEMPI
   endif
   !
   ! process p has now the intptvector (hopefully!)
   ! send it now to the master with a simpel send,
   ! master will receive this
   ! This is not necessary if master has everything
   !
   ! p is the rank in xmpi_comm
   ! what is the rank in xmpi_ocomm?
   ! for the time being, assume that that is
   ! p(xmpi_comm) + 1
   !
   p = p + 1
   if (p .ne. xmpi_omaster) then
      if (xmpi_orank .eq. p) then
         call MPI_Send(intpvector, numvars, MPI_DOUBLE_PRECISION, xmpi_omaster, &
         311, xmpi_ocomm,                    ierr)
      endif
      if (xomaster) then
         call MPI_Recv(intpvector, numvars, MPI_DOUBLE_PRECISION, p, &
         311, xmpi_ocomm, MPI_STATUS_IGNORE, ierr)
      endif
      ! this barrier is really needed:
      call xmpi_barrier(toall)
   endif
#endif
   end subroutine makeintpvector






   subroutine outarray_r0(index,x)
      implicit none
      integer                       :: index
      real*8,intent(in)             :: x
      integer                       :: unit,reclen,jtg

      inquire(iolength=reclen) CONVREAL(x)
      call checkfile(index,unit,reclen,jtg)

      write(unit,rec=jtg) CONVREAL(x)

   end subroutine outarray_r0

   subroutine outarray_r1(index,x)
      implicit none
      integer                       :: index
      real*8, dimension(:), pointer :: x
      integer                       :: unit,reclen,jtg
      character*20                  :: mnem
      real*8, parameter             :: pi = 4*atan(1.0d0)

      inquire(iolength=reclen) CONVREAL(x)
      call checkfile(index,unit,reclen,jtg)

      mnem = mnemonics(outnumbers(index))
      if (mnem .eq. mnem_theta .or. &
      mnem .eq. mnem_theta0 ) then
         write(unit,rec=jtg) CONVREAL(270-(x*(180/pi)))
      else
         write(unit,rec=jtg) CONVREAL(x)
      endif

   end subroutine outarray_r1

   subroutine outarray_r2(s,index,x)
      use spaceparams
      implicit none
      type(spacepars),intent(in)    :: s
      integer                       :: index
      real*8, dimension(:,:), pointer :: x
      integer                       :: unit,reclen,jtg
      character*20                  :: mnem
      real*8, parameter             :: pi = 4*atan(1.0d0)

      inquire(iolength=reclen) CONVREAL(x)
      call checkfile(index,unit,reclen,jtg)

      mnem = mnemonics(outnumbers(index))

      select case(mnem)
       case(mnem_thetamean)
         write(unit,rec=jtg)CONVREAL(270-((x)*(180/pi)))
       case(mnem_Fx)
         write(unit,rec=jtg)CONVREAL(x*cos(s%alfaz)-s%Fy*sin(s%alfaz))
       case(mnem_Fy)
         write(unit,rec=jtg)CONVREAL(s%Fx*sin(s%alfaz)+x*cos(s%alfaz))
       case(mnem_u)
         write(unit,rec=jtg)CONVREAL(x*cos(s%alfaz)-s%v*sin(s%alfaz))
       case(mnem_gwu)
         write(unit,rec=jtg)CONVREAL(x*cos(s%alfaz)-s%gwv*sin(s%alfaz))
       case(mnem_v)
         write(unit,rec=jtg)CONVREAL(s%u*sin(s%alfaz)+x*cos(s%alfaz))
       case(mnem_gwv)
         write(unit,rec=jtg)CONVREAL(s%gwu*sin(s%alfaz)+x*cos(s%alfaz))
       case(mnem_ue)
         write(unit,rec=jtg)CONVREAL(x*cos(s%alfaz)-s%ve*sin(s%alfaz))
       case(mnem_ve)
         write(unit,rec=jtg)CONVREAL(s%ue*sin(s%alfaz)+x*cos(s%alfaz))
       case(mnem_ui)
         write(unit,rec=jtg)CONVREAL(x*cos(s%alfaz)-s%vi*sin(s%alfaz))
       case(mnem_vi)
         write(unit,rec=jtg)CONVREAL(s%ui*sin(s%alfaz)+x*cos(s%alfaz))
       case(mnem_umean)
         write(unit,rec=jtg)CONVREAL(x*cos(s%alfaz)-s%vmean*sin(s%alfaz))
       case(mnem_vmean)
         write(unit,rec=jtg)CONVREAL(s%umean*sin(s%alfaz)+x*cos(s%alfaz))
       case(mnem_uwf)
         write(unit,rec=jtg)CONVREAL(x*cos(s%alfaz)-s%vwf*sin(s%alfaz))
       case(mnem_vwf)
         write(unit,rec=jtg)CONVREAL(s%uwf*sin(s%alfaz)+x*cos(s%alfaz))
       case(mnem_Sutot)
         write(unit,rec=jtg)CONVREAL((sum(s%Subg,DIM=3)+sum(s%Susg,DIM=3))*cos(s%alfaz) - (sum(s%Svbg,DIM=3)+sum(s%Svsg,DIM=3))*sin(s%alfaz))
       case(mnem_Svtot)
         write(unit,rec=jtg)CONVREAL((sum(s%Subg,DIM=3)+sum(s%Susg,DIM=3))*sin(s%alfaz) + (sum(s%Svbg,DIM=3)+sum(s%Svsg,DIM=3))*cos(s%alfaz))
       case(mnem_cctot)
         write(unit,rec=jtg)CONVREAL(sum(s%ccg,DIM=3))
       case default
         write(unit,rec=jtg) CONVREAL(x)
      end select

   end subroutine outarray_r2

   subroutine outarray_r3(s,index,x)
      use spaceparams
      implicit none
      type(spacepars),intent(in)    :: s
      integer index
      real*8, dimension(:,:,:), pointer :: x
      integer                       :: unit,reclen,jtg
      character*(dimnamelen)        :: mnem
      real*8,parameter              :: pi = 4*atan(1.0d0)

      inquire(iolength=reclen) CONVREAL(x)
      call checkfile(index,unit,reclen,jtg)

      mnem = mnemonics(outnumbers(index))

      !!
      !! Dano: need to find elegant way to multiply 3d array woth 2d alfaz
      !  select case(mnem)
      !  case(mnem_cgx)
      !     write(unit,rec=jtg)x*cos(s%alfaz)-s%cgy*sin(s%alfaz)
      !  case(mnem_cgy)
      !     write(unit,rec=jtg)s%cgx*sin(s%alfaz)+x*cos(s%alfaz)
      !  case(mnem_cx)
      !     write(unit,rec=jtg)x*cos(s%alfaz)-s%cy*sin(s%alfaz)
      !  case(mnem_cy)
      !     write(unit,rec=jtg)s%cx*sin(s%alfaz)+x*cos(s%alfaz)
      !  case(mnem_thet)
      !     write(unit,rec=jtg)270-((s%thet+s%alfaz)*(180/pi))
      !  case(mnem_Susg)
      !     write(unit,rec=jtg)x*cos(s%alfaz)-s%Svsg*sin(s%alfaz)
      !  case(mnem_Svsg)
      !     write(unit,rec=jtg)s%Susg*sin(s%alfaz)+x*cos(s%alfaz)
      !  case(mnem_Subg)
      !     write(unit,rec=jtg)x*cos(s%alfaz)-s%Svbg*sin(s%alfaz)
      !  case(mnem_Svbg)
      !     write(unit,rec=jtg)s%Subg*sin(s%alfaz)+x*cos(s%alfaz)
      !  case default
      write(unit,rec=jtg) CONVREAL(x)
      !  end select
      ! wwvv to avoid warning about unused parameter s:
      if (s%nx .eq. -1) return
   end subroutine outarray_r3

   subroutine outarray_r4(index,x)
      implicit none
      integer index
      real*8, dimension(:,:,:,:), pointer :: x
      integer                       :: unit,reclen,jtg

      inquire(iolength=reclen) CONVREAL(x)
      call checkfile(index,unit,reclen,jtg)
      write(unit,rec=jtg) CONVREAL(x)
   end subroutine outarray_r4

   subroutine outarray_i0(index,x)
      implicit none
      integer                       :: index
      integer,intent(in)            :: x
      integer                       :: unit,reclen,jtg

      inquire(iolength=reclen) x
      call checkfile(index,unit,reclen,jtg)

      write(unit,rec=jtg) x

   end subroutine outarray_i0

   subroutine outarray_i1(index,x)
      implicit none
      integer index
      integer, dimension(:), pointer :: x
      integer                       :: unit,reclen,jtg

      inquire(iolength=reclen) x
      call checkfile(index,unit,reclen,jtg)
      write(unit,rec=jtg) x
   end subroutine outarray_i1

   subroutine outarray_i2(index,x)
      implicit none
      integer index
      integer, dimension(:,:), pointer :: x
      integer                       :: unit,reclen,jtg

      inquire(iolength=reclen) x
      call checkfile(index,unit,reclen,jtg)
      write(unit,rec=jtg) x
   end subroutine outarray_i2

   subroutine outarray_i3(index,x)
      implicit none
      integer index
      integer, dimension(:,:,:), pointer :: x
      integer                       :: unit,reclen,jtg

      inquire(iolength=reclen) x
      call checkfile(index,unit,reclen,jtg)
      write(unit,rec=jtg) x
   end subroutine outarray_i3

   subroutine outarray_i4(index,x)
      implicit none
      integer index
      integer, dimension(:,:,:,:), pointer :: x
      integer                       :: unit,reclen,jtg

      inquire(iolength=reclen) x
      call checkfile(index,unit,reclen,jtg)
      write(unit,rec=jtg) x
   end subroutine outarray_i4

   subroutine checkfile(index,unit,reclen,jtg)
      implicit none
      integer, intent(in)  :: index,reclen
      integer, intent(out) :: unit,jtg
      logical              :: lopen
      character(len=1000)  :: filename

      unit = indextoglobalunit(index)
      inquire(unit=unit, opened=lopen)
      if ( .not. lopen ) then
         filename = trim(mnemonics(outnumbers(index)))//'.dat'
         open(unit, file=trim(filename),form='unformatted',&
         access='direct',recl=reclen)
      endif
      inquire(unit=unit,nextrec=jtg)
   end subroutine checkfile

   integer function indextoglobalunit(index)
      implicit none
      integer, intent(in) :: index
      indextoglobalunit = 100+index
   end function indextoglobalunit

   integer function indextomeanunit(index)
      implicit none
      integer, intent(in) :: index
      indextomeanunit = 100+numvars+index
   end function indextomeanunit

   integer function indextopointsunit(index)
      implicit none
      integer, intent(in) :: index
      indextopointsunit = 100+20*numvars+index
   end function indextopointsunit

   integer function indextocrossunit(index)
      implicit none
      integer, intent(in) :: index
      indextocrossunit = 100+30*numvars+index
   end function indextocrossunit

   integer function indextominunit(index)
      implicit none
      integer, intent(in) :: index
      indextominunit = 100+40*numvars+index
   end function indextominunit

   integer function indextomaxunit(index)
      implicit none
      integer, intent(in) :: index
      indextomaxunit = 100+50*numvars+index
   end function indextomaxunit

   integer function indextovarunit(index)
      implicit none
      integer, intent(in) :: index
      indextovarunit = 100+60*numvars+index
   end function indextovarunit

   integer function indextodrifterunit(index)
      implicit none
      integer, intent(in) :: index
      indextodrifterunit = 700+index
   end function indextodrifterunit

end module fortoutput_module
! matching #if 0 at the start:
#endif
