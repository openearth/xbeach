module waveparams
   implicit none
   save
   type waveparameters

      integer                                 :: K, Npy, Nr
      integer*4                               :: listline
      integer, dimension(:), pointer          :: index_vector

      real*8                                  :: mainang,dang, scoeff                    !scoeff is now a wp
      real*8                                  :: h0t0
      real*8                                  :: hm0gew, df
      real*8                                  :: Ly, dt, rt
      real*8,dimension(:),pointer             :: S0, dthetafin, fgen, theta0
      real*8,dimension(:),pointer             :: Sf, Dd, f, theta, window
      real*8,dimension(:,:),pointer           :: S_array

      double complex,dimension(:,:),pointer   :: CompFn

   end type waveparameters

contains

   ! --------------------------------------------------------------
   ! ------------- Sorting for calling functions ------------------
   ! --------------------------------------------------------------

   !> Creates energy and flux boundary condition files (*.bcf).
   !! @param   par struct with values from params.txt
   !! @param   s   struct with grid information
   !! @param   wp  struct with wave information
   !! @return  void
   !! @todo    clean up
   subroutine makebcf(par,s,wp)

      use params
      use spaceparams
      use readkey_module
      use xmpi_module
      use typesandkinds
      IMPLICIT NONE

      ! Input / output variables
      type(parameters), INTENT(INOUT)             :: par
      type(spacepars), intent(IN)                 :: s

      type(waveparameters), intent(inout)         :: wp
      real*8,save                                 :: bcendtime
      character(len=slen)                         :: fname,Ebcfname,qbcfname
      character*8                                 :: testc
      logical                                     :: makefile
      integer,save                                :: reuse  ! = 0 used to be in code
      integer,save                                :: counter

      ! Flag that determines whether new BCF files are created or not
      makefile=.false.

      ! Check whether this is the first time step
      if (abs(par%t-par%dt)<1.d-6) then
         bcendtime=0
         wp%listline=0
         counter=0
         wp%scoeff=-1

         if(xmaster) then
            open(74,file=par%bcfile,form='formatted')
            if (par%wbctype/=WBCTYPE_JONS_TABLE) read(74,*)testc
         endif

         if(xmaster) then
            ! Create new BCF files, if in master process
            open(53,file='ebcflist.bcf',form='formatted',status='replace')
            open(54,file='qbcflist.bcf',form='formatted',status='replace')
            close(53)
            close(54)
         endif

         ! Check whether BCF files should be reused
         if (testc=='FILELIST' .or. par%wbctype==WBCTYPE_JONS_TABLE) then
            reuse=0
         else
            reuse=1
            if(xmaster) then
               close(74) ! only close if this is not the list of files
            endif
         end if
      end if

      ! Check whether this is the last time step
      if (par%t>(par%tstop-par%dt)) then
         ! Prevent recalculation of boundary conditions for last timestep
         return
         close(74)
      end if

      ! wp%listline is not increased, therefore, first line of current bcf
      ! is used (i.e. all zeros)

      ! Check whether BCF files should be reused
      if (reuse==0) then

         ! Read rt and dt values
         if(xmaster) then
            if (par%wbctype/=WBCTYPE_JONS_TABLE) then
               read(74,*)wp%rt,wp%dt,fname
               if (par%morfacopt==1) wp%rt = wp%rt / max(par%morfac,1.d0)
            endif
         endif

#ifdef USEMPI
         !Dano call xmpi_bcast(wp%rt)
         !Dano call xmpi_bcast(wp%dt)
         !Dano call xmpi_bcast(fname)
#endif

         ! Create filenames for BCF files
         if (par%wbctype/=WBCTYPE_JONS_TABLE) then
            Ebcfname='E_'//trim(fname)
            qbcfname='q_'//trim(fname)
         else
            counter=counter+1
            write(Ebcfname, '(A,I0.6,A)') 'Ejonsw', counter, '.bcf'
            write(qbcfname, '(A,I0.6,A)') 'qjonsw', counter, '.bcf'
         endif
      else

         ! Read rt and dt values in first timestep
         if (abs(par%t-par%dt)<1.d-6) then
            wp%rt = readkey_dbl ('params.txt', 'rt', 3600.d0, 1200.d0, 7200.d0, bcast=.false.)
            if (par%morfacopt==1) wp%rt = wp%rt / max(par%morfac,1.d0)
            wp%dt = readkey_dbl ('params.txt', 'dtbc', 0.1d0, 0.01d0, 1.0d0, bcast=.false.)
         end if

         ! Create filenames for BCF files
         fname=par%bcfile
         Ebcfname='E_reuse.bcf'
         qbcfname='q_reuse.bcf'
      end if

      ! (Re)create BCF files if this is the first time step or it is explicitly
      ! requested
      if (abs(par%t-par%dt)<1.d-6) then
         makefile=.true.
      else
         if (reuse==0) then
            makefile=.true.
         end if
      end if

      ! (Re)create BCF files, if requested
      if (makefile) then

         ! Calculate weighed average water depth at offshore boundary
         if (s%ny>0) then
            wp%h0t0=sum(s%hh(1,2:s%ny)*s%dnv(1,2:s%ny))/sum(s%dnv(1,2:s%ny))
         else
            wp%h0t0 = s%hh(1,1)
         endif

         if (par%wbctype==WBCTYPE_PARAMETRIC .or. par%wbctype==WBCTYPE_JONS_TABLE) then
            ! Use JONSWAP spectrum
            call build_jonswap(par,s,wp,fname)
            call build_etdir(par,s,wp,Ebcfname)
            call build_boundw(par,s,wp,qbcfname)
         elseif (par%wbctype==WBCTYPE_SWAN) then
            ! Use SWAN data
            call swanreader(par,s,wp,fname)
            call build_etdir(par,s,wp,Ebcfname)
            call build_boundw(par,s,wp,qbcfname)
         elseif (par%wbctype==WBCTYPE_VARDENS) then
            ! Use variance density spectrum
            call vardensreader(par,s,wp,fname)
            call build_etdir(par,s,wp,Ebcfname)
            call build_boundw(par,s,wp,qbcfname)
         endif
      end if

      ! Define line counter for boundaryconditions.f90
      wp%listline=wp%listline+1

      ! Define endtime for boundary conditions, boundaryconditions.f90 should
      ! recreate BCF files after this moment
      bcendtime=bcendtime+wp%rt

      if(xmaster) then

         ! Write lists with BCF information
         open(53,file='ebcflist.bcf',form='formatted',position='append')
         open(54,file='qbcflist.bcf',form='formatted',position='append')
         write(53,'(f12.3,a,f12.3,a,f9.3,a,f9.5,a,f11.5,a)') &
         & bcendtime,' ',wp%rt,' ',wp%dt,' ',par%Trep,' ',wp%mainang,' '//trim(Ebcfname)
         write(54,'(f12.3,a,f12.3,a,f9.3,a,f9.5,a,f11.5,a)') &
         & bcendtime,' ',wp%rt,' ',wp%dt,' ',par%Trep,' ',wp%mainang,' '//trim(qbcfname)
         close(53)
         close(54)
      endif

   end subroutine makebcf

   ! --------------------------------------------------------------
   ! ----------------------Read spectrum files --------------------
   ! --------------------------------------------------------------
   subroutine build_jonswap(par,s,wp,fname)

      use readkey_module
      use params
      use spaceparams
      use xmpi_module
      use logging_module

      IMPLICIT NONE

      ! Input / output variables
      type(parameters), INTENT(INout)         :: par
      type(spacepars), intent(IN)             :: s
      type(waveparameters), INTENT(INOUT)     :: wp
      character(len=*)                        :: fname

      ! Internal variables
      integer                                 :: i=0,ii,nang,nfreq,ier
      integer                                 :: firstp, lastp
      real*8,dimension(:),allocatable         :: temp, x, y
      real*8                                  :: dfj, fnyq, fp
      real*8                                  :: gam
      character(len=80)                       :: dummystring

      fnyq = -123
      dfj  = -123
      ! Read JONSWAP data
      if(xmaster) then

         ! Check whether spectrum characteristics or table should be used
         if (par%wbctype /= WBCTYPE_JONS_TABLE) then

            ! Use spectrum characteristics
            call writelog('sl','','waveparams: Reading from ',trim(fname),' ...')
         else

            ! Use spectrum table
            call readkey('params.txt','bcfile',fname)
            call writelog('sl','','waveparams: Reading from table',fname,' ...')
            read(74,*,iostat=ier)wp%hm0gew,fp,wp%mainang,gam,wp%scoeff,wp%rt,wp%dt
            if (par%morfacopt==1) wp%rt = wp%rt/max(par%morfac,1.d0)

            ! Set extra parameters
            fp=1.d0/fp
            fnyq=3.d0*fp
            dfj=fp/20.d0
         endif
      endif

      ! Read JONSWAP characteristics
      if (par%wbctype/=WBCTYPE_JONS_TABLE) then
         !                                      Input file   Keyword     Default     Minimum     Maximum
         wp%hm0gew           = readkey_dbl (fname,       'Hm0',      0.0d0,      0.00d0,     5.0d0,      bcast=.false. )
         fp                  = readkey_dbl (fname,       'fp',       0.08d0,     0.0625d0,   0.4d0,      bcast=.false. )
         fnyq                = readkey_dbl (fname,       'fnyq',     max(0.3d0,3.d0*fp),    0.2d0,      1.0d0, bcast=.false. )
         dfj                 = readkey_dbl (fname,       'dfj',      fnyq/200,   fnyq/1000,  fnyq/20,    bcast=.false. )
         gam                 = readkey_dbl (fname,       'gammajsp', 3.3d0,      1.0d0,      5.0d0,      bcast=.false. )
         wp%scoeff           = readkey_dbl (fname,       's',        10.0d0,     1.0d0,      1000.0d0,   bcast=.false. )
         wp%mainang          = readkey_dbl (fname,       'mainang',  270.0d0,    0.0d0,      360.0d0,    bcast=.false. )
         !
         if(xmaster) then
            call readkey(fname,'checkparams',dummystring)
         endif
      endif

      ! Determine number of nodes in y-direction
      wp%Npy=s%ny+1

      ! Print JONSWAP characteristics to screen
      call writelog('sl','(a,f0.3,a,f0.3,a,f0.3,a,f0.3)','Input checked: Hm0 = ',wp%hm0gew,' Tp = ',1.d0/fp, &
      ' dir = ',wp%mainang,' duration = ',wp%rt)

      ! approximation from Coastal Engineering: Processes, Theory and Design Practice
      ! Dominic Reeve, Andrew Chadwick 2004
      ! par%Trep=0.8345d0*(1/fp)

      ! Make the sample resoltion depending on the time record length, thus
      ! discarding the input parameter dfj. The base frequency now perfectly fits the
      ! given time record and the calculation of the number of wave components will
      ! result in an integer value.
      ! dfj = 1/wp%rt

      ! Define number of frequency bins by defining an array of the necessary length
      ! using the Nyquist frequency and frequency step size
      allocate(temp(ceiling((fnyq-dfj)/dfj)))
      temp=(/(i,i=1,size(temp))/)

      ! Define array with actual eqidistant frequency bins
      allocate(wp%f(size(temp)))
      wp%f=temp*dfj
      deallocate (temp)

      ! Determine frequency bins relative to peak frequency
      allocate(x(size(wp%f)))
      x=wp%f/fp

      ! Calculate unscaled and non-directional JONSWAP spectrum using
      ! peak-enhancement factor and pre-determined frequency bins
      allocate(y(size(wp%f)))
      call jonswapgk(x,gam,y)

      ! Determine scaled and non-directional JONSWAP spectrum using the JONSWAP
      ! characteristics
      y=(wp%hm0gew/(4.d0*sqrt(sum(y)*dfj)))**2*y
      deallocate (x)

      ! Define 200 directions relative to main angle running from -pi to pi
      allocate(temp(201))
      allocate(wp%theta(201))
      temp=(/(i,i=0,200)/)
      wp%theta=temp*((2*par%px)/200.d0)-par%px
      deallocate (temp)

      ! Define pi/2
      !t1=-(par%px)/2.d0

      ! Define 100 directions relative to main angle running from -pi/2 to pi/2
      !allocate(temp(101))
      !allocate(wp%theta(101))
      !temp=(/(i,i=0,100)/)
      !wp%theta=temp*((par%px)/100.d0)+t1
      !deallocate (temp)

      ! Determine directional step size: pi/200
      wp%dang=wp%theta(2)-wp%theta(1)

      ! Define 200 directional bins accordingly
      allocate (wp%Dd(size(wp%theta)))

      ! Convert main angle from degrees to radians and from nautical convention to
      ! internal grid
      wp%mainang=(1.5d0*par%px)-wp%mainang*par%px/180.d0

      ! Make sure the main angle is defined between 0 and 2*pi
      !do while (wp%mainang>2*par%px .or. wp%mainang<0.d0)
      !    if (wp%mainang>2*par%px) then
      !        wp%mainang=wp%mainang-2*par%px
      !    elseif (wp%mainang<0.d0*par%px) then
      !        wp%mainang=wp%mainang+2*par%px
      !    endif
      !enddo

      do while (wp%mainang>par%px .or. wp%mainang<-par%px) !Robert en Ap
         if (wp%mainang>par%px) then
            wp%mainang=wp%mainang-2*par%px
         elseif (wp%mainang<-par%px) then
            wp%mainang=wp%mainang+2*par%px
         endif
      enddo

      ! Convert 200 directions relative to main angle to directions relative to
      ! internal grid                                                                ! Bas: apparently division by 2 for cosine law happens already here
      allocate(temp(size(wp%theta)))
      temp = (wp%theta-wp%mainang)/2.d0

      ! Make sure all directions around the main angle are defined between 0 and 2*pi
      do while (any(temp>2*par%px) .or. any(temp<0.d0))
         where (temp>2*par%px)
            temp=temp-2*par%px
         elsewhere (temp<0.d0*par%px)
            temp=temp+2*par%px
         endwhere
      enddo

      ! Calculate directional spreading based on cosine law
      wp%Dd = dcos(temp)**(2*nint(wp%scoeff))                                            ! Robert: apparently nint is needed here, else MATH domain error
      deallocate(temp)

      ! Scale directional spreading to have a surface of unity by dividing by it's
      ! own surface
      wp%Dd = wp%Dd / (sum(wp%Dd)*wp%dang)

      ! Define number of directional and frequency bins
      nang=size(wp%theta)
      nfreq=size(y)

      ! Define two-dimensional variance density spectrum array and distribute
      ! variance density for each frequency over directional bins
      allocate(wp%S_array(nfreq,nang))
      do i=1,nang
         do ii=1,nfreq
            wp%S_array(ii,i)=y(ii)*wp%Dd(i)
         end do
      end do
      deallocate (y)

      ! Back integrate two-dimensional variance density spectrum over directions
      allocate(wp%Sf(size(wp%f)))
      wp%Sf = sum(wp%S_array, DIM = 2)*wp%dang

      ! Calculate mean wave period based on one-dimensional non-directional variance
      ! density spectrum and factor trepfac
      call tpDcalc(wp%Sf,wp%f,par%Trep,par%trepfac,par%Tm01switch)
      call writelog('sl','(a,f0.3)','Derived Trep = ',par%Trep)
      ! par%Trep=1.d0/par%Trep
      ! Jaap try Tm-1,0
      ! par%Trep = 1/fp/1.1d0

      ! Determine frequencies around peak frequency of one-dimensional
      ! non-directional variance density spectrum, based on factor sprdthr, which
      ! should be included in the determination of the wave boundary conditions
      allocate(temp(size(wp%Sf)))
      temp=0.d0
      call frange(par,wp%Sf,firstp,lastp,temp)
      deallocate (temp)

      ! Calculate number of wave components to be included in determination of the
      ! wave boundary conditions based on the wave record length and width of the
      ! wave frequency range
      wp%K=ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1)                               ! ja/ap: changed wp%K=max(100,ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1))
      !wp%K=max(256*s%ntheta,wp%K)

      return

   end subroutine build_jonswap

   ! --------------------------------------------------------------
   ! ----------------------Read SWAN files ------------------------
   ! --------------------------------------------------------------
   subroutine swanreader(par,s,wp,fname)

      use params
      use spaceparams
      use readkey_module
      use math_tools
      use xmpi_module
      use logging_module

      IMPLICIT NONE

      ! Input / output variables
      type(parameters), INTENT(INout)         :: par
      type(spacepars), intent(IN)             :: s
      type(waveparameters), INTENT(INOUT)     :: wp
      character(len=*), INTENT(IN)            :: fname

      ! Internal variables
      character(6)                            :: rtext
      real*8                                  :: factor,exc,m0,dthetaS_XB
      integer                                 :: nfreq, ndir, switch, i, flipped
      integer                                 :: firstp,lastp,nt,Ashift
      real*8, dimension(:),allocatable        :: temp, findline
      real*8, dimension(:,:),allocatable      :: tempA

      dthetaS_XB = readkey_dbl ('params.txt','dthetaS_XB', 0.0d0, -360.0d0, 360.0d0,bcast=.false.)

      flipped=0
      wp%Npy=s%ny+1

      switch = 0
      if(xmaster) then
         call writelog('sl','','Reading from SWAN file: ',fname,' ...')
         open(44,file=fname,form='formatted',status='old')


         ! Read file until RFREQ or AFREQ is found

         do while (switch==0)
            read(44,'(a)')rtext
            if (rtext == 'RFREQ ') then
               switch = 1
            elseif (rtext == 'AFREQ ') then
               switch = 2
            end if
         end do

         ! Read nfreq and f

         read(44,*)nfreq
      endif
#ifdef USEMPI
      !Dano call xmpi_bcast(nfreq)
      !Dano call xmpi_bcast(switch)
#endif
      allocate(wp%f(nfreq))
      if(xmaster) then
         do i=1,nfreq
            read(44,*)wp%f(i)
         end do
      endif
#ifdef USEMPI
      !Dano call xmpi_bcast(wp%f)
#endif

      ! Convert to absolute frequencies

      if (switch == 1) then
         wp%f = wp%f
      else
         wp%f = wp%f
      end if

      ! Read CDIR or NDIR

      if(xmaster) then
         read(44,'(a)')rtext
         if (rtext == 'NDIR  ') then
            switch = 1
         elseif (rtext == 'CDIR  ') then
            switch = 2
         else
            call writelog('els','', 'SWAN directional bins keyword not found')
#ifdef USEMPI
            call halt_program
#else
            stop
#endif
         end if

         ! Read ndir, theta

         read(44,*)ndir
      endif
#ifdef USEMPI
      !Dano call xmpi_bcast(ndir)
#endif
      allocate(wp%theta(ndir))

      if(xmaster) then
         do i=1,ndir
            read(44,*)wp%theta(i)
         end do
      endif

#ifdef USEMPI
      !Dano call xmpi_bcast(wp%theta)
#endif

      ! Convert angles to XBeach angles and radians

      if (switch == 1) then
         wp%theta = 270-wp%theta
      else
         wp%theta = wp%theta-dthetaS_XB                  ! dthetaS_XB is the angle in the degrees to rotate the x-axis in SWAN to the
         ! x-axis in XBeach (in Cartesian terms) (Have Fun :-))
      end if

      ! Ensure angles are increasing instead of decreasing
      if ((wp%theta(2)-wp%theta(1))<0) then
         call flipv(wp%theta,size(wp%theta))
         flipped=1
      end if

      nt = 0
      Ashift = 0
      ! Make sure that all angles are in -180 to 180 degrees
      if(minval(wp%theta)<-180)then
         allocate (temp(ndir))
         Ashift=-1
         temp=0
         do i=1,ndir
            if (wp%theta(i)<-180) then
               wp%theta(i)=wp%theta(i)+360.0d0
               nt = nt+1
            endif
         enddo
         !temp(1:ndir-nt) = wp%theta(nt+1:ndir)
         !temp(ndir-nt+1:ndir) = wp%theta(1:nt)
         !  temp(1:nt) = wp%theta(ndir-nt+1:ndir)
         !  temp(nt+1:ndir) = wp%theta(1:ndir-nt)
         temp(1:ndir-nt)=wp%theta(nt+1:ndir)
         temp(ndir-nt+1:ndir)=wp%theta(1:nt)
         wp%theta=temp
         deallocate(temp)
      elseif(maxval(wp%theta)>180.0d0)then
         allocate(temp(ndir))
         Ashift=1
         temp=0
         do i=1,ndir
            if (wp%theta(i)>180.0d0) then
               wp%theta(i)=wp%theta(i)-360.0d0
               nt = nt+1
            endif
         enddo
         !temp(1:ndir-nt) = wp%theta(nt+1:ndir)
         !temp(ndir-nt+1:ndir) = wp%theta(1:nt)
         !temp(1:nt) = wp%theta(ndir-nt+1:ndir)
         !temp(nt+1:ndir) = wp%theta(1:ndir-nt)
         temp(nt+1:ndir)=wp%theta(1:ndir-nt)
         temp(1:nt)=wp%theta(ndir-nt+1:ndir)
         wp%theta=temp
         deallocate(temp)
      endif

      wp%theta=wp%theta*par%px/180
      wp%dang=wp%theta(2)-wp%theta(1)

      ! Skip Quant, next line, read VaDens or EnDens
      if(xmaster) then
         read(44,'(a)')rtext
         read(44,'(a)')rtext
         read(44,'(a)')rtext
         if (rtext == 'VaDens') then
            switch = 1
         elseif (rtext == 'EnDens') then
            switch = 2
         else
            call writelog('sle','', 'SWAN VaDens/EnDens keyword not found')
#ifdef USEMPI
            call halt_program
#else
            stop
#endif
         end if
      endif
#ifdef USEMPI
      !Dano call xmpi_bcast(switch)
#endif

      if(xmaster) then
         read(44,'(a)')rtext
         read(44,*)exc
      endif

#ifdef USEMPI
      !Dano call xmpi_bcast(rtext)
      !Dano call xmpi_bcast(exc)
#endif

      if(xmaster) then
         i=0
         ! Find FACTOR keyword
         do while (i==0)
            read(44,'(a)')rtext
            if (rtext == 'FACTOR') then
               i=1
            elseif (rtext == 'ZERO  ') then
               call writelog('lse','','Zero energy density input for this point')
#ifdef USEMPI
               call halt_program
#else
               stop
#endif
            elseif (rtext == 'NODATA') then
               call writelog('lse','','SWAN file has no data for this point')
#ifdef USEMPI
               call halt_program
#else
               stop
#endif
            end if
         end do
         read(44,*)factor
      endif

#ifdef USEMPI
      !Dano call xmpi_bcast(factor)
#endif
      ! Read S_array
      allocate(wp%S_array(nfreq,ndir))

      if(xmaster) then
         do i=1,nfreq
            read(44,*)wp%S_array(i,:)
         end do
      endif

#ifdef USEMPI
      !Dano call xmpi_bcast(wp%S_array)
#endif

      where (wp%S_array == exc)
         wp%S_array =0
      endwhere

      ! If angles were decreasing, flip S_array as also dir is flipped
      if (flipped == 1) then
         flipped=2
         call flipa(wp%S_array,nfreq,ndir,flipped)
      end if


      ! Make sure that all wave variance is between -180 to 180 degrees range
      if(Ashift==-1)then
         allocate(tempA(nfreq,ndir))
         tempA=0
         !   tempA(:,ndir-nt+1:ndir) = wp%S_array(:,1:nt)
         !   tempA(:,1:ndir-nt) = wp%S_array(:,nt+1:ndir)
         !    tempA(:,1:nt) = wp%S_array(:,ndir-nt+1:ndir)
         !    tempA(:,nt+1:ndir) = wp%S_array(:,1:ndir-nt)
         tempA(:,1:ndir-nt)=wp%S_array(:,nt+1:ndir)
         tempA(:,ndir-nt+1:ndir)=wp%S_array(:,1:nt)
         wp%S_array=tempA
         deallocate(tempA)
      elseif (Ashift==1) then
         allocate(tempA(nfreq,ndir))
         tempA=0
         !   tempA(:,ndir-nt+1:ndir) = wp%S_array(:,1:nt)
         !   tempA(:,1:ndir-nt) = wp%S_array(:,nt+1:ndir)
         !    tempA(:,1:nt) = wp%S_array(:,ndir-nt+1:ndir)
         !    tempA(:,nt+1:ndir) = wp%S_array(:,1:ndir-nt)
         !        tempA(:,nt+1:ndir)=wp%S_array(:,1:ndir-nt)
         !    tempA(:,ndir-nt+1:ndir)=wp%S_array(:,1:nt)
         tempA(:,nt+1:ndir)=wp%S_array(:,1:ndir-nt)
         tempA(:,1:nt)=wp%S_array(:,ndir-nt+1:ndir)
         wp%S_array=tempA
         deallocate(tempA)
      endif

      wp%S_array=wp%S_array*factor

      if(xmaster) then
         close(44)                               ! Finished reading file
      endif

      ! Convert to m2/Hz/rad

      wp%S_array=wp%S_array*180/par%px

      ! Convert from energy density to variance density

      if (switch == 2) then
         wp%S_array=wp%S_array/(par%rho*par%g)
      end if

      allocate(wp%Sf(nfreq))
      wp%Sf = sum(wp%S_array, DIM = 2)*wp%dang

      ! Find main wave direction
      allocate (temp(ndir))
      temp=sum(wp%S_array, DIM = 1)
      i=maxval(maxloc(temp))
      wp%mainang=wp%theta(i)
      deallocate(temp)

      allocate (temp(nfreq+1))
      temp(1)=0
      temp(2:nfreq)=0.5*wp%f(1:nfreq-1)+0.5*wp%f(2:nfreq)
      temp(nfreq+1)=wp%f(nfreq)
      ! Calculate zero-order moment
      m0=0
      !m1=0
      do i=1,nfreq
         m0=m0+wp%Sf(i)*(temp(i+1)-temp(i))
         !    m1=m1+wp%f(i)*wp%Sf(i)*(temp(i+1)-temp(i))
      end do
      deallocate (temp)

      wp%hm0gew=4*sqrt(m0)
      !par%Trep=1/(m1/m0)

      call tpDcalc(wp%Sf,wp%f,par%Trep,par%trepfac,par%Tm01switch)
      call writelog('sl','','Trep = ',par%Trep,'s')
      call writelog('sl','','Hm0  = ',wp%hm0gew,'m')
      call writelog('sl','','Peak dir  = ',wp%mainang*180/par%px,' Cartesian degrees relative to x-axis')
      call writelog('sl','','Peak dir  = ',270-wp%mainang*180/par%px,' Nautical degrees')

      allocate(findline(size(wp%Sf)))

      firstp=0
      lastp=0
      call frange(par, wp%Sf,firstp,lastp,findline)
      deallocate(findline)


      !!!!! ja/ap wp%K=max(100,ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1))

      wp%K=ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1)  !!! this has changed
      !wp%K=max(2048,wp%K)

      allocate(wp%Dd(ndir))

      wp%Dd=sum(wp%S_array, DIM = 1)

   end subroutine swanreader

   ! --------------------------------------------------------------
   ! -----------------Read variance density files -----------------
   ! --------------------------------------------------------------
   subroutine vardensreader(par,s,wp,fname)

      use params
      use spaceparams
      use xmpi_module
      use logging_module
      use readkey_module

      IMPLICIT NONE

      ! Input / output variables
      type(parameters), INTENT(INout)         :: par
      type(spacepars), intent(IN)             :: s
      type(waveparameters), INTENT(INOUT)     :: wp
      character(len=*), INTENT(IN)            :: fname

      ! Internal variables
      real*8                                  :: m0
      integer                                 :: nfreq, ndir,i
      integer                                 :: firstp,lastp
      real*8, dimension(:),allocatable        :: temp, findline

      wp%Npy=s%ny+1

      if(xmaster) then
         call writelog('ls','','Reading from VarDens file ',fname,' ...')
         open(44,file=fname,form='formatted',status='old')

         read(44,*)nfreq
      endif
#ifdef USEMPI
      !Dano call xmpi_bcast(nfreq)
#endif
      allocate(wp%f(nfreq))

      if(xmaster) then
         do i=1,nfreq
            read(44,*)wp%f(i)
         end do

         read(44,*)ndir
      endif
#ifdef USEMPI
      !Dano call xmpi_bcast(wp%f)
      !Dano call xmpi_bcast(ndir)
#endif
      allocate(wp%theta(ndir))

      if(xmaster) then
         do i=1,ndir
            read(44,*)wp%theta(i)
         end do
      endif
#ifdef USEMPI
      !Dano call xmpi_bcast(wp%theta)
#endif

      wp%theta=wp%theta*par%px/180
      wp%dang=wp%theta(2)-wp%theta(1)

      ! Read S_array
      allocate(wp%S_array(nfreq,ndir))

      if(xmaster) then
         do i=1,nfreq
            read(44,*)wp%S_array(i,:)
         end do

         close(44)                               ! Finished reading file
      endif

#ifdef USEMPI
      !Dano call xmpi_bcast(wp%S_array)
#endif

      ! Convert to m2/Hz/rad

      wp%S_array=wp%S_array*180/par%px

      allocate(wp%Sf(nfreq))
      wp%Sf = sum(wp%S_array, DIM = 2)*wp%dang

      allocate (temp(nfreq+1))
      temp(1)=0
      temp(2:nfreq)=0.5*wp%f(1:nfreq-1)+0.5*wp%f(2:nfreq)
      temp(nfreq+1)=wp%f(nfreq)
      ! Calculate zero-order moment
      m0=0.0d0
      !m1=0.0d0
      do i=1,nfreq
         m0=m0+wp%Sf(i)*(temp(i+1)-temp(i))
         !    m1=m1+wp%f(i)*wp%Sf(i)*(temp(i+1)-temp(i))
      end do
      deallocate (temp)

      ! Find main wave direction
      allocate (temp(ndir))
      temp=sum(wp%S_array, DIM = 1)
      i=maxval(maxloc(temp))
      wp%mainang=wp%theta(i)
      deallocate(temp)

      wp%hm0gew=4*sqrt(m0)
      !par%Trep=1/(m1/m0)
      call tpDcalc(wp%Sf,wp%f,par%Trep,par%trepfac,par%Tm01switch)
      call writelog('sl','','Trep = ',par%Trep,'s')
      call writelog('sl','','Hm0  = ',wp%hm0gew,'m')
      call writelog('sl','','Peak dir  = ',wp%mainang*180/par%px,' Cartesian degrees relative to x-axis')
      !
      allocate(findline(size(wp%Sf)))
      firstp=0
      lastp=0
      call frange(par, wp%Sf,firstp,lastp,findline)
      deallocate(findline)

      !!!!! ja/ap wp%K=max(100,ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1))

      wp%K=ceiling(wp%rt*(wp%f(lastp)-wp%f(firstp))+1)  !!! this has changed
      !wp%K=max(2048,wp%K)


      allocate(wp%Dd(ndir))

      wp%Dd=sum(wp%S_array, DIM = 1)

   end subroutine vardensreader

   ! --------------------------------------------------------------
   ! ---------------------- Build E_tdir file ---------------------
   ! --------------------------------------------------------------
   subroutine build_etdir(par,s,wp,Ebcfname)

      use params
      use math_tools
      use spaceparams
      use interp
      use xmpi_module
      use logging_module

      IMPLICIT NONE

      ! Input / output variables
      type(parameters), INTENT(IN)            :: par
      type(waveparameters), INTENT(INOUT)     :: wp
      type(spacepars), INTENT(IN)             :: s

      character(len=*), INTENT(IN)            :: Ebcfname

      ! Internal variables

      ! Help integers
      integer                                 :: Ns ! Number of theta bins
      integer                                 :: reclen
      logical                                 :: disptext

      ! Counters
      integer                                 :: i, ii, iii, stepf, stepang, index2
      integer                                 :: firstp, lastp, M

      ! Nothings
      integer                                 :: F2
      real*8                                  :: pp

      ! Help single variables with meaning
      real*8                                  :: TT, kmax
      real*8                                  :: hm0now, s1, s2, modf, modang
      real*8                                  :: stdeta,stdzeta

      ! Help vectors
      integer, dimension(wp%K)                :: Nbin
      real*8,dimension(size(wp%Sf))           :: findline
      real*8,dimension(size(wp%Dd))           :: Dmean, P
      real*8,dimension(wp%K)                  :: P0, k,phase, Sf0, A, Sf0org,S0org
      real*8,dimension(wp%K*2)                                :: randummy
      real*8,dimension(:),allocatable         :: temp, temp2, t, Nbox,rD
      real*8,dimension(1:401)                 :: ktemp, ftemp

      ! Help arrays
      real*8,dimension(:,:), allocatable      :: D
      real*8,dimension(:,:,:), allocatable    :: zeta, Ampzeta, E_tdir
      real*8,dimension(:,:), allocatable      :: eta, Amp

      ! Complex help variables
      ! double complex                          :: ctemp
      ! wwvv double complex,dimension(:),allocatable :: Gn, Comptemp
      complex(fftkind),dimension(:),allocatable   :: Gn, Comptemp

      ! Check whether maximum frequency is smaller than the Nyquist frequency,
      ! otherwise limit the time step to fit this frequency
      pp=maxval(wp%f)*2.d0
      if (wp%dt>(1.d0/pp)) then
         wp%dt=1.d0/pp
         if (xmaster) then
            call writelog('ls','(a)','Changing dtbc in wave boundary conditions to satisfy Nyquist condition:')
            call writelog('ls','(a,f0.4,a)','New dtbc = ',wp%dt,' s.')
         endif
      endif

      ! Print message to screen
      if(xmaster) then
         call writelog('ls','','Calculating wave energy at boundary')
      endif

      ! Determine frequencies around peak frequency of one-dimensional
      ! non-directional variance density spectrum, based on factor sprdthr, which
      ! should be included in the determination of the wave boundary conditions
      findline=0.0d0
      call frange(par,wp%Sf,firstp,lastp,findline)

      ! Determine number of frequencies in discrete variance density spectrum to be
      ! included (not equal to wp%K)
      M=int(sum(findline))

      ! Define one-dimensional non-directional variance density spectrum array       ! Bas: I guess this is unnecessary
      allocate (temp(size(wp%Sf)))
      temp=1.0d0

      ! Define two-dimensional spectrum array to be used to fill with directional
      ! spreading information in the next step
      allocate (D(size(findline),size(wp%Dd)))
      D=0.0d0

      ! Copy directional spreading array into two-dimensional spectrum array for each
      ! frequency
      do i=1,size(wp%Dd)
         do ii=1,size(findline)
            D(ii,i)=temp(ii)*wp%Dd(i)                                              ! Bas: temp(ii) is always unity here
         end do
      end do                                                                         ! Bas: D is filled with copies of the same row, which is wp%Dd

      deallocate (temp)

      ! Discard directional spreading information for frequencies that are outside
      ! the range around the peak frequency determined by sprdthr
      do i=1,size(wp%Dd)
         D(:,i)=D(:,i)*findline
      end do                                                                         ! Bas: D is still filled with copies of the same row, but about 90% of the values is zero now

      ! Determine the average contribution of a certain frequency for each direction
      Dmean=sum(D, DIM=1)/M                                                          ! Bas: D is still filled with copies of the same row, averaging in the first dimension thus results in that specific row again, still containing 90% zeros... so, Dmean=wp%Dd*findline

      ! Define number of wave components to be used
      allocate (temp(wp%K))
      temp=(/(i,i=0,wp%K-1)/)

      ! Select equidistant wave components between the earlier selected range of
      ! frequencies around the peak frequency based on sprdthr
      allocate(wp%fgen(wp%K))
      wp%fgen=temp*((wp%f(lastp)-wp%f(firstp))/(wp%K-1))+wp%f(firstp)
      deallocate(temp)

      ! Determine equidistant frequency step size
      wp%df=(wp%fgen(wp%K)-wp%fgen(1))/(dble(wp%K)-1.d0)
      ! Avoid leakage
      ! Pieter, Jaap and Ap 28/5: wp%fgen=floor(wp%fgen/wp%df)*wp%df taken out because this gives double values of fgen

      ! Determine equidistant period step size, which approximately equals wp%rt due
      ! to the dependence of wp%K on wp%rt
      TT=1/wp%df                                                                     ! Bas: because of ceiling statement in wp%K definition TT can deviate from wp%rt

      ! Determine maximum wave number based on maximum frequency and dispersion
      ! relation w^2 = g*k*tanh(k*d) and max(tanh(k*d))=1
      kmax=((2*(par%px)*wp%f(lastp))**2)/par%g

      ! Determine frequency array with 400 frequencies corresponding to wave numbers
      ! running from 0 to 2*kmax using the dispersion relation w^2 = g*k*tanh(k*d)
      allocate(temp(401))
      temp=(/(i,i=0,400)/)
      ktemp=temp*(kmax/200)
      ftemp=sqrt((par%g)*ktemp*tanh(ktemp*wp%h0t0))/(2*par%px)
      deallocate (temp)

      ! Determine all wave numbers of the selected frequencies around the peak
      ! frequency by linear interpolation
      do i=1,size(wp%fgen)
         call LINEAR_INTERP(ftemp,ktemp,401,wp%fgen(i),pp,F2)
         k(i)=pp
      end do

      ! Define normalization factor for wave variance
      pp=1/(sum(Dmean)*wp%dang)

      ! Determine normalized wave variance for each directional bin to be used as
      ! probability density function, so surface is equal to unity
      do i=1,size(wp%theta)
         ! P(i)=(sum(Dmean(1:i))-Dmean(i)/2)*wp%dang*pp                                            ! Bas: this is equal to P(i)=sum(Dmean(1:i))/sum(Dmean)
         P(i)=sum(Dmean(1:i))*wp%dang*pp
      end do

      ! Update random seed, if requested
      if (par%random==1) CALL init_seed

      ! Define random number between 0.025 and 0975 for each wave component
      !call random_number(randummy)
      ! wwvvrandom
      do i=1,wp%K*2
         randummy(i) = random(0)
      enddo

      P0=randummy(1:wp%K)
      !P0=0.95*P0+0.05/2                                                             ! Bas: do not crop cdf, only needed in Matlab to ensure monotonicity

      ! Define direction for each wave component based on random number and linear
      ! interpolation of the probability density function
      allocate(wp%theta0(wp%K))

      if (wp%scoeff >= 1000.d0) then
         wp%theta0 = wp%mainang
         !Ap longcrested waves
      else
         do i=1,size(P0)
            !call LINEAR_INTERP(P(1:size(P)-1),wp%theta(1:size(P)-1),size(P)-1,P0(i),pp,F2)
            call LINEAR_INTERP(P(1:size(P)),wp%theta(1:size(P)),size(P),P0(i),pp,F2)   ! Bas: do not crop cdf, only needed in Matlab to ensure monotonicity
            wp%theta0(i)=pp
         end do
      end if

      ! Determine number of time steps in wave record and make it even
      F2=nint(TT/wp%dt)                                                              ! Bas: why not simply use nint(wp%rt/wp%dt) ??
      if (mod(F2,2)/=0) then
         F2=F2+1
      end if

      ! Define time axis based on time steps
      allocate(t(F2))
      do i=1,F2
         t(i)=wp%dt*i
      end do

      ! Determine number of time steps in wave record and make it even               ! Bas: redundant with above, use wp%Nr = F2
      wp%Nr=nint(TT/wp%dt)
      if (mod(wp%Nr,2)/=0) then
         wp%Nr=wp%Nr+1
      end if

      ! Define a random phase for each wave component based on a subsequent series of
      ! random numbers
      phase=randummy(wp%K+1:2*wp%K)
      phase=2*phase*par%px

      ! Determine variance density spectrum values for all relevant wave components
      ! around the peak frequency by interpolation of two-dimensional variance
      ! density spectrum array. This is done by looping through the corresponding
      ! frequencies and find for each component the two-dimensional
      ! frequency/directional bin where it is located
      allocate(wp%S0(wp%K))
      do i=1,size(wp%fgen)

         ! Define frequency indices of frequencies that are equal or larger than the
         ! currently selected component around the peak frequency
         allocate(temp(size(wp%f)))
         allocate(temp2(size(wp%f)))
         temp2=(/(ii,ii=1,size(wp%f))/)
         temp=1
         where (wp%f < wp%fgen(i))
            temp=0
         endwhere
         temp=temp*temp2

         ! Check whether any indices are defined. If so, select the first selected
         ! index. Otherwise, select the last but one from all frequencies
         if (sum(temp)==0) then
            stepf=size(wp%f)-1
         else
            stepf=max(nint(minval(temp, MASK = temp .gt. 0)-1),1)
         end if

         ! Determine relative location of the currently selected component in the
         ! selected frequency bin
         modf=(wp%fgen(i)-wp%f(stepf))/(wp%f(stepf+1)-wp%f(stepf))
         deallocate(temp,temp2)

         ! Define directional indices of directions that are equal or larger than
         ! the currently selected component around the peak frequency
         allocate(temp(size(wp%theta)))
         allocate(temp2(size(wp%theta)))
         temp2=(/(ii,ii=1,size(wp%theta))/)
         temp=1
         where (wp%theta < wp%theta0(i) )
            temp=0
         endwhere
         temp=temp*temp2

         ! Check whether any indices are defined. If so, select the first selected
         ! index. Otherwise, select the first from all directions
         if (wp%theta0(i)==wp%theta(1)) then
            stepang=1
         else
            stepang=nint(minval(temp, MASK = temp .gt. 0) -1)
         end if

         ! Determine relative location of the currently selected component in the
         ! selected directional bin
         modang=(wp%theta0(i)-wp%theta(stepang))/(wp%theta(stepang+1)-wp%theta(stepang))
         deallocate(temp,temp2)

         ! Determine variance density spectrum values at frequency boundaries of
         ! selected two-dimensional bin by linear interpolation in the directional
         ! dimension along these two boundaries
         s1=(1.d0-modang)*wp%S_array(stepf,stepang)+modang*wp%S_array(stepf,stepang+1)
         s2=(1.d0-modang)*wp%S_array(stepf+1,stepang)+modang*wp%S_array(stepf+1,stepang+1)

         ! Determine variance density spectrum value at currently selected component
         ! around peak frequency by linear interpolation of variance density
         ! spectrum values at frequency boundaries in frequency direction
         wp%S0(i)=max(tiny(0.d0),0.00000001d0,(1.d0-modf)*s1+modf*s2)               ! Robert: limit to positive values only in case no energy is drawn
      end do

      ! Determine the variance density spectrum values for all relevant wave
      ! components around the peak frequency again, but now using the one-dimensional
      ! non-directional variance density spectrum only
      do i=1,size(wp%fgen)
         call Linear_interp(wp%f,wp%Sf,size(wp%f),wp%fgen(i),pp,F2)
         Sf0(i)=pp
      end do

      ! Determine significant wave height using Hm0 = 4*sqrt(m0) using the
      ! one-dimensional non-directional variance density spectrum
      hm0now = 4*sqrt(sum(Sf0)*wp%df)

      ! Backup original spectra just calculated
      Sf0org=Sf0
      S0org=wp%S0

      ! Correct spectra for wave height
      if (par%correctHm0 == 1) then
         wp%S0 = (wp%hm0gew/hm0now)**2*wp%S0                                            ! Robert: back on ?
         Sf0 = (wp%hm0gew/hm0now)**2*Sf0                                                ! Robert: back on ?
      endif

      ! Determine directional step size
      allocate(wp%dthetafin(wp%K))
      wp%dthetafin = Sf0/wp%S0

      ! Determine amplitude of components from two-dimensional variance density
      ! spectrum using Var = 1/2*a^2
      A = sqrt(2*wp%S0*wp%df*wp%dthetafin)

      ! Restore original spectra just calculated
      Sf0=Sf0org
      wp%S0=S0org

      ! Allocate Fourier coefficients for each y-position at the seaside border and
      ! each time step
      allocate(wp%CompFn(wp%Npy,wp%Nr))
      wp%CompFn=0.d0

      ! Determine indices of relevant wave components around peak frequencies in
      ! Fourier transform result
      allocate(wp%index_vector(wp%K))
      wp%index_vector = floor(wp%f(firstp)/wp%df)+1+nint((wp%fgen-wp%f(firstp))/wp%df)    ! Bas: too complex ?? 1+nint((wp%fgen-wp%f(firstp))/wp%df) equals 1:size(wp%fgen)

      ! Determine first half of complex Fourier coefficients of relevant wave
      ! components for first y-coordinate using random phase and amplitudes from
      ! sampled spectrum until Nyquist frequency. The amplitudes are represented in a
      ! two-sided spectrum, which results in the factor 1/2.
      do i=1,wp%K
         wp%CompFn(1,wp%index_vector(i)) = A(i)/2*exp(par%compi*phase(i))           ! Bas: wp%index_vector used in time dimension because t=j*dt in frequency space
      enddo

      ! Determine Fourier coefficients beyond Nyquist frequency (equal to
      ! coefficients at negative frequency axis) of relevant wave components for
      ! first y-coordinate by mirroring
      allocate(Comptemp(size(wp%CompFn(1,wp%Nr/2+2:wp%Nr))))                         ! Bas: too complex ?? size(wp%CompFn(1,wp%Nr/2+2:wp%Nr)) equals (1,size(wp%Nr)/2-2)
      Comptemp = conjg(wp%CompFn(1,2:wp%Nr/2))
      call flipiv(Comptemp,size(Comptemp))
      wp%CompFn(1,wp%Nr/2+2:wp%Nr)=Comptemp

      ! Determine Fourier coefficients for all other y-coordinates along seaside
      ! border in the same manner
      do index2=2,wp%Npy
         wp%CompFn(index2,wp%index_vector)=wp%CompFn(1,wp%index_vector)* &
         exp(-par%compi*k*(dsin(wp%theta0)*(s%yz(1,index2)-s%yz(1,1)) &
         +dcos(wp%theta0)*(s%xz(1,index2)-s%xz(1,1))) )

         Comptemp = conjg(wp%CompFn(index2,2:wp%Nr/2))
         call flipiv(Comptemp,size(Comptemp))
         wp%CompFn(index2,wp%Nr/2+2:wp%Nr)=Comptemp
      end do

      deallocate(Comptemp)

      ! Determine directional bins in computation (not spectrum) ensuring that
      ! s%thetamax is included in the range
      Ns=s%ntheta
      !allocate(temp(Ns+1))
      !temp=(/(i,i=0,Ns)/)
      !temp(Ns+1)=temp(Ns+1)+epsilon(1.d0)
      allocate(Nbox(Ns+1))
      allocate(rD(Ns))
      do i=1,Ns+1
         Nbox(i)=s%thetamin+(i-1)*s%dtheta
      enddo

      if (par%wbctype==WBCTYPE_PARAMETRIC .or. par%wbctype==WBCTYPE_JONS_TABLE) then
         rD = dcos(0.5d0*(Nbox(1:Ns)+Nbox(2:Ns+1))-wp%mainang)**(2*nint(wp%scoeff))
         rD = rD/sum(rD)
      endif


      !deallocate (temp)

      ! try to put all wave directions between thetamax and thetamin
      do i=1,size(wp%theta0)
         if (wp%theta0(i)>s%thetamax) then
            F2=1
         elseif (wp%theta0(i)<s%thetamin) then
            F2=-1
         else
            F2=0
         endif
         ! now turn over 2pi
         if (F2==1) then
            do while (F2==1)
               wp%theta0(i)=wp%theta0(i)-2*par%px
               if (wp%theta0(i)<=s%thetamax) then
                  F2=0
               endif
            enddo
         elseif (F2==-1) then
            do while (F2==-1)
               wp%theta0(i)=wp%theta0(i)+2*par%px
               if (wp%theta0(i)>=s%thetamin) then
                  F2=0
               endif
            enddo
         endif
      enddo

      ! Determine computational directional bin for each wave component
      do i=1,size(wp%theta0)
         do ii=1,Ns
            if (wp%theta0(i)>=Nbox(ii) .and. wp%theta0(i)<Nbox(ii+1)) then
               Nbin(i)=ii
            endif
         enddo
      enddo

      !    Nbin(i)=ceiling((wp%theta0(i)-Nbox(1))/(par%dtheta*par%px/180.d0))
      !
      !    ! Ensure lower bin boundaries to be part of succeeding bin
      !    if (mod((wp%theta0(i)-Nbox(1)),(par%dtheta*par%px/180.d0))==0) then
      !        Nbin(i)=Nbin(i)+1
      !    end if
      !end do

      ! Determine highest bin containing energy
      i=(maxval(Nbin))                                                               ! Bas: not used

      ! Add wave components outside computational directional bins to outer bins, if
      ! nspr parameter is set to one
      if (par%nspr==1) then
         do i=1,wp%K
            if (Nbin(i)<=0) then
               Nbin(i)=1
            elseif (Nbin(i)>Ns) then
               Nbin(i)=Ns
            endif
            wp%theta0(i)=s%theta(Nbin(i))
         enddo
      endif

      ! deallocate(Nbox)

      ! Define time window that gradually increases and decreases energy input over
      ! the wave time record according to tanh(fc*t/max(t))^2*tanh(fc*(1-t/max(t)))^2
      allocate(wp%window(size(t)))
      allocate(temp(size(t)))
      temp=t
      where (t>wp%rt)
         temp=0                                                         ! Bas: skip time steps that exceed wp%rt, might not be necessary when wp%Nr, F2 and t are defined well, see above
      endwhere
      wp%window=1
      wp%window=wp%window*(tanh(192.d0*temp/maxval(temp))**2)*(tanh(192.d0*(1.d0-temp/maxval(temp)))**2)
      !wp%window=wp%window*(tanh(192.d0*t/maxval(t))**2)*(tanh(192.d0*(1.d0-t/maxval(t)))**2)      ! Bas: array t matches wp%rt, so truncating via temp is not necessary anymore
      deallocate(temp)

      ! Allocate variables for water level exitation and amplitude with and without
      ! directional spreading dependent envelope
      allocate(zeta(wp%Npy,wp%Nr,Ns))
      allocate(Ampzeta(wp%Npy,wp%Nr,Ns))
      zeta=0
      Ampzeta=0

      allocate(eta(wp%Npy,wp%Nr))
      allocate(Amp(wp%Npy,wp%Nr))
      eta=0
      Amp=0

      ! Calculate wave energy for each computational directional bin
      do ii=1,Ns

         ! Print message to screen
         if(xmaster) then
            call writelog('ls','(A,I0,A,I0)','Calculating wave energy for theta bin ',ii,' of ',Ns)
         endif

         ! Calculate wave energy for each y-coordinate along seaside boundary for
         ! current computational directional bin
         disptext = .true.
         do index2=1,wp%Npy

            ! Select wave components that are in the current computational
            ! directional bin
            allocate(Gn(wp%Nr))
            Gn=0
            allocate(temp(size(Nbin)))
            temp=0
            where (Nbin==ii)
               temp=1.
            end where

            ! Determine number of wave components that are in the current
            ! computational directional bin
            F2=nint(sum(temp))

            ! Check whether any wave components are in the current computational
            ! directional bin
            if (F2/=0) then

               ! Determine for each wave component in the current computational
               ! directional bin it's index in the Fourier coefficients array
               ! ordered from hight to low frequency
               allocate(temp2(F2))
               temp2=0
               do i=1,F2
                  iii=maxval(maxloc(temp))
                  temp(iii)=0
                  temp2(i)=wp%index_vector(iii)
               end do

               ! Determine Fourier coefficients of all wave components for current
               ! y-coordinate in the current computational directional bin
               Gn(int(temp2))=wp%CompFn(index2,int(temp2))
               deallocate(temp2)
               allocate(Comptemp(size(Gn(wp%Nr/2+2:wp%Nr))))
               Comptemp = conjg(Gn(2:wp%Nr/2))
               call flipiv(Comptemp,size(Comptemp))
               Gn(wp%Nr/2+2:wp%Nr)=Comptemp
               deallocate(Comptemp)

               ! Inverse Discrete Fourier transformation to transform back to time
               ! domain from frequency domain
               allocate(Comptemp(size(Gn)))
               Comptemp=Gn
               F2=0
               Comptemp=fft(Comptemp,inv=.true.,stat=F2)

               ! Scale result
               Comptemp=Comptemp/sqrt(real(size(Comptemp)))

               ! Superimpose gradual increase and decrease of energy input for
               ! current y-coordinate and computational diretional bin on
               ! instantaneous water level excitation
               zeta(index2,:,ii)=dble(Comptemp*wp%Nr)*wp%window
               Comptemp(:)=zeta(index2,:,ii)
               !
               ! Hilbert tranformation to determine envelope for each directional bin seperately
               call hilbert(Comptemp,size(Comptemp))

               ! Integrate instantaneous water level excitation of wave
               ! components over directions
               eta(index2,:) = sum(zeta(index2,:,:),2)
               Comptemp=eta(index2,:)
               !
               ! Hilbert transformation to determine envelope of all total non-directional wave components
               call hilbert(Comptemp,size(Comptemp))
               !
               ! Determine amplitude of water level envelope by calculating
               ! the absolute value of the complex wave envelope descriptions
               Amp(index2,:)=abs(Comptemp)
               !
               ! Calculate standard deviations of directional and
               ! non-directional instantaneous water level excitation of all
               ! wave components to be used as weighing factor
               stdzeta = sqrt(sum(zeta(index2,:,ii)**2)/(size(zeta(index2,:,ii)-1)))
               stdeta = sqrt(sum(eta(index2,:)**2)/(size(eta(index2,:)-1)))
               !
               ! Calculate amplitude of directional wave envelope
               Ampzeta(index2,:,ii)= Amp(index2,:)*stdzeta/stdeta
               !
               ! Print status message to screen
               !
               if(xmaster) then
                  if (F2/=0) then
                     call writelog('ls','(A,I0,A,I0,A,I0)','Y-point ',index2,' of ',wp%Npy,' done. Error code: ',F2)
                  else
                     call writelog('s','(A,I0,A,I0,A)','Y-point ',index2,' of ',wp%Npy,' done.')
                  end if
               endif

               deallocate(Comptemp)
            else
               ! Current computational directional bin does not contain any wave
               ! components, so print message to screen
               if(xmaster) then
                  if (disptext) then
                     call writelog('ls','(A,I0,A)','Theta bin ',ii,' empty at this point. Continuing to next point')
                     disptext = .false.
                  endif
               endif
            end if

            deallocate(temp)
            deallocate(Gn)
         end do
      end do



      ! Print message to screen
      if(xmaster) then
         call writelog('ls','','writing wave energy to ',trim(Ebcfname),' ...')
      endif

      ! Determine energy distribution over y-coordinates, computational directional
      ! bins and time using E = 1/2*rho*g*a^2
      allocate(E_tdir(wp%Npy,wp%Nr,Ns))

      ! Jaap: apply symmetric distribution in case of jons or jons_table
      ! REMARK: in this case printing messages to screen can be erroneous
      if (par%wbctype==WBCTYPE_PARAMS .or. par%wbctype==WBCTYPE_JONS_TABLE) then
         call writelog('ls','','Apply symmetric energy distribution w.r.t mean wave direction')
         do ii=1,Ns
            do index2=1,wp%Npy
               E_tdir(index2,:,ii)=0.0d0
               E_tdir(index2,:,ii)=0.5d0*(par%rho)*(par%g)*Amp(index2,:)**2*rD(ii)
            enddo
         enddo
      else
         E_tdir=0.0d0
         E_tdir=0.5d0*(par%rho)*(par%g)*Ampzeta**2
      endif
      E_tdir=E_tdir/s%dtheta


      deallocate(Nbox)
      deallocate(rD)

      ! Write energy distribution to BCF file
      if(xmaster) then
         inquire(iolength=reclen) 1.d0
         reclen=reclen*(wp%Npy)*(Ns)
         write(*,*) 'Opening', trim(Ebcfname)
         open(12,file=trim(Ebcfname),form='unformatted',access='direct',recl=reclen)
         do i=1,wp%Nr+4                                                             ! Bas: why add 4 ??
            write(12,rec=i)E_tdir(:,min(i,wp%Nr),:)
         end do
         close(12)
         call writelog('sl','','file done')
      endif

      deallocate (D,t,zeta,Ampzeta,E_tdir, Amp, eta)

      return

   end subroutine build_etdir

   ! --------------------------------------------------------------
   ! ----------------------- Bound long wave ----------------------
   ! --------------------------------------------------------------
   subroutine build_boundw(par,s,wp,qbcfname)

      use params
      use spaceparams
      use math_tools
      use xmpi_module
      use logging_module

      IMPLICIT NONE


      ! Input / output variables

      type(parameters), INTENT(IN)            :: par
      type(spacepars), INTENT(IN)             :: s
      type(waveparameters), INTENT(INOUT)     :: wp
      character(len=*), INTENT(IN)            :: qbcfname

      ! Internal variables
      integer                                 :: K, m, index1, Npy, Nr, i=0, jj
      integer                                 :: reclen
      integer,dimension(:),allocatable        :: index2

      logical                                 :: firsttime

      real*8                                  :: g
      real*8                                  :: df, deltaf
      real*8,dimension(:), allocatable        :: w1, k1
      real*8,dimension(:), allocatable        :: Ebnd
      real*8,dimension(:), allocatable        :: term1, term2, term2new, dif, chk1, chk2
      real*8,dimension(:,:),allocatable       :: Eforc, D, deltheta, KKx, KKy, theta3
      real*8,dimension(:,:),allocatable       :: dphi3, k3, cg3, Abnd
      real*8,dimension(:,:,:),allocatable     :: q

      double complex,dimension(:),allocatable       :: Comptemp, Comptemp2
      complex(fftkind),dimension(:,:),  allocatable :: Gn, Ftemp2
      complex(fftkind),dimension(:,:,:),allocatable :: Ftemp
      character(4),dimension(4)               :: qstr

      g=par%g
      K=wp%K
      df=wp%df
      index1=wp%index_vector(1)
      Npy=wp%Npy
      Nr=wp%Nr

      ! Print message to screen
      if(xmaster) then
         call writelog('sl','', 'Calculating flux at boundary')
      endif

      ! Allocate two-dimensional variables for all combinations of interacting wave
      ! components to be filled triangular
      allocate(Eforc(K-1,K),D(K-1,K),deltheta(K-1,K),KKx(K-1,K),KKy(K-1,K))
      allocate(dphi3(K-1,K),k3(K-1,K),cg3(K-1,K))

      ! Initialize variables as zero
      Eforc = 0
      D = 0
      deltheta = 0
      KKx = 0
      KKy = 0
      dphi3 = 0
      k3 = 0
      cg3 = 0

      ! Allocate variables for angular velocity and wave numbers for wave components
      allocate(w1(size(wp%fgen)),k1(size(wp%fgen)))

      ! Initialize variables as zero
      w1=0
      k1=0

      ! Determine for each wave component interactions with all other wave components
      ! as far as not processed yet (each loop step the number of interactions thus
      ! decrease with one)

      ! First time is set true for each time new wave bc are generated
      firsttime = .true.

      do m=1,K-1

         ! Determine difference frequency
         deltaf=m*df

         ! Determine angular velocity of primary waves
         w1=2*par%px*wp%fgen

         ! Determine wave numbers of primary waves
         call bc_disper(k1,w1,size(w1),wp%h0t0,g)

         ! Determine difference angles (pi already added)
         deltheta(m,1:K-m) = abs(wp%theta0(m+1:K)-wp%theta0(1:K-m))+par%px

         ! Determine x- and y-components of wave numbers of difference waves
         KKy(m,1:K-m)=k1(m+1:K)*dsin(wp%theta0(m+1:K))-k1(1:K-m)*dsin(wp%theta0(1:K-m))
         KKx(m,1:K-m)=k1(m+1:K)*dcos(wp%theta0(m+1:K))-k1(1:K-m)*dcos(wp%theta0(1:K-m))

         ! Determine difference wave numbers according to Van Dongeren et al. 2003
         ! eq. 19
         k3(m,1:K-m) =sqrt(k1(1:K-m)**2+k1(m+1:K)**2+2*k1(1:K-m)*k1(m+1:K)*dcos(deltheta(m,1:K-m)))

         ! Determine group velocity of difference waves
         cg3(m,1:K-m)= 2.d0*par%px*deltaf/k3(m,1:K-m)

         ! Ideetje Robert Jaap: make sure that we don't blow up bound long wave
         !                      when offshore boundary is too close to shore
         ! cg3 = min(cg3,par%nmax*sqrt(par%g*wp%h0t0))
         cg3(m,1:K-m) = min(cg3(m,1:K-m),par%nmax*sqrt(g/k3(m,1:K-m)*tanh(k3(m,1:K-m)*wp%h0t0)))

         ! Determine difference-interaction coefficient according to Herbers 1994
         ! eq. A5
         allocate(term1(K-m),term2(K-m),term2new(K-m),dif(K-m),chk1(K-m),chk2(K-m))

         term1 = (-w1(1:K-m))*w1(m+1:K)
         term2 = (-w1(1:K-m))+w1(m+1:K)
         term2new = cg3(m,1:K-m)*k3(m,1:K-m)
         dif = (abs(term2-term2new))
         if (any(dif>0.01*term2 .and. firsttime .eqv. .true.)) then
            firsttime = .false.
            call writelog('lws','','Warning: shallow water so long wave variance is reduced using par%nmax');
         endif

         chk1  = cosh(k1(1:K-m)*wp%h0t0)
         chk2  = cosh(k1(m+1:K)*wp%h0t0)

         D(m,1:K-m) = -g*k1(1:K-m)*k1(m+1:K)*dcos(deltheta(m,1:K-m))/2.d0/term1+g*term2*(chk1*chk2)/ &
         ((g*k3(m,1:K-m)*tanh(k3(m,1:K-m)*wp%h0t0)-(term2new)**2)*term1*cosh(k3(m,1:K-m)*wp%h0t0))* &
         (term2*((term1)**2/g/g - k1(1:K-m)*k1(m+1:K)*dcos(deltheta(m,1:K-m))) &
         - 0.50d0*((-w1(1:K-m))*k1(m+1:K)**2/(chk2**2)+w1(m+1:K)*k1(1:K-m)**2/(chk1**2)))

         deallocate(term1,term2,term2new,dif,chk1,chk2)

         ! Correct for surface elevation input and output instead of bottom pressure
         ! so it is consistent with Van Dongeren et al 2003 eq. 18
         D(m,1:K-m) = D(m,1:K-m)*cosh(k3(m,1:K-m)*wp%h0t0)/(cosh(k1(1:K-m)*wp%h0t0)*cosh(k1(m+1:K)*wp%h0t0))

         ! Exclude interactions with components smaller than or equal to current
         ! component according to lower limit Herbers 1994 eq. 1
         where(wp%fgen<=m*df)
            D(m,:)=0.d0                                           ! Bas: redundant with initial determination of D ??
         endwhere

         ! Exclude interactions with components that are cut-off by the fcutoff
         ! parameter
         if (m*df<=par%fcutoff) D(m,:)=0.d0

         ! Determine energy of bound long wave according to Herbers 1994 eq. 1 based
         ! on difference-interaction coefficient and energy density spectra of
         ! primary waves
         Eforc(m,1:K-m) = 2*D(m,1:K-m)**2*wp%S0(1:K-m)*wp%S0(m+1:K)*wp%dthetafin(1:K-m)*wp%dthetafin(m+1:K)*df

         ! Determine phase of bound long wave assuming a local equilibrium with
         ! forcing of interacting primary waves according to Van Dongeren et al.
         ! 2003 eq. 21 (the angle is the imaginary part of the natural log of a
         ! complex number as long as the complex number is not zero)
         allocate(Comptemp(K-m),Comptemp2(K-m))
         Comptemp=conjg(wp%CompFn(1,index1+m:index1+K-1))
         Comptemp2=conjg(wp%CompFn(1,index1:index1+K-m-1))
         dphi3(m,1:K-m) = par%px+imag(log(Comptemp))-imag(log(Comptemp2))
         deallocate (Comptemp,Comptemp2)

      end do

      ! Determine angle of bound long wave according to Van Dongeren et al. 2003 eq. 22
      allocate(theta3(K-1,K))
      where (abs(KKx)>0.00001d0)
         theta3 = atan(KKy/KKx)
      elsewhere
         theta3 = atan(KKy/sign(0.00001d0,KKx))
      endwhere

      ! Allocate variables for amplitude and Fourier coefficients of bound long wave
      allocate(Gn(Npy,Nr))
      allocate(Abnd(K-1,K))
      allocate(Ebnd(K-1))

      Ebnd = sum(Eforc,2)

      Abnd = sqrt(2*Eforc*df)

      allocate(index2(K-1))
      index2=(/(i,i=1,K-1)/)

      ! Determine complex description of bound long wave per interaction pair of
      ! primary waves for first y-coordinate along seaside boundary
      allocate(Ftemp(K-1,K,4)) ! Jaap qx, qy qtot
      Ftemp(:,:,1) = Abnd/2*exp(-1*par%compi*dphi3)*cg3*dcos(theta3) ! qx
      Ftemp(:,:,2) = Abnd/2*exp(-1*par%compi*dphi3)*cg3*dsin(theta3) ! qy
      Ftemp(:,:,3) = Abnd/2*exp(-1*par%compi*dphi3)*cg3              ! qtot
      Ftemp(:,:,4) = Abnd/2*exp(-1*par%compi*dphi3)                  ! eta

      allocate(q(Npy,Nr,4))           ! qx qy qtot eta
      allocate(Comptemp(Nr/2-1))
      allocate(Comptemp2(Nr))         ! Allocate mass flux as function of time
      allocate(Ftemp2(K-1,K))

      q=0.0d0

      qstr = (/'qx  ','qy  ','qtot','eta '/)

      do m=1,4
         ! Determine complex description of bound long wave per primary wave component
         ! for first y-coordinate along seaside boundary
         Gn=0
         Gn(1,index2+1)=(sum(Ftemp(:,:,m),DIM=2))

         ! Determine Fourier coefficients

         Comptemp=conjg(Gn(1,2:Nr/2))
         call flipiv(Comptemp,Nr/2-1)
         Gn(1,Nr/2+2:Nr)=Comptemp

         ! Fourier transformation
         Comptemp2=fft(Gn(1,:),inv=.true.)

         ! Determine mass flux as function of time and let the flux gradually increase
         ! and decrease in and out the wave time record using the earlier specified
         ! window
         Comptemp2=Comptemp2/sqrt(real(Nr))
         q(1,:,m)=real(Comptemp2*Nr)*wp%window

         ! Determine mass flux of bound long wave for every y-coordinate at the seaside
         ! boundary

         do jj=2,Npy

            ! Determine phase shift due to y-coordinate and y-component wave number
            ! with respect to first y-coordinate alogn seaside boundary
            Ftemp2 = Ftemp(:,:,m)*exp(-1*par%compi*(KKy*(s%yz(1,jj)-s%yz(1,1))+KKx*(s%xz(1,jj)-s%xz(1,1))))

            ! Determine Fourier coefficients
            Gn(jj,index2+1) = (sum(Ftemp2,DIM=2))
            Comptemp = conjg(Gn(jj,2:Nr/2))
            call flipiv(Comptemp,Nr/2-1)
            Gn(jj,Nr/2+2:Nr) = Comptemp

            ! Inverse Discrete Fourier transformation to transform back to time space
            ! from frequency space
            Comptemp2=fft(Gn(jj,:),inv=.true.)

            ! Determine mass flux as function of time and let the flux gradually
            ! increase and decrease in and out the wave time record using the earlier
            ! specified window
            Comptemp2=Comptemp2/sqrt(real(Nr))    ! u-direction
            q(jj,:,m)=real(Comptemp2*Nr)*wp%window

         end do !jj loop
      end do !m loop

      ! Jaap and Bas: Fix for curvi-grids
      ! REMARK: Need to do something about Ftemp?
      do jj=2,Npy
         q(jj,:,1) = dcos(datan(q(jj,:,2)/max(q(jj,:,1),par%eps))-s%alfaz(1,jj))*q(jj,:,3)
         q(jj,:,2) = dsin(datan(q(jj,:,2)/max(q(jj,:,1),par%eps))-s%alfaz(1,jj))*q(jj,:,3)
      end do

      deallocate(Comptemp)
      deallocate(Comptemp2)
      deallocate(Ftemp)
      deallocate(Ftemp2)
      deallocate(index2)

      ! Write bound long wave flux to BCF file
      if(xmaster) then
         call writelog('ls','','writing long wave mass flux to ',trim(qbcfname),' ...')
         inquire(iolength=reclen) 1.d0
         reclen=reclen*(Npy*4)
         open(21,file=qbcfname,form='unformatted',access='direct',recl=reclen)
         do i=1,wp%Nr+4                                                             ! Bas: why add 4 ??
            write(21,rec=i)q(:,min(i,wp%Nr),:)
         end do
         close(21)
         call writelog('sl','','file done')
      endif

      deallocate(wp%index_vector,wp%S0,wp%dthetafin,wp%fgen,wp%theta0,wp%window)
      deallocate(wp%Sf,wp%Dd,wp%f,wp%theta,wp%S_array,wp%CompFn)

   end subroutine build_boundw




   ! -----------------------------------------------------------
   ! --------- JONSWAP  unscaled JONSWAP spectrum --------------
   ! ----------------(used by build_jonswap)--------------------
   subroutine jonswapgk(x,gam,y)

      IMPLICIT NONE
      ! Required input: - x           : nondimensional frequency, divided by the peak frequency
      !                 - gam         : peak enhancement factor, optional parameter (DEFAULT 3.3)
      !                 - y is output : nondimensional relative spectral density, equal to one at the peak

      real*8, INTENT(IN)                  :: gam
      real*8,dimension(:), INTENT(IN)     :: x
      real*8,dimension(:), INTENT(INOUT)  :: y

      ! Internal variables
      real*8,dimension(size(x))           :: xa, sigma, fac1, fac2, fac3, temp

      xa=abs(x)

      where (xa==0)
         xa=1e-20
      end where

      sigma=xa

      where (sigma<1.)
         sigma=0.07
      end where

      where (sigma>=1.)
         sigma=0.09
      end where

      temp=0*xa+1

      fac1=xa**(-5)
      fac2=exp(-1.25*(xa**(-4)))
      fac3=(gam*temp)**(exp(-((xa-1)**2)/(2*(sigma**2))))

      y=fac1*fac2*fac3
      y=y/maxval(y)

      return

   end subroutine jonswapgk

   ! -----------------------------------------------------------
   ! ---- Small subroutine to determine f-range round peak -----
   ! ----(used by build_jonswap, swanreader, vardensreader)-----
   subroutine frange(par,Sf,firstp,lastp,findlineout)

      use params

      implicit none

      type(parameters)                        :: par


      real*8, dimension(:), intent(in)        :: Sf
      integer, intent(out)                    :: firstp, lastp

      real*8, dimension(:), intent(out)       :: findlineout
      real*8, dimension(:),allocatable        :: temp, findline
      integer                                 :: i = 0

      allocate(findline(size(Sf)))
      findline=0*Sf                           ! find frequency range around peak

      where (Sf>par%sprdthr*maxval(Sf))
         findline=1
      end where


      firstp=maxval(maxloc(findline))         ! Picks the first "1" in temp

      allocate (temp(size(findline)))
      temp=(/(i,i=1,size(findline))/)
      lastp=maxval(maxloc(temp*findline))     ! Picks the last "1" in temp

      findlineout=findline
      deallocate(temp, findline)

   end subroutine frange


   ! -----------------------------------------------------------
   ! ----------- Small subroutine to determine tpD -------------
   ! ----(used by build_jonswap, swanreader, vardensreader)-----
   subroutine tpDcalc(Sf,f,Trep,trepfac,switch)

      implicit none

      real*8, dimension(:), intent(in)        :: Sf, f
      real*8, intent(out)                     :: Trep
      real*8, intent(in)                      :: trepfac
      integer, intent(in)                     :: switch

      real*8, dimension(:),allocatable        :: temp

      allocate(temp(size(Sf)))
      temp=0.d0
      where (Sf>=trepfac*maxval(Sf))
         temp=1.d0
      end where

      if (switch == 1) then
         Trep=sum(temp*Sf)/sum(temp*Sf*f)    ! Tm01
      else
         Trep = sum(temp*Sf/f)/sum(temp*Sf)    ! Tm-1,0
      endif

   end subroutine tpDcalc


   ! --------------------------------------------------------------
   ! --------------------- Dispersion relation --------------------
   ! ----------------- (used only by build_boundw) ----------------
   subroutine bc_disper(k1,w1,m,h,g)
      !          k  = wave number             (2 * pi / wave length)
      !          w  = wave angular frequency  (2 * pi / wave period)
      !          m  = size k and w vectors
      !          h  = water depth
      !          g  = gravitational acceleration constant, optional (DEFAULT 9.81)
      !
      !          absolute error in k*h < 5.0e-16 for all k*h
      !
      !
      !          original Matlab code by: G. Klopman, Delft Hydraulics, 6 Dec 1994

      integer, intent(in)                     :: m
      real*8,dimension(m),intent(in)          :: w1
      real*8,dimension(m),intent(out)         :: k1
      real*8, intent(in)                      :: h, g

      ! internal variables

      real*8,dimension(m)                     :: w2,q,thq,thq2,a,b,c,arg,sign
      integer                                 :: j
      real*8                                  :: hu

      w2 = w1**2*(h/g)
      q = w2/(1.0d0-exp(-(w2**(5.0d0/4.0d0))))**(2.0d0/5.0d0)

      do j=1,4
         thq  = tanh(q)
         thq2 = 1.0d0-thq**2
         a    = (1.0d0-q*thq)*thq2
         b    = thq + q*thq2
         c    = q*thq-w2
         where (abs(a*c)<(b**2*1.0e-8))
            arg = -c/b
         elsewhere
            arg  = (b**2)-4.0d0*a*c
            arg  = (-b + sqrt(arg))/(2.0d0*a)
         endwhere
         q    = q+arg
      end do

      where (w1>0.0d0)
         sign=1.0d0
      endwhere

      where (w1==0.0d0)
         sign=0.0d0
      endwhere

      where (w1<0.0d0)
         sign=-1.0d0
      endwhere

      k1 = sign*q/h

      where (k1==huge(hu))
         k1=0.0d0
      endwhere

      where (k1==-1.0d0*huge(hu))
         k1=0.0d0
      endwhere

      return

   end subroutine bc_disper



   subroutine init_seed   ! based on usage of random  wwvv
      use math_tools
      integer clock
      real*8 dummy
      call system_clock(count=clock)
      dummy = random(clock)
      return
   end subroutine init_seed

end module waveparams
