! Everything in this module should be handled by xmaster, so only call using xmaster
module spectral_wave_bc_module
  use typesandkinds
  implicit none

  type spectrum                                         ! These are related to input spectra
     real*8,dimension(:,:),pointer          :: S        ! 2D variance density spectrum
     real*8,dimension(:),pointer            :: f,ang    ! 1D frequency and direction vectors for S
     real*8,dimension(:),pointer            :: Sf       ! S integrated over all directions
     integer                                :: nf,nang  ! number of frequencies and angles
     real*8                                 :: df,dang  ! frequency and angle step size
     real*8                                 :: hm0,fp,dir0,scoeff  ! imposed significant wave height, peak frequency,
     ! main wave angle and spreading coefficient
     real*8                                 :: trep,dirm ! representative period and mean wave direction
  endtype spectrum
  type shortspectrum
     real*8,dimension(:),pointer            :: Sf         ! S integrated over all directions
     !   real*8,dimension(:),pointer            :: Srep       ! S at frequency/direction locations of wave train components
     !   real*8                                 :: dangrep    ! representative angle step size for Srep
  endtype shortspectrum
  type waveparamsnew                                      ! These are related to the generated bc series
     real*8                                 :: h0         ! average water depth on offshore boundary
     integer                                :: K          ! number of wave train components
     real*8                                 :: rtbc,dtbc  ! duration and time step to be written to boundary condition file
     real*8                                 :: rtin,dtin  ! duration and time step for the internal time axis, based on Fourier
     ! limitations and the number of wave train components (wp%K)
     real*8,dimension(:),pointer            :: tin        ! internal time axis
     real*8,dimension(:),pointer            :: taperf,taperw  ! internal taper function for flow and waves
     integer                                :: tslen      ! internal time axis length (is always even, not odd)
     integer                                :: tslenbc    ! time axis length for boundary condition file
     logical                                :: dtchanged  ! quick check to see if dtin == dtbc (useful for interpolation purpose)
     real*8,dimension(:),pointer            :: fgen,thetagen,phigen,kgen,wgen ! frequency, angle, phase, wave number, radian frequency
     ! of wavetrain components for boundary signal
     real*8                                 :: dfgen      ! frequency grid size in the generated components
     type(shortspectrum),dimension(:),pointer :: vargen   ! This is where the variance for each wave train at each spectrum location
     ! is stored
     type(shortspectrum),dimension(:),pointer :: vargenq  ! This is where the variance for each wave train at each spectrum location
     ! is stored, which is not scaled and is used for the generation of bound waves
     !   real*8,dimension(:),pointer            :: danggen    ! Representative integration angle for Srep to Sf, per offshore grid point
     real*8,dimension(:,:),pointer          :: A          ! Amplitude, per wave train component, per offshore grid point A(ny+1,K)
     real*8,dimension(:,:),pointer          :: Sfinterp   ! S integrated over all directions at frequency locations of fgen,
     ! per offshore grid point Sfinterp(ny+1,K)
     real*8,dimension(:,:),pointer          :: Sfinterpq  ! S integrated over all directions at frequency locations of fgen,
     ! per offshore grid point Sfinterpq(ny+1,K), uncorrected for generation of bound waves
     real*8,dimension(:),pointer            :: Hm0interp  ! Hm0 per offshore point, based on intergration of Sfinterp and used to scale
     ! final time series
     integer, dimension(:),pointer          :: Findex     ! Index of wave train component locations on frequency/Fourier axis
     integer, dimension(:),pointer          :: WDindex    ! Index of wave train component locations on wave directional bin axis
     double complex,dimension(:,:),pointer  :: CompFn     ! Fourier components of the wave trains
     character(slen)                        :: Efilename,qfilename,nhfilename
     real*8,dimension(:,:),pointer          :: zsits      ! time series of total surface elevation for nonhspectrum==1
     real*8,dimension(:,:),pointer          :: uits       ! time series of depth-averaged horizontal velocity nonhspectrum==1
     real*8,dimension(:,:),pointer          :: wits       ! time series of depth-averaged vertical velocity for nonhspectrum==1  ??
  endtype waveparamsnew
  type filenames                                      ! Place to store multiple file names
     character(slen)                        :: fname  ! file name of boundary condition file
     integer                                :: listline ! read position in FILELIST files
     logical                                :: reuse = .false. ! indicate to reuse this file every rtbc cycle
  endtype filenames

  ! These are for administration purposes and are initialized in initialize.F90
  integer,dimension(:),allocatable,save         :: n_index_loc     ! y-index locations of all input spectra, set in init spectrum
  integer,save                                  :: nspectra        ! number of input spectrs, set in init spectrum
  type(filenames),dimension(:),allocatable,save :: bcfiles         ! input wave spectrum files
  logical,save                                  :: reuseall        ! switch to reuse all of the wave boundary conditions
  integer,save                                  :: bccount         ! number of times boundary conditions have been generated, set in init spectrum
  real*8,save                                   :: spectrumendtime ! end time of boundary condition written to administration file
  real*8,dimension(:,:),allocatable,save        :: lastwaveelevation ! wave height at the end of the last spectrum
  integer,save                                  :: ind_end_taper   ! index of where the taper function equals rtbc
  ! These parameters control a lot how the spectra are handled. They could be put in params.txt,
  ! but most users will want to keep these at their default values anyway
  integer,parameter,private                 :: nfint = 801   ! size of standard 2D spectrum in frequency dimension
  integer,parameter,private                 :: naint = 401   ! size of standard 2D spectrum in angular dimension
  integer,parameter,private                 :: Kmin  = 200   ! minimum number of wave train components
  real*8,parameter,private                  :: wdmax = 5.d0  ! maximum depth*reliable angular wave frequency that can be resolved by
                                                             ! nonhydrostatic wave model. All frequencies above this are removed
                                                             ! from nonhspectrum generation

contains

  subroutine spectral_wave_bc(s,par,curline)
    use params
    use spaceparams
    use logging_module
    use filefunctions, only: create_new_fid
    use xmpi_module

    implicit none
    type(spacepars),intent(in)  :: s
    type(parameters),intent(inout) :: par
    integer,intent(out)            :: curline
    type(spectrum),dimension(:),allocatable         :: specin,specinterp
    type(spectrum)              :: combspec
    type(waveparamsnew)         :: wp             ! Most will be deallocated, but some variables useful to keep?
    integer                     :: iloc
    integer                     :: fidelist,fidqlist,fidnhlist ! file identifiers for ebcflist and qbcflist
    integer                     :: iostat
    real*8                      :: spectrumendtimeold,fmax

    ! Updat the number of boundary conditions generated by this module
    bccount = bccount+1

  if (xmaster) then
    if (reuseall) then
       ! Wave boundary conditions have already been computed
       ! Filename, rtbc, dtbc, etc. still kept in memory.
    else
       call writelog('l','','--------------------------------')
       call writelog('ls','','Calculating spectral wave boundary conditions ')
       call writelog('l','','--------------------------------')

       ! allocate temporary input storage
       if (.not. allocated(specin)) then
          allocate(specin(nspectra))
          allocate(specinterp(nspectra))
       endif

       ! Read through input spectra files
       fmax = 1.d0   ! assume 1Hz as maximum frequency. Increase in loop below if needed.
       do iloc = 1,nspectra
          call writelog('sl','(a,i0)','Reading spectrum at location ',iloc)
          ! Read input file 
          call read_spectrum_input(par,wp,bcfiles(iloc),specin(iloc))
          fmax = max(fmax,maxval(specin(iloc)%f))
       enddo
       do iloc = 1,nspectra
          call writelog('sl','(a,i0)','Interpreting spectrum at location ',iloc)
          ! Interpolate input 2D spectrum to standard 2D spectrum
          call interpolate_spectrum(specin(iloc),specinterp(iloc),par,fmax)
          call writelog('sl','','Values calculated from interpolated spectrum:')
          call writelog('sl','(a,f0.2,a)','Hm0       = ',specinterp(iloc)%hm0,' m')
          call writelog('sl','(a,f0.2,a)','Trep      = ',specinterp(iloc)%trep,' s')
          call writelog('sl','(a,f0.2,a)','Mean dir  = ',specinterp(iloc)%dirm,' degN')
       enddo

       ! Determine whether all the spectra are to be reused, which implies that the global reuseall should be
       ! set to true (no further computations required in future calls)
       call set_reuseall

       ! Determine the file names of the output boundary condition time series
       call set_bcfilenames(wp)

       ! calculate the mean combined spectra (used for combined Trep, determination of wave components, etc.)
       ! now still uses simple averaging, but could be improved to use weighting for distance etc.
       call generate_combined_spectrum(specinterp,combspec,par%px)
       par%Trep = combspec%trep
       call writelog('sl','(a,f0.2,a)','Overall Trep from all spectra calculated: ',par%Trep,' s')

       ! calculate the average offshore water depth, which is used in various
       ! computation in this module
       wp%h0 = calculate_average_water_depth(s)

       ! Wave trains that are used by XBeach. The number of wave trains, their frequencies and directions
       ! are based on the combined spectra of all the locations to ensure all wave conditions are
       ! represented in the XBeach model
       call generate_wavetrain_components(combspec,wp,par)

       ! We can now apply a correction to the wave train components if necessary. This section can be
       ! improved later
       if (nspectra==1 .and. specin(1)%scoeff>1000.d0) then
          ! this can be used both for Jonswap and vardens input
          wp%thetagen = specin(1)%dir0
       endif

       ! Set up time axis, including the time axis for output to boundary condition files and an
       ! internal time axis, which may differ in length to the output time axis
       call generate_wave_time_axis(par,wp)

       ! Determine the variance for each wave train component, at every spectrum location point
       call generate_wave_train_variance(wp,specinterp,par)

       ! Determine the amplitude of each wave train component, at every point along the
       ! offshore boundary
       call generate_wave_train_properties_per_offshore_point(wp,s)

       ! Generate Fourier components for all wave train component, at every point along
       ! the offshore boundary
       call generate_wave_train_Fourier(wp,s,par)

       ! time series of short wave energy or surface elevation
       if (par%nonhspectrum==0) then
          ! Distribute all wave train components among the wave direction bins. Also rearrage
          ! the randomly drawn wave directions to match the centres of the wave bins if the
          ! user-defined nspr is set on.
          call distribute_wave_train_directions(wp,s,par%px,par%nspr,.false.)

          ! Calculate the wave energy envelope per offshore grid point and write to output file
          call generate_ebcf(wp,s,par)
       else
          ! If user set nspr on, then force all short wave components to head along computational
          ! x-axis.
          call distribute_wave_train_directions(wp,s,par%px,par%nspr,.true.)    
          ! Generate time series of surface elevation, horizontal velocity and vertical velocity
          call generate_swts(wp,s,par)
       endif

       ! Calculate the bound long wave from the wave train components and write to output file
       call generate_qbcf(wp,s,par)

       ! Write non-hydrostatic time series of combined short and long waves if necessary
       if (par%nonhspectrum==1) then
          call generate_nhtimeseries_file(wp,par)
       endif

       ! Deallocate a lot of memory
       deallocate(specin,specinterp)
       deallocate(wp%tin,wp%taperf,wp%taperw)
       deallocate(wp%fgen,wp%thetagen,wp%phigen,wp%kgen,wp%wgen)
       deallocate(wp%vargen)
       deallocate(wp%vargenq)
       deallocate(wp%Sfinterp)
       deallocate(wp%Sfinterpq)
       deallocate(wp%Hm0interp)
       deallocate(wp%A)
       deallocate(wp%Findex)
       deallocate(wp%CompFn)
       if (par%nonhspectrum==0) then
          deallocate(wp%WDindex)
       else
          deallocate(wp%zsits)
          deallocate(wp%uits)
       endif
       ! Send message to screen and log
       call writelog('l','','--------------------------------')
       call writelog('ls','','Spectral wave boundary conditions complete ')
       call writelog('l','','--------------------------------')
    endif ! reuseall
    !
    ! The end time for the boundary conditions calculated in the last call
    ! is equal to the old end time, plus the duration of the boundary conditions
    ! created at the last call
    spectrumendtimeold = spectrumendtime ! Robert: need this in case new run with
                                         ! instat='reuse' does not have exact same
                                         ! time steps as original simulation
    spectrumendtime = spectrumendtime + wp%rtbc

    ! Collect new file identifiers for administration list files
    if (bccount==1) then
       if (par%nonhspectrum==0) then
          fidelist = create_new_fid()
          open(fidelist,file='ebcflist.bcf',form='formatted',position='append',status='replace')
          fidqlist = create_new_fid()
          open(fidqlist,file='qbcflist.bcf',form='formatted',position='append',status='replace')
       else
          fidnhlist = create_new_fid()
          open(fidnhlist,file='nhbcflist.bcf',form='formatted',status='replace',iostat=iostat)
       endif
    else
       if (par%nonhspectrum==0) then
          fidelist = create_new_fid()
          open(fidelist,file='ebcflist.bcf',form='formatted',position='append')
          fidqlist = create_new_fid()
          open(fidqlist,file='qbcflist.bcf',form='formatted',position='append')
       else
          fidnhlist = create_new_fid()
          open(fidnhlist,file='nhbcflist.bcf',form='formatted',position='append',iostat=iostat)
       endif
    endif
    ! Write new line
    if (par%nonhspectrum==0) then
       write(fidelist,'(f12.3,a,f12.3,a,f9.3,a,f9.5,a,f11.5,a)') &
            & spectrumendtime,' ',wp%rtbc,' ',wp%dtbc,' ',par%Trep,' ',combspec%dir0,' '//trim(wp%Efilename)
       write(fidqlist,'(f12.3,a,f12.3,a,f9.3,a,f9.5,a,f11.5,a)') &
            & spectrumendtime,' ',wp%rtbc,' ',wp%dtbc,' ',par%Trep,' ',combspec%dir0,' '//trim(wp%qfilename)
       ! Close administation files
       close(fidelist)
       close(fidqlist)
    else
       write(fidnhlist,'(f12.3,a,f12.3,a)',iostat=iostat) &
           & spectrumendtimeold,'   ',spectrumendtime,' '//trim(wp%nhfilename)
       close(fidnhlist)
    endif

    ! Set output of the correct line
    curline = bccount
  endif

  endsubroutine spectral_wave_bc

  ! --------------------------------------------------------------
  ! ---------------- Read input spectra files --------------------
  ! --------------------------------------------------------------
  subroutine read_spectrum_input(par,wp,fn,specin)
    use params
    use filefunctions
    use logging_module

    implicit none
    ! Interface
    type(parameters),intent(in)    :: par
    type(waveparamsnew),intent(inout) :: wp
    type(filenames),intent(inout)  :: fn
    type(spectrum),intent(inout)   :: specin
    ! internal
    integer                     :: fid
    character(8)                :: testline
    character(slen)             :: readfile
    logical                     :: filelist
    integer                     :: i,ier

    ! check the first line of the boundary condition file for FILELIST keyword
    fid = create_new_fid()
    open(fid,file=fn%fname,status='old',form='formatted')
    read(fid,*,iostat=ier)testline
    if (ier .ne. 0) then
       call report_file_read_error(fn%fname)
    endif
    if (trim(testline)=='FILELIST') then
       filelist = .true.
       ! move listline off its default position of zero to the first row
       if (fn%listline==0) then
          fn%listline = 1
       endif
       fn%reuse = .false.
    else
       filelist = .false.
       if (trim(par%instat) /= 'jons_table') then
          fn%reuse = .true.
       else
          fn%reuse = .false.
       endif
    endif
    close(fid)
    ! If file has list of spectra, read through the lines to find the correct
    ! spectrum file name and rtbc and dtbc. Else the filename and rtbc, dtbc are in
    ! params.txt
    if (filelist) then
       fid = create_new_fid()
       open(fid,file=fn%fname,status='old',form='formatted')
       do i=1,fn%listline
          read(fid,*)testline  ! old stuff, not needed anymore
       enddo
       read(fid,*,iostat=ier)wp%rtbc,wp%dtbc,readfile  ! new boundary condition
       if (ier .ne. 0) then
          call report_file_read_error(fn%fname)
       endif
       ! we have to adjust this to morphological time, as done in params.txt
       if (par%morfacopt==1) then
          wp%rtbc = wp%rtbc/max(par%morfac,1.d0)
       endif
       fn%listline = fn%listline + 1   ! move one on from the last time we opened this file
       close(fid)
    else
       wp%rtbc = par%rt  ! already set to morphological time
       wp%dtbc = par%dtbc
       readfile = fn%fname
    endif

    ! based on the value of instat, we need to read either Jonswap, Swan or vardens files
    ! note: jons_table is also handeled by read_jonswap_file subroutine
    select case (par%instat(1:4))
    case ('jons')
       ! wp type sent in to receive rtbc and dtbc from jons_table file
       ! fn%listline sent in to find correct row in jons_table file
       ! pfff..
       call read_jonswap_file(par,wp,readfile,fn%listline,specin)
    case ('swan')
       call read_swan_file(par,readfile,specin)
    case ('vard')
       call read_vardens_file(par,readfile,specin)
    endselect

  endsubroutine read_spectrum_input

  ! --------------------------------------------------------------
  ! ------------------- Read JONSWAP files -----------------------
  ! --------------------------------------------------------------
  subroutine read_jonswap_file(par,wp,readfile,listline,specin)

    use readkey_module
    use params
    use logging_module
    use filefunctions

    IMPLICIT NONE

    ! Input / output variables
    type(parameters), INTENT(IN)            :: par
    type(waveparamsnew),intent(inout)       :: wp
    character(len=*), intent(IN)            :: readfile
    integer, intent(INOUT)                  :: listline
    type(spectrum),intent(inout)            :: specin

    ! Internal variables
    integer                                 :: i,ii,ier,ip,ind
    integer                                 :: nmodal
    integer                                 :: fid
    integer                                 :: forcepartition
    integer,dimension(2)                    :: indvec
    real*8,dimension(:),allocatable         :: x, y, Dd, tempdir
    real*8,dimension(:),allocatable         :: Hm0,fp,gam,mainang,scoeff
    real*8                                  :: dfj, fnyq
    real*8                                  :: gammajsp,Tp
    character(len=80)                       :: dummystring
    type(spectrum),dimension(:),allocatable :: multinomalspec,scaledspec
    logical                                 :: cont
    real*8,dimension(:),allocatable         :: scalefac1,scalefac2,tempmax,avgscale
    real*8,dimension(:),allocatable         :: oldvariance,newvariance
    real*8                                  :: newconv,oldconv

    ! First part: read JONSWAP parameter data

    ! Check whether spectrum characteristics or table should be used
    if (trim(par%instat) /= 'jons_table') then
       ! Use spectrum characteristics
       call writelog('sl','','waveparams: Reading from ',trim(readfile),' ...')
       !
       ! First read if the spectrum is multinodal, and how many peaks there should be
       !
       nmodal     = readkey_int (readfile, 'nmodal',  1,  1, 4, bcast=.false.)
       if (nmodal<1) then
          call writelog('lswe','','Error: number of spectral partions may not be less than 1 in ',trim(readfile))
          call halt_program
       endif
       !
       ! Allocate space for all spectral parameters
       !
       allocate(Hm0(nmodal))
       allocate(fp(nmodal))
       allocate(gam(nmodal))
       allocate(mainang(nmodal))
       allocate(scoeff(nmodal))
       !
       ! Read the spectral parameters for all spectrum components
       !
       ! Wave height (required)
       Hm0    = readkey_dblvec(readfile, 'Hm0',nmodal,nmodal, 0.0d0, 0.0d0, 5.0d0, bcast=.false.,required=.true. )
       !
       ! Wave period (required)
       ! allow both Tp and fp specification to bring in line with params.txt
       if (isSetParameter(readfile,'Tp',bcast=.false.) .and. .not. isSetParameter(readfile,'fp',bcast=.false.)) then
          fp     = 1.d0/readkey_dblvec(readfile, 'Tp',nmodal,nmodal, 12.5d0, 2.5d0, 20.0d0, bcast=.false.)
       elseif (isSetParameter(readfile,'fp',bcast=.false.) .and. .not. isSetParameter(readfile,'Tp',bcast=.false.)) then
          fp     = readkey_dblvec(readfile, 'fp',nmodal,nmodal, 0.08d0,0.0625d0,   0.4d0, bcast=.false.)
       elseif (.not. isSetParameter(readfile,'fp',bcast=.false.) .and. .not. isSetParameter(readfile,'Tp',bcast=.false.)) then
          call writelog('lswe','','Error: missing required value for parameter ''Tp'' or ''fp'' in ',trim(readfile))
          call halt_program
       else
          fp     = 1.d0/readkey_dblvec(readfile, 'Tp',nmodal,nmodal, 12.5d0, 2.5d0, 20.0d0, bcast=.false.)
          call writelog('lsw','','Warning: selecting to read peak period (Tp) instead of frequency (fp) in ',trim(readfile))
       endif
       !
       ! Wave spreading in frequency domain (peakedness)
       !
       gam    = readkey_dblvec(readfile, 'gammajsp',nmodal,nmodal, 3.3d0, 1.0d0, 5.0d0, bcast=.false.)
       !
       ! Wave spreading in directional domain
       !
       scoeff = readkey_dblvec(readfile, 's',nmodal,nmodal, 10.0d0, 1.0d0, 1000.0d0, bcast=.false.)
       !
       ! Main wave direction
       ! allow both mainang and dir0 specification to bring in line with params.txt
       if (isSetParameter(readfile,'mainang',bcast=.false.) .and. .not. isSetParameter(readfile,'dir0',bcast=.false.)) then
          mainang = readkey_dblvec(readfile, 'mainang',nmodal,nmodal, 270.0d0, 0.0d0, 360.0d0, bcast=.false.)
       elseif (isSetParameter(readfile,'dir0',bcast=.false.) .and. .not. isSetParameter(readfile,'mainang',bcast=.false.)) then
          mainang = readkey_dblvec(readfile, 'dir0',nmodal,nmodal, 270.0d0, 0.0d0, 360.0d0, bcast=.false.)
       elseif (.not. isSetParameter(readfile,'dir0',bcast=.false.) .and. .not. isSetParameter(readfile,'mainang',bcast=.false.)) then
          mainang = 270.d0
       else
          mainang = readkey_dblvec(readfile, 'mainang',nmodal,nmodal, 270.0d0, 0.0d0, 360.0d0, bcast=.false.)
          call writelog('lsw','','Warning: selecting to read ''mainang'' instead of ''dir0'' in ',trim(readfile))
       endif
       !
       ! Nyquist parameters used only in this subroutine
       ! are not read individually for each spectrum partition
       if (par%oldnyq==1) then
          fnyq    = readkey_dbl (readfile, 'fnyq',       0.3d0,    0.2d0,      1.0d0,      bcast=.false. )
       else
         fnyq    = readkey_dbl (readfile, 'fnyq',max(0.3d0,3.d0*maxval(fp)), 0.2d0, 1.0d0, bcast=.false. )
       endif
       dfj        = readkey_dbl (readfile, 'dfj',      fnyq/200,   fnyq/1000,  fnyq/20,    bcast=.false. )
       !
       ! Finally check if XBeach should accept even the most stupid partioning (sets error level to warning
       ! level when computing partition overlap
       if (nmodal>1) then
          forcepartition = readkey_int (readfile, 'forcepartition',  0,  0, 1, bcast=.false.)
       endif
       ! check for other strange values in this file
       call readkey(readfile,'checkparams',dummystring)
    else
       nmodal = 1
       allocate(Hm0(nmodal))
       allocate(fp(nmodal))
       allocate(gam(nmodal))
       allocate(mainang(nmodal))
       allocate(scoeff(nmodal))
       ! Use spectrum table
       fid = create_new_fid()
       call writelog('sl','','waveparams: Reading from table ',trim(readfile),' ...')
       open(fid,file=readfile,status='old',form='formatted')
       ! read junk up to the correct line in the file
       do i=1,listline
          read(fid,*,iostat=ier)dummystring
          if (ier .ne. 0) then
             call report_file_read_error(readfile)
          endif
       enddo
       read(fid,*,iostat=ier)Hm0(1),Tp,mainang(1),gam(1),scoeff(1),wp%rtbc,wp%dtbc
       if (ier .ne. 0) then
          call writelog('lswe','','Error reading file ',trim(readfile))
          close(fid)
          call halt_program
       endif
       ! move the line pointer in the file
       listline = listline+1
       ! convert to morphological time
       if (par%morfacopt==1) then
          wp%rtbc = wp%rtbc/max(par%morfac,1.d0)
       endif
       fp(1)=1.d0/Tp
       fnyq=3.d0*fp(1)
       dfj=fp(1)/50
       close(fid)
    endif
    !
    ! Second part: generate 2D spectrum from input parameters
    !
    ! Define number of frequency bins by defining an array of the necessary length
    ! using the Nyquist frequency and frequency step size
    specin%nf = ceiling((fnyq-dfj)/dfj)
    specin%df = dfj
        !
    ! Define array with actual eqidistant frequency bins
    allocate(specin%f(specin%nf))
    do i=1,specin%nf
       specin%f(i)=i*dfj
    enddo
    !
    ! we need a normalised frequency and variance vector for JONSWAP generation
    allocate(x(size(specin%f)))
    allocate(y(size(specin%f)))
        !
    ! Define 200 directions relative to main angle running from 0 to 2*pi
    ! we also need a temporary vector for direction distribution
    specin%nang = naint
    allocate(tempdir(specin%nang))
    allocate(Dd(specin%nang))
    allocate(specin%ang(specin%nang))
    specin%dang = 2*par%px/(naint-1)
    do i=1,specin%nang
       specin%ang(i)=(i-1)*specin%dang
    enddo
    !
    ! Generate density spectrum for each spectrum partition
    !
    allocate(multinomalspec(nmodal))
    !
    do ip=1,nmodal
       ! relative frequenct vector
       x=specin%f/fp(ip)
       ! Calculate unscaled and non-directional JONSWAP spectrum using
       ! peak-enhancement factor and pre-determined frequency bins
       call jonswapgk(x,gam(ip),y)
       ! Determine scaled and non-directional JONSWAP spectrum using the JONSWAP
       ! characteristics
       y=(Hm0(ip)/(4.d0*sqrt(sum(y)*dfj)))**2*y
       ! Convert main angle from degrees to radians and from nautical convention to
       ! internal grid
       mainang(ip)=(1.5d0*par%px)-mainang(ip)*par%px/180
       ! Make sure the main angle is defined between 0 and 2*pi
       do while (mainang(ip)>2*par%px .or. mainang(ip)<0.d0) !Robert en Ap
          if (mainang(ip)>2*par%px) then
             mainang(ip)=mainang(ip)-2*par%px
          elseif (mainang(ip)<0.d0) then
             mainang(ip)=mainang(ip)+2*par%px
          endif
       enddo
       ! Convert 200 directions relative to main angle to directions relative to
       ! internal grid
       ! Bas: apparently division by 2 for cosine law happens already here
       tempdir = (specin%ang-mainang(ip))/2
       ! Make sure all directions around the main angle are defined between 0 and 2*pi
       do while (any(tempdir>2*par%px) .or. any(tempdir<0.d0))
          where (tempdir>2*par%px)
             tempdir=tempdir-2*par%px
          elsewhere (tempdir<0.d0)
             tempdir=tempdir+2*par%px
          endwhere
       enddo
       ! Calculate directional spreading based on cosine law
       Dd = dcos(tempdir)**(2*nint(scoeff(ip)))  ! Robert: apparently nint is needed here, else MATH error
       ! Scale directional spreading to have a surface of unity by dividing by its
       ! own surface
       Dd = Dd / (sum(Dd)*specin%dang)
       ! Define two-dimensional variance density spectrum array and distribute
       ! variance density for each frequency over directional bins
       allocate(multinomalspec(ip)%S(specin%nf,specin%nang))
       do i=1,specin%nang
          do ii=1,specin%nf
             multinomalspec(ip)%S(ii,i)=y(ii)*Dd(i)
          end do
       end do
     enddo  ! 1,nmodal
     do ip=1,nmodal
        write(*,*)'Hm0 = ',4*sqrt(sum(multinomalspec(ip)%S)*specin%dang*specin%df)
     enddo
     !
     ! Combine spectrum partitions so that the total spectrum is correct
     !
     if (nmodal==1) then
        ! Set all the useful parameters and arrays
        specin%Hm0 = Hm0(1)
        specin%fp = fp(1)
        specin%dir0 = mainang(1)
        specin%scoeff = scoeff(1)
        allocate(specin%S(specin%nf,specin%nang))
        specin%S = multinomalspec(1)%S
     else
        specin%Hm0 = sqrt(sum(Hm0**2))  ! total wave height
        ind = maxval(maxloc(Hm0))       ! spectrum that is largest
        specin%fp = fp(ind)             ! not really used in further calculation
        specin%dir0 = mainang(ind)      ! again not really used in further calculation
        ! if all scoeff>1000 then all waves should be in the same direction exactly
        if (all(scoeff>=1024.d0) .and. all(mainang==mainang(ind))) then
           specin%scoeff = 1024.d0
        else
           specin%scoeff = min(scoeff(ind),999.d0)
        endif
        ! Now we have to loop over all partitioned spectra. Where two or more spectra
        ! overlap, only the largest is counted, and all others are set to zero. Afterwards
        ! all spectra are scaled so that the total energy is maintained. Since the scaling
        ! affects where spectra overlap, this loop is repeated until only minor changes
        ! in the scaling occur. Warnings or errors are given if the spectrum is scaled too
        ! much, i.e. spectra overlap too much.
        !
        ! Allocate space
        allocate(scaledspec(nmodal)) ! used to store scaled density spectra
        do ip=1,nmodal
           allocate(scaledspec(ip)%S(specin%nf,specin%nang))
           scaledspec(ip)%S = multinomalspec(ip)%S
        enddo
        allocate(scalefac1(nmodal))  ! this is the scaling factor required to maintain Hm0
                                     ! including that parts of the spectrum are set to zero
        scalefac1 = 1.d0
        allocate(scalefac2(nmodal))  ! this is the scaling factor required to maintain Hm0
                                     ! in the previous iteration
        scalefac2 = 1.d0
        allocate(avgscale(nmodal))
        avgscale = 1.d0
        ! these are convergence criteria
        newconv = 0.d0
        oldconv = huge(0.d0)
        cont = .true.
        allocate(tempmax(nmodal))     ! used to store maximum value in f,theta space
        allocate(oldvariance(nmodal)) ! used to store the sum of variance in the original partitions
        do ip=1,nmodal
           oldvariance(ip) = sum(multinomalspec(ip)%S)
        enddo
        allocate(newvariance(nmodal)) ! used to store the sum of variance in the new partitions
        !
        ! Start convergence loop
        do while (cont)
           avgscale = (avgscale+scalefac1)/2
           ! First scale the spectra
           do ip=1,nmodal
              ! scale the spectrum using less than half the additional scale factor, this to
              ! ensure the method does not become unstable
              ! scaledspec(ip)%S = multinomalspec(ip)%S*(0.51d0+scalefac1(ip)*0.49d0)
              ! scaledspec(ip)%S = multinomalspec(ip)%S*(scalefac1(ip)+scalefac2(ip))/2
              scaledspec(ip)%S = multinomalspec(ip)%S*avgscale(ip)

           enddo
           do i=1,specin%nang
              do ii=1,specin%nf
                 ! vector of variance densities at this point in f,theta space
                 do ip=1,nmodal
                    tempmax(ip) = scaledspec(ip)%S(ii,i)
                 end do
                 ! All spectra other than the one that is greatest in this point
                 ! is set to zero. In case multiple spectra are equal largest,
                 ! the first spectrum in the list is chosen (minval(maxloc))
                 ind = minval(maxloc(tempmax))
                 do ip=1,nmodal
                    if (ip/=ind) then
                       scaledspec(ip)%S(ii,i) = 0.d0
                    endif
                 enddo
              enddo
           end do
           !
           ! Now rescale adjusted partition spectra so that they match the incident Hm0
           scalefac2 = scalefac1  ! keep previous results
           do ip=1,nmodal
              newvariance(ip) = sum(scaledspec(ip)%S)
              if (newvariance(ip)>0.01d0*oldvariance(ip)) then
                 scalefac1(ip) = oldvariance(ip)/newvariance(ip) ! want to maximise to a factor 2
                                                                    ! else can generate rediculous results
              else
                 scalefac1(ip) = 0.d0                            ! completely remove this spectrum
              endif
           enddo
           !
           ! check convergence criteria (if error is increasing, we have passed best fit)
           newconv = maxval(abs(scalefac2-scalefac1)/scalefac1)
           if (newconv<0.0001d0 .or. abs(newconv-oldconv)<0.0001d0) then
              cont = .false.
           endif
           oldconv = newconv
        end do
        ! Ensure full energy conservation by scaling now according to full scalefac, not just the
        ! half used in the iteration loop
        do ip=1,nmodal
           scaledspec(ip)%S = multinomalspec(ip)%S*scalefac1(ip)
        enddo
        !
        ! Now check the total scaling that has taken place across the spectrum, also accounting
        ! for the fact that some parts of the spectra are set to zero. This differs from the
        ! scalefac1 and scalefac2 vectors in the loop that only adjust the part of the spectrum
        ! that is not zero.
        do ip=1,nmodal
           write(*,*)'Hm0m = ',4*sqrt(sum(scaledspec(ip)%S)*specin%dang*specin%df)
        enddo
        do ip=1,nmodal
           indvec = maxloc(scaledspec(ip)%S)
           scalefac1(ip) = scaledspec(ip)%S(indvec(1),indvec(2)) / &
                           multinomalspec(ip)%S(indvec(1),indvec(2))-1.d0
        enddo
        !
        ! Warning and/or error criteria here if spectra overlap each other too much
        do ip=1,nmodal
           if(scalefac1(ip)>0.5d0) then
              if (forcepartition==1) then
                 call writelog('lsw','(a,f0.0,a,i0,a,a,a)', &
                               'Warning: ',scalefac1(ip)*100,'% of energy in spectrum partition ''',ip, &
                               ''' in  ', trim(readfile),' is overlapped by other partitions')
                 call writelog('lsw','',' Check spectral partitioning in ',trim(readfile))
              else
                 call writelog('lswe','(a,f0.0,a,i0,a,a,a)', &
                               'Error: ',scalefac1(ip)*100,'% of energy in spectrum partition ''',ip, &
                               ''' in  ', trim(readfile),' is overlapped by other partitions')
                 call writelog('lswe','','This spectrum overlaps too much with another spectrum partition.', &
                                         ' Check spectral partitioning in ',trim(readfile))
                 call writelog('lswe','','If the partitioning should be carried out regardles of energy loss, ', &
                                          'set ''forcepartition = 1'' in ',trim(readfile))
                 call halt_program
              endif
           elseif (scalefac1(ip)>0.2d0 .and. scalefac1(ip)<=0.5d0) then
              call writelog('lsw','(a,f0.0,a,i0,a,a,a)', &
                            'Warning: ',scalefac1(ip)*100,'% of energy in spectrum partition ''',ip, &
                            ''' in  ', trim(readfile),' is overlapped by other partitions')
              call writelog('lsw','',' Check spectral partitioning in ',trim(readfile))
           elseif (scalefac1(ip)<0.d0) then
              call writelog('lsw','(a,i0,a,a,a)','Warning: spectrum partition ''',ip,''' in  ',trim(readfile), &
                                   ' has been removed')
              call writelog('lsw','','This spectrum is entirely overlapped by another spectrum partition.',&
                                     ' Check spectral partitioning in ',trim(readfile))
           endif
        enddo
        !
        ! Now set total spectrum
        allocate(specin%S(specin%nf,specin%nang))
        specin%S = 0.d0
        do ip=1,nmodal
           specin%S = specin%S+scaledspec(ip)%S
        enddo

        deallocate(scaledspec)
        deallocate(tempmax)
        deallocate(scalefac1,scalefac2)
        deallocate(newvariance,oldvariance)
     endif
     
    ! We need frequency spectrum to ensure Sf remains correct between interpoloation
    ! routines
    allocate(specin%Sf(specin%nf))
    specin%Sf=0.d0
    do i=1,specin%nf
       do ii=1,specin%nang
          if (ii==1) then
             specin%Sf(i) = specin%Sf(i)+specin%S(i,ii)*abs(specin%ang(2)-specin%ang(1))
          elseif (ii==specin%nang) then
             specin%Sf(i) = specin%Sf(i)+specin%S(i,ii)*abs(specin%ang(ii)-specin%ang(ii-1))
          else
             specin%Sf(i) = specin%Sf(i)+specin%S(i,ii)*abs(specin%ang(ii+1)-specin%ang(ii-1))/2
          endif
       enddo
    enddo

    deallocate (Dd)
    deallocate(Hm0,fp,gam,mainang,scoeff)
    deallocate (x,y)
    deallocate(tempdir)
    ! TODO dereference pointers first....
    do ip=1,nmodal
       deallocate(multinomalspec(ip)%S)
    end do
    deallocate(multinomalspec)


  end subroutine read_jonswap_file



  ! -----------------------------------------------------------
  ! --------- JONSWAP  unscaled JONSWAP spectrum --------------
  ! -------------(used by read_jonswap_files)------------------
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


  ! --------------------------------------------------------------
  ! ----------------------Read SWAN files ------------------------
  ! --------------------------------------------------------------
  subroutine read_swan_file(par,readfile,specin)

    use params
    use logging_module
    use filefunctions
    use xmpi_module, only: Halt_Program
    use math_tools, only: flipv,flipa

    IMPLICIT NONE

    ! Input / output variables
    type(parameters), INTENT(IN)            :: par
    character(len=*), intent(IN)            :: readfile
    type(spectrum),intent(inout)            :: specin

    ! Internal variables
    character(6)                            :: rtext
    real*8                                  :: factor,exc
    integer                                 :: fid,switch
    integer                                 :: i,ii,ier,ier2,ier3
    logical                                 :: flipped
    integer                                 :: nt,Ashift
    real*8, dimension(:),allocatable        :: temp
    real*8, dimension(:,:),allocatable      :: tempA

    flipped =.false.
    switch  = 0

    call writelog('sl','','Reading from SWAN file ',trim(readfile),' ...')
    fid = create_new_fid()
    open(fid,file=readfile,form='formatted',status='old')

    ! Read file until RFREQ or AFREQ is found
    do while (switch==0)
       read(fid,'(a)',iostat=ier)rtext
       if (ier .ne. 0) then
          call report_file_read_error(readfile)
       endif
       if (rtext == 'RFREQ ') then
          switch = 1
       elseif (rtext == 'AFREQ ') then
          switch = 2
       end if
    end do

    ! Read nfreq and f
    ! Note f is not monotonically increasing in most simulations
    read(fid,*,iostat=ier)specin%nf
    if (ier .ne. 0) then
       call report_file_read_error(readfile)
    endif
    allocate(specin%f(specin%nf))
    do i=1,specin%nf
       read(fid,*,iostat=ier)specin%f(i)
       if (ier .ne. 0) then
          call report_file_read_error(readfile)
       endif
    end do

    ! Convert to absolute frequencies:
    ! STILL TO BE DONE
    if (switch == 1) then
       specin%f = specin%f
    else
       specin%f = specin%f
    end if

    ! Read CDIR or NDIR
    read(fid,'(a)',iostat=ier)rtext
    if (ier .ne. 0) then
       call report_file_read_error(readfile)
    endif
    if (rtext == 'NDIR  ') then
       switch = 1
    elseif (rtext == 'CDIR  ') then
       switch = 2
    else
       call writelog('ewls','', 'SWAN directional bins keyword not found')
       call halt_program
    endif

    ! Read ndir, theta
    read(fid,*,iostat=ier)specin%nang
    if (ier .ne. 0) then
       call report_file_read_error(readfile)
    endif
    allocate(specin%ang(specin%nang))
    do i=1,specin%nang
       read(fid,*,iostat=ier)specin%ang(i)
       if (ier .ne. 0) then
          call report_file_read_error(readfile)
       endif
    end do

    ! Convert angles to cartesian degrees relative to East
    if (switch == 1) then
       ! nautical to cartesian East
       specin%ang = 270.d0-specin%ang
    else
       ! cartesian to cartesian East
       specin%ang = specin%ang-par%dthetaS_XB
       ! dthetaS_XB is the (counter-clockwise) angle in the degrees to rotate from the x-axis in SWAN to the
       ! x-axis pointing East
    end if

    ! Ensure angles are increasing instead of decreasing
    if (specin%ang(2)<specin%ang(1)) then
       call flipv(specin%ang,size(specin%ang))
       flipped=.true.
    end if

    nt = 0
    Ashift = 0
    ! Make sure that all angles are in range of 0 to 360 degrees
    if(minval(specin%ang)<0.d0)then
       allocate (temp(specin%nang))
       Ashift=-1
       temp=0.d0
       do i=1,specin%nang
          if (specin%ang(i)<0.d0) then
             specin%ang(i)=specin%ang(i)+360.0d0
             nt = nt+1
          endif
       enddo
       temp(1:specin%nang-nt)=specin%ang(nt+1:specin%nang)
       temp(specin%nang-nt+1:specin%nang)=specin%ang(1:nt)
       specin%ang=temp
       deallocate(temp)
    elseif(maxval(specin%ang)>360.0d0)then
       allocate (temp(specin%nang))
       Ashift=1
       temp=0.d0
       do i=1,specin%nang
          if (specin%ang(i)>360.d0) then
             specin%ang(i)=specin%ang(i)-360.0d0
             nt = nt+1
          endif
       enddo
       temp(nt+1:specin%nang)=specin%ang(1:specin%nang-nt)
       temp(1:nt)=specin%ang(specin%nang-nt+1:specin%nang)
       specin%ang=temp
       deallocate(temp)
    endif

    ! convert to radians
    specin%ang=specin%ang*par%px/180
    specin%dang=specin%ang(2)-specin%ang(1)

    ! Skip Quant, next line, read VaDens or EnDens
    read(fid,'(a)',iostat=ier)rtext
    read(fid,'(a)',iostat=ier2)rtext
    read(fid,'(a)',iostat=ier3)rtext
    if (ier+ier2+ier3 .ne. 0) then
       call report_file_read_error(readfile)
    endif
    if (rtext == 'VaDens') then
       switch = 1
    elseif (rtext == 'EnDens') then
       switch = 2
    else
       call writelog('slwe','', 'SWAN VaDens/EnDens keyword not found')
       call halt_program
    end if
    read(fid,'(a)',iostat=ier)rtext
    read(fid,*,iostat=ier2)exc
    if (ier+ier2 .ne. 0) then
       call report_file_read_error(readfile)
    endif

    i=0
    ! Find FACTOR keyword
    do while (i==0)
       read(fid,'(a)',iostat=ier)rtext
       if (ier .ne. 0) then
          call report_file_read_error(readfile)
       endif
       if (rtext == 'FACTOR') then
          i=1
       elseif (rtext == 'ZERO  ') then
          call writelog('lswe','','Zero energy density input for this point')
          call halt_program
       elseif (rtext == 'NODATA') then
          call writelog('lwse','','SWAN file has no data for this point')
          call halt_program
       end if
    end do
    read(fid,*,iostat=ier)factor
    if (ier .ne. 0) then
       call report_file_read_error(readfile)
    endif

    ! Read 2D S array
    allocate(specin%S(specin%nf,specin%nang))
    do i=1,specin%nf
       read(fid,*,iostat=ier)(specin%S(i,ii),ii=1,specin%nang)
       if (ier .ne. 0) then
          call report_file_read_error(readfile)
       endif
    end do

    ! Finished reading file
    close(fid)

    ! Replace exception value
    where (specin%S == exc)
       specin%S = 0.d0
    endwhere

    ! If angles were decreasing, flip S_array as also dir is flipped
    if (flipped) then
       call flipa(specin%S,specin%nf,specin%nang,2)
    end if

    ! If the order of the angles in specin%ang was reordered, so the same in
    ! specin%S array
    if(Ashift==-1)then
       allocate(tempA(specin%nf,specin%nang))
       tempA=0
       tempA(:,1:specin%nang-nt)=specin%S(:,nt+1:specin%nang)
       tempA(:,specin%nang-nt+1:specin%nang)=specin%S(:,1:nt)
       specin%S=tempA
       deallocate(tempA)
    elseif (Ashift==1) then
       allocate(tempA(specin%nf,specin%nang))
       tempA=0
       tempA(:,nt+1:specin%nang)=specin%S(:,1:specin%nang-nt)
       tempA(:,1:nt)=specin%S(:,specin%nang-nt+1:specin%nang)
       specin%S=tempA
       deallocate(tempA)
    endif

    ! multiply by SWAN output factor
    specin%S=specin%S*factor

    ! Convert from energy density to variance density
    if (switch == 2) then
       specin%S=specin%S/(par%rho*par%g)
    end if

    ! Convert to m2/Hz/rad
    specin%S=specin%S*180/par%px

    ! We need a value for spreading. The assumption is that it is less than 1000
    ! This way, wp%fgen will not be set to just one angle.
    specin%scoeff = -1.d0
    ! We need to know if hm0 was set explicitly, not the case for Swan files
    specin%hm0 = -1.d0
    
    ! We need frequency spectrum to ensure Sf remains correct between interpoloation
    ! routines
    allocate(specin%Sf(specin%nf))
    specin%Sf=0.d0
    do i=1,specin%nf
       do ii=1,specin%nang
          if (ii==1) then
             specin%Sf(i) = specin%Sf(i)+specin%S(i,ii)*abs(specin%ang(2)-specin%ang(1))
          elseif (ii==specin%nang) then
             specin%Sf(i) = specin%Sf(i)+specin%S(i,ii)*abs(specin%ang(ii)-specin%ang(ii-1))
          else
             specin%Sf(i) = specin%Sf(i)+specin%S(i,ii)*abs(specin%ang(ii+1)-specin%ang(ii-1))/2
          endif
       enddo
    enddo

  end subroutine read_swan_file


  ! --------------------------------------------------------------
  ! -----------------Read variance density files -----------------
  ! --------------------------------------------------------------
  subroutine read_vardens_file(par,readfile,specin)

    use params
    use logging_module
    use filefunctions
    use xmpi_module, only: Halt_Program
    use math_tools, only: flipv,flipa

    IMPLICIT NONE

    ! Input / output variables
    type(parameters), INTENT(IN)            :: par
    character(len=*), intent(IN)            :: readfile
    type(spectrum),intent(inout)            :: specin

    ! Internal variables
    integer                                 :: fid,i,ii,nnz,ier,nt,Ashift
    real*8,dimension(:),allocatable         :: Sd,temp
    real*8, dimension(:,:),allocatable      :: tempA

    ! Open file to start read
    call writelog('sl','','Reading from vardens file ',trim(readfile),' ...')
    fid = create_new_fid()
    open(fid,file=readfile,form='formatted',status='old')

    ! Read number of frequencies and frequency vector
    read(fid,*,iostat=ier)specin%nf
    if (ier .ne. 0) then
       call report_file_read_error(readfile)
    endif
    allocate(specin%f(specin%nf))
    do i=1,specin%nf
       read(fid,*,iostat=ier)specin%f(i)
       if (ier .ne. 0) then
          call report_file_read_error(readfile)
       endif
    end do

    ! Read number of angles and angles vector
    read(fid,*,iostat=ier)specin%nang
    if (ier .ne. 0) then
       call report_file_read_error(readfile)
    endif
    allocate(specin%ang(specin%nang))
    do i=1,specin%nang
       read(fid,*,iostat=ier)specin%ang(i)
       if (ier .ne. 0) then
          call report_file_read_error(readfile)
       endif
    end do
    ! Directions should span 0-360 degrees, not negative
    nt = 0
    Ashift = 0
    ! Make sure that all angles are in range of 0 to 360 degrees
    if(minval(specin%ang)<0.d0)then
       allocate (temp(specin%nang))
       Ashift=-1
       temp=0.d0
       do i=1,specin%nang
          if (specin%ang(i)<0.d0) then
             specin%ang(i)=specin%ang(i)+360.0d0
             nt = nt+1
          endif
       enddo
       temp(1:specin%nang-nt)=specin%ang(nt+1:specin%nang)
       temp(specin%nang-nt+1:specin%nang)=specin%ang(1:nt)
       specin%ang=temp
       deallocate(temp)
    elseif(maxval(specin%ang)>360.0d0)then
       allocate (temp(specin%nang))
       Ashift=1
       temp=0.d0
       do i=1,specin%nang
          if (specin%ang(i)>360.d0) then
             specin%ang(i)=specin%ang(i)-360.0d0
             nt = nt+1
          endif
       enddo
       temp(nt+1:specin%nang)=specin%ang(1:specin%nang-nt)
       temp(1:nt)=specin%ang(specin%nang-nt+1:specin%nang)
       specin%ang=temp
       deallocate(temp)
    endif

    ! Convert from degrees to rad
    specin%ang=specin%ang*par%px/180
    specin%dang=specin%ang(2)-specin%ang(1)

    ! Read 2D S array
    allocate(specin%S(specin%nf,specin%nang))
    do i=1,specin%nf
       read(fid,*,iostat=ier)(specin%S(i,ii),ii=1,specin%nang)
       if (ier .ne. 0) then
          call report_file_read_error(readfile)
       endif
    end do

    ! If the order of the angles in specin%ang was reordered, so the same in
    ! specin%S array
    if(Ashift==-1)then
       allocate(tempA(specin%nf,specin%nang))
       tempA=0
       tempA(:,1:specin%nang-nt)=specin%S(:,nt+1:specin%nang)
       tempA(:,specin%nang-nt+1:specin%nang)=specin%S(:,1:nt)
       specin%S=tempA
       deallocate(tempA)
    elseif (Ashift==1) then
       allocate(tempA(specin%nf,specin%nang))
       tempA=0
       tempA(:,nt+1:specin%nang)=specin%S(:,1:specin%nang-nt)
       tempA(:,1:nt)=specin%S(:,specin%nang-nt+1:specin%nang)
       specin%S=tempA
       deallocate(tempA)
    endif

    ! Finished reading file
    close(fid)

    ! Convert to m2/Hz/rad
    specin%S=specin%S*180/par%px

    ! We need a value for spreading. The assumption is that it is less than 1000
    ! if more than one direction column has variance. If only one column has variance
    ! we assume the model requires long crested waves, and wp%thetagen should be exactly
    ! equal to the input direction. We then also need to find the dominant angle.
    ! First flatten spectrum to direction only, then determine if all but one are zero,
    ! use the direction of the non-zero as the main wave direction, set spreading to
    ! a high value (greater than 1000). Else set spreading to a negative value.
    allocate(Sd(specin%nang))
    Sd = sum(specin%S,DIM = 1)
    nnz = count(Sd>0.d0)
    if (nnz == 1) then
       specin%scoeff = 1024.d0
       specin%dir0 = specin%ang(minval(maxloc(Sd)))
    else
       specin%scoeff = -1.d0
    endif
    
    ! We need frequency spectrum to ensure Sf remains correct between interpoloation
    ! routines
    allocate(specin%Sf(specin%nf))
    specin%Sf=0.d0
    do i=1,specin%nf
       do ii=1,specin%nang
          if (ii==1) then
             specin%Sf(i) = specin%Sf(i)+specin%S(i,ii)*abs(specin%ang(2)-specin%ang(1))
          elseif (ii==specin%nang) then
             specin%Sf(i) = specin%Sf(i)+specin%S(i,ii)*abs(specin%ang(ii)-specin%ang(ii-1))
          else
             specin%Sf(i) = specin%Sf(i)+specin%S(i,ii)*abs(specin%ang(ii+1)-specin%ang(ii-1))/2
          endif
       enddo
    enddo
    

    ! We need to know if hm0 was set explicitly, not the case for vardens files
    specin%hm0 = -1.d0

    ! Free memory
    deallocate(Sd)
  end subroutine read_vardens_file


  ! --------------------------------------------------------------
  ! ------------- Interpolate to standard spectrum ---------------
  ! --------------------------------------------------------------
  subroutine interpolate_spectrum(specin,specinterp,par,fmax)

    use interp
    use params

    IMPLICIT NONE

    ! input/output
    type(spectrum),intent(in)            :: specin
    type(spectrum),intent(inout)         :: specinterp
    type(parameters),intent(in)          :: par
    real*8, intent(in)                   :: fmax
    ! internal
    integer                              :: i,j,dummy
    real*8                               :: m0,df,dang
    real*8,dimension(naint)              :: Sd
    real*8                               :: hm0pre,hm0post,Sfnow,factor,tempt0

    ! allocate size of f,ang,Sf and S arrays in specinterp
    allocate(specinterp%f(nfint))
    allocate(specinterp%ang(naint))
    allocate(specinterp%Sf(nfint))
    allocate(specinterp%S(nfint,naint))

    ! fill f and ang arrays
    specinterp%nf = nfint
    specinterp%df = fmax/(nfint-1)     ! this usually gives a range of 0-1Hz, 
                                       ! with 1/800Hz resolution
    do i=1,nfint
       specinterp%f(i)=(i-1)*specinterp%df
    enddo
    specinterp%nang = naint
    specinterp%dang = 2*par%px/(naint-1) ! this is exactly the same as in the JONSWAP construction
    do i=1,specinterp%nang
       specinterp%ang(i)=(i-1)*specinterp%dang
    enddo

    ! If hm0 was set explicitly, then use that, else calculate hm0
    if (specin%hm0>0.d0) then
       hm0pre = specin%hm0
    else
       ! pre-interpolation hm0 value (can be on a non-monotonic f,ang grid)
       m0 = 0
       do j=1,specin%nang
          do i=1,specin%nf
             df = specin%f(max(2,i))-specin%f(max(1,i-1))
             dang = specin%ang(max(2,j))-specin%ang(max(1,j-1))
             m0 = m0 + specin%S(i,j)*df*dang
          enddo
       enddo
       hm0pre = 4*sqrt(m0)
    endif

    ! interpolation (no extrapolation) of input 2D spectrum to standard 2D spectrum
    do j=1,naint
       do i=1,nfint
          call linear_interp_2d(specin%f,specin%nf,specin%ang,specin%nang,specin%S, &
               specinterp%f(i),specinterp%ang(j),specinterp%S(i,j),'interp',0.d0)
       enddo
    enddo
    ! hm0 post (is always on a monotonic f,ang grid
    m0 = sum(specinterp%S)*specinterp%df*specinterp%dang
    hm0post = 4*sqrt(m0)

    ! correct the wave energy to ensure hm0post == hm0pre
    ! specinterp%S = (hm0pre/hm0post)**2*specinterp%S
    ! This is done later using Sf
    
    ! calculate 1D spectrum, summed over directions
    specinterp%Sf = sum(specinterp%S, DIM = 2)*specinterp%dang
    
    ! correct the wave energy by setting Sfpost == Sfpre
    do i=1,nfint
       if (specinterp%f(i)>=minval(specin%f) .and. specinterp%f(i)<=maxval(specin%f)) then
          call LINEAR_INTERP(specin%f,specin%Sf,specin%nf,specinterp%f(i),Sfnow,dummy)
          if (specinterp%Sf(i)>0.d0 .and. Sfnow>0.d0) then
             factor = Sfnow/specinterp%Sf(i)
             specinterp%Sf(i)  = specinterp%Sf(i)*factor
             specinterp%S(i,:) = specinterp%S(i,:)*factor
          elseif (Sfnow==0.d0) then
             specinterp%Sf(i)  = 0.d0
             specinterp%S(i,:) = 0.d0
          else
             specinterp%Sf(i)  = Sfnow
             dummy = maxval(maxloc(specin%S(i,:)))
             tempt0 = specin%ang(dummy)
             dummy = maxval(minloc(abs(specinterp%ang-tempt0)))
             specinterp%S(i,dummy) = Sfnow/specinterp%dang
          endif
       endif
    enddo

    ! calculate other wave statistics from interpolated spectrum
    m0 = sum(specinterp%Sf)*specinterp%df
    specinterp%hm0 = 4*sqrt(m0)
    Sd = sum(specinterp%S, DIM = 1)*specinterp%df
    i  = maxval(maxloc(Sd))
    specinterp%dir0 = 270.d0 - specinterp%ang(i)*180/par%px  ! converted back into nautical degrees
    i  = maxval(maxloc(specinterp%Sf))
    specinterp%fp = specinterp%f(i)
    call tpDcalc(specinterp%Sf,specinterp%f,specinterp%trep,par%trepfac,par%Tm01switch)
    specinterp%dirm = 270.d0-180.d0/par%px*atan2( sum(sin(specinterp%ang)*Sd)/sum(Sd),&
         sum(cos(specinterp%ang)*Sd)/sum(Sd) )

  end subroutine interpolate_spectrum

  ! --------------------------------------------------------------
  ! ----------- Small subroutine to determine if the  ------------
  ! ---------- global reuseall should be true or false -----------
  subroutine set_reuseall

    implicit none
    ! internal
    integer                                           :: i

    ! set default to reuse all boundary conditions
    reuseall = .true.
    ! find any spectrum that cannot be reused, then change global
    ! reuseall to false
    do i=1,nspectra
       if (.not. bcfiles(i)%reuse) then
          reuseall = .false.
       endif
    enddo

  end subroutine set_reuseall

  ! --------------------------------------------------------------
  ! ----------- Small subroutine to set the filenames ------------
  ! ----------- of the boundary condition output files -----------
  subroutine set_bcfilenames(wp)

    implicit none
    ! input/output
    type(waveparamsnew),intent(inout)            :: wp
    ! internal
    integer                                      :: i1,i2,i3,i4,i5

    if (reuseall) then
       wp%Efilename = 'E_reuse.bcf'
       wp%qfilename = 'q_reuse.bcf'
       wp%nhfilename = 'nh_reuse.bcf'
    else
       i1=floor(real(bccount)/10000)
       i2=floor(real(bccount-i1*10000)/1000)
       i3=floor(real(bccount-i1*10000-i2*1000)/100)
       i4=floor(real(bccount-i1*10000-i2*1000-i3*100)/10)
       i5=bccount-i1*10000-i2*1000-i3*100-i4*10
       wp%Efilename='E_series'//char(48+i1)//char(48+i2)//char(48+i3)//char(48+i4)//char(48+i5)//'.bcf'
       wp%qfilename='q_series'//char(48+i1)//char(48+i2)//char(48+i3)//char(48+i4)//char(48+i5)//'.bcf'
       wp%nhfilename='nh_series'//char(48+i1)//char(48+i2)//char(48+i3)//char(48+i4)//char(48+i5)//'.bcf'
    endif

  end subroutine set_bcfilenames

  ! --------------------------------------------------------------
  ! ----------- Merge all separate spectra into one --------------
  ! -------------- average spectrum for other use ----------------
  subroutine generate_combined_spectrum(specinterp,combspec,px)

    implicit none
    type(spectrum),dimension(nspectra),intent(in)   :: specinterp
    type(spectrum),intent(inout)                    :: combspec
    real*8, intent(in)                              :: px
    integer                                         :: iloc
    real*8,dimension(naint)                         :: Sd
    real*8,dimension(3)                             :: peakSd,peakang
    real*8,dimension(:),allocatable                 :: tempSf

    allocate(combspec%f(nfint))
    allocate(combspec%Sf(nfint))
    allocate(combspec%ang(naint))
    allocate(combspec%S(nfint,naint))
    combspec%f=specinterp(1)%f
    combspec%nf=nfint
    combspec%df=specinterp(1)%df
    combspec%ang=specinterp(1)%ang
    combspec%nang=naint
    combspec%dang=specinterp(1)%dang
    combspec%trep = 0.d0
    combspec%S=0.d0
    combspec%Sf=0.d0

    do iloc = 1,nspectra
       combspec%trep = combspec%trep+specinterp(iloc)%trep/nspectra
       combspec%S    = combspec%S   +specinterp(iloc)%S   /nspectra
       combspec%Sf   = combspec%Sf  +specinterp(iloc)%Sf  /nspectra
    enddo
    !
    ! Calculate peak wave angle
    !
    ! frequency integrated variance array
    Sd = sum(combspec%S,DIM=1)
    ! peak location of f-int array
    iloc = maxval(maxloc(Sd))
    ! pick two neighbouring directional bins, including effect of closing circle at
    ! 0 and 2pi rad
    if (iloc>1 .and. iloc<naint) then
       peakSd = (/Sd(iloc-1),Sd(iloc),Sd(iloc+1)/)
       peakang = (/combspec%ang(iloc-1),combspec%ang(iloc),combspec%ang(iloc+1)/)
    elseif (iloc==1) then
       peakSd = (/Sd(naint),Sd(1),Sd(2)/)
       peakang = (/combspec%ang(naint)-2*px,combspec%ang(1),combspec%ang(2)/)
    elseif (iloc==naint) then
       peakSd = (/Sd(naint-1),Sd(naint),Sd(1)/)
       peakang = (/combspec%ang(naint-1),combspec%ang(naint),combspec%ang(1)+2*px/)
    endif
    ! dir0 calculated as mean over peak and two neighbouring cells
    combspec%dir0 = sum(peakSd*peakang)/sum(peakSd)
    ! return to 0<=dir0<=2*px
    if (combspec%dir0>2*px) then
       combspec%dir0 = combspec%dir0-2*px
    elseif (combspec%dir0<0) then
       combspec%dir0 = combspec%dir0+2*px
    endif
    !
    ! Compute peak wave frequency
    !
    allocate(tempSf(size(combspec%Sf)))
    tempSf=0*combspec%Sf
    ! find frequency range of 95% wave energy around peak
    where (combspec%Sf>0.95d0*maxval(combspec%Sf))
       tempSf=combspec%Sf
    end where
    ! smoothed peak is weighted mean frequency of 95% range
    combspec%fp = sum(tempSf*combspec%f)/sum(tempSf)
  end subroutine generate_combined_spectrum


  ! --------------------------------------------------------------
  ! -------------- Calculate average water depth -----------------
  ! ------------------ on offshore boundary ----------------------
  pure function calculate_average_water_depth(s) result(h)

    use spaceparams

    implicit none

    ! input
    type(spacepars),intent(in)     :: s
    ! output
    real*8                         :: h
    ! internal

    h = sum((s%zs0(1,:)-s%zb(1,:))*s%dnc(1,:))/sum(s%dnc(1,:))

  end function calculate_average_water_depth


  ! --------------------------------------------------------------
  ! ----------- Choose wave train components based on ------------
  ! --------------------- combined spectrum ----------------------
  subroutine generate_wavetrain_components(combspec,wp,par)

    use logging_module
    use params
    use interp
    use wave_functions_module, only: iteratedispersion

    implicit none
    ! input/output
    type(spectrum),intent(in)                    :: combspec
    type(waveparamsnew),intent(inout)            :: wp
    type(parameters),intent(in)                  :: par
    ! internal
    integer                                      :: i,ii
    integer                                      :: ind1,ind2,dummy
    real*8,dimension(:),allocatable              :: randnums,pdflocal,cdflocal
    real*8                                       :: L0,L,kmax,fmax


    ! If we are running non-hydrostatic boundary conditions, we want to remove
    ! very high frequency components, as these will not be computed well anyway
    if (par%nonhspectrum==1) then
       kmax = wdmax/wp%h0
       fmax = par%g*kmax*tanh(kmax*wp%h0)/2/par%px
    else
       ! this is really already taken into account
       ! by specifying par%sprdthr
       fmax = 2*maxval(combspec%f)
    endif

    ! Determine frequencies around peak frequency of one-dimensional
    ! non-directional variance density spectrum, based on factor sprdthr, which
    ! should be included in the determination of the wave boundary conditions
    call frange(par,combspec%f,combspec%Sf,fmax,ind1,ind2)

    ! Calculate number of wave components to be included in determination of the
    ! wave boundary conditions based on the wave record length and width of the
    ! wave frequency range
    wp%K = ceiling(wp%rtbc*(combspec%f(ind2)-combspec%f(ind1))+1)
    ! also include minimum number of components
    wp%K = max(wp%K,Kmin)

    ! Allocate space in waveparams for all wave train components
    allocate(wp%fgen(wp%K))
    allocate(wp%thetagen(wp%K))
    allocate(wp%phigen(wp%K))
    allocate(wp%kgen(wp%K))
    allocate(wp%wgen(wp%K))

    ! Select equidistant wave components between the earlier selected range of
    ! frequencies around the peak frequency in the frange subfunction
    wp%dfgen = (combspec%f(ind2)-combspec%f(ind1))/(wp%K-1)
    do i=1,wp%K
       wp%fgen(i)=combspec%f(ind1)+(i-1)*wp%dfgen
    enddo

    ! This subroutine needs to generate random phase / directions for the individual
    ! wave trains. Due to some strange Fortran properties, it is better to select
    ! one long vector with random numbers for both the phase and the direction, than
    ! one vector for each.
    ! Update random seed, if requested
    if (par%random==1) CALL init_seed
    allocate(randnums(2*wp%K))
    call random_number(randnums)

    ! Determine a random phase for each wave train component between 0 and 2pi
    wp%phigen=randnums(1:wp%K)*2*par%px

    ! Determine random directions for each wave train component, based on the CDF of
    ! the directional spectrum. For each wave train we will interpolate the directional
    ! distribution at that frequency, generate a CDF, and then interpolate the wave
    ! train direction from a random number draw and the CDF.
    allocate(pdflocal(naint))
    allocate(cdflocal(naint))
    do i=1,wp%K
       ! interpolate spectrum at this frequency per angle in the spectrum
       do ii=1,naint
          call LINEAR_INTERP(combspec%f,combspec%S(:,ii),combspec%nf, wp%fgen(i),pdflocal(ii),dummy)
       enddo
       ! convert to pdf by ensuring total integral == 1, assuming constant directional bin size
       pdflocal = pdflocal/sum(pdflocal)
       ! convert to cdf by trapezoidal integration, which is not biased for integration
       ! direction.
       ! Boundary condition, which may be nonzero:
       cdflocal(1) = pdflocal(1)
       do ii=2,naint
          cdflocal(ii) = cdflocal(ii-1) + (pdflocal(ii)+pdflocal(ii-1))/2
          ! Note: this only works if the directional
          ! bins are constant in size. Assumed multiplication
          ! by one.
       enddo
       ! interpolate random number draw across the cdf. Note that cdf(1) is not assumed to be zero and
       ! there is a posibility of drawing less than cdf(1). Therefore, if the random number is .lt. cdf(1),
       ! interpolation should take place across the back of the angle spectrum
       ! (i.e., from theta(end)-2pi : theta(1)).
       ! Note that we are using the second half of randnums now (K+1:end)
       if (randnums(wp%K+i)>=cdflocal(1)) then
          call LINEAR_INTERP(cdflocal,combspec%ang,naint,randnums(wp%K+i), wp%thetagen(i),dummy)
       else
          call LINEAR_INTERP((/0.d0,cdflocal(1)/),(/combspec%ang(naint)-2*par%px,combspec%ang(1)/), &
               2,randnums(wp%K+i),wp%thetagen(i),dummy)
          if (wp%thetagen(i)<0.d0) then
             wp%thetagen(i)=wp%thetagen(i)+2*par%px
          endif
       endif
    enddo

    ! determine wave number for each wave train component, using standard dispersion relation
    ! solver from wave_functions module. This function returns a negative wave length if the
    ! solver did not converge.
    do i=1,wp%K
       L0 = par%g*(1/wp%fgen(i))**2/2/par%px    ! deep water wave length
       L = iteratedispersion(L0,L0,par%px,wp%h0)
       if (L<0.d0) then
          call writelog('lsw','','No dispersion convergence found for wave train ',i, &
               ' in boundary condition generation')
          L = -L
       endif
       wp%kgen(i) = 2*par%px/L
    enddo

    ! Angular frequency
    wp%wgen = 2*par%px*wp%fgen
       
    ! Free memory
    deallocate(pdflocal,cdflocal,randnums)

  end subroutine generate_wavetrain_components

  ! -----------------------------------------------------------
  ! ----------- Small subroutine to determine -----------------
  ! ------------ representative wave period -------------------
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
       Trep = sum(temp*Sf/max(f,0.001d0))/sum(temp*Sf)    ! Tm-1,0
    endif

    deallocate(temp)

  end subroutine tpDcalc


  ! -----------------------------------------------------------
  ! ---- Small subroutine to determine f-range round peak -----
  ! -----------------------------------------------------------
  subroutine frange(par,f,Sf,fmax,firstp,lastp,findlineout)

    use params

    implicit none

    type(parameters),intent(in)             :: par
    real*8, dimension(:), intent(in)        :: f,Sf
    real*8,intent(in)                       :: fmax

    integer, intent(out)                    :: firstp, lastp

    real*8, dimension(:), intent(out),allocatable,optional  :: findlineout
    real*8, dimension(:),allocatable        :: temp, findline
    integer                                 :: i = 0

    ! find frequency range around peak
    allocate(findline(size(Sf)))
    findline=0*Sf
    where (Sf>par%sprdthr*maxval(Sf))
       findline=1
    end where

    ! find frequency range below fmax
    where(f>fmax)
       findline = 0
    endwhere

    firstp=maxval(maxloc(findline))         ! Picks the first "1" in temp

    allocate (temp(size(findline)))
    temp=(/(i,i=1,size(findline))/)
    lastp=maxval(maxloc(temp*findline))     ! Picks the last "1" in temp

    if (present(findlineout)) then
       allocate(findlineout(size(Sf)))
       findlineout=findline
    endif
    deallocate(temp, findline)

  end subroutine frange

  ! -----------------------------------------------------------
  ! --- Small subroutine to reseed random number generator ----
  ! -----------------------------------------------------------
  subroutine init_seed
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    ! RANDOM_SEED does not result in a random seed for each compiler, so we better make sure that we have a pseudo random seed.
    ! I am not sure what is a good n
    n = 40
    i = 1
    ! not sure what size means here
    ! first call with size
    CALL RANDOM_SEED(size = n)
    ! define a seed array of size n
    ALLOCATE(seed(n))
    ! what time is it
    CALL SYSTEM_CLOCK(COUNT=clock)
    ! define the seed vector based on a prime, the clock and the set of integers
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    ! if mpi do we need a different seed on each node or the same???
    ! if we do need different seeds on each node
    ! seed *= some big prime * rank ?
    ! now use the seed
    CALL RANDOM_SEED(PUT = seed)
  end subroutine init_seed


  ! --------------------------------------------------------------
  ! ---------------- Setup time axis for waves -------------------
  ! --------------------------------------------------------------
  subroutine generate_wave_time_axis(par,wp)

    use params

    implicit none
    ! input/output
    type(waveparamsnew),intent(inout)            :: wp
    type(parameters),intent(in)                  :: par
    ! internal
    integer                                      :: ntaper, indend
    integer                                      :: i



    ! First assume that internal and bc-writing time step is the same
    wp%dtin = wp%dtbc
    wp%dtchanged = .false.
    wp%tslenbc = nint(wp%rtbc/wp%dtbc)+1

    ! Check whether the internal frequency is high enough to describe the highest frequency
    ! wave train returned from frange (which can be used in the boundary conditions)
    if (par%nonhspectrum==0) then
       if (wp%dtin>0.5d0/wp%fgen(wp%K)) then
          wp%dtin = 0.5d0/wp%fgen(wp%K)
          wp%dtchanged = .true.
       endif
    else
       if (wp%dtin>0.1d0/wp%fgen(wp%K)) then
          wp%dtin = 0.1d0/wp%fgen(wp%K)
          wp%dtchanged = .true.
       endif
    endif

    ! The length of the internal time axis should be even (for Fourier transform) and
    ! depends on the internal time step needed and the internal duration (~1/dfgen):
    wp%tslen = ceiling(1/wp%dfgen/wp%dtin)+1
    if (mod(wp%tslen,2)/=0) then
       wp%tslen = wp%tslen +1
    end if

    ! Now we can make the internal time axis
    wp%rtin = wp%tslen * wp%dtin
    allocate(wp%tin(wp%tslen))
    do i=1,wp%tslen
       wp%tin(i) = (i-1)*wp%dtin
    enddo

    ! Make a taper function to slowly increase and decrease the boundary condition forcing
    ! at the start and the end of the boundary condition file (including any time beyond
    ! the external rtbc
    allocate(wp%taperf(wp%tslen))
    allocate(wp%taperw(wp%tslen))
    ! fill majority with unity
    wp%taperf = 1.d0
    wp%taperw = 1.d0
    if (par%nonhspectrum==1) then
       ! begin taper by building up the wave conditions over 2 wave periods
       ntaper = nint((2.0d0*par%Trep)/wp%dtin)
    else
       ntaper = nint((5.d0*par%Trep)/wp%dtin)
    endif
    do i=1,min(ntaper,size(wp%taperf))
       wp%taperf(i) = tanh(5.d0*i/ntaper)     ! multiplied by five because tanh(5)=~1
       wp%taperw(i) = tanh(5.d0*i/ntaper)
    enddo
    ! We do not want to taperw the end anymore. Instead we pass the wave height at the end of rtbc to
    ! the next wave generation iteration.
    ! end taper by finding where tin=rtbc, taper before that and set everything to zero after
    ! that.
    if (wp%tin(wp%tslen)>wp%rtbc) then
       indend = minval(minloc(wp%tin,MASK = wp%tin>=wp%rtbc))
    else
       indend = wp%tslen
    endif
    do i=1,min(ntaper,indend)
       wp%taperf(indend+1-i) = min(wp%taperf(indend+1-i),tanh(5.d0*i/ntaper))
    enddo
    wp%taperf(indend:wp%tslen) = 0.d0
    ind_end_taper = indend
  end subroutine generate_wave_time_axis


  ! --------------------------------------------------------------
  ! -------- Calculate variance at each spectrum location --------
  ! --------------------------------------------------------------
  subroutine generate_wave_train_variance(wp,specinterp,par)

    use params
    use interp

    implicit none
    ! input/output
    type(waveparamsnew),intent(inout)            :: wp
    type(spectrum),dimension(nspectra),intent(in):: specinterp
    type(parameters),intent(in)                  :: par
    ! internal
    integer                                      :: i,ii, dummy
    real*8                                       :: hm0post

    ! allocate space for the variance arrays
    allocate(wp%vargen(nspectra))
    allocate(wp%vargenq(nspectra))

    ! Determine variance at each spectrum location
    do i=1,nspectra
       allocate(wp%vargen(i)%Sf(wp%K))
       allocate(wp%vargenq(i)%Sf(wp%K))
       do ii=1,wp%K
          ! In order to maintain energy density per frequency, an interpolation
          ! is carried out on the input directional variance density spectrum
          ! at the frequency and direction locations in fgen. This must be done
          ! over the 2D spectrum, not 1D spectrum.
!          call linear_interp_2d(specinterp(i)%f,nfint, &   ! input frequency, size
!               specinterp(i)%ang,naint, &                  ! input angles, size
!               specinterp(i)%S*specinterp(i)%dang, &       ! input variance density times grid angle resolution
!               wp%fgen(ii),wp%thetagen(ii), &              ! output frquency/angle
!               wp%vargen(i)%Sf(ii), &                      ! output variance density
!               'interp',0.d0)                              ! method and exception return value
          call linear_interp(specinterp(i)%f, &
                             specinterp(i)%Sf, &
                             nfint, &
                             wp%fgen(ii), &
                             wp%vargen(i)%Sf(ii),dummy)
       enddo
       ! Correct variance to ensure the total variance remains the same as the total variance in the
       ! (interpolated) input spectrum
       hm0post = 4*sqrt(sum(wp%vargen(i)%Sf)*wp%dfgen)
       wp%vargen(i)%Sf   = (specinterp(i)%hm0/hm0post)**2*wp%vargen(i)%Sf
       ! For the generation of long waves we cannot use wp%vargen%Sf, because it contains an overestimation
       ! of energy in the peak frequencies. We can also not use the standard directionally-integrated spectrum
       ! because this stores all energy at fgen(K), where it is possible that for this current spectrum, there
       ! is no energy at S(f(K),theta(K0))
       ! The current solution is to take the minimum of both methods
       do ii=1,wp%K
          ! Map Sf input to fgen
          call LINEAR_INTERP(specinterp(i)%f,specinterp(i)%Sf,specinterp(i)%nf,wp%fgen(ii),wp%vargenq(i)%Sf(ii),dummy)
       enddo
       wp%vargenq(i)%Sf = min(wp%vargen(i)%Sf,wp%vargenq(i)%Sf)
    enddo
 
  end subroutine generate_wave_train_variance

  ! --------------------------------------------------------------
  ! ------ Calculate amplitudes at each spectrum location --------
  ! ------------ for every wave train component ------------------
  subroutine generate_wave_train_properties_per_offshore_point(wp,s)

    use spaceparams
    use interp

    implicit none
    ! input/output
    type(waveparamsnew),intent(inout)            :: wp
    type(spacepars),intent(in)                   :: s
    ! internal
    integer                                      :: i,ii,dummy
    integer,dimension(2)                         :: interpindex
    real*8,dimension(2)                          :: positions,sf
    logical                                      :: interpcase
    real*8                                       :: here,sfnow,sfimp
    integer,dimension(size(n_index_loc))         :: temp_index_loc

    ! allocate space for the amplitude array and representative integration angle
    allocate(wp%A(s%ny+1,wp%K))
    !    allocate(wp%danggen(s%ny+1))
    allocate(wp%Sfinterp(s%ny+1,wp%K))
    allocate(wp%Sfinterpq(s%ny+1,wp%K))
    allocate(wp%Hm0interp(s%ny+1))

    ! where necessary, interpolate Sf of each spectrum location
    ! to the current grid cell, and use this to calculate A
    do i=1,s%ny+1
       ! Four possibilities:
       ! 1: this exact point has an input spectrum
       ! 2: this point lies between two input spectra
       ! 3: this point lies before any input spectrum
       ! 4: this point lies beyond any input spectrum
       if (any(n_index_loc==i)) then  ! case 1
          interpcase = .false.
          do ii=1,nspectra
             if(n_index_loc(ii)==i) then
                interpindex = ii
             endif
          enddo
       elseif(i<minval(n_index_loc)) then   ! case 3
          interpcase = .false.
          do ii=1,nspectra
             if(n_index_loc(ii)==minval(n_index_loc)) then
                interpindex = ii
             endif
          enddo
       elseif(i>maxval(n_index_loc)) then   ! case 4
          interpcase = .false.
          do ii=1,nspectra
             if(n_index_loc(ii)==maxval(n_index_loc)) then
                interpindex = ii
             endif
          enddo
       else   ! case 2
          interpcase = .true.
          temp_index_loc = n_index_loc
          where(n_index_loc>i)
             temp_index_loc = -huge(0)
          endwhere
          interpindex(1) = minval(maxloc(temp_index_loc))

          temp_index_loc = n_index_loc
          where(n_index_loc<i)
             temp_index_loc = huge(0)
          endwhere
          interpindex(2) = minval(minloc(temp_index_loc))
       endif

       ! interpolate Sf
       if (interpcase) then
          positions(1) = s%ndist(1,n_index_loc(interpindex(1)))
          positions(2) = s%ndist(1,n_index_loc(interpindex(2)))
          here = s%ndist(1,i)
          do ii=1,wp%K
             sf(1) = wp%vargen(interpindex(1))%Sf(ii)
             sf(2) = wp%vargen(interpindex(2))%Sf(ii)
             call LINEAR_INTERP(positions,sf,2,here,wp%Sfinterp(i,ii),dummy)
             sf(1) = wp%vargenq(interpindex(1))%Sf(ii)
             sf(2) = wp%vargenq(interpindex(2))%Sf(ii)
             call LINEAR_INTERP(positions,sf,2,here,wp%Sfinterpq(i,ii),dummy)
          enddo
          ! Ensure that the total amount of variance has not been reduced/increased
          ! during interpolation for short wave energy only. This is not done for
          ! Sf for bound long wave generation
          sfnow = sum(wp%Sfinterp(i,:))  ! integration over fixed bin df size unnecessary
          sf(1) = sum(wp%vargen(interpindex(1))%Sf)
          sf(2) = sum(wp%vargen(interpindex(2))%Sf)
          call LINEAR_INTERP(positions,sf,2,here,sfimp,dummy)
          wp%Sfinterp(i,:)=wp%Sfinterp(i,:)*sfimp/sfnow
       else
          wp%Sfinterp(i,:) = wp%vargen(interpindex(1))%Sf
          wp%Sfinterpq(i,:) = wp%vargenq(interpindex(1))%Sf
       endif

       wp%A(i,:) = sqrt(2*wp%Sfinterp(i,:)*wp%dfgen)

       wp%Hm0interp(i) = 4*sqrt(sum(wp%Sfinterp(i,:))*wp%dfgen)

    enddo ! i=1,s%ny+1

  end subroutine generate_wave_train_properties_per_offshore_point

  ! --------------------------------------------------------------
  ! --------- Calculate Fourier componets for each wave ----------
  ! ---------- train component on the offshore boundary ----------
  subroutine generate_wave_train_Fourier(wp,s,par)

    use spaceparams
    use interp
    use math_tools
    use params
    use logging_module

    implicit none
    ! input/output
    type(waveparamsnew),intent(inout)            :: wp
    type(spacepars),intent(in)                   :: s
    type(parameters),intent(in)                  :: par
    ! internal
    integer                                      :: i,ii,tempi
    complex(fftkind),dimension(:),allocatable    :: tempcmplx

    ! Determine indices of wave train components in frequency axis and
    ! Fourier transform result
    allocate(wp%Findex(wp%K))
    tempi = floor(wp%fgen(1)/wp%dfgen)
    do i=1,wp%K
       wp%Findex(i) = tempi + i
    enddo

    ! Allocate Fourier coefficients for each y-position at the offshore boundary and
    ! each time step
    allocate(wp%CompFn(s%ny+1,wp%tslen))
    wp%CompFn=0.d0
    allocate(tempcmplx(wp%tslen/2-1))
    call writelog('ls','','Calculating Fourier components')
    call progress_indicator(.true.,0.d0,5.d0,2.d0)
    do i=1,wp%K
       call progress_indicator(.false.,dble(i)/wp%K*100,5.d0,2.d0)
       do ii=1,s%ny+1
          ! Determine first half of complex Fourier coefficients of wave train
          ! components using random phase and amplitudes from sampled spectrum
          ! until Nyquist frequency. The amplitudes are represented in a
          ! two-sided spectrum, which results in the factor 1/2.
          ! Unroll Fourier components along offshore boundary, assuming all wave trains
          ! start at x(1,1), y(1,1).
          wp%CompFn(ii,wp%Findex(i)) = wp%A(ii,i)/2*exp(par%compi*wp%phigen(i))* &    ! Bas: wp%Findex used in time dimension because t=j*dt in frequency space
               exp(-par%compi*wp%kgen(i)*(dsin(wp%thetagen(i))*(s%yz(1,ii)-s%yz(1,1)) &
               +dcos(wp%thetagen(i))*(s%xz(1,ii)-s%xz(1,1))) )

          ! Determine Fourier coefficients beyond Nyquist frequency (equal to
          ! coefficients at negative frequency axis) of relevant wave components for
          ! first y-coordinate by mirroring
          tempcmplx = conjg(wp%CompFn(ii,2:wp%tslen/2))
          call flipiv(tempcmplx,size(tempcmplx))
          wp%CompFn(ii,wp%tslen/2+2:wp%tslen)=tempcmplx
       enddo
    enddo

    ! Free memory
    deallocate(tempcmplx)

  end subroutine generate_wave_train_Fourier

  ! --------------------------------------------------------------
  ! --------- Calculate in which computational wave bin ----------
  ! ------------- each wave train component belongs --------------
  subroutine distribute_wave_train_directions(wp,s,px,nspr,nonhspec)

    use spaceparams
    use logging_module

    implicit none

    ! input/output
    type(waveparamsnew),intent(inout)            :: wp
    type(spacepars),intent(in)                   :: s
    real*8,intent(in)                            :: px
    integer,intent(in)                           :: nspr
    logical,intent(in)                           :: nonhspec
    ! internal
    integer                                      :: i,ii
    real*8,dimension(:),allocatable              :: binedges
    logical                                      :: toosmall,toolarge
    real*8                                       :: lostvar,keptvar,perclost
    integer                                      :: lntheta              ! local copies that can be changed without
    real*8                                       :: lthetamin,ldtheta    ! damage to the rest of the model

    ! Set basic parameters for comparison if using nonh spectrum
    if (nonhspec) then
       lntheta = 1
       lthetamin = -0.5d0*px
       ldtheta = px
    else
       lntheta = s%ntheta
       lthetamin = s%thetamin
       ldtheta = s%dtheta
    endif

    ! Calculate the bin edges of all the computational wave bins in the
    ! XBeach model (not the input spectrum)
    allocate(binedges(lntheta+1))
    do i=1,lntheta+1
       binedges(i) = lthetamin+(i-1)*ldtheta
    enddo

    ! Try to fit as many wave train components as possible between the
    ! minimum and maximum bin edges by adding/subtracting 2pi rad. This
    ! solves the problem of XBeach not including energy in the computation
    ! than should be in there. For instance XBeach does not realize that
    ! 90 degrees == -270 degrees == 450 degrees.
    ! Note: this does not ensure all energy is included in the wave bins,
    ! as wave energy may still fall outside the computational domain.
    do i=1,wp%K
       if (wp%thetagen(i)>maxval(binedges)) then
          toosmall = .false.
          toolarge = .true.
       elseif (wp%thetagen(i)<minval(binedges)) then
          toosmall = .true.
          toolarge = .false.
       else
          toosmall = .false.
          toolarge = .false.
       endif
       ! Continue to reduce/increase wave direction until
       ! condition is not met.
       do while (toosmall)
          wp%thetagen(i) = wp%thetagen(i) + 2*px
          toosmall = wp%thetagen(i) < minval(binedges)
       enddo
       do while (toolarge)
          wp%thetagen(i) = wp%thetagen(i) - 2*px
          toolarge = wp%thetagen(i) > maxval(binedges)
       enddo
    enddo

    ! distribute the wave train directions among the wave directional bins,
    ! store the corresponding bin per wave train in WDindex.
    allocate(wp%WDindex(wp%K))
    wp%WDindex = 0
    do i=1,wp%K
       if (wp%thetagen(i) < minval(binedges)) then
          ! remove from computation if wave direction still not within computational
          ! bounds
          wp%WDindex(i) = 0
       elseif (wp%thetagen(i) > maxval(binedges)) then
          wp%WDindex(i) = lntheta+1
       elseif (wp%thetagen(i) == binedges(lntheta+1) ) then
          ! catch upper boundary, and add to upper bin
          wp%WDindex(i) = lntheta
       else
          ! fit this wave train component into one of the lower bins
          do ii=1,lntheta
             if (wp%thetagen(i)>=binedges(ii) .and. wp%thetagen(i)<binedges(ii+1)) then
                wp%WDindex(i) = ii
             endif
          enddo
       endif
    enddo

    ! If the user has set nspr == 1 then the randomly drawn wave directions
    ! should be set to the centres of the wave directional bins.
    ! Also move all wave energy falling outside the computational bins, into
    ! the computational domain (in the outer wave direction bins)
    if (nspr==1) then
       if (nonhspec) then
          ! We only really want to include all energy between -90 and +90 degrees.
          ! Waves outside this direction are sent out of the model
          where(wp%WDindex==1)
             wp%thetagen = 0.d0
          elsewhere
             wp%thetagen = -px
          endwhere
       else
          do i=1,wp%K
             if (wp%WDindex(i)==0) then
                wp%WDindex(i)=1
             elseif (wp%WDindex(i)>lntheta) then
                wp%WDindex(i)=lntheta
             endif
             ! reset the direction of this wave train to the centre of the bin
             wp%thetagen(i)=s%theta(wp%WDindex(i))
          enddo
       endif
    endif

    ! Check the amount of energy lost to wave trains falling outside the computational
    ! domain
    lostvar = 0.d0
    keptvar = 0.d0
    do i=1,wp%K
       if (wp%WDindex(i)==0) then
          lostvar = lostvar + sum(wp%A(:,i)**2)
       else
          keptvar = keptvar + sum(wp%A(:,i)**2)
       endif
    enddo
    perclost = 100*(lostvar/(lostvar+keptvar))
    if (perclost>5.0d0) then
       call writelog('lsw','(a,f0.1,a)','Large amounts of energy (',perclost, &
            '%) fall outside computational domain at the offshore boundary')
       call writelog('lsw','','Check specification of input wave angles and wave directional grid')
    else
       call writelog('ls','(a,f0.1,a)','Wave energy outside computational domain at offshore boundary: ',perclost,'%')
    endif

    ! Free memory
    deallocate(binedges)

  end subroutine distribute_wave_train_directions

  ! --------------------------------------------------------------
  ! --------- Calculate energy envelope time series from ---------
  ! -------- Fourier components, and write to output file --------
  subroutine generate_ebcf(wp,s,par)

    use params
    use spaceparams
    use math_tools
    use interp
    use logging_module
    use filefunctions, only: create_new_fid

    implicit none

    ! input/output
    type(waveparamsnew),intent(inout)            :: wp
    type(spacepars),intent(in)                   :: s
    type(parameters),intent(in)                  :: par
    ! internal
    integer                                      :: itheta,iy,iwc,it,irec
    integer                                      :: index,status
    integer                                      :: reclen,fid
    integer,dimension(:),allocatable             :: nwc
    real*8,dimension(:,:,:), allocatable         :: zeta, Ampzeta, E_tdir, E_interp
    real*8,dimension(:,:), allocatable           :: eta, Amp
    complex(fftkind),dimension(:),allocatable    :: Gn, tempcmplx,tempcmplxhalf
    integer,dimension(:),allocatable             :: tempindex,tempinclude
    real*8                                       :: stdzeta,stdeta,etot,perc


    ! Allocate variables for water level exitation and amplitude with and without
    ! directional spreading dependent envelope
    allocate(zeta(s%ny+1,wp%tslen,s%ntheta))
    allocate(Ampzeta(s%ny+1,wp%tslen,s%ntheta))
    zeta=0.d0
    Ampzeta=0.d0

    allocate(eta(s%ny+1,wp%tslen))
    allocate(Amp(s%ny+1,wp%tslen))
    eta=0.d0
    Amp=0.d0

    ! Calculate wave energy for each y-coordinate along seaside boundary for
    ! current computational directional bin
    allocate(nwc(s%ntheta))
    allocate(Gn(wp%tslen))
    allocate(tempinclude(wp%K))
    allocate(tempcmplx(wp%tslen))
    allocate(tempcmplxhalf(size(Gn(wp%tslen/2+2:wp%tslen))))
    do itheta=1,s%ntheta
       call writelog('ls','(A,I0,A,I0)','Calculating short wave time series for theta bin ',itheta,' of ',s%ntheta)
       ! Select wave components that are in the current computational
       ! directional bin
       tempinclude=0
       where (wp%WDindex==itheta)
          tempinclude=1
       end where

       ! Determine number of wave components that are in the current
       ! computational directional bin
       nwc(itheta)=sum(tempinclude)

       ! Determine for each wave component in the current computational
       ! directional bin its index in the Fourier coefficients array
       ! ordered from high to low frequency
       allocate(tempindex(nwc(itheta)))
       tempindex=0
       do iwc=1,nwc(itheta)
          ! find highest
          index=maxval(maxloc(tempinclude))
          ! reset that one so that the next highest in found in next iteration
          tempinclude(index)=0
          tempindex(iwc)=wp%Findex(index)
       end do

       ! Check whether any wave components are in the current computational
       ! directional bin
       if (nwc(itheta)>0) then
          do iy=1,s%ny+1
             ! Reset
             Gn=0

             ! Determine Fourier coefficients of all wave components for current
             ! y-coordinate in the current computational directional bin
             Gn(tempindex)=wp%CompFn(iy,tempindex)
             tempcmplxhalf = conjg(Gn(2:wp%tslen/2))
             call flipiv(tempcmplxhalf,size(tempcmplxhalf))
             Gn(wp%tslen/2+2:wp%tslen)=tempcmplxhalf

             ! Inverse Discrete Fourier transformation to transform back to time
             ! domain from frequency domain
             tempcmplx=Gn
             status=0
             tempcmplx=fft(tempcmplx,inv=.true.,stat=status)

             ! Scale result
             tempcmplx=tempcmplx/sqrt(dble(size(tempcmplx)))

             ! Superimpose gradual increase and decrease of energy input for
             ! current y-coordinate and computational diretional bin on
             ! instantaneous water level excitation
             ! zeta(iy,:,itheta)=dble(tempcmplx*wp%tslen)*wp%taperw
             ! Robert: use final wave elevation from last iteration to startup
             ! this boundary condition
             zeta(iy,:,itheta)=dble(tempcmplx*wp%tslen)
             zeta(iy,:,itheta)=zeta(iy,:,itheta)*wp%taperw+(1-wp%taperw)*lastwaveelevation(iy,itheta)
             ! note that taperw is only zero at the start, not the end of the time series, so we
             ! can safely copy zeta at the point in time where t=rtbc to lastwaveelevation for use
             ! in the generation of the next time series
             lastwaveelevation(:,itheta) = zeta(iy,ind_end_taper,itheta)
          enddo ! iy=1,ny+1
       endif  ! nwc>0
       deallocate(tempindex)
    enddo ! itheta = 1,ntheta
    deallocate(tempinclude)
    deallocate(Gn)
    deallocate(tempcmplxhalf)
    !
    ! Calculate energy envelope amplitude

    ! Note Robert: oldwbc should be deprecated for wbcversion>2
    if (par%oldwbc==1) then
       do itheta=1,s%ntheta
          ! Print message to screen
          call writelog('ls','(A,I0,A,I0)','Calculating wave energy for theta bin ',itheta,' of ',s%ntheta)
          if (nwc(itheta)>0) then
             do iy=1,s%ny+1
                tempcmplx=zeta(iy,:,itheta)
                ! Hilbert tranformation to determine envelope for each directional
                ! bin seperately
                call hilbert(tempcmplx,size(tempcmplx))
                ! Determine amplitude of water level envelope by calculating
                ! the absolute value of the complex wave envelope descriptions
                Ampzeta(iy,:,itheta)=abs(tempcmplx)
                ! Print status message to screen
                call writelog('ls','(A,I0,A,I0,A)','Y-point ',iy,' of ',s%ny+1,' done.')
             end do   ! 1:ny+1
          else  !  nwc==0
             ! Current computational directional bin does not contain any wave
             ! components, so print message to screen
             Ampzeta(iy,:,itheta)=0.d0
          endif  ! nwc>0
       end do   ! 1:ntheta
    else
       do iy=1,s%ny+1
          ! Integrate instantaneous water level excitation of wave
          ! components over directions
          eta(iy,:) = sum(zeta(iy,:,:),2)
          tempcmplx=eta(iy,:)

          ! Hilbert transformation to determine envelope of all total
          ! non-directional wave components
          call hilbert(tempcmplx,size(tempcmplx))

          ! Determine amplitude of water level envelope by calculating
          ! the absolute value of the complex wave envelope descriptions
          Amp(iy,:)=abs(tempcmplx)

          ! Calculate standard deviation of non-directional
          ! instantaneous water level excitation of all
          ! wave components to be used as weighing factor
          stdeta = sqrt(sum(eta(iy,:)**2)/(size(eta(iy,:))-1))
          do itheta=1,s%ntheta
             if (nwc(itheta)>0) then
                ! Calculate standard deviations of directional
                ! instantaneous water level excitation of all
                ! wave components to be used as weighing factor
                stdzeta = sqrt(sum(zeta(iy,:,itheta)**2)/(size(zeta(iy,:,itheta))-1))

                ! Calculate amplitude of directional wave envelope
                Ampzeta(iy,:,itheta)= Amp(iy,:)*stdzeta/stdeta
             else  !  nwc==0
                ! Current computational directional bin does not contain any wave
                ! components, so print message to screen
                Ampzeta(iy,:,itheta)=0.d0
             endif  ! nwc>0
          end do   ! 1:ntheta
          ! Print status message to screen
          call writelog('ls','(A,I0,A,I0,A)','Y-point ',iy,' of ',s%ny+1,' done.')
       end do   ! 1:ny+1
    endif ! oldwbc
    !

    ! free memory
    deallocate(tempcmplx)
    deallocate(nwc)

    ! Allocate memory for energy time series
    allocate(E_tdir(s%ny+1,wp%tslen,s%ntheta))
    E_tdir=0.0d0
    E_tdir=0.5d0*par%rho*par%g*Ampzeta**2
    E_tdir=E_tdir/s%dtheta

    !    ! Ensure we scale back to the correct Hm0
    !    do iy=1,s%ny+1
    !       stdeta = sum(E_tdir(iy,:,:))*s%dtheta   ! sum energy
    !       stdeta = stdeta/wp%tslen       ! mean energy
    !
    !       stdzeta = (wp%Hm0interp(iy)/sqrt(2.d0))**2 * par%rhog8
    !
    !       E_tdir(iy,:,:) = E_tdir(iy,:,:)*stdzeta/stdeta
    !    enddo

    ! Print directional energy distribution to screen
    etot = sum(E_tdir)
    do itheta=1,s%ntheta
       perc = sum(E_tdir(:,:,itheta))/etot*100
       call writelog('ls','(a,i0,a,f0.2,a)','Wave bin ',itheta,' contains ',perc,'% of total energy')
    enddo


    ! Open file for storage
    call writelog('ls','','Writing wave energy to ',trim(wp%Efilename),' ...')
    inquire(iolength=reclen) 1.d0
    reclen=reclen*(s%ny+1)*(s%ntheta)
    fid = create_new_fid()
    open(fid,file=trim(wp%Efilename),form='unformatted',access='direct',recl=reclen,status='REPLACE')

    ! Write to external file
    if (wp%dtchanged) then
       ! Interpolate from internal time axis to output time axis
       allocate(E_interp(s%ny+1,wp%tslenbc,s%ntheta))
       do itheta=1,s%ntheta
          do it=1,wp%tslenbc
             do iy=1,s%ny+1
                call linear_interp(wp%tin,E_tdir(iy,:,itheta),wp%tslen,(it-1)*wp%dtbc,E_interp(iy,it,itheta),status)
             enddo
          enddo
       enddo
       ! write to file
       do irec=1,wp%tslenbc+1
          write(fid,rec=irec)E_interp(:,min(irec,wp%tslenbc),:)
       end do
       deallocate(E_interp)
    else
       ! no need for interpolation
       do irec=1,wp%tslenbc+1
          write(fid,rec=irec)E_tdir(:,min(irec,wp%tslenbc),:)
       end do
    endif
    close(fid)
    call writelog('sl','','file done')

    ! Free memory
    deallocate(zeta,Ampzeta,E_tdir, Amp, eta)

  end subroutine generate_ebcf


  ! --------------------------------------------------------------
  ! ------ Calculate time series of short wave at offshore -------
  ! --------------------------- boundary -------------------------
  subroutine generate_swts(wp,s,par)
    use params
    use spaceparams
    use logging_module
    use filefunctions, only: create_new_fid

    implicit none

    ! input/output
    type(waveparamsnew),intent(inout)            :: wp
    type(spacepars),intent(in)                   :: s
    type(parameters),intent(in)                  :: par
    ! internal
    integer                                      :: j,it,ik,iz

    ! allocate memory for time series of data
    allocate(wp%zsits(s%ny+1,wp%tslen))
    allocate(wp%uits(s%ny+1,wp%tslen))
    wp%zsits=0.d0
    wp%uits=0.d0
    ! total surface elevation
    call writelog('ls','','Calculating short wave elevation time series')
    call progress_indicator(.true.,0.d0,5.d0,2.d0)
    do it=1,wp%tslen
       call progress_indicator(.false.,dble(it)/wp%tslen*100,5.d0,2.d0)
       do ik=1,wp%K
          do j=1,s%ny+1
             wp%zsits(j,it)=wp%zsits(j,it)+wp%A(j,ik)*dsin( &
                  +wp%wgen(ik)*wp%tin(it)&
                  -wp%kgen(ik)*( dsin(wp%thetagen(ik))*(s%yz(1,j)-s%yz(1,1)) &
                  +dcos(wp%thetagen(ik))*(s%xz(1,j)-s%xz(1,1))) &
                  +wp%phigen(ik) &
                  )
          enddo
       enddo
    enddo
    ! depth-averaged velocity
    call writelog('ls','','Calculating short wave velocity time series')
    call progress_indicator(.true.,0.d0,5.d0,2.d0)
    do it=1,wp%tslen
       call progress_indicator(.false.,dble(it)/wp%tslen*100,5.d0,2.d0)
       do ik=1,wp%K
          do j=1,s%ny+1
             wp%uits(j,it) = wp%uits(j,it) + &
                             1.d0/wp%h0*wp%wgen(ik)*wp%A(j,ik)/sinh(wp%kgen(ik)*wp%h0) * &
                                                dsin( &
                                                      +wp%wgen(ik)*wp%tin(it)&
                                                      -wp%kgen(ik)*( dsin(wp%thetagen(ik))*(s%yz(1,j)-s%yz(1,1)) &
                                                                    +dcos(wp%thetagen(ik))*(s%xz(1,j)-s%xz(1,1))) &
                                                      +wp%phigen(ik) &
                                                     ) * &
                                                1.d0/wp%kgen(ik)*sinh(wp%kgen(ik)*wp%h0)

          enddo
       enddo
    enddo
    !
    !
    ! Apply tapering to time series
    do j=1,s%ny+1
       wp%uits(j,:)=wp%uits(j,:)*wp%taperf
       wp%zsits(j,:)=wp%zsits(j,:)*wp%taperf
    enddo
  end subroutine generate_swts


  ! --------------------------------------------------------------
  ! ----------------------- Bound long wave ----------------------
  ! --------------------------------------------------------------
  subroutine generate_qbcf(wp,s,par)

    use params
    use spaceparams
    use logging_module
    use math_tools
    use interp
    use filefunctions, only: create_new_fid

    implicit none

    ! input/output
    type(waveparamsnew),intent(inout)            :: wp
    type(spacepars),intent(in)                   :: s
    type(parameters),intent(in)                  :: par
    ! internal
    integer                                      :: j,m,iq,irec,it  ! counters
    integer                                      :: K               ! copy of K
    integer                                      :: halflen,reclen,fid,status
    logical                                      :: firsttime    ! used for output message only
    real*8                                       :: deltaf       ! difference frequency
    real*8,dimension(:), allocatable             :: term1,term2,term2new,dif,chk1,chk2
    real*8,dimension(:,:),allocatable            :: Eforc,D,deltheta,KKx,KKy,dphi3,k3,cg3,theta3,Abnd
    real*8,dimension(:,:,:),allocatable          :: q,qinterp
    complex(fftkind),dimension(:),allocatable    :: Comptemp,Comptemp2,Gn
    complex(fftkind),dimension(:,:,:),allocatable:: Ftemp



    ! This function has changed with respect to previous versions of XBeach, in that
    ! the bound long wave has to be calculated separately at each longshore point,
    ! owing to longshore varying incident wave spectra

    ! shortcut variable
    K = wp%K

    ! Print message to screen
    call writelog('sl','', 'Calculating primary wave interaction')

    ! Allocate two-dimensional variables for all combinations of interacting wave
    ! components to be filled triangular
    allocate(Eforc(K-1,K))
    allocate(D(K-1,K))
    allocate(deltheta(K-1,K))
    allocate(KKx(K-1,K))
    allocate(KKy(K-1,K))
    allocate(dphi3(K-1,K))
    allocate(k3(K-1,K))
    allocate(cg3(K-1,K))
    allocate(theta3(K-1,K))
    ! Allocate variables for amplitude and Fourier coefficients of bound long wave
    allocate(Gn(wp%tslen))
    allocate(Abnd(K-1,K))
    allocate(Ftemp(K-1,K,4)) ! Jaap qx, qy qtot, zeta
    ! Storage for output discharge
    allocate(q(s%ny+1,wp%tslen,4))   ! qx qy qtot, zeta
    !
    ! Initialize variables as zero
    Eforc = 0
    D = 0
    deltheta = 0
    KKx = 0
    KKy = 0
    dphi3 = 0
    k3 = 0
    cg3 = 0
    theta3 = 0
    Gn = 0*par%compi
    Abnd = 0
    Ftemp = 0*par%compi
    q = 0

    ! First time is set true for each time new wave bc are generated
    firsttime = .true.

    ! upper half of frequency space
    halflen = wp%tslen/2

    ! Run loop over wave-wave interaction components
    call progress_indicator(.true.,0.d0,5.d0,2.d0)
    do m=1,K-1
       !       call writelog('ls','(a,i0,a,i0)','Wave component ',m,' of ',K-1)
       call progress_indicator(.false.,dble(m)/(K-1)*100,5.d0,2.d0)
       ! Allocate memory
       allocate(term1(K-m),term2(K-m),term2new(K-m),dif(K-m),chk1(K-m),chk2(K-m))

       ! Determine difference frequency
       deltaf=m*wp%dfgen

       ! Determine difference angles (pi already added)
       deltheta(m,1:K-m) = abs(wp%thetagen(m+1:K)-wp%thetagen(1:K-m))+par%px

       ! Determine x- and y-components of wave numbers of difference waves
       KKy(m,1:K-m)=wp%kgen(m+1:K)*dsin(wp%thetagen(m+1:K))-wp%kgen(1:K-m)*dsin(wp%thetagen(1:K-m))
       KKx(m,1:K-m)=wp%kgen(m+1:K)*dcos(wp%thetagen(m+1:K))-wp%kgen(1:K-m)*dcos(wp%thetagen(1:K-m))

       ! Determine difference wave numbers according to Van Dongeren et al. 2003
       ! eq. 19
       k3(m,1:K-m) =sqrt(wp%kgen(1:K-m)**2+wp%kgen(m+1:K)**2+ &
            2*wp%kgen(1:K-m)*wp%kgen(m+1:K)*dcos(deltheta(m,1:K-m)))

       ! Determine group velocity of difference waves
       cg3(m,1:K-m)= 2.d0*par%px*deltaf/k3(m,1:K-m)

       ! Modification Robert + Jaap: make sure that the bound long wave amplitude does not
       !                             explode when offshore boundary is too close to shore,
       !                             by limiting the interaction group velocity
       cg3(m,1:K-m) = min(cg3(m,1:K-m),par%nmax*sqrt(par%g/k3(m,1:K-m)*tanh(k3(m,1:K-m)*wp%h0)))

       ! Determine difference-interaction coefficient according to Herbers 1994
       ! eq. A5
       term1 = (-wp%wgen(1:K-m))*wp%wgen(m+1:K)
       term2 = (-wp%wgen(1:K-m))+wp%wgen(m+1:K)
       term2new = cg3(m,1:K-m)*k3(m,1:K-m)
       dif = (abs(term2-term2new))
       if (any(dif>0.01*term2) .and. firsttime) then
          firsttime = .false.
          !          call writelog('lws','','Warning: shallow water so long wave variance is reduced using par%nmax')
       endif
       chk1  = cosh(wp%kgen(1:K-m)*wp%h0)
       chk2  = cosh(wp%kgen(m+1:K)*wp%h0)

       D(m,1:K-m) = -par%g*wp%kgen(1:K-m)*wp%kgen(m+1:K)*dcos(deltheta(m,1:K-m))/2.d0/term1+par%g*term2*(chk1*chk2)/ &
            ((par%g*k3(m,1:K-m)*tanh(k3(m,1:K-m)*wp%h0)-(term2new)**2)*term1*cosh(k3(m,1:K-m)*wp%h0))* &
            (term2*((term1)**2/par%g/par%g - wp%kgen(1:K-m)*wp%kgen(m+1:K)*dcos(deltheta(m,1:K-m))) &
            - 0.50d0*((-wp%wgen(1:K-m))*wp%kgen(m+1:K)**2/(chk2**2)+wp%wgen(m+1:K)*wp%kgen(1:K-m)**2/(chk1**2)))

       ! Correct for surface elevation input and output instead of bottom pressure
       ! so it is consistent with Van Dongeren et al 2003 eq. 18
       D(m,1:K-m) = D(m,1:K-m)*cosh(k3(m,1:K-m)*wp%h0)/(cosh(wp%kgen(1:K-m)*wp%h0)*cosh(wp%kgen(m+1:K)*wp%h0))

       ! Exclude interactions with components smaller than or equal to current
       ! component according to lower limit Herbers 1994 eq. 1
       where(wp%fgen<=deltaf) D(m,:)=0.d0

       ! Exclude interactions with components that are cut-off by the fcutoff
       ! parameter
       if (deltaf<=par%fcutoff) D(m,:)=0.d0

       ! Determine phase of bound long wave assuming a local equilibrium with
       ! forcing of interacting primary waves according to Van Dongeren et al.
       ! 2003 eq. 21 (the angle is the imaginary part of the natural log of a
       ! complex number as long as the complex number is not zero)
       allocate(Comptemp(K-m),Comptemp2(K-m))
       Comptemp=conjg(wp%CompFn(1,wp%Findex(1)+m:wp%Findex(1)+K-1))
       Comptemp2=conjg(wp%CompFn(1,wp%Findex(1):wp%Findex(1)+K-m-1))
       dphi3(m,1:K-m) = par%px+imag(log(Comptemp))-imag(log(Comptemp2))
       deallocate (Comptemp,Comptemp2)
       !
       ! Determine angle of bound long wave according to Van Dongeren et al. 2003 eq. 22
       theta3 = atan2(KKy,KKx)
       !
       ! free memory
       deallocate(term1,term2,term2new,dif,chk1,chk2)
    enddo ! m=1,K-1
    !
    ! Output to screen
    if (.not. firsttime) then
       call writelog('lws','','Warning: shallow water so long wave variance is reduced using par%nmax')
    endif
    call writelog('sl','', 'Calculating flux at boundary')
    !
    ! Allocate temporary arrays for upcoming loop
    allocate(Comptemp(halflen-1))
    allocate(Comptemp2(wp%tslen))
    !
    ! Run a loop over the offshore boundary
    do j=1,s%ny+1
       ! Determine energy of bound long wave according to Herbers 1994 eq. 1 based
       ! on difference-interaction coefficient and energy density spectra of
       ! primary waves
       ! Robert: E = 2*D**2*S**2*dtheta**2*df can be rewritten as
       !         E = 2*D**2*Sf**2*df
       Eforc = 0
       do m=1,K-1
          Eforc(m,1:K-m) = 2*D(m,1:K-m)**2*wp%Sfinterpq(j,1:K-m)*wp%Sfinterpq(j,m+1:K)*wp%dfgen
       enddo
       !
       ! Calculate bound wave amplitude for this offshore grid point
       Abnd = sqrt(2*Eforc*wp%dfgen)
       !
       ! Determine complex description of bound long wave per interaction pair of
       ! primary waves for first y-coordinate along seaside boundary
       Ftemp(:,:,1) = Abnd/2*exp(-1*par%compi*dphi3)*cg3*dcos(theta3) ! qx
       Ftemp(:,:,2) = Abnd/2*exp(-1*par%compi*dphi3)*cg3*dsin(theta3) ! qy
       Ftemp(:,:,3) = Abnd/2*exp(-1*par%compi*dphi3)*cg3              ! qtot
       Ftemp(:,:,4) = Abnd/2*exp(-1*par%compi*dphi3)                  ! eta
       !
       ! loop over qx,qy and qtot
       do iq=1,4
          ! Unroll wave component to correct place along the offshore boundary
          Ftemp(:,:,iq) = Ftemp(:,:,iq)* &
               exp(-1*par%compi*(KKy*(s%yz(1,j)-s%yz(1,1))+KKx*(s%xz(1,j)-s%xz(1,1))))
          ! Determine Fourier coefficients
          Gn(2:K) = sum(Ftemp(:,:,iq),DIM=2)
          Comptemp = conjg(Gn(2:halflen))
          call flipiv(Comptemp,halflen-1)
          Gn(halflen+2:wp%tslen) = Comptemp
          !
          ! Print status message to screen
          if (iq==3) then
             call writelog('ls','(A,I0,A,I0)','Flux ',j,' of ',s%ny+1)
          endif
          !
          ! Inverse Discrete Fourier transformation to transform back to time space
          ! from frequency space
          Comptemp2=fft(Gn,inv=.true.)
          !
          ! Determine mass flux as function of time and let the flux gradually
          ! increase and decrease in and out the wave time record using the earlier
          ! specified window
          Comptemp2=Comptemp2/sqrt(dble(wp%tslen))
          q(j,:,iq)=dreal(Comptemp2*wp%tslen)*wp%taperf
       enddo ! iq=1,3
    enddo ! j=1,s%ny+1
    !
    ! free memory
    deallocate(Comptemp,Comptemp2,Ftemp)
    !
    if (par%nonhspectrum==0) then
       ! Open file for storage of bound long wave flux
       call writelog('ls','','Writing long wave mass flux to ',trim(wp%qfilename),' ...')
       inquire(iolength=reclen) 1.d0
       reclen=reclen*((s%ny+1)*4)
       fid = create_new_fid()
       open(fid,file=trim(wp%qfilename),form='unformatted',access='direct',recl=reclen,status='REPLACE')
       !
       ! Write to external file
       if (wp%dtchanged) then
          ! Interpolate from internal time axis to output time axis
          allocate(qinterp(s%ny+1,wp%tslenbc,3))
          do iq=1,3
             do it=1,wp%tslenbc
                do j=1,s%ny+1
                   call linear_interp(wp%tin,q(j,:,iq),wp%tslen,(it-1)*wp%dtbc,qinterp(j,it,iq),status)
                enddo
             enddo
          enddo
          ! write to file
          do irec=1,wp%tslenbc+1
             write(fid,rec=irec)qinterp(:,min(irec,wp%tslenbc),:)
          end do
          deallocate(qinterp)
       else
          ! no need for interpolation
          do irec=1,wp%tslenbc+1
             write(fid,rec=irec)q(:,min(irec,wp%tslenbc),:)
          end do
       endif
       close(fid)
       call writelog('sl','','file done')
    else
       do j=1,s%ny+1
          ! add to velocity time series
          wp%uits(j,:)=wp%uits(j,:)+q(j,:,1)/wp%h0
          ! add to surface elevation time series
          wp%zsits(j,:)=wp%zsits(j,:)+q(j,:,4)
       enddo
    endif ! par%nonhspectrum==1

    ! Free memory
    deallocate(Eforc,D,deltheta,KKx,KKy,dphi3,k3,cg3,theta3,Gn,Abnd,q)

  end subroutine generate_qbcf


  ! --------------------------------------------------------------
  ! --------------- Non-hydrostatic wave time --------------------
  ! ---------------- series file generation ----------------------
  subroutine generate_nhtimeseries_file(wp,par)

    use params
    use filefunctions, only: create_new_fid
    use logging_module

    implicit none

    ! input/output
    type(waveparamsnew),intent(inout)            :: wp
    type(parameters),intent(in)                  :: par
    ! internal
    integer                                      :: fid
    integer                                      :: it,j


    call writelog('ls','','Writing short wave time series to ',wp%nhfilename)

    fid = create_new_fid()
    open(fid,file=trim(wp%nhfilename),status='REPLACE')
    write(fid,'(a)')'vector'
    write(fid,'(a)')'3'
    write(fid,'(a)')'t,U,Z'
    do it=1,wp%tslen
       write(fid,*)wp%tin(it)+par%t-par%dt,wp%uits(:,it),wp%zsits(:,it)
       if (wp%tin(it)>wp%rtbc) exit
    enddo
    close(fid)

  end subroutine generate_nhtimeseries_file


end module spectral_wave_bc_module
