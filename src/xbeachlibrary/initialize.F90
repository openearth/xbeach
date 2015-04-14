module initialize 
  use typesandkinds
  implicit none
  save
  integer imin_ee,imax_ee,jmin_ee,jmax_ee
  integer imin_uu,imax_uu,jmin_uu,jmax_uu
  integer imin_vv,imax_vv,jmin_vv,jmax_vv
  integer imin_zs,imax_zs,jmin_zs,jmax_zs

contains
  subroutine grid_bathy(s,par)                                         

    use params   
    use spaceparams                       
    use xmpi_module
    use general_mpi_module
    use readkey_module
    use logging_module
    use paramsconst
#ifdef USEMPI
    use mpi
#endif

    implicit none                                                            

    type(spacepars)                :: s               
    type(parameters)               :: par               

    integer                        :: i,ier,idum,n,m,ier2,ier3,dum
    integer                        :: j
    integer                        :: itheta
    real*8                         :: degrad
    character(slen)                :: line

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

    if(.not. (xmaster .or. xomaster)) return
    if(.true.) then
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

    call writelog('l','' ,'------------------------------------')
    call writelog('ls','','Building Grid and Bathymetry')
    call writelog('l','', '------------------------------------')

    !
    ! Create grid and bathymetry
    !
    ! Jaap make switch here to read XBeach or Delft3D format respectively 
    !
    ! wv in the following select case construct, s%zb, s%x and s%y are determined
    !  and they are read from file
    ! we let xmaster read these entities and send them to xomaster.
    !  because xomaster has also a need for these entities
    !

    if (xmaster) then
      select case(par%gridform)
       case(GRIDFORM_XBEACH)
         select case(s%vardx)
           case(0)
            if (par%setbathy .ne. 1) then
             open(31,file=par%depfile)
             do j=1,s%ny+1
                read(31,*,iostat=ier)(s%zb(i,j),i=1,s%nx+1)
                if (ier .ne. 0) then
                   call report_file_read_error(par%depfile)
                endif
             end do
             close(31)
             endif
             do j=1,s%ny+1
                do i=1,s%nx+1
                   s%x(i,j)=(i-1)*s%dx
                   s%y(i,j)=(j-1)*s%dy
                end do
             end do
             case (1)   ! Robert keep vardx == 1 for backwards compatibility??
             if (par%setbathy .ne. 1) then
             open (31,file=par%depfile)
             read (31,*,iostat=ier)((s%zb(i,j),i=1,s%nx+1),j=1,s%ny+1)
             if (ier .ne. 0) then
                call report_file_read_error(par%depfile)
             endif
             close(31)
             endif

             open (32,file=par%xfile)
             read (32,*,iostat=ier)((s%x(i,j),i=1,s%nx+1),j=1,s%ny+1)
             if (ier .ne. 0) then
                call report_file_read_error(par%xfile)
             endif
             close(32)

             if (s%ny>0 .and. par%yfile/=' ') then
                open (33,file=par%yfile)
                read (33,*,iostat=ier)((s%y(i,j),i=1,s%nx+1),j=1,s%ny+1)
                if (ier .ne. 0) then
                   call report_file_read_error(par%yfile)
                endif
                close(33)
             elseif (s%ny==0 .and. par%yfile/=' ') then
                open (33,file=par%yfile)
                read (33,*,iostat=ier)((s%y(i,j),i=1,s%nx+1),j=1,s%ny+1)
                if (ier .ne. 0) then
                   call report_file_read_error(par%yfile)
                endif
                close(33)
             else
                s%y = 0.d0
             end if
               !endif
             case default
          call writelog('esl','','Invalid value for vardx: ',par%vardx)
          call halt_program
            end select

          case (GRIDFORM_DELFT3D)
          ! 
          ! Gridfile
          !
          open(31,file=par%xyfile,status='old',iostat=ier)  
            if (ier .ne. 0) then
               call report_file_read_error(par%xyfile)
            endif
            ! http://oss.deltares.nl/documents/183920/185723/Delft3D-FLOW_User_Manual.pdf section A.2.2
          ! skip comment text in file...
            do
             read(31,'(a)',iostat=ier)line
             if (line(1:1)/='*') then 
                  exit
             endif
          enddo
          read(31,*,iostat=ier2) idum,idum
          ! new grid format now specifies missing value, so catch this error
          if (ier2 .ne. 0) then
             read(31,*,iostat=ier2) idum,idum
          endif
          read(31,*,iostat=ier3) idum,idum,idum
          ! if any iostat is still /= 0 then there is an error reading the file
          if (ier+ier2+ier3 .ne. 0) then
             call report_file_read_error(par%xyfile)
          endif
          
          read(31,*,iostat=ier) &
              (line,dum,(s%x(m,n),m=1,s%nx+1),n=1,s%ny+1), &
              (line,dum,(s%y(m,n),m=1,s%nx+1),n=1,s%ny+1) 
            if (ier .ne. 0) then
               call report_file_read_error(par%xyfile)
            endif

            close(31)
          !
          ! Depfile
          !
          if (par%setbathy .ne. 1) then
          open(33,file=par%depfile,status='old')
          do n=1,s%ny+1
             read(33,*,iostat=ier)(s%zb(m,n),m=1,s%nx+1)
             if (ier .ne. 0) then
                call report_file_read_error(par%depfile)
             endif
          enddo
             close(33)
          endif
         end select
      endif
      ! xmaster

    degrad=par%px/180.d0
      ! send zb, x, y to xomaster

#ifdef USEMPI
      ! wwvv todo  use xmpi_send
    if(xmaster) then
         ! MPI_SEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
         ! <type>    BUF(*)
         ! INTEGER    COUNT, DATATYPE, DEST, TAG, COMM, IERROR
         call MPI_Send(s%zb, size(s%zb), MPI_DOUBLE_PRECISION, xmpi_omaster, 1, xmpi_ocomm, ier)
         call MPI_Send(s%x,  size(s%x),  MPI_DOUBLE_PRECISION, xmpi_omaster, 2, xmpi_ocomm, ier)
         call MPI_Send(s%y,  size(s%y),  MPI_DOUBLE_PRECISION, xmpi_omaster, 3, xmpi_ocomm, ier)
      else
         ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)
         ! <type>    BUF(*)
         ! INTEGER    COUNT, DATATYPE, SOURCE, TAG, COMM
         ! INTEGER    STATUS(MPI_STATUS_SIZE), IERROR
         call MPI_Recv(s%zb, size(s%zb), MPI_DOUBLE_PRECISION, xmpi_imaster, 1, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
         call MPI_Recv(s%x,  size(s%x),  MPI_DOUBLE_PRECISION, xmpi_imaster, 2, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
         call MPI_Recv(s%y,  size(s%y),  MPI_DOUBLE_PRECISION, xmpi_imaster, 3, xmpi_ocomm, MPI_STATUS_IGNORE, ier)
      endif
#endif

      if(.true.) then

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

       !degrad=par%px/180.d0

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
      
       if (par%swave==1) then
          s%dtheta=par%dtheta*degrad
          s%ntheta=nint((s%thetamax-s%thetamin)/s%dtheta)
       else
          s%dtheta=par%dtheta*degrad
          s%ntheta = 1
       endif

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

       if (par%single_dir==1) then
          s%dtheta_s=par%dtheta_s*degrad
          s%ntheta_s=nint((s%thetamax-s%thetamin)/s%dtheta_s)

          allocate(s%theta_s(1:s%ntheta_s))
          allocate(s%thet_s(1:s%nx+1,1:s%ny+1,1:s%ntheta_s))
          allocate(s%costh_s(1:s%nx+1,1:s%ny+1,1:s%ntheta_s))
          allocate(s%sinth_s(1:s%nx+1,1:s%ny+1,1:s%ntheta_s))

          do itheta=1,s%ntheta_s
             s%theta_s(itheta)=s%thetamin+s%dtheta_s/2+s%dtheta_s*(itheta-1)
          end do

          do itheta=1,s%ntheta_s
             do j=1,s%ny+1
                do i=1,s%nx+1
                   s%thet_s(i,j,itheta) = s%theta_s(itheta)
                   s%costh_s(i,j,itheta)=cos(s%theta_s(itheta)-s%alfaz(i,j))
                   s%sinth_s(i,j,itheta)=sin(s%theta_s(itheta)-s%alfaz(i,j))
                enddo
             enddo
          enddo
      else
          s%ntheta_s = 0
          allocate(s%theta_s(s%ntheta_s))
          allocate(s%thet_s (s%nx+1,s%ny+1,s%ntheta_s))
          allocate(s%costh_s(s%nx+1,s%ny+1,s%ntheta_s))
          allocate(s%sinth_s(s%nx+1,s%ny+1,s%ntheta_s))
      endif
       
    endif

      !if (xmaster) then
      if(.true.) then
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

  subroutine setbathy_init(s,par)
    use params
    use spaceparams
    use filefunctions
    use logging_module
    use interp
    
    type(spacepars)                     :: s
    type(parameters)                    :: par
    
    integer                             :: i,j,it
    integer                             :: ier,fid,dummy
    
    if(.not. xmaster) return

    if (par%setbathy==1) then
       allocate(s%setbathy(s%nx+1,s%ny+1,par%nsetbathy))
       allocate(s%tsetbathy(par%nsetbathy))
       ! start file read
       fid = create_new_fid()
       open (fid,file=par%setbathyfile)
       do it=1,par%nsetbathy
          read(fid,*,iostat=ier)s%tsetbathy(it)
          if (ier .ne. 0) then
             call report_file_read_error(par%setbathyfile)
          endif
          do j=1,s%ny+1
             read(fid,*,iostat=ier)(s%setbathy(i,j,it),i=1,s%nx+1)
             if (ier .ne. 0) then
                call report_file_read_error(par%setbathyfile)
             endif
          enddo
       enddo
       close(fid)
       ! Interpolate initial bathymetry
       do j=1,s%ny+1
          do i=1,s%nx+1
             call LINEAR_INTERP(s%tsetbathy,s%setbathy(i,j,:),par%nsetbathy, &
                                 0.d0,s%zb(i,j),dummy)  
          enddo
       enddo 
    else
       ! give MPI bcast a memory address 
       allocate(s%setbathy(0,0,0))
       allocate(s%tsetbathy(0))
    endif
  end subroutine setbathy_init

  subroutine wave_init (s,par)
    use params
    use spaceparams
    use readkey_module
    use xmpi_module
    use wave_functions_module
    use paramsconst

    IMPLICIT NONE

    type(spacepars),target              :: s
    type(parameters)                    :: par

    integer                             :: itheta


    if(.not. xmaster) return

    allocate(s%thetamean(1:s%nx+1,1:s%ny+1))
    allocate(s%Fx(1:s%nx+1,1:s%ny+1))
    allocate(s%Fy(1:s%nx+1,1:s%ny+1))
    allocate(s%Sxx(1:s%nx+1,1:s%ny+1))
    allocate(s%Sxy(1:s%nx+1,1:s%ny+1))
    allocate(s%Syy(1:s%nx+1,1:s%ny+1))
    allocate(s%n(1:s%nx+1,1:s%ny+1))
    allocate(s%H(1:s%nx+1,1:s%ny+1))
    allocate(s%cgx(1:s%nx+1,1:s%ny+1,1:s%ntheta))
    allocate(s%cgy(1:s%nx+1,1:s%ny+1,1:s%ntheta))
    allocate(s%cx(1:s%nx+1,1:s%ny+1,1:s%ntheta))
    allocate(s%cy(1:s%nx+1,1:s%ny+1,1:s%ntheta))
    allocate(s%ctheta(1:s%nx+1,1:s%ny+1,1:s%ntheta))
    allocate(s%sigt(1:s%nx+1,1:s%ny+1,1:s%ntheta))
    allocate(s%ee(1:s%nx+1,1:s%ny+1,1:s%ntheta))
    allocate(s%rr(1:s%nx+1,1:s%ny+1,1:s%ntheta))
    if (par%single_dir==1) then
       allocate(s%cgx_s(1:s%nx+1,1:s%ny+1,1:s%ntheta_s))
       allocate(s%cgy_s(1:s%nx+1,1:s%ny+1,1:s%ntheta_s))
       allocate(s%ctheta_s(1:s%nx+1,1:s%ny+1,1:s%ntheta_s))
       allocate(s%ee_s(1:s%nx+1,1:s%ny+1,1:s%ntheta_s))
    endif
    allocate(s%sigm(1:s%nx+1,1:s%ny+1))
    allocate(s%c(1:s%nx+1,1:s%ny+1))
    allocate(s%cg(1:s%nx+1,1:s%ny+1))
    allocate(s%k(1:s%nx+1,1:s%ny+1))
    allocate(s%ui(1:s%nx+1,1:s%ny+1))
    allocate(s%vi(1:s%nx+1,1:s%ny+1))
    allocate(s%E(1:s%nx+1,1:s%ny+1)) 
    allocate(s%R(1:s%nx+1,1:s%ny+1)) 
    allocate(s%urms(1:s%nx+1,1:s%ny+1)) 
    allocate(s%D(1:s%nx+1,1:s%ny+1)) 
    allocate(s%Df(1:s%nx+1,1:s%ny+1)) 
    allocate(s%Dveg(1:s%nx+1,1:s%ny+1))
    allocate(s%Fvegu(1:s%nx+1,1:s%ny+1))
    allocate(s%Fvegv(1:s%nx+1,1:s%ny+1))
    allocate(s%Dp(1:s%nx+1,1:s%ny+1)) 
    allocate(s%Qb(1:s%nx+1,1:s%ny+1)) 
    allocate(s%ust(1:s%nx+1,1:s%ny+1)) 
    allocate(s%tm(1:s%nx+1,1:s%ny+1)) 
    allocate(s%uwf(1:s%nx+1,1:s%ny+1)) 
    allocate(s%vwf(1:s%nx+1,1:s%ny+1)) 
    allocate(s%ustr(1:s%nx+1,1:s%ny+1)) 
    allocate(s%usd(1:s%nx+1,1:s%ny+1))
    allocate(s%bi(1:s%ny+1))
    allocate(s%DR(1:s%nx+1,1:s%ny+1)) 
    allocate(s%umwci       (1:s%nx+1,1:s%ny+1))
    allocate(s%vmwci       (1:s%nx+1,1:s%ny+1))
    allocate(s%zswci       (1:s%nx+1,1:s%ny+1))
    allocate(s%BR(1:s%nx+1,1:s%ny+1))
    !
    ! Initial condition
    !
    s%ee        = 0.d0 
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
    if (par%single_dir==1) then
       s%ee_s      = 0.d0
       s%cgx_s     = 0.d0
       s%cgy_s     = 0.d0
       s%ctheta_s  = 0.d0
    endif
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
    s%Df        = 0.d0
    s%Dveg      = 0.d0
    s%Fvegu     = 0.d0
    s%Fvegv     = 0.d0
    s%BR        = par%Beta


    ! introduce intrinsic frequencies for wave action
    if ( par%instat==INSTAT_JONS .or. &
         par%instat==INSTAT_JONS_TABLE .or. &
         par%instat==INSTAT_SWAN .or. &
         par%instat==INSTAT_VARDENS .or. &
         par%instat==INSTAT_REUSE .or. &
         par%instat==INSTAT_TS_NONH &
         ) par%Trep=10.d0 
    !Robert
    ! incorrect values are computed below for instat = 4/5/6/7
    ! in this case right values are computed in wave params.f90
    if ( par%instat==INSTAT_JONS .or. &
         par%instat==INSTAT_JONS_TABLE .or. &
         par%instat==INSTAT_SWAN .or. &
         par%instat==INSTAT_VARDENS) then 
       if(xmaster) call spectral_wave_init (s,par)  ! only used by xmaster
    endif
    do itheta=1,s%ntheta
       s%sigt(:,:,itheta) = 2*par%px/par%Trep
    end do
    s%sigm = sum(s%sigt,3)/s%ntheta
    call dispersion(par,s,.false.)


  end subroutine wave_init


  subroutine spectral_wave_init(s,par)
    use params
    use spaceparams
    use filefunctions
    use logging_module
    use spectral_wave_bc_module
  
    type(spacepars)             :: s
    type(parameters)            :: par
  
    integer                     :: fid,err
    integer                     :: i
    integer,dimension(1)        :: minlocation
    character(slen)             :: testline
    real*8,dimension(:),allocatable :: xspec,yspec,mindist
    real*8                      :: mindistr
  
    call writelog('l','','--------------------------------')
    call writelog('l','','Initializing spectral wave boundary conditions ')
    ! Initialize that wave boundary conditions need to be calculated (first time at least)
    ! Stored and defined in spectral_wave_bc_module
    reuseall = .false.
    ! Initialize the number of times wave boundary conditions have been generated.
    ! Stored and defined in spectral_wave_bc_module
    bccount  = 0
    ! Initialize bcendtime to zero.
    ! Stored and defined in spectral_wave_bc_module
    spectrumendtime = 0.d0
    ! Initialise lastwaveheight to zero
    ! Stored and defined in spectral_wave_bc_module
    allocate(lastwaveelevation(s%ny+1,s%ntheta))
    lastwaveelevation = 0.d0
  
    if (par%nspectrumloc<1) then
       call writelog('ewls','','number of boundary spectra (''nspectrumloc'') may not be less than 1')
       call halt_program
    endif
  
    ! open location list file
    fid = create_new_fid()
    open(fid,file=par%bcfile,status='old',form='formatted')
    ! check for LOCLIST
    read(fid,*)testline
    if (trim(testline)=='LOCLIST') then
       allocate(n_index_loc(par%nspectrumloc)) ! stored and defined in spectral_wave_bc_module
       allocate(bcfiles(par%nspectrumloc))     ! stored and defined in spectral_wave_bc_module
       allocate(xspec(par%nspectrumloc))
       allocate(yspec(par%nspectrumloc))
       allocate(mindist(s%ny+1))
       do i=1,par%nspectrumloc
          ! read x,y and file name per location
          read(fid,*,IOSTAT=err)xspec(i),yspec(i),bcfiles(i)%fname
          bcfiles(i)%listline = 0
          if (err /= 0) then
             ! something has gone wrong during the read of this file
             call writelog('ewls','a,i0,a,a)','error reading line ',i+1,' of file ',par%bcfile)
             call writelog('ewls','','check file for format errors and ensure the number of  ',&
                  'lines is equal to nspectrumloc')
             call halt_program
          endif
       enddo
       ! convert x,y locations of the spectra to the nearest grid points on s=1 boundary
       do i=1,par%nspectrumloc
          mindist=sqrt((xspec(i)-s%xz(1,:))**2+(yspec(i)-s%yz(1,:))**2)
          ! look up the location of the found minimum
          minlocation=minloc(mindist)
          ! minimum distance
          mindistr = mindist(minlocation(1))
          ! store location
          n_index_loc(i) = minlocation(1)
          call writelog('ls','(a,i0,a,i0)','Spectrum ',i,' placed at n = ',minlocation(1))
          call writelog('ls','(a,f0.3)','Distance spectrum to grid: ',mindistr)
       enddo
       nspectra = par%nspectrumloc     ! stored and defined in spectral_wave_bc_module
       deallocate(xspec,yspec,mindist)
    else
       if (par%nspectrumloc==1) then
          allocate(n_index_loc(1)) ! stored and defined in spectral_wave_bc_module
          allocate(bcfiles(1))     ! stored and defined in spectral_wave_bc_module
          n_index_loc = 1
          bcfiles(1)%fname = par%bcfile
          bcfiles(1)%listline = 0      ! for files that have multiple lines, set listline to 0
          nspectra = 1                 ! stored and defined in spectral_wave_bc_module
       else
          call writelog('ewls','','if nspectrumloc>1 then bcfile should contain spectra locations with LOCLIST header')
          close(fid)
          call halt_program
       endif
    endif
    close(fid)
  
    call writelog('l','','--------------------------------')
  
  end subroutine spectral_wave_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine flow_init (s,par)
    use params
    use spaceparams
    use readkey_module
    use logging_module
    use interp
    use xmpi_module
    use paramsconst

    IMPLICIT NONE

    type(spacepars),target                  :: s
    type(parameters), intent(in)            :: par

    integer*4                               :: i,j,ig,indt
    logical                                 :: exists
    logical                                 :: offshoreregime
    integer                                 :: indoff,indbay,ier

    real*8,dimension(:),allocatable         :: xzs0,yzs0,szs0
    real*8,dimension(:,:),allocatable       :: vmagvold,vmaguold
    real*8                                  :: flowerr


    if(.not. xmaster) return

    allocate(s%zs(1:s%nx+1,1:s%ny+1))
    allocate(s%dzsdt(1:s%nx+1,1:s%ny+1))
    allocate(s%dzsdx(1:s%nx+1,1:s%ny+1))
    allocate(s%dzsdy(1:s%nx+1,1:s%ny+1))
    allocate(s%dzbdt(1:s%nx+1,1:s%ny+1))
    allocate(s%dzbnow(1:s%nx+1,1:s%ny+1))
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
    allocate(s%wscrit(1:s%nx+1,1:s%ny+1))
    allocate(s%wb(1:s%nx+1,1:s%ny+1))
    allocate(s%nuh(1:s%nx+1,1:s%ny+1))  
    allocate(s%pres(1:s%nx+1,1:s%ny+1))
    allocate(s%ph(1:s%nx+1,1:s%ny+1))
    allocate(s%wi(2,1:s%ny+1))
    allocate(s%zi(2,1:s%ny+1))
    allocate(s%bedfriccoef(1:s%nx+1,1:s%ny+1))
    allocate(s%cf(1:s%nx+1,1:s%ny+1))
    allocate(s%cfu(1:s%nx+1,1:s%ny+1))
    allocate(s%cfv(1:s%nx+1,1:s%ny+1))
    allocate(s%zs0(1:s%nx+1,1:s%ny+1)) 
    allocate(s%zs0fac(1:s%nx+1,1:s%ny+1,2))
    allocate(s%wm(1:s%nx+1,1:s%ny+1))
    allocate(s%umean(1:s%nx+1,1:s%ny+1))
    allocate(s%vmean(1:s%nx+1,1:s%ny+1))
    allocate(s%ur(1:s%nx+1,1:s%ny+1))
    allocate(s%xyzs01(2))
    allocate(s%xyzs02(2))
    allocate(s%xyzs03(2))
    allocate(s%xyzs04(2))

    allocate(szs0(1:2)) 
    allocate(xzs0(1:2)) 
    allocate(yzs0(1:2)) 
    allocate(vmagvold(1:s%nx+1,1:s%ny+1)) 
    allocate(vmaguold(1:s%nx+1,1:s%ny+1)) 

    allocate(s%ududx(1:s%nx+1,1:s%ny+1))
    allocate(s%vdvdy(1:s%nx+1,1:s%ny+1))
    allocate(s%udvdx(1:s%nx+1,1:s%ny+1))
    allocate(s%vdudy(1:s%nx+1,1:s%ny+1))
    allocate(s%viscu(1:s%nx+1,1:s%ny+1))
    allocate(s%viscv(1:s%nx+1,1:s%ny+1))
    ! nonh hydrostatic flow initialisation
    allocate(s%breaking(1:s%nx+1,1:s%ny+1))    

    ! Just to be sure!
    s%zs = 0.0d0
    s%dzsdt = 0.0d0
    s%dzsdx = 0.0d0
    s%dzsdy = 0.0d0
    s%dzbdt = 0.0d0
    s%dzbnow = 0.d0
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
    s%ph = 0.0d0 
    s%wi = 0.0d0 
    s%zi = 0.0d0 
    s%cf = 0.0d0
    s%zs0 = 0.0d0
    s%zs0fac = 0.0d0
    s%wm = 0.0d0
    s%umean = 0.0d0
    s%vmean = 0.0d0
    s%ur = 0.0d0
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
    s%wscrit = 0.d0
    s%wb   = 0.0d0
    s%pres = 0.0d0
    s%zi   = 0.0d0
    s%wi   = 0.0d0
    s%nuh  = 0.0d0
    s%cf   = par%cf
    s%wm   =0.d0
    s%zs0fac = 0.0d0

    s%ududx = 0.d0
    s%vdvdy = 0.d0
    s%udvdx = 0.d0
    s%vdudy = 0.d0
    s%viscu = 0.d0
    s%viscv = 0.d0   

    if (par%nonh==1) then
       s%breaking = 0
    endif

    !
    ! set-up tide and surge waterlevels
    s%zs01=par%zs0
    ! I: read zs0 at model corners using zs0file
    if (par%tideloc>0) then

       ! Need to interpolate to the correct moment in time. First point in tidal 
       ! record not necessarily == 0.0
       call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,1),s%tidelen,0.d0, s%zs01, indt)

       if(par%tideloc.eq.1) s%zs02=s%zs01

       if(par%tideloc.eq.2 .and. par%paulrevere==PAULREVERE_LAND) then
          call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,2),s%tidelen,0.d0, s%zs03, indt)
          s%zs02=s%zs01
          s%zs04=s%zs03
       endif

       if(par%tideloc.eq.2 .and. par%paulrevere==PAULREVERE_SEA) then
          call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,2),s%tidelen,0.d0, s%zs02, indt)
          s%zs03=0.d0
          s%zs04=0.d0
       endif

       if(par%tideloc.eq.4) then
          call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,2),s%tidelen,0.d0, s%zs02, indt)
          call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,3),s%tidelen,0.d0, s%zs03, indt)
          call LINEAR_INTERP(s%tideinpt,s%tideinpz(:,4),s%tidelen,0.d0, s%zs04, indt)
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

       if(par%tideloc.eq.2 .and. par%paulrevere==PAULREVERE_SEA) then
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

       if(par%tideloc.eq.2 .and. par%paulrevere==PAULREVERE_LAND) then
          yzs0(1)=s%sdist(1,1)
          yzs0(2)=s%sdist(s%nx+1,1)
          szs0(1)=s%zs01
          szs0(2)=s%zs04
          s%zs0(1,:)=s%zs01
          s%zs0(s%nx+1,:)=s%zs03
          do j = 1,s%ny+1 
             indoff = s%nx+1
             indbay = 1
             ! look for intersect of bed with zs offshore
             do i = 2,s%nx+1
                if (s%zb(i,j).gt.s%zs0(1,j)+par%eps) then
                   indoff = i-1
                   exit
                endif
             enddo
             ! look for intersect of bed with zs bay
             do i = s%nx,2,-1
                if (s%zb(i,j).gt.s%zs0(s%nx+1,j)+par%eps) then
                   indbay = i+1
                   exit
                endif
             enddo
             ! do both intersects exist?
             ! apply two water levels in the domain
             if (indoff<s%nx+1 .and. indbay>1) then
                s%zs0(2:indoff,j) = s%zs0(1,j)
                s%zs0(indbay:s%nx,j) = s%zs0(s%nx+1,j)
                ! linear interpolation between intersects
                ! maximize to bed level
                do ig=indoff+1,indbay-1
                   xzs0(1)=s%sdist(indoff,j)
                   xzs0(2)=s%sdist(indbay,j)
                   szs0(1)=s%zs0(indoff,j)
                   szs0(2)=s%zs0(indbay,j)
                   call LINEAR_INTERP(xzs0, szs0, 2, s%sdist(ig,j), s%zs0(ig,j), indt)
                   s%zs0(ig,j)=min(s%zs0(ig,j),s%zb(ig,j))
                enddo
                ! only bay intersect exists -> all land below offshore sea level
             elseif (indoff==s%nx+1 .and. indbay>1) then
                s%zs0(:,j) = s%zs0(1,j)             
                ! only offshore intersect exists -> all land below bay sea level
             elseif (indoff<s%nx+1 .and. indbay==1) then
                s%zs0(:,j) = s%zs0(s%nx+1,j)
                ! no intersects exist -> all land below bay and offshore sea level
                ! linear interpolation between offshore and bay sea level
             else
                do ig=1,s%nx+1
                   xzs0(1)=s%sdist(1,j)
                   xzs0(2)=s%sdist(s%nx+1,j)
                   szs0(1)=s%zs0(1,j)
                   szs0(2)=s%zs0(s%nx+1,j)
                   call LINEAR_INTERP(xzs0, szs0, 2, s%sdist(ig,j), s%zs0(ig,j), indt)
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
             indoff = s%nx+1
             indbay = 1
             ! look for intersect of bed with zs offshore
             do i = 2,s%nx+1
                if (s%zb(i,j).gt.s%zs0(1,j)+par%eps) then
                   indoff = i-1
                   exit
                endif
             enddo
             ! look for intersect of bed with zs bay
             do i = s%nx,2,-1
                if (s%zb(i,j).gt.s%zs0(s%nx+1,j)+par%eps) then
                   indbay = i+1
                   exit
                endif
             enddo
             ! do both intersects exist?
             ! apply two water levels in the domain
             if (indoff<s%nx+1 .and. indbay>1) then
                s%zs0(2:indoff,j) = s%zs0(1,j)
                s%zs0(indbay:s%nx,j) = s%zs0(s%nx+1,j)
                ! linear interpolation between intersects
                ! maximize to bed level
                do ig=indoff+1,indbay-1
                   xzs0(1)=s%sdist(indoff,j)
                   xzs0(2)=s%sdist(indbay,j)
                   szs0(1)=s%zs0(indoff,j)
                   szs0(2)=s%zs0(indbay,j)
                   call LINEAR_INTERP(xzs0, szs0, 2, s%sdist(ig,j), s%zs0(ig,j), indt)
                   s%zs0(ig,j)=min(s%zs0(ig,j),s%zb(ig,j))
                enddo
                ! only bay intersect exists -> all land below offshore sea level
             elseif (indoff==s%nx+1 .and. indbay>1) then
                s%zs0(:,j) = s%zs0(1,j)             
                ! only offshore intersect exists -> all land below bay sea level
             elseif (indoff<s%nx+1 .and. indbay==1) then
                s%zs0(:,j) = s%zs0(s%nx+1,j)
                ! no intersects exist -> all land below bay and offshore sea level
                ! linear interpolation between offshore and bay sea level
             else
                do ig=1,s%nx+1
                   xzs0(1)=s%sdist(1,j)
                   xzs0(2)=s%sdist(s%nx+1,j)
                   szs0(1)=s%zs0(1,j)
                   szs0(2)=s%zs0(s%nx+1,j)
                   call LINEAR_INTERP(xzs0, szs0, 2, s%sdist(ig,j), s%zs0(ig,j), indt)
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
          read(723,*,iostat=ier)(s%zs0(i,j),i=1,s%nx+1)
          if (ier .ne. 0) then
             call report_file_read_error(par%zsinitfile)
          endif
       enddo
       close(723)
    endif

    inquire(file=par%bedfricfile,exist=exists)
    if ((exists)) then
       open(723,file=par%bedfricfile)
       do j=1,s%ny+1
          read(723,*,iostat=ier)(s%bedfriccoef(i,j),i=1,s%nx+1)
          if (ier .ne. 0) then
             call report_file_read_error(par%bedfricfile)
          endif
       enddo
       close(723)
    else
       s%bedfriccoef=par%bedfriccoef
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

    ! R+L: Why are these variable reinitialise? 
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
    ! Start with initial velocities based on balance between water level gradient
    ! and bed friction term
    if (par%hotstartflow==1) then
       call writelog('ls','','Calculating stationary flow field for initial condition')
       ! water level gradients
       do j=1,s%ny+1
          do i=2,s%nx
             s%dzsdx(i,j)=(s%zs0(i+1,j)-s%zs0(i,j))/s%dsu(i,j)
          end do
       end do
       do j=1,s%ny ! Dano need to get correct slope on boundary y=0
          do i=1,s%nx+1
             s%dzsdy(i,j)=(s%zs0(i,j+1)-s%zs0(i,j))/s%dnv(i,j)
          end do
       end do
       ! Water depth in u,v-points mean for wind forcing
       do j=1,s%ny+1
          do i=1,s%nx+1 !Ap
             s%hum(i,j)=max(.5d0*(s%hh(i,j)+s%hh(min(s%nx,i)+1,j)),par%eps) 
          end do
       end do
       do j=1,s%ny+1
          do i=1,s%nx+1
             s%hvm(i,j)=max(.5d0*(s%hh(i,j)+s%hh(i,min(s%ny,j)+1)),par%eps)
          end do
       end do
       ! residual error
       flowerr = huge(0.d0)
       vmagvold = 0.d0
       vmaguold = 0.d0
       do while (flowerr > 0.00001d0)
          !
          ! Balance in longshore
          vmagvold=0.5d0*(vmagvold+sqrt(s%uv**2+s%vv**2))   ! mean needed for convergence
          ! solve v-balance of pressure gradient, wind forcing and bed friction
          where (s%wetv==1)
             s%vv = s%hv/s%cf/max(vmagvold,0.000001d0) &
                  *(-par%g*s%dzsdy+par%rhoa*par%Cd*s%windnv**2/(par%rho*s%hvm))
          elsewhere
             s%vv = 0.d0
          endwhere
          ! update vmagev
          ! u velocity in v points
          if (s%ny>0) then
             s%uv(2:s%nx,1:s%ny)= .25d0*(s%uu(1:s%nx-1,1:s%ny)+s%uu(2:s%nx,1:s%ny)+ &
                  s%uu(1:s%nx-1,2:s%ny+1)+s%uu(2:s%nx,2:s%ny+1))
             ! boundaries?
             ! wwvv and what about uv(:,1) ?
             if(xmpi_isright) then
                s%uv(:,s%ny+1) = s%uv(:,s%ny)
             endif
          else
             s%uv(2:s%nx,1)= .5d0*(s%uu(1:s%nx-1,1)+s%uu(2:s%nx,1))
          endif !s%ny>0
          s%vmagev = sqrt(s%uv**2+s%vv**2)
          !
          ! Balance in cross shore
          vmaguold= 0.5d0*(vmaguold+sqrt(s%uu**2+s%vu**2)) ! mean needed for convergence
          ! Solve balance of forces
          where (s%wetu==1)
             s%uu = s%hu/s%cf/max(vmaguold,0.000001d0) &
                  *(-par%g*s%dzsdx+par%rhoa*par%Cd*s%windsu**2/(par%rho*s%hum))
          elsewhere
             s%uu = 0.d0
          endwhere
          ! update vmageu
          if (s%ny>0) then
             s%vu(1:s%nx,2:s%ny)= 0.25d0*(s%vv(1:s%nx,1:s%ny-1)+s%vv(1:s%nx,2:s%ny)+ &
                  s%vv(2:s%nx+1,1:s%ny-1)+s%vv(2:s%nx+1,2:s%ny))
             if(xmpi_isleft) then
                s%vu(:,1) = s%vu(:,2)
             endif
             if(xmpi_isright) then
                s%vu(:,s%ny+1) = s%vu(:,s%ny)
             endif
          else 
             s%vu(1:s%nx,1)= 0.5d0*(s%vv(1:s%nx,1)+s%vv(2:s%nx+1,1))
          endif !s%ny>0
          s%vmageu = sqrt(s%uu**2+s%vu**2)
          !
          ! Check residual error
          flowerr = max(maxval(abs(s%vmagev-vmagvold)),maxval(abs(s%vmageu-vmaguold)))          
       enddo
       ! calculate all derivatives
       s%ueu=s%uu
       s%vev=s%vv
       s%qx=s%uu*s%hu
       s%qy=s%vv*s%hv
       s%vmagu=s%vmageu
       s%vmagv=s%vmagev
       s%u(2:s%nx,:)=0.5d0*(s%uu(1:s%nx-1,:)+s%uu(2:s%nx,:))
       if(xmpi_istop) then
          s%u(1,:)=s%uu(1,:)
       endif
       if(xmpi_isbot) then
          s%u(s%nx+1,:)=s%u(s%nx,:)
       endif
       if (s%ny>0) then
          s%v(:,2:s%ny)=0.5d0*(s%vv(:,1:s%ny-1)+s%vv(:,2:s%ny))
          if(xmpi_isleft) then
             s%v(:,1)=s%vv(:,1)
          endif
          if(xmpi_isright) then
             s%v(:,s%ny+1)=s%v(:,s%ny)
          endif
          s%v(s%nx+1,:)=s%v(s%nx,:)
       else ! Dano
          s%v=s%vv
       endif !s%ny>0
    endif
    !
    ! Initialize for tide instant boundary condition
    !
    if (par%tidetype==TIDETYPE_INSTANT) then
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

          if (indbay==1 .and. indoff==s%nx+1) then ! in case of completely wet arrays linear interpolation for s%zs0fac
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

    type(spacepars),target              :: s
    type(parameters)                    :: par

    integer                             :: i,j,m,jg,start,ier
    character(slen)                     :: fnameg
    character(len=4)                    :: tempc
    real*8                              :: tempr


    if(.not. xmaster) return

    allocate(s%ccg(1:s%nx+1,1:s%ny+1,par%ngd))
    allocate(s%dcbdy(1:s%nx+1,1:s%ny+1))
    allocate(s%dcbdx(1:s%nx+1,1:s%ny+1))
    allocate(s%dcsdy(1:s%nx+1,1:s%ny+1))
    allocate(s%dcsdx(1:s%nx+1,1:s%ny+1))
    allocate(s%Tsg(1:s%nx+1,1:s%ny+1,par%ngd)) 
    allocate(s%Susg(1:s%nx+1,1:s%ny+1,par%ngd)) 
    allocate(s%Svsg(1:s%nx+1,1:s%ny+1,par%ngd)) 
    allocate(s%Subg(1:s%nx+1,1:s%ny+1,par%ngd)) 
    allocate(s%Svbg(1:s%nx+1,1:s%ny+1,par%ngd)) 
    allocate(s%vmag(1:s%nx+1,1:s%ny+1)) 
    allocate(s%ceqsg(1:s%nx+1,1:s%ny+1,par%ngd))
    allocate(s%ceqbg(1:s%nx+1,1:s%ny+1,par%ngd))
    !if (par%dilatancy==1) allocate(s%D15(1:par%ngd)) ! Lodewijk
    !if (par%dilatancy==1) then
    !  allocate(s%D15(1:par%ngd)) ! Lodewijk
    !else
    !  allocate(s%D15(1)) ! wwvv: s%D15 will be distributed, so it must exist
    !endif
    allocate(s%D15(1:par%ngd)) ! Robert: this is distributed according to stencil in spaceparams.tmpl
    allocate(s%D50(1:par%ngd))
    allocate(s%D90(1:par%ngd))
    allocate(s%D50top(1:s%nx+1,1:s%ny+1))
    allocate(s%D90top(1:s%nx+1,1:s%ny+1))
    allocate(s%sedcal(1:par%ngd))
    allocate(s%ucrcal(1:par%ngd))
    allocate(s%nd(1:s%nx+1,1:s%ny+1))
    allocate(s%dzbed(1:s%nx+1,1:s%ny+1,1:max(par%nd,3))) 
    allocate(s%pbbed(1:s%nx+1,1:s%ny+1,1:max(par%nd,3),1:par%ngd)) 
    allocate(s%z0bed(1:s%nx+1,1:s%ny+1))
    allocate(s%ureps(1:s%nx+1,1:s%ny+1))
    allocate(s%urepb(1:s%nx+1,1:s%ny+1))
    allocate(s%vreps(1:s%nx+1,1:s%ny+1))
    allocate(s%vrepb(1:s%nx+1,1:s%ny+1))
    allocate(s%ero(1:s%nx+1,1:s%ny+1,1:par%ngd))
    allocate(s%depo_ex(1:s%nx+1,1:s%ny+1,1:par%ngd))
    allocate(s%depo_im(1:s%nx+1,1:s%ny+1,1:par%ngd))
    allocate(s%kb(1:s%nx+1,1:s%ny+1))
    allocate(s%Tbore(1:s%nx+1,1:s%ny+1))
    allocate(s%ua(1:s%nx+1,1:s%ny+1))  
    allocate(s%dzav(1:s%nx+1,1:s%ny+1))  
    allocate(s%Sk(1:s%nx+1,1:s%ny+1))
    allocate(s%As(1:s%nx+1,1:s%ny+1))
    allocate(s%kturb(1:s%nx+1,1:s%ny+1))
    allocate(s%rolthick(1:s%nx+1,1:s%ny+1))
    allocate(s%Sutot(1:s%nx+1,1:s%ny+1))     ! Only really for easy output 
    allocate(s%Svtot(1:s%nx+1,1:s%ny+1))     ! Only really for easy output
    allocate(s%cctot(1:s%nx+1,1:s%ny+1))     ! Only really for easy output
    allocate(s%runup(1:s%ny+1))
    allocate(s%Hrunup(1:s%ny+1))
    allocate(s%xHrunup(1:s%ny+1))
    allocate(s%istruct(1:s%ny+1))
    allocate(s%iwl(1:s%ny+1))
    allocate(s%strucslope(1:s%ny+1))
    allocate(s%Dc(1:s%nx+1,1:s%ny+1))

    ! Initialize so structures can be implemented more easily
    s%pbbed = 0.d0
    !
    ! Set grain size(s)
    !
    if (par%dilatancy==1) s%D15 = par%D15(1:par%ngd)
    s%D50 = par%D50(1:par%ngd)
    s%D90 = par%D90(1:par%ngd)
    s%sedcal = par%sedcal(1:par%ngd)
    s%ucrcal = par%ucrcal(1:par%ngd)


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
             do j=1,s%ny+1
                read(31,*,iostat=ier)(s%pbbed(i,j,m,jg),i=1,s%nx+1)
                if (ier .ne. 0) then
                   call report_file_read_error(fnameg)
                endif
             enddo
          enddo
          close(31)
       enddo
       ! Rework pbbed so that sum fractions = 1
       do m=1,par%nd
          do j=1,s%ny+1     !Jaap instead of 2:s%ny
             do i=1,s%nx+1 !Jaap instead of 2:s%nx

                tempr=sum(s%pbbed(i,j,m,1:par%ngd))
                if (abs(1.d0-tempr)>0.d0) then
                   ! Maybe fix this warning if in combination with structures
                   call writelog('lws','ai0ai0ai0a',' Warning: Resetting sum of sediment fractions in point (',&
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
    do j=1,s%ny+1
       do i=1,s%nx+1
          s%D50top(i,j) =  sum(s%pbbed(i,j,1,:)*s%D50)
          s%D90top(i,j) =  sum(s%pbbed(i,j,1,:)*s%D90)
       enddo
    enddo
    ! 
    ! Set non-erodable layer
    !
    allocate(s%structdepth(s%nx+1,s%ny+1))

    s%structdepth = 100.d0

    if (par%struct==1) then
       !call readkey('params.txt','ne_layer',fnameh)
       !open(31,file=fnameh)
       open(31,file=par%ne_layer)

       do j=1,s%ny+1
          read(31,*,iostat=ier)(s%structdepth(i,j),i=1,s%nx+1)
          if (ier .ne. 0) then
             call report_file_read_error(par%ne_layer)
          endif
       end do

       close(31)

    endif

    ! bottom of sediment model
    s%z0bed = s%zb - sum(s%dzbed,DIM=3)

    s%nd = max(par%nd,2)

    s%ureps      = 0.d0
    s%vreps      = 0.d0
    s%ccg        = 0.d0
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
    s%runup      = 0.d0
    s%Hrunup     = 0.d0
    s%xHrunup    = 0.d0
    s%istruct    = s%nx+1
    s%strucslope = 0.d0
    s%Dc         = 0.d0

    ! Initialize dzbdx, dzbdy
    do j=1,s%ny+1
       do i=1,s%nx
          s%dzbdx(i,j)=(s%zb(i+1,j)-s%zb(i,j))/s%dsu(i,j)
       enddo
    enddo
    ! dummy, needed to keep compiler happy
    s%dzbdx(s%nx+1,:)=s%dzbdx(s%nx,:)

    if (s%ny>0) then
       do j=1,s%ny
          do i=1,s%nx+1
             s%dzbdy(i,j)=(s%zb(i,j+1)-s%zb(i,j))/s%dnv(i,j)
          enddo
       enddo
       s%dzbdy(:,s%ny+1)=s%dzbdy(:,s%ny)
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

    type(spacepars),target                  :: s
    type(parameters)                        :: par

    integer                                 :: i,j
    integer                                 :: io
    integer                                 :: m1,m2,n1,n2
    real*8, dimension(:),allocatable        :: xdb,ydb,xde,yde
    integer,dimension(2)                    :: mnb,mne


    if(.not. xmaster) return

    io          = 0

    allocate(xdb            (par%ndischarge)      )
    allocate(ydb            (par%ndischarge)      )
    allocate(xde            (par%ndischarge)      )
    allocate(yde            (par%ndischarge)      )

    allocate(s%pntdisch     (1:par%ndischarge)                      )
    allocate(s%pdisch       (1:par%ndischarge   , 1:4)              )
    allocate(s%tdisch       (1:par%ntdischarge)                     )
    allocate(s%qdisch       (1:par%ntdischarge  , 1:par%ndischarge) )

    s%pntdisch  = 0
    s%tdisch    = 0.d0
    s%pdisch    = 0
    s%qdisch    = 0.d0

    !if (xmaster) then
    if(.true.) then
       if (par%ndischarge>0) then

          ! read discharge locations
          open(10,file=par%disch_loc_file)
          do i=1,par%ndischarge
             read(10,*,IOSTAT=io) xdb(i),ydb(i),xde(i),yde(i)
             if (io .ne. 0) then
                call report_file_read_error(par%disch_loc_file)
             endif
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
                if (io .ne. 0) then
                   call report_file_read_error(par%disch_timeseries_file)
                endif
             enddo
             close(10)
          endif
       endif
    endif

    !if (xmaster) then
    if(.true.) then

       ! initialise each discharge location
       do i=1,par%ndischarge

!          dxd = abs(xde(i)-xdb(i))
!          dyd = abs(yde(i)-ydb(i))
          mnb = minloc(sqrt((s%xz-xdb(i))**2+(s%yz-ydb(i))**2))
          mne = minloc(sqrt((s%xz-xde(i))**2+(s%yz-yde(i))**2))

          ! convert discharge location to cell indices depending on type of discharge:
          !     point discharge, in v-direction or in u-direction
          if (s%pntdisch(i).eq.1) then

             ! point discharge (no orientation, no added momentum, just mass)

 !            mnb = minloc(sqrt((s%xz-xdb(i))**2+(s%yz-ydb(i))**2))
 !            mne = minloc(sqrt((s%xz-xde(i))**2+(s%yz-yde(i))**2))

             s%pdisch(i,:) = (/mnb(1),mnb(2),0,0/)
!          elseif (dxd.gt.dyd) then
          elseif (mnb(1).ne.mne(1)) then

             ! discharge through v-points

!             mnb = minloc(sqrt((s%xv-xdb(i))**2+(s%yv-ydb(i))**2))
!             mne = minloc(sqrt((s%xv-xde(i))**2+(s%yv-yde(i))**2))

             m1 = minval((/mnb(1),mne(1)/))
             m2 = maxval((/mnb(1),mne(1)/))
             n1 = nint(0.5*(mnb(2)+mne(2)))

             if (n1.lt.1)    n1 = 1
             if (n1.gt.s%ny) n1 = s%ny

             s%pdisch(i,:) = (/m1,n1,m2,n1/)
          else

             ! discharge through u-points

!             mnb = minloc(sqrt((s%xu-xdb(i))**2+(s%yu-ydb(i))**2))
!             mne = minloc(sqrt((s%xu-xde(i))**2+(s%yu-yde(i))**2))

             m1 = nint(0.5*(mnb(1)+mne(1)))
             n1 = minval((/mnb(2),mne(2)/))
             n2 = maxval((/mnb(2),mne(2)/))

             if (m1.lt.1)    m1 = 1
             if (m1.gt.s%nx) m1 = s%nx

             s%pdisch(i,:) = (/m1,n1,m1,n2/)
          endif
       enddo

       ! incorporate morfac
       if (par%morfacopt == 1) then
          s%tdisch = s%tdisch/max(par%morfac,1.d0)
       endif
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

    type(spacepars),target                  :: s
    type(parameters)                        :: par

    integer                                 :: i,ier
    real*8                                  :: xdrift,ydrift
    real*8                                  :: ds,dn
    integer,dimension(2)                    :: mn


    if(.not. (xmaster .or. xomaster)) return

    allocate(s%idrift   (par%ndrifter))
    allocate(s%jdrift   (par%ndrifter))
    allocate(s%tdriftb  (par%ndrifter))
    allocate(s%tdrifte  (par%ndrifter))

    if (par%ndrifter>0) then

          ! read drifter file
       if(xmaster) then
          open(10,file=par%drifterfile)
          do i=1,par%ndrifter
             read(10,*,iostat=ier)xdrift,ydrift,s%tdriftb(i),s%tdrifte(i)
             if (ier .ne. 0) then
                call report_file_read_error(par%drifterfile)
             endif

             mn          = minloc(sqrt((s%xz-xdrift)**2+(s%yz-ydrift)**2))

             ds          =  (xdrift - s%xz(mn(1),mn(2)))*cos(s%alfaz(mn(1),mn(2))) &
                  +(ydrift - s%yz(mn(1),mn(2)))*sin(s%alfaz(mn(1),mn(2)))
             dn          = -(xdrift - s%xz(mn(1),mn(2)))*sin(s%alfaz(mn(1),mn(2))) &
                  +(ydrift - s%yz(mn(1),mn(2)))*cos(s%alfaz(mn(1),mn(2)))

             s%idrift(i) = mn(1) + ds/s%dsu(mn(1),mn(2))
             s%jdrift(i) = mn(2) + dn/s%dnv(mn(1),mn(2))
          enddo
          close(10)

          ! incorporate morfac
          if (par%morfacopt == 1) then
             s%tdriftb   = s%tdriftb/max(par%morfac,1.d0)
             s%tdrifte   = s%tdrifte/max(par%morfac,1.d0)
          endif
       endif

#ifdef USEMPI
       call xmpi_send(xmpi_imaster, xmpi_omaster,s%idrift)
       call xmpi_send(xmpi_imaster, xmpi_omaster,s%jdrift)
       call xmpi_send(xmpi_imaster, xmpi_omaster,s%tdriftb)
       call xmpi_send(xmpi_imaster, xmpi_omaster,s%tdrifte)
#endif

    endif
  end subroutine drifter_init

end module initialize
