module beachwizard_module
   implicit none
   save
  type beachwiz

     real*8, dimension(1000,8)          :: timesnew          !time in minutes when remote sensed sources apply - eight sources max
     integer, dimension(8)                       :: ntimesnew = 0     !number of maps per source
     integer                                         :: nflagsigpr = 0    !flags
     integer                                         :: nflagsetobs = 0
     real*8, dimension(1000,8)          :: timvalnew         !duration in minutes of each map
     character(232), dimension(1000,8)  :: fname             !filename of the map
     character(232), dimension(1000,8)  :: fnameerr          !filename of the uncertainties in the map
     character(232), dimension(:), &
          &  pointer            :: fnamim
     character(232), dimension(:), &
          &  pointer            :: fnamzb
     character(232), dimension(:), &
          &  pointer            :: fnamcx
     character(232), dimension(:), &
          &  pointer            :: fnamsh
     integer, dimension(8)                       :: status = 0        !status of the map (1=read, 0=not read/does not exist)
     integer, dimension(8)                       :: statuserr = 0     !status of the error map
     integer, dimension(8)              :: itimnew           !counter of the time "timesnew"
     integer, dimension(8)              :: itimrdnew         !flag to indicate time has been read
     integer                            :: itim
     integer                            :: itimrd
     integer                            :: itim2
     integer                            :: itimrd2
     integer                            :: itim3
     integer                            :: itimrd3
     integer                            :: itim4
     integer                            :: itimrd4
     real*8, dimension(:,:), pointer    :: fobs1             !value of the observed quantity in the coordinate system of the file
     real*8, dimension(:,:), pointer    :: fobs1err          !same, for error file
     real*8, dimension(:,:), pointer    :: disobs 
     real*8, dimension(:,:), pointer    :: zbobs1
     real*8, dimension(:,:), pointer    :: shobs1
     !     real*8,dimension(:,:),allocatable  :: dobs              !observed    dissipation
     real*8,dimension(:,:),allocatable  :: fobs
     real*8,dimension(:,:),allocatable  :: cobs              !observed  celerity
     real*8,dimension(:,:),allocatable  :: ccom              !computed  celerity
     real*8,dimension(:,:),allocatable  :: kobs              !observed  wave number
     !     real*8,dimension(:,:),allocatable  :: shobs             !observed    shoreline
     !     real*8,dimension(:,:),allocatable  :: zbobs             !observed bathymetry
     real*8,dimension(:,:),allocatable  :: dcmdo             !computed dissipation minus observed dissipation
     real*8,dimension(:,:),allocatable  :: ccmco             !computed celerity minus observed celerity
     real*8,dimension(:,:),allocatable  :: scmso             !computed shoreline minus observed shoreline
     real*8,dimension(:,:),allocatable  :: dassim        
     real*8,dimension(:,:),allocatable  :: sigD              !standard dev. of dissipation
     real*8,dimension(:,:),allocatable  :: sigC              !standard dev. of celerity
     real*8,dimension(:,:),allocatable  :: sigK              !standard dev. of wave number
     real*8,dimension(:,:),allocatable  :: sigZ              !standard dev. of bathymetry       
     real*8,dimension(:,:),allocatable  :: sigS              !standard dev. of shoreline        
     !    real*8,dimension(:,:),allocatable  :: sig2prior
     real*8                             :: xll               !origin of the coordinate system used in the map
     real*8                             :: yll               !idem
     real*8                             :: dx                !dx in that c.s.
     real*8                             :: dy                !dy in that c.s.
     real*8                             :: angle             !angle in degrees
     real*8                             :: tradar            !period used in the wavenumber or celerity file
     integer                            :: nx                !number of points in the x-direction in the c.s. of the map
     integer                            :: ny                !same in y
     integer                                         :: ix                !counter
     integer                                         :: iy                !counter
     real*8                                          :: sigDdef           !default value of the dissipation uncertainty
     real*8                                          :: sigCdef           !default value of the celerity uncertainty
     real*8                                          :: sigZdef           !default value of the bathymetry uncertainty
     real*8                                          :: sigSdef           !default value of the shoreline uncertainty


  end type beachwiz

  ! Type of beachwizard, stored local. Please store persistent info in s, par.
   type(beachwiz) :: bw

contains      

  subroutine bwinit(s)
    use params
    use spaceparams
    use xmpi_module

      implicit none

    type(spacepars), target             :: s 

      if(.not. xmaster) return

      if(.true.) then
       allocate(s%dobs(1:s%nx+1,1:s%ny+1))    
       allocate(s%zbobs(1:s%nx+1,1:s%ny+1)) 
       allocate(s%sig2prior(1:s%nx+1,1:s%ny+1))
       allocate(s%shobs(1:s%nx+1,1:s%ny+1))
       allocate(s%bwalpha(1:s%nx+1,1:s%ny+1))
       allocate(s%dcmdo(1:s%nx+1,1:s%ny+1))
       allocate(s%dassim(1:s%nx+1,1:s%ny+1))
       allocate(s%cobs(1:s%nx+1,1:s%ny+1))
      endif


  end subroutine bwinit

  subroutine assim(s,par)

    use params
    use spaceparams



      implicit none

    type(spacepars), target             :: s 
    type(parameters)                    :: par

    character*80                        :: infofile
    integer                                             :: im,in
    integer                                             :: srcnr
    real*8,  external                   :: readkey_dbl
    integer, external                   :: readkey_int
    logical                                             :: exists
    real*8, allocatable, dimension(:,:) :: h1

    integer                             :: jj
    real*8                              :: om
    real*8, allocatable, dimension(:,:) :: ome2
    real*8, allocatable, dimension(:,:) :: num
    real*8, allocatable, dimension(:,:) :: den
    real*8, allocatable, dimension(:,:) :: rk
    real*8, allocatable, dimension(:,:) :: Hrms
    real*8, allocatable, dimension(:,:) :: E
    real*8, allocatable, dimension(:,:) :: hh

    real*8, parameter                   :: a1 = 5.060219360721177D-01
    real*8, parameter                   :: a2 = 2.663457535068147D-01
    real*8, parameter                   :: a3 = 1.108728659243231D-01
    real*8, parameter                   :: a4 = 4.197392043833136D-02
    real*8, parameter                   :: a5 = 8.670877524768146D-03
    real*8, parameter                   :: a6 = 4.890806291366061D-03
    real*8, parameter                   :: b1 = 1.727544632667079D-01
    real*8, parameter                   :: b2 = 1.191224998569728D-01
    real*8, parameter                   :: b3 = 4.165097693766726D-02
    real*8, parameter                   :: b4 = 8.674993032204639D-03 



    if (.not. allocated(bw%fobs)) then
       !          allocate(s%dobs(1:s%nx+1,1:s%ny+1))
       allocate(bw%fobs(1:s%nx+1,1:s%ny+1))
       allocate(bw%cobs(1:s%nx+1,1:s%ny+1))
       allocate(bw%ccom(1:s%nx+1,1:s%ny+1))
       allocate(bw%kobs(1:s%nx+1,1:s%ny+1))
       !          allocate(s%shobs(1:s%nx+1,1:s%ny+1))
       !          allocate(s%zbobs(1:s%nx+1,1:s%ny+1))
       allocate(bw%dcmdo(1:s%nx+1,1:s%ny+1))
       allocate(bw%ccmco(1:s%nx+1,1:s%ny+1))
       allocate(bw%scmso(1:s%nx+1,1:s%ny+1))      
       allocate(bw%dassim(1:s%nx+1,1:s%ny+1))
       allocate(bw%sigD(1:s%nx+1,1:s%ny+1))
       allocate(bw%sigC(1:s%nx+1,1:s%ny+1))
       allocate(bw%sigK(1:s%nx+1,1:s%ny+1))      
       allocate(bw%sigZ(1:s%nx+1,1:s%ny+1))
       allocate(bw%sigS(1:s%nx+1,1:s%ny+1))
       !          allocate(s%sig2prior(1:s%nx+1,1:s%ny+1))
    end if

    allocate(h1(1:s%nx+1,1:s%ny+1))
    allocate(hh(1:s%nx+1,1:s%ny+1))
    allocate(ome2(1:s%nx+1,1:s%ny+1))
    allocate(num(1:s%nx+1,1:s%ny+1))
    allocate(rk(1:s%nx+1,1:s%ny+1))
    allocate(E(1:s%nx+1,1:s%ny+1))
    allocate(Hrms(1:s%nx+1,1:s%ny+1))     

    ! wwvv somebody forgot:
    allocate(den(s%nx+1,s%ny+1))

    bw%dcmdo=0.
    bw%ccmco=0.
    bw%scmso=0.

    if (bw%nflagsetobs==0) then !do this only once.

       s%dobs=-999.
       bw%cobs=-999.
       s%shobs=-999.
       s%zbobs = -999.

       bw%nflagsetobs = 1

    end if

    if (bw%nflagsigpr==0) then  !read prior uncertainty in terms of variance

       inquire (file='sigpr.dep', exist=exists)
       if  (exists) then
          open(444,file='sigpr.dep',err=999)
          do im = 1,s%ny+1
             read(444,*)(s%sig2prior(in,im),in =1,s%nx+1)  !Leo: nx and ny --> nx+1 ny+1
          enddo
          close(444)
       else
          s%sig2prior = 1**2 !if the file does not exist, set to hardcoded 1**2
       end if

       inquire (file='bwdefaults.inp', exist=exists) !see if defaults file exists
       if  (exists) then
          open(444,file='bwdefaults.inp',err=999)
          read(444,*) bw%sigDdef
          read(444,*) bw%sigCdef 
          read(444,*) bw%sigZdef
          read(444,*) bw%sigSdef
          close(444)
          bw%sigD=bw%sigDdef
          bw%sigC=bw%sigCdef
          bw%sigZ=bw%sigZdef
          bw%sigS=bw%sigSdef
       else
          bw%sigCdef=1. !set to hardcoded defaults
          bw%sigDdef=20.
          bw%sigZdef=0.5
          bw%sigSdef=0.5
       end if
       bw%nflagsigpr=1

    end if

999 continue


    inquire (file='imageinfoimap',exist=exists) !read file for dissipation 
    if (exists) then
       infofile='imageinfoimap'
       srcnr = 1
       call assim_rd(bw, s, par, infofile,   srcnr)
       if (bw%status(srcnr) ==1) then 
          do in=1,s%nx
             do im=1,s%ny
                if ((s%wetu(in,im)>0).and.(s%dobs(in,im)>-5.)) then
                   bw%dcmdo(in,im) = s%Dr(in,im) - s%dobs(in,im)
                end if
             end do
          end do
       end if
       if (bw%statuserr(srcnr) ==0) then 
          h1 = max(s%hh+par%delta*s%H,par%hmin)
         !! bw%sigD=bw%sigDdef+(999.-bw%sigDdef)*(0.5+0.5*tanh((s%dobs/h1-bw%sigDdef)/10.)) ! change 50 --> bw%sigDdef 
          bw%sigD = bw%sigDdef               !!using only default dissipation uncertainty for now...     
     
       end if
    end if

    inquire (file='imageinfozb',exist=exists) !read file for bathymetry
    if (exists) then
       infofile='imageinfozb'
       srcnr = 2
       call assim_rd(bw, s, par, infofile,   srcnr)
       if (bw%statuserr(srcnr) ==0) then 
          bw%sigZ=bw%sigZDef
       end if
    end if

    inquire (file='imageinfocx',exist=exists) !read file for celerity
    if (exists) then
       infofile='imageinfocx'
       srcnr = 3
       call assim_rd(bw, s, par, infofile,   srcnr)
       if (bw%status(srcnr) ==1) then
          do in=1,s%nx+1
             do im=1,s%ny+1
                om = 2.*3.1415/bw%tradar             
                ome2(in,im) = om**2*s%hh(in,im)/par%g
                num(in,im) = 1.0 +   &
                     & ome2(in,im)*(a1 + ome2(in,im)*(a2 + ome2(in,im)*(a3 + ome2(in,im)   &
                     & *(a4 + ome2(in,im)*(a5 + ome2(in,im)*a6)))))
                den(in,im) = 1.0 + ome2(in,im)*(b1 + ome2(in,im)*(b2 + ome2(in,im)*(b3 + &
                     & ome2(in,im)*(b4 + ome2(in,im)*a6))))
                rk(in,im) = sqrt(ome2(in,im)*num(in,im)/den(in,im))/s%hh(in,im)
                Hrms(in,im) = max(0.01d0,sqrt(8*max(0.01d0,s%E(in,im))/par%rho/par%g))
                bw%ccom(in,im) = sqrt(par%g/rk(in,im)*tanh(rk(in,im)*s%hh(in,im)+Hrms(in,im)))
                if ((s%wetu(in,im)>0).and.(bw%cobs(in,im)>-990.)) then
                   bw%ccmco(in,im) = bw%ccom(in,im) - bw%cobs(in,im)
                end if
             enddo
          enddo
       end if
       if (bw%statuserr(srcnr) ==0) then 
          bw%sigC=bw%sigCdef
       endif
    end if

    inquire (file='imageinfokabs',exist=exists) !do the same as above for wavenumbers
    if (exists) then
       infofile='imageinfokabs'
       srcnr = 4
       call assim_rd(bw, s, par, infofile,   srcnr)
       if (bw%status(srcnr) ==1) then 
          do in = 1,s%nx+1
             do im = 1,s%ny+1
                om = 2.*3.1415/bw%tradar             
                ome2(in,im) = om**2*s%hh(in,im)/par%g
                num(in,im) = 1.0 +   &
                     & ome2(in,im)*(a1 + ome2(in,im)*(a2 + ome2(in,im)*(a3 + ome2(in,im)   &
                     & *(a4 + ome2(in,im)*(a5 + ome2(in,im)*a6)))))
                den(in,im) = 1.0 + ome2(in,im)*(b1 + ome2(in,im)*(b2 + ome2(in,im)*(b3 + &
                     & ome2(in,im)*(b4 + ome2(in,im)*a6))))
                rk(in,im) = sqrt(ome2(in,im)*num(in,im)/den(in,im))/s%hh(in,im)
                Hrms(in,im) = max(0.01d0,sqrt(8*max(0.01d0,s%E(in,im))/par%rho/par%g))
                bw%ccom(in,im) = sqrt(par%g/rk(in,im)*tanh(rk(in,im)*s%hh(in,im)+Hrms(in,im)))
                bw%cobs(in,im) = om/max(bw%kobs(in,im),0.01d0)    ! here we revert back to celerities.   
                if (bw%kobs(in,im)<-990.) then
                   bw%cobs(in,im)=-999.                
                end if
                if ((s%wetu(in,im)>0).and.(bw%kobs(in,im)>-990.)) then
                   bw%ccmco(in,im) = bw%ccom(in,im) - bw%cobs(in,im)
                end if
             end do
          end do
       end if
       if (bw%statuserr(srcnr) ==1) then
          do in = 1,s%nx+1
             do im = 1,s%ny+1
                om = 2.*3.1415/bw%tradar             
                ome2(in,im) = om**2*s%hh(in,im)/par%g
                num(in,im) = 1.0 +   &
                     & ome2(in,im)*(a1 + ome2(in,im)*(a2 + ome2(in,im)*(a3 + ome2(in,im)   &
                     & *(a4 + ome2(in,im)*(a5 + ome2(in,im)*a6)))))
                den(in,im) = 1.0 + ome2(in,im)*(b1 + ome2(in,im)*(b2 + ome2(in,im)*(b3 + &
                     & ome2(in,im)*(b4 + ome2(in,im)*a6))))
                rk(in,im) = sqrt(ome2(in,im)*num(in,im)/den(in,im))/s%hh(in,im)
                bw%sigC(in,im) = om/rk(in,im)*(bw%sigK(in,im)/(rk(in,im)+bw%sigK(in,im))) ! express sigK in terms of sigK, 
                ! from c+sigC = om/(k+sigK)
             end do
          end do
       end if
       if (bw%statuserr(srcnr) ==0) then 
          bw%sigC = bw%sigCdef 
       end if
    end if

    inquire (file='imageinfoibathy',exist=exists) !read file for shorelines
    if (exists) then
       infofile='imageinfoibathy'
       srcnr = 5
       call assim_rd(bw, s, par, infofile,   srcnr)
       if (bw%status(srcnr) ==1) then

          do im=1,s%ny
             do in = 1,s%nx
                if (abs(s%shobs(in,im)).le.5.) then
            !                bw%scmso(in,im) = s%hh(in,im) + s%shobs(in,im)  !!plus?
                   bw%scmso(in,im) = s%zb(in,im) - s%shobs(in,im) 
                end if
             end do
          end do


          do jj=1,2
             do in = 2,s%nx
                do im = 2,s%ny
                   bw%scmso(in,im) = (bw%scmso(in-1,im-1)+bw%scmso(in,im-1)+bw%scmso(in+1,im-1)+ &
                        &          bw%scmso(in-1,im)+  bw%scmso(in,im)+  bw%scmso(in+1,im)+     &
                        &          bw%scmso(in-1,im+1)+bw%scmso(in,im+1)+bw%scmso(in+1,im+1))/9. 
                end do
             end do
          end do

       end if !if status
       if (bw%statuserr(srcnr) ==0) then 
          bw%sigS=bw%sigSdef
       end if
    end if !if exists
    
    !expose variables to output module
    
    s%dcmdo = bw%dcmdo

    call comp_depchg(bw, s, par)

  end subroutine assim

  subroutine assim_rd(bw, s, par, infofile,   srcnr)
  
    

    use params
    use spaceparams
    implicit none
    !     Global variables
    type(spacepars), target             :: s 
    type(parameters)                    :: par
    type(beachwiz)                                              :: bw
    character*80                        :: infofile
    integer                                                             :: im,in,iy,ix
    integer                                                             :: srcnr

    !     Local variables
    real*8                           :: degrad
    real*8                           :: a
    real*8                           :: b
    real*8                           :: w11
    real*8                           :: w21
    real*8                           :: w12
    real*8                           :: w22
    real*8                           :: x1
    real*8                           :: y1
    real*8                           :: xs
    real*8                           :: ys
    real*8                           :: cs
    real*8                           :: sn
    real*8                           :: x1max
    real*8                           :: y1max
    integer                            :: ixp1
    integer                            :: iyp1
    integer                            :: i
    integer                            :: istat
    logical                            :: exists

    !        Actual time based on timhr relative to itdate
    !        


    if (bw%ntimesnew(srcnr)==0) then
       inquire (file=infofile,exist=exists)


       if (exists) then

          open(31,file=infofile,err=999)
          read(31,*)bw%ntimesnew(srcnr)

          do i=1,bw%ntimesnew(srcnr)
             read(31,*)bw%timesnew(i,srcnr),bw%fname(i,srcnr),bw%fnameerr(i,srcnr), &
                  & bw%timvalnew(i,srcnr)
             write(*,*) bw%timesnew(i,srcnr), bw%fname(i,srcnr), bw%fnameerr(i,srcnr), &
                  & bw%timvalnew(i,srcnr)
          enddo
          !pause
          close(31)
          bw%timesnew(bw%ntimesnew(srcnr)+1:1000,srcnr)=1.e10


       else
          stop ' Argus assimilation active and no imageinfoimap file'
       endif
       bw%itimnew(srcnr)=1
       bw%itimrdnew(srcnr)=0

    endif

    bw%status(srcnr) = 0


    !        Check if time matches time of an image
    if (par%t<60.*bw%timesnew(1,srcnr)) then
       !
       !           Time less than time first image; leave dcmdo at 0
    else
10     if (par%t>=60.*bw%timesnew(bw%itimnew(srcnr)+1,srcnr)) then
          bw%itimnew(srcnr)=bw%itimnew(srcnr)+1
          goto 10
       endif
       !

       if (par%t-60.*bw%timesnew(bw%itimnew(srcnr),srcnr)<=60*bw%timvalnew(bw%itimnew(srcnr),srcnr)) then
          !              Read image if not read before
          !               write(*,*) 'smaller than timval',timmin,timesnew(itimnew(srcnr),srcnr), &
          !&                       timvalnew(itimnew(srcnr),srcnr)             
          if (bw%itimnew(srcnr)/=bw%itimrdnew(srcnr)) then
             !
             !                 Read header file that defines argus grid
             inquire(file=bw%fname(bw%itimnew(srcnr),srcnr),exist=exists)
             if (exists) then
                open(31,file=bw%fname(bw%itimnew(srcnr),srcnr),err=999)
                read(31,*,err=999)bw%xll
                read(31,*,err=999)bw%yll
                read(31,*,err=999)bw%dx
                read(31,*,err=999)bw%dy
                read(31,*,err=999)bw%nx
                read(31,*,err=999)bw%ny
                read(31,*,err=999)bw%angle 

                if (associated(bw%fobs1)) deallocate (bw%fobs1, STAT = istat)
                allocate (bw%fobs1(bw%nx,bw%ny))
                do iy=1,bw%ny
                   read(31,*)(bw%fobs1(ix,iy),ix=1,bw%nx)
                enddo
                close(31)
                bw%itimrdnew(srcnr)=bw%itimnew(srcnr)
             else
                goto 999
             endif
             !    inquire(file=fnameerr(itimnew(srcnr),srcnr),exist=exists)
             !    write(*,*) fnameerr(itimnew(srcnr),srcnr)
             !pause
             !     if (exists) then
             !       31=NEWLUN(GDP) 
             !        open(31,file=fnameerr(itimnew(srcnr),srcnr),err=999)
             !        read(31,*,err=999)xll
             !        read(31,*,err=999)yll
             !        read(31,*,err=999)dx
             !        read(31,*,err=999)dy
             !        read(31,*,err=999)nx
             !        read(31,*,err=999)ny
             !        read(31,*,err=999)angle 
             !                   write(*,*) xll, yll, dx, dy, nx, ny, angle
             !pause
             !            if (associated(fobs1err)) deallocate (fobs1err, STAT = istat)
             !           allocate (fobs1err(nx,ny))
             !          do iy=1,ny
             !            read(31,*)(fobs1err(ix,iy),ix=1,nx)
             !        enddo
             !       close(31)
             !           statuserr(srcnr)=1
             !           write(*,*) 'here read ', statuserr(srcnr), srcnr
             !pause
             !     else
             !       fobs1err=-999
             !           statuserr(srcnr)=0
             !           write(*,*) 'here not read ', statuserr(srcnr), srcnr
             !           !pause
             !     endif 
             degrad=atan(1.)/45.
             cs=cos(bw%angle*degrad)
             sn=sin(bw%angle*degrad)
             x1max=(bw%nx-2)*bw%dx
             y1max=(bw%ny-2)*bw%dy

             do in=1,s%nx+1                  !Leo and Roberto changed nx --> nx+1
                do im=1,s%ny+1               !Leo and Roberto changed ny --> ny+1
                   bw%fobs(in,im)=-999.d0 !Changed by Ap 31/5
                   !             fobserr(n,m) = 999.d0 ! changed by Ap 31/5
                   !   if (kfs(n,m)>0) then !Ap
                   xs = s%xz(in,im) - bw%xll     ! Jaap changed x --> xz
                   ys = s%yz(in,im) - bw%yll
                   x1 = xs*cs + ys*sn
                   y1 =-xs*sn + ys*cs
                   x1 = min(max(x1,bw%dx),x1max)
                   y1 = min(max(y1,bw%dy),y1max)
                   ix = int(x1/bw%dx)+1
                   iy = int(y1/bw%dy)+1
                   ixp1 = max(min(ix+1,bw%nx-1),2)
                   iyp1 = max(min(iy+1,bw%ny-1),2)
                   ix   = max(min(ix,bw%nx-1),2)
                   iy   = max(min(iy,bw%ny-1),2)
                   a    = mod(x1,bw%dx)/bw%dx
                   b    = mod(y1,bw%dy)/bw%dy
                   w11  = (1.-b)*(1.-a)
                   w21  = (1.-b)*    a
                   w12  =     b *(1.-a)
                   w22  =     b *    a
                   !                       if (kfs(n,m)>0) then
                   bw%fobs(in,im)      = w11*bw%fobs1(ix  ,iy  ) + &
                        &                              w21*bw%fobs1(ixp1,iy  ) + &
                        &                              w12*bw%fobs1(ix  ,iyp1) + &
                        &                              w22*bw%fobs1(ixp1,iyp1)

                   !                       endif
                   !                                            if (statuserr(srcnr)==1) then
                   !
                   !                        fobserr(n,m)      = w11*fobs1err(ix  ,iy  ) + &
                   !     &                              w21*fobs1err(ixp1,iy  ) + &
                   !     &                              w12*fobs1err(ix  ,iyp1) + &
                   !     &                              w22*fobs1err(ixp1,iyp1)
                   !                                    end if

                enddo
             enddo
             if (srcnr==1) s%dobs=bw%fobs
             if (srcnr==2) s%zbobs=bw%fobs
             if (srcnr==3) bw%cobs=bw%fobs
             if (srcnr==4) bw%cobs=bw%fobs
             if (srcnr==5) s%shobs=bw%fobs

          endif

          bw%status(srcnr)=1

       else
          bw%fobs=0.
       endif
    endif

999 continue

  end subroutine assim_rd

  subroutine assim_update(s, par)
    use params
    use spaceparams
    type(spacepars), target             :: s 
    type(parameters)                    :: par

    ! These are doing the same thing now
    if (par%bchwiz .eq. 1) then
       ! s%dzbdt(i,j) = s%dzbdt(i,j)-bw%dassim(i,j)/par%dt
       ! Times dt*morfac gives 
       s%zb = s%zb - bw%dassim
    endif
    if (par%bchwiz .eq. 2) then
       s%zb = s%zb - bw%dassim  !if we have only beachwizard
    end if
  end subroutine assim_update
  subroutine comp_depchg(bw, s, par)
    
    use params
    use spaceparams
    implicit none
    !     Global variables
    type(spacepars), target             :: s 
    type(parameters)                    :: par
    type(beachwiz)                                              :: bw


    !     Local variables
    integer  :: in,i
    integer  :: im,j 
    integer  :: ind
    real*8 :: Hrmsmax 
    real*8 :: mxdDdh
    real*8 :: alp_as
    real*8 :: Nassim
    real*8, allocatable, dimension(:,:)  :: errD
    real*8, allocatable, dimension(:,:)  :: errC
    real*8, allocatable, dimension(:,:)  :: errS
    real*8, allocatable, dimension(:,:)  :: Hb
    real*8, allocatable, dimension(:,:)  :: kh
    real*8, allocatable, dimension(:,:)  :: gambal
    real*8, allocatable, dimension(:,:)  :: Ga 
    real*8, allocatable, dimension(:,:)  :: dDdGa 
    real*8, allocatable, dimension(:,:)  :: dGadHb 
    real*8, allocatable, dimension(:,:)  :: dHbdh 
    real*8, allocatable, dimension(:,:)  :: h1
    real*8, allocatable, dimension(:,:)  :: dDdh
    real*8, allocatable, dimension(:,:)  :: sig2obs 
      !    real*8, allocatable, dimension(:,:)  :: alpha
    real*8, allocatable, dimension(:,:)  :: Hrms
    real*8, allocatable, dimension(:,:)  :: alphafac
    !
    real*8                              :: backdis,disfac
    integer                             :: index
    real*8, allocatable, dimension(:,:)  :: h2
    if (.not.(allocated(errD))) then
       allocate (errD(1:s%nx+1,1:s%ny+1))
       allocate (errC(1:s%nx+1,1:s%ny+1))
       allocate (errS(1:s%nx+1,1:s%ny+1))
       allocate (Hb(1:s%nx+1,1:s%ny+1))
       allocate (kh(1:s%nx+1,1:s%ny+1))
       allocate (gambal(1:s%nx+1,1:s%ny+1))
       allocate (Ga(1:s%nx+1,1:s%ny+1))
       allocate (dGadHb(1:s%nx+1,1:s%ny+1))
       allocate (dHbdh(1:s%nx+1,1:s%ny+1))
       allocate (dDdGa(1:s%nx+1,1:s%ny+1))
       allocate (h1(1:s%nx+1,1:s%ny+1))
       allocate (dDdh(1:s%nx+1,1:s%ny+1))
       allocate (sig2obs(1:s%nx+1,1:s%ny+1))
      ! allocate (alpha(1:s%nx+1,1:s%ny+1))
       allocate (Hrms(1:s%nx+1,1:s%ny+1))
       allocate (h2(1:s%nx+1,1:s%ny+1))
       allocate (alphafac(1:s%nx+1,1:s%ny+1))
    end if


    ! general
    Hrms = max(s%H,0.01d0);
    h2 = 0.d0           ! modified wave length, initially set to L1
     do j = 2,s%ny
        do i = 2,s%nx+1
           index = i       ! start index
           backdis = 0.d0  ! relative distance backward
           do while (backdis<1.d0)
              ! disfac = s%dsc(index,j)/(par%facsd*s%L1(index,j))
              ! use average wavelength over distance dsc
              disfac = s%dsc(index,j)/(0.5d0*(s%L1(index,j)+s%L1(max(index-1,1),j)))
              disfac = min(disfac,1.d0-backdis)
              
              h2(i,j) = h2(i,j) + disfac*s%hh(index,j)
              backdis = backdis+disfac
              
              index = max(index-1,1)
           enddo
        enddo
     enddo
     h2(:,1)=h2(:,2)
     h2(:,s%ny+1)=h2(:,s%ny)
    h1 = max(s%hh+par%delta*Hrms,par%hmin)     !total water depth
    kh = min(10.0, s%k*h1) 
    
                                ! ! The following equation numbers and variables follow the Beach Wizard Paper 
                                                !!       (Van Dongeren et al., 2008)

    ! derivatives for dissipation assimilation
    gambal = 0.29+0.76*kh
    Hb = 0.88/s%k*tanh(gambal*kh/0.88)          ! eq. (A4)
    Hrmsmax = maxval(Hrms)
    Ga = (Hb/Hrms)**2        

    dDdGa = -0.25*par%rho*par%g/par%Trep*(Hrms**2)*Ga*exp(-Ga)   !eq. (A5)
    dGadHb = 2*Hb/(Hrms**2)
    
    dHbdh = ((0.29+2*0.76*kh)*(1.-(kh/(sinh(kh)*cosh(kh)+kh)))) &      ! eq. (A6)
         & /(cosh((0.29*kh+0.76*kh**2)/0.88)**2) + &
         & 0.88*tanh((0.29*kh+0.76*kh**2)/0.88)/(sinh(kh)*cosh(kh)+kh)

    dDdh = dDdGa*dGadHb*dHbdh                 ! eq. (A3), derivative of dissipation with respect to water depth

                                                
    mxdDdh = maxval(dDdh**2)    


   !! bw%sigD =  bw%sigD+(100.d0-bw%sigD)*(1-tanh((3.7*Hrmsmax/h1)**20))  !measurement error for dissipatiopn, used in eq. (6)
    !bw%sigD above is commented out for now to use default value of measurement error (Roberto)
    
        
    bw%sigD = 0.15*maxval(s%dobs)!! dissipation measurement error equal to a percentage (15%) of the maximum observed Dissipation.
    
  !!  do in=1,s%nx+1
  !!    do im=1,s%ny+1  !! here the measurement error is spatially varied, increasing for shallow waters where sandy areas or persistent foam may be. (Roberto)
  !!      if (bw%sigD(in,im)/(h1(in,im)+1.e-16)>bw%sigD(in,im)) bw%sigD(in,im)=bw%sigD(in,im)/(h1(in,im)+1.e-16)
  !!    end do
  !!  end do
                                                    
    alp_as = 0.005                                                  
    errD = ((bw%dcmdo**2+bw%sigD**2)/(dDdh**2+alp_as*mxdDdh+1.e-16)) !eq.(6)  uncertainty in observed data (Dissipation)
                                                                      
    do in=1,s%nx+1
       do im=1,s%ny+1
          if (s%dobs(in,im)<-990) errD(in,im)=999
       end do
    end do

    errC=999 !!! for now, if we want c as a source this needs the change to 
    ! errC(nm) = ((ccmco(nm)**2+sigC(nm)**2)/(dcdh(nm)**2+1.e-16)) without indices
    do in=1,s%nx+1
       do im=1,s%ny+1
          if (bw%cobs(in,im)<-990) errC(in,im)=999
       end do
    end do

    ! for output 
    s%cobs = bw%cobs
    
    errS = ((bw%scmso**2+bw%sigS**2)/(1.d0))                       !eq.(6)  uncertainty in observed data (ibathy)
    do in=1,s%nx+1                             
       do im=1,s%ny+1
          if (s%shobs(in,im)<-990) errS(in,im)=999
       end do
    end do



    Nassim = par%tstop/par%dt    !! changed par%tstop/par%dt -->par%tstop/par%wavint

    sig2obs = 1.0d0/(1./errD+1./errS)
    s%bwalpha = s%sig2prior/(s%sig2prior+Nassim*sig2obs)  ! eq.(2) optimal weighting of prior and observed estimates
    !ap2 s%bwalpha = s%bwalpha*tanh(h1**20)                   ! tanh(h1**20) included to reduce effect at shallow waters...

     
     
    do in=1,s%nx+1
       do im=1,s%ny+1
         alphafac(in,im) = (cosh(s%sdist(in,im)/100.d0-0.65*par%t*s%sdist(s%nx+1,im)/par%tstop/100.d0-2.d0))**(-10.d0)
       enddo
    enddo
    ! max left hand to 1
    do im=1,s%ny+1
       ind = minval(maxloc(alphafac(:,im)))
       alphafac(1:ind,im)=1.d0
    enddo
        
      !    s%bwalpha = s%bwalpha*alphafac
            
    bw%dassim = -s%bwalpha*tanh((h1/0.85)**5)*((dDdh-sqrt(alp_as*mxdDdh))/(dDdh**2+alp_as*mxdDdh)*bw%dcmdo  )*alphafac !eq.(5) depth change
    !ap2                        
    !& +dcdh(nm)/(dcdh(nm)**2+alp_as*dcdhmean**2)*ccmco(nm))

    do in=1,s%nx+1
       do im=1,s%ny+1
          if (bw%dassim(in,im)>0.) then
             bw%dassim(in,im) = min(bw%dassim(in,im),0.1*(h1(in,im)))
          elseif (bw%dassim(in,im)<0.) then
             bw%dassim(in,im) = max(bw%dassim(in,im),-0.1*(h1(in,im))) 
          else
             bw%dassim(in,im) = 0.
          end if
       end do
    end do


    bw%dassim = bw%dassim + s%bwalpha*bw%scmso !!! in xbeach zb=zb-bw%dassim, so zb=zb-s%bwalpha*(zbcomp-zbobs)
    !  bw%dassim = bw%dassim - 0.0001*bw%scmso

    !expose to output
    s%dassim = bw%dassim

    do in=1,s%nx+1
       do im=1,s%ny+1                                               !!Roberto changed tanh term below
          if (s%dobs(in,im)>-990) s%sig2prior(in,im) = s%bwalpha(in,im)*tanh((h1(in,im)/0.85)**5)*Nassim*sig2obs(in,im) 
          !!  if (s%dobs(in,im)>-990) s%sig2prior(in,im) = s%bwalpha(in,im)*tanh((h2(in,im)/0.85)**5)*Nassim*sig2obs(in,im) 
       end do
    end do

  end subroutine comp_depchg

end module beachwizard_module
