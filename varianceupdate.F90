module means_module
  ! This module uses the On-line algorithm according to Knuth (1998) to
  ! determine the variance on the fly
  ! Note: this approximation can be improved.
  use xmpi_module
  use mnemmodule

  type meanspars
     ! Name
     character(maxnamelen)             :: name
     ! Rank
     integer                           :: rank

     ! Array type
     type(arraytype)                   :: t

     ! Keep time average variables
     real*8,dimension(:,:), pointer:: mean2d   => NULL()
     real*8,dimension(:,:,:), pointer:: mean3d   => NULL()
     real*8,dimension(:,:,:,:), pointer:: mean4d=> NULL()
     ! Keep variance of variables   
     real*8,dimension(:,:), pointer:: variance2d     => NULL()
     real*8,dimension(:,:,:), pointer:: variance3d     => NULL()
     real*8,dimension(:,:,:,:), pointer:: variance4d     => NULL()
     ! Needed for variance calculation
     real*8,dimension(:,:), pointer:: variancecrossterm2d  => NULL()
     real*8,dimension(:,:,:), pointer:: variancecrossterm3d  => NULL()
     real*8,dimension(:,:,:,:), pointer:: variancecrossterm4d=> NULL()
     ! Needed for variance calculation  
     real*8,dimension(:,:), pointer:: variancesquareterm2d=> NULL()
     real*8,dimension(:,:,:), pointer:: variancesquareterm3d => NULL()
     real*8,dimension(:,:,:,:), pointer:: variancesquareterm4d => NULL()
     ! Keep time min variables
     real*8,dimension(:,:), pointer:: min2d   => NULL()
     real*8,dimension(:,:,:), pointer:: min3d   => NULL()
     real*8,dimension(:,:,:,:), pointer:: min4d   => NULL()
     ! Keep time max variables
     real*8,dimension(:,:), pointer:: max2d   => NULL()
     real*8,dimension(:,:,:), pointer:: max3d   => NULL()
     real*8,dimension(:,:,:,:), pointer:: max4d   => NULL()
  end type meanspars

  type(meanspars),dimension(:),allocatable  :: meansparsglobal
  type(meanspars),dimension(:),allocatable  :: meansparslocal

contains

  subroutine means_init(sg,sl,par)
    use spaceparams
    use params
    use mnemmodule

    implicit none

    type(spacepars),intent(in)        :: sg,sl
    type(parameters),intent(in)       :: par
    type(arraytype)                   :: t
    integer                           :: i,index,d1,d2,d3,d4

    if (par%nmeanvar>0) then
    
        allocate(meansparsglobal(par%nmeanvar))
        allocate(meansparslocal(par%nmeanvar))
        
        do i=1,par%nmeanvar
           index=chartoindex(trim(par%meanvars(i)))
           call indextos(sg,index,t)
           meansparsglobal(i)%name=t%name
           meansparsglobal(i)%rank=t%rank
           meansparsglobal(i)%t = t
           select case (t%rank)
           case (2)
              if (t%type == 'r') then
                 d1 = size(t%r2,1)
                 d2 = size(t%r2,2)
              elseif (t%type == 'i') then
                 d1 = size(t%i2,1)
                 d2 = size(t%i2,2)
              else
                 call halt_program
              end if
              allocate(meansparsglobal(i)%mean2d(d1,d2))
              allocate(meansparsglobal(i)%variance2d(d1,d2))
              allocate(meansparsglobal(i)%variancecrossterm2d(d1,d2))
              allocate(meansparsglobal(i)%variancesquareterm2d(d1,d2))
              allocate(meansparsglobal(i)%min2d(d1,d2))
              allocate(meansparsglobal(i)%max2d(d1,d2))
              meansparsglobal(i)%mean2d = 0.d0
              meansparsglobal(i)%variance2d = 0.d0
              meansparsglobal(i)%variancecrossterm2d = 0.d0
              meansparsglobal(i)%variancesquareterm2d = 0.d0
              meansparsglobal(i)%min2d = huge(0.d0)
              meansparsglobal(i)%max2d = -1.d0*huge(0.d0)
           case (3)
              if (t%type == 'r') then
                 d1 = size(t%r3,1)
                 d2 = size(t%r3,2)
                 d3 = size(t%r3,3)
              elseif (t%type == 'i') then
                 d1 = size(t%r3,1)
                 d2 = size(t%r3,2)
                 d3 = size(t%r3,3)
              else
                 call halt_program
              end if
              allocate(meansparsglobal(i)%mean3d(d1,d2,d3))
              allocate(meansparsglobal(i)%variance3d(d1,d2,d3))
              allocate(meansparsglobal(i)%variancecrossterm3d(d1,d2,d3))
              allocate(meansparsglobal(i)%variancesquareterm3d(d1,d2,d3))
              allocate(meansparsglobal(i)%min3d(d1,d2,d3))
              allocate(meansparsglobal(i)%max3d(d1,d2,d3))
              meansparsglobal(i)%mean3d = 0.d0
              meansparsglobal(i)%variance3d = 0.d0
              meansparsglobal(i)%variancecrossterm3d = 0.d0
              meansparsglobal(i)%variancesquareterm3d = 0.d0
              meansparsglobal(i)%min3d = huge(0.d0)
              meansparsglobal(i)%max3d = -1.d0*huge(0.d0)
           case (4)
              if (t%type == 'r') then
                 d1 = size(t%r4,1)
                 d2 = size(t%r4,2)
                 d3 = size(t%r4,3)
                 d4 = size(t%r4,4)
              elseif (t%type == 'i') then
                 d1 = size(t%r4,1)
                 d2 = size(t%r4,2)
                 d3 = size(t%r4,3)
                 d4 = size(t%r4,4)
              else
                 call halt_program
              end if
              allocate(meansparsglobal(i)%mean4d(d1,d2,d3,d4))
              allocate(meansparsglobal(i)%variance4d(d1,d2,d3,d4))
              allocate(meansparsglobal(i)%variancecrossterm4d(d1,d2,d3,d4))
              allocate(meansparsglobal(i)%variancesquareterm4d(d1,d2,d3,d4))
              allocate(meansparsglobal(i)%min4d(d1,d2,d3,d4))
              allocate(meansparsglobal(i)%max4d(d1,d2,d3,d4))
              meansparsglobal(i)%mean4d = 0.d0
              meansparsglobal(i)%variance4d = 0.d0
              meansparsglobal(i)%variancecrossterm4d = 0.d0
              meansparsglobal(i)%variancesquareterm4d = 0.d0
              meansparsglobal(i)%min4d = huge(0.d0)
              meansparsglobal(i)%max4d = -1.d0*huge(0.d0)
           end select
#ifdef USEMPI
           call indextos(sl,index,t)
#else
           ! do nothing, we use s global here, no local s available
#endif
           meansparslocal(i)%name=t%name
           meansparslocal(i)%rank=t%rank
           meansparslocal(i)%t = t
           select case (t%rank)
           case (2)
              if (t%type == 'r') then
                 d1 = size(t%r2,1)
                 d2 = size(t%r2,2)
              elseif (t%type == 'i') then
                 d1 = size(t%r2,1)
                 d2 = size(t%r2,2)
              else
                 call halt_program
              end if
              allocate(meansparslocal(i)%mean2d(d1,d2))
              allocate(meansparslocal(i)%variance2d(d1,d2))
              allocate(meansparslocal(i)%variancecrossterm2d(d1,d2))
              allocate(meansparslocal(i)%variancesquareterm2d(d1,d2))
              allocate(meansparslocal(i)%min2d(d1,d2))
              allocate(meansparslocal(i)%max2d(d1,d2))
              meansparslocal(i)%mean2d = 0.d0
              meansparslocal(i)%variance2d = 0.d0
              meansparslocal(i)%variancecrossterm2d = 0.d0
              meansparslocal(i)%variancesquareterm2d = 0.d0
              meansparslocal(i)%min2d = huge(0.d0)
              meansparslocal(i)%max2d = -1.d0*huge(0.d0)
           case (3)
              if (t%type == 'r') then
                 d1 = size(t%r3,1)
                 d2 = size(t%r3,2)
                 d3 = size(t%r3,3)
              elseif (t%type == 'i') then
                 d1 = size(t%r3,1)
                 d2 = size(t%r3,2)
                 d3 = size(t%r3,3)
              else
                 call halt_program
              end if
              allocate(meansparslocal(i)%mean3d(d1,d2,d3))
              allocate(meansparslocal(i)%variance3d(d1,d2,d3))
              allocate(meansparslocal(i)%variancecrossterm3d(d1,d2,d3))
              allocate(meansparslocal(i)%variancesquareterm3d(d1,d2,d3))
              allocate(meansparslocal(i)%min3d(d1,d2,d3))
              allocate(meansparslocal(i)%max3d(d1,d2,d3))
              meansparslocal(i)%mean3d = 0.d0
              meansparslocal(i)%variance3d = 0.d0
              meansparslocal(i)%variancecrossterm3d = 0.d0
              meansparslocal(i)%variancesquareterm3d = 0.d0
              meansparslocal(i)%min3d = huge(0.d0)
              meansparslocal(i)%max3d = -1.d0*huge(0.d0)
           case (4)
              if (t%type == 'r') then
                 d1 = size(t%r4,1)
                 d2 = size(t%r4,2)
                 d3 = size(t%r4,3)
                 d4 = size(t%r4,4)
              elseif (t%type == 'i') then
                 d1 = size(t%r4,1)
                 d2 = size(t%r4,2)
                 d3 = size(t%r4,3)
                 d4 = size(t%r4,4)
              else
                 call halt_program
              end if
              allocate(meansparslocal(i)%mean4d(d1,d2,d3,d4))
              allocate(meansparslocal(i)%variance4d(d1,d2,d3,d4))
              allocate(meansparslocal(i)%variancecrossterm4d(d1,d2,d3,d4))
              allocate(meansparslocal(i)%variancesquareterm4d(d1,d2,d3,d4))
              allocate(meansparslocal(i)%min4d(d1,d2,d3,d4))
              allocate(meansparslocal(i)%max4d(d1,d2,d3,d4))
              meansparslocal(i)%mean4d = 0.d0
              meansparslocal(i)%variance4d = 0.d0
              meansparslocal(i)%variancecrossterm4d = 0.d0
              meansparslocal(i)%variancesquareterm4d = 0.d0
              meansparslocal(i)%min4d = huge(0.d0)
              meansparslocal(i)%max4d = -1.d0*huge(0.d0)
           end select

        enddo
    endif
  end subroutine means_init


  subroutine makeaverage(sl,par)

    use params
    use spaceparams
    use mnemmodule
    use logging_module
    use postprocessmod
    use timestep_module

    IMPLICIT NONE

    type(parameters),   intent(in)                      :: par
    type(spacepars),    intent(in)                      :: sl

    ! keep track of which mean variables are used
    integer                                             :: index 
    integer                                             :: i
    real*8                                              :: mult
    
    type(arraytype)                                     :: t
    
    integer,dimension(sl%nx+1,sl%ny+1)                  :: tvar2di
    
    real*8, dimension(sl%nx+1,sl%ny+1)                  :: oldmean2d,tvar2d
    real*8, dimension(:,:,:),           allocatable     :: oldmean3d,tvar3d
    real*8, dimension(:,:,:,:),         allocatable     :: oldmean4d,tvar4d



    !	avgtime = tpar%tpm(tpar%itm+1)-tpar%tpm(tpar%itm) 
    ! updated in varoutput, not needed to calculate here
    mult = max(par%dt/par%tintm,0.d0) ! Catch initialization at t=0

    do i=1,par%nmeanvar
       index=chartoindex(par%meanvars(i))
       call indextos(sl,index,t)
       select case (t%rank)
       case (2)
          oldmean2d=meansparslocal(i)%mean2d
          ! Robert: One where, elsewhere, endwhere statement leads to memory leak. 
          ! maybe ifort bug ?
          where (oldmean2d<tiny(0.d0) .and. oldmean2d>=0.d0)
             oldmean2d=tiny(0.d0)
          endwhere
          where (oldmean2d>-1.d0*tiny(0.d0) .and. oldmean2d<0.d0)
             oldmean2d=-1.d0*tiny(0.d0)
          endwhere
          ! Some variables (vectors) are rotated to N-S and E-W direction
          if (t%type=='i') then 
             if (par%rotate .eq. 1) call gridrotate(sl,t,tvar2di)
             tvar2d=dble(tvar2di)
          else
             if (par%rotate .eq. 1) call gridrotate(sl,t,tvar2d)
          endif
          if (par%meanvars(i)=='thetamean') then
            meansparslocal(i)%mean2d = meansparslocal(i)%mean2d +   &
                nint((1+mult*sin(tvar2d))*1e6)*1e3              +   &
                nint((1+mult*cos(tvar2d))*1e6)/1e9
          else
            meansparslocal(i)%mean2d = meansparslocal(i)%mean2d + mult*tvar2d
          endif
          meansparslocal(i)%variancecrossterm2d = &
               meansparslocal(i)%variancecrossterm2d/oldmean2d*meansparslocal(i)%mean2d + &
               mult*2.d0*tvar2d*meansparslocal(i)%mean2d
          meansparslocal(i)%variancesquareterm2d = &
               meansparslocal(i)%variancesquareterm2d+mult*(tvar2d)**2
          meansparslocal(i)%variance2d = &
               meansparslocal(i)%variancesquareterm2d-meansparslocal(i)%variancecrossterm2d + &
               meansparslocal(i)%mean2d**2
          meansparslocal(i)%max2d = max(meansparslocal(i)%max2d,tvar2d)
          meansparslocal(i)%min2d = min(meansparslocal(i)%min2d,tvar2d)
       case (3)
          allocate(oldmean3d(size(t%r3,1),size(t%r3,2),size(t%r3,3)))
          allocate(tvar3d(size(t%r3,1),size(t%r3,2),size(t%r3,3)))
          oldmean3d=meansparslocal(i)%mean3d
          ! bug in elsewere --> memory leak see commet Robert above
          where (oldmean3d<tiny(0.d0) .and. oldmean3d>=0.d0)
             oldmean3d=tiny(0.d0)
          endwhere
          where (oldmean3d>-1.d0*tiny(0.d0) .and. oldmean3d<0.d0)
             oldmean3d=-1.d0*tiny(0.d0)
          endwhere
          if (par%rotate .eq. 1) call gridrotate(sl,t,tvar3d)
          meansparslocal(i)%mean3d = meansparslocal(i)%mean3d + mult*tvar3d
          meansparslocal(i)%variancecrossterm3d = &
               meansparslocal(i)%variancecrossterm3d/oldmean3d*meansparslocal(i)%mean3d + &
               mult*2.d0*tvar3d*meansparslocal(i)%mean3d
          meansparslocal(i)%variancesquareterm3d = &
               meansparslocal(i)%variancesquareterm3d+mult*(tvar3d)**2
          meansparslocal(i)%variance3d = &
               meansparslocal(i)%variancesquareterm3d-meansparslocal(i)%variancecrossterm3d + &
               meansparslocal(i)%mean3d**2
          meansparslocal(i)%max3d = max(meansparslocal(i)%max3d,tvar3d)
          meansparslocal(i)%min3d = min(meansparslocal(i)%min3d,tvar3d)
          deallocate(oldmean3d,tvar3d)
       case (4)
          allocate(oldmean4d(size(t%r3,1),size(t%r3,2),size(t%r3,3),size(t%r4,4)))
          allocate(tvar4d(size(t%r3,1),size(t%r3,2),size(t%r3,3),size(t%r4,4)))
          oldmean4d=meansparslocal(i)%mean4d
          where (oldmean4d<tiny(0.d0) .and. oldmean4d>=0.d0)
             oldmean4d=tiny(0.d0)
          endwhere
          where (oldmean4d>-1.d0*tiny(0.d0) .and. oldmean4d<0.d0)
             oldmean4d=-1.d0*tiny(0.d0)
          endwhere
          if (par%rotate .eq. 1) call gridrotate(sl,t,tvar4d)
          meansparslocal(i)%mean4d = meansparslocal(i)%mean4d + mult*tvar4d
          meansparslocal(i)%variancecrossterm4d = &
               meansparslocal(i)%variancecrossterm4d/oldmean4d*meansparslocal(i)%mean4d + &
               mult*2.d0*tvar4d*meansparslocal(i)%mean4d
          meansparslocal(i)%variancesquareterm4d = &
               meansparslocal(i)%variancesquareterm4d+mult*(tvar4d)**2
          meansparslocal(i)%variance4d = &
               meansparslocal(i)%variancesquareterm4d-meansparslocal(i)%variancecrossterm4d + &
               meansparslocal(i)%mean4d**2
          meansparslocal(i)%max4d = max(meansparslocal(i)%max4d,tvar4d)
          meansparslocal(i)%min4d = min(meansparslocal(i)%min4d,tvar4d)
          deallocate(oldmean4d,tvar4d)
       end select
    enddo ! par%nmeanvar

  end subroutine makeaverage
  
  ! clear averages for new averaging period
  subroutine clearaverage(par)
  
    use params
    
    implicit none

    integer                         :: i
    type(parameters), intent(in)    :: par
    
    do i=1,par%nmeanvar
      select case (meansparsglobal(i)%rank)
        case (2)
          meansparslocal(i)%mean2d = 0.d0
          meansparslocal(i)%variance2d = 0.d0
          meansparslocal(i)%variancecrossterm2d = 0.d0
          meansparslocal(i)%variancesquareterm2d = 0.d0
          meansparslocal(i)%min2d = huge(0.d0)
          meansparslocal(i)%max2d = -1.d0*huge(0.d0)
        case (3)
          meansparslocal(i)%mean3d = 0.d0
          meansparslocal(i)%variance3d = 0.d0
          meansparslocal(i)%variancecrossterm3d = 0.d0
          meansparslocal(i)%variancesquareterm3d = 0.d0
          meansparslocal(i)%min3d = huge(0.d0)
          meansparslocal(i)%max3d = -1.d0*huge(0.d0)
        case (4)
          meansparslocal(i)%mean4d = 0.d0
          meansparslocal(i)%variance4d = 0.d0
          meansparslocal(i)%variancecrossterm4d = 0.d0
          meansparslocal(i)%variancesquareterm4d = 0.d0
          meansparslocal(i)%min4d = huge(0.d0)
          meansparslocal(i)%max4d = -1.d0*huge(0.d0)
      end select
    enddo
  end subroutine clearaverage

  ! Subroutine make time series file name

  subroutine makeaveragenames(counter,fnamemean,fnamevar,fnamemin&
       &,fnamemax)
    use mnemmodule

    implicit none

    character(99)      :: fnamemean,fnamevar,fnamemin,fnamemax
    integer            :: counter

    fnamemean = trim(mnemonics(counter))//'_mean.dat'
    fnamevar  = trim(mnemonics(counter))//'_var.dat'
    fnamemin  = trim(mnemonics(counter))//'_min.dat'
    fnamemax  = trim(mnemonics(counter))//'_max.dat'

  end subroutine makeaveragenames

  subroutine makecrossvector(s,sl,crossvararray,nvar,varindexvec,mg&
       &,cstype)
    use params
    use spaceparams
    use mnemmodule

    IMPLICIT NONE

    type(spacepars), intent(in)       :: s,sl    
    real*8,dimension(:,:)             :: crossvararray
    integer,intent(in)                :: mg,cstype,nvar
    integer,dimension(nvar),intent(in):: varindexvec
    type(arraytype)                   :: t

    integer                           :: i

    crossvararray=-999.d0
    do i=1,nvar
#ifdef USEMPI
       call space_collect_index(s,sl,varindexvec(i))
#endif  
       if (xmaster) then
          call indextos(s,varindexvec(i),t)
          select case(t%type)
          case('r')
             if(cstype==0) then
                crossvararray(i,:)=t%r2(:,mg)
             else
                crossvararray(i,:)=t%r2(mg,:)
             endif
          case('i')
             if(cstype==0) then
                crossvararray(i,:)=t%i2(:,mg)
             else
                crossvararray(i,:)=t%i2(mg,:)
             endif
          end select
       endif
    enddo
#ifdef USEMPI
    call xmpi_bcast(crossvararray)
#endif 
  end subroutine makecrossvector

#ifdef USEMPI
  subroutine means_collect(sl,a,b)
    use spaceparams

    implicit none

    type(spacepars),intent(in)  :: sl
    type(meanspars),intent(inout) :: a
    type(meanspars),intent(in)    :: b

    select case (b%rank)
    case (2)
       call space_collect(sl,a%mean2d,b%mean2d)
       call space_collect(sl,a%variance2d,b%variance2d)
       call space_collect(sl,a%max2d,b%max2d)
       call space_collect(sl,a%min2d,b%min2d)
       ! for H and urms :
       call space_collect(sl,a%variancesquareterm2d,b%variancesquareterm2d)
    case (3)
       call space_collect(sl,a%mean3d,b%mean3d)
       call space_collect(sl,a%variance3d,b%variance3d)
       call space_collect(sl,a%max3d,b%max3d)
       call space_collect(sl,a%min3d,b%min3d)
    case (4)
       call space_collect(sl,a%mean4d,b%mean4d)
       call space_collect(sl,a%variance4d,b%variance4d)
       call space_collect(sl,a%max4d,b%max4d)
       call space_collect(sl,a%min4d,b%min4d)
    end select

  end subroutine means_collect
#endif


end module means_module
