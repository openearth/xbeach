!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   MODULE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module outputmod
use xmpi_module
use mnemmodule
use varianceupdate_module

implicit none
private
public init_output, var_output

! Robert: Add choice of output variables
!logical,dimension(999)              :: outputindex ! [-]      tracks which  global variables are to be outputted.
integer*4                           :: nglobalvar  ! number of global output variables
integer*4                           :: npoints     ! number of output points
integer*4                           :: nrugauge    ! number of runup gauges
integer*4                           :: ncross      ! number of cross section profiles
integer                             :: timings     ! 0: no timings, 1 timings output
integer*4,dimension(:),allocatable  :: pointtype   ! 0 = point output, 1 = runup gauge
integer*4,dimension(:),allocatable  :: crosstype   ! 0 = cross shore (x), 1 = longshore (y)
integer*4,dimension(:),allocatable  :: xpoints     ! model x-coordinate of output points
integer*4,dimension(:),allocatable  :: ypoints     ! model y-coordinate of output points
integer*4,dimension(:),allocatable  :: xcross      ! model x-coordinate of output cross sections
integer*4,dimension(:),allocatable  :: ycross      ! model y-coordinate of output cross sections
integer*4,dimension(:),allocatable  :: nvarpoint   ! vector with number of output variable per output point
integer*4,dimension(:),allocatable  :: nvarcross   ! vector with number of output variable per output cross section
integer*4,dimension(:,:),allocatable:: Avarpoint   ! Array with associated index of output variables per point
integer*4,dimension(:,:),allocatable:: Avarcross   ! Array with associated index of output variables per cross section
                                                   ! Only alive at xmaster
integer*4                           :: nmeanvar    ! number of time-average variables
integer*4,dimension(:),allocatable  :: meanvec     ! keep track of which mean variables are used
real*8,dimension(:),allocatable     :: tpg,tpp,tpm,tpc,tpw ! output time points
integer*4                           :: stpm        ! size of tpm

integer                             :: noutnumbers = 0  ! the number of outnumbers

integer, dimension(numvars)         :: outnumbers  ! numbers, corrsponding to mnemonics, which are to be output
integer                             :: ndt,itg,itp,itc,itm,day,ot,itprev
real*8,dimension(10)                :: tlast10
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



subroutine init_output(s,sl,par,it)
  use params
  use spaceparams
  use readkey_module

  IMPLICIT NONE

  type(spacepars),intent(inout)       :: s,sl
  type(parameters)                    :: par

  integer                             :: id,ic,icold,i,ii,index,it
  integer,dimension(2)                :: minlocation
  character(80)                       :: line, keyword
  character(80)                       :: var, fname
  integer, dimension(:,:),allocatable :: temparray
  real*8,dimension(s%nx+1,s%ny+1)	  :: mindist
  real*8                              :: tg1,tp1,tm1,tw1,tc1
  real*8,dimension(:),allocatable     :: xpointsw,ypointsw
  real*8,dimension(:),allocatable     :: xcrossw,ycrossw
  character(1)								  :: char
  

  !!!!! OUTPUT POINTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  npoints  = readkey_int ('params.txt','npoints',      0,       0,     50)
  nrugauge = readkey_int ('params.txt','nrugauge',     0,       0,     50)

  timings  = readkey_int ('params.txt','timings',      1,       0,      1)

  if ((npoints+nrugauge)>0) then 
      allocate(pointtype(npoints+nrugauge))
      allocate(xpoints(npoints+nrugauge))
      allocate(ypoints(npoints+nrugauge))
	  allocate(xpointsw(npoints+nrugauge))
      allocate(ypointsw(npoints+nrugauge))
      allocate(nvarpoint(npoints+nrugauge))
      allocate(temparray(npoints+nrugauge,99))
      pointtype(1:npoints)=0
      pointtype(npoints+1:npoints+nrugauge)=1
      temparray=-1

      if (xmaster) then
        if (npoints>0) then
          id=0
          ! Look for keyword npoints in params.txt
          open(10,file='params.txt')
          do while (id == 0)
            read(10,'(a)')line
            ic=scan(line,'=')
            if (ic>0) then
               keyword=adjustl(line(1:ic-1))
               if (keyword == 'npoints') id=1
            endif
          enddo

          do i=1,npoints
              read(10,*) xpointsw(i),ypointsw(i),nvarpoint(i),line
			  write(*,'(a,i0)') ' Output point ',i
              write(*,'(a,f0.2,a,f0.2)') ' xpoint: ',xpointsw(i),'   ypoint: ',ypointsw(i)
              ! Convert world coordinates of points to nearest (lsm) grid point
			     mindist=sqrt((xpointsw(i)-s%xw)**2+(ypointsw(i)-s%yw)**2)
			     minlocation=minloc(mindist)
			     xpoints(i)=minlocation(1)
			     ypoints(i)=minlocation(2)
			     write(*,'(a,i0,a,i0,a,f0.2,a)')' Distance output point to nearest grid point ('&
			                        ,minlocation(1),',',minlocation(2),') is '&
			                        ,mindist(minlocation(1),minlocation(2)), ' meters'

              icold=0
              do ii =1,nvarpoint(i)
                  ic=scan(line(icold+1:80),'#')
                  ic=ic+icold
                  var=line(icold+1:ic-1)
                  index = chartoindex(var)

                  if (index/=-1) then
                      temparray(i,ii)=index
                      write(*,*)' Output type: "'//trim(var)//'"'
                  else
                      write(*,*)' Unknown point output type: "'//trim(var)//'"'
                      call halt_program
                  endif
                  icold=ic
              enddo
			 
          enddo
          close(10)
        endif ! npoints>0  

        if (nrugauge>0) then
          id=0
          ! Look for keyword nrugauge in params.txt
          open(10,file='params.txt')
          do while (id == 0)
            read(10,'(a)')line
            ic=scan(line,'=')
            if (ic>0) then
               keyword=adjustl(line(1:ic-1))
               if (keyword == 'nrugauge') id=1
            endif
          enddo

          do i=1+npoints,nrugauge+npoints
              read(10,*) xpointsw(i),ypointsw(i),nvarpoint(i),line
			  write(*,'(a,i0)') ' Output runup gauge ',i-npoints
              write(*,'(a,f0.2,a,f0.2)') ' xpoint: ',xpointsw(i),'   ypoint: ',ypointsw(i)
				  ! Convert world coordinates of points to nearest (lsm) grid row
			     mindist=sqrt((xpointsw(i)-s%xw)**2+(ypointsw(i)-s%yw)**2)
			     minlocation=minloc(mindist)
			     xpoints(i)=1
			     ypoints(i)=minlocation(2)
			     write(*,'(a,i0)')'Runup gauge at grid line iy=',ypoints(i)
              icold=0
              do ii =1,nvarpoint(i)
                  ic=scan(line(icold+1:80),'#')
                  ic=ic+icold
                  var=line(icold+1:ic-1)
                  index = chartoindex(var)

                  if (index/=-1) then
                      temparray(i,ii)=index
                      write(*,*)' Output type:"'//trim(var)//'"'
                  else
                      write(*,*)' Unknown point output type: "'//trim(var)//'"'
                      call halt_program
                  endif
                  icold=ic
              enddo
			  
          enddo
          close(10)
        endif  ! nrugauge > 0	    
      endif ! xmaster

	 

#ifdef USEMPI
      call xmpi_bcast(nvarpoint)
      call xmpi_bcast(pointtype)
      call xmpi_bcast(xpoints)
      call xmpi_bcast(ypoints)
      call xmpi_bcast(nvarpoint)
      call xmpi_bcast(temparray)   
#endif
      ! Tidy up information
      allocate(Avarpoint(npoints+nrugauge,maxval(nvarpoint)))
      Avarpoint(:,:)=temparray(:,1:maxval(nvarpoint))
      deallocate (temparray)
  endif ! npoints + nrugauge > 0


  !!!!! CROSS SECTION OUTPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ncross  = readkey_int ('params.txt','ncross',      0,       0,     50)

  if (ncross>0) then
     allocate(crosstype(ncross))
     allocate(xcross(ncross))
     allocate(ycross(ncross))
	  allocate(xcrossw(ncross))
     allocate(ycrossw(ncross))
     allocate(nvarcross(ncross))
     allocate(temparray(ncross,99))
     temparray=-1

	  if (xmaster) then
	     id=0
        ! Look for keyword ncross in params.txt
        open(10,file='params.txt')
        do while (id == 0)
          read(10,'(a)')line
          ic=scan(line,'=')
          if (ic>0) then
            keyword=adjustl(line(1:ic-1))
            if (keyword == 'ncross') id=1
          endif
        enddo
        
		  do i=1,ncross
!          read(10,*) xcrossw(i),ycrossw(i),crosstype(i),nvarcross(i),line
          read(10,*) xcrossw(i),ycrossw(i),char,nvarcross(i),line
			 write(*,'(a,i0)') ' Cross section ',i
			 if(char=='x'.or.char=='X') then
			    crosstype(i)=0
			 elseif(char=='y'.or.char=='Y') then
			    crosstype(i)=1
			 else
                write(*,*)' Unknown cross section type: ',char
                call halt_program
			 endif
			 ! Convert world coordinates of points to nearest (lsm) grid row
			 mindist=(xcrossw(i)-s%xw)**2+(ycrossw(i)-s%yw)**2
			 minlocation=minloc(mindist)
			 xcross(i)=minlocation(1)
			 ycross(i)=minlocation(2)
          write(*,'(a,f0.2,a,f0.2)') ' positioned at x-coordinate: ',xcrossw(i),&
			                                  ' y-coordinate: ',ycrossw(i)			 
			 if (crosstype(i)==0) then
			   write(*,'(a,i0)')' placed in x-direction at row iy = ',ycross(i)
			 else
			   write(*,'(a,i0)')' placed in y-direction at column ix = ',xcross(i) 
			 endif
          icold=0
          do ii =1,nvarcross(i)
            ic=scan(line(icold+1:80),'#')
            ic=ic+icold
            var=line(icold+1:ic-1)
            index = chartoindex(var)
            if (index/=-1) then
              temparray(i,ii)=index
				  call indextos(s,index,At)
		        if (At%rank.ne.2) then
		          write(*,'(a,i0,a,a,a)')' Cross section output not designed for rank ',At%rank,' array "',trim(At%name),'"'
			       call halt_program
		        endif
              write(*,*)' Output type: "'//trim(var)//'"'
            else
              write(*,*)' Unknown point output type: "'//trim(var)//'"'
              call halt_program
            endif
            icold=ic
          enddo
        enddo
        close(10)
     endif !xmaster

#ifdef USEMPI
     call xmpi_bcast(nvarcross)
     call xmpi_bcast(crosstype)
     call xmpi_bcast(xcross)
     call xmpi_bcast(ycross)
     call xmpi_bcast(nvarcross)
     call xmpi_bcast(temparray)   
#endif     

     ! Tidy up information
     allocate(Avarcross(ncross,maxval(nvarcross)))
     Avarcross(:,:)=temparray(:,1:maxval(nvarcross))
     deallocate (temparray)
  endif ! ncross > 0



  !!!!! GLOBAL OUTPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! outputindex is set in this subroutine. Row number corresponds to line in s.inp.
  ! Default choice = H,zs,zb,u,v,cc,Su,Sv, dims
  ! dims is given number 999 !!!!!

  nglobalvar  = readkey_int ('params.txt','nglobalvar',   -1,       -1,     20)

  !outputindex = .false.

  !if (nglobalvar == 0) then                       ! User wants no output
  !        outputindex = .false.
  if (nglobalvar == -1)       then    ! Default output
          call add_outmnem('H')
          call add_outmnem('zs')
          call add_outmnem('zs0')
          call add_outmnem('zb')
          call add_outmnem('hh')
          call add_outmnem('u')
          call add_outmnem('v')
          call add_outmnem('ue')
          call add_outmnem('ve')
          call add_outmnem('urms')
          call add_outmnem('Fx')
          call add_outmnem('Fy')
          call add_outmnem('ccg')    ! todo wwvv: cc does not exist
          call add_outmnem('ceqsg')
		  call add_outmnem('ceqbg')
          call add_outmnem('Sug')
          call add_outmnem('Svg')
          call add_outmnem('E')
          call add_outmnem('R')
          call add_outmnem('D')
          call add_outmnem('DR')
  elseif (nglobalvar == 999) then ! Output all
          do i=1,numvars
            call add_outnumber(i)
          enddo
  else                                                            ! User specified output
          ! Look for keyword nglobalvar in params.txt
          id=0
          if (xmaster) then
            open(10,file='params.txt')
            do while (id == 0)
              read(10,'(a)')line
              ic=scan(line,'=')
              if (ic>0) then
                keyword=adjustl(line(1:ic-1))
                if (keyword == 'nglobalvar') id=1
              endif
            enddo
          ! Read through the variables lines
            do i=1,nglobalvar
                    read(10,'(a)')line
                    line=trim(line)        ! useless   wwvv
                    call add_outmnem(line)
            end do
            close(10)
          endif  ! xmaster
#ifdef USEMPI
          call xmpi_bcast(noutnumbers)
          call xmpi_bcast(outnumbers)
#endif

  endif


  !!!!! TIME-AVEARGE, VARIANCE and MIN-MAX ARRAYS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  nmeanvar = readkey_int ('params.txt','nmeanvar',0,0,15)
  if (nmeanvar>0) then
      allocate(meanvec(nmeanvar))
      id=0
	  ! Look for keyword nmeanvar in params.txt
	  if(xmaster) then
        open(10,file='params.txt')
        do while (id == 0)
          read(10,'(a)')line
          ic=scan(line,'=')
          if (ic>0) then
             keyword=adjustl(line(1:ic-1))
             if (keyword == 'nmeanvar') id=1
          endif
        enddo
		! Read through the variables lines
        do i=1,nmeanvar
            read(10,'(a)')line
            index = chartoindex(trim(line))
            if (index/=-1) then
				call indextos(s,index,At)
		        if (At%rank.ne.2) then
		           write(*,'(a,f2.0,a,a,a)')' Time-average output not designed for rank ',real(At%rank),&
						                                       ' array "',trim(At%name),'"'
			        call halt_program
		        endif
                meanvec(i)=index
                write(*,'(a)')' Will generate mean, min, max and variance output for variable: '//trim(mnemonics(index))
            else
                write(*,*)' Unknown output type ',trim(line)
                call halt_program
            endif
        end do
        close(10)
		call initialize_mean_arrays(s%nx,s%ny,nmeanvar)
      endif  ! xmaster

#ifdef USEMPI
      call xmpi_bcast(meanvec)
#endif

  endif ! nmeanvar > 0

  !!!!! OUTPUT TIME POINTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! If working in instat 0 or instat 40 we want time to coincide at wavint timesteps
  if (par%instat==0 .or. par%instat==40) then
	  if(xmaster) then
	     ii=floor(par%tstop/par%wavint)
	     allocate(tpw(ii))
	     do i=1,ii
	        tpw(i)=real(i)*par%wavint
	     enddo
      endif
#ifdef USEMPI
      call xmpi_bcast(ii)
      if (.not.xmaster) then
        allocate(tpw(ii))
      endif
      call xmpi_bcast(tpw)
#endif
  endif    

  ! If we want global output then
  if (nglobalvar/=0) then
      if(xmaster) then
        call readkey('params.txt','tsglobal',fname)
        if (fname/=' ') then
            open(10,file=fname)
            read(10,*)ii
            allocate(tpg(ii))
            do i=1,ii
                    read(10,*)tpg(i)
            enddo
				tpg=tpg/max(par%morfac,1.d0)
            close(10)
        else
            ii=floor((par%tstop-par%tstart)/par%tintg)+1
            allocate(tpg(ii))
            do i=1,ii
                    tpg(i)=par%tstart+par%tintg*real(i-1)
            enddo
        endif
      endif  !xmaster
#ifdef USEMPI
      call xmpi_bcast(ii)
      if (.not.xmaster) then
        allocate(tpg(ii))
      endif
      call xmpi_bcast(tpg)
#endif
  endif ! nglobalvar /=0

  ! If we want point output then
  if ((npoints+nrugauge)>0) then 
      if (xmaster) then
        call readkey('params.txt','tspoints',fname)
        if (fname/=' ') then
            open(10,file=fname)
            read(10,*)ii
            allocate(tpp(ii))
            do i=1,ii
                    read(10,*)tpp(i)
            enddo
				tpp=tpp/max(par%morfac,1.d0)
            close(10)
        else
            ii=floor((par%tstop-par%tstart)/par%tintp)+1
            allocate(tpp(ii))
            do i=1,ii
                    tpp(i)=par%tstart+par%tintp*real(i-1)
            enddo
        endif
      endif    ! xmaster
#ifdef USEMPI
      call xmpi_bcast(ii)
      if (.not.xmaster) then
        allocate(tpp(ii))
      endif
      call xmpi_bcast(tpp)
#endif
  endif ! (npoints+nrugauge)>0 

! If we want cross section output then
  if ((ncross)>0) then 
      if (xmaster) then
        call readkey('params.txt','tscross',fname)
        if (fname/=' ') then
            open(10,file=fname)
            read(10,*)ii
            allocate(tpc(ii))
            do i=1,ii
               read(10,*)tpc(i)
            enddo
            close(10)
        else
            ii=floor((par%tstop-par%tstart)/par%tintc)+1
            allocate(tpc(ii))
            do i=1,ii
               tpc(i)=par%tstart+par%tintc*real(i-1)
            enddo
        endif
      endif    ! xmaster
#ifdef USEMPI
      call xmpi_bcast(ii)
      if (.not.xmaster) then
        allocate(tpc(ii))
      endif
      call xmpi_bcast(tpc)
#endif
  endif ! (ncross)>0 

  ! If we want time-average output then
  if (nmeanvar>0) then
      if(xmaster) then
        call readkey('params.txt','tsmean',fname)
        if (fname/=' ') then
            open(10,file=fname)
            read(10,*)ii
            allocate(tpm(ii))
            do i=1,ii
                    read(10,*)tpm(i)
            enddo
				tpm=tpm/max(par%morfac,1.d0)
            close(10)
        else
            ii=floor((par%tstop-par%tstart)/par%tintm)+1
			if (ii<=1) then
				write(*,*)'Tintm is larger than output simulation time '
				call halt_program
			endif
            allocate(tpm(ii))
            do i=1,ii
                    tpm(i)=par%tstart+par%tintm*real(i-1)
            enddo
        endif
      endif
#ifdef USEMPI
      call xmpi_bcast(ii)
      if (.not.xmaster) then
        allocate(tpm(ii))
      endif
      call xmpi_bcast(tpm)
#endif
  endif  ! nmeanvar > 0

  ! If tp series not defined, then no output wanted, so large timestep
  if (.not. allocated(tpw)) then 
      allocate(tpw(1))
	  tpw=par%tstop
  endif

  if (.not. allocated(tpg)) then
      allocate(tpg(1))
      tpg=par%tstop
  endif
  if (.not. allocated(tpp)) then
      allocate(tpp(1))
      tpp=par%tstop
  endif
  if (.not. allocated(tpc)) then
      allocate(tpc(1))
      tpc=par%tstop
  endif
  if (.not. allocated(tpm)) then  ! Need minimum two in this array
      allocate(tpm(2))
      tpm(1)=par%tstop
      tpm(2)=par%tstop+1.d0
  endif

  tg1=minval(tpg)
  tp1=minval(tpp)
  tc1=minval(tpc)
  tm1=minval(tpm)
  tw1=minval(tpw)
  par%tintm=tpm(2)-tpm(1)
  stpm=size(tpm)

  if (min(tg1,tp1,tm1,tc1)>0.d0) then    ! No output wanted at initialization
      par%tnext=min(tg1,tp1,tm1,tc1,tw1)
  else
      it=1
      tg1=minval(tpg,MASK=tpg .gt. 0.d0)
      tp1=minval(tpp,MASK=tpp .gt. 0.d0)
		tc1=minval(tpc,MASK=tpc .gt. 0.d0)
      tm1=minval(tpm,MASK=tpm .gt. 0.d0)
		tw1=minval(tpw,MASK=tpw .gt. 0.d0)
      par%tnext=min(tg1,tp1,tm1,tc1,tw1)
  endif


end subroutine init_output

subroutine add_outnumber(number)
  use xmpi_module
  implicit none
  integer, intent(in)  :: number ! to add
  if (noutnumbers .ge. numvars) then
    write(*,'(a,i0)')'Too many outnumbers asked, max is ',numvars
    write(*,*)'Program will stop'
    call halt_program
  endif

  noutnumbers = noutnumbers+1
  outnumbers(noutnumbers) = number

end subroutine add_outnumber

subroutine add_outmnem(mnem)
  implicit none
  character(len=*), intent(in) :: mnem
  integer              :: i
  i = chartoindex(mnem)
  if (i .lt. 1 .or. i .gt. numvars) then
    if(xmaster) then
      write(*,*)'Warning: cannot locate variable "',trim(mnem),'", no output for this one'
    endif
    return
  endif
  call add_outnumber(i)
  if(xmaster) then
    write(*,*)'Will generate global output for variable "',trim(mnem),'"'
  endif
  return
end subroutine add_outmnem



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!   OUTPUT AT EVERY TIMESTEP    !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine var_output(it,s,sl,par)
  use params
  use spaceparams

  IMPLICIT NONE

  type(spacepars)                         :: s,sl
  type(parameters)                        :: par

  integer                                 :: i,ii,tt,k
  integer                                 :: it,i1,i2,i3
  integer                                 :: wordsize, idum
  integer                                 :: reclen,reclen2,reclenc,reclenp
  real*8                                  :: dtimestep,tpredicted,tnow, t1,t2,t3,t4,t5
  character(12)                           :: fname
  character(99)                           :: fnamemean,fnamevar,fnamemin,fnamemax
  real*8,dimension(numvars)               :: intpvector
  real*8,dimension(numvars,s%nx+1)        :: crossvararray0
  real*8,dimension(numvars,s%ny+1)        :: crossvararray1
  integer,dimension(:),allocatable        :: tempvectori
  real*8,dimension(:),allocatable         :: tempvectorr, temp
  integer,dimension(8)                    :: datetime
  real*8,dimension(size(tpg)+size(tpp)+size(tpc)+size(tpm)) :: outputtimes
  type(arraytype)                         :: t
  real*8,dimension(:,:),pointer				:: MI,MA
  
  inquire(iolength=wordsize) 1.d0
  reclen=wordsize*(s%nx+1)*(s%ny+1)
  reclen2=wordsize*(s%nx+1)*(s%ny+1)*(par%ngd)*(par%nd)
  
  ! Complicated check, so only carry out if it = 0 or it = 1
  if (it<2) then 
  ! Do only at very first time step
	  if ((abs(par%t-par%dt)<1.d-6.and.it==0).or.(par%t<1.d-6.and.it==1)) then
		  ndt=0
		  tlast10=0.d0
		  day=0
		  ot=0
		  itprev=0
	  endif
  endif

  ndt=ndt+1                                       ! Number of calculation time steps per output time step

  if (nmeanvar/=0) then
    if (par%t>tpm(1) .and. par%t<=tpm(stpm)) then
	     call makeaverage(s,sl,par,meanvec, nmeanvar)                ! Make averages and min-max every flow timestep
	 endif
  endif


  if (par%t>=(dble(ot)*par%tstop/500.d0)) then
      ot=ot+1
      ! Predicted time
      !call cpu_time(tnow)
      call date_and_time(VALUES=datetime)
      tnow=day*24.d0*3600.d0+datetime(5)*3600.d0+60.d0*datetime(6)+1.d0*datetime(7)+0.001d0*datetime(8)

      if (tnow<minval(tlast10)) then
          day=day+1
          tnow=day*24.d0*3600.d0+datetime(5)*3600.d0+60.d0*datetime(6)+1.d0*datetime(7)+0.001d0*datetime(8)
      endif

      tlast10(minloc(tlast10))=tnow
      if (ot==1) then
          dtimestep=0.d0
      elseif (ot<=9) then
          tt=count(tlast10 .gt. 0.d0)
          dtimestep=(maxval(tlast10)-minval(tlast10,MASK=tlast10 .gt. (0.d0)))/real(max(1,tt-1))
      else
          dtimestep=(maxval(tlast10)-minval(tlast10))/9.d0
      endif

      !       tpredicted=((par%tstop/par%tint)-it)*dtimestep
      if(timings .ne. 0) then
        if(xmaster) then
          write(*,fmt='(a,f5.1,a)')'Simulation ',100.d0*par%t/par%tstop,&
             ' percent complete'
        endif
      endif

      tpredicted=(par%tstop-par%t)/(par%tstop/500.d0)*dtimestep
      if (timings .ne. 0) then
        if (dtimestep<1 .and. tpredicted<120 .or. .not. xmaster) then 
                ! Write nothing
        elseif (tpredicted>=3600) then 
            write(*,fmt='(a,I3,a,I2,a)')'Time remaining ',&
              floor(tpredicted/3600.0d0),' hours and ',&
              nint((tpredicted-3600.0d0*floor(tpredicted/3600.0d0))/60.0d0),&
              ' minutes'
        elseif (tpredicted>=600) then
            write(*,fmt='(a,I2,a)')'Time remaining ',&
              floor(tpredicted/60.0d0),' minutes'
        elseif (tpredicted>=60) then
            write(*,fmt='(a,I2,a,I2,a)')'Time remaining ',&
              floor(tpredicted/60.0d0),' minutes and ',&
              nint((tpredicted-60.0d0*floor(tpredicted/60.0d0))),' seconds'
        else
          write(*,fmt='(a,I2,a)')'Time remaining ',nint(tpredicted),' seconds'
        endif
      endif
  endif


  ! Determine if this is an output timestep
  if (it>itprev) then
          itprev=it

     !!!!!!! Start writing output

     if (it==1) then   !!! Only first time round: initialisation of files

          itp=0
          if (xmaster) then
            !! First time file opening for point output
            if (npoints+nrugauge>0) then
                do i=1,npoints+nrugauge
                    if (pointtype(i)==0) then
                        fname(1:5)='point'
                        i1=floor(real(i)/100.d0)
                        i2=floor(real(i-i1*100)/10.d0)
                        i3=i-i1*100-i2*10
                    else
                        fname(1:5)='rugau'
                        i1=floor(real(i-npoints)/100.d0)
                        i2=floor(real((i-npoints)-i1*100)/10.d0)
                        i3=(i-npoints)-i1*100-i2*10
                    endif
                    fname(6:6)=char(48+i1)
                    fname(7:7)=char(48+i2)
                    fname(8:8)=char(48+i3)
                    fname(9:12)='.dat'
                    reclenp=wordsize*(nvarpoint(i)+1)*1
                    open(indextopointsunit(i),file=fname,&
                        form='unformatted',access='direct',recl=reclenp)
                enddo
            endif
          endif  ! xmaster


          itc=0
          if (xmaster) then
            !! First time file opening for cross section output
            if (ncross>0) then
                do i=1,ncross
                    fname(1:5)='cross'
                    i1=floor(real(i)/100.d0)
                    i2=floor(real(i-i1*100)/10.d0)
                    i3=i-i1*100-i2*10
                    fname(6:6)=char(48+i1)
                    fname(7:7)=char(48+i2)
                    fname(8:8)=char(48+i3)
                    fname(9:12)='.dat'
						  if (crosstype(i)==0) then
						    reclenc=wordsize*(s%nx+1)*(nvarcross(i))
                    else
						    reclenc=wordsize*(s%ny+1)*(nvarcross(i))
						  endif 
                    open(indextocrossunit(i),file=fname,&
                        form='unformatted',access='direct',recl=reclenc)
                enddo
            endif
          endif  ! xmaster

          !! First time file opening for time-average output
          itm=0
          if(xmaster) then
            if (nmeanvar>0) then
                do i=1,nmeanvar
                    call makeaveragenames(meanvec(i),fnamemean,fnamevar,fnamemin,fnamemax)
                    fnamemean=trim(fnamemean)
						  fnamevar =trim(fnamevar)
						  fnamemin =trim(fnamemin)
						  fnamemax =trim(fnamemax)
                    open(indextomeanunit(i),file=fnamemean, &
                     form='unformatted',access='direct',recl=reclen)
						  open(indextovarunit(i),file=fnamevar, &
                     form='unformatted',access='direct',recl=reclen)
						  open(indextominunit(i),file=fnamemin, &
                     form='unformatted',access='direct',recl=reclen)
						  open(indextomaxunit(i),file=fnamemax, &
                     form='unformatted',access='direct',recl=reclen)
                enddo
            endif
          endif ! xmaster


          !! First time file opening for global output
          itg=0
          if (xmaster) then
            if (nglobalvar/=0) then
                open(100,file='xy.dat',form='unformatted',access='direct',recl=reclen)
                write(100,rec=1)s%xw
                write(100,rec=2)s%yw
                write(100,rec=3)s%x
                write(100,rec=4)s%y
                close(100)
            endif
          endif

!                do i=1,noutnumbers
!                  j = outnumbers(i)
!                  open (indextoglobalunit(j),file=mnemonics(i)//'.dat',form='unformatted',access='direct',recl=reclen) 
!                enddo
     endif                !! End first time opening


    !!! Write at every output timestep


    !!! Write point variables
    !!! Only write if it is output time for points
    if (npoints+nrugauge>0) then
        if (any(abs(par%t-tpp) .le. 0.0000001d0)) then
            itp=itp+1
            ! wwvv check if there is any pointtype.eq. 1
            ! In that case, we need the values in wetz
            ! Probably, this test can be coded somewhat smarter
#ifdef USEMPI
            do i=1,npoints+nrugauge
              if (pointtype(i) .eq. 1) then
                call space_collect(sl,s%wetz,sl%wetz)
                exit
              endif
            enddo
#endif
            do i=1,npoints+nrugauge
                !!! Make vector of all s% values at n,m grid coordinate
                if (pointtype(i)==1) then
                    if (xmaster) then
                      allocate (temp(1:s%nx+1))
                      temp=(/(k,k=1,s%nx+1)/)
                      ! wwvv wetz ok?, probably not
                      ! idum = maxval(maxloc(temp*s%wetz(1:s%nx+1,ypoints(i)))) ! Picks the last "1" in temp
                      idum =  max(maxval(minloc(s%wetz(:,ypoints(i))))-1,1)
                      deallocate(temp)
                    endif
#ifdef USEMPI
                    call xmpi_bcast(idum)
#endif
                    call makeintpvector(par,sl,intpvector,idum,ypoints(i))
                    ! need sl here, because of availability
                    ! of is,js,lm,ln
                else
                    call makeintpvector(par,sl,intpvector,xpoints(i),ypoints(i))
                endif
                if (xmaster) then
                  allocate(tempvectori(nvarpoint(i)))
                  allocate(tempvectorr(nvarpoint(i)+1))
                  tempvectori=Avarpoint(i,1:nvarpoint(i))
                  tempvectorr(1)=par%t*max(par%morfac,1.d0)
                  do ii=1,nvarpoint(i)
                          tempvectorr(ii+1)=intpvector(tempvectori(ii))
                  enddo
                  write(indextopointsunit(i),rec=itp)tempvectorr
                  deallocate(tempvectori)
                  deallocate(tempvectorr)
                endif
            enddo
        endif
    endif

    !!! Write cross section variables
    !!! Only write if it is output time for cross sections
    if (ncross>0) then
        if (any(abs(par%t-tpc) .le. 0.0000001d0)) then
            itc=itc+1
            do i=1,ncross
                !!! Make array of all s% values at n or m grid line
					 if (crosstype(i)==0) then
					   !    makecrossvector(par, s, local s, output array, no of variables, index of variables in output, m or n coordinate, cross section type) 
                  call makecrossvector(s,sl,crossvararray0,nvarcross(i),Avarcross(i,1:nvarcross(i)),ycross(i),0)
						if (xmaster) write(indextocrossunit(i),rec=itc)crossvararray0(1:nvarcross(i),:)
					 else
					   call makecrossvector(s,sl,crossvararray1,nvarcross(i),Avarcross(i,1:nvarcross(i)),xcross(i),1)
						if (xmaster) write(indextocrossunit(i),rec=itc)crossvararray1(1:nvarcross(i),:)
					 endif
            enddo
        endif
    endif

    !!! Write average variables
    if(xmaster) then
      if (nmeanvar>0) then
          ! Not at the first in tpm as this is the start of averaging. Only output after second in tpm
          if (any(abs(par%t-tpm(2:stpm)) .le. 0.0000001d0)) then
              itm=itm+1
              do i=1,nmeanvar
                  if (mnemonics(meanvec(i))=='H') then                ! Hrms changed to H
                      meanarrays(:,:,i)=sqrt(meanarrays(:,:,i))
                      write(indextomeanunit(i),rec=itm)meanarrays(:,:,i)
                  elseif (mnemonics(meanvec(i))=='urms') then    ! urms
                      meanarrays(:,:,i)=sqrt(meanarrays(:,:,i))
                      write(indextomeanunit(i),rec=itm)meanarrays(:,:,i)
                  else                                                    ! non-rms variables
                      write(indextomeanunit(i),rec=itm)meanarrays(:,:,i)
                  endif
						write(indextovarunit(i),rec=itm)variancearrays(:,:,i)
						MI=>minarrays(:,:,i)
						MA=>maxarrays(:,:,i)
						where(MI>0.99d0*huge(0.d0))
						    MI=-999.d0
						endwhere
						where(MA<-0.99d0*huge(0.d0))
						    MA=-999.d0
						endwhere
						write(indextominunit(i),rec=itm)MI
                  write(indextomaxunit(i),rec=itm)MA
              enddo
              meanarrays=0.0d0
				  variancearrays=0.d0
				  variancecrossterm=0.d0
				  variancesquareterm=0.d0
				  minarrays=huge(0.d0)
				  maxarrays=-1.d0*huge(0.d0)
              par%tintm=tpm(min(itm+2,stpm))-tpm(itm+1)  ! Next averaging period (min to stop array out of bounds)
              par%tintm=max(par%tintm,tiny(0.d0))        ! to prevent par%tintm=0 after last output
          endif
      endif
    endif ! xmaster

    !!! Write global variables
    if (nglobalvar/=0) then
      if (any(abs(par%t-tpg) .le. 0.000001d0)) then
          itg=itg+1      
          do i = 1,noutnumbers
#ifdef USEMPI
            call space_collect_index(s,sl,outnumbers(i))
#endif
            if(xmaster) then
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
            endif ! xmaster
          enddo   ! outnumbers loop
      endif
    endif  ! end global file writing


    if(xmaster) then
      outputtimes=-999.d0
      outputtimes(1:itg)=tpg(1:itg)
      outputtimes(itg+1:itg+itp)=tpp(1:itp)
		outputtimes(itg+itp+1:itg+itp+itc)=tpc(1:itc)
      outputtimes(itg+itp+itc+1:itg+itp+itc+itm)=tpm(2:itm+1)          ! mean output always shifted by 1
	   outputtimes=outputtimes*max(par%morfac,1.d0)
      open(999,file='dims.dat',form='unformatted',access='direct',recl=wordsize*(10+size(outputtimes)))
      write(999,rec=1)		 itg*1.d0,&
							 s%nx*1.d0,&
							 s%ny*1.d0,&
                             s%ntheta*1.d0,&
                             par%kmax*1.d0,&
							 par%ngd*1.d0,&
							 par%nd*1.d0, &
                             itp*1.d0,&
							 itc*1.d0,&
							 itm*1.d0,&
							 outputtimes
      ! Just output for MICORE for backwards compat
!      open(999,file='dims.dat',form='unformatted',access='direct',recl=wordsize*(7+size(outputtimes)))
!      write(999,rec=1)itg*1.d0,&
!							 s%nx*1.d0,&
!							 s%ny*1.d0,&
!                      par%ngd*1.d0,&
!							 par%nd*1.d0, &
!                      itp*1.d0,&
!							 itm*1.d0,&
!							 outputtimes
!		close(999)
    endif  !xmaster

    ! Determine next time step
    t1=minval(tpg,MASK=tpg .gt. par%t)
    t2=minval(tpp,MASK=tpp .gt. par%t)
    t3=minval(tpm,MASK=tpm .gt. par%t)
	t4=minval(tpw,MASK=tpw .gt. par%t)
	t5=minval(tpc,MASK=tpc .gt. par%t)
    par%tnext=min(t1,t2,t3,t4,t5)
  end if  ! it > itpref

  !!! Close files

  if(xmaster) then
    if(par%t>=par%tstop) then
        do i=1,npoints+nrugauge
          close(indextopointsunit(i))
        enddo

        do i=1,nmeanvar
          close(indextomeanunit(i))
			 close(indextovarunit(i))
          close(indextominunit(i))
			 close(indextomaxunit(i))
        enddo

        do i=1,noutnumbers
          close (indextoglobalunit(i))
        enddo
    end if
  endif

end subroutine var_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!   INTERNAL SUBROUTINE   !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine makeintpvector(par,s,intpvector,mg,ng)

  use params
  use spaceparams

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
    if (mg .ge. s%is(i) .and. mg .le. s%is(i) + s%lm(i) .and. &
        ng .ge. s%js(i) .and. ng .le. s%js(i) + s%ln(i)) then
      m = mg - s%is(i) + 1  ! m and n now the indices for process p
      n = ng - s%js(i) + 1
      p = i-1
      exit
    endif
  enddo
  if (p .lt. 0) then
    write(*,*)'Catastrophic error in makeintpvector'
    write(*,*)'Cannot find processor for indexes',mg,ng
    call halt_program
  endif


  if (xmpi_rank .eq. p) then
#endif
    do j=1,numvars
      call indextos(s,j,t)
      value = -999   ! default
      select case (t%name)
        case(mnem_thetamean)
          value = 270-((s%thetamean(m,n)+s%alfa)*(180/par%px))  ! iwwvv thetamean: is that
                                                                ! different on each
                                                                ! process?
        case(mnem_Fx) 
          value = s%Fx(m,n)*cos(s%alfa)-s%Fy(m,n)*sin(s%alfa)
        case(mnem_Fy)
          value = s%Fx(m,n)*sin(s%alfa)+s%Fy(m,n)*cos(s%alfa)
        case(mnem_u)
          value = s%u(m,n)*cos(s%alfa)-s%v(m,n)*sin(s%alfa)
        case(mnem_v)
          value = s%u(m,n)*sin(s%alfa)+s%v(m,n)*cos(s%alfa)
        case(mnem_ue)
          value = s%ue(m,n)*cos(s%alfa)-s%ve(m,n)*sin(s%alfa)
        case(mnem_ve)
          value = s%ue(m,n)*sin(s%alfa)+s%ve(m,n)*cos(s%alfa)
        case(mnem_uwf)
          value = s%uwf(m,n)*cos(s%alfa)-s%vwf(m,n)*sin(s%alfa)
        case(mnem_vwf)
          value = s%uwf(m,n)*sin(s%alfa)+s%vwf(m,n)*cos(s%alfa)
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
                  case (mnem_yz,mnem_yv)
                    value = t%r1(n)
                  case (mnem_xz,mnem_xu)
                    value = t%r1(m)
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
  if (p .ne. xmpi_master) then
    if (xmpi_rank .eq. p) then
      call MPI_Send(intpvector, numvars, MPI_DOUBLE_PRECISION, xmpi_master, &
                    311, xmpi_comm,                    ierr)
    endif
    if (xmaster) then
      call MPI_Recv(intpvector, numvars, MPI_DOUBLE_PRECISION, p, &
                    311, xmpi_comm, MPI_STATUS_IGNORE, ierr)
    endif
  ! this barrier is really needed:
    call xmpi_barrier
  endif
#endif
end subroutine makeintpvector






subroutine outarray_r0(index,x)
  implicit none
  integer                       :: index
  real*8,intent(in)             :: x
  integer                       :: unit,reclen,jtg

  if(xmaster) then
    inquire(iolength=reclen) x
    call checkfile(index,unit,reclen,jtg)

    write(unit,rec=jtg) x
  endif

end subroutine outarray_r0

subroutine outarray_r1(index,x)
  implicit none
  integer                       :: index
  real*8, dimension(:), pointer :: x
  integer                       :: unit,reclen,jtg
  character*20                  :: mnem
  real*8                        :: pi

  pi = 4*atan(1.0d0)

  if(xmaster) then
    inquire(iolength=reclen) x
    call checkfile(index,unit,reclen,jtg)
  endif

  mnem = mnemonics(outnumbers(index))
  ! wwvv todo theta0 is scalar
  if (mnem .eq. 'theta' .or. &
      mnem .eq. 'theta0' ) then
    write(unit,rec=jtg) 270-(x*(180/pi))
  else
    write(unit,rec=jtg) x
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
  real*8                        :: pi

  pi = 4*atan(1.0d0)

  inquire(iolength=reclen) x
  call checkfile(index,unit,reclen,jtg)

  mnem = mnemonics(outnumbers(index))

  select case(mnem)
      case(mnem_thetamean)
        write(unit,rec=jtg)270-((x+s%alfa)*(180/pi))
      case(mnem_Fx)
        write(unit,rec=jtg)x*cos(s%alfa)-s%Fy*sin(s%alfa)
      case(mnem_Fy)
        write(unit,rec=jtg)s%Fx*sin(s%alfa)+x*cos(s%alfa)
      case(mnem_u)
        write(unit,rec=jtg)x*cos(s%alfa)-s%v*sin(s%alfa)
      case(mnem_v)
        write(unit,rec=jtg)s%u*sin(s%alfa)+x*cos(s%alfa)
      case(mnem_ue)
        write(unit,rec=jtg)x*cos(s%alfa)-s%ve*sin(s%alfa)
      case(mnem_ve)
        write(unit,rec=jtg)s%ue*sin(s%alfa)+x*cos(s%alfa)
      case(mnem_uwf)
        write(unit,rec=jtg)x*cos(s%alfa)-s%vwf*sin(s%alfa)
      case(mnem_vwf)
        write(unit,rec=jtg)s%uwf*sin(s%alfa)+x*cos(s%alfa)
      case default
    write(unit,rec=jtg) x
  end select

end subroutine outarray_r2

subroutine outarray_r3(s,index,x)
  use spaceparams
  implicit none
  type(spacepars),intent(in)    :: s
  integer index
  real*8, dimension(:,:,:), pointer :: x
  integer                       :: unit,reclen,jtg
  character*20                  :: mnem
  real*8                        :: pi

  pi=4*atan(1.0d0)

  inquire(iolength=reclen) x
  call checkfile(index,unit,reclen,jtg)

  mnem = mnemonics(outnumbers(index))
  select case(mnem)
    case(mnem_cgx)
      write(unit,rec=jtg)x*cos(s%alfa)-s%cgy*sin(s%alfa)
    case(mnem_cgy)
      write(unit,rec=jtg)s%cgx*sin(s%alfa)+x*cos(s%alfa)
    case(mnem_cx)
      write(unit,rec=jtg)x*cos(s%alfa)-s%cy*sin(s%alfa)
    case(mnem_cy)
      write(unit,rec=jtg)s%cx*sin(s%alfa)+x*cos(s%alfa)
    case(mnem_thet)
      write(unit,rec=jtg)270-((s%thet+s%alfa)*(180/pi))
    case(mnem_Susg)
      write(unit,rec=jtg)x*cos(s%alfa)-s%Svsg*sin(s%alfa)
    case(mnem_Svsg)
      write(unit,rec=jtg)s%Susg*sin(s%alfa)+x*cos(s%alfa)
	case(mnem_Subg)
      write(unit,rec=jtg)x*cos(s%alfa)-s%Svbg*sin(s%alfa)
    case(mnem_Svbg)
      write(unit,rec=jtg)s%Subg*sin(s%alfa)+x*cos(s%alfa)
    case default
      write(unit,rec=jtg) x
  end select
end subroutine outarray_r3

subroutine outarray_r4(index,x)
  implicit none
  integer index
  real*8, dimension(:,:,:,:), pointer :: x
  integer                       :: unit,reclen,jtg

  inquire(iolength=reclen) x
  call checkfile(index,unit,reclen,jtg)
  write(unit,rec=jtg) x
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

end module outputmod
