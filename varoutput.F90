!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!   MODULE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module outputmod
implicit none

! Robert: Add choice of output variables
logical,dimension(999)                          :: outputindex ! [-]      tracks which  global variables are to be outputted.
integer*4                                                       :: nglobalvar  ! number of global output variables
integer*4                                                       :: npoints     ! number of output points
integer*4                                                       :: nrugauge    ! number of runup gauges
integer*4,dimension(:),allocatable      :: pointtype   ! 0 = point output, 1 = runup gauge
integer*4,dimension(:),allocatable  :: xpoints     ! model x-coordinate of output points
integer*4,dimension(:),allocatable  :: ypoints     ! model y-coordinate of output points
integer*4,dimension(:),allocatable  :: nassocvar   ! vector with number of output variable per output point
integer*4,dimension(:,:),allocatable:: arrayassocvar ! Array with associated index of output variables per point
real*8,dimension(:,:,:), allocatable:: meanarrays  ! Keep time average variables
integer*4                                                       :: nmeanvar    ! number of time-average variables
integer*4,dimension(:),allocatable  :: meanvec     ! keep track of which mean variables are used
real*8,dimension(:),allocatable         :: tpg,tpp,tpm ! output time points
integer*4                                                       :: stpm            ! size of tpm





contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!   INITIALISE OUTPUT    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine init_output(s,par)
use params
use spaceparams
use readkey_module

IMPLICIT NONE

type(spacepars)                     :: s
type(parameters)                    :: par

integer                                                 :: id,ic,icold,i,ii,index
character(80)                                           :: line, keyword
character(80)                                           :: var, fname
integer, dimension(:,:),allocatable :: temparray
real*8                                                          :: tg1,tp1,tm1




!!!!! OUTPUT POINTS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


npoints  = readkey_int ('params.txt','npoints',      0,       0,     50)
nrugauge = readkey_int ('params.txt','nrugauge',     0,       0,     50)

if ((npoints+nrugauge)>0) then 
        allocate(pointtype(npoints+nrugauge))
        allocate(xpoints(npoints+nrugauge))
        allocate(ypoints(npoints+nrugauge))
        allocate(nassocvar(npoints+nrugauge))
        allocate(temparray(npoints+nrugauge,99))
        pointtype(1:npoints)=0
        pointtype(npoints+1:npoints+nrugauge)=1
        temparray=-1

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
                        read(10,*) xpoints(i),ypoints(i),nassocvar(i),line

                        icold=0
                        do ii =1,nassocvar(i)
                                ic=scan(line(icold+1:80),'#')
                                ic=ic+icold
                                var=line(icold+1:ic-1)
                                call chartoindex(var,index)

                                if (index/=-1) then
                                        temparray(i,ii)=index
                                else
                                        write(*,*)'Unknown point output type ',var
                                        stop
                                endif
                                icold=ic
                        enddo
                enddo
                close(10)
        endif

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
 		        	    xpoints(i)=1
 			            read(10,*) ypoints(i),nassocvar(i),line
                        icold=0
                        do ii =1,nassocvar(i)
                                ic=scan(line(icold+1:80),'#')
                                ic=ic+icold
                                var=line(icold+1:ic-1)
                                call chartoindex(var,index)

                                if (index/=-1) then
                                        temparray(i,ii)=index
                                else
                                        write(*,*)'Unknown point output type ',var
                                        stop
                                endif
                                icold=ic
                        enddo
                enddo
                close(10)
        endif


        ! Tidy up information
        allocate(arrayassocvar(npoints+nrugauge,maxval(nassocvar)))
        arrayassocvar(:,:)=temparray(:,1:maxval(nassocvar))
        deallocate (temparray)

endif





!!!!! GLOBAL OUTPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! outputindex is set in this subroutine. Row number corresponds to line in s.inp.
! Default choice = H,zs,zb,u,v,cc,Su,Sv, dims
! dims is given number 999 !!!!!

nglobalvar  = readkey_int ('params.txt','nglobalvar',   -1,       -1,     20)

outputindex = .false.

if (nglobalvar == 0) then                       ! User wants no output
        outputindex = .false.
elseif (nglobalvar == -1)       then    ! Default output
        outputindex(26)=.true.          ! H
        outputindex(43)=.true.          ! zs
        outputindex(44)=.true.          ! zs0
        outputindex(11)=.true.          ! zb
        outputindex(42)=.true.          ! hh
        outputindex(67)=.true.          ! u
        outputindex(68)=.true.          ! v
        outputindex(69)=.true.          ! ue
        outputindex(70)=.true.          ! ve
        outputindex(59)=.true.          ! urms
        outputindex(20)=.true.          ! Fx
        outputindex(21)=.true.          ! Fy
        outputindex(84)=.true.          ! cc
        outputindex(79)=.true.          ! ceq
        outputindex(81)=.true.          ! Su
        outputindex(82)=.true.          ! Sv
        outputindex(57)=.true.          ! E
        outputindex(58)=.true.          ! R
        outputindex(60)=.true.          ! D
        outputindex(90)=.true.          ! DR
        outputindex(999)=.true.         ! dims
elseif (nglobalvar == 999) then ! Output all
        outputindex = .true.
else                                                            ! User specified output
        ! Look for keyword nglobalvar in params.txt
        id=0
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
                call chartoindex(trim(line),index)

                if (index/=-1) then
                        outputindex(index)=.true.
                else
                        write(*,*)'Unknown global output type ',line
                        stop
                endif
        end do
        close(10)
endif


!!!!! TIME-AVEARGE ARRAYS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


nmeanvar = readkey_int ('params.txt','nmeanvar',0,0,15)
if (nmeanvar>0) then
        allocate(meanarrays(s%nx+1,s%ny+1,nmeanvar))
        allocate(meanvec(nmeanvar))
        id=0
        ! Look for keyword nmeanvar in params.txt
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
                call chartoindex(trim(line),index)
                if (index/=-1) then
                        meanvec(i)=index
                else
                        write(*,*)'Unknown time-average output type ',line
                        stop
                endif
        end do
        close(10)
meanarrays=0.d0
endif

!!!!! OUTPUT TIME POINTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! If we want global output then
if (nglobalvar/=0) then
        call readkey('params.txt','tsglobal',fname)
        if (fname/=' ') then
                open(10,file=fname)
                read(10,*)ii
                allocate(tpg(ii))
                do i=1,ii
                        read(10,*)tpg(i)
                enddo
                close(10)
        else
                ii=floor((par%tstop-par%tstart)/par%tintg)+1
                allocate(tpg(ii))
                do i=1,ii
                        tpg(i)=par%tstart+par%tintg*real(i-1)
                enddo
        endif
endif

! If we want point output then
if ((npoints+nrugauge)>0) then 
        call readkey('params.txt','tspoints',fname)
        if (fname/=' ') then
                open(10,file=fname)
                read(10,*)ii
                allocate(tpp(ii))
                do i=1,ii
                        read(10,*)tpp(i)
                enddo
                close(10)
        else
                ii=floor((par%tstop-par%tstart)/par%tintp)+1
                allocate(tpp(ii))
                do i=1,ii
                        tpp(i)=par%tstart+par%tintp*real(i-1)
                enddo
        endif
endif

if (nmeanvar>0) then
        call readkey('params.txt','tsmean',fname)
        if (fname/=' ') then
                open(10,file=fname)
                read(10,*)ii
                allocate(tpm(ii))
                do i=1,ii
                        read(10,*)tpm(i)
                enddo
                close(10)
        else
                ii=floor((par%tstop-par%tstart)/par%tintm)+1
                allocate(tpm(ii))
                do i=1,ii
                        tpm(i)=par%tstart+par%tintm*real(i-1)
                enddo
        endif
endif

! If tp series not defined, then no output wanted, so large timestep
if (.not. allocated(tpg)) then
        allocate(tpg(1))
        tpg=par%tstop
endif
if (.not. allocated(tpp)) then
        allocate(tpp(1))
        tpp=par%tstop
endif
if (.not. allocated(tpm)) then  ! Need minimum two in this array
        allocate(tpm(2))
        tpm(1)=par%tstop
        tpm(2)=par%tstop+1.d0
endif

tg1=minval(tpg,MASK=tpg .gt. par%t)
tp1=minval(tpp,MASK=tpp .gt. par%t)
tm1=minval(tpm,MASK=tpm .gt. par%t)
par%tnext=par%t+min(tg1,tp1,tm1)
par%tintm=tpm(2)-tpm(1)
stpm=size(tpm)
end subroutine init_output




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!   OUTPUT AT EVERY TIMESTEP    !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine var_output(it,s,par)
use params
use spaceparams

IMPLICIT NONE

type(spacepars)                         :: s
type(parameters)                        :: par

!logical,dimension(93)                              :: outputindex
integer                                 :: i,ii,tt
integer,save                            :: ndt,itg,itp,itm,day,ot,itprev
integer                                 :: it,i1,i2,i3
integer                                 :: wordsize
integer                                 :: reclen,reclen2,reclenp, idum
real*8                                                                  :: dtimestep,tpredicted,tnow, t1,t2,t3
!real*8,save                                                            :: told
real*8,dimension(10),save                               :: tlast10
character(12)                                                   :: fname
character(99)                                                   :: fname2
real*8,dimension(98)                                    :: intpvector
integer,dimension(:),allocatable                :: tempvectori
real*8,dimension(:),allocatable             :: tempvectorr, temp
integer,dimension(8)                                    :: datetime
real*8,dimension(size(tpg)+size(tpp)+size(tpm)) :: outputtimes


inquire(iolength=wordsize) 1.d0
reclen=wordsize*(s%nx+1)*(s%ny+1)
reclen2=wordsize*(s%nx+1)*(s%ny+1)*(par%ngd)*(par%nd)

! Do only at very first time step
if (par%t==par%dt) then
        ndt=0
        tlast10=0.d0
        day=0
        ot=0
        itprev=0
endif

ndt=ndt+1                                       ! Number of calculation time steps per output time step

if (nmeanvar/=0) then
        if (par%t>=tpm(1) .and. par%t<=tpm(stpm)) call makeaverage(s,par)       ! Make averages every flow timestep
endif


if (par%t>=(real(ot)*par%tstop/500.d0)) then
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
        write(*,fmt='(a,f5.1,a)')'Simulation ',100.d0*par%t/par%tstop,' percent complete'
        tpredicted=(par%tstop-par%t)/(par%tstop/500.d0)*dtimestep
        if (dtimestep<1 .and. tpredicted<120) then
                ! Write nothing
        elseif (tpredicted>=3600) then 
                write(*,fmt='(a,I3,a,I2,a)')'Time remaining ',floor(tpredicted/3600.0d0),' hours and ',&
                      nint((tpredicted-3600.0d0*floor(tpredicted/3600.0d0))/60.0d0),' minutes'
        elseif (tpredicted>=600) then
                write(*,fmt='(a,I2,a)')'Time remaining ',floor(tpredicted/60.0d0),' minutes'
        elseif (tpredicted>=60) then
                write(*,fmt='(a,I2,a,I2,a)')'Time remaining ',floor(tpredicted/60.0d0),' minutes and ',&
                        nint((tpredicted-60.0d0*floor(tpredicted/60.0d0))),' seconds'
        else
                write(*,fmt='(a,I2,a)')'Time remaining ',nint(tpredicted),' seconds'
        endif

endif



! Determine if this is an output timestep
if (it>itprev) then
        itprev=it


   !!!!!!! Start writing output


   if (it==1) then   !!! Only first time round: initialisation of files

                itp=0
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
                                reclenp=wordsize*(nassocvar(i)+1)*1
                                open(300+i,file=fname,form='unformatted',access='direct',recl=reclenp)
                        enddo
                endif


                !! First time file opening for time-average output
                itm=0
                if (nmeanvar>0) then
                        do i=1,nmeanvar
                                call makeaveragenames(meanvec(i),fname2)
                                fname2=trim(fname2)
                                open(400+i,file=fname2,form='unformatted',access='direct',recl=reclen)
                        enddo
                endif


                !! First time file opening for global output
                itg=0
                if (nglobalvar/=0) then
                        open(100,file='xy.dat',form='unformatted',access='direct',recl=reclen)
                        write(100,rec=1)s%xw
                        write(100,rec=2)s%yw
                        write(100,rec=3)s%x
                        write(100,rec=4)s%y
                        close(100)
                endif

        if (outputindex(1)) open(101,file='x.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(2)) open(102,file='y.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(3)) open(103,file='dx.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(4)) open(104,file='dy.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(5)) open(105,file='xz.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(6)) open(106,file='yz.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(7)) open(107,file='xu.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(8)) open(108,file='yv.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(9)) open(109,file='nx.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(10)) open(110,file='ny.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(11)) open(111,file='zb.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(12)) open(112,file='zb0.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(13)) open(113,file='theta.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(14)) open(114,file='ntheta.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(15)) open(115,file='dtheta.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(16)) open(116,file='theta0.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(17)) open(117,file='cxsth.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(18)) open(118,file='sxnth.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(19)) open(119,file='thetamean.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(20)) open(120,file='Fx.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(21)) open(121,file='Fy.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(22)) open(122,file='Sxy.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(23)) open(123,file='Syy.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(24)) open(124,file='Sxx.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(25)) open(125,file='n.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(26)) open(126,file='Hrms.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(27)) open(127,file='cgx.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(28)) open(128,file='cgy.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(29)) open(129,file='cx.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(30)) open(130,file='cy.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(31)) open(131,file='ctheta.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(32)) open(132,file='ee.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(33)) open(133,file='rr.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(34)) open(134,file='thet.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(35)) open(135,file='costhet.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(36)) open(136,file='sinthet.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(37)) open(137,file='sigt.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(38)) open(138,file='k.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(39)) open(139,file='c.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(40)) open(140,file='cg.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(41)) open(141,file='sigm.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(42)) open(142,file='hh.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(43)) open(143,file='zs.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(44)) open(144,file='zs0.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(45)) open(145,file='tideinpt.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(46)) open(146,file='tideinpz.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(47)) open(147,file='dzsdt.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(48)) open(148,file='dzbdt.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(49)) open(149,file='uu.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(50)) open(150,file='vv.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(51)) open(151,file='qx.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(52)) open(152,file='qy.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(53)) open(153,file='sedero.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(54)) open(154,file='dcdx.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(55)) open(155,file='dcdy.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(56)) open(156,file='ui.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(57)) open(157,file='E.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(58)) open(158,file='R.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(59)) open(159,file='urms.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(60)) open(160,file='D.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(61)) open(161,file='ust.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(62)) open(162,file='tm.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(63)) open(163,file='ueu.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(64)) open(164,file='vev.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(65)) open(165,file='vmagu.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(66)) open(166,file='vmagv.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(67)) open(167,file='u.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(68)) open(168,file='v.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(69)) open(169,file='ue.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(70)) open(170,file='ve.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(71)) open(171,file='hold.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(72)) open(172,file='wetu.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(73)) open(173,file='wetv.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(74)) open(174,file='wetz.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(75)) open(175,file='hu.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(76)) open(176,file='hv.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(77)) open(177,file='hum.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(78)) open(178,file='hvm.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(79)) open(179,file='ceq.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(80)) open(180,file='vmag.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(81)) open(181,file='Su.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(82)) open(182,file='Sv.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(83)) open(183,file='Ts.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(84)) open(184,file='cc.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(85)) open(185,file='uwf.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(86)) open(186,file='vwf.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(87)) open(187,file='ustr.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(88)) open(188,file='usd.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(89)) open(189,file='bi.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(90)) open(190,file='DR.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(91)) open(191,file='vardx.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(92)) open(192,file='vu.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(93)) open(193,file='Beta.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(94)) open(194,file='kb.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(95)) open(195,file='uon.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(96)) open(196,file='uoff.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(97)) open(197,file='Tbore.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(98)) open(198,file='dzav.dat',form='unformatted',access='direct',recl=reclen)
        if (outputindex(99)) open(199,file='maxzs.dat',form='unformatted',access='direct',recl=reclen)
        if (outputindex(100)) open(200,file='minzs.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(101)) open(201,file='uxbgrid.dat',form='unformatted',access='direct',recl=reclen)
                if (outputindex(102)) open(202,file='vxbgrid.dat',form='unformatted',access='direct',recl=reclen)  
   endif                !! End first time opening


        !!! Write at every output timestep


        !!! Write point variables
        !!! Only write if it is output time for points
        if (npoints+nrugauge>0) then
                if (any(abs(par%t-tpp) .le. 0.000001d0)) then
                        itp=itp+1
                        do i=1,npoints+nrugauge
                                !!! Make vector of all s% values at n,m grid coordinate
                                if (pointtype(i)==1) then
                                        if (.not.allocated(temp)) allocate (temp(1:s%nx+1))
                                        temp=(/(i,i=1,s%nx+1)/)
                                        idum = maxval(maxloc(temp*s%wetz(1:s%nx+1,ypoints(i)))) ! Picks the last "1" in temp
                                        call makeintpvector(par,s,intpvector,idum,ypoints(i))
                                else
                                        call makeintpvector(par,s,intpvector,xpoints(i),ypoints(i))
                                endif
                                allocate(tempvectori(nassocvar(i)))
                                allocate(tempvectorr(nassocvar(i)+1))
                                tempvectori=arrayassocvar(i,1:nassocvar(i))
                                tempvectorr(1)=par%t
                                do ii=1,nassocvar(i)
                                        tempvectorr(ii+1)=intpvector(tempvectori(ii))
                                enddo
                                write(300+i,rec=itp)tempvectorr
                                deallocate(tempvectori)
                                deallocate(tempvectorr)
                        enddo
                endif
        endif

        !!! Write average variables
        if (nmeanvar>0) then
                ! Not at the first in tpm as this is the start of averaging. Only output after second in tpm
                if (any(abs(par%t-tpm(2:stpm)) .le. 0.000001d0)) then
                        itm=itm+1
                        do i=1,nmeanvar
                                if (meanvec(i)==26) then                ! Hrms
                                        meanarrays(:,:,i)=sqrt(meanarrays(:,:,i))
                                        write(400+i,rec=itm)meanarrays(:,:,i)
                                elseif (meanvec(i)==59) then    ! urms
                                        meanarrays(:,:,i)=sqrt(meanarrays(:,:,i))
                                        write(400+i,rec=itm)meanarrays(:,:,i)
                                else                                                    ! non-rms variables
                                        write(400+i,rec=itm)meanarrays(:,:,i)
                                endif
                        enddo
                        meanarrays=0.0d0
                        par%tintm=tpm(min(itm+2,stpm))-tpm(itm+1)               ! Next averaging period (min to stop array out of bounds)
                        par%tintm=max(par%tintm,tiny(0.d0))                             ! to prevent par%tintm=0 after last output
                endif
        endif

        !!! Write global variables
        if (nglobalvar/=0) then
                if (any(abs(par%t-tpg) .le. 0.000001d0)) then
                itg=itg+1
                if (outputindex(1)) write(101,rec=itg)s%x
                if (outputindex(2)) write(102,rec=itg)s%y
                if (outputindex(3)) write(103,rec=itg)s%dx
                if (outputindex(4)) write(104,rec=itg)s%dy
                if (outputindex(5)) write(105,rec=itg)s%xz
                if (outputindex(6)) write(106,rec=itg)s%yz
                if (outputindex(7)) write(107,rec=itg)s%xu
                if (outputindex(8)) write(108,rec=itg)s%yv
                if (outputindex(9)) write(109,rec=itg)s%nx
                if (outputindex(10)) write(110,rec=itg)s%ny
                if (outputindex(11)) write(111,rec=itg)s%zb
                if (outputindex(12)) write(112,rec=itg)s%zb0
                if (outputindex(13)) write(113,rec=itg)270-(s%theta*(180/par%px))
                if (outputindex(14)) write(114,rec=itg)s%ntheta
                if (outputindex(15)) write(115,rec=itg)s%dtheta
                if (outputindex(16)) write(116,rec=itg)270-(s%theta0*(180/par%px))
                if (outputindex(17)) write(117,rec=itg)s%cxsth
                if (outputindex(18)) write(118,rec=itg)s%sxnth
                if (outputindex(19)) write(119,rec=itg)270-((s%thetamean+s%alfa)*(180/par%px))
                if (outputindex(20)) write(120,rec=itg)s%Fx*cos(s%alfa)-s%Fy*sin(s%alfa)
                if (outputindex(21)) write(121,rec=itg)s%Fx*sin(s%alfa)+s%Fy*cos(s%alfa)
                if (outputindex(22)) write(122,rec=itg)s%Sxy
                if (outputindex(23)) write(123,rec=itg)s%Syy
                if (outputindex(24)) write(124,rec=itg)s%Sxx
                if (outputindex(25)) write(125,rec=itg)s%n
                if (outputindex(26)) write(126,rec=itg)s%H
                if (outputindex(27)) write(127,rec=itg)s%cgx*cos(s%alfa)-s%cgy*sin(s%alfa)
                if (outputindex(28)) write(128,rec=itg)s%cgx*sin(s%alfa)+s%cgy*cos(s%alfa)
                if (outputindex(29)) write(129,rec=itg)s%cx*cos(s%alfa)-s%cy*sin(s%alfa)
                if (outputindex(30)) write(130,rec=itg)s%cx*sin(s%alfa)+s%cy*cos(s%alfa)
                if (outputindex(31)) write(131,rec=itg)s%ctheta
                if (outputindex(32)) write(132,rec=itg)s%ee
                if (outputindex(33)) write(133,rec=itg)s%rr
                if (outputindex(34)) write(134,rec=itg)270-((s%thet+s%alfa)*(180/par%px))
                if (outputindex(35)) write(135,rec=itg)s%costhet
                if (outputindex(36)) write(136,rec=itg)s%sinthet
                if (outputindex(37)) write(137,rec=itg)s%sigt
                if (outputindex(38)) write(138,rec=itg)s%k
                if (outputindex(39)) write(139,rec=itg)s%c
                if (outputindex(40)) write(140,rec=itg)s%cg
                if (outputindex(41)) write(141,rec=itg)s%sigm
                if (outputindex(42)) write(142,rec=itg)s%hh
                if (outputindex(43)) write(143,rec=itg)s%zs
                if (outputindex(44)) write(144,rec=itg)s%zs0
                if (outputindex(45)) write(145,rec=itg)s%tideinpt
                if (outputindex(46)) write(146,rec=itg)s%tideinpz
                if (outputindex(47)) write(147,rec=itg)s%dzsdt
                if (outputindex(48)) write(148,rec=itg)s%dzbdt
                if (outputindex(49)) write(149,rec=itg)s%uu
                if (outputindex(50)) write(150,rec=itg)s%vv
                if (outputindex(51)) write(151,rec=itg)s%qx
                if (outputindex(52)) write(152,rec=itg)s%qy
                if (outputindex(53)) write(153,rec=itg)s%sedero
                if (outputindex(54)) write(154,rec=itg)s%dcdx
                if (outputindex(55)) write(155,rec=itg)s%dcdy
                if (outputindex(56)) write(156,rec=itg)s%ui
                if (outputindex(57)) write(157,rec=itg)s%E
                if (outputindex(58)) write(158,rec=itg)s%R
                if (outputindex(59)) write(159,rec=itg)s%urms
                if (outputindex(60)) write(160,rec=itg)s%D
                if (outputindex(61)) write(161,rec=itg)s%ust
                if (outputindex(62)) write(162,rec=itg)s%tm
                if (outputindex(63)) write(163,rec=itg)s%ueu
                if (outputindex(64)) write(164,rec=itg)s%vev
                if (outputindex(65)) write(165,rec=itg)s%vmagu
                if (outputindex(66)) write(166,rec=itg)s%vmagv
                if (outputindex(67)) write(167,rec=itg)s%u*cos(s%alfa)-s%v*sin(s%alfa)
                if (outputindex(68)) write(168,rec=itg)s%u*sin(s%alfa)+s%v*cos(s%alfa)
                if (outputindex(69)) write(169,rec=itg)s%ue*cos(s%alfa)-s%ve*sin(s%alfa)
                if (outputindex(70)) write(170,rec=itg)s%ue*sin(s%alfa)+s%ve*cos(s%alfa)
                if (outputindex(71)) write(171,rec=itg)s%hold
                if (outputindex(72)) write(172,rec=itg)s%wetu
                if (outputindex(73)) write(173,rec=itg)s%wetv
                if (outputindex(74)) write(174,rec=itg)s%wetz
                if (outputindex(75)) write(175,rec=itg)s%hu
                if (outputindex(76)) write(176,rec=itg)s%hv
                if (outputindex(77)) write(177,rec=itg)s%hum
                if (outputindex(78)) write(178,rec=itg)s%hvm
                if (outputindex(79)) write(179,rec=itg)s%ceqg
                if (outputindex(80)) write(180,rec=itg)s%vmag
                if (outputindex(81)) write(181,rec=itg)s%Sug*cos(s%alfa)-s%Svg*sin(s%alfa)
                if (outputindex(82)) write(182,rec=itg)s%Sug*sin(s%alfa)+s%Svg*cos(s%alfa)
                if (outputindex(83)) write(183,rec=itg)s%Tsg
                if (outputindex(84)) write(184,rec=itg)s%ccg
                if (outputindex(85)) write(185,rec=itg)s%uwf*cos(s%alfa)-s%vwf*sin(s%alfa)
                if (outputindex(86)) write(186,rec=itg)s%uwf*sin(s%alfa)+s%vwf*cos(s%alfa)
                if (outputindex(87)) write(187,rec=itg)s%ustr
                if (outputindex(88)) write(188,rec=itg)s%usd
                if (outputindex(89)) write(189,rec=itg)s%bi
                if (outputindex(90)) write(190,rec=itg)s%DR
                if (outputindex(91)) write(191,rec=itg)s%vardx
                if (outputindex(92)) write(192,rec=itg)s%vu
                if (outputindex(93)) write(193,rec=itg)s%BR
                if (outputindex(94)) write(194,rec=itg)s%kb
                if (outputindex(95)) write(195,rec=itg)s%uon
                if (outputindex(96)) write(196,rec=itg)s%uoff
                if (outputindex(97)) write(197,rec=itg)s%Tbore
                if (outputindex(98)) write(198,rec=itg)s%dzav
                if (outputindex(99)) write(199,rec=itg)s%maxzs
                if (outputindex(100)) write(200,rec=itg)s%minzs
                if (outputindex(101)) write(201,rec=itg)s%u
                if (outputindex(102)) write(202,rec=itg)s%v

          endif
        endif  ! end global file writing


outputtimes=-999.d0
outputtimes(1:itg)=tpg(1:itg)
outputtimes(itg+1:itg+itp)=tpp(1:itp)
outputtimes(itg+itp+1:itg+itp+itm)=tpm(2:itm+1)                 ! mean output always shifted by 1
open(999,file='dims.dat',form='unformatted',access='direct',recl=wordsize*(7+size(outputtimes)))
                write(999,rec=1)itg*1.d0,s%nx*1.d0,s%ny*1.d0,par%ngd*1.d0,par%nd*1.d0, &
                                    itp*1.d0,itm*1.d0,outputtimes
close(999)




! Determine next time step
t1=minval(tpg,MASK=tpg .gt. par%t)
t2=minval(tpp,MASK=tpp .gt. par%t)
t3=minval(tpm,MASK=tpm .gt. par%t)
par%tnext=min(t1,t2,t3)
end if

!!! Close files


if(par%t>=par%tstop) then

        if (npoints+nrugauge>0) then
                do i=1,npoints+nrugauge
                        close(300+i)
                enddo
        endif

        if (nmeanvar>0) then
                do i=1,nmeanvar
                        close(400+i)
                enddo
        endif


   if (outputindex(1)) close(101)
        if (outputindex(2)) close(102)
        if (outputindex(3)) close(103)
        if (outputindex(4)) close(104)
        if (outputindex(5)) close(105)
        if (outputindex(6)) close(106)
        if (outputindex(7)) close(107)
        if (outputindex(8)) close(108)
        if (outputindex(9)) close(109)
        if (outputindex(10)) close(110)
        if (outputindex(11)) close(111)
        if (outputindex(12)) close(112)
        if (outputindex(13)) close(113)
        if (outputindex(14)) close(114)
        if (outputindex(15)) close(115)
        if (outputindex(16)) close(116)
        if (outputindex(17)) close(117)
        if (outputindex(18)) close(118)
        if (outputindex(19)) close(119)
        if (outputindex(20)) close(120)
        if (outputindex(21)) close(121)
        if (outputindex(22)) close(122)
        if (outputindex(23)) close(123)
        if (outputindex(24)) close(124)
        if (outputindex(25)) close(125)
        if (outputindex(26)) close(126)
        if (outputindex(27)) close(127)
        if (outputindex(28)) close(128)
        if (outputindex(29)) close(129)
        if (outputindex(30)) close(130)
        if (outputindex(31)) close(131)
        if (outputindex(32)) close(132)
        if (outputindex(33)) close(133)
        if (outputindex(34)) close(134)
        if (outputindex(35)) close(135)
        if (outputindex(36)) close(136)
        if (outputindex(37)) close(137)
        if (outputindex(38)) close(138)
        if (outputindex(39)) close(139)
        if (outputindex(40)) close(140)
        if (outputindex(41)) close(141)
        if (outputindex(42)) close(142)
        if (outputindex(43)) close(143)
        if (outputindex(44)) close(144)
        if (outputindex(45)) close(145)
        if (outputindex(46)) close(146)
        if (outputindex(47)) close(147)
        if (outputindex(48)) close(148)
        if (outputindex(49)) close(149)
        if (outputindex(50)) close(150)
        if (outputindex(51)) close(151)
        if (outputindex(52)) close(152)
        if (outputindex(53)) close(153)
        if (outputindex(54)) close(154)
        if (outputindex(55)) close(155)
        if (outputindex(56)) close(156)
        if (outputindex(57)) close(157)
        if (outputindex(58)) close(158)
        if (outputindex(59)) close(159)
        if (outputindex(60)) close(160)
        if (outputindex(61)) close(161)
        if (outputindex(62)) close(162)
        if (outputindex(63)) close(163)
        if (outputindex(64)) close(164)
        if (outputindex(65)) close(165)
        if (outputindex(66)) close(166)
        if (outputindex(67)) close(167)
        if (outputindex(68)) close(168)
        if (outputindex(69)) close(169)
        if (outputindex(70)) close(170)
        if (outputindex(71)) close(171)
        if (outputindex(72)) close(172)
        if (outputindex(73)) close(173)
        if (outputindex(74)) close(174)
        if (outputindex(75)) close(175)
        if (outputindex(76)) close(176)
        if (outputindex(77)) close(177)
        if (outputindex(78)) close(178)
        if (outputindex(79)) close(179)
        if (outputindex(80)) close(180)
        if (outputindex(81)) close(181)
        if (outputindex(82)) close(182)
        if (outputindex(83)) close(183)
        if (outputindex(84)) close(184)
        if (outputindex(85)) close(185)
        if (outputindex(86)) close(186)
        if (outputindex(87)) close(187)
        if (outputindex(88)) close(188)
        if (outputindex(89)) close(189)
        if (outputindex(90)) close(190)
        if (outputindex(91)) close(191)
        if (outputindex(92)) close(192)
        if (outputindex(93)) close(193)
        if (outputindex(94)) close(194)
        if (outputindex(95)) close(195)
        if (outputindex(96)) close(196)
        if (outputindex(97)) close(197)
        if (outputindex(98)) close(198)
    if (outputindex(99)) close(199)
    if (outputindex(100)) close(200)
        if (outputindex(101)) close(201)
        if (outputindex(102)) close(202)
end if


end subroutine var_output





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!   INTERNAL SUBROUTINE   !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine chartoindex(line,index)

IMPLICIT NONE

character(len=*)                :: line
integer                                 :: index

        index=-1
                if (line == 'x')   index=1
                if (line == 'y')   index=2
                if (line == 'dx')  index=3
                if (line == 'dy')  index=4
                if (line == 'xz')  index=5
                if (line == 'yz')  index=6
                if (line == 'xu')  index=7
                if (line == 'yu')  index=8
                if (line == 'nx')  index=9
                if (line == 'ny')  index=10
                if (line == 'zb')  index=11
                if (line == 'zb0')  index=12
                if (line == 'theta')  index=13
                if (line == 'ntheta')  index=14
                if (line == 'dtheta')  index=15
                if (line == 'theta0')  index=16
                if (line == 'cxsth')  index=17
                if (line == 'sxnth')  index=18
                if (line == 'thetamean')  index=19
                if (line == 'Fx')  index=20
                if (line == 'Fy')  index=21
                if (line == 'Sxy')  index=22
                if (line == 'Syy')  index=23
                if (line == 'Sxx')  index=24
                if (line == 'n')  index=25
                if (line == 'H')  index=26
                if (line == 'cgx')  index=27
                if (line == 'cgy')  index=28
                if (line == 'cx')  index=29
                if (line == 'cy')  index=30
                if (line == 'ctheta')  index=31
                if (line == 'ee')  index=32
                if (line == 'rr')  index=33
                if (line == 'thet')  index=34
                if (line == 'costhet')  index=35
                if (line == 'sinthet')  index=36
                if (line == 'sigt')  index=37
                if (line == 'k')  index=38
                if (line == 'c')  index=39
                if (line == 'cg')  index=40
                if (line == 'sigm')  index=41
                if (line == 'hh')  index=42
                if (line == 'zs')  index=43
                if (line == 'zs0')  index=44
                if (line == 'tideinpt')  index=45
                if (line == 'tideinpz')  index=46
                if (line == 'dzsdt')  index=47
                if (line == 'dzbdt')  index=48
                if (line == 'uu')  index=49
                if (line == 'vv')  index=50
                if (line == 'qx')  index=51
                if (line == 'qy')  index=52
                if (line == 'sedero')  index=53
                if (line == 'dcdx')  index=54
                if (line == 'dcdy')  index=55
                if (line == 'ui')  index=56
                if (line == 'E')  index=57
                if (line == 'R')  index=58
                if (line == 'urms')  index=59
                if (line == 'D')  index=60
                if (line == 'ust')  index=61
                if (line == 'tm')  index=62
                if (line == 'ueu')  index=63
                if (line == 'vev')  index=64
                if (line == 'vmagu')  index=65
                if (line == 'vmagv')  index=66
                if (line == 'u')  index=67
                if (line == 'v')  index=68
                if (line == 'ue')  index=69
                if (line == 've')  index=70
                if (line == 'hold')  index=71
                if (line == 'wetu')  index=72
                if (line == 'wetv')  index=73
                if (line == 'wetz')  index=74
                if (line == 'hu')  index=75
                if (line == 'hv')  index=76
                if (line == 'hum')  index=77
                if (line == 'hvm')  index=78
                if (line == 'ceq')  index=79
                if (line == 'vmag')  index=80
                if (line == 'Su')  index=81
                if (line == 'Sv')  index=82
                if (line == 'Ts')  index=83
                if (line == 'cc')  index=84
                if (line == 'uwf')  index=85
                if (line == 'vwf')  index=86
                if (line == 'ustr')  index=87
                if (line == 'usd')  index=88
                if (line == 'bi')  index=89
                if (line == 'DR')  index=90
                if (line == 'vardx')  index=91
                if (line == 'vu')  index=92
                if (line == 'Beta')  index=93
                if (line == 'kb')  index=94
                if (line == 'uon')  index=95
                if (line == 'uoff')  index=96
                if (line == 'Tbore')  index=97
                if (line == 'dzav')  index=98
                if (line == 'maxzs')  index=99
                if (line == 'minzs')  index=100
                if (line == 'uxbgrid') index=101
                if (line == 'vxbgrid') index=102
                if (line == 'dims')  index=999

end subroutine chartoindex


subroutine makeintpvector(par,s,intpvector,m,n)

use params
use spaceparams

IMPLICIT NONE

type(parameters), intent(IN)            :: par
type(spacepars), intent(IN)             :: s
real*8,dimension(98)                                    :: intpvector
integer                                                                 :: m,n

intpvector(1    )=              s%x(m,n)
intpvector(2    )=              s%y     (m,n)
intpvector(3    )=              s%dx
intpvector(4    )=              s%dy
intpvector(5    )=              s%xz    (m)
intpvector(6    )=              s%yz    (n)
intpvector(7    )=              s%xu    (m)
intpvector(8    )=              s%yv    (n)
intpvector(9    )=              real(s%nx)
intpvector(10   )=              real(s%ny)
intpvector(11   )=              s%zb    (m,n)
intpvector(12   )=              s%zb0   (m,n)
intpvector(13   )=              -999.d0
intpvector(14   )=              real(s%ntheta)
intpvector(15   )=              s%dtheta
intpvector(16   )=              270-(s%theta0*(180/par%px))
intpvector(17   )=              -999.d0
intpvector(18   )=              -999.d0
intpvector(19   )=              270-((s%thetamean(m,n)+s%alfa)*(180/par%px))
intpvector(20   )=              s%Fx(m,n)*cos(s%alfa)-s%Fy(m,n)*sin(s%alfa)
intpvector(21   )=              s%Fx(m,n)*sin(s%alfa)+s%Fy(m,n)*cos(s%alfa)
intpvector(22   )=              s%Sxy   (m,n)
intpvector(23   )=              s%Syy   (m,n)
intpvector(24   )=              s%Sxx   (m,n)
intpvector(25   )=              s%n     (m,n)
intpvector(26   )=              s%H     (m,n)
intpvector(27   )=              -999.d0
intpvector(28   )=              -999.d0
intpvector(29   )=              -999.d0
intpvector(30   )=              -999.d0
intpvector(31   )=              -999.d0
intpvector(32   )=              -999.d0
intpvector(33   )=              -999.d0
intpvector(34   )=              -999.d0
intpvector(35   )=              -999.d0
intpvector(36   )=              -999.d0
intpvector(37   )=              -999.d0
intpvector(38   )=              s%k     (m,n)
intpvector(39   )=              s%c     (m,n)
intpvector(40   )=              s%cg    (m,n)
intpvector(41   )=              s%sigm  (m,n)
intpvector(42   )=              s%hh    (m,n)
intpvector(43   )=              s%zs    (m,n)
intpvector(44   )=              s%zs0   (m,n)
intpvector(45   )=              -999.d0
intpvector(46   )=          -999.d0
intpvector(47   )=              s%dzsdt         (m,n)
intpvector(48   )=              s%dzbdt         (m,n)
intpvector(49   )=              s%uu            (m,n)
intpvector(50   )=              s%vv            (m,n)
intpvector(51   )=              s%qx            (m,n)
intpvector(52   )=              s%qy            (m,n)
intpvector(53   )=              s%sedero        (m,n)
intpvector(54   )=              s%dcdx          (m,n)
intpvector(55   )=              s%dcdy          (m,n)
intpvector(56   )=              s%ui            (m,n)
intpvector(57   )=              s%E             (m,n)
intpvector(58   )=              s%R             (m,n)
intpvector(59   )=              s%urms          (m,n)
intpvector(60   )=              s%D             (m,n)
intpvector(61   )=              s%ust           (m,n)
intpvector(62   )=              s%tm            (m,n)
intpvector(63   )=              s%ueu           (m,n)
intpvector(64   )=              s%vev           (m,n)
intpvector(65   )=              s%vmagu         (m,n)
intpvector(66   )=              s%vmagv         (m,n)
intpvector(67   )=              s%u(m,n)*cos(s%alfa)-s%v(m,n)*sin(s%alfa)
intpvector(68   )=              s%u(m,n)*sin(s%alfa)+s%v(m,n)*cos(s%alfa)
intpvector(69   )=              s%ue(m,n)*cos(s%alfa)-s%ve(m,n)*sin(s%alfa)
intpvector(70   )=              s%ue(m,n)*sin(s%alfa)+s%ve(m,n)*cos(s%alfa)
intpvector(71   )=              s%hold          (m,n)
intpvector(72   )=              s%wetu          (m,n)
intpvector(73   )=              s%wetv          (m,n)
intpvector(74   )=              s%wetz          (m,n)
intpvector(75   )=              s%hu            (m,n)
intpvector(76   )=              s%hv            (m,n)
intpvector(77   )=              s%hum           (m,n)
intpvector(78   )=              s%hvm           (m,n)
intpvector(79   )=              -999.d0
intpvector(80   )=              s%vmag          (m,n)
intpvector(81   )=              -999.d0
intpvector(82   )=              -999.d0
intpvector(83   )=              -999.d0
intpvector(84   )=              -999.d0
intpvector(85   )=              s%uwf(m,n)*cos(s%alfa)-s%vwf(m,n)*sin(s%alfa)
intpvector(86   )=              s%uwf(m,n)*sin(s%alfa)+s%vwf(m,n)*cos(s%alfa)
intpvector(87   )=              s%ustr          (m,n)
intpvector(88   )=              s%usd   (m,n)
intpvector(89   )=              -999.d0
intpvector(90   )=              s%DR    (m,n)
intpvector(91   )=              -999.d0
intpvector(92   )=              s%vu    (m,n)
intpvector(93   )=              s%BR    (m,n)
intpvector(94   )=              s%kb    (m,n)
intpvector(95   )=              s%uon   (m,n)
intpvector(96   )=              s%uoff  (m,n)
intpvector(97   )=              s%Tbore (m,n)
intpvector(98   )=              s%dzav  (m,n)

end subroutine makeintpvector

!!!! Routine to calculate time-average variables

subroutine makeaverage(s,par)

use params
use spaceparams

IMPLICIT NONE

type(parameters), intent(IN)            :: par
type(spacepars), intent(IN)             :: s
integer                                                                 :: i

do i=1,nmeanvar

        if (meanvec(i)==1 ) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%x
        if (meanvec(i)==2 ) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%y
        if (meanvec(i)==3 ) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%dx
        if (meanvec(i)==4 ) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%dy
        if (meanvec(i)==11) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%zb
        if (meanvec(i)==12) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%zb0
        if (meanvec(i)==16) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*(270-(s%theta0*(180/par%px)))
        if (meanvec(i)==19) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*(270-((s%thetamean+s%alfa) *(180/par%px)))
        if (meanvec(i)==20) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*(s%Fx *cos(s%alfa)-s%Fy *sin(s%alfa))
        if (meanvec(i)==21) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*(s%Fx *sin(s%alfa)+s%Fy *cos(s%alfa))
        if (meanvec(i)==22) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%Sxy
        if (meanvec(i)==23) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%Syy
        if (meanvec(i)==24) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%Sxx
        if (meanvec(i)==25) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%n
        if (meanvec(i)==26) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%H**2
        if (meanvec(i)==38) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%k
        if (meanvec(i)==39) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%c
        if (meanvec(i)==40) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%cg
        if (meanvec(i)==41) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%sigm
        if (meanvec(i)==42) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%hh
        if (meanvec(i)==43) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%zs
        if (meanvec(i)==44) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%zs0
        if (meanvec(i)==47) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%dzsdt
        if (meanvec(i)==48) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%dzbdt
        if (meanvec(i)==49) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%uu
        if (meanvec(i)==50) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%vv
        if (meanvec(i)==51) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%qx
        if (meanvec(i)==52) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%qy
        if (meanvec(i)==53) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%sedero
        if (meanvec(i)==54) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%dcdx
        if (meanvec(i)==55) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%dcdy
        if (meanvec(i)==56) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%ui
        if (meanvec(i)==57) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%E
        if (meanvec(i)==58) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%R
        if (meanvec(i)==59) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%urms**2
        if (meanvec(i)==60) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%D
        if (meanvec(i)==61) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%ust
        if (meanvec(i)==62) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%tm
        if (meanvec(i)==63) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%ueu
        if (meanvec(i)==64) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%vev
        if (meanvec(i)==65) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%vmagu
        if (meanvec(i)==66) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%vmagv
        if (meanvec(i)==67) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*(s%u *cos(s%alfa)-s%v *sin(s%alfa))
        if (meanvec(i)==68) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*(s%u *sin(s%alfa)+s%v *cos(s%alfa))
        if (meanvec(i)==69) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*(s%ue *cos(s%alfa)-s%ve *sin(s%alfa))
        if (meanvec(i)==70) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*(s%ue *sin(s%alfa)+s%ve *cos(s%alfa))
        if (meanvec(i)==71) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%hold
        if (meanvec(i)==72) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%wetu
        if (meanvec(i)==73) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%wetv
        if (meanvec(i)==74) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%wetz
        if (meanvec(i)==75) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%hu
        if (meanvec(i)==76) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%hv
        if (meanvec(i)==77) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%hum
        if (meanvec(i)==78) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%hvm
        if (meanvec(i)==80) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%vmag
        if (meanvec(i)==85) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*(s%uwf *cos(s%alfa)-s%vwf *sin(s%alfa))
        if (meanvec(i)==86) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*(s%uwf *sin(s%alfa)+s%vwf *cos(s%alfa))
        if (meanvec(i)==87) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%ustr
        if (meanvec(i)==88) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%usd
        if (meanvec(i)==90) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%DR
        if (meanvec(i)==92) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%vu
        if (meanvec(i)==93) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%BR
        if (meanvec(i)==94) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%kb
        if (meanvec(i)==95) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%uon
        if (meanvec(i)==96) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%uoff
        if (meanvec(i)==97) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%Tbore
        if (meanvec(i)==98) meanarrays(:,:,i)=  meanarrays(:,:,i) + (par%dt/par%tintm)*s%dzav
enddo

end subroutine makeaverage

! Subroutine make time series file name

subroutine makeaveragenames(counter,name)


IMPLICIT NONE

character(99)                                                   :: name
integer                                                                 :: counter

name='unknown_mean.dat'

if (counter==1 ) name ='x_mean.dat'
if (counter==2 ) name ='y_mean.dat'
if (counter==3 ) name ='dx_mean.dat'
if (counter==4 ) name ='dy_mean.dat'
if (counter==11) name ='zb_mean.dat'
if (counter==12) name ='zb0_mean.dat'
if (counter==16) name ='theta0_mean.dat'
if (counter==19) name ='thetamean_mean.dat'
if (counter==20) name ='Fx_mean.dat'
if (counter==21) name ='Fy_mean.dat'
if (counter==22) name ='Sxy_mean.dat'
if (counter==23) name ='Syy_mean.dat'
if (counter==24) name ='Sxx_mean.dat'
if (counter==25) name ='n_mean.dat'
if (counter==26) name ='H_mean.dat'
if (counter==38) name ='k_mean.dat'
if (counter==39) name ='c_mean.dat'
if (counter==40) name ='cg_mean.dat'
if (counter==41) name ='sigm_mean.dat'
if (counter==42) name ='hh_mean.dat'
if (counter==43) name ='zs_mean.dat'
if (counter==44) name ='zs0_mean.dat'
if (counter==47) name ='dzsdt_mean.dat'
if (counter==48) name ='dzbdt_mean.dat'
if (counter==49) name ='uu_mean.dat'
if (counter==50) name ='vv_mean.dat'
if (counter==51) name ='qx_mean.dat'
if (counter==52) name ='qy_mean.dat'
if (counter==53) name ='sedero_mean.dat'
if (counter==54) name ='dcdx_mean.dat'
if (counter==55) name ='dcdy_mean.dat'
if (counter==56) name ='ui_mean.dat'
if (counter==57) name ='E_mean.dat'
if (counter==58) name ='R_mean.dat'
if (counter==59) name ='urms_mean.dat'
if (counter==60) name ='D_mean.dat'
if (counter==61) name ='ust_mean.dat'
if (counter==62) name ='tm_mean.dat'
if (counter==63) name ='ueu_mean.dat'
if (counter==64) name ='vev_mean.dat'
if (counter==65) name ='vmagu_mean.dat'
if (counter==66) name ='vmagv_mean.dat'
if (counter==67) name ='u_mean.dat'
if (counter==68) name ='v_mean.dat'
if (counter==69) name ='ue_mean.dat'
if (counter==70) name ='ve_mean.dat'
if (counter==71) name ='hold_mean.dat'
if (counter==72) name ='wetu_mean.dat'
if (counter==73) name ='wetv_mean.dat'
if (counter==74) name ='wetz_mean.dat'
if (counter==75) name ='hu_mean.dat'
if (counter==76) name ='hv_mean.dat'
if (counter==77) name ='hum_mean.dat'
if (counter==78) name ='hvm_mean.dat'
if (counter==80) name ='vmag_mean.dat'
if (counter==85) name ='uwf_mean.dat'
if (counter==86) name ='vwf_mean.dat'
if (counter==87) name ='ustr_mean.dat'
if (counter==88) name ='usd_mean.dat'
if (counter==90) name ='DR_mean.dat'
if (counter==92) name ='vu_mean.dat'
if (counter==93) name ='BR_mean.dat'
if (counter==94) name ='kb_mean.dat'
if (counter==95) name ='uon_mean.dat'
if (counter==96) name ='uoff_mean.dat'
if (counter==97) name ='Tbore_mean.dat'
if (counter==98) name ='dzav_mean.dat'

end subroutine makeaveragenames

end module outputmod
