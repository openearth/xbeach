subroutine output(it,s,par)
use params
use spaceparams

IMPLICIT NONE

type(spacepars)                         :: s
type(parameters)                        :: par

integer                                 :: i,nt
integer,save                            :: ndt
integer                                 :: j,ig
integer                                 :: it
!integer                                 :: wordsize=8   ! Gfortran under Cygwin
integer                                 :: wordsize=2   ! Compaq/Windows or Ifort Linux
integer                                 :: reclen,reclen2

reclen=wordsize*(s%nx+1)*(s%ny+1)
reclen2=wordsize*(s%nx+1)*(s%ny+1)*(par%ngd)*(par%nd)

if (par%t==par%dt) ndt=0
ndt=ndt+1					! Number of time steps per second
if (MODULO(par%t,par%tint)==0) then
    if(par%t>=par%tstart) then
       it=it+1
       if (it==1) then
          open(100,file='xy.dat',form='unformatted',access='direct',recl=reclen)
          write(100,rec=1)s%xw
          write(100,rec=2)s%yw
          write(100,rec=3)s%x
          write(100,rec=4)s%y
          close(100)
          open(101,file='zs.dat'  ,form='unformatted',access='direct',recl=reclen)
          open(102,file='u.dat'   ,form='unformatted',access='direct',recl=reclen)
          open(103,file='v.dat'   ,form='unformatted',access='direct',recl=reclen)
          open(104,file='ue.dat'  ,form='unformatted',access='direct',recl=reclen)
          open(105,file='ve.dat'  ,form='unformatted',access='direct',recl=reclen)
          open(106,file='Hrms.dat',form='unformatted',access='direct',recl=reclen)
          open(107,file='urms.dat',form='unformatted',access='direct',recl=reclen)
          open(108,file='zb.dat'  ,form='unformatted',access='direct',recl=reclen)
          open(109,file='hh.dat'  ,form='unformatted',access='direct',recl=reclen)
          open(110,file='Fx.dat'  ,form='unformatted',access='direct',recl=reclen)
          open(111,file='Fy.dat'  ,form='unformatted',access='direct',recl=reclen)
          open(112,file='E.dat'   ,form='unformatted',access='direct',recl=reclen)
          open(113,file='R.dat'   ,form='unformatted',access='direct',recl=reclen)
          open(114,file='D.dat'   ,form='unformatted',access='direct',recl=reclen)
          open(115,file='cc.dat'  ,form='unformatted',access='direct',recl=reclen)
          open(116,file='Su.dat'  ,form='unformatted',access='direct',recl=reclen)
          open(117,file='Sv.dat'  ,form='unformatted',access='direct',recl=reclen)
          open(118,file='Gd.dat'  ,form='unformatted',access='direct',recl=reclen2)
       endif
       write(*,*) 'it=', it,' used ',ndt,' time steps'
	   ndt=0
       write(101,rec=it)s%zs
       write(102,rec=it)s%u*cos(s%alfa)-s%v*sin(s%alfa)
       write(103,rec=it)s%u*sin(s%alfa)+s%v*cos(s%alfa)
       write(104,rec=it)s%ue*cos(s%alfa)-s%ve*sin(s%alfa)
       write(105,rec=it)s%ue*sin(s%alfa)+s%ve*cos(s%alfa)
       write(106,rec=it)s%H
       write(107,rec=it)s%urms
       write(108,rec=it)s%zb
       write(109,rec=it)s%hh
       write(110,rec=it)s%Fx*cos(s%alfa)-s%Fy*sin(s%alfa)
       write(111,rec=it)s%Fx*sin(s%alfa)+s%Fy*cos(s%alfa)
       write(112,rec=it)s%E
       write(113,rec=it)s%R
       write(114,rec=it)s%D
       write(115,rec=it)sum(s%ccg,dim=3)
       write(116,rec=it)sum(s%Sug,dim=3)*cos(s%alfa)-sum(s%Svg,dim=3)*sin(s%alfa)
       write(117,rec=it)sum(s%Sug,dim=3)*sin(s%alfa)+sum(s%Svg,dim=3)*cos(s%alfa)
       write(118,rec=it)s%graindistr
       open(99,file='dims.dat',form='unformatted',access='direct',recl=wordsize*5)
       write(99,rec=1)it*1.d0,s%nx*1.d0,s%ny*1.d0,par%ngd*1.d0,par%nd*1.d0
       close(99)
    end if
end if

if(par%t>=par%tstop) then
   close(101)
   close(102)
   close(103)
   close(104)
   close(105)
   close(106)
   close(107)
   close(108)
   close(109)
   close(110)
   close(111)
   close(112)
   close(113)
   close(114)
   close(115)
   close(116)
   close(117)
end if


end subroutine

