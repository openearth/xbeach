module readkey_module

contains
real*8 function readkey_dbl(fname,key,defval,mnval,mxval)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Copyright (C) 2007 UNESCO-IHE, WL|Delft Hydraulics and Delft University !
! Dano Roelvink, Ap van Dongeren, Ad Reniers, Jamie Lescinski,            !
! Jaap van Thiel de Vries, Robert McCall                                  !       
!                                                                         !
! d.roelvink@unesco-ihe.org                                               !
! UNESCO-IHE Institute for Water Education                                !
! P.O. Box 3015                                                           !
! 2601 DA Delft                                                           !
! The Netherlands                                                         !
!                                                                         !
! This library is free software; you can redistribute it and/or           !
! modify it under the terms of the GNU Lesser General Public              !
! License as published by the Free Software Foundation; either            !
! version 2.1 of the License, or (at your option) any later version.      !
!                                                                         !
! This library is distributed in the hope that it will be useful,         !
! but WITHOUT ANY WARRANTY; without even the implied warranty of          !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        !
! Lesser General Public License for more details.                         !
!                                                                         !
! You should have received a copy of the GNU Lesser General Public        !
! License along with this library; if not, write to the Free Software     !
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307     !
! USA                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if USEMPI then the master process will read the parameter,
! this value is subsequently broadcasted to the other processes

use xmpi_module
implicit none
character(len=*)  :: fname,key
real*8            :: defval,mnval,mxval

character*80   :: value
real*8         :: value_dbl
if (xmaster) then
  call readkey(fname,key,value)

  if (value/=' ') then
     read(value,'(f10.0)')value_dbl
     if (value_dbl>mxval) then
        write(*,*)'Warning: variable ',trim(key),value_dbl,' > recommended value of ',mxval 
     elseif (value_dbl<mnval) then
        write(*,*)'Warning: variable ',trim(key),value_dbl,' < recommended value of ',mnval 
     else
        write(*,*)trim(key),' = ',value_dbl
     endif
  else
     value_dbl=defval
     if(xmaster) then
       write(*,*)trim(key),' = ',value_dbl,' (no record found, default value used)'
     endif
  endif
endif

#ifdef USEMPI
  if (xmpi_bckey) call xmpi_bcast(value_dbl)
#endif

readkey_dbl=value_dbl
end function readkey_dbl

function readkey_int(fname,key,defval,mnval,mxval) result (value_int)
use xmpi_module
implicit none
character*(*)  :: fname,key
character*80   :: value
integer*4      :: value_int
integer*4      :: defval,mnval,mxval

if (xmaster) then
  call readkey(fname,key,value)

  if (value/=' ') then
     read(value,'(i80)')value_int
     if (value_int>mxval) then
        write(*,*)'Warning: variable ',trim(key),value_int,' > recommended value of ',mxval 
     elseif (value_int<mnval) then
        write(*,*)'Warning: variable ',trim(key),value_int,' < recommended value of ',mnval 
     else
        write(*,*)trim(key),' = ',value_int
     endif
  else
     value_int=defval
     if(xmaster) then
       write(*,*)trim(key),' = ',value_int,' (no record found, default value used)'
     endif
  endif
endif
#ifdef USEMPI
  if (xmpi_bckey) call xmpi_bcast(value_int)
#endif

end function readkey_int

!
!  readkey is only to be called from master, ie:
!  if(xmaster) then
!    call readkey(....)
!
subroutine readkey(fname,key,value)
integer                                     :: lun,i,ier,nlines,ic,ikey
character*1                                 :: ch
character(len=*), intent(in)                :: fname,key
character(len=*), intent(out)               :: value
character*80, dimension(:),allocatable,save :: keyword,values
character*80                                :: line
logical, save                               :: first=.true.
integer, save                               :: nkeys
character*80, save                          :: fnameold='first_time.exe'
integer, dimension(:),allocatable,save          :: readindex

if (fname/=fnameold) then                   ! Open new file if fname changes
    if (fnameold/='first_time.exe') then    ! only if not the first time older versions
        deallocate(keyword)
        deallocate(values)
                deallocate(readindex)
    end if
    first=.true.
    fnameold=fname
    nkeys=0
    ier=0
end if


if (first) then
   write(*,*)'readkey: Reading from ',trim(fname),' ...........'
   first=.false.
   lun=99
   i=0
   open(lun,file=fname)
   do while (ier==0)
      read(lun,'(a)',iostat=ier)ch
      if (ier==0)i=i+1
   enddo
   close(lun)
   nlines=i

   allocate(keyword(nlines))
   allocate(values(nlines))

   open(lun,file=fname)
   ikey=0
   do i=1,nlines
      read(lun,'(a)')line
      ic=scan(line,'=')
      if (ic>0) then
         ikey=ikey+1
         keyword(ikey)=adjustl(line(1:ic-1))
         values(ikey)=adjustl(line(ic+1:80))
      endif
   enddo
   nkeys=ikey
   close(lun)
   allocate(readindex(nkeys))
   readindex=0
endif

value=' '
do ikey=1,nkeys
   if (key.eq.keyword(ikey)) then
      value=values(ikey)
          readindex(ikey)=1
   endif
enddo


! If required, do a check whether params are not used or unknown
if (key .eq. 'checkparams') then
        do ikey=1,nkeys
                if (readindex(ikey)==0) then
                        write(*,*) 'Unknown, unused or multiple statements of parameter ',trim(keyword(ikey)),' in ',trim(fname)
                endif
        enddo
endif


end subroutine readkey

end module readkey_module
