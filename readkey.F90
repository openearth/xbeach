module readkey_module

contains
  real*8 function readkey_dbl(fname,key,defval,mnval,mxval,bcast,required)
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
    use logging_module
    implicit none
    character(len=*)  :: fname,key
    character(24)     :: printkey
    real*8            :: defval,mnval,mxval
    logical, intent(in), optional :: bcast,required

    character(len=256)   :: value
    real*8         :: value_dbl
    logical        :: lbcast,lrequired
    character(24)  :: fmt

    fmt = '(a,a,a,f0.4,a,f0.4)'

    if (present(bcast)) then
       lbcast = bcast
    else
       lbcast = .true.
    endif

    if (present(required)) then
       lrequired = required
    else
       lrequired = .false.
    endif

    !printkey=key
    printkey(2:24)=key
    printkey(1:1)=' '

    if (xmaster) then
       call readkey(fname,key,value)

       if (value/=' ') then
          read(value,'(f10.0)')value_dbl
          if (value_dbl>mxval) then
             call writelog('l','(a,a,f0.4,a,f0.4)',(printkey),' = ',value_dbl,' Warning: value > recommended value of ',mxval)
             call writelog('s','(a,a,a,f0.4)','Warning: ',trim(printkey),' > recommended value of ',mxval)
          elseif (value_dbl<mnval) then
             call writelog('l','(a,a,f0.4,a,f0.4)',(printkey),' = ',value_dbl,' Warning: value < recommended value of ',mnval)
             call writelog('s','(a,a,a,f0.4)','Warning: ',trim(printkey),' < recommended value of ',mnval)
          else
             call writelog('l','(a,a,f0.4)',(printkey),' = ',value_dbl)
          endif
       else
          if (lrequired) then
             call writelog('lse','','Error: missing required value for parameter ',printkey)
             call halt_program
          else
             value_dbl=defval
             call writelog('l','(a,a,f0.4,a)',(printkey),' = ',value_dbl,' (no record found, default value used)')
          endif
       endif
       ! write to basic params data file
       !  write(pardatfileid,*)'f ',printkey,' ',value_dbl
    endif

#ifdef USEMPI
    if (lbcast) then
       call xmpi_bcast(value_dbl)
    endif
#endif

    readkey_dbl=value_dbl
  end function readkey_dbl

  function readkey_int(fname,key,defval,mnval,mxval,bcast,required) result (value_int)
    use xmpi_module
    use logging_module
    implicit none
    character*(*)  :: fname,key
    character(24)  :: printkey
    character*256   :: value
    integer*4      :: value_int
    integer*4      :: defval,mnval,mxval
    logical, intent(in), optional :: bcast, required
    logical        :: lbcast,lrequired
    character(24)  :: fmt

    fmt = '(a,a,a,i0,a,i0)'

    if (present(bcast)) then
       lbcast = bcast
    else
       lbcast = .true.
    endif

    if (present(required)) then
       lrequired = required
    else
       lrequired = .false.
    endif

    printkey(2:24)=key
    printkey(1:1)=' '
    if (xmaster) then
       call readkey(fname,key,value)

       if (value/=' ') then
          read(value,'(i256)')value_int
          if (value_int>mxval) then
             call writelog('l',fmt,'Warning: variable ',(printkey),' ',value_int,' > recommended value of ',mxval)
             call writelog('s','(a,a,a,i0)','Warning: ',trim(printkey),' > recommended value of ',mxval)
          elseif (value_int<mnval) then
             call writelog('l',fmt,'Warning: variable ',(printkey),' ',value_int,' < recommended value of ',mnval)
             call writelog('s','(a,a,a,i0)','Warning: ',trim(printkey),' < recommended value of ',mnval)
          else
             call writelog('l','(a,a,i0)',(printkey),' = ',value_int)
          endif
       else
          if (lrequired) then
             call writelog('lse','','Error: missing required value for parameter ',printkey)
             call halt_program
          else
             value_int=defval
             call writelog('l','(a,a,i0,a)',(printkey),' = ',value_int,' (no record found, default value used)')
          endif
       endif
       ! write to basic params data file
       !  write(pardatfileid,*)'i ',printkey,' ',value_int
    endif
#ifdef USEMPI
    if (lbcast) then
       call xmpi_bcast(value_int)
    endif
#endif

  end function readkey_int

  function readkey_str(fname,key,defval,nv,nov,allowed,old,bcast,required) result (value_str)
    use xmpi_module
    use logging_module
    implicit none
    character*(*)  :: fname,key,defval
    character(24)  :: value_str
    character(24)   :: value
    integer*4      :: nv,nov,i,j
    character(24),dimension(nv) :: allowed
    character(24),dimension(nov):: old
    logical, intent(in), optional :: bcast,required
    logical        :: lbcast,lrequired,passed
    character(24)  :: printkey

    printkey(2:24)=key
    printkey(1:1)=' '

    if (present(bcast)) then
       lbcast = bcast
    else
       lbcast = .true.
    endif

    if (present(required)) then
       lrequired = required
    else
       lrequired = .false.
    endif

    passed = .false.
    if (xmaster) then
       call readkey(fname,key,value)
       ! Change to lowercase
       call lowercase(value)
       if (value == ' ') then
          if (lrequired) then
             call writelog('lse','','Error: missing required value for parameter ',printkey)
             call halt_program
          else 
             value_str=defval
             call writelog('l','(a,a,a,a)',(printkey),' = ',trim(value_str),' (no record found, default value used)')
          endif
       else
          value=adjustl(value)
          do i=1,nv
             if (trim(value)==trim(allowed(i))) then
                passed = .true.
                value_str = value
             endif
          enddo
          do j=1,nov
             if (trim(value)==trim(old(j))) then
                passed = .true.
                value_str = allowed(j)
             endif
          enddo
          if (passed) then
             call writelog('l','(a,a,a)',printkey,' = ',trim(value_str))
          else
             call writelog('sle','(a,a,a,a)','Invalid option for ',trim(printkey),' : ',trim(value))
             call halt_program
          endif
       endif
       ! write to basic params data file
       !  write(pardatfileid,*)'c ',printkey,' ',value_str
    endif
#ifdef USEMPI
    if (lbcast) then
       call xmpi_bcast(value_str)
    endif
#endif   
  end function readkey_str


  function readkey_name(fname,key,bcast,required) result (value_str)
    use xmpi_module
    use logging_module
    implicit none
    character*(*)  :: fname,key
    character(256)  :: value_str
    character*256   :: value
    logical, intent(in), optional :: bcast,required
    logical        :: lbcast,lrequired
    character(24)  :: printkey

    printkey(2:24)=key
    printkey(1:1)=' '

    if (present(bcast)) then
       lbcast = bcast
    else
       lbcast = .true.
    endif

    if (present(required)) then
       lrequired = required
    else
       lrequired = .false.
    endif

    if (xmaster) then
       call readkey(fname,key,value)
       if (value == ' ') then
          if (lrequired) then
             call writelog('lse','','Error: missing required value for parameter ',printkey)
             call halt_program
          else 
             value_str=' '
             call writelog('l',' (a,a)'    ,printkey,' = None specified')
             ! write to basic params data file
             !    write(pardatfileid,*)'c ',key,' ','none'
          endif
       else
          value_str=adjustl(value)
          call writelog('l','(a,a,a)',printkey,' = ',trim(value_str))
          ! write to basic params data file
          !    write(pardatfileid,*)'c ',printkey,' ',value_str
       endif
    endif
#ifdef USEMPI
    if (lbcast) then
       call xmpi_bcast(value_str)
    endif
#endif   
  end function readkey_name



  !
  !  readkey is only to be called from master, ie:
  !  if(xmaster) then
  !    call readkey(....)
  !  No need to cache these results. 
  !
  subroutine readkey(fname, key, value)
    use logging_module
    character(len=*), intent(in)                :: fname,key
    character(len=*), intent(out)               :: value

    integer                                     :: stat,ic
    integer                                     :: fh
    character(len=256)                            :: line
    character(len=24)                           :: foundkey
    logical                                     :: found=.true.

    open(fh, file=fname)
    found = .false.
    value = ' '
    line = 'x'
    do 
       read(fh, '(a)', iostat=stat) line
       if (stat /= 0) exit
       ic = scan(line,'=')
       if (ic>0) then
          ! remove leading and trailing spaces
          foundkey=trim(adjustl(line(1:ic-1)))
          if (trim(key).eq.trim(foundkey)) then
             ! remove leading and trailing spaces
             value=trim(adjustl(line(ic+1:256)))
             found = .true.
             exit
          end if
       endif
    enddo
    close(fh)
       
  end subroutine readkey


  ! The following code is taken from program "CHCASE" @ http://www.davidgsimpson.com/software/chcase_f90.txt:
  !  Programmer:   Dr. David G. Simpson
  !                NASA Goddard Space Flight Center
  !                Greenbelt, Maryland  20771
  !
  !  Date:         January 24, 2003
  !
  !  Language:     Fortran-90
  !
  !  Version:      1.00a
  !

  SUBROUTINE UPPERCASE(STR)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN OUT) :: STR
    INTEGER :: I, DEL


    DEL = IACHAR('a') - IACHAR('A')

    DO I = 1, LEN_TRIM(STR)
       IF (LGE(STR(I:I),'a') .AND. LLE(STR(I:I),'z')) THEN
          STR(I:I) = ACHAR(IACHAR(STR(I:I)) - DEL)
       END IF
    END DO

    RETURN

  END SUBROUTINE UPPERCASE
  !
  !  LOWERCASE
  !
  SUBROUTINE LOWERCASE(STR)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN OUT) :: STR
    INTEGER :: I, DEL


    DEL = IACHAR('a') - IACHAR('A')

    DO I = 1, LEN_TRIM(STR)
       IF (LGE(STR(I:I),'A') .AND. LLE(STR(I:I),'Z')) THEN
          STR(I:I) = ACHAR(IACHAR(STR(I:I)) + DEL)
       END IF
    END DO

    RETURN

  END SUBROUTINE LOWERCASE

! End of code taken from CHCASE




end module readkey_module
