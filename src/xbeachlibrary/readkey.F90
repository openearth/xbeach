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


module readkey_module
  use typesandkinds
  implicit none

  integer, parameter, private                   :: maxnames = 20
  character(slen), dimension(maxnames), private :: allowednames
  character(slen), dimension(maxnames), private :: oldnames
  character(slen), private                      :: varname
  integer,         dimension(maxnames), private :: intvalues
  integer, private                              :: numallowednames
  integer, private                              :: numoldnames

  interface read_v
    module procedure read_v_array
    module procedure read_v_9
  end interface read_v

contains
  real*8 function readkey_dbl(fname,key,defval,mnval,mxval,bcast,required)
    ! if USEMPI then the master process will read the parameter,
    ! this value is subsequently broadcasted to the other processes

    use xmpi_module
    use logging_module
    implicit none
    character(len=*)  :: fname,key
    character(slen)     :: printkey
    real*8            :: defval,mnval,mxval
    logical, intent(in), optional :: bcast,required

    character(slen)   :: value,tempout
    real*8         :: value_dbl
    logical        :: lbcast,lrequired
    character(slen)  :: fmt
    integer          :: ier

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
    printkey = ' '
    printkey(2:24)=trim(key)
    printkey(1:1)=' '

    if (xmaster) then
       call readkey(fname,key,value)

       if (value/=' ') then
          read(value,'(f10.0)',iostat=ier)value_dbl
          if (ier .ne. 0) then
             tempout = trim(fname)//' (value of '''//trim(printkey)//''' cannot be interpreted)'
             call report_file_read_error(tempout)
          endif
          if (value_dbl>mxval) then
             call writelog('lw','(a24,a,f0.4,a,f0.4)',(printkey),' = ',value_dbl,' Warning: value > recommended value of ',mxval)
             call writelog('s','(a24,a,a,f0.4)','Warning: ',trim(printkey),' > recommended value of ',mxval)
          elseif (value_dbl<mnval) then
             call writelog('lw','(a24,a,f0.4,a,f0.4)',(printkey),' = ',value_dbl,' Warning: value < recommended value of ',mnval)
             call writelog('s','(a24,a,a,f0.4)','Warning: ',trim(printkey),' < recommended value of ',mnval)
          else
             call writelog('l','(a24,a,f0.4)',(printkey),' = ',value_dbl)
          endif
       else
          if (lrequired) then
             call writelog('lse','','Error: missing required value for parameter ',printkey)
             call halt_program
          else
             value_dbl=defval
             call writelog('l','(a24,a,f0.4,a)',(printkey),' = ',value_dbl,' (no record found, default value used)')
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
    character(slen)  :: printkey
    character(slen)  :: value
    integer*4      :: value_int
    integer*4      :: defval,mnval,mxval,ier
    logical, intent(in), optional :: bcast, required
    logical        :: lbcast,lrequired
    character(slen)  :: fmt,tempout

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
    printkey = ' '
    printkey(2:24)=trim(key)
    printkey(1:1)=' '
    if (xmaster) then
       call readkey(fname,key,value)

       if (value/=' ') then
          read(value,'(i256)',iostat=ier)value_int
          if (ier .ne. 0) then
             tempout = trim(fname)//' (value of '''//trim(printkey)//''' cannot be interpreted)'
             call report_file_read_error(tempout)
          endif
          if (value_int>mxval) then
             call writelog('lw',fmt,'Warning: variable ',(printkey),' ',value_int,' > recommended value of ',mxval)
             call writelog('s','(a24,a,a,i0)','Warning: ',trim(printkey),' > recommended value of ',mxval)
          elseif (value_int<mnval) then
             call writelog('lw',fmt,'Warning: variable ',(printkey),' ',value_int,' < recommended value of ',mnval)
             call writelog('s','(a24,a,a,i0)','Warning: ',trim(printkey),' < recommended value of ',mnval)
          else
             call writelog('l','(a24,a,i0)',(printkey),' = ',value_int)
          endif
       else
          if (lrequired) then
             call writelog('lse','','Error: missing required value for parameter ',printkey)
             call halt_program
          else
             value_int=defval
             call writelog('l','(a24,a,i0,a)',(printkey),' = ',value_int,' (no record found, default value used)')
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
    character(slen)  :: value_str
    character(slen)   :: value
    integer*4      :: nv,nov,i,j
    character(slen),dimension(nv) :: allowed
    character(slen),dimension(nov):: old
    logical, intent(in), optional :: bcast,required
    logical        :: lbcast,lrequired,passed
    character(slen)  :: printkey

    printkey(2:slen)=key
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
             call writelog('l','(a24,a,a,a)',(printkey),' = ',trim(value_str),' (no record found, default value used)')
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
             call writelog('l','(a24,a,a)',printkey,' = ',trim(value_str))
          else
             call writelog('sle','(a24,a,a,a)','Invalid option for ',trim(printkey),' : ',trim(value))
             call writelog('sle','(a24,a,a)','Valid options for ',trim(printkey),' are:')
             do i=1,nv
                call writelog('sle','(a24)',trim(allowed(i)))
             enddo
             do j=1,nov
                call writelog('sle','(a24)',trim(old(j)))
             enddo
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
    character(slen)  :: value_str
    character(slen)   :: value
    logical, intent(in), optional :: bcast,required
    logical        :: lbcast,lrequired
    character(slen)  :: printkey

    printkey(2:slen)=key
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
             call writelog('l',' (a24,a)'    ,printkey,' = None specified')
             ! write to basic params data file
             !    write(pardatfileid,*)'c ',key,' ','none'
          endif
       else
          value_str=adjustl(value)
          call writelog('l','(a24,a,a)',printkey,' = ',trim(value_str))
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

  function readkey_dblvec(fname,key,vlength,tlength,defval,mnval,mxval,bcast,required) result (value_vec)
    use xmpi_module
    use logging_module
    implicit none
    character*(*)  :: fname,key
    integer, intent(in) :: vlength,tlength
    real*8,dimension(tlength)  :: value_vec
    real*8            :: defval,mnval,mxval
    logical, intent(in), optional :: bcast,required
    logical        :: lbcast,lrequired

    integer          :: i, ioerr
    character(slen)   :: value
    character(slen)  :: printkey,tempout

    printkey(2:slen)=key
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
       if (value/=' ') then
          read(value,*,IOSTAT=ioerr)value_vec(1:vlength)
          if (ioerr < 0) then
             call writelog('lse','','Error reading value for parameter ',printkey)
             call writelog('lse','','Check whether parameter is given sufficient number of input values')
             call halt_program
          endif
          do i=1,vlength
             if (value_vec(i)>mxval) then
                call writelog('lw','(a24,a,f0.4,a,f0.4)',(printkey),' = ',value_vec(i), &
                     ' Warning: value > recommended value of ',mxval)
                call writelog('s','(a24,a,a,f0.4)','Warning: ',trim(printkey),' > recommended value of ',mxval)
             elseif (value_vec(i)<mnval) then
                call writelog('lw','(a24,a,f0.4,a,f0.4)',(printkey),' = ',value_vec(i), &
                     ' Warning: value < recommended value of ',mnval)
                call writelog('s','(a24,a,a,f0.4)','Warning: ',trim(printkey),' < recommended value of ',mnval)
             else
                call writelog('l','(a24,a,f0.4)',(printkey),' = ',value_vec(i))
             endif
          enddo
       else
          if (lrequired) then
             call writelog('lse','','Error: missing required value for parameter ',printkey)
             call halt_program
          else
             value_vec(1:vlength)=defval
             do i=1,vlength
                call writelog('l','(a,a,f0.4,a)',(printkey),' = ',value_vec(i),' (no record found, default value used)')
             enddo
          endif
       endif
    endif

#ifdef USEMPI
    if (lbcast) then
       do i=1,vlength
          call xmpi_bcast(value_vec(i))
       enddo
    endif
#endif
  end function readkey_dblvec

  function isSetParameter(fname,key,bcast) result (isSet)
    ! Function return logical true if the keyword is specified in file,
    ! or logical false if the keyword is not specified in the file.
    use xmpi_module
    implicit none
    character*(*)   :: fname,key
    logical, intent(in), optional :: bcast
    logical         :: isSet
    character(slen)   :: value
    logical         :: lbcast

    isSet = .false.
    if (present(bcast)) then
       lbcast = bcast
    else
       lbcast = .true.
    endif

    if (xmaster) then
       call readkey(fname,key,value)
       if (value == ' ') then
          isSet = .false.
       else
          isSet = .true.
       endif
    endif
#ifdef USEMPI
    if (lbcast) then
       call xmpi_bcast(isSet)
    endif
#endif  
  end function isSetParameter

  !
  !  readkey is only to be called from master, ie:
  !  if(xmaster) then
  !    call readkey(....)
  !  No need to cache these results. 
  !
  subroutine readkey(fname,key,value)
    ! Reads through input file (fname) looking for key = value combinations
    ! Return value as string
    ! Subroutine also used to keep track of which lines have been succesfully read
    ! If called by readkey('params.txt','checkparams'), will output unsuccesful key = value
    ! combinations in params.txt
    use xmpi_module
    use logging_module
    integer                                     :: lun,i,ier,nlines,ic,ikey,itab
    character*1                                 :: ch
    character(len=*), intent(in)                :: fname,key
    character(len=*), intent(out)               :: value
    character(slen), dimension(1024),save          :: keyword,values
    character(slen)                                :: line,lineWithoutSpecials
    integer, save                               :: nkeys
    character(slen), save                          :: fnameold=''
    integer, dimension(:),allocatable,save      :: readindex

    ! If the file name of the input file changes, the file should be reread
    if (fname/=fnameold) then
       ! Make sure this reset only recurs when the input file name changes
       fnameold=fname
       nkeys=0
       ier=0
       ! Read the file for all lines with "=" 
       call writelog('ls','','XBeach reading from ',trim(fname))
       lun=99
       i=0
       open(lun,file=fname)
       do while (ier==0)
          read(lun,'(a)',iostat=ier)ch
          if (ier==0)i=i+1
       enddo
       close(lun)
       nlines=i
       ! reset keyword values and readindex 
       keyword = ''
       values = '' 
       if (allocated(readindex)) deallocate(readindex)
       ! Read through the file to fill all the keyword = value combinations
       open(lun,file=fname)
       ikey=0
       do i=1,nlines
          read(lun,'(a)')line
          do itab=1,slen
              if (ichar(line(itab:itab))<32 .or. ichar(line(itab:itab))>126) then  ! this is anything not in standard
                                                                                   ! alphanumeric
                  lineWithoutSpecials(itab:itab) = ' '
              else
                  lineWithoutSpecials(itab:itab) = line(itab:itab)
              endif
          enddo
          line = lineWithoutSpecials
          ic=scan(line,'=')
          if (ic>0) then
             ikey=ikey+1
             keyword(ikey)=adjustl(line(1:ic-1))
             values(ikey)=adjustl(line(ic+1:slen))
          endif
       enddo
       nkeys=ikey
       close(lun)
       ! allocate index vector that stores which values have succesfully been called to be read
       allocate(readindex(nkeys))
       readindex=0
    endif

    ! Compare the input key with any keyword stored in the keyword vector and return the value.
    ! A succesful key - keyword match is recorded in readindex with a value "1"
    ! Note: in case more than one keyword matches the key, the first keyword - value combination is returned
    value=' '
    do ikey=1,nkeys
       if (key.eq.keyword(ikey)) then
          value=values(ikey)
          readindex(ikey)=1
          exit
       endif
    enddo

    ! Easter egg!
    ! With call for key "checkparams", the subroutine searches readindex for keyword - value combinations that
    ! have not yet been read. It returns a warning to screen and log file for each unsuccesful keyword.
    if (key .eq. 'checkparams') then
       do ikey=1,nkeys
          if (readindex(ikey)==0) then
             call writelog('slw','','Unknown, unused or multiple statements of parameter ', &
                                    trim(uppercase(keyword(ikey))),' in ',trim(fname))
          endif
       enddo
    endif

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
  !                1.1 : Modified uppercase into function form by R.T. McCall 23/7/2013
  !

  pure function UPPERCASE(STR) result(upperstr)

    IMPLICIT NONE

    CHARACTER(LEN=*),intent(in) :: STR
    character(slen)             :: upperstr
    INTEGER                     :: I, DEL

    upperstr = STR

    DEL = IACHAR('a') - IACHAR('A')
    
    DO I = 1, LEN_TRIM(upperstr)
       IF (LGE(upperstr(I:I),'a') .AND. LLE(upperstr(I:I),'z')) THEN
          upperstr(I:I) = ACHAR(IACHAR(upperstr(I:I)) - DEL)
       END IF
    END DO

  end function UPPERCASE
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



! The following
! defines the following method to read a string parameter:
!
!  for example:

!  integer NAME1, NAME2, NAME3
!  call setallowednames('name1',NAME1,'name2',NAME2,'name3',NAME3)
!   this defines the allowed names ('name1', 'name2', 'name3') and
!   the to-be associated integer values (NAME1, NAME2, NAME3)
!  call setoldname('0','1')
!   this defines alternate allowed names for 'name1' and 'name2'
!  integer intvalue
!  character(slen) stringvalue
!  call parmapply('gridform',2,intvalue[,stringvalue])
!                            |
!   this searches for 'gridform' in file 'params.txt'. 
!   let us assume that there is a line
!   gridform = name3
!   then intvalue becomes NAME2, and stringvalue becomes 'name3'
!   If no line defining gridform is found, the default numer allowed name
!   is used, in this case number 2: NAME2 and 'name2'
!
!   Notes:
!     setoldnames must be called after setallowednames and is optional
!     the 3rd parameter in parmapply is optional
!


    subroutine parmapply(vname,idefname,parm,parm_str,bcast,required)
      use typesandkinds
      use xmpi_module
      implicit none
      character(*), intent(in)            :: vname
      integer,      intent(in)            :: idefname
      integer,      intent(out)           :: parm
      character(*), intent(out), optional :: parm_str
      logical,      intent(in), optional  :: bcast,required

      character(slen)                     :: d
      integer                             :: i
      logical                             :: lbcast

      d = readkey_str('params.txt',vname,allowednames(idefname), &
                      numallowednames,numoldnames,allowednames,oldnames, &
                      bcast, required)

      if (present(bcast)) then
        lbcast = bcast
      else
        lbcast = .true.
      endif

      if (xmaster) then
        do i=1,numallowednames
          if (d .eq. allowednames(i)) then
            parm = intvalues(i)
            if (present(parm_str)) then
              parm_str = d
            endif
            exit
          endif
        enddo
      endif

#ifdef USEMPI
      if (lbcast) then
        call xmpi_bcast(parm)
      endif
#endif

    end subroutine parmapply
     
    subroutine setallowednames(a1,v1,a2,v2,a3,v3,a4,v4,a5,v5,a6,v6,a7,v7,a8,v8, &
                               a9,v9,a10,v10,a11,v11,a12,v12,a13,v13,a14,v14,   &
                               a15,v15,a16,v16,a17,v17,a18,v18,a19,v19,a20,v20)
    character(*), intent(in) :: a1
    character(*), intent(in), optional :: a2,a3,a4,a5,a6,a7,a8,a9,a10 &
                    ,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20
    integer   ,   intent(in) :: v1
    integer   ,   intent(in), optional :: v2,v3,v4,v5,v6,v7,v8,v9,v10 &
                    ,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20
    numoldnames = 0
    allowednames(1) = a1
    intvalues(1)    = v1
    numallowednames = 1
    if (present(a2)) then
      allowednames(2) = a2
      intvalues(2)    = v2
      numallowednames = 2
    endif
    if (present(a3)) then
      allowednames(3) = a3
      intvalues(3)    = v3
      numallowednames = 3
    endif
    if (present(a4)) then
      allowednames(4) = a4
      intvalues(4)    = v4
      numallowednames = 4
    endif
    if (present(a5)) then
      allowednames(5) = a5
      intvalues(5)    = v5
      numallowednames = 5
    endif
    if (present(a6)) then
      allowednames(6) = a6
      intvalues(6)    = v6
      numallowednames = 6
    endif
    if (present(a7)) then
      allowednames(7) = a7
      intvalues(7)    = v7
      numallowednames = 7
    endif
    if (present(a8)) then
      allowednames(8) = a8
      intvalues(8)    = v8
      numallowednames = 8
    endif
    if (present(a9)) then
      allowednames(9) = a9
      intvalues(9)    = v9
      numallowednames = 9
    endif
    if (present(a10)) then
      allowednames(10) = a10
      intvalues(10)    = v10
      numallowednames = 10
    endif
    if (present(a11)) then
      allowednames(11) = a11
      intvalues(11)    = v11
      numallowednames = 11
    endif
    if (present(a12)) then
      allowednames(12) = a12
      intvalues(12)    = v12
      numallowednames = 12
    endif
    if (present(a13)) then
      allowednames(13) = a13
      intvalues(13)    = v13
      numallowednames = 13
    endif
    if (present(a14)) then
      allowednames(14) = a14
      intvalues(14)    = v14
      numallowednames = 14
    endif
    if (present(a15)) then
      allowednames(15) = a15
      intvalues(15)    = v15
      numallowednames = 15
    endif
    if (present(a16)) then
      allowednames(16) = a16
      intvalues(16)    = v16
      numallowednames = 16
    endif
    if (present(a17)) then
      allowednames(17) = a17
      intvalues(17)    = v17
      numallowednames = 17
    endif
    if (present(a18)) then
      allowednames(18) = a18
      intvalues(18)    = v18
      numallowednames = 18
    endif
    if (present(a19)) then
      allowednames(19) = a19
      intvalues(19)    = v19
      numallowednames = 19
    endif
    if (present(a20)) then
      allowednames(20) = a20
      intvalues(20)    = v20
      numallowednames = 20
    endif

    end subroutine setallowednames

    subroutine setoldnames(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10 &
                       ,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)
    character(*), intent(in) :: a1
    character(*), intent(in), optional :: a2,a3,a4,a5,a6,a7,a8,a9,a10 &
                    ,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20
    oldnames(1) = a1
    numoldnames = 1
    if (present(a2)) then
      oldnames(2) = a2
      numoldnames = 2
    endif
    if (present(a3)) then
      oldnames(3) = a3
      numoldnames = 3
    endif
    if (present(a4)) then
      oldnames(4) = a4
      numoldnames = 4
    endif
    if (present(a5)) then
      oldnames(5) = a5
      numoldnames = 5
    endif
    if (present(a6)) then
      oldnames(6) = a6
      numoldnames = 6
    endif
    if (present(a7)) then
      oldnames(7) = a7
      numoldnames = 7
    endif
    if (present(a8)) then
      oldnames(8) = a8
      numoldnames = 8
    endif
    if (present(a9)) then
      oldnames(9) = a9
      numoldnames = 9
    endif
    if (present(a10)) then
      oldnames(10) = a10
      numoldnames = 10
    endif
    if (present(a11)) then
      oldnames(11) = a11
      numoldnames = 11
    endif
    if (present(a12)) then
      oldnames(12) = a12
      numoldnames = 12
    endif
    if (present(a13)) then
      oldnames(13) = a13
      numoldnames = 13
    endif
    if (present(a14)) then
      oldnames(14) = a14
      numoldnames = 14
    endif
    if (present(a15)) then
      oldnames(15) = a15
      numoldnames = 15
    endif
    if (present(a16)) then
      oldnames(16) = a16
      numoldnames = 16
    endif
    if (present(a17)) then
      oldnames(17) = a17
      numoldnames = 17
    endif
    if (present(a18)) then
      oldnames(18) = a18
      numoldnames = 18
    endif
    if (present(a19)) then
      oldnames(19) = a19
      numoldnames = 19
    endif
    if (present(a20)) then
      oldnames(20) = a20
      numoldnames = 20
    endif

    end subroutine setoldnames

    subroutine read_v_array(fid,a)
    use xmpi_module
    !
    ! reads array from unit fid
    ! to be called by everyone
    ! only xmaster reads
    ! fid is only needed on xmaster
    !
    integer, intent(in)               :: fid
    real*8, dimension(:), intent(out) :: a
    if (xmaster) then
      read(fid,*) a
    endif
#ifdef USEMPI
    call xmpi_bcast(a)
#endif
    end subroutine read_v_array

    subroutine read_v_9(fid,a,a1,a2,a3,a4,a5,a6,a7,a8)
    use xmpi_module
    !
    ! reads doubles from unit fid
    ! to be called by everyone
    ! only xmaster reads
    ! fid is only needed on xmaster
    !
    integer, intent(in) :: fid
    real*8, intent(out) :: a
    real*8, optional, intent(out) :: a1,a2,a3,a4,a5,a6,a7,a8
    if (present(a8)) then
      if(xmaster) read(fid,*) a,a1,a2,a3,a4,a5,a6,a7,a8
#ifdef USEMPI
      call xmpi_bcast(a1)
      call xmpi_bcast(a2)
      call xmpi_bcast(a3)
      call xmpi_bcast(a4)
      call xmpi_bcast(a5)
      call xmpi_bcast(a6)
      call xmpi_bcast(a7)
      call xmpi_bcast(a8)
#endif
    elseif (present(a7)) then
      if(xmaster) read(fid,*) a,a1,a2,a3,a4,a5,a6,a7
#ifdef USEMPI
      call xmpi_bcast(a1)
      call xmpi_bcast(a2)
      call xmpi_bcast(a3)
      call xmpi_bcast(a4)
      call xmpi_bcast(a5)
      call xmpi_bcast(a6)
      call xmpi_bcast(a7)
#endif
    elseif (present(a6)) then
      if(xmaster) read(fid,*) a,a1,a2,a3,a4,a5,a6
#ifdef USEMPI
      call xmpi_bcast(a1)
      call xmpi_bcast(a2)
      call xmpi_bcast(a3)
      call xmpi_bcast(a4)
      call xmpi_bcast(a5)
      call xmpi_bcast(a6)
#endif
    elseif (present(a5)) then
      if(xmaster) read(fid,*) a,a1,a2,a3,a4,a5
#ifdef USEMPI
      call xmpi_bcast(a1)
      call xmpi_bcast(a2)
      call xmpi_bcast(a3)
      call xmpi_bcast(a4)
      call xmpi_bcast(a5)
#endif
    elseif (present(a4)) then
      if(xmaster) read(fid,*) a,a1,a2,a3,a4
#ifdef USEMPI
      call xmpi_bcast(a1)
      call xmpi_bcast(a2)
      call xmpi_bcast(a3)
      call xmpi_bcast(a4)
#endif
    elseif (present(a3)) then
      if(xmaster) read(fid,*) a,a1,a2,a3
#ifdef USEMPI
      call xmpi_bcast(a1)
      call xmpi_bcast(a2)
      call xmpi_bcast(a3)
#endif
    elseif (present(a2)) then
      if(xmaster) read(fid,*) a,a1,a2
#ifdef USEMPI
      call xmpi_bcast(a1)
      call xmpi_bcast(a2)
#endif
    elseif (present(a1)) then
      if(xmaster) read(fid,*) a,a1
#ifdef USEMPI
      call xmpi_bcast(a1)
#endif
    else 
      if(xmaster) read(fid,*) a
    endif
#ifdef USEMPI
    call xmpi_bcast(a)
#endif
    end subroutine read_v_9

    integer function count_lines(f)
    !
    ! returns number of lines in file f
    ! to be called by all, xmaster will count
    ! result is broadcasted to all
    !
    use xmpi_module
    use filefunctions
    character(len=*), intent(in) :: f
    integer lines,fid,ierr
    fid = create_new_fid()
    if(xmaster) then
      open(fid,file=f)
      lines = 0
      ierr = 0
      do while(ierr == 0)
        read(fid,*,iostat=ierr)
        if (ierr == 0) lines = lines+1
      enddo
      close(fid)
    endif
#ifdef USEMPI
    call xmpi_bcast(lines)
#endif
    count_lines = lines

    end function count_lines
end module readkey_module
