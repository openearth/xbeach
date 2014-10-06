!  test version  wwvv, do not use in production version of xbeach
!
! this program makes includefiles from spaceparams.tmpl
! indextos.gen            is for spaceparams.F90
! mnemonic.gen            is for mnemonic.F90
! chartoindex.gen         is for mnemonic.F90
! space_alloc_arrays.gen  is for spaceparams.F90
! spacedecl.gen           is for spaceparams.F90
! space_ind.gen           is for s.ind
! space_inp.gen           is for s.inp
! spaceparams.h           definition of macro's: #define nx s%nx
! spaceparamsu.h          undefinition of macro's: #undef nx
! space.sed               sed commands
! parameters.inc          is for xbeach.F90 after params.F90
! space_define.h          defines macros like #define x s%X
! space_undef.h           undefines what space_define.h defines
!
! and also parameters.inc and getkey.gen, both from params.F90
!
! see also README.includefiles
!          spaceparams.tmpl
!
module makemodule
! To avoid referencing typesandkinds...
  integer, parameter         :: slen = 1024
  integer, parameter         :: dimnamelen = 20
!
  character(slen), parameter :: templatefilename = 'spaceparams.tmpl'
  character(slen)            :: outputfilename
  integer             :: infile = 1
  integer             :: outfile = 2
  character(slen), parameter :: spacedeclname           = 'spacedecl.gen'
  character(slen), parameter :: mnemonicname            = 'mnemonic.gen'
  character(slen), parameter :: indextosname            = 'indextos.gen'
  character(slen), parameter :: space_alloc_arraysname  = 'space_alloc_arrays.gen'
  character(slen), parameter :: space_indname           = 'space_ind.gen'
  character(slen), parameter :: space_inpname           = 'space_inp.gen'
  character(slen), parameter :: spaceparams_hname       = 'spaceparams.h'
  character(slen), parameter :: uspaceparams_hname      = 'spaceparamsu.h'
  character(slen), parameter :: space_sedname           = 'space.sed'
  character(slen), parameter :: chartoindexname         = 'chartoindex.gen'
  character(slen), parameter :: paramsincname           = 'parameters.inc'
  character(slen), parameter :: getkeygenname           = 'getkey.gen'
  character(slen), parameter :: space_definename        = 'space_define.h'
  character(slen), parameter :: space_undefname         = 'space_undef.h'

  integer                        :: numvars, maxnamelen, inumvars0, rnumvars0
  character(slen)            :: type, name, comment, line, broadcast, description, units, standardname
  integer                        :: rank
  logical                        :: varfound
  logical                        :: endfound
  logical                        :: dimensioned
  character(dimnamelen), dimension(10)   :: dims       ! dimensions are maximum 20 length statement, maximum 10 different dimensions
  integer                        :: maxdimlen, maxrank
  character(6),parameter         :: sp = '      '

contains

  Pure Function caseswitch (str) Result (string)

!   ==============================
!   switches case of string: AbC -> aBc
!   ==============================

    Implicit None
    Character(*), Intent(In) :: str
    Character(LEN(str))      :: string

    Integer :: ic, i

    Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, LEN_TRIM(str)
        ic = INDEX(low, str(i:i))
        if (ic > 0) then
          string(i:i) = cap(ic:ic)
        endif
        ic = INDEX(cap, str(i:i))
        if (ic > 0) then
          string(i:i) = low(ic:ic)
        endif
    end do

  End Function caseswitch

  !
  !  reads next line
  !
  subroutine getline
    implicit none
    ! character(slen) :: s
    integer         :: i
    do
       read(infile,'(a)',end=100) line
       ! some compilers don't think carriage return is white space.
       ! we convert all cr, nl, tab into space
       do i=1,len(line)
          select case(ichar(line(i:i)))
          case (9,10,13)
             line(i:i) = ' '
          end select
       enddo

       line = adjustl(line)
       if (line(1:2) .ne. '!*') then
          exit
       endif
    enddo
    endfound = .false.
    return
100 continue
    endfound = .true.
  end subroutine getline
  !
  !  opens inputfile
  !
  subroutine openinfile
    implicit none
    open(infile, file=templatefilename)
  end subroutine openinfile
  !
  !  opens outputfile
  !
  subroutine openoutput
    implicit none
    open (outfile,file=trim(adjustl(outputfilename)),recl=slen)
  end subroutine openoutput
  !
  !  extract items from an input line
  !  type, rank, varfound, comment, name
  !
  subroutine getitems
    implicit none

    character(slen)  :: w
    integer              :: j,l,i

    l=1
    type = ''
    name = ''
    units = ''
    standardname = ''
    description = ''
    comment = ''

    call getnext(line,l,type)
    if (len(trim(type)) .eq. 0) then
       varfound = .false.
       return
    endif

    if (type(1:1) .eq. '!' ) then
       comment = line
       type = ''
       varfound = .false.
       return
    endif

    varfound = .true.
    call getnext(line,l,w)
    read (w,'(i20)') rank

    call getnext(line,l,name)

    do i=1,rank
       call getnext(line,l,dims(i))
       maxdimlen = max(maxdimlen,(len_trim(dims(i))))
    enddo

    call getnext(line,l,broadcast)

    call getnext(line,l,units)

    call getnext(line,l,standardname)

    ! Get the first word of the description
    call getnext(line,l,description)

    ! search for a comment
    j = scan(line,'!')
    if (j .ne. 0) then
       comment = line(j:)
       ! store the line up to comment
       description = trim(description) // trim(line(l:j-1))
    else
       ! store the rest of the line in the description
       description = trim(description) // trim(line(l:))
    endif


  end subroutine getitems
  !
  !  reads next word from characterstring
  !
  subroutine getnext(a,l,w)
    implicit none
    character(len=*), intent(in)  :: a    ! string to read from
    character(len=*), intent(out) :: w    ! word found
    integer, intent(inout)        :: l    ! index, where to start in a
    ! normally, at first call of this routine, 
    ! l should be one. l is updated.

    integer                       :: j,i,k
    integer                       :: la 
    character                     :: tab 

    tab = char(9)

    la = len(a)
    j = la
    do i = l, la
       if (a(i:i) .ne. ' ' .and. a(i:i) .ne. tab) then
          j = i
          exit
       endif
    enddo
    l = j
    w = ''
    j = la
    k = 1
    do i = l, la
       if (a(i:i) .eq. ' ' .or. a(i:i) .eq. tab) then 
          j = i
          exit
       endif
       w(k:k) = a(i:i)
       k = k + 1
    enddo
    l = j
  end subroutine getnext
  !
  !  emits a waring on the outputfile
  !
  subroutine warning
    implicit none
    write(outfile,'(a)') '!  DO NOT EDIT THIS FILE'
    write(outfile,'(a)') '!  But edit '//trim(templatefilename)//' and/or makeincludes.F90'
    write(outfile,'(a)') '!  compile makeincludes.F90 and run the compilate'
    write(outfile,'(a)') '!  Compiling and running is taken care of by the Makefile'

  end subroutine warning
  !
  !  emits directions for vim editor
  !
  subroutine format_fortran
    implicit none
    write(outfile,'(a)') '!directions for vi vim: filetype=fortran : syntax=fortran'
  end subroutine format_fortran
  !
  !  just counting
  !
  subroutine counting
    implicit none
    call openinfile

    !
    !      just counting
    !
    numvars = 0
    maxnamelen = 0
    maxrank = 0
    maxdimlen = 0
    inumvars0 = 0
    rnumvars0 = 0
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          numvars = numvars+1
          maxnamelen = max(maxnamelen,len_trim(name))
          maxrank = max(maxrank, rank)
          if (rank .eq. 0) then
             if (type(1:1) .eq. 'i' .or. type(1:1) .eq. 'I') then
                inumvars0=inumvars0+1
             else
                rnumvars0=rnumvars0+1
             endif
          endif
       endif
    enddo
    close(infile)
  end subroutine counting
  !
  !  construction of spacedecl.gen
  !
  subroutine makespacedecl
    implicit none
    character(slen) rankstr
    call openinfile
    call openoutput
    call warning
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          select case(rank)
          case(0)
             rankstr = '                                '
          case(1)
             rankstr = ',dimension(:),       allocatable'
          case(2)
             rankstr = ',dimension(:,:),     allocatable'
          case(3)
             rankstr = ',dimension(:,:,:),   allocatable'
          case(4)
             rankstr = ',dimension(:,:,:,:), allocatable'
          case default
             rankstr = 'invalid number of dimensions'
          end select
          write(outfile,'(a)') sp//trim(type)//trim(rankstr)//' :: '//trim(name)// &
               '  '//trim(comment)
       else
          write(outfile,'(a)') trim(comment)
       endif
    enddo
    call format_fortran
    close (outfile)
    close(infile)
  end subroutine makespacedecl
  !
  !  construction of mnemonic.gen
  !
  subroutine makemnemonic
    implicit none
    integer :: j
    call openinfile
    call openoutput
    call warning
    write (outfile,'(a,i4)') sp//'  integer, parameter :: numvars    = ', numvars
    write (outfile,'(a,i4)') sp//'  integer, parameter :: maxnamelen = ', maxnamelen
    write (outfile,'(a,i4)') sp//'  integer, parameter :: maxrank = ', maxrank

    do
       call getline
       if(endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          write(outfile,'(a)') sp//'  character(len=maxnamelen),parameter ::  mnem_'// &
               name(1:maxnamelen)//" = '"//name(1:maxnamelen)//"'"
       endif
    enddo

    rewind infile
    write (outfile,'(a)') sp//&
         '  character(len=maxnamelen),dimension(numvars),parameter :: mnemonics= (/ &'

    j=1
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          if ( j .ne. numvars) then
             write(outfile,'(a,i4)') sp//'   mnem_'//name(1:maxnamelen)//', & ! ',j
          else
             write(outfile,'(a,i4)') sp//'   mnem_'//name(1:maxnamelen)//'  & ! ',j
          endif
          j = j + 1
       endif
    enddo
    write(outfile,'(a)') sp//"  /)"

    call format_fortran
    close(outfile)
    close(infile)
  end subroutine makemnemonic
  !
  ! construction of chartoindex.gen
  !
  subroutine makechartoindex
    implicit none
    integer :: j
    call openinfile
    call openoutput
    call warning

    write(outfile,'(a)') sp//'  select case(line)'
    j = 0
    do
       call getline
       if(endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          j = j+1
          write(outfile,'(a)') sp//'    case(mnem_'//name(1:maxnamelen)//')'
          write(outfile,'(a,i5)') sp//'             chartoindex =',j
       endif
    enddo

    write(outfile,'(a)') sp//'  end select'

    call format_fortran
    close(outfile)
    close(infile)
  end subroutine makechartoindex
  !
  !  construction of indextos.gen
  !
  subroutine makeindextos
    implicit none
    integer             :: i,j
    character           :: ctype
    character(slen) :: rankstr,rankstr1
    character(slen) :: extrastr
    call openinfile
    call openoutput
    call warning

    j=1
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          if (type(1:1) .eq. 'i' .or. type(1:1) .eq. 'I') then
             ctype = 'i'
          else
             ctype = 'r'
          endif
          select case(ctype)
          case('r')
             select case(rank)
             case(0)
                rankstr = 't%r0'
                rankstr1 = 'r0'
             case(1)
                rankstr = 't%r1'
                rankstr1 = 'r1'
             case(2)
                rankstr = 't%r2'
                rankstr1 = 'r2'
             case(3)
                rankstr = 't%r3'
                rankstr1 = 'r3'
             case(4)
                rankstr = 't%r4'
                rankstr1 = 'r4'
             case default
                rankstr = 'wrong dimensionality'
             end select
          case ('i')
             select case(rank)
             case(0)
                rankstr = 't%i0'
                rankstr1 = 'i0'
             case(1)
                rankstr = 't%i1'
                rankstr1 = 'i1'
             case(2)
                rankstr = 't%i2'
                rankstr1 = 'i2'
             case(3)
                rankstr = 't%i3'
                rankstr1 = 'i3'
             case(4)
                rankstr = 't%i4'
                rankstr1 = 'i4'
             case default
                rankstr = 'wrong dimensionality'
             end select
          end select
          write(outfile,'(a,i3,a)') sp//'  case(',j,')'
          write(outfile,'(a)')      sp//'    point'//trim(rankstr1)//'  => s%'//trim(name)
          write(outfile,'(a,i3)')   sp//'    t%rank = ',rank
          write(outfile,'(a)')      sp//"    t%type = '"//ctype//"'"
          write(outfile,'(a)')      sp//"    t%name = '"//trim(name(1:maxnamelen))//"'"
          write(outfile,'(a)')      sp//"    t%btype = '"//broadcast(1:1)//"'"
          ! remove the [ ] around the units, this should be easier....
          units = trim(units) // '' ! make sure we trim spaces first
          if (units(1:1) .eq. '[') then
             units=units(2:) // '' ! avoid length issues
          endif
          if (index(units,']') > 0) then
             units=units(1:index(units,']')-1) // '' ! make sure we add an extra ''
          endif
          write(outfile,'(a)')      sp//"    t%units= '"//trim(units)//"'"
          ! write standard names
          if (standardname(1:1) .eq. '[') then
             standardname=standardname(2:) // '' ! avoid length issues
          endif
          if (index(standardname,']') > 0) then
             standardname=standardname(1:index(standardname,']')-1) // '' ! make sure we add an extra ''
          endif
          write(outfile,'(a)')      sp//"    t%standardname= '"//trim(standardname)//"'"
          write(outfile,'(a)')      sp//"    t%description= '"//trim(description)//"'"
          if (rank>0) then
             ! Stre the dimensions
             write(outfile,'(a,i3,a)', advance='no')      sp//"    t%dimensions(1:", rank,") = (/ " 
             do i=1,rank
                if (i < rank) then
                   extrastr = ", "
                else
                   extrastr = "/)"
                endif
                write(outfile,'(a)',advance='no')   "'"
                write(outfile,'(a)', advance='no'), dims(i) ! Don't trim this one, we need the spaces...
                write(outfile,'(a)',advance='no'),  "'" // trim(extrastr)
             enddo
             write(outfile, *) ! empty line
          endif
          j = j + 1
       endif
    enddo
    call format_fortran
    close (outfile)
    close( infile)
  end subroutine makeindextos

  !
  !  make allocation statements 
  !
  subroutine makespace_alloc_arrays
    implicit none
    integer              :: j

    call openinfile
    call openoutput
    call warning

    j=1
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          select case(rank)
          case (0)
          case (1)
             write(outfile,'(a)') sp//'  allocate(s%'//name(1:maxnamelen)// &
                  '('//dims(1)(1:maxdimlen)//repeat(' ',3*maxdimlen+3)//'))'
          case (2)
             write(outfile,'(a)') sp//'  allocate(s%'//name(1:maxnamelen)// &
                  '('//dims(1)(1:maxdimlen)//','//dims(2)(1:maxdimlen)//repeat(' ',2*maxdimlen+2)//'))'
          case (3)
             write(outfile,'(a)') sp//'  allocate(s%'//name(1:maxnamelen)// &
                  '('//dims(1)(1:maxdimlen)//','//dims(2)(1:maxdimlen)//','//dims(3)(1:maxdimlen)// &
                  repeat(' ',maxdimlen+1)//'))'
          case (4)
             write(outfile,'(a)') sp//'  allocate(s%'//name(1:maxnamelen)// &
                  '('//dims(1)(1:maxdimlen)//','//dims(2)(1:maxdimlen)//','//dims(3)(1:maxdimlen)// &
                  ','//dims(4)(1:maxdimlen)//'))'
          case default
             write(outfile,'(a)') sp//'wrong line in '//trim(templatefilename)
             write(outfile,'(a)') sp//trim(line)
          end select
       endif
    enddo

    call format_fortran
    close(outfile)
    close(infile)

  end subroutine makespace_alloc_arrays

  subroutine makespaceind
    implicit none
    character(slen) rankstr
    call openinfile
    call openoutput
    call warning
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          select case(rank)
          case(0)
             rankstr = ',pointer                   '
          case(1)
             rankstr = ',dimension(:),      pointer'
          case(2)
             rankstr = ',dimension(:,:),    pointer'
          case(3)
             rankstr = ',dimension(:,:,:),  pointer'
          case(4)
             rankstr = ',dimension(:,:,:,:),pointer'
          case default
             rankstr = 'invalid number of dimensions'
          end select
          write(outfile,'(a)') sp//trim(type)//trim(rankstr)//' :: '//trim(name)
       endif
    enddo
    call format_fortran
    close (outfile)
    close(infile)
  end subroutine makespaceind

  subroutine makespaceinp
    implicit none
    call openinfile
    call openoutput
    call warning
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          write(outfile,'(a)') sp//trim(name)//' => '//'s%'//trim(name)
       endif
    enddo
    call format_fortran
    close (outfile)
    close(infile)
  end subroutine makespaceinp

  subroutine makespaceparamsh
    implicit none
    call openinfile
    call openoutput
    call warning
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          write(outfile,'(a)') "#define "//trim(name)//' '//'s%'//trim(name)
       endif
    enddo
    close (outfile)
    close(infile)
  end subroutine makespaceparamsh

  subroutine umakespaceparamsh
    implicit none
    call openinfile
    call openoutput
    call warning
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          write(outfile,'(a)') "#undef "//trim(name)
       endif
    enddo
    close (outfile)
    close(infile)
  end subroutine umakespaceparamsh

  subroutine makespacesed
    implicit none
    call openinfile
    call openoutput
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          write(outfile,'(a)') &
              "/\<s\.ind\>/d;/\<s\.inp\>/d;"//&
              "/^[ 	]*!/b;s/%"//trim(name)//"\>/___"//trim(name)//&
              "/g;s/\<"//trim(name)//"\>/s%"//trim(name)//"/g;s/___"//&
              trim(name)//"\>/%"//trim(name)//"/g"
       endif
    enddo
    close (outfile)
    close(infile)
  end subroutine makespacesed

  ! just for testing some things, not used in production version wwvv
  subroutine makespacedefine
    implicit none
    call openinfile
    call openoutput
    call warning
    write(outfile,'(a)'),'#ifndef SPACE_DEFINE_H'
    write(outfile,'(a)'),'#define SPACE_DEFINE_H'
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          write(outfile,'(a)') '#define '//trim(name)//' s%'//caseswitch(trim(name))
       endif
    enddo
    write(outfile,'(a)'),'#endif'
    close (outfile)
    close(infile)
  end subroutine makespacedefine

  subroutine makespaceundef
    implicit none
    call openinfile
    call openoutput
    call warning
    write(outfile,'(a)'),'#ifdef SPACE_DEFINE_H'
    write(outfile,'(a)'),'#undef SPACE_DEFINE_H'
    do
       call getline
       if (endfound) then
          exit
       endif
       call getitems
       if (varfound) then
          write(outfile,'(a)') '#undef '//trim(name)
       endif
    enddo
    write(outfile,'(a)'),'#endif'
    close (outfile)
    close(infile)
  end subroutine makespaceundef

  subroutine makeparamsinc
    implicit none
    integer             :: startpar,endpar,i
    logical             :: readtype
    character(slen)     :: parname

    infile   = 10
    outfile   = 11
    endfound = .false.
    readtype = .false.
    open(infile,file='params.F90')
    open(outfile,file='parameters.inc')
    write(outfile,'(a)')'!! This code is generated automatocally by makeincludes from params.F90'
    write(outfile,'(a)')'!! Do not change manually'
    write(outfile,'(a)')'subroutine outputparameters(par)'
    write(outfile,'(a)')'   use filefunctions'
    write(outfile,'(a)')'   implicit none'
    write(outfile,'(a)')'   type(parameters)   :: par'
    write(outfile,'(a)')'   integer            :: fid,i'
    write(outfile,'(a)')' '
    write(outfile,'(a)')'   fid=create_new_fid()'
    write(outfile,'(a)')'   open(fid,file=''params.dat'',status=''replace'', action=''write'')'
    do
       call getline
       if (endfound) then
          exit
       endif
       if (.not. readtype) then
          if (line(1:15) == 'type parameters') then
             readtype = .true.
          endif
       else
          ! look for end of type declaration
          if (line(1:19) == 'end type parameters') then
             readtype = .false.
             exit
          else  ! we're still reading in parameters
             if (line(1:1) .ne. '!') then ! not reading comment lines
                ! look for '::' and '='
                startpar = 0
                endpar = 0
                line=adjustl(line)
                do i=1,slen-1
                   if (line(i:i+1) == '::') then
                      startpar = i+2
                      exit
                   endif
                enddo
                do i=1,slen-1
                   if ((line(i:i) == '=') .and. (i>startpar) )then
                      endpar = i-1
                      exit
                   endif
                enddo
                if ((startpar .ne. 0).and.(endpar.ne.0)) then ! we have a valid parameter
                   dimensioned = .false.
                   parname = trim(adjustl(line(startpar:endpar)))
                   do i=1,startpar-1-9
                      if (line(i:i+8)=='dimension') then
                         dimensioned=.true.
                         exit
                      endif
                   enddo
                   if (line(1:9) =='character') then
                      if (dimensioned) then
                         if (parname(1:10)=='globalvars') then
                            write(outfile,'(a)')'   do i=1,par%nglobalvar'
                            write(outfile,'(a)')'      write(fid,*)trim(adjustl(par%globalvars(i)))'
                            write(outfile,'(a)')'   enddo'   
                         elseif(parname(1:9)=='meanvars') then
                            write(outfile,'(a)')'   do i=1,par%nmeanvar'
                            write(outfile,'(a)')'      write(fid,*)trim(adjustl(par%meanvars(i)))'
                            write(outfile,'(a)')'   enddo'   
                         elseif(parname(1:9)=='pointvars') then
                            write(outfile,'(a)')'   do i=1,par%npointvar'
                            write(outfile,'(a)')'      write(fid,*)trim(adjustl(par%pointvars(i)))'
                            write(outfile,'(a)')'   enddo' 
                         endif
                      else
                         write(outfile,'(a)')'   write(fid,*)'''//trim(parname)//'='',trim(adjustl(par%'//trim(parname)//'))'
                      endif
                   else  ! not character, i.e. real/integer/logical
                      if (dimensioned) then
                         if(parname(1:3)=='D50') then
                            write(outfile,'(a)')'   write(fid,*)'''//trim(parname)//'='',par%'//trim(parname)//'(1:par%ngd)'
                         elseif(parname(1:3)=='D90') then
                            write(outfile,'(a)')'   write(fid,*)'''//trim(parname)//'='',par%'//trim(parname)//'(1:par%ngd)'
                         elseif(parname(1:6)=='sedcal') then
                            write(outfile,'(a)')'   write(fid,*)'''//trim(parname)//'='',par%'//trim(parname)//'(1:par%ngd)'
                         elseif(parname(1:6)=='ucrcal') then
                            write(outfile,'(a)')'   write(fid,*)'''//trim(parname)//'='',par%'//trim(parname)//'(1:par%ngd)'
                         elseif(parname(1:6)=='ucrcal') then
                            write(outfile,'(a)')'   write(fid,*)'''//trim(parname)//'='',par%'//trim(parname)//'(1:par%ngd)'
                         endif
                      else
                         if(parname(1:7)=='npoints') then
                            write(outfile,'(a)')'   write(fid,*)''npoints='',trim(adjustl(par%npoints))'
                            write(outfile,'(a)')'   do i=1,par%npoints'
                            write(outfile,'(a)')'      write(fid,*)xpointsw(i),''  '',ypointsw(i)'
                            write(outfile,'(a)')'   enddo'
                         elseif(parname(1:8)=='nrugauge') then
                            write(outfile,'(a)')'   write(fid,*)''nrugauge='',trim(adjustl(par%nrugauge))'
                            write(outfile,'(a)')'   do i=1,par%nrugauge'
                            write(outfile,'(a)')'      write(fid,*)xpointsw(i+par%npoints),''  '',ypointsw(i+par%npoints)'
                            write(outfile,'(a)')'   enddo'
                         else   

                            write(outfile,'(a)')'   write(fid,*)'''//trim(parname)//'='',par%'//trim(parname)
                         endif
                      endif
                   endif
                endif
             endif
          endif
       endif
    enddo
    write(outfile,'(a)')'   close(fid)'
    write(outfile,'(a)')'end subroutine'
    close(infile)
    close(outfile)

  end subroutine makeparamsinc

  subroutine makegetkeygen
    implicit none
    integer             :: startpar,endpar,i, ncharacterkeys, nintegerkeys, nrealkeys
    logical             :: readtype
    character(slen)     :: parname
    integer, parameter  :: maxnvar = 1024
    character(slen), dimension(maxnvar) :: characterkeys, integerkeys, realkeys

    ncharacterkeys = 0
    nrealkeys      = 0
    nintegerkeys   = 0
    infile         = 10
    outfile        = 11
    endfound       = .false.
    readtype       = .false.
    open(infile,file='params.F90')
    do
       call getline
       if (endfound) then
          exit
       endif
       if (line(1:15) == 'type parameters') then
          readtype = .true.
          ! look for end of type declaration
       elseif (line(1:19) == 'end type parameters') then
          readtype = .false.
          exit
       end if
       ! If we're not reading a type, just continue to the next line
       if (.not.readtype) then
          cycle
       end if
       if (line(1:1) .eq. '!') then ! not reading comment lines
          cycle
       end if


       ! look for '::' and '='. In between is the parameter name.
       startpar = 0
       endpar = 0
       line=adjustl(line)

       do i=1,slen-1
          if (line(i:i+1) == '::') then
             startpar = i+2
             exit
          endif
       enddo
       do i=1,slen-1
          if ((line(i:i) == '=') .and. (i>startpar) )then
             endpar = i-1
             exit
          endif
       enddo

       if ((startpar .eq. 0).or.(endpar.eq.0)) then ! we have an invalid parameter
          ! Let's return....
          
          cycle
       end if


       dimensioned = .false.
       parname = trim(adjustl(line(startpar:endpar)))
       do i=1,startpar-1-9
          if (line(i:i+8)=='dimension') then
             dimensioned=.true.
             exit
          endif
       enddo


       if (dimensioned) then
          if (line(1:9) =='character') then
             if (parname(1:10)=='globalvars') then
             elseif(parname(1:9)=='meanvars') then
             elseif(parname(1:9)=='pointvars') then
             endif
          else
             if(parname(1:3)=='D50') then
             elseif(parname(1:3)=='D90') then
             elseif(parname(1:6)=='ucrcal') then
             elseif(parname(1:6)=='sedcal') then
             elseif(parname(1:7)=='npoints') then
             elseif(parname(1:8)=='nrugauge') then
             end if
          endif
       else  ! not character, i.e. real/integer/logical
          if (line(1:9) == 'character') then
             ncharacterkeys = ncharacterkeys + 1
             characterkeys(ncharacterkeys) = trim(parname)
          elseif (line(1:4) ==  'real') then
             nrealkeys = nrealkeys + 1
             realkeys(nrealkeys) = trim(parname)
          elseif (line(1:7) == 'integer') then
             nintegerkeys = nintegerkeys + 1
             integerkeys(nintegerkeys) = trim(parname)
          end if

       endif


    enddo
    close(infile)

    ! Write the getkey.gen
    open(outfile,file='getkey.gen')

    write(outfile,*) 'integer, parameter :: ncharacterkeys=', ncharacterkeys
    write(outfile,*) 'integer, parameter :: nrealkeys=', nrealkeys
    write(outfile,*) 'integer, parameter :: nintegerkeys=', nintegerkeys
    ! key definitions
    write(outfile,*) 'character(slen), dimension(ncharacterkeys) :: characterkeys'
    write(outfile,*) 'character(slen), dimension(nrealkeys) :: realkeys'
    write(outfile,*) 'character(slen), dimension(nintegerkeys) :: integerkeys'
    ! value definitions
    write(outfile,*) 'character(slen), dimension(ncharacterkeys) :: charactervalues'
    write(outfile,*) 'real*8, dimension(nrealkeys) :: realvalues'
    write(outfile,*) 'integer*4, dimension(nintegerkeys) :: integervalues'

    do i=1,ncharacterkeys
       write(outfile,*) 'characterkeys(', i, ') = "', trim(characterkeys(i)),  '"'
    end do
    do i=1,ncharacterkeys
       write(outfile,*) 'charactervalues(', i, ') = par%', trim(characterkeys(i))
    end do
    do i=1,nrealkeys
       write(outfile,*) 'realkeys(', i, ') = "', trim(realkeys(i)),  '"'
    end do
    do i=1,nrealkeys
       write(outfile,*) 'realvalues(',i, ') = par%', trim(realkeys(i))
    end do
    do i=1,nintegerkeys
       write(outfile,*) 'integerkeys(', i, ') = "', trim(integerkeys(i)),  '"'
    end do
    do i=1,nintegerkeys
       write(outfile,*) 'integervalues(', i, ') = par%', trim(integerkeys(i))
    end do
    close(outfile)

  end subroutine makegetkeygen


end module makemodule

program makeincludes
  use makemodule
  implicit none
  character(slen) :: command
  !
  ! count number of variables
  !
  call counting

  read (*,'(a)',end=100) command
100 continue
  if (trim(command) .eq. ' ' .or. ichar(command(1:1)) .eq. 0) then
    command = trim(spacedeclname)//' '// &
         trim(mnemonicname)//' '// &
         trim(indextosname)//' '// &
  !       trim(space_alloc_scalarsname)//' '// &
         trim(space_alloc_arraysname)//' '// &
         trim(space_indname)//' '// &
         trim(chartoindexname)//' '// &
         trim(space_inpname)//' '// &
         trim(spaceparams_hname)//' '// &
         trim(uspaceparams_hname)//' '// &
         trim(space_sedname)//' '// &
         trim(paramsincname)//' '// &
         trim(getkeygenname)
   endif

  if (index(command,trim(spacedeclname)) .ne. 0) then
     write(*,*)'Making '//trim(spacedeclname)
     outputfilename = spacedeclname
     call makespacedecl
  endif

  if (index(command,trim(mnemonicname)) .ne. 0) then
     write(*,*)'Making '//trim(mnemonicname)
     outputfilename = mnemonicname
     call makemnemonic
  endif

  if (index(command,trim(indextosname)) .ne. 0) then
     write(*,*)'Making '//trim(indextosname)
     outputfilename = indextosname
     call makeindextos
  endif

  if (index(command,trim(space_alloc_arraysname)) .ne. 0) then
     write(*,*)'Making '//trim(space_alloc_arraysname)
     outputfilename = space_alloc_arraysname
     call makespace_alloc_arrays
  endif

  if (index(command,trim(space_indname)) .ne. 0) then
     write(*,*)'Making '//trim(space_indname)
     outputfilename = space_indname
     call makespaceind
  endif

  if (index(command,trim(space_inpname)) .ne. 0) then
     write(*,*)'Making '//trim(space_inpname)
     outputfilename = space_inpname
     call makespaceinp
  endif

  if (index(command,trim(spaceparams_hname)) .ne. 0) then
     write(*,*)'Making '//trim(spaceparams_hname)
     outputfilename = spaceparams_hname
     call makespaceparamsh
  endif

  if (index(command,trim(uspaceparams_hname)) .ne. 0) then
     write(*,*)'Making '//trim(uspaceparams_hname)
     outputfilename = uspaceparams_hname
     call umakespaceparamsh
  endif

  if (index(command,trim(space_sedname)) .ne. 0) then
     write(*,*)'Making '//trim(space_sedname)
     outputfilename = space_sedname
     call makespacesed
  endif

  if (index(command,trim(chartoindexname)) .ne. 0) then
     write(*,*)'Making '//trim(chartoindexname)
     outputfilename = chartoindexname
     call makechartoindex
  endif

  if (index(command,trim(paramsincname)) .ne. 0) then
     write(*,*)'Making '//trim(paramsincname)
     outputfilename = paramsincname
     call makeparamsinc
  endif

  if (index(command,trim(getkeygenname)) .ne. 0) then
     write(*,*)'Making '//trim(getkeygenname)
     outputfilename = getkeygenname
     call makegetkeygen
  endif

  if (index(command,trim(space_definename)) .ne. 0) then
     write(*,*)'Making '//trim(space_definename)
     outputfilename = space_definename
     call makespacedefine
  endif

  if (index(command,trim(space_undefname)) .ne. 0) then
     write(*,*)'Making '//trim(space_undefname)
     outputfilename = space_undefname
     call makespaceundef
  endif

end program makeincludes

