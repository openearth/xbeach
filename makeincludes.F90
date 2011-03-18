! this program makes includefiles from spaceparams.tmpl
! indextos.gen            is for spaceparams.F90
! mnemonic.gen            is for mnemonic.F90
! chartoindex.gen         is for mnemonic.F90
! space_alloc_arrays.gen  is for spaceparams.F90
! space_alloc_scalars.gen is for spaceparams.F90
! spacedecl.gen           is for spaceparams.F90
! space_ind.gen           is for s.ind
! space_inp.gen           is for s.inp
! parameters.inc          is for xbeach.F90 after params.F90
!
! see also README.includefiles
!          spaceparams.tmpl
!
module makemodule
integer, parameter             :: slen = 1024
character(len=slen), parameter :: templatefilename = 'spaceparams.tmpl'
character(len=slen)            :: outputfilename
integer             :: infile = 1
integer             :: outfile = 2
character(len=slen), parameter :: spacedeclname = 'spacedecl.gen'
character(len=slen), parameter :: mnemonicname = 'mnemonic.gen'
character(len=slen), parameter :: indextosname = 'indextos.gen'
character(len=slen), parameter :: space_alloc_scalarsname = 'space_alloc_scalars.gen'
character(len=slen), parameter :: space_alloc_arraysname = 'space_alloc_arrays.gen'
character(len=slen), parameter :: space_indname = 'space_ind.gen'
character(len=slen), parameter :: space_inpname = 'space_inp.gen'
character(len=slen), parameter :: chartoindexname = 'chartoindex.gen'
integer                        :: numvars, maxnamelen, inumvars0, rnumvars0
character(len=slen)            :: type, name, comment, line, broadcast, description, units
integer                        :: rank
logical                        :: varfound
logical                        :: endfound
logical                        :: dimensioned
character(20), dimension(10) :: dims       ! dimensions are maximum 10 length statement, maximum 10 different dimensions
integer                        :: maxdimlen, maxrank
character(6),parameter         :: sp = '      '

contains
!
!  reads next line
!
subroutine getline
  implicit none
  character(slen) :: s
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

    s = adjustl(line)
    if (s(1:2) .ne. '!*') then
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

  character(len=slen)  :: w
  integer              :: j,l,i

  l=1
  type = ''
  name = ''
  units = ''
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
  character(len=slen) nullstr,rankstr
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
      nullstr = ' =>NULL() '
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
      write(outfile,'(a)') sp//trim(type)//trim(rankstr)//' :: '//trim(name)// &
           trim(nullstr)//'  '//trim(comment)
    else
      write(outfile,'(a)') trim(comment)
    endif
  enddo
  call format_fortran
  close (outfile)
  close(infile)
end subroutine makespacedecl
!
! construct subroutine to initialize the scalars
!
subroutine makespace_alloc_scalars
  implicit none
  integer ii,rr
  ii=1
  rr=1
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
      if (rank .eq. 0) then
        write(outfile,'(a)') sp//'allocate(s%'//trim(name)//')'
      endif
    endif
  enddo

  call format_fortran
  close(outfile)
  close(infile)

end subroutine makespace_alloc_scalars
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
  character(len=slen) :: rankstr
  character(len=10) :: extrastr
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
            case(1)
              rankstr = 't%r1'
            case(2)
              rankstr = 't%r2'
            case(3)
              rankstr = 't%r3'
            case(4)
              rankstr = 't%r4'
            case default
              rankstr = 'wrong dimensionality'
          end select
          case ('i')
          select case(rank)
            case(0)
              rankstr = 't%i0'
            case(1)
              rankstr = 't%i1'
            case(2)
              rankstr = 't%i2'
            case(3)
              rankstr = 't%i3'
            case(4)
              rankstr = 't%i4'
            case default
              rankstr = 'wrong dimensionality'
          end select
         end select
        write(outfile,'(a,i3,a)') sp//'  case(',j,')'
        write(outfile,'(a)')      sp//'    '//trim(rankstr)//'  => s%'//trim(name)
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
              write(outfile,'(a)',advance='no')   "'"//(dims(i))//"'"//trim(extrastr)
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
  character(len=slen) rankstr
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

subroutine makeparamsinc
   implicit none
   integer             :: startpar,endpar,i
   logical             :: readtype
   character(len=slen) :: parname

   infile = 10
   outfile = 11
   endfound = .false.
   readtype = .false.
   open(infile,file='params.F90')
   open(outfile,file='parameters.inc')
   write(outfile,'(a)')'!! This code is generated automatocally by makeincludes'
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



end module makemodule

program makeincludes
  use makemodule
  implicit none
  character(len=slen) :: command
  !
  ! count number of variables
  !
   call counting

  read (*,'(a)',end=100,err=100) command
  goto 110
  100 continue
      command = trim(spacedeclname)//' '// &
                trim(mnemonicname)//' '// &
                trim(indextosname)//' '// &
                trim(space_alloc_scalarsname)//' '// &
                trim(space_alloc_arraysname)//' '// &
                trim(space_indname)//' '// &
                trim(chartoindexname)//' '// &
                trim(space_inpname)
  110 continue
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

   if (index(command,trim(space_alloc_scalarsname)) .ne. 0) then
     write(*,*)'Making '//trim(space_alloc_scalarsname)
     outputfilename = space_alloc_scalarsname
     call makespace_alloc_scalars
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

   if (index(command,trim(chartoindexname)) .ne. 0) then
     write(*,*)'Making '//trim(chartoindexname)
     outputfilename = chartoindexname
     call makechartoindex
   endif
   
   call makeparamsinc
   
end program makeincludes

