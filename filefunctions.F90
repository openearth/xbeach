module filefunctions

  interface check_file_length
     module procedure check_file_length_1D
     module procedure check_file_length_2D
     module procedure check_file_length_3D
  end interface check_file_length

contains

  integer function create_new_fid()
    use xmpi_module
    use logging_module

    integer    :: fileunit

    fileunit = -1 ! temporary
    if (xmaster) then 
       fileunit = create_new_fid_generic()
       if (fileunit==-1) then
          call writelog('les','','Serious problem: not enough free unit ids to create new file')
          call halt_program
       endif
    endif
    create_new_fid = fileunit   
  end function create_new_fid

  subroutine check_file_exist(filename)
    use xmpi_module
    use logging_module

    implicit none

    character(*)               :: filename
    integer                    :: error

    error = 0
    if (xmaster) call check_file_exist_generic(filename,error)

    if (error==1) then
       if (xmaster) then
          call writelog('sle','','File ''',trim(filename),''' not found. Terminating simulation')
       endif
	   call halt_program
    endif
  end subroutine check_file_exist



  subroutine check_file_length_1D(fname,d1)
    use xmpi_module
    use logging_module

    implicit none
    character(*)                   ::  fname
    integer, intent(in)            ::  d1
    integer                        ::  fid,iost
    integer                        ::  i
    real,dimension(:),allocatable  ::  dat
    character(256)                   ::  msg
    
    if (xmaster) then
	   allocate(dat(d1))
       fid = create_new_fid()
	   open(fid,file=trim(fname))
       read(fid,*,iostat=iost,iomsg=msg)(dat(i),i=1,d1)
       if (iost .ne. 0) then
          call writelog('sle','','Error processing file ''',trim(fname),'''. File may be too short or contains invalid values.', & 
               ' Terminating simulation' )
          call writelog('sle','', msg)
          call halt_program
       endif
       close(fid)
	   deallocate(dat)
	endif
    
  end subroutine check_file_length_1D

  subroutine check_file_length_2D(fname,d1,d2)
    use xmpi_module
    use logging_module

    implicit none
    character(*)                     :: fname
    integer, intent(in)              :: d1,d2
    integer                          :: fid,iost
    integer                          :: i,j
    real,dimension(:,:),allocatable  :: dat
    character(256)                   ::  msg

   
    if (xmaster) then 
	   allocate(dat(d1,d2))
       fid = create_new_fid()
       open(fid,file=trim(fname))
       read(fid,*,iostat=iost, iomsg=msg)((dat(i,j),i=1,d1),j=1,d2)
       if (iost .ne. 0) then
          call writelog('sle','','Error processing file ''',trim(fname),'''. File may be too short or contains invalid values.', & 
               ' Terminating simulation')
          call writelog('sle','',msg)
          call halt_program
       endif
       close(fid)
	   deallocate(dat)
    endif
  end subroutine check_file_length_2D

  subroutine check_file_length_3D(fname,d1,d2,d3)
    use xmpi_module
    use logging_module

    implicit none
    character(*)                       ::  fname
    integer, intent(in)                ::  d1,d2,d3
    integer                            ::  fid,iost
    integer                            ::  i,j,k
    real,dimension(:,:,:),allocatable  ::  dat
    character(256)                   ::  msg

    
    if (xmaster) then 
	   allocate(dat(d1,d2,d3))
       fid = create_new_fid()
       open(fid,file=trim(fname))
       read(fid,*,iostat=iost,iomsg=msg)(((dat(i,j,k),i=1,d1),j=1,d2),k=1,d3)
       if (iost .ne. 0) then
          call writelog('esl','', msg)
          call halt_program
       endif
	   close(fid)
	   deallocate(dat)
    endif    
  end subroutine check_file_length_3D

  subroutine checkbcfilelength(tstop,instat,filename,filetype)
    use logging_module
    use xmpi_module

    IMPLICIT NONE
    real*8, intent(in) :: tstop
    character(24), intent(in):: instat
    character*80      :: filename,dummy
    character*8       :: testc
    character*1       :: ch
    integer           :: i,ier=0,nlines,filetype,fid
    real*8            :: t,dt,total,d1,d2,d3,d4,d5


    if (xmaster) then 
       ier = 0
       fid = create_new_fid()
       open(fid,file=trim(filename))
       i=0
       do while (ier==0)
          read(fid,'(a)',iostat=ier)ch
          if (ier==0)i=i+1
       enddo
       nlines=i
       rewind(fid)    

       if (trim(instat)=='jons' .or. trim(instat)=='swan' .or. trim(instat)=='vardens') then 
          read(fid,*)testc
          if (testc=='FILELIST') then
             filetype = 1
             nlines=nlines-1
          else
             filetype = 0
          endif
       elseif (trim(instat)=='stat_table' .or. trim(instat)=='jons_table') then
          filetype = 2
       endif

       total=0.d0
       i=0
       select case (filetype)
       case(0)
          total=2.d0*tstop
       case(1)
          do while (total<tstop .and. i<nlines)
             read(fid,*)t,dt,dummy
             total=total+t
             i=i+1
          enddo
       case(2)
          do while (total<tstop .and. i<nlines)
             read(fid,*)d1,d2,d3,d4,d5,t,dt
             total=total+t
             i=i+1
          enddo
       end select
       close(fid)
       if (total<tstop) then
          call writelog('sle',' ','Error: Wave boundary condition time series too short in ',trim(filename))
          call writelog('sle','(a,f0.2,a,f0.2)',' Total wave condition time series is ',total, ' but simulation length is ',tstop)
          call writelog('sle',' ','Stopping calculation')
          call halt_program
       endif
    endif ! xmaster

  end subroutine checkbcfilelength

  subroutine check_file_exist_generic(filename,error)
    implicit none

    character(*)               :: filename
    integer                    :: error
    logical                    :: file_exists

    inquire(file=filename,exist=file_exists)

    error = 0

    if (.not. file_exists) then
       error = 1
    endif

  end subroutine check_file_exist_generic


  integer function create_new_fid_generic()
    integer    :: tryunit = 98
    logical    :: fileopen

    fileopen = .true.    
    do while (fileopen)
       inquire(tryunit,OPENED=fileopen)
       if (fileopen) then
          tryunit=tryunit-1
       endif
       if (tryunit<=10) then 
          tryunit = -1
          fileopen = .false.
       endif
    enddo
    create_new_fid_generic = tryunit   
  end function create_new_fid_generic

end module filefunctions
