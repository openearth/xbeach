module filefunctions
   use typesandkinds
   use paramsconst
   implicit none
   save
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
      fileunit = create_new_fid_generic()
      if (fileunit==-1) then
         call writelog('les','','Serious problem: not enough free unit ids to create new file')
         call halt_program
      endif
      create_new_fid = fileunit
   end function create_new_fid

   subroutine check_file_exist(filename,exist,forceclose)
      use xmpi_module
      use logging_module

      implicit none

      character(*)               :: filename
      logical,intent(out),optional :: exist
      logical,intent(in), optional :: forceclose
      logical                    :: endsim
      integer                    :: error

      if (present(forceclose)) then
         endsim = forceclose
      else
         endsim = .true.
      endif

      error = 0
      if (xmaster) call check_file_exist_generic(filename,error)

      if (error==1 .and. endsim) then
         if (xmaster) then
            call writelog('sle','','File ''',trim(filename),''' not found. Terminating simulation')
         endif
         call halt_program
      endif

      if (present(exist)) then
         if (error==1) then
            exist = .false.
         else
            exist = .true.
         endif
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

      if (xmaster) then
         allocate(dat(d1))
         fid = create_new_fid()
         open(fid,file=trim(fname))
         read(fid,*,iostat=iost)(dat(i),i=1,d1)
         if (iost .ne. 0) then
            call writelog('sle','','Error processing file ''',trim(fname),'''. File may be too short or contains invalid values.', &
            ' Terminating simulation' )
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


      if (xmaster) then
         allocate(dat(d1,d2))
         fid = create_new_fid()
         open(fid,file=trim(fname))
         read(fid,*,iostat=iost)((dat(i,j),i=1,d1),j=1,d2)
         if (iost .ne. 0) then
            call writelog('sle','','Error processing file ''',trim(fname),'''. File may be too short or contains invalid values.',&
            ' Terminating simulation')
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


      if (xmaster) then
         allocate(dat(d1,d2,d3))
         fid = create_new_fid()
         open(fid,file=trim(fname))
         read(fid,*,iostat=iost)(((dat(i,j,k),i=1,d1),j=1,d2),k=1,d3)
         if (iost .ne. 0) then
            call writelog('esl','Error processing file ''',trim(fname),'''. File may be too short or contains invalid values.', &
            ' Terminating simulation')
            call halt_program
         endif
         close(fid)
         deallocate(dat)
      endif
   end subroutine check_file_length_3D

   subroutine checkbcfilelength(tstop,instat,filename,filetype,nonh)
      use logging_module
      use xmpi_module

      IMPLICIT NONE
      type fileinfo
         character(slen)  :: fname
         integer          :: nlines
      end type

      real*8, intent(in) :: tstop
      integer, intent(in):: instat
      character(slen)     :: filename,dummy
      character(slen)     :: testc
      character(len=1)    :: ch
      integer           :: i,ier=0,nlines,filetype,fid,nlocs,ifid,fid2
      real*8            :: t,dt,total,d1,d2,d3,d4,d5
      type(fileinfo),dimension(:),allocatable :: bcfiles
      logical,intent(in),optional :: nonh
      logical                     :: lnonh

      if (present(nonh)) then
         lnonh=nonh
      else
         lnonh = .false.
      endif

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

         ! test for multiple locations setting
         read(fid,*,iostat=ier)testc
         if (ier .ne. 0) then
            call report_file_read_error(filename)
         endif
         ! wwvv fid2 was not initialized, so:
         fid2=create_new_fid()
         if (trim(testc)=='LOCLIST') then
            nlocs = nlines-1
            allocate(bcfiles(nlocs))
            do ifid = 1,nlocs
               read(fid,*,iostat=ier)d1,d2,bcfiles(ifid)%fname
               if (ier .ne. 0) then
                  call report_file_read_error(filename)
               endif
               call check_file_exist(trim(bcfiles(ifid)%fname))
               open(fid2,file=trim(bcfiles(ifid)%fname))
               i=0
               ier = 0
               do while (ier==0)
                  read(fid2,'(a)',iostat=ier)ch
                  if (ier==0)i=i+1
               enddo
               close(fid2)
               bcfiles(ifid)%nlines=i
            enddo
         else
            nlocs = 1
            allocate(bcfiles(1))
            bcfiles(1)%fname = filename
            bcfiles(1)%nlines = nlines
         endif
         close(fid)

         do ifid=1,nlocs
            fid = create_new_fid()
            open(fid,file=trim(bcfiles(ifid)%fname))
            if (instat==INSTAT_JONS .or. instat==INSTAT_SWAN .or. instat==INSTAT_VARDENS) then
               read(fid,*,iostat=ier)testc
               if (ier .ne. 0) then
                  call report_file_read_error(bcfiles(ifid)%fname)
               endif
               if (trim(testc)=='FILELIST') then
                  filetype = 1
                  bcfiles(ifid)%nlines=bcfiles(ifid)%nlines-1
               else
                  filetype = 0
               endif
            elseif (instat==INSTAT_STAT_TABLE .or. instat==INSTAT_JONS_TABLE) then
               filetype = 2
            elseif (instat==INSTAT_REUSE) then
               filetype = 3
            endif

            total=0.d0
            i=0
            select case (filetype)
             case(0)
               total=2.d0*tstop
             case(1)
               do while (total<tstop .and. i<bcfiles(ifid)%nlines)
                  read(fid,*,iostat=ier)t,dt,dummy
                  if (ier .ne. 0) then
                     call report_file_read_error(bcfiles(ifid)%fname)
                  endif
                  total=total+t
                  i=i+1
                  call check_file_exist(trim(dummy))
               enddo
             case(2)
               do while (total<tstop .and. i<bcfiles(ifid)%nlines)
                  read(fid,*,iostat=ier)d1,d2,d3,d4,d5,t,dt
                  if (ier .ne. 0) then
                     call report_file_read_error(bcfiles(ifid)%fname)
                  endif
                  total=total+t
                  i=i+1
               enddo
             case (3)
               do while (total<tstop .and. i<bcfiles(ifid)%nlines)
                  if (lnonh) then
                     read(fid,*,iostat=ier)d1,total,d2,dummy
                  else
                     read(fid,*,iostat=ier)total,d2,d3,d4,d5,dummy
                  endif
                  if (ier .ne. 0) then
                     call report_file_read_error(bcfiles(ifid)%fname)
                  endif
                  call check_file_exist(trim(dummy))
                  i=i+1
               enddo
            end select
            close(fid)
            if (total<tstop) then
               call writelog('sle',' ','Error: Wave boundary condition time series too short in ',trim(bcfiles(ifid)%fname))
               call writelog('sle','(a,f0.2,a,f0.2)',' Total wave condition time series is ',total, &
               ' but simulation length is ',tstop)
               call writelog('sle',' ','Stopping calculation')
               call halt_program
            endif
         enddo ! nlocs
      endif ! xmaster

   end subroutine checkbcfilelength

   function get_file_length(filename) result (n)

      implicit none

      character(slen), intent(in)             :: filename
      integer                                 :: n
      integer                                 :: io, error
      real*8                                  :: temp

      n   = 0
      io  = 0

      if (filename==' ') then
         n = 0
      else
         call check_file_exist_generic(filename, error)

         if (error == 1) then
            n = 0
         else
            open(11,file=filename)
            do while (io==0)
               n = n + 1
               read(11,*,IOSTAT=io) temp
            enddo
            close(11)
            n = n - 1
         endif
      endif

   end function get_file_length

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
      integer    :: tryunit = 9999
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
