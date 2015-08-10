module spaceparams
   use spaceparamsdef
   use xmpi_module
   use mnemmodule
   implicit none
   save
   !

#ifdef USEMPI

   ! This interface an cause complains by memcache because the it is conditional on uninitialized variables.
   ! I don't think this is an error
   interface space_distribute
      module procedure space_distribute_matrix_real8
      module procedure space_distribute_matrix_integer
      module procedure space_distribute_block_real8
      module procedure space_distribute_block_integer
      module procedure space_distribute_block4_real8
      module procedure space_distribute_vector
      module procedure space_distribute_block_vector
   end interface space_distribute

   interface space_shift_borders
      module procedure space_shift_borders_matrix_real8
      module procedure space_shift_borders_block_real8
   end interface space_shift_borders

   interface space_collect
      module procedure space_collect_block_real8
      module procedure space_collect_block_integer
      module procedure space_collect_block4_real8
      module procedure space_collect_block4_integer
      module procedure space_collect_matrix_real8
      module procedure space_collect_matrix_integer
   end interface space_collect

#endif

contains

  subroutine indextos(s,index,t)
    !
    ! given s and index, this subroutine returns in t
    ! a pointer to the requested array
    !
    use mnemmodule
    use logging_module
    use spaceparamsdef
    implicit none
    type (spacepars), intent(in),target :: s
    integer, intent(in)                 :: index
    type(arraytype), intent(out)        :: t
      
    if (index .lt. 1 .or. index .gt. numvars) then
       call writelog('els','(a,i3,a)','invalid index ',index,' in indextos. Program will stop')
       call halt_program
    endif
    
    select case(index)
       include 'indextos.inc'
    end select
    
  end subroutine indextos
  
  subroutine index_allocate(s,par,index,choice)
    ! allocates, deallocates reallocates in type s, based on index
    ! choice:
    !   'a': allocate
    !   'd': deallocate
    !   'r': reallocate
    use params
    use logging_module
    use spaceparamsdef
    implicit none
    type (spacepars), intent(inout) :: s
    type(parameters), intent(in)    :: par
    integer, intent(in)             :: index
    character(*)                    :: choice
    
    ! the .inc files contain code for allocatable entities
    ! scalars will be skipped silently
    select case(choice(1:1))
    case('a','A')
       select case(index)
       case default
          continue
          include 'index_allocate.inc'
       end select
    case ('d','D')
       select case(index)
       case default
          continue
          include 'index_deallocate.inc'
       end select
    case ('r','R')
       select case(index)
       case default
          continue
          include 'index_reallocate.inc'
       end select
    end select
  end subroutine index_allocate
  
  logical function index_allocated(s,index)
    use logging_module
    use spaceparamsdef
    implicit none
    type (spacepars), intent(in) :: s
    integer, intent(in)          :: index
    
    logical                      :: r
    
    r = .true. ! for scalars: function will return always .true.
    select case(index)
    case default
       continue
       include 'index_allocated.inc'
    end select
    
    index_allocated = r
    
  end function index_allocated

  ! Generated subroutine to allocate the scalars in s
  subroutine space_alloc_scalars(s)
    use mnemmodule
    implicit none
    type(spacepars),intent(inout)  :: s

    include 'space_alloc_scalars.inc'

  end subroutine space_alloc_scalars

  ! Generated subroutine to allocate all arrays in s
  subroutine space_alloc_arrays(s,par)
    use mnemmodule
    use params
    implicit none
    type(spacepars),intent(inout)  :: s
    type(parameters),intent(in)    :: par
    
    include 'space_alloc_arrays.inc'
    
  end subroutine space_alloc_arrays
  
#ifdef USEMPI

   ! copies scalars from sg to sl on xmaster, and distributes
   ! them
   subroutine space_copy_scalars(sg,sl)
      use mnemmodule
      implicit none
      type(spacepars),intent(inout)  :: sg,sl

      type(arraytype)                :: tg,tl
      integer                        :: j

      do j = 1,numvars
         call indextos(sg,j,tg)
         if (tg%rank .eq. 0) then
            call indextos(sl,j,tl)
            select case (tg%type)
             case('i')
               tl%i0 = tg%i0
             case('r')
               tl%r0 = tg%r0
            end select
         endif
      enddo
   end subroutine space_copy_scalars

   ! The scalars are needed for allocating the arrays,
   ! we distribute them here:
   subroutine space_distribute_scalars(sl)
      use mnemmodule
      use xmpi_module
      type (spacepars) :: sl

      type (arraytype) :: tl
      integer          :: i
      logical, parameter :: toall = .true.

      do i=1,numvars
         call indextos(sl,i,tl)
         if (tl%rank .eq. 0) then
            if (tl%type .eq. 'i') then
               call xmpi_bcast(tl%i0,toall)
            else
               call xmpi_bcast(tl%r0,toall)
            endif
         endif
      enddo

   end subroutine space_distribute_scalars
   !________________________________________________________________________
   subroutine space_distribute_matrix_real8(sl,a,b,toall)
      use xmpi_module
      use general_mpi_module
      implicit none
      type (spacepars), intent(inout)     :: sl
      real*8, dimension(:,:), intent(in)  :: a
      real*8, dimension(:,:), intent(out) :: b
      include 'space_distribute.inc'
      call matrix_distr(a,b,sl%is,sl%lm,sl%js,sl%ln,src,comm)
   end subroutine space_distribute_matrix_real8

   !________________________________________________________________________
   subroutine space_distribute_matrix_integer(sl,a,b,toall)
      use xmpi_module
      use general_mpi_module
      implicit none
      type (spacepars), intent(inout)      :: sl
      integer, dimension(:,:), intent(in)  :: a
      integer, dimension(:,:), intent(out) :: b
      include 'space_distribute.inc'

      call matrix_distr(a,b,sl%is,sl%lm,sl%js,sl%ln,src,comm)

   end subroutine space_distribute_matrix_integer

   !________________________________________________________________________
   subroutine space_distribute_block_real8(sl,a,b,toall)
      use xmpi_module
      use general_mpi_module
      implicit none
      type (spacepars), intent(inout)                :: sl
      real*8, dimension(:,:,:), intent(in)           :: a
      real*8, dimension(:,:,:), intent(out)          :: b

      real*8, dimension(1,1)                 :: dum
      integer                                        :: i
      include 'space_distribute.inc'

      do i=1,size(b,3)   ! assuming that b is allocated on all processes
         if (iamsrc) then
            call matrix_distr(a(:,:,i),b(:,:,i),sl%is,sl%lm,sl%js,sl%ln,src,comm)
         else
            call matrix_distr(dum,     b(:,:,i),sl%is,sl%lm,sl%js,sl%ln,src,comm)
         endif
      enddo

   end subroutine space_distribute_block_real8
   !________________________________________________________________________

   subroutine space_distribute_block_integer(sl,a,b,toall)
      use xmpi_module
      use general_mpi_module
      implicit none
      type (spacepars), intent(inout)                :: sl
      integer, dimension(:,:,:), intent(in)          :: a
      integer, dimension(:,:,:), intent(out)         :: b

      integer, dimension(1,1)                :: dum
      integer                                        :: i
      include 'space_distribute.inc'

      do i=1,size(b,3)   ! assuming that b is allocated on all processes
         if(iamsrc) then
            call matrix_distr(a(:,:,i),b(:,:,i),sl%is,sl%lm,sl%js,sl%ln,src,comm)
         else
            call matrix_distr(dum,     b(:,:,i),sl%is,sl%lm,sl%js,sl%ln,src,comm)
         endif

      enddo

   end subroutine space_distribute_block_integer
   !________________________________________________________________________

   subroutine space_distribute_block4_real8(sl,a,b,toall)
      use xmpi_module
      use general_mpi_module
      implicit none
      type (spacepars), intent(inout)                :: sl
      real*8, dimension(:,:,:,:), intent(in)         :: a
      real*8, dimension(:,:,:,:), intent(out)        :: b

      integer                                        :: i,j
      real*8, dimension(1,1)                         :: dum
      include 'space_distribute.inc'

      do i=1,size(b,3)
         do j=1,size(b,4)
            if(iamsrc) then
               call matrix_distr(a(:,:,i,j),b(:,:,i,j),sl%is,sl%lm,sl%js,sl%ln,src,comm)
            else
               call matrix_distr(dum,       b(:,:,i,j),sl%is,sl%lm,sl%js,sl%ln,src,comm)
            endif
         enddo
      enddo

   end subroutine space_distribute_block4_real8
   !________________________________________________________________________

   subroutine space_distribute_vector(xy,sl,a,b,toall)
      use xmpi_module
      use general_mpi_module
      implicit none
      ! Not sure what this xy  is doing here.
      ! Ok, the master process holds a vector of length ny+1 or nx+1
      ! This vector will be distributed among all processes.
      ! if xy .eq. 'x', the distribution is done according
      ! to the values in is and lm, i.e. the global vector
      ! has a size equal to the size of the first dimension of the
      ! global area: nx+1.
      ! In the other case, the distribution is done according to
      ! the values in js and ln, i.e. the global vector has a size
      ! equal to the second dimension of the global area: ny+1
      !
      character, intent(in)                     :: xy
      type (spacepars), intent(inout), target   :: sl
      real*8, dimension(:), intent(in)          :: a
      real*8, dimension(:), intent(out)         :: b

      integer, dimension(:), pointer    :: ijs,lmn
      include 'space_distribute.inc'

      if(xy .eq.'x') then
         ijs => sl%is
         lmn => sl%lm
      else
         ijs => sl%js
         lmn => sl%ln
      endif
      call vector_distr(a,b,ijs,lmn,src,comm)

   end subroutine space_distribute_vector
   !________________________________________________________________________

   subroutine space_distribute_block_vector(xy,sl,a,b,toall)
      use xmpi_module
      use general_mpi_module
      implicit none
      character, intent(in)               :: xy
      type (spacepars), intent(in)        :: sl
      real*8, dimension(:,:), intent(in)  :: a
      real*8, dimension(:,:), intent(out) :: b

      !    integer                             :: i
      include 'space_distribute.inc'

      !DF  do i=1,sl%ntheta
      !    do i=1,size(b,2)
      !       call space_distribute(xy,sl,a(:,i),b(:,i))
      !    enddo

      select case(xy)
       case('y')
         call vector_distr(a,b,sl%js,sl%ln,src,comm)
       case default
         print *,'error in space_distribute_block_vector, other than "y" not tested'
         call xmpi_abort()
      end select

   end subroutine space_distribute_block_vector
   !________________________________________________________________________


   subroutine space_distribute_space(sg,sl,par)
      use xmpi_module
      use logging_module
      use general_mpi_module
      use params
      use mnemmodule

      implicit none
      type(spacepars), intent(inout)  :: sg
      type(spacepars), intent(inout)  :: sl
      type(parameters)                :: par

      integer                         :: i,j,lid,eid, wid
      real*8, pointer, dimension(:)   :: vectorg, vectorl
      type (arraytype)                :: tg, tl

      logical, parameter :: toall = .true.
      integer, parameter :: nbord = 2
      character(1000)    :: txt

      !
      ! This subroutine takes care that all contents of the global
      ! space is distributed to the local space.
      !

      ! allocate scalars

      !call space_alloc_scalars(sl)

      ! copy scalars to sl, only on master
      ! distributing will take place later


      if(xmaster) then
         call space_copy_scalars(sg,sl)
      endif

      ! copy scalars to all processes, nx and ny will be adapted later
      call space_distribute_scalars(sl)
      call space_distribute_scalars(sg)

      ! Also, the isleft, isright, istop and isbot logicals from
      ! the xmpi module are put in sg and sl.
      !

      !
      ! Distribute is,js,ln,lm,isleft,isright,istop,isbot
      !

      if(xmaster) then
         allocate(sg%is(xmpi_size))
         allocate(sg%js(xmpi_size))
         allocate(sg%lm(xmpi_size))
         allocate(sg%ln(xmpi_size))
         allocate(sg%isleft(xmpi_size))
         allocate(sg%isright(xmpi_size))
         allocate(sg%istop(xmpi_size))
         allocate(sg%isbot(xmpi_size))
      endif


      allocate(sl%is(xmpi_size))
      allocate(sl%js(xmpi_size))
      allocate(sl%lm(xmpi_size))
      allocate(sl%ln(xmpi_size))
      allocate(sl%isleft(xmpi_size))
      allocate(sl%isright(xmpi_size))
      allocate(sl%istop(xmpi_size))
      allocate(sl%isbot(xmpi_size))

      if(xmaster) then
         call det_submatrices(sg%nx+1, sg%ny+1, xmpi_m, xmpi_n, &
         sg%is, sg%lm, sg%js, sg%ln, &
         sg%isleft, sg%isright, sg%istop, sg%isbot)
         call writelog('l','','--------------------------------')
         call writelog('l','','MPI implementation: ')
         call writelog('sl','','Distribution of matrix on processors')
         call writelog('sl','',' proc   is   lm   js   ln')
         do i=1,xmpi_size
            call writelog('sl','(i5,i5,i5,i5,i5)',i-1,sg%is(i),sg%lm(i),sg%js(i),sg%ln(i))
         enddo
         call writelog('ls','',' proc   left right top bot')
         do i=1,xmpi_size
            call writelog('ls','',i-1,sg%isleft(i),sg%isright(i),sg%istop(i),sg%isbot(i))
         enddo
         call writelog('l','','--------------------------------')
      endif

      ! determine the computational regions, i.e. the parts of the matrices
      ! that are computed in the processes.
      ! In general, the 2 borders left, right, top and bottom are not
      ! computed but received from the neighbours.
      ! Exceptions for the rows and colomns of left, right, top or
      ! bottom processes.
      !
      ! icgs, icge: start and end indices in global coordinates 1st dimension
      ! jcgs: jcge: start and end indices in global coordinates 2nd dimension
      ! icls, icle: start and end indices in local coordinates 1st dimension
      ! jcls, jcle: start and end indices in local coordinates 2nd dimension

      if (xmaster) then
         allocate (sg%icgs(xmpi_size))
         allocate (sg%icge(xmpi_size))
         allocate (sg%jcgs(xmpi_size))
         allocate (sg%jcge(xmpi_size))
         allocate (sg%icls(xmpi_size))
         allocate (sg%icle(xmpi_size))
         allocate (sg%jcls(xmpi_size))
         allocate (sg%jcle(xmpi_size))
      endif

      allocate (sl%icgs(xmpi_size))
      allocate (sl%icge(xmpi_size))
      allocate (sl%jcgs(xmpi_size))
      allocate (sl%jcge(xmpi_size))
      allocate (sl%icls(xmpi_size))
      allocate (sl%icle(xmpi_size))
      allocate (sl%jcls(xmpi_size))
      allocate (sl%jcle(xmpi_size))

      if (xmaster) then
         do i = 1,xmpi_size
            sg%icgs(i) = sg%is(i) + nbord
            sg%icge(i) = sg%is(i) + sg%lm(i) - 1 - nbord
            sg%jcgs(i) = sg%js(i) + nbord
            sg%jcge(i) = sg%js(i) + sg%ln(i) - 1 - nbord
            sg%icls(i) = nbord + 1
            sg%icle(i) = sg%lm(i) - nbord
            sg%jcls(i) = nbord + 1
            sg%jcle(i) = sg%ln(i) - nbord

            if (sg%istop(i)) then
               sg%icgs(i) = sg%icgs(i) - nbord
               sg%icls(i) = sg%icls(i) - nbord
            endif

            if (sg%isbot(i)) then
               sg%icge(i) = sg%icge(i) + nbord
               sg%icle(i) = sg%icle(i) + nbord
            endif

            if (sg%isleft(i)) then
               sg%jcgs(i) = sg%jcgs(i) - nbord
               sg%jcls(i) = sg%jcls(i) - nbord
            endif

            if (sg%isright(i)) then
               sg%jcge(i) = sg%jcge(i) + nbord
               sg%jcle(i) = sg%jcle(i) + nbord
            endif

         enddo
      endif

      if (xmaster) then
         sl%is       = sg%is
         sl%js       = sg%js
         sl%lm       = sg%lm
         sl%ln       = sg%ln
         sl%isleft   = sg%isleft
         sl%isright  = sg%isright
         sl%istop    = sg%istop
         sl%isbot    = sg%isbot

         sl%icgs     = sg%icgs
         sl%icge     = sg%icge
         sl%jcgs     = sg%jcgs
         sl%jcge     = sg%jcge
         sl%icls     = sg%icls
         sl%icle     = sg%icle
         sl%jcls     = sg%jcls
         sl%jcle     = sg%jcle
      endif

      call xmpi_bcast(sl%is,toall)
      call xmpi_bcast(sl%js,toall)
      call xmpi_bcast(sl%lm,toall)
      call xmpi_bcast(sl%ln,toall)
      call xmpi_bcast(sl%isleft,toall)
      call xmpi_bcast(sl%isright,toall)
      call xmpi_bcast(sl%istop,toall)
      call xmpi_bcast(sl%isbot,toall)

      call xmpi_bcast(sl%icgs,toall)
      call xmpi_bcast(sl%icge,toall)
      call xmpi_bcast(sl%jcgs,toall)
      call xmpi_bcast(sl%jcge,toall)
      call xmpi_bcast(sl%icls,toall)
      call xmpi_bcast(sl%icle,toall)
      call xmpi_bcast(sl%jcls,toall)
      call xmpi_bcast(sl%jcle,toall)

      if (xmaster) then
         call writelog('l','','--------------------------------')
         call writelog('sl','','computational domains on processors')
         call writelog('sl','','   proc   icgs   icge   jcgs   jcge   icls   icle   jcls   jcle')
         do i=1,xmpi_size
            write(txt,'(9i7)')i-1,sg%icgs(i),sg%icge(i),sg%jcgs(i),sg%jcge(i), &
            sg%icls(i),sg%icle(i),sg%jcls(i),sg%jcle(i)
            call writelog('sl','',txt)
         enddo
         call writelog('l','','--------------------------------')
      endif

      !
      ! compute the values for local nx and ny
      !

      if(xcompute) then
         sl%nx = sl%lm(xmpi_rank+1) - 1
         sl%ny = sl%ln(xmpi_rank+1) - 1
      endif


      !
      ! allocate all arrays in sl
      !

      call space_alloc_arrays(sl,par)

      ! for each variable in sg, find out how to distribute it
      ! and distribute
      !

      do i = 1,numvars

         call indextos(sg,i,tg)
         call indextos(sl,i,tl)
         select case (tl%btype)
          case('b')           ! have to broadcast this
            select case (tl%type)
             case('i')
               select case(tl%rank)
                case(0) ! do nothing, scalars are already in place
                  ! Have to be prudent here. all scalars are broadcasted, exept nx
                  ! and ny
                  !if (tl%name .ne. mnem_nx .and. tl%name .ne. mnem_ny) then
                  !  if(xmaster) then
                  !    tl%i0 = tg%i0
                  !  endif
                  !  call xmpi_bcast(tl%i0)
                  !endif
                case(1)
                  if(xmaster) then
                     tl%i1 = tg%i1
                  endif
                  call xmpi_bcast(tl%i1,toall)
                case(2)   ! wwvv : try to fix error with integer pdish
                  if(xmaster) then
                     tl%i2 = tg%i2
                  endif
                  call xmpi_bcast(tl%i2,toall)
                case default
                  goto 100
               end select   ! rank
             case('r')
               select case(tl%rank)
                case(0) ! do nothing here, scalars are already in place
                  ! Have to be prudent here. all scalars are broadcasted, exept nx
                  ! and ny
                  ! here only real*8 is broadcasted, so no check necessary here
                  !if(xmaster) then
                  !  tl%r0 = tg%r0
                  !endif
                  !call xmpi_bcast(tl%r0)
                case(1)
                  if(xmaster) then
                     tl%r1 = tg%r1
                  endif
                  call xmpi_bcast(tl%r1,toall)
                case(2)
                  if(xmaster) then
                     tl%r2 = tg%r2
                  endif
                  call xmpi_bcast(tl%r2,toall)
                case default
                  goto 100
               end select   ! rank
             case default
               goto 100
            end select     ! type
          case('d')        ! have to distribute this
            select case(tl%type)
             case ('i')
               select case(tl%rank)
                case(1) ! wwvv 2013: superfluous, there is a default case
                  goto 100
                case(2)
                  call space_distribute(sl,tg%i2,tl%i2)
                case(3)
                  call space_distribute(sl,tg%i3,tl%i3)
                case default
                  goto 100
               end select  ! rank
             case ('r')
               select case(tl%rank)
                case(1)
                  !         Robert: these don't exist anymore?
                  !            case(mnem_xz,mnem_xu)
                  !              call space_distribute_vector('x',sl,tg%r1,tl%r1)
                  !            case(mnem_yz, mnem_yv, mnem_bi)
                  select case(tl%name)
                   case(mnem_bi)
                     call space_distribute_vector('y',sl,tg%r1,tl%r1)
                   case(mnem_runup)
                     ! Not sure why vector expects a name....
                     ! This name is x or y, it relates to the length of the vector.
                     call space_distribute_vector('y',sl,tg%r1,tl%r1)
                   case(mnem_hrunup)
                     ! Not sure why vector expects a name....
                     ! This name is x or y, it relates to the length of the vector.
                     call space_distribute_vector('y',sl,tg%r1,tl%r1)
                   case(mnem_xhrunup)
                     ! Not sure why vector expects a name....
                     ! This name is x or y, it relates to the length of the vector.
                     call space_distribute_vector('y',sl,tg%r1,tl%r1)
                   case(mnem_strucslope)
                     call space_distribute_vector('y',sl,tg%r1,tl%r1)
                   case(mnem_istruct)
                     call space_distribute_vector('y',sl,tg%r1,tl%r1)
                   case(mnem_iwl)
                     call space_distribute_vector('y',sl,tg%r1,tl%r1)
                   case default
                     goto 100
                  end select
                case(2)
                  call space_distribute(sl,tg%r2,tl%r2)
                case(3)
                  call space_distribute(sl,tg%r3,tl%r3)
                case(4)
                  call space_distribute(sl,tg%r4,tl%r4)
                case default
                  continue
               end select  ! rank
             case default
               goto 100
            end select       ! type
          case('2')
            ! the umean case
            ! dimension 2,s%ny+1
            if(xmaster) then
               allocate(vectorg(sg%ny+1))
            else
               allocate(vectorg(1))
               vectorg = 0
            endif
            allocate(vectorl(sl%ny+1))
            ! Let's keep this 2 case explicit. (because tg%r2 is only known on master)
            ! This one is distributed as 2 vectors.
            do j=1,2
               if (xmaster) then
                  vectorg = tg%r2(j,:)
               endif
               call space_distribute("y", sl,vectorg,vectorl)
               tl%r2(j,:) = vectorl
            enddo
            if(xmaster) then
               deallocate(vectorg)
            endif
            deallocate(vectorl)
          case default
            goto 100
         end select  ! btype
      enddo  ! numvars
      return

100   continue
      call writelog('sel','',xmpi_rank,': Error in space_distribute_space, trying to distribute:')
      call get_logfileid(lid,eid, wid)
      call printvar(tl,lid,eid, wid)
      call halt_program
      return

   end subroutine space_distribute_space

   subroutine space_who_has(sl,i,j,p)
      ! determine which process contains element(i,j)
      ! sl: local spacepars
      ! i,j: global indices
      ! p:   process in xmpi_ocomm which contains the element

      ! icgs, icge: start and end indices in global coordinates 1st dimension
      ! jcgs: jcge: start and end indices in global coordinates 2nd dimension
      implicit none
      type(spacepars), intent(in) :: sl
      integer, intent(in)         :: i,j
      integer, intent(out)        :: p

      integer k

      p = -123

      do k = 1,xmpi_size
         if (i .ge. sl%icgs(k) .and. i .le. sl%icge(k) .and. &
         &   j .ge. sl%jcgs(k) .and. j .le. sl%jcge(k)) then
            p = xmpi_rank_to_orank(k-1)
            exit
         endif
      enddo

   end subroutine space_who_has

   subroutine space_global_to_local(sl,ig,jg,il,jl)
      ! given global coordinates (ig,jg), compute local coordinates (il,jl)
      implicit none
      type(spacepars), intent(in) :: sl
      integer, intent(in)         :: ig,jg
      integer, intent(out)        :: il,jl

      il = ig - sl%is(xmpi_rank+1) + 1
      jl = jg - sl%js(xmpi_rank+1) + 1

   end subroutine space_global_to_local

   subroutine space_local_to_global(sl,il,jl,ig,jg)
      ! given local coordinates (il,jl), compute global coordinates (ig,jg)
      implicit none
      type(spacepars), intent(in) :: sl
      integer, intent(in)         :: il,jl
      integer, intent(out)        :: ig,jg

      ig = il + sl%is(xmpi_rank+1) - 1
      jg = jl + sl%js(xmpi_rank+1) - 1

   end subroutine space_local_to_global

   subroutine space_shift_borders_matrix_real8(a)
      use general_mpi_module
      use xmpi_module
      implicit none
      real*8, intent(inout),dimension(:,:) :: a
      call shift_borders_matrix_real8(a,xmpi_left,xmpi_right, &
      xmpi_top,xmpi_bot,xmpi_comm)
   end subroutine space_shift_borders_matrix_real8

   subroutine space_shift_borders_block_real8(a)
      use general_mpi_module
      use xmpi_module
      implicit none
      real*8, intent(inout),dimension(:,:,:) :: a

      integer i

      do i=1,size(a,3)
         call shift_borders_matrix_real8(a(:,:,i),xmpi_left,xmpi_right,&
         xmpi_top,xmpi_bot,xmpi_comm)
      enddo
   end subroutine space_shift_borders_block_real8

   !
   ! wwvv a subtle point with the collect subroutines: the second
   ! argument: the matrix wherein the submatrices are to be collected,
   ! does not have to be available on the non-master processes, so
   ! the dimensions are not defined.
   ! The third argument is always defined, on master and non-master
   ! processes, so its dimensions (notably the 3rd in the block subroutines
   ! are available
   !
   ! parameters of the space_collect routines:
   ! s: spacepars: LOCAL s
   ! a: output: in this matrix the submatrices are collected
   ! b: input:  the local submatrix

   subroutine space_collect_block_real8(s,a,b)
      use general_mpi_module
      use xmpi_module
      implicit none
      type(spacepars), intent(in)            :: s
      real*8, dimension(:,:,:), intent(out)  :: a
      real*8, dimension(:,:,:), intent(in)   :: b

      integer i
      real*8, dimension(1,1) :: dum

      do i = 1,size(b,3)
         if (xomaster) then
            call matrix_coll(a(:,:,i),b(:,:,i),s%is,s%lm,s%js,s%ln, &
            s%isleft,s%isright,s%istop,s%isbot, &
            xmpi_omaster,xmpi_ocomm)
         else
            call matrix_coll(dum,     b(:,:,i),s%is,s%lm,s%js,s%ln, &
            s%isleft,s%isright,s%istop,s%isbot, &
            xmpi_omaster,xmpi_ocomm)
         endif
      enddo

   end subroutine space_collect_block_real8

   subroutine space_collect_block_integer(s,a,b)
      use general_mpi_module
      use xmpi_module
      implicit none
      type(spacepars), intent(in)            :: s
      integer, dimension(:,:,:), intent(out)  :: a
      integer, dimension(:,:,:), intent(in)   :: b

      integer i

      !real*8, dimension(:,:,:), allocatable :: ra,rb
      ! wwvv changed this into
      ! wwvv huh?  TODO
      integer, dimension(:,:,:), allocatable :: ra,rb
      integer                             :: m,n,o
      integer, dimension(1,1)                :: dum

      m = size(b,1)
      n = size(b,2)
      o = size(b,3)

      allocate(rb(m,n,o))

      if (xomaster) then
         m = size(a,1)
         n = size(a,2)
         o = size(a,3)
         allocate(ra(m,n,o))
      else
         allocate(ra(1,1,1))
      endif

      rb = b

      do i = 1,o
         if(xomaster) then
            call matrix_coll(ra(:,:,i),rb(:,:,i),s%is,s%lm,s%js,s%ln, &
            s%isleft,s%isright,s%istop,s%isbot, &
            xmpi_omaster,xmpi_ocomm)
         else
            call matrix_coll(dum,      rb(:,:,i),s%is,s%lm,s%js,s%ln, &
            s%isleft,s%isright,s%istop,s%isbot, &
            xmpi_omaster,xmpi_ocomm)
         endif
      enddo

      if (xomaster) then
         a = ra
      endif

      deallocate(ra,rb)

   end subroutine space_collect_block_integer

   subroutine space_collect_block4_real8(s,a,b)
      use general_mpi_module
      use xmpi_module
      implicit none
      type(spacepars), intent(in)              :: s
      real*8, dimension(:,:,:,:), intent(out)  :: a
      real*8, dimension(:,:,:,:), intent(in)   :: b

      integer i,j
      real*8, dimension(1,1) :: dum

      do j = 1,size(b,4)
         do i = 1,size(b,3)
            if (xomaster) then
               call matrix_coll(a(:,:,i,j),b(:,:,i,j),s%is,s%lm,s%js,s%ln, &
               s%isleft,s%isright,s%istop,s%isbot, &
               xmpi_omaster,xmpi_ocomm)
            else
               call matrix_coll(dum,       b(:,:,i,j),s%is,s%lm,s%js,s%ln, &
               s%isleft,s%isright,s%istop,s%isbot, &
               xmpi_omaster,xmpi_ocomm)
            endif
         enddo
      enddo

   end subroutine space_collect_block4_real8

   subroutine space_collect_block4_integer(s,a,b)
      use general_mpi_module
      use xmpi_module
      implicit none
      type(spacepars), intent(in)              :: s
      integer, dimension(:,:,:,:), intent(out)  :: a
      integer, dimension(:,:,:,:), intent(in)   :: b

      integer i,j
      integer :: dum(1,1)

      do j = 1,size(b,4)
         do i = 1,size(b,3)
            if(xomaster) then
               call matrix_coll(a(:,:,i,j),b(:,:,i,j),s%is,s%lm,s%js,s%ln, &
               s%isleft,s%isright,s%istop,s%isbot, &
               xmpi_omaster,xmpi_ocomm)
            else
               call matrix_coll(dum,       b(:,:,i,j),s%is,s%lm,s%js,s%ln, &
               s%isleft,s%isright,s%istop,s%isbot, &
               xmpi_omaster,xmpi_ocomm)
            endif
         enddo
      enddo

   end subroutine space_collect_block4_integer

   subroutine space_collect_matrix_real8(s,a,b)
      use general_mpi_module
      use xmpi_module
      implicit none
      type(spacepars), intent(in)          :: s
      real*8, dimension(:,:), intent(out)  :: a
      real*8, dimension(:,:), intent(in)   :: b

      call matrix_coll(a,b,s%is,s%lm,s%js,s%ln, &
      s%isleft,s%isright,s%istop,s%isbot, &
      xmpi_omaster,xmpi_ocomm)

   end subroutine space_collect_matrix_real8

   subroutine space_collect_matrix_integer(s,a,b)
      use general_mpi_module
      use xmpi_module
      implicit none
      type(spacepars), intent(in)           :: s
      integer, dimension(:,:), intent(out)  :: a
      integer, dimension(:,:), intent(in)   :: b

      ! not used often, so we convert the integers to real*8,
      ! collect and convert back.
      ! if this routine becomes heavily used, than a special
      ! matrix_coll has to be made. (Now we understand why
      ! C++ has templates)

      !real*8, dimension(:,:), allocatable :: ra,rb
      ! wwvv changed this into
      integer, dimension(:,:), allocatable :: ra,rb

      integer                             :: m,n

      m = size(b,1)
      n = size(b,2)

      allocate(rb(m,n))

      if (xomaster) then
         m = size(a,1)
         n = size(a,2)
         allocate(ra(m,n))
      else
         allocate(ra(1,1))
      endif

      rb = b

      call matrix_coll(ra,rb,s%is,s%lm,s%js,s%ln, &
      s%isleft,s%isright,s%istop,s%isbot, &
      xmpi_omaster,xmpi_ocomm)

      if (xomaster) then
         a = ra
      endif

      deallocate(ra,rb)

   end subroutine space_collect_matrix_integer

#endif


#ifdef USEMPI
   !
   !  collects data from processes in master
   !  using the index number of the variable to be
   !  collected
   !
   subroutine space_collect_index(sg,sl,par,index)
      use xmpi_module
      use mnemmodule
      use logging_module
      use params
      type(spacepars)                 :: sg
      type(spacepars), intent(in)     :: sl
      type(parameters)                :: par
      integer, intent(in)             :: index

      integer                         :: lid,eid, wid,ierr
      type(arraytype)                 :: tg,tl
      logical                         :: iscollected

#ifdef USEMPE
      call MPE_Log_event(event_coll_start,0,'cstart')
#endif

      if (xomaster) then
         iscollected = sg%collected(index)
      endif

      call MPI_Bcast(iscollected,1,MPI_LOGICAL,xmpi_omaster,xmpi_ocomm,ierr)

      ! return if this alreay has been collected
      if(iscollected) then
         return
      endif

#ifdef USEMPI
      if (xomaster) then
         call index_allocate(sg,par,index,'r')
      endif
#endif

      call indextos(sl,index,tl)
      call indextos(sg,index,tg)

      select case(tl%type)
       case('i')
         select case(tl%rank)
          case(0)             ! nothing to do
          case(2)
            call space_collect(sl, tg%i2, tl%i2)
          case(3)
            call space_collect(sl, tg%i3, tl%i3)
          case default     ! case 1 and 4 are not handled
            goto 100
         end select   ! rank
       case('r')
         select case(tl%rank)
          case(0)             ! nothing to do
          case(2)
            !if (tl%name .eq. mnem_umean) then
            !  goto 100
            !endif
            call space_collect(sl,tg%r2,tl%r2)
          case(3)
            call space_collect(sl,tg%r3,tl%r3)
          case(4)
            call space_collect(sl,tg%r4,tl%r4)
          case default
         end select   ! rank
       case default
      end select   ! type

#ifdef USEMPE
      call MPE_Log_event(event_coll_start,0,'cend')
#endif

      sg%collected(index) = .true.
      return

100   continue
      call writelog('lse','','Problem in space_collect_index with variable ',trim(tg%name ) )
      call writelog('lse','','Don''t know how to collect that on the masternode')
      call get_logfileid(lid,eid, wid)
      call printvar(tl,lid,eid, wid)
      call halt_program

   end subroutine space_collect_index

   subroutine space_collect_mnem(sg,sl,par,mnem)
      use params
      type(spacepars)                 :: sg
      type(spacepars), intent(in)     :: sl
      type(parameters)                :: par
      character(*)                    :: mnem

      integer index
      index = chartoindex(mnem)
      call space_collect_index(sg,sl,par,index)

   end subroutine space_collect_mnem
#endif

   subroutine gridprops (s)

      IMPLICIT NONE

      type(spacepars),target                  :: s
      !  Temporary local arrays and variables
      integer                                 :: i, j
      real*8,dimension(:,:),allocatable       :: xc      ! x-coordinate c-points
      real*8,dimension(:,:),allocatable       :: yc      ! y-coordinate c-points
      real*8                                  :: dsdnu   ! surface of cell centered around u-point
      real*8                                  :: dsdnv   ! surface of cell centered around v-point
      real*8                                  :: dsdnz   ! surface of cell centered around z-point
      real*8                                  :: x1,y1,x2,y2,x3,y3,x4,y4

      allocate (xc(s%nx+1,s%ny+1))
      allocate (yc(s%nx+1,s%ny+1))

      ! x and y were read in grid_bathy.f90 (initialize.f90) and can be either (cartesian) world coordinates or
      ! XBeach coordinates (xori, yori and alfa are nonzero). Here x and y are transformed to cartesian world coordinates.
      ! XBeach performs all computations (after implementation of curvi-lineair option) world coordinate grid.

      ! world coordinates of z-points
      s%xz=s%xori+s%x*cos(s%alfa)-s%y*sin(s%alfa)
      s%yz=s%yori+s%x*sin(s%alfa)+s%y*cos(s%alfa)

      ! world coordinates of u-points
      do j=1,s%ny+1
         do i=1,s%nx
            s%xu(i,j)=.5d0*(s%xz(i,j)+s%xz(i+1,j))
            s%yu(i,j)=.5d0*(s%yz(i,j)+s%yz(i+1,j))
         enddo
         s%xu(s%nx+1,j)=1.5d0*s%xz(s%nx+1,j)-0.5d0*s%xz(s%nx,j)
         s%yu(s%nx+1,j)=1.5d0*s%yz(s%nx+1,j)-0.5d0*s%yz(s%nx,j)
      enddo

      ! world coordinates of v-points
      if (s%ny>0) then
         do i=1,s%nx+1
            do j=1,s%ny
               s%xv(i,j)=.5d0*(s%xz(i,j)+s%xz(i,j+1))
               s%yv(i,j)=.5d0*(s%yz(i,j)+s%yz(i,j+1))
            enddo
            s%xv(i,s%ny+1)=1.5d0*s%xz(i,s%ny+1)-0.5d0*s%xz(i,s%ny)
            s%yv(i,s%ny+1)=1.5d0*s%yz(i,s%ny+1)-0.5d0*s%yz(i,s%ny)
         enddo
      else
         s%xv=s%xz
         s%yv=s%yz
      endif

      ! world coordinates of corner points
      if (s%ny>0) then
         do j=1,s%ny
            do i=1,s%nx
               xc(i,j)=.25d0*(s%xz(i,j)+s%xz(i+1,j)+s%xz(i,j+1)+s%xz(i+1,j+1))
               yc(i,j)=.25d0*(s%yz(i,j)+s%yz(i+1,j)+s%yz(i,j+1)+s%yz(i+1,j+1))
            enddo
            xc(s%nx+1,j)=0.5d0*(s%xu(s%nx+1,j)+s%xu(s%nx+1,j+1))
            yc(s%nx+1,j)=0.5d0*(s%yu(s%nx+1,j)+s%yu(s%nx+1,j+1))
         enddo
         do i=1,s%nx
            xc(i,s%ny+1)=0.5d0*(s%xv(i,s%ny+1)+s%xv(i+1,s%ny+1))
            yc(i,s%ny+1)=0.5d0*(s%yv(i,s%ny+1)+s%yv(i+1,s%ny+1))
         enddo
         xc(s%nx+1,s%ny+1)=1.5d0*s%xu(s%nx+1,s%ny+1)-0.5*s%xu(s%nx,s%ny+1)
         yc(s%nx+1,s%ny+1)=1.5d0*s%yu(s%nx+1,s%ny+1)-0.5*s%yu(s%nx,s%ny+1)
      else
         xc=s%xu
         yc=s%yu
      endif

      ! s%dsu

      do j=1,s%ny+1
         do i=1,s%nx
            s%dsu(i,j)=((s%xz(i+1,j)-s%xz(i,j))**2+(s%yz(i+1,j)-s%yz(i,j))**2)**(0.5d0)
         enddo
         s%dsu(s%nx+1,j)=s%dsu(s%nx,j)
      enddo

      ! s%dsz

      do j=1,s%ny+1
         do i=2,s%nx+1
            s%dsz(i,j)=((s%xu(i,j)-s%xu(i-1,j))**2+(s%yu(i,j)-s%yu(i-1,j))**2)**(0.5d0)
         enddo
         s%dsz(1,j)=s%dsz(2,j)        ! Robert -> was i=s%nx+1 (calculated in loop), so i=1 more likely
      enddo

      ! s%dsv

      if (s%ny>0) then
         do j=1,s%ny+1
            do i=2,s%nx+1
               s%dsv(i,j)=((xc(i,j)-xc(i-1,j))**2+(yc(i,j)-yc(i-1,j))**2)**(0.5d0)
            enddo
            s%dsv(1,j)=s%dsv(2,j)      ! Robert, need to have this for wave advec
            s%dsv(s%nx+1,j)=s%dsv(s%nx,j)
         enddo
      else
         s%dsv=s%dsz
      endif

      ! s%dsc

      if (s%ny>0) then
         do j=1,s%ny+1
            do i=1,s%nx
               s%dsc(i,j)=((s%xv(i+1,j)-s%xv(i,j))**2+(s%yv(i+1,j)-s%yv(i,j))**2)**(0.5d0)
            enddo
            s%dsc(s%nx+1,j)=s%dsc(s%nx,j)
         enddo
      else
         s%dsc=s%dsu
      endif

      ! s%dnu

      if (s%ny>0) then
         do j=2,s%ny+1
            do i=1,s%nx+1
               s%dnu(i,j)=((xc(i,j)-xc(i,j-1))**2+(yc(i,j)-yc(i,j-1))**2)**(0.5d0)
            enddo
         enddo
         s%dnu(:,1)=s%dnu(:,2)
      else
         s%dnu=100.d0
      endif

      ! s%dnz

      if (s%ny>0) then
         do j=2,s%ny+1
            do i=1,s%nx+1
               s%dnz(i,j)=((s%xv(i,j)-s%xv(i,j-1))**2+(s%yv(i,j)-s%yv(i,j-1))**2)**(0.5d0)
            enddo
         enddo
         s%dnz(:,1)=s%dnz(:,2)
      else
         s%dnz=100.d0
      endif

      ! s%dnv

      if (s%ny>0) then
         do j=1,s%ny
            do i=1,s%nx+1
               s%dnv(i,j)=((s%xz(i,j+1)-s%xz(i,j))**2+(s%yz(i,j+1)-s%yz(i,j))**2)**(0.5d0)
            enddo
         enddo
         s%dnv(:,s%ny+1)=s%dnv(:,s%ny)
      else
         s%dnv=100.d0
      endif

      ! s%dnc

      if (s%ny>0) then
         do j=1,s%ny
            do i=1,s%nx+1
               s%dnc(i,j)=((s%xu(i,j+1)-s%xu(i,j))**2+(s%yu(i,j+1)-s%yu(i,j))**2)**(0.5d0)
            enddo
         enddo
         s%dnc(:,s%ny+1)=s%dnc(:,s%ny)
      else
         s%dnc=100.d0
      endif


      if (s%ny>0) then

         ! dsdnu

         do j=2,s%ny+1
            do i=1,s%nx
               x1=s%xv(i  ,j  ) - s%xv(i  ,j-1)
               x3=s%xv(i+1,j-1) - s%xv(i  ,j-1)
               x2=s%xv(i+1,j  ) - s%xv(i+1,j-1)
               x4=s%xv(i+1,j  ) - s%xv(i  ,j  )
               y1=s%yv(i  ,j  ) - s%yv(i  ,j-1)
               y3=s%yv(i+1,j-1) - s%yv(i  ,j-1)
               y2=s%yv(i+1,j  ) - s%yv(i+1,j-1)
               y4=s%yv(i+1,j  ) - s%yv(i  ,j  )
               dsdnu=0.5d0*(abs(x1*y3-x3*y1)+abs(x2*y4-x4*y2))
               s%dsdnui(i,j)=1.d0/dsdnu
            enddo
         enddo
         s%dsdnui(:,1)=s%dsdnui(:,2)
         s%dsdnui(s%nx+1,:)=s%dsdnui(s%nx,:)

         ! dsdnv

         do j=1,s%ny
            do i=2,s%nx+1
               x1=s%xu(i-1,j+1) - s%xu(i-1,j  )
               x3=s%xu(i  ,j  ) - s%xu(i-1,j  )
               x2=s%xu(i  ,j+1) - s%xu(i  ,j  )
               x4=s%xu(i  ,j+1) - s%xu(i-1,j+1)
               y1=s%yu(i-1,j+1) - s%yu(i-1,j  )
               y3=s%yu(i  ,j  ) - s%yu(i-1,j  )
               y2=s%yu(i  ,j+1) - s%yu(i  ,j  )
               y4=s%yu(i  ,j+1) - s%yu(i-1,j+1)
               dsdnv=0.5d0*(abs(x1*y3-x3*y1)+abs(x2*y4-x4*y2))
               s%dsdnvi(i,j)=1.d0/dsdnv
            enddo
         enddo
         s%dsdnvi(:,s%ny+1)=s%dsdnvi(:,s%ny)
         s%dsdnvi(1,:)=s%dsdnvi(2,:)

         ! dsdnz

         do j=2,s%ny+1
            do i=2,s%nx+1
               x1=xc(i-1,j  ) - xc(i-1,j-1)
               x3=xc(i  ,j-1) - xc(i-1,j-1)
               x2=xc(i  ,j  ) - xc(i  ,j-1)
               x4=xc(i  ,j  ) - xc(i-1,j  )
               y1=yc(i-1,j  ) - yc(i-1,j-1)
               y3=yc(i  ,j-1) - yc(i-1,j-1)
               y2=yc(i  ,j  ) - yc(i  ,j-1)
               y4=yc(i  ,j  ) - yc(i-1,j  )
               dsdnz=0.5d0*(abs(x1*y3-x3*y1)+abs(x2*y4-x4*y2))
               s%dsdnzi(i,j)=1.d0/dsdnz
            enddo
         enddo
         s%dsdnzi(:,1)=s%dsdnzi(:,2)
         s%dsdnzi(1,:)=s%dsdnzi(2,:)

      else

         s%dsdnui=1.d0/(s%dsu*s%dnu)
         s%dsdnvi=1.d0/(s%dsv*s%dnv)
         s%dsdnzi=1.d0/(s%dsz*s%dnz)

      endif

      ! s%alfaz, grid orientation in z-points

      do j=1,s%ny+1
         do i=2,s%nx
            s%alfaz(i,j)=atan2(s%yz(i+1,j)-s%yz(i-1,j),s%xz(i+1,j)-s%xz(i-1,j))
         enddo
         s%alfaz(1,j)=s%alfaz(2,j)
         s%alfaz(s%nx+1,j)=s%alfaz(s%nx,j)
      enddo

      ! s%alfau, grid orientation in u-points

      do j=1,s%ny+1
         do i=1,s%nx
            s%alfau(i,j)=atan2(s%yz(i+1,j)-s%yz(i,j),s%xz(i+1,j)-s%xz(i,j))
         enddo
         s%alfau(s%nx+1,j)=s%alfau(s%nx,j)
      enddo

      ! s%alfav, grid orientation in v-points

      if (s%ny>0) then
         do i=1,s%nx+1
            do j=1,s%ny
               s%alfav(i,j)=atan2(s%yz(i,j+1)-s%yz(i,j),s%xz(i,j+1)-s%xz(i,j))
            enddo
            s%alfav(i,s%ny+1)=s%alfav(i,s%ny)
         enddo
      else
         s%alfav=s%alfaz
      endif

      do j=1,s%ny+1
         s%sdist(1,j)=0
         do i=2,s%nx+1
            s%sdist(i,j)=s%sdist(i-1,j)+s%dsu(i-1,j)
         enddo
      enddo

      do i=1,s%nx+1
         s%ndist(i,1)=0
         do j=2,s%ny+1
            s%ndist(i,j)=s%ndist(i,j-1)+s%dnv(i,j-1)
         enddo
      enddo

      deallocate (xc)
      deallocate (yc)

   end subroutine gridprops

   subroutine ranges_init(s)
      use xmpi_module
      implicit none
      type(spacepars)    :: s

      imin_ee = 2
      imax_ee = s%nx
      jmin_ee = 1
      jmax_ee = s%ny+1

      imin_uu = 2
      imax_uu = s%nx-1
      if (s%ny>0) then
         jmin_uu = 2
         jmax_uu = s%ny
      else
         jmin_uu = 1
         jmax_uu = 1
      endif

      imin_vv = 2
      imax_vv = s%nx
      if (s%ny>0) then
         jmin_vv = 2
         jmax_vv = s%ny-1
      else
         jmin_vv = 1
         jmax_vv = 1
      endif

      imin_zs = 2
      imax_zs = s%nx
      if (s%ny>0) then
         jmin_zs = 2
         jmax_zs = s%ny
      else
         jmin_zs = 1
         jmax_zs = 1
      endif
#ifdef USEMPI

      if (.not. xmpi_istop) then
         imin_ee = 3
         imin_uu = 2
         imin_vv = 3
         imin_zs = 3
      endif

      if (.not. xmpi_isbot) then
         imax_ee = s%nx-1
         imax_uu = s%nx-2
         imax_vv = s%nx-1
         imax_zs = s%nx-1
      endif

      if (.not. xmpi_isleft) then
         jmin_ee = 3
         jmin_uu = 3
         jmin_vv = 2
         jmin_zs = 3
      endif

      if (.not. xmpi_isright) then
         jmax_ee = s%ny-1
         jmax_uu = s%ny-1
         jmax_vv = s%ny-2
         jmax_zs = s%ny-1
      endif

#endif
   end subroutine ranges_init

end module spaceparams
