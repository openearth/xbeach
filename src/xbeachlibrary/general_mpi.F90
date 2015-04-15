module general_mpi_module
   implicit none
   save

   !#define COLLECT_TIMER
   !#define DISTRIBUTE_TIMER

#ifdef USEMPI

  !
  !  NOTE: the matrices are considered to be bordered by a columns
  !   left and right and rows up and down. So, n general, 
  !   only b(2:mb-1,2:nb-1) is send to the master process. 
  !   HOWEVER: if matrix b is 
  !   on the left: b(:,1) is also sent 
  !   on the right: b(:,nb) is also sent
  !   on the top:   b(1,:)  is also sent
  !   on the bottom: b(mb,:) is also sent
  ! The master process receives a(1:ma,1:na). 
  !   dimension of b (1:mb,1:nb). 
  !   dimension of a (1:ma,1:na)
  !
  ! 

  interface matrix_distr
      module procedure matrix_distr_real8
      module procedure matrix_distr_integer
  end interface matrix_distr

  interface matrix_coll
      module procedure matrix_coll_real8
      module procedure matrix_coll_integer
  end interface matrix_coll

  interface shift_borders
     module procedure shift_borders_matrix_real8
  end interface shift_borders

   interface vector_distr
      module procedure vector_distr_real8
      module procedure block_vector_distr_real8
   end interface vector_distr

contains
   !_____________________________________________________________________________
   subroutine vector_distr_real8(a,b,is,lm,root,comm)
    !
      ! Distribute vector a on non-root processes
    ! Process P receives a(is(P):is(P)+lm(P)-1)
    !
    ! Parameters:
    ! real*8 a (in): vector to distribute, only used at root
    ! real*8 b (out): local vector
    ! integer is(number-of-processes) (in): see above: needed on all processes
    ! integer lm(number-of-processes) (in): see above: needed on all processes
    ! integer root (in): process that holds a, is, lm
    ! integer comm (in): MPI communicator
    !
      ! complication: a and b are not necessarily contiguous,
      !  but calling mpi_alltoallw with a and b, the compiler
      !  will generate contiguous arrays.
      !
      ! root process must have rank .eq. 0
    !
    use mpi
    implicit none
    real*8, dimension(:), intent(in)  :: a
    real*8, dimension(:), intent(out) :: b
    integer, dimension(:), intent(in) :: is
    integer, dimension(:), intent(in) :: lm
    integer, intent(in)               :: root
    integer, intent(in)               :: comm

      integer                          :: ier,ra,sz,basic_type,i
      integer,allocatable,dimension(:) :: recvtypes,sendtypes,recvcounts,sendcounts,sdispls,rdispls
      integer,dimension(1)             :: sizes,subsizes,starts

      basic_type = MPI_DOUBLE_PRECISION
      call MPI_Comm_size(comm, sz, ier)
      call MPI_Comm_rank(comm, ra, ier)

      if (root /= 0) then
         print *,'Error in vector_distr_real8: root must be 0, but is:',root
         call MPI_Abort(MPI_COMM_WORLD,1,ier)
          endif

      allocate(recvtypes(sz))
      allocate(sendtypes(sz))
      allocate(recvcounts(sz))
      allocate(sendcounts(sz))
      allocate(sdispls(sz))
      allocate(rdispls(sz))

      sdispls    = 0
      rdispls    = 0
      recvtypes  = MPI_CHARACTER
      sendtypes  = MPI_CHARACTER
      recvcounts = 0
      sendcounts = 0

      !
      ! Create MPI types
      !

      ! MPI_TYPE_CREATE_SUBARRAY(NDIMS, ARRAY_OF_SIZES, ARRAY_OF_SUBSIZES,
      !   ARRAY_OF_STARTS, ORDER, OLDTYPE, NEWTYPE, IERROR)
      !   INTEGER    NDIMS, ARRAY_OF_SIZES(*), ARRAY_OF_SUBSIZES(*),
      !   ARRAY_OF_STARTS(*), ORDER, OLDTYPE, NEWTYPE, IERROR

      ! determine mpi_types for the receive matrices
      ! all processes will receive only from root
      sizes     = shape(b)
      subsizes  = lm(ra+1)
      starts    = 0
      call MPI_Type_create_subarray(1,sizes,subsizes,starts, &
      MPI_ORDER_FORTRAN, basic_type,recvtypes(1), ier)
      call MPI_Type_commit(recvtypes(1),ier)
      recvcounts(1) = 1


      ! determine mpi types for the senders
      if(ra == root) then
         ! root sends to everyone

         do i=1,sz
            sizes      = shape(a)
            subsizes   = lm(i)
            starts     = is(i) - 1
            sendcounts = 1
            call MPI_Type_create_subarray(1,sizes,subsizes,starts,  &
            MPI_ORDER_FORTRAN, basic_type,sendtypes(i),ier)
            call MPI_Type_commit(sendtypes(i),ier)
       enddo
    endif

      call MPI_Alltoallw(a,sendcounts,sdispls,sendtypes, &
      b,recvcounts,rdispls,recvtypes,comm,ier)

      do i=1,sz
         if (sendtypes(i) /= MPI_CHARACTER) then
            call MPI_Type_free(sendtypes(i),ier)
         endif
         if (recvtypes(i) /= MPI_CHARACTER) then
            call MPI_Type_free(recvtypes(i),ier)
         endif
      enddo

   end subroutine vector_distr_real8

   !_____________________________________________________________________________

   subroutine block_vector_distr_real8(a,b,is,lm,root,comm)
      !
      ! Distribute matrix a on non-root processes
      ! Process P receives a(is(P):is(P)+lm(P)-1,:)
      !
      ! Parameters:
      ! real*8 a(:,:) (in): matrix to distribute, only used at root
      ! real*8 b(:,:) (out): local matrix
      ! integer is(number-of-processes) (in): see above: needed on all processes
      ! integer lm(number-of-processes) (in): see above: needed on all processes
      ! integer root (in): process that holds a, is, lm
      ! integer comm (in): MPI communicator
      !
      ! complication: a and b are not necessarily contiguous,
      !  but calling mpi_alltoallw with a and b, the compiler
      !  will generate contiguous arrays.
      !
      ! root process must have rank .eq. 0
      !
    use mpi
    implicit none
    real*8,  dimension(:,:), intent(in)  :: a
    real*8,  dimension(:,:), intent(out) :: b
      integer, dimension(:), intent(in) :: is
      integer, dimension(:), intent(in) :: lm
      integer, intent(in)               :: root
      integer, intent(in)               :: comm

      integer                          :: ier,ra,sz,basic_type,i
      integer,allocatable,dimension(:) :: recvtypes,sendtypes,recvcounts,sendcounts,sdispls,rdispls
      integer,dimension(2)             :: sizes,subsizes,starts

      basic_type = MPI_DOUBLE_PRECISION
      call MPI_Comm_size(comm, sz, ier)
      call MPI_Comm_rank(comm, ra, ier)

      if (root /= 0) then
         print *,'Error in block_vector_distr_real8: root must be 0, but is:',root
         call MPI_Abort(MPI_COMM_WORLD,1,ier)
      endif

      allocate(recvtypes(sz))
      allocate(sendtypes(sz))
      allocate(recvcounts(sz))
      allocate(sendcounts(sz))
      allocate(sdispls(sz))
      allocate(rdispls(sz))

      sdispls    = 0
      rdispls    = 0
      recvtypes  = MPI_CHARACTER
      sendtypes  = MPI_CHARACTER
      recvcounts = 0
      sendcounts = 0

      !
      ! Create MPI types
      !

      ! MPI_TYPE_CREATE_SUBARRAY(NDIMS, ARRAY_OF_SIZES, ARRAY_OF_SUBSIZES,
      !   ARRAY_OF_STARTS, ORDER, OLDTYPE, NEWTYPE, IERROR)
      !   INTEGER    NDIMS, ARRAY_OF_SIZES(*), ARRAY_OF_SUBSIZES(*),
      !   ARRAY_OF_STARTS(*), ORDER, OLDTYPE, NEWTYPE, IERROR

      ! determine mpi_types for the receive matrices

      ! all processes will receive only from root
      sizes     = shape(b)
      subsizes  = (/ lm(ra+1), size(b,2) /)
      starts    = 0
      call MPI_Type_create_subarray(2,sizes,subsizes,starts, &
      MPI_ORDER_FORTRAN, basic_type,recvtypes(1), ier)
      call MPI_Type_commit(recvtypes(1),ier)
      recvcounts(1) = 1

      ! determine mpi types for the senders
      if (ra == root) then
         ! root sends to everyone

         do i=1,sz
            sizes      = shape(a)
            subsizes   = (/ lm(i),     size(a,2) /)
            starts     = (/ is(i) - 1, 0         /)
            sendcounts = 1
            call MPI_Type_create_subarray(2,sizes,subsizes,starts,  &
            MPI_ORDER_FORTRAN, basic_type,sendtypes(i),ier)
            call MPI_Type_commit(sendtypes(i),ier)
        enddo
    endif

      call MPI_Alltoallw(a,sendcounts,sdispls,sendtypes, &
      b,recvcounts,rdispls,recvtypes,comm,ier)

      do i=1,sz
         if (sendtypes(i) /= MPI_CHARACTER) then
            call MPI_Type_free(sendtypes(i),ier)
         endif
         if (recvtypes(i) /= MPI_CHARACTER) then
            call MPI_Type_free(recvtypes(i),ier)
         endif
      enddo

   end subroutine block_vector_distr_real8

   !___________________________________________________________________________

   subroutine matrix_distr_real8(a,b,is,lm,js,ln,  &
   root,comm)
      !
      ! This subroutine distributes matrix on root to the non-root
      ! processes.
      ! Root will not send to itself.
      ! Root MUST be process #0 in comm
      !
      ! a: matrix to be distributed
      ! b: receive matrix, process p
      !    (1:lm(p),1:ln(p)) will be filled (p=1..size of comm)
      !    Normally: the whole of b
      ! is,js: on process p, b(1,1) coincides with a(is(p),js(p)),
      !        p=1..size of comm
      ! root: root process, must be zero
      ! comm: communicator
      !
      ! Note: root process does not send or copy anything to itself,
      !       'there is no b on root'
      !

      use mpi
      implicit none
      real*8, intent(in),  dimension(:,:)  :: a
      real*8, intent(out), dimension(:,:)  :: b
      integer, intent(in), dimension(:)    :: is, lm
      integer, intent(in), dimension(:)    :: js, ln
      integer, intent(in)                  :: root, comm

      integer              :: sz,ier,ra,i
      integer, allocatable :: recvtypes(:), sendtypes(:), recvcounts(:),sendcounts(:)
      integer, allocatable :: sdispls(:), rdispls(:)
      integer sizes(2),subsizes(2),starts(2)
      integer basic_type

      basic_type = MPI_DOUBLE_PRECISION

      include 'genmpi_distr.inc'

   end subroutine matrix_distr_real8

   !___________________________________________________________________________

   subroutine matrix_distr_integer(a,b,is,lm,js,ln,  &
   root,comm)
      !
      ! This subroutine distributes matrix on root to the non-root
      ! processes.
      ! Root will not send to itself.
      ! Root MUST be process #0 in comm
      !
      ! a: matrix to be distributed
      ! b: receive matrix, process p
      !    (1:lm(p),1:ln(p)) will be filled (p=1..size of comm)
      !    Normally: the whole of b
      ! is,js: on process p, b(1,1) coincides with a(is(p),js(p)),
      !        p=1..size of comm
      ! root: root process, must be zero
      ! comm: communicator
      !
      ! Note: root process does not send or copy anything to itself,
      !       'there is no b on root'
      !

      use mpi
      implicit none
      integer, intent(in),  dimension(:,:)  :: a
      integer, intent(out), dimension(:,:)  :: b
      integer, intent(in), dimension(:)    :: is, lm
      integer, intent(in), dimension(:)    :: js, ln
      integer, intent(in)                  :: root, comm

      integer              :: sz,ier,ra,i
      integer, allocatable :: recvtypes(:), sendtypes(:), recvcounts(:),sendcounts(:)
      integer, allocatable :: sdispls(:), rdispls(:)
      integer sizes(2),subsizes(2),starts(2)
      integer basic_type

      basic_type = MPI_INTEGER

      include 'genmpi_distr.inc'

   end subroutine matrix_distr_integer

   !___________________________________________________________________________

   subroutine matrix_coll_real8(a,b,is,lm,js,ln,&
   isleft,isright,istop,isbot,root,comm)
      ! collect matrix on root process
      !  NOTE: the collect is done from the non-root processes only
    !
    ! parameters:
      ! real*8 a(:,:) (out): matrix to be collected.
    !                      This matrix has to be available, it is not
    !                      allocated by this subroutine.
      ! real*8 b(:,:) (in): local scattered matrix
      ! integer is,lm,js,ln(number of non-root processes) (in): description of
    !   the location of the submatrices b in matrix a:
    !   Element b(1,1) coincides with a(is(rank+1),js(rank+1)), where rank is
    !   the rank of the process (0..procs-1)
    !   The length of the 1st dimension of b is given by lm(rank+1)
    !   The length of the 2nd dimension of b is given by ln(rank+1)
      ! logical isleft, isright, istop, isbot (in): is this submatrix at the
      !  left, right top or bottom of the general matrix.
    ! integer root (in): the rank of the process where a is available
    ! integer comm (in): the MPI communicator to be used.
    !
    ! Note: a,is,js are used at the root process only
    !    
      ! Note: the root process must have rank .eq. 0
    !
    use mpi
    implicit none
      real*8, intent(out),  dimension(:,:)  :: a
      real*8, intent(in), dimension(:,:)  :: b
    integer, intent(in), dimension(:)    :: is,lm
    integer, intent(in), dimension(:)    :: js,ln
      logical, intent(in), dimension(:)    :: isleft, isright, istop, isbot
    integer, intent(in)                  :: root,comm

      integer sz,ier,ra,i
      integer, allocatable :: recvtypes(:), sendtypes(:), recvcounts(:),sendcounts(:)
      integer, allocatable :: sdispls(:), rdispls(:)
      integer sizes(2),subsizes(2),starts(2)
      integer basic_type
      ! nover is number of overlapping rows/columns of submatrices.
      ! nbord is number of borders.

      integer, parameter                   :: nbord = 2
      integer, parameter                   :: nover = 2*nbord

      basic_type = MPI_DOUBLE_PRECISION

      include 'genmpi_coll.inc'

   end subroutine matrix_coll_real8
   !___________________________________________________________________________

   subroutine matrix_coll_integer(a,b,is,lm,js,ln,&
   isleft,isright,istop,isbot,root,comm)
      ! for description, see matrix-coll_real8
      use mpi
      implicit none
      integer, intent(out),  dimension(:,:) :: a
      integer, intent(in), dimension(:,:)   :: b
      integer, intent(in), dimension(:)     :: is, lm
      integer, intent(in), dimension(:)     :: js, ln
      logical, intent(in), dimension(:)     :: isleft, isright, istop, isbot
      integer, intent(in)                   :: root, comm

      integer sz,ier,ra,i
      integer, allocatable :: recvtypes(:), sendtypes(:), recvcounts(:),sendcounts(:)
      integer, allocatable :: sdispls(:), rdispls(:)
      integer sizes(2),subsizes(2),starts(2)
      integer basic_type
      ! nover is number of overlapping rows/columns of submatrices.
      ! nbord is number of borders.

      integer, parameter                   :: nbord = 2
      integer, parameter                   :: nover = 2*nbord

      basic_type = MPI_INTEGER

      include 'genmpi_coll.inc'

   end subroutine matrix_coll_integer
   !___________________________________________________________________________



  subroutine matrix_distr_scatter_integer(a,b,is,lm,js,ln,root,comm)
    ! distribute matrix on processes
    !
    ! parameters:
    ! integer a(:,:) (in): matrix to be scattered. 
    ! integer b(:,:) (out): local matrix to scatter to.
    !                      This matrix has to be available, it is not
    !                      allocated by this subroutine.
    ! integer is,lm,js,ln(number of processes) (in): description of
    !   the location of the submatrices b in matrix a:
    !   Element b(1,1) coincides with a(is(rank+1),js(rank+1)), where rank is
    !   the rank of the process (0..procs-1)
    !   The length of the 1st dimension of b is given by lm(rank+1)
    !   The length of the 2nd dimension of b is given by ln(rank+1)
    ! integer root (in): the rank of the process where a is available
    ! integer comm (in): the MPI communicator to be used.
    !
    ! Note: a,is,js are used at the root process only
    !    
    ! Method
    !
    !  Distribute in a number of scatterv's the colums of the submatrices.
    !  Per scatterv, distribute (number of processes) columns, so
    !  we need maxval(ln) scatters.
    !
    use mpi
    implicit none
    integer, intent(in),  dimension(:,:)  :: a
    integer, intent(out), dimension(:,:)  :: b
    integer, intent(in), dimension(:)    :: is,lm
    integer, intent(in), dimension(:)    :: js,ln
    integer, intent(in)                  :: root,comm

    integer                              :: j,p,rank,ierror,procs,jj,ma,jmax
    integer                              :: nlocal,mlocal,rcnt
    integer, allocatable, dimension(:)   :: displs, cnts
    integer, allocatable, dimension(:,:) :: aa

    integer, parameter                   :: tag = 1

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    nlocal = ln(rank+1)
    mlocal = lm(rank+1)

    if (rank .eq. root) then
       ma = size(a,1)
       allocate(displs(procs), cnts(procs))
       allocate(aa(size(a,1),size(a,2)))
       aa = a
    else
       ma = 1                       ! to quiet the compiler
       allocate(displs(1), cnts(1)) ! these are not needed in non-root process, but
       ! it is safer to have valid addresses
       allocate(aa(1,1))
    endif

    jmax = maxval(ln(1:procs))     ! number of scatters needed

    do j=1,jmax ! scatter all j'th columns 
       if (rank .eq. root) then
          do p=1,procs
             if (ln(p) .ge. j) then ! compute displacement of the start of the j-th 
                ! column of submatrix(p)
                displs(p) = (js(p) -1 + j -1)*ma + is(p) - 1
                cnts(p) = lm(p)      ! number of elements in this column
             else                   ! for submatrix(p), we are out of columns
                displs(p) = 0        ! put in some valid value for the displacement
                cnts(p) = 0          ! and use zero for the number of elements in
                ! this column
             endif
          enddo
       endif
       if (nlocal .ge. j) then    ! this is for all processes
          jj=j                     ! jj equals the column number of the local
          ! matrix
          rcnt = mlocal            ! number of elements in this column
       else
          jj=1                     ! we are out of columns, give jj a valid number
          rcnt = 0                 ! and set the number of elements to receive to zero
       endif

       call MPI_Scatterv(aa(1,1),  cnts(1), displs(1), MPI_INTEGER,  &
            b(:,jj), rcnt,         MPI_INTEGER,  &
            root,    comm, ierror)
    enddo

    deallocate(displs,cnts,aa)

  end subroutine matrix_distr_scatter_integer


  subroutine decomp(n, numprocs, myid, s, e)

    !  From the mpich-distribution: MPE_Decomp1d
    !
    !  n        : integer, input (for example: number of tasks or 
    !                               dimension of an array)
    !  numprocs : integer, input (number of processes)
    !  myid     : integer, input (MPI id of this process)
    !  s        : integer, output
    !  e        : integer, output
    !  
    !  The purpose of this subroutine is to get an as equal as possible
    !  division of the integers 1..n among numprocs intervals.
    !  A common use is to divide work among MPI processes.
    !  
    !  example:
    !  
    !  integer numprocs, ier, myid, n, s, e, i
    !  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ier)
    !  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ier)
    !  n = 1000
    !  call decomp(n, nprocs, myid, s, e)
    !  do i = s,e
    !    ... 
    !  enddo
    !
    integer, intent(in) :: n, numprocs, myid
    integer, intent(out):: s, e
    integer nlocal
    integer deficit

    nlocal  = n / numprocs
    s       = myid * nlocal + 1
    deficit = mod(n,numprocs)
    s       = s + min(myid,deficit)
    if (myid .lt. deficit) then
       nlocal = nlocal + 1
    endif
    e = s + nlocal - 1
    if (e .gt. n .or. myid .eq. numprocs-1) then
       e = n
    endif

  end subroutine decomp

  subroutine det_submatrices(ma,na,mp,np,is,lm,js,ln,isleft,isright,istop,isbot)
    !
    ! determine a division of a ma*na matrix on mp*np processes
    ! ma: integer(in): 1st dimension of matrix to divide
    ! na: integer(in): 2nd dimension of this matrix
    ! mp: integer(in): 1st dimension of processorgrid
    ! np: integer(in): 2nd dimension of processorgrid
    ! is(mp*np): integer(out): is(i) will be equal to the
    !                          start row of the i-th submatrix
    ! lm(mp*np): integer(out): lm(i) will be equal to the
    !                          first dimension of the i-th submatrix
    ! js(mp*np): integer(out): js(i) will be equal to the
    !                          start column of the i-th submatrix
    ! ln(mp*np): integer(out): ln(i) will be equal to the
    !                          2nd dimension of the i-th submatrix
    ! isleft(mp*np): logical(out): isleft(i) is true if the i-th submatrix 
    !                              shares the first column with the global matrix
    !                              a(:,1)
    ! isright(mp*np): logical(out): isright(i) is true if the i-th submatrix 
    !                              shares the last column with the global matrix
    !                              a(ma,:)
    ! istop(mp*np): logical(out): istop(i) is true if the i-th submatrix 
    !                              shares the first row with the global matrix
    !                              a(1,:)
    ! isbot(mp*np): logical(out): isbot(i) is true if the i-th submatrix 
    !                              shares the last row with the global matrix
    !                              a(na,:)
    ! 
    ! The submatrices overlap:
    ! 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
    ! 1 2 3 4 5 6 7 8 9           1 2 3 4 5 6 7 8 9 0 1
    !           1 2 3 4 5 6 7 8 9 0 1 2 3 
    !
    ! + + + + + x x x x + + + + + x x x x + + + + + + +
    ! + + + + + x x x x + + + + + x x x x + + + + + + +
    ! + + + + + x x x x + + + + + x x x x + + + + + + +
    ! + + + + + x x x x + + + + + x x x x + + + + + + +
    ! + + + + + x x x x + + + + + x x x x + + + + + + +
    ! + + + + + x x x x + + + + + x x x x + + + + + + +
    ! + + + + + x x x x + + + + + x x x x + + + + + + +
    !
    ! This picture is for the overlap in horizontal direction,
    ! The overlap is 4.
    ! + -> matrix element
    ! x -> overlapping matrix element
    ! the order of the submatrices is 9,13,11.
    ! The order of the global matrix is 25
    !  
    implicit none
    integer, intent(in)                :: ma,na,mp,np
    integer, dimension(:), intent(out) :: is,lm,js,ln
    logical, dimension(:), intent(out) :: isleft,isright,istop,isbot

    integer i,j,s,e,k
    integer, parameter :: nover = 4

    isleft = .false.
    isright = .false.
    istop = .false.
    isbot = .false.

    do i = 1, mp
       call decomp(ma-nover, mp, i - 1, s, e)
       k = 0
       do j = 1, np
          if (j .eq. 1) then
             isleft(i+k) = .true.
          endif
          if (j .eq. np) then
             isright(i+k) = .true.
          endif
          is(i + k) = s
          lm(i + k) = e - s + nover + 1
          k = k + mp
       enddo
    enddo
    k = 0
    do j=1, np
       call decomp(na-nover, np, j - 1, s, e)
       do i = 1, mp
          if (i .eq. 1) then
             istop(i+k) = .true.
          endif
          if (i .eq. mp) then
             isbot(i+k) = .true.
          endif
          js(i + k) = s
          ln(i + k) = e - s + nover + 1
       enddo
       k = k + mp
    enddo

  end subroutine det_submatrices

  subroutine shift_borders_matrix_real8(a,left,right,top,bot,comm)
    !
    ! shift borders to and from neighbours.
    ! 
    use mpi
    implicit none
    real*8, dimension(:,:), intent(inout) :: a
    integer, intent(in)                   :: left,right,top,bot,comm

    integer, parameter                    :: datatype = MPI_DOUBLE_PRECISION
    integer ierror,ma,na

    ma = ubound(a,1)
    na = ubound(a,2)

    ! send to left, receive from right
    call MPI_Sendrecv(a(2:ma-1,2), ma-2, datatype,       &
         left, 1,                           &
         a(2:ma-1,na), ma-2, datatype,      &
         right, 1,                          &
         comm, MPI_STATUS_IGNORE, ierror)
    ! send to right, receive from left
    call MPI_Sendrecv(a(2:ma-1,na-1), ma-2, datatype,    &
         right, 2,                          &
         a(2:ma-1,1), ma-2, datatype,       &
         left, 2,                           &
         comm, MPI_STATUS_IGNORE, ierror)

    ! send to bot, receive from top
    call MPI_Sendrecv(a(ma-1,:), na, datatype,           &
         bot, 3,                            &
         a(1,:), na, datatype,              &
         top, 3,                            &
         comm, MPI_STATUS_IGNORE, ierror)
    ! send to top, receive from bot
    call MPI_Sendrecv(a(2,:), na, datatype,              &
         top, 4,                            &
         a(ma,:), na, datatype,             &
         bot, 4,                            &
         comm, MPI_STATUS_IGNORE, ierror)

  end subroutine shift_borders_matrix_real8

#endif
end module general_mpi_module
