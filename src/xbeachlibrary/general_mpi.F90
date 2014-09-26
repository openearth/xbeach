module general_mpi_module
#ifdef USEMPI

  ! This module contains three variants of subroutines to scatter
  ! a matrix from the root process to the other processes.
  ! matrix_distr_scatter: uses MPI_Scatterv, scattering groups
  !                       of columns 
  ! matrix_distr_send:    uses MPI_Isend and MPI_Irecv, scattering
  !                       groups od columns
  ! matrix_distr_sendmat: uses MPI_Send and MPI_Recv, scattering
  !                       submatrices as a whole
  !                       
  ! and the opposite operation:
  ! matrix_coll_gather:   uses MPI_Gatherv, gathering groups
  !                       of columns
  ! matrix_coll_recv:     uses MPI_Isend and MPI_Recv, gathering groups
  !                       of columns
  ! matrix_coll_matrecv:  uses MPI_Send and MPI_Recv, gathering 
  !                       submatrices as a whole
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

  ! The parameter lists are the same.
  ! On the lisa system at sara (Xeon, infiniband), the send variety
  ! is about twice as fast as the scatter variety.
  ! 
  ! Most probable cause is that in the scatter variety, a number
  ! of MPI_Scatterv's are used. They will run after the previous
  ! has been completed, while in the other, as many MPI_isend's and
  ! MPI_Irecv's are executed in parallel, without waiting.
  ! There is no non-waiting MPI_Scatterv available, alas.
  !
  ! On the same system, there is no significant difference in 
  ! performance between matrix_coll_gather and matrix_coll_recv.
  !
  ! The matrecv and matsend varieties are even faster than the
  ! recv and send varieties. Moreover, the code is simpler, so
  ! it is bet to use these, see interfaces matrix_distr and matrix_coll.

  ! Most of the subroutines end with a MPI_Barrier, this was
  ! necessary for usage with openmpi, at least for the
  ! distbribution routines. Didn't understand why, but now
  ! I know: this is a bug in the shared memory interface of
  ! openmpi.
  ! 


  interface matrix_distr
     !  module procedure matrix_distr_send_real8
     !  module procedure matrix_distr_send_integer
     !  module procedure matrix_distr_scatter_real8
     !  module procedure matrix_distr_scatter_integer
     module procedure matrix_distr_sendmat_real8
     module procedure matrix_distr_sendmat_integer
  end interface matrix_distr

  interface matrix_coll
     !  module procedure matrix_coll_gather_real8
     !  module procedure matrix_coll_recv_real8
     module procedure matrix_coll_recvmat_int
     module procedure matrix_coll_recvmat_real8
  end interface matrix_coll

  interface shift_borders
     module procedure shift_borders_matrix_real8
  end interface shift_borders

contains
  subroutine vector_distr_send(a,b,is,lm,root,comm)
    !
    ! Distribute vector a on processes
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
    ! Method: root process sends, using MPI_Isend parts of a to
    !         the other processes. On the root process a copy is
    !         made.
    !
    use mpi
    implicit none
    real*8, dimension(:), intent(in)  :: a
    real*8, dimension(:), intent(out) :: b
    integer, dimension(:), intent(in) :: is
    integer, dimension(:), intent(in) :: lm
    integer, intent(in)               :: root
    integer, intent(in)               :: comm

    integer                             :: ierror,p,procs,rank,pp,mlocal
    integer, dimension(:), allocatable  :: req
    real*8, dimension(:), allocatable   :: aa

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)
    allocate(req(procs-1))

    mlocal = lm(rank+1)

    pp=0
    if ( rank .eq. root ) then
       allocate(aa(size(a)))
       aa = a               ! must have contiguous memory for a because
       ! of mpi_isend
       do p = 1,procs
          if ( p-1 .eq. rank ) then
             b(1:lm(p)) = a(is(p) : is(p)+lm(p)-1)
          else
             pp = pp + 1
             call MPI_Isend(aa(is(p)), lm(p), MPI_DOUBLE_PRECISION, p-1, &
                  11, comm, req(pp), ierror)
          endif
       enddo
       call MPI_Waitall(pp,req,MPI_STATUSES_IGNORE,ierror)
       deallocate(aa)
    else
       call MPI_Recv(b, mlocal, MPI_DOUBLE_PRECISION, root, 11,&
            comm, MPI_STATUS_IGNORE, ierror) 
    endif

    deallocate(req)

  end subroutine vector_distr_send

  subroutine block_vector_distr_y(a,b,is,lm,root,comm)
    !distribute matrix on processes, the first dimension is 
    ! divided among the processes.
    ! parameters:
    ! real*8 a(:,:) (in): matrix to be scattered
    ! real*8 b(:,:) (out): local matrix to scatter to
    !                      This matrix has to be available, it is not
    !                      allocated by this subroutine.
    ! is, lm(number of processes) (in): description
    !   of the location of the submatrices b in matrix a:
    !   b(1,1) coincides with a(is(rank+1),1)
    !   The size of the first dimension of b is given by lm(rank+1)
    !   The size of the second dimension of b is equal to size(a,2)
    ! integer root (in): the rank of the process where a is available
    ! integer comm (in): the MPI communicator to be used
    ! Note: a,is are used at the root process only

    use mpi
    implicit none

    real*8,  dimension(:,:), intent(in)  :: a
    real*8,  dimension(:,:), intent(out) :: b
    integer, dimension(:),   intent(in)  :: is,lm
    integer,                 intent(in)  :: root,comm

    integer                              :: ierror,p,procs,rank,pp
    integer                              :: i,j,l,buflen
    real*8,  dimension(:), allocatable   :: buf
    integer, dimension(:), allocatable   :: counts, displs

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    allocate(counts(procs))
    allocate(displs(procs))

    pp=0
    if ( rank .eq. root ) then
      ! create buffer to hold the submatrices to scatter.
      ! Note that this buffer is larger that the matrix a because
      ! of overlap of the submatrices
      ! also compute the displacemets and counts needed
      ! for mpi_scatterv

      buflen = 0
      do p=1,procs
        buflen = buflen + lm(p)*size(a,2)

        counts(p) = lm(p)*size(a,2)

        if (p .eq. 1) then
          displs(p) = 0
        else
          displs(p) = displs(p-1)+counts(p-1)
        endif
      enddo

      allocate(buf(buflen))

      ! copy the submatrices into the buffer

      l = 0
      do p=1,procs
        do j=1,size(a,2)
          do i=is(p),is(p)+lm(p)-1
            l = l+1
            buf(l) = a(i,j)
          enddo
        enddo
      enddo
    else
      allocate(buf(1))  ! to get a valid address for buf on non-root processes
    endif

    ! scatter the submatrices

    call MPI_Scatterv(buf, counts, displs, MPI_DOUBLE_PRECISION, b,  &
                      size(b), MPI_DOUBLE_PRECISION, root, comm, ierror)

    deallocate(displs)
    deallocate(counts)
    deallocate(buf)

  end subroutine block_vector_distr_y

  subroutine matrix_distr_scatter_real8(a,b,is,lm,js,ln,root,comm)
    ! distribute matrix on processes
    !
    ! parameters:
    ! real*8 a(:,:) (in): matrix to be scattered. 
    ! real*8 b(:,:) (out): local matrix to scatter to.
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
    real*8, intent(in),  dimension(:,:)  :: a
    real*8, intent(out), dimension(:,:)  :: b
    integer, intent(in), dimension(:)    :: is,lm
    integer, intent(in), dimension(:)    :: js,ln
    integer, intent(in)                  :: root,comm

    integer                              :: j,p,rank,ierror,procs,jj,ma,jmax
    integer                              :: nlocal,mlocal,rcnt
    integer, allocatable, dimension(:)   :: displs, cnts
    real*8, allocatable, dimension(:,:)  :: aa

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    nlocal = ln(rank+1)
    mlocal = lm(rank+1)

    if (rank .eq. root) then
       ma = size(a,1)
       allocate(displs(procs), cnts(procs))
       ! because the counts and displacements are only valid for a 
       ! contiguous matrix, we take care of that:

       allocate(aa(size(a,1),size(a,2)))
       aa = a

    else
       ma = 1                       ! to quiet the compiler
       allocate(displs(1), cnts(1)) ! these are not needed in non-root process, but
       allocate(aa(1,1))
       displs(1) = 0
       cnts(1)   = 0
       ! it is safer to have valid addresses
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

       call MPI_Scatterv(aa,  cnts, displs, MPI_DOUBLE_PRECISION,  &
            b(:,jj), rcnt,         MPI_DOUBLE_PRECISION,  &
            root,    comm, ierror)
    enddo

    deallocate(displs,cnts,aa)

  end subroutine matrix_distr_scatter_real8

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

       call MPI_Scatterv(aa,  cnts, displs, MPI_INTEGER,  &
            b(:,jj), rcnt,         MPI_INTEGER,  &
            root,    comm, ierror)
    enddo

    deallocate(displs,cnts,aa)

  end subroutine matrix_distr_scatter_integer

  subroutine matrix_distr_send_real8(a,b,is,lm,js,ln,root,comm)
    ! distribute matrix on processes
    !
    ! parameters:
    ! real*8 a(:,:) (in): matrix to be scatterd. 
    ! real*8 b(:,:) (out): local matrix to scatter to.
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
    ! Sends colums of submatrices to the processes using MPI_Isend. 
    ! This routine tries to send as many columns in parallel
    ! to different processes as possible, see the usage of req(:)
    ! and maxreq below.
    ! The columns are received using MPI_Irecv.
    ! 
    use mpi
    implicit none
    real*8, intent(in),  dimension(:,:)  :: a
    real*8, intent(out), dimension(:,:)  :: b
    integer, intent(in), dimension(:)    :: is,lm
    integer, intent(in), dimension(:)    :: js,ln
    integer, intent(in)                  :: root,comm

    integer                              :: j,rank,ierror,l,procs

    integer                              :: sends, sends_todo, length
    integer                              :: p,mlocal,nlocal

    logical                              :: done

    integer                              :: maxreq
    integer, allocatable, dimension(:)   :: req

    integer, parameter                   :: tag = 1

    integer, parameter                   :: datatype = MPI_DOUBLE_PRECISION
    real*8, allocatable, dimension(:,:)  :: aa,bb

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    nlocal = ln(rank+1)
    mlocal = lm(rank+1)

    ! calculate total number of sends to be done: that is equal to 
    ! the total number of columns in the submatrices

    if (rank .eq. root) then
       sends_todo = sum(ln(1:procs))
       ! maxreq is the maximum number of outstanding Isend or Irecv requests.
       ! Experiments learned that the value of maxreq is not critical,
       ! a value of 10 gave on the lisa system (Xeon, infiniband) 
       ! only somewhat worse results. Probably, the slower the network,
       ! the higher this value can be set. Maximum usable is the total number
       ! of columns (sends_todo, see above)
       maxreq = 1000
    else
       sends_todo = 0    ! to quiet the compiler
       maxreq = nlocal   ! on non-root processes, we need only the number of local
       ! columns
    endif

    ! distribute matrix. 
    ! to ensure that every process is as busy as possible, we do the
    ! loops on the processes in the inner loop.

    allocate(req(maxreq))
    req=MPI_REQUEST_NULL    ! initialize req with a null-value

    allocate(bb(size(b,1),size(b,2)))
    if ( rank .eq. root) then
       ! we need a contiguous matrix a, because
       !  we are going to use isend.
       allocate(aa(size(a,1),size(a,2)))
       aa = a

       sends=0
       l = 0
       j = 0
       loop: do 
          j=j+1                 ! j runs over the columns
          do p=0, procs - 1
             if (ln(p + 1) .ge. j) then ! is j in the range of submatrix(p)
                length = lm(p + 1)
                if (p .eq. rank) then    ! submtarix(p) is on root, just copy the column
                   b(1:length,j) = a(is(rank+1):is(rank+1)+length-1,j+js(rank+1)-1)
                else
                   ! search a free request holder:
                   done = .false.
                   do while( .not. done) 
                      l=l+1
                      if (l .gt. maxreq) l=1
                      if (req(l) .eq. MPI_REQUEST_NULL) then
                         done = .true.
                      else
                         !
                         ! MPI_Test gives .true. if req(l).eq.MPI_REQUEST_NULL,
                         ! but the standard says nothing about this case,
                         ! so we test first if req(l) is equal to MPI_REQUEST_NULL
                         !
                         call MPI_Test(req(l),done,MPI_STATUS_IGNORE,ierror)
                      endif
                   enddo ! now we have a free request holder req(l)
                   call MPI_Isend(aa(is(p + 1), j + js(p + 1) - 1), length, &
                        datatype, p, tag, comm, req(l), ierror)
                endif
                sends = sends + 1
                if(sends .ge. sends_todo) exit loop ! finish if the number of sends is
                ! equal to the total number of sends 
                ! needed
             endif
          enddo
       enddo loop
       deallocate(aa)
    else   ! below the code for non-root:
       ! because of irecv, we must be sure that we have a contiguous matrix
       do j=1,nlocal
          call MPI_Irecv(bb(1,j), mlocal, datatype, root, tag, comm, req(j), ierror)
       enddo
    endif

    call MPI_Waitall(maxreq,req,MPI_STATUSES_IGNORE,ierror)

    if (rank .ne. root) then
       b = bb
    endif

    deallocate(bb)
    deallocate(req)

  end subroutine matrix_distr_send_real8

  subroutine matrix_distr_send_integer(a,b,is,lm,js,ln,root,comm)
    ! distribute matrix on processes
    !
    ! parameters:
    ! integer a(:,:) (in): matrix to be scatterd. 
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
    ! Sends colums of submatrices to the processes using MPI_Isend. 
    ! This routine tries to send as many columns in parallel
    ! to different processes as possible, see the usage of req(:)
    ! and maxreq below.
    ! The columns are received using MPI_Irecv.
    ! 
    use mpi
    implicit none
    integer, intent(in),  dimension(:,:)  :: a
    integer, intent(out), dimension(:,:)  :: b
    integer, intent(in), dimension(:)    :: is,lm
    integer, intent(in), dimension(:)    :: js,ln
    integer, intent(in)                  :: root,comm

    integer                              :: j,rank,ierror,l,procs

    integer                              :: sends, sends_todo, length
    integer                              :: p,mlocal,nlocal

    logical                              :: done

    integer                              :: maxreq
    integer, allocatable, dimension(:)   :: req

    integer, parameter                   :: tag = 1

    integer, parameter                   :: datatype = MPI_INTEGER
    integer, allocatable, dimension(:,:) :: aa,bb


    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    nlocal = ln(rank+1)
    mlocal = lm(rank+1)

    ! calculate total number of sends to be done: that is equal to 
    ! the total number of columns in the submatrices

    if (rank .eq. root) then
       sends_todo = sum(ln(1:procs))
       ! maxreq is the maximum number of outstanding Isend or Irecv requests.
       ! Experiments learned that the value of maxreq is not critical,
       ! a value of 10 gave on the lisa system (Xeon, infiniband) 
       ! only somewhat worse results. Probably, the slower the network,
       ! the higher this value can be set. Maximum usable is the total number
       ! of columns (sends_todo, see above)
       maxreq = 1000
    else
       sends_todo = 0    ! to quiet the compiler
       maxreq = nlocal   ! on non-root processes, we need only the number of local
       ! columns
    endif

    ! distribute matrix. 
    ! to ensure that every process is as busy as possible, we do the
    ! loops on the processes in the inner loop.

    allocate(req(maxreq))
    req=MPI_REQUEST_NULL    ! initialize req with a null-value

    allocate(bb(size(b,1),size(b,2)))
    if ( rank .eq. root) then
       allocate(aa(size(a,1),size(a,2)))
       aa = a
       sends=0
       l = 0
       j = 0
       loop: do 
          j=j+1                 ! j runs over the columns
          do p=0, procs - 1
             if (ln(p + 1) .ge. j) then ! is j in the range of submatrix(p)
                length = lm(p + 1)
                if (p .eq. rank) then    ! submtarix(p) is on root, just copy the column
                   b(1:length,j) = a(is(rank+1):is(rank+1)+length-1,j+js(rank+1)-1)
                else
                   ! search a free request holder:
                   done = .false.
                   do while( .not. done) 
                      l=l+1
                      if (l .gt. maxreq) l=1
                      if (req(l) .eq. MPI_REQUEST_NULL) then
                         done = .true.
                      else
                         !
                         ! MPI_Test gives .true. if req(l).eq.MPI_REQUEST_NULL,
                         ! but the standard says nothing about this case,
                         ! so we test first if req(l) is equal to MPI_REQUEST_NULL
                         !
                         call MPI_Test(req(l),done,MPI_STATUS_IGNORE,ierror)
                      endif
                   enddo ! now we have a free request holder req(l)
                   call MPI_Isend(aa(is(p + 1), j + js(p + 1) - 1), length, &
                        datatype, p, tag, comm, req(l), ierror)
                endif
                sends = sends + 1
                if(sends .ge. sends_todo) exit loop ! finish if the number of sends is
                ! equal to the total number of sends 
                ! needed
             endif
          enddo
       enddo loop
       deallocate(aa)
    else   ! below the code for non-root:
       do j=1,nlocal
          call MPI_Irecv(bb(1,j), mlocal, datatype, root, tag, comm, req(j), ierror)
       enddo
    endif

    do l=1,maxreq   ! wait until all requests are done
       if (req(l) .ne. MPI_REQUEST_NULL) then
          call MPI_Wait(req(l), MPI_STATUS_IGNORE, ierror)
       endif
    enddo

    if (rank .ne. root) then
       b = bb
    endif
    deallocate(bb)

    deallocate(req)

  end subroutine matrix_distr_send_integer

  subroutine matrix_distr_sendmat_real8(a,b,is,lm,js,ln,&
       root,comm)

    ! distribute matrix on root process to all processes
    !
    ! parameters:
    ! real*8 a(:,:) (in):  matrix to be scattered
    ! real*8 b(:,:) (out): local scattered matrix
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
    !  root receives submatrices as a whole, one by one. Local
    !  submatrix is copied. 
    ! 

    use mpi
    implicit none
    real*8, intent(in),  dimension(:,:)  :: a
    real*8, intent(out), dimension(:,:)  :: b
    integer, intent(in), dimension(:)    :: is, lm
    integer, intent(in), dimension(:)    :: js, ln
    integer, intent(in)                  :: root, comm

    integer                              :: i, rank, ierror, procs, p
    integer                              :: mlocal, nlocal
    integer                              :: acolstart, acolend, arowstart, arowend
    integer                              :: mp, np, k, ia, ja
    integer, parameter                   :: tag1 = 125
    integer, parameter                   :: tag2 = 126

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    nlocal = ln(rank+1)
    mlocal = lm(rank+1)

    ! define start-end cols and rows for the general case:
    ! computations on arow.. and acol.. are only meaningful at root

    ! See also spaceparams comments:
    !    b concides with 
    !       a(is(p+1):is(p+1)+lm(p+1)-1,js(p+1):js(p+1)+ln(p+1)-1)

    arowstart  = is(root + 1)
    arowend    = is(root + 1) + mlocal - 1
    acolstart  = js(root + 1)
    acolend    = js(root + 1) + nlocal - 1
    if (rank .eq. root) then
       ! copy of the the relevant part of A to local submtrix B:
       b(:,:) =  a(arowstart:arowend,acolstart:acolend)

       ! now, send to each non-root process 
       ! it's part of matrix a

       ploop: do p = 1, procs
          if ( p - 1 .eq. root ) then
             cycle ploop
          endif
          ! howmany elements to send to process p-1  (= mp*np)
          ! and where to put this in a  (ia:ia+mp-1, ja:ja+np-1)
          !

          mp = lm(p)
          np = ln(p)
          ia = is(p)
          ja = js(p)

          ! to avoid deadlocks, master first sends a
          ! short message to p. After receiving that message,
          ! p will start receiving its submatrix

          call MPI_Send(mp*np, 1, MPI_INTEGER, p - 1, tag1, comm, ierror)

          ! send a submatrix to process p-1:

          call MPI_Send(a(ia:ia+mp-1,ja:ja+np-1), mp*np, MPI_DOUBLE_PRECISION, &
               p - 1, tag2, comm, ierror)

       enddo  ploop ! loop over non-root processes

    else ! above is the root code, now the non-root code
       ! receive the short message from root, so I can start receiving b
       call MPI_Recv(i, 1, MPI_INTEGER,&
            root, tag1, comm, MPI_STATUS_IGNORE, ierror)

       ! sanity check

       k = size(b)
       if (i .ne. k) then
          write(*,*)'process ',rank,' has a problem in matrix_distr_sendmat_real8: '
          write(*,*)'number of words to receive:',k
          write(*,*)'root expects:              ',i
          call MPI_Abort(MPI_COMM_WORLD,1,ierror)
       endif

       call MPI_Recv(b, k, MPI_DOUBLE_PRECISION, root, &
            tag2, comm, MPI_STATUS_IGNORE, ierror)
    endif

  end subroutine matrix_distr_sendmat_real8

  subroutine matrix_distr_sendmat_integer(a,b,is,lm,js,ln,&
       root,comm)

    ! distribute matrix on root process to all processes
    !
    ! parameters:
    ! integer a(:,:) (in):  matrix to be scattered
    ! integer b(:,:) (out): local scattered matrix
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
    !  root receives submatrices as a whole, one by one. Local
    !  submatrix is copied. 
    ! 

    use mpi
    implicit none
    integer, intent(in),  dimension(:,:) :: a
    integer, intent(out), dimension(:,:) :: b
    integer, intent(in), dimension(:)    :: is, lm
    integer, intent(in), dimension(:)    :: js, ln
    integer, intent(in)                  :: root, comm
    integer                              :: i, rank, ierror, procs, p
    integer                              :: mlocal, nlocal
    integer                              :: acolstart, acolend, arowstart, arowend
    integer                              :: mp, np, k, ia, ja
    integer, parameter                   :: tag1 = 125
    integer, parameter                   :: tag2 = 126

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    nlocal = ln(rank+1)
    mlocal = lm(rank+1)

    ! define start-end cols and rows for the general case:
    ! computations on arow.. and acol.. are only meaningful at root

    arowstart  = is(root + 1)
    arowend    = is(root + 1) + mlocal - 1 
    acolstart  = js(root + 1)
    acolend    = js(root + 1) + nlocal - 1

    if (rank .eq. root) then

       ! copy of the the relevant part of A to local submtrix B:

       b(:,:) =  a(arowstart:arowend,acolstart:acolend)

       ! now, send to each non-root process 
       ! it's part of matrix a

       ploop: do p = 1, procs
          if ( p - 1 .eq. root ) then
             cycle ploop
          endif

          ! howmany elements to send to process p-1  (= mp*np)
          ! and where to put this in a  (ia:ia+mp-1, ja:ja+np-1)
          !

          mp = lm(p)
          np = ln(p)
          ia = is(p)
          ja = js(p)

          ! to avoid deadlocks, master first sends a
          ! short message to p - 1. After receiving that message,
          ! p - 1 will start receiving its submatrix

          call MPI_Send(mp*np, 1, MPI_INTEGER, p - 1, tag1, comm, ierror)

          ! send a submatrix to process p-1

          call MPI_Send(a(ia:ia+mp-1,ja:ja+np-1), mp*np, MPI_INTEGER, &
               p - 1, tag2, comm, ierror)

       enddo  ploop ! loop over non-root processes

    else ! above is the root code, below the non-root code

       ! receive the short message from root, so I can start receiving b

       call MPI_Recv(i, 1, MPI_INTEGER,&
            root, tag1, comm, MPI_STATUS_IGNORE, ierror)

       ! sanity check

       k = size(b)
       if (i .ne. k) then
          write(*,*)'process ',rank,' has a problem in matrix_distr_sendmat_integer: '
          write(*,*)'number of words to receive:',k
          write(*,*)'root expects:              ',i
          call MPI_Abort(MPI_COMM_WORLD,1,ierror)
       endif

       call MPI_Recv(b, k, MPI_INTEGER, root, &
            tag2, comm, MPI_STATUS_IGNORE, ierror)

    endif

  end subroutine matrix_distr_sendmat_integer

  subroutine matrix_coll_gather_real8(a,b,is,lm,js,ln, &
       isleft,isright,istop,isbot,root,comm)
    ! collect matrix on root process
    !
    ! parameters:
    ! real*8 a(:,:) (out): matrix to be collected. 
    !                      This matrix has to be available, it is not
    !                      allocated by this subroutine.
    ! real*8 b(:,:) (in): local scattered matrix
    ! integer is,lm,js,ln(number of processes) (in): description of
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
    ! Method
    !
    ! Use gatherv. Each gatherv call collects one column of
    !  all submatrices
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
    use mpi
    implicit none
    real*8, intent(out),  dimension(:,:) :: a
    real*8, intent(in), dimension(:,:)   :: b
    integer, intent(in), dimension(:)    :: is,lm
    integer, intent(in), dimension(:)    :: js,ln
    logical, intent(in), dimension(:)    :: isleft, isright, istop, isbot
    integer, intent(in)                  :: root,comm

    integer                              :: j,p,rank,ierror,procs,jj,ma,jmax,k
    integer                              :: nlocal,mlocal,scnt
    integer, allocatable, dimension(:)   :: displs, cnts
    logical                              :: isleftl, isrightl, istopl, isbotl
    real*8, allocatable, dimension(:,:)  :: aa

#ifdef USEMPE
    call MPE_Log_event(event_coll_start,0,'cstart')
#endif
    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    nlocal   = ln(rank+1)
    mlocal   = lm(rank+1)
    isleftl  = isleft(rank+1)
    isrightl = isright(rank+1)
    istopl   = istop(rank+1)
    isbotl   = isbot(rank+1)

    ma = 1 ! to quiet the compiler
    if (rank .eq. root) then
       ma = size(a,1)
    endif

    if (rank .eq. root) then
       allocate(displs(procs), cnts(procs))
       ! need a contiguous matrix, otherwise the counts and 
       ! displacements will be in error
       allocate(aa(size(a,1),size(a,2)))
    else
       allocate(displs(1),cnts(1))    ! need valid addresses on non-root
       allocate(aa(1,1))
    endif

    nlocal = ln(rank+1)
    mlocal = lm(rank+1)

    jmax = maxval(ln(1:procs))     ! number of gathers needed

    do j=1,jmax
       if (rank .eq. root) then  ! displacements and counts only needed on root
          do p=1,procs
             if (ln(p) .ge. j) then
                displs(p) = (js(p) -1 + j -1)*ma + is(p)     ! compute displacement for
                ! column j of submatrix p
                ! note that first and last elements
                ! of this columns are borders, not
                ! to be gathered
                cnts(p) = lm(p)-2               ! number of non-border elements in
                ! this column
                ! Normally, the first and last rows and columns are not to
                ! be send to root, but we make an exception for the 
                ! rows and columns that are the first rows and columns
                ! from the global matrix.
                ! 
                displs(p) = (js(p) -1 + j - 1)*ma + is(p) - 1
                cnts(p) = lm(p)
                if ( .not. isleft(p) ) then
                   ! this is a submatrix, not on the left.
                   if ( j .eq. 1 ) then
                      ! The first column is  not to be sent, 
                      ! so make cnts(p) = 0, and give displs(p) a
                      ! valid value
                      displs(p) = 0
                      cnts(p)   = 0
                   endif
                endif
                if ( .not. isright(p) ) then
                   ! this is a submatrix not on the right
                   if ( j .eq. ln(p) ) then
                      ! The last column is not to be sent,
                      ! so make cnts(p) = 0, and give displs(p) a
                      ! valid value
                      displs(p) = 0
                      cnts(p)   = 0
                   endif
                endif
                if ( .not. istop(p) ) then
                   ! this is a submatrix, not on the top
                   ! so we don't send the first element of column j
                   ! adapt displs(p) and cnts(p) accordingly:
                   displs(p) = displs(p)+1
                   cnts(p) = max(cnts(p)-1,0)
                endif
                if ( .not. isbot(p) ) then
                   ! this is a submatrix, not on the bottom
                   ! so we don't send the last element of column j
                   ! adapt cnts(p) accordingly:
                   cnts(p) = max(cnts(p)-1,0)
                endif
             else
                displs(p) = 0  ! we are out of columns, use a valid value for
                ! the displacement
                cnts(p) = 0    ! and set number of elements to receive for 
                ! this submatrix to zero
             endif  !ln(p) .ge. j
          enddo  ! do p = 1,procs
       endif  ! rank .eq. root
       jj   = j            ! column of submatrix to send
       scnt = mlocal       ! number of elements to send
       k    = 1            ! starting element of column j to send
       if (nlocal .ge. j) then 
          if (.not. isleftl ) then
             if ( j .eq. 1 ) then
                jj   = 1
                scnt = 0
             endif
          endif
          if (.not. isrightl ) then
             if ( j .eq. nlocal ) then
                jj = 1
                scnt = 0
             endif
          endif
          if (.not. istopl ) then
             k = k+1
             scnt = max(scnt-1, 0)
          endif
          if (.not. isbotl ) then
             scnt = max(scnt-1, 0)
          endif
       else
          jj   = 1              ! out of columns, use valid value for jj
          scnt = 0              ! and put number of elements to send to zero
       endif

       call MPI_Gatherv(b(k:k+scnt-1,jj), scnt, MPI_DOUBLE_PRECISION,  &
            aa(1,1), cnts, displs, MPI_DOUBLE_PRECISION,  &
            root, comm, ierror)
    enddo

    if(rank .eq. root) then
       a = aa
    endif

    deallocate(displs,cnts,aa)

  end subroutine matrix_coll_gather_real8

  subroutine matrix_coll_recv_real8    (a,b,is,lm,js,ln,&
       isleft,isright,istop,isbot,root,comm)
    ! collect matrix on root process
    !
    ! parameters:
    ! real*8 a(:,:) (out): matrix to be collected. 
    !                      This matrix has to be available, it is not
    !                      allocated by this subroutine.
    ! real*8 b(:,:) (in): local scattered matrix
    ! integer is,lm,js,ln(number of processes) (in): description of
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
    ! Method
    !
    ! Receive colums of submatrices using MPI_Probe and MPI_Recv. 
    ! The non-root processes use MPI_Isend to send the columns
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
    use mpi
    implicit none
    real*8, intent(out), dimension(:,:)  :: a
    real*8, intent(in),  dimension(:,:)  :: b
    integer, intent(in), dimension(:)    :: is,lm
    integer, intent(in), dimension(:)    :: js,ln
    logical, intent(in), dimension(:)    :: isleft,isright,istop,isbot
    integer, intent(in)                  :: root,comm

    integer                              :: i,j,rank,ierror,l,procs,ii,jj
    integer                              :: jstart, jend, p, m, n

    integer                              :: recvs, recvs_todo
    integer                              :: mlocal,nlocal,sender
    logical                              :: isleftl, isrightl, istopl, isbotl

    integer, allocatable, dimension(:)   :: col
    integer, allocatable, dimension(:)   :: req

    integer, parameter                   :: tag = 2
    integer, dimension(MPI_STATUS_SIZE)  :: mpistatus
    real*8, allocatable, dimension(:,:)  :: bb

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    nlocal   = ln(rank+1)
    mlocal   = lm(rank+1)
    isleftl  = isleft(rank+1)
    isrightl = isright(rank+1)
    istopl   = istop(rank+1)
    isbotl   = isbot(rank+1)

    if (rank .eq. root) then
       ! calculate total number of recvs to be done: that is equal to 
       ! the total number of columns in the submatrices 
       recvs_todo = sum(ln(1:procs))

       recvs_todo =  recvs_todo - nlocal  ! correct for columns local available

       do p=1,procs                       ! correct for borders
          if (p .eq. root+1) then
             cycle
          endif
          if ( .not. isleft(p)) then
             recvs_todo = recvs_todo -1
          endif
          if ( .not. isright(p)) then
             recvs_todo = recvs_todo -1
          endif
       enddo
       allocate(col(procs)) ! to remember which was the last column
       ! received for submatrix(p)
       col   = 0
       recvs = 0
       ! copy submatrix on root:
       i = 1
       ii = is(rank+1)
       m = mlocal
       j = 1
       jj = js(rank+1)
       n = nlocal
       ! take care of borders:
       if ( .not. isleftl ) then
          j  = j  + 1
          jj = jj + 1
       endif
       if ( .not. isrightl ) then
          n = n - 1
       endif
       if ( .not. istopl ) then
          i  = i  + 1
          ii = ii + 1
          m  = m  - 1
       endif
       if ( .not. isbotl ) then
          m = m - 1
       endif
       a(ii:ii+m-1,jj:jj+n-1) = b(i:i+m-1,j:j+n-1)

       do recvs=1,recvs_todo  ! loop over the total number of receives needed
          ! wait until anybody did send a column to root:

          call MPI_Probe(MPI_ANY_SOURCE, tag, comm, mpistatus, ierror)

          ! find out who was the sender:

          sender = mpistatus(MPI_SOURCE)

          ! update col, and compute i and j: the startlocation of 
          ! to be received column

          p      = sender + 1
          col(p) = col(p) + 1
          i      = is(p)
          j      = col(p)+js(p)-1
          m      = lm(p)

          ! perform the actual receive:
          ! and correct for borders

          if ( .not. isleft(p)) then
             j = j + 1
          endif
          if ( .not. isright(p)) then
             ! nothing to correct here
          endif
          if ( .not. istop(p)) then
             i = i + 1
             m = m - 1
          endif
          if ( .not. isbot(p)) then
             m = m - 1
          endif

          call MPI_Recv(a(i:i+m-1,j),m,MPI_DOUBLE_PRECISION, &
               sender,tag,comm,MPI_STATUS_IGNORE,ierror)
       enddo
       deallocate(col)
       !    deallocate(isleft, isright, istop, isbot)
    else  ! non-root code
       allocate(req(nlocal))  ! for holding the requests
       allocate(bb(size(b,1),size(b,2)))
       bb = b     ! want to be sure that the matrix to send is
       ! contiguous, because isend's are used
       req = MPI_REQUEST_NULL

       jstart = 1
       jend   = nlocal
       i      = 1
       m      = mlocal

       ! take care of borders
       if ( .not. isleftl ) then
          jstart = jstart + 1
       endif
       if ( .not. isrightl ) then
          jend = jend - 1
       endif
       if ( .not. istopl ) then
          i = i + 1
          m = m - 1
       endif
       if ( .not. isbotl ) then
          m = m - 1
       endif

       do j = jstart, jend
          call MPI_Isend(bb(i,j), m, MPI_DOUBLE_PRECISION, root, tag, comm, req(j), ierror)
       enddo

       do l=1,nlocal ! wait until all requests have been handled
          if (req(l) .ne. MPI_REQUEST_NULL) then
             call MPI_Wait(req(l),MPI_STATUS_IGNORE,ierror)
          endif
       enddo
       deallocate(req,bb)
    endif

  end subroutine matrix_coll_recv_real8

  subroutine matrix_coll_recvmat_real8(a,b,is,lm,js,ln,&
       isleft,isright,istop,isbot,root,comm)
    ! collect matrix on root process
    !
    ! parameters:
    ! real*8 a(:,:) (out): matrix to be collected. 
    !                      This matrix has to be available, it is not
    !                      allocated by this subroutine.
    ! real*8 b(:,:) (in): local scattered matrix
    ! integer is,lm,js,ln(number of processes) (in): description of
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
    ! Method
    !
    !  root receives submatrices as a whole, one by one. Local
    !  submatrix is copied. 
    ! 
    !  NOTE: the matrices are considered to be bordered by columns
    !   left and right and rows up and down. 
    !   The number of border columns/rows is nbord.
    !   The number of overlapping columns nover = 2*nbord.
    !   So, in general, 
    !   only b(nbord+1:mb-nbord,nbord+1:nb-nbord) is send to the master process. 
    !   HOWEVER: if matrix b is 
    !   on the left: b(:,1:nbord) is also sent 
    !   on the right: b(:,nb-nbord+1:nb) is also sent
    !   on the top:   b(1:nbord,:)  is also sent
    !   on the bottom: b(mb-nbord+1:mb,:) is also sent
    ! The master process receives a(1:ma,1:na). 
    !   dimension of b (1:mb,1:nb). 
    !   dimension of a (1:ma,1:na)
    !
    use mpi
    implicit none
    real*8, intent(out), dimension(:,:)  :: a
    real*8, intent(in),  dimension(:,:)  :: b
    integer, intent(in), dimension(:)    :: is, lm
    integer, intent(in), dimension(:)    :: js, ln
    logical, intent(in), dimension(:)    :: isleft, isright, istop, isbot
    integer, intent(in)                  :: root, comm

    integer                              :: i, rank, ierror, procs, p

    integer                              :: mlocal, nlocal

    !logical, allocatable, dimension(:)   :: isleft, isright, istop, isbot

    integer, parameter                   :: tag = 2

    integer                              :: acolstart, acolend, arowstart, arowend
    integer                              :: bcolstart, bcolend, browstart, browend
    integer                              :: mp, np, k, ia, ja
    integer, parameter                   :: tag1 = 123
    integer, parameter                   :: tag2 = 124
    logical                              :: isleftl, isrightl, istopl, isbotl

    ! nover is number of overlapping rows/columns of submatrices.
    ! nbord is number of borders.

    integer, parameter                   :: nbord = 2
    integer, parameter                   :: nover = 2*nbord

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    nlocal   = ln(rank+1)
    mlocal   = lm(rank+1)
    isleftl  = isleft(rank+1)
    isrightl = isright(rank+1)
    istopl   = istop(rank+1)
    isbotl   = isbot(rank+1)

    ! define start-end cols and rows for the general case:
    ! computations on arow.. and acol.. are only meaningful at root
    browstart  = nbord  + 1
    browend    = mlocal - nbord
    bcolstart  = nbord + 1
    bcolend    = nlocal - nbord
    arowstart  = is(root + 1) + nbord
    arowend    = is(root + 1) + mlocal - 1 - nbord
    acolstart  = js(root + 1) + nbord
    acolend    = js(root + 1) + nlocal - 1 - nbord

    ! correct for the border cases:

    if (isleftl) then
       bcolstart = bcolstart - nbord
       acolstart = acolstart - nbord
    endif
    if (isrightl) then
       bcolend = bcolend + nbord
       acolend = acolend + nbord
    endif
    if (istopl) then
       browstart = browstart - nbord
       arowstart = arowstart - nbord
    endif
    if (isbotl) then
       browend = browend + nbord
       arowend = arowend + nbord
    endif

    if (rank .eq. root) then

       ! copy local submatrix 

       ! coy of the local submtrix B to A:
       a(arowstart:arowend,acolstart:acolend) = &
            b(browstart:browend,bcolstart:bcolend)

       ! now, collect from each non-root process 
       ! it's matrix b, and put it in the proper place in a

       ploop: do p = 1, procs
          if ( p - 1 .eq. root ) then
             cycle ploop
          endif
          ! howmany elements to expect from process p-1  (= mp*np)
          ! and where to put this in a  (ia:ia+mp-1, ja:ja+np-1)
          !

          mp = lm(p) - nover
          np = ln(p) - nover
          ia = is(p) + nbord
          ja = js(p) + nbord

          ! correction if b is at a border:

          if (isleft(p)) then
             ja = ja - nbord
             np = np + nbord
          endif
          if (isright(p)) then
             np = np + nbord
          endif
          if (istop(p)) then
             ia = ia - nbord
             mp = mp + nbord
          endif
          if (isbot(p)) then
             mp = mp + nbord
          endif

          ! to avoid deadlocks, master first sends a
          ! short message to p. After receiving that message,
          ! p will start sending its submatrix

          call MPI_Send(mp*np, 1, MPI_INTEGER, p - 1, tag1, comm, ierror)

          ! receive a submatrix from process p-1

          call MPI_Recv(a(ia:ia+mp-1,ja:ja+np-1), mp*np, MPI_DOUBLE_PRECISION, &
               p - 1, tag2, comm, MPI_STATUS_IGNORE, ierror)

       enddo  ploop ! loop over non-root processes
    else ! above is the root code, now the non-root code
       ! receive the short message from root, so I can start sending b
       call MPI_Recv(i, 1, MPI_INTEGER,&
            root, tag1, comm, MPI_STATUS_IGNORE, ierror)

       ! which part of b to send: b(browstart:browend,bcolstart:bcolend) 

       ! sanity check
       k = (browend - browstart + 1)*(bcolend - bcolstart + 1)
       if (i .ne. k) then
          write(*,*)'process ',rank,' has a problem in matrix_coll_recvmat_real8: '
          write(*,*)'number of words to send:',k
          write(*,*)'root expects:           ',i
          do i=1,procs
             write(*,*)i-1,isleft(i),isright(i),istop(i),isbot(i)
          enddo
          write(*,*)isleftl,isrightl,istopl,isbotl
          call MPI_Abort(MPI_COMM_WORLD,1,ierror)
       endif

       call MPI_Send(b(browstart:browend,bcolstart:bcolend), k, &
            MPI_DOUBLE_PRECISION, root, tag2, comm, ierror)
    endif

  end subroutine matrix_coll_recvmat_real8

  subroutine matrix_coll_recvmat_int(a,b,is,lm,js,ln,&
       isleft,isright,istop,isbot,root,comm)
    ! This is a copy of matrix_coll_recvmat_real8 for integers
    use mpi
    implicit none
    integer, intent(out), dimension(:,:)  :: a
    integer, intent(in),  dimension(:,:)  :: b
    integer, intent(in), dimension(:)    :: is, lm
    integer, intent(in), dimension(:)    :: js, ln
    logical, intent(in), dimension(:)    :: isleft, isright, istop, isbot
    integer, intent(in)                  :: root, comm

    integer                              :: i, rank, ierror, procs, p

    integer                              :: mlocal, nlocal

    !logical, allocatable, dimension(:)   :: isleft, isright, istop, isbot

    integer, parameter                   :: tag = 2

    integer                              :: acolstart, acolend, arowstart, arowend
    integer                              :: bcolstart, bcolend, browstart, browend
    integer                              :: mp, np, k, ia, ja
    integer, parameter                   :: tag1 = 123
    integer, parameter                   :: tag2 = 124
    logical                              :: isleftl, isrightl, istopl, isbotl

    ! nover is number of overlapping rows/columns of submatrices.
    ! nbord is number of borders.

    integer, parameter                   :: nover=4
    integer, parameter                   :: nbord=nover/2
    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)

    nlocal   = ln(rank+1)
    mlocal   = lm(rank+1)
    isleftl  = isleft(rank+1)
    isrightl = isright(rank+1)
    istopl   = istop(rank+1)
    isbotl   = isbot(rank+1)

    ! define start-end cols and rows for the general case:
    ! computations on arow.. and acol.. are only meaningful at root
    browstart  = nbord  + 1
    browend    = mlocal - nbord
    bcolstart  = nbord + 1
    bcolend    = nlocal - nbord
    arowstart  = is(root + 1) + nbord
    arowend    = is(root + 1) + mlocal - 1 - nbord
    acolstart  = js(root + 1) + nbord
    acolend    = js(root + 1) + nlocal - 1 - nbord

    ! correct for the border cases:

    if (isleftl) then
       bcolstart = bcolstart - nbord
       acolstart = acolstart - nbord
    endif
    if (isrightl) then
       bcolend = bcolend + nbord
       acolend = acolend + nbord
    endif
    if (istopl) then
       browstart = browstart - nbord
       arowstart = arowstart - nbord
    endif
    if (isbotl) then
       browend = browend + nbord
       arowend = arowend + nbord
    endif

    if (rank .eq. root) then

       ! copy local submatrix 

       ! coy of the local submtrix B to A:

       a(arowstart:arowend,acolstart:acolend) = &
            b(browstart:browend,bcolstart:bcolend)

       ! now, collect from each non-root process 
       ! it's matrix b, and put it in the proper place in a

       ploop: do p = 1, procs
          if ( p - 1 .eq. root ) then
             cycle ploop
          endif
          ! howmany elements to expect from process p-1  (= mp*np)
          ! and where to put this in a  (ia:ia+mp-1, ja:ja+np-1)
          !

          mp = lm(p) - nover
          np = ln(p) - nover
          ia = is(p) + nbord
          ja = js(p) + nbord

          ! correction if b is at a border:

          if (isleft(p)) then
             ja = ja - nbord
             np = np + nbord
          endif
          if (isright(p)) then
             np = np + nbord
          endif
          if (istop(p)) then
             ia = ia - nbord
             mp = mp + nbord
          endif
          if (isbot(p)) then
             mp = mp + nbord
          endif

          ! to avoid deadlocks, master first sends a
          ! short message to p. After receiving that message,
          ! p will start sending its submatrix

          call MPI_Send(mp*np, 1, MPI_INTEGER, p - 1, tag1, comm, ierror)

          ! receive a submatrix from process p-1

          call MPI_Recv(a(ia:ia+mp-1,ja:ja+np-1), mp*np, MPI_INTEGER, &
               p - 1, tag2, comm, MPI_STATUS_IGNORE, ierror)

       enddo  ploop ! loop over non-root processes
    else ! above is the root code, now the non-root code
       ! receive the short message from root, so I can start sending b
       call MPI_Recv(i, 1, MPI_INTEGER,&
            root, tag1, comm, MPI_STATUS_IGNORE, ierror)

       ! which part of b to send: b(browstart:browend,bcolstart:bcolend) 

       ! sanity check
       k = (browend - browstart + 1)*(bcolend - bcolstart + 1)
       if (i .ne. k) then
          write(*,*)'process ',rank,' has a problem in matrix_coll_recvmat_int: '
          write(*,*)'number of words to send:',k
          write(*,*)'root expects:           ',i
          do i=1,procs
             write(*,*)i-1,isleft(i),isright(i),istop(i),isbot(i)
          enddo
          write(*,*)isleftl,isrightl,istopl,isbotl
          call MPI_Abort(MPI_COMM_WORLD,1,ierror)
       endif

       call MPI_Send(b(browstart:browend,bcolstart:bcolend), k, &
            MPI_INTEGER, root, tag2, comm, ierror)
    endif

  end subroutine matrix_coll_recvmat_int

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
