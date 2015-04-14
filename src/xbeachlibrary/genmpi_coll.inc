!
call MPI_Comm_size(comm, sz, ier)
call MPI_Comm_rank(comm, ra, ier)

if (root /= 0) then
  print *,'Error in matrix_distr_real8: root must be 0, but is:',root
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

if(ra == 0) then
  ! the root process will not receive from itself
  ! but will receive from everybody else
  do i=2,sz
    sizes    = shape(a)
    subsizes = (/ lm(i-1) - 2*nbord   , ln(i-1) - 2*nbord    /)
    starts   = (/ is(i-1) +   nbord -1, js(i-1) +   nbord -1 /)
    if (istop(i-1)) then
      subsizes(1) = subsizes(1) + nbord
      starts(1)   = starts(1)   - nbord
    endif
    if (isbot(i-1)) then
      subsizes(1) = subsizes(1) + nbord
    endif
    if (isleft(i-1)) then
      subsizes(2) = subsizes(2) + nbord
      starts(2)   = starts(2)   - nbord
    endif
    if (isright(i-1)) then
      subsizes(2) = subsizes(2) + nbord
    endif
    call MPI_Type_create_subarray(2,sizes,subsizes,starts,  &
        MPI_ORDER_FORTRAN, basic_type,recvtypes(i),ier)
    call MPI_Type_commit(recvtypes(i),ier)
    recvcounts(i) = 1
  enddo
endif

! determine mpi types for the senders
! root does not send anythink
if(ra /= 0) then
  ! non-root processes send only to root
  sizes      = shape(b)
  subsizes   = (/ lm(ra) - 2*nbord    , ln(ra) - 2*nbord     /)
  starts     = (/ nbord, nbord/)
  if (istop(ra)) then
    subsizes(1) = subsizes(1) + nbord
    starts(1)   = starts(1)   - nbord
  endif
  if (isbot(ra)) then
    subsizes(1) = subsizes(1) + nbord
  endif
  if (isleft(ra)) then
    subsizes(2) = subsizes(2) + nbord
    starts(2)   = starts(2)   - nbord
  endif
  if (isright(ra)) then
    subsizes(2) = subsizes(2) + nbord
  endif
  call MPI_Type_create_subarray(2,sizes,subsizes,starts,  &
      MPI_ORDER_FORTRAN, basic_type,sendtypes(1),ier)
  call MPI_Type_commit(sendtypes(1),ier)
  sendcounts(1) = 1
endif

call MPI_Alltoallw(b,sendcounts,sdispls,sendtypes, &
    a,recvcounts,rdispls,recvtypes,comm,ier)

do i=1,sz
  if (sendtypes(i) /= MPI_CHARACTER) then
    call MPI_Type_free(sendtypes(i),ier)
  endif
  if (recvtypes(i) /= MPI_CHARACTER) then
    call MPI_Type_free(recvtypes(i),ier)
  endif
enddo
!directions for vi vim: filetype=fortran
