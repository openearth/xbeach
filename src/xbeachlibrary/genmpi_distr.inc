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
!
! every process receives from root

! all processes will receive only from root

sizes     = shape(b)
subsizes  = (/ lm(ra+1), ln(ra+1) /)
starts    = 0
call MPI_Type_create_subarray(2,sizes,subsizes,starts, &
    MPI_ORDER_FORTRAN, basic_type,recvtypes(1), ier)
call MPI_Type_commit(recvtypes(1),ier)
recvcounts(1) = 1

! determine mpi types for the senders
!
! only root will send
!
if(ra == 0) then
  do i=1,sz
    sizes    = shape(a)
    subsizes = (/ lm(i)     ,ln(i)     /)
    starts   = (/ is(i) - 1 ,js(i) - 1 /)
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

!directions for vi vim: filetype=fortran : syntax=fortran
