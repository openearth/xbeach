module debugging_module
#if 0
! these routines need a thorough rewrite
implicit none
#ifdef USEMPI
  interface compare
     module procedure comparei2
     module procedure comparer2
     module procedure comparer3
  end interface compare
#endif

  interface printsum
     module procedure printsum0
     module procedure printsum1
     module procedure printsum2
     module procedure printsum3
     module procedure printsum4
     module procedure printsumi0
     module procedure printsumi1
     module procedure printsumi2
     module procedure printsumi3
     module procedure printsumi4
  end interface printsum

contains                                         

#ifdef USEMPI
  !
  ! consistency check: in general, all distributed matrices and blocks
  ! should have the following properties:
  !
  !  a  is distributed matrix in this process
  !  al is distributed matrix in left neighbour process
  !  ar is distributed matrix in right neighbour process
  !  at is distributed matrix in top neighbour process
  !  ab is distributed matrix in bottom neighbour process
  !
  !  a(:,2) = al(:,ny+1)
  !  a(2,:) = at(nx+1,:)
  !  except the first and last elements of these arrays
  !  
  !
  !  When all these checks are ok on every process, automatically care 
  ! has been taken for the other neighbours

  !
  ! consistency check for mnem
  ! if mnem = 'ALL' then all
  subroutine space_consistency(s,mnem,verbose)
    use mnemmodule
    use spaceparamsdef
    implicit none
    type(spacepars)                        :: s
    character(len=*)                       :: mnem
    character(len=*), intent(in), optional :: verbose
    integer          :: j,jmin,jmax
    type(arraytype)  :: t

    if(mnem .eq. 'ALL') then
       jmin = 1
       jmax = numvars
    else
       jmin = chartoindex(mnem)
       jmax = jmin
    endif
    do j=jmin,jmax
       call indextos(s,j,t)
       select case(t%type)
       case('r')
          select case (t%rank)
          case(0)
             !call compare(t%r0,t%name)
          case(1)
             !call compare(t%r1,t%name)
          case(2)
             call compare(t%r2,t%name,verbose)
          case(3)
             call compare(t%r3,t%name,verbose)
          case(4)
             !call compare(t%r4,t%name)
          end select ! rank
       case('i')
          select case (t%rank)
          case(0)
             !call compare(t%i0,t%name)
          case(1)
             !call compare(t%i1,t%name)
          case(2)
             call compare(t%i2,t%name,verbose)
          case(3)
             !call compare(t%i3,t%name)
          case(4)
             !call compare(t%i4,t%name)
          end select ! rank
       end select ! type
    enddo

  end subroutine space_consistency

  subroutine comparer2(x,s,verbose)
    use xmpi_module
    use mnemmodule
    implicit none
    real*8, dimension(:,:)    :: x
    character(len=*)          :: s
    character(len=*),optional :: verbose

    real*8, parameter                   :: eps=1.0d-60
    integer                             :: m
    integer                             :: n
    real*8, dimension(:,:), allocatable :: c
    real*8, dimension(:,:), allocatable :: r
    real*8, dimension(2)                :: dif,difmax
    character*100                       :: warning
    integer i,j

    select case(s)
    case (mnem_tideinpz)
       return
    end select
    m=size(x,1)
    n=size(x,2)
    allocate(c(m,2))
    allocate(r(2,n))
    !c = x(:,1)
    c = x(:,1:2)
    !call xmpi_shift(x,':1')
    call xmpi_shift(x,SHIFT_Y_R,1,2)

    !dif(1) = sum(abs(c(2:m-1)-x(2:m-1,1)))
    !dif(2) = sum(abs(c(1:m  )-x(1:m  ,1)))
    dif(1) = sum(abs(c(2:m-1,:)-x(2:m-1,1:2)))
    dif(2) = sum(abs(c(1:m  ,:)-x(1:m  ,1:2)))

    call xmpi_allreduce(dif,MPI_SUM)
    difmax = dif

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (:,1:2) '//trim(s)//': ',difmax,trim(warning)
    endif
    if (present(verbose)) then
      if (verbose .eq. 'verbose') then
        if (sum(difmax) .gt. eps) then
          do i=1,m
            do j=1,2
              !if(abs(c(i,j)-x(i,j)) .gt. 0.0001) then
              if(c(i,j) .ne. 0 .or. x(i,j) .ne. 0) then
                print *,xmpi_rank,i,j,c(i,j),x(i,j),abs(c(i,j)-x(i,j))
              endif
            enddo
          enddo
        endif
      endif
    endif

    x(:,1:2) = c

    !c = x(:,n)
    c = x(:,n-1:n)
    !call xmpi_shift(x,':n')
    call xmpi_shift(x,SHIFT_Y_L,3,4)

    !dif(1) = sum(abs(c(2:m-1)-x(2:m-1,n)))
    !dif(2) = sum(abs(c(1:m  )-x(1:m  ,n)))
    dif(1) = sum(abs(c(2:m-1,:)-x(2:m-1,n-1:n)))
    dif(2) = sum(abs(c(1:m  ,:)-x(1:m  ,n-1:n)))
    !x(:,n) = c
    x(:,n-1:n) = c

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (:,n-1:n) '//trim(s)//': ',difmax,trim(warning)
    endif

    !r = x(1,:)
    r = x(1:2,:)
    !call xmpi_shift(x,'1:')
    call xmpi_shift(x,SHIFT_X_D,1,2)

    !dif(1) = sum(abs(r(2:n-1)-x(1,2:n-1)))
    !dif(2) = sum(abs(r(1:n  )-x(1,1:n  )))
    dif(1) = sum(abs(r(:,2:n-1)-x(1:2,2:n-1)))
    dif(2) = sum(abs(r(:,1:n  )-x(1:2,1:n  )))
    x(1:2,:) = r

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (1:2,:) '//trim(s)//': ',difmax,trim(warning)
    endif

    !r = x(m,:)
    r = x(m-1:m,:)
    !call xmpi_shift(x,'m:')
    call xmpi_shift(x,SHIFT_X_U,3,4)

    !dif(1) = sum(abs(r(2:n-1)-x(m,2:n-1)))
    !dif(2) = sum(abs(r(1:n  )-x(m,1:n  )))
    dif(1) = sum(abs(r(:,2:n-1)-x(m-1:m,2:n-1)))
    dif(2) = sum(abs(r(:,1:n  )-x(m-1:m,1:n  )))
    x(m-1:m,:) = r

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (m-1:m,:) '//trim(s)//': ',difmax,trim(warning)
    endif
  end subroutine comparer2

  subroutine comparei2(x,s,verbose)
    use xmpi_module
    implicit none
    integer, dimension(:,:)   :: x
    character(len=*)          :: s
    character(len=*),optional :: verbose

    integer, parameter    :: eps=0
    integer               :: m
    integer               :: n
    integer, dimension(:), allocatable :: c
    integer, dimension(:), allocatable :: r
    integer, dimension(2) :: dif,difmax
    character*100         :: warning

    print *,'comparei2 not adapted for double border scheme'
    return

    if (present(verbose)) then
      print *,'comparei2: verbose not implemented'
    endif

    m=size(x,1)
    n=size(x,2)
    allocate(c(m))
    allocate(r(n))
    c = x(:,1)
    call xmpi_shift(x,':1')

    dif(1) = sum(abs(c(2:m-1)-x(2:m-1,1)))
    dif(2) = sum(abs(c(1:m  )-x(1:m  ,1)))
    x(:,1) = c

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (:,1) '//trim(s)//': ',difmax,trim(warning)
    endif

    c = x(:,n)
    call xmpi_shift(x,':n')

    dif(1) = sum(abs(c(2:m-1)-x(2:m-1,n)))
    dif(2) = sum(abs(c(1:m  )-x(1:m  ,n)))
    x(:,n) = c

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (:,n) '//trim(s)//': ',difmax,trim(warning)
    endif
    r = x(1,:)
    call xmpi_shift(x,'1:')

    dif(1) = sum(abs(r(2:n-1)-x(1,2:n-1)))
    dif(2) = sum(abs(r(1:n  )-x(1,1:n  )))
    x(1,:) = r

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (1,:) '//trim(s)//': ',difmax,trim(warning)
    endif

    r = x(m,:)
    call xmpi_shift(x,'m:')

    dif(1) = sum(abs(r(2:n-1)-x(m,2:n-1)))
    dif(2) = sum(abs(r(1:n  )-x(m,1:n  )))
    x(m,:) = r

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (m,:) '//trim(s)//': ',difmax,trim(warning)
    endif
  end subroutine comparei2

  subroutine comparer3(x,s,verbose)
    use xmpi_module
    implicit none
    real*8, dimension(:,:,:)  :: x
    character(len=*)          :: s
    character(len=*),optional :: verbose

    real*8, parameter     :: eps=1.0d-60
    integer               :: m,n,l
    real*8, dimension(:,:), allocatable :: c
    real*8, dimension(:,:), allocatable :: r
    real*8, dimension(2)  :: dif,difmax
    character*100         :: warning

    print *,'comparei2 not adapted for double border scheme'
    return
    if (present(verbose)) then
      print *,'comparer3: verbose not implemented'
    endif
    m=size(x,1)
    n=size(x,2)
    l=size(x,3)

    allocate(c(m,l))
    allocate(r(n,l))
    c = x(:,1,:)
    call xmpi_shift(x,':1')

    dif(1) = sum(abs(c(2:m-1,:)-x(2:m-1,1,:)))
    dif(2) = sum(abs(c(1:m,:)  -x(1:m  ,1,:)))
    x(:,1,:) = c

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (:,1) '//trim(s)//': ',difmax,trim(warning)
    endif

    c = x(:,n,:)
    call xmpi_shift(x,':n')

    dif(1) = sum(abs(c(2:m-1,:)-x(2:m-1,n,:)))
    dif(2) = sum(abs(c(1:m  ,:)-x(1:m  ,n,:)))
    x(:,n,:) = c

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (:,n) '//trim(s)//': ',difmax,trim(warning)
    endif

    r = x(1,:,:)
    call xmpi_shift(x,'1:')

    dif(1) = sum(abs(r(2:n-1,:)-x(1,2:n-1,:)))
    dif(2) = sum(abs(r(1:n  ,:)-x(1,1:n  ,:)))
    x(1,:,:) = r

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (1,:) '//trim(s)//': ',difmax,trim(warning)
    endif

    r = x(m,:,:)
    call xmpi_shift(x,'m:')

    dif(1) = sum(abs(r(2:n-1,:)-x(m,2:n-1,:)))
    dif(2) = sum(abs(r(1:n  ,:)-x(m,1:n  ,:)))
    x(m,:,:) = r

    call xmpi_reduce(dif,difmax,MPI_SUM)

    if(xmaster) then
       warning=' '
       if (sum(difmax) .gt. eps) then
          warning = '<===++++++++++++++++++++++'
       endif
       write (*,*) 'compare (m,:) '//trim(s)//': ',difmax,trim(warning)
    endif
  end subroutine comparer3
  !
#endif

  subroutine printsum0(f,str,id,val)
    implicit none
    integer, intent(in)     :: f
    character(*),intent(in) :: str
    integer, intent(in)     :: id
    real*8, intent(in)      :: val
    write(f,*) 'printsum ',id,' ',str,':',val
  end subroutine printsum0

  subroutine printsum1(f,str,id,val)
    implicit none
    integer, intent(in) :: f
    character(*), intent(in) :: str
    integer, intent(in) :: id
    real*8, pointer, dimension(:) :: val
    if (associated(val) ) then
       write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
    else
       write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
    endif
  end subroutine printsum1

  subroutine printsum2(f,str,id,val)
    implicit none
    integer, intent(in) :: f
    character(*), intent(in) :: str
    integer, intent(in) :: id
    real*8, pointer, dimension(:,:) :: val
    if (associated(val) ) then
       write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
    else
       write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
    endif

  end subroutine printsum2

  subroutine printsum3(f,str,id,val)
    implicit none
    integer, intent(in) :: f
    character(*), intent(in) :: str
    integer, intent(in) :: id
    real*8, pointer, dimension(:,:,:) :: val
    if (associated(val) ) then
       write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
    else
       write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
    endif
  end subroutine printsum3

  subroutine printsum4(f,str,id,val)
    implicit none
    integer, intent(in) :: f
    character(*), intent(in) :: str
    integer, intent(in) :: id
    real*8, pointer, dimension(:,:,:,:) :: val
    if (associated(val) ) then
       write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
    else
       write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
    endif
  end subroutine printsum4

  subroutine printsumi0(f,str,id,val)
    implicit none
    integer, intent(in) :: f
    character(*),intent(in) :: str
    integer, intent(in) :: id
    integer, intent(in) :: val
    write(f,*) 'printsum ',id,' ',str,':',val
  end subroutine printsumi0

  subroutine printsumi1(f,str,id,val)
    implicit none
    integer, intent(in) :: f
    character(*), intent(in) :: str
    integer, intent(in) :: id
    integer, pointer, dimension(:) :: val
    if (associated(val) ) then
       write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
    else
       write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
    endif
  end subroutine printsumi1

  subroutine printsumi2(f,str,id,val)
    implicit none
    integer, intent(in) :: f
    character(*), intent(in) :: str
    integer, intent(in) :: id
    integer, pointer, dimension(:,:) :: val
    if (associated(val) ) then
       write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
    else
       write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
    endif

  end subroutine printsumi2

  subroutine printsumi3(f,str,id,val)
    implicit none
    integer, intent(in) :: f
    character(*), intent(in) :: str
    integer, intent(in) :: id
    integer, pointer, dimension(:,:,:) :: val
    if (associated(val) ) then
       write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
    else
       write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
    endif
  end subroutine printsumi3

  subroutine printsumi4(f,str,id,val)
    implicit none
    integer, intent(in) :: f
    character(*), intent(in) :: str
    integer, intent(in) :: id
    integer, allocatable, dimension(:,:,:,:) :: val
    if (allocated(val) ) then
       write(f,*) 'printsum ',id,' ',str,':',sum(val),shape(val)
    else
       write(f,*) 'printsum ',id,' ',str,':',' Not allocated'
    endif
  end subroutine printsumi4

  subroutine printssums(s,str)
#ifdef USEMPI
    use xmpi_module
#endif
    use mnemmodule
    use spaceparams
    type (spacepars), intent(in) :: s
    character(*), intent(in) :: str
    integer :: id, f,i
    type(arraytype) :: t

#ifdef USEMPI
    id = xmpi_rank
    if (id .gt. 0 ) then
       return
    endif
#else
    id=0
#endif

    f = 100+4*numvars+id+1

    write(f, *) 'printsum: ',id,'Start of printssums ',str
    do i=1,numvars
       call indextos(s,i,t)
       select case (t%rank)
       case (0)
          if (t%type .eq. 'i') then
             call printsum(f,t%name,id,t%i0)
          else
             call printsum(f,t%name,id,t%r0)
          endif
       case (1)
          if (t%type .eq. 'i') then
             call printsum(f,t%name,id,t%i1)
          else
             call printsum(f,t%name,id,t%r1)
          endif
       case (2)
          if (t%type .eq. 'i') then
             call printsum(f,t%name,id,t%i2)
          else
             call printsum(f,t%name,id,t%r2)
          endif
       case (3)
          if (t%type .eq. 'i') then
             call printsum(f,t%name,id,t%i3)
          else
             call printsum(f,t%name,id,t%r3)
          endif
       case (4)
          if (t%type .eq. 'i') then
             call printsum(f,t%name,id,t%i4)
          else
             call printsum(f,t%name,id,t%r4)
          endif
       end select
    enddo
#ifdef USEMPI
    call printsum(f,'s%is',id,s%is) 
    call printsum(f,'s%js',id,s%js) 
    call printsum(f,'s%lm',id,s%lm) 
    call printsum(f,'s%ln',id,s%ln) 
#endif
  end subroutine printssums

  subroutine printssumso(s)
    use spaceparams
    implicit none
    type (spacepars), intent(in) :: s

    write(*,*)'Start of printssumso'

    write(*,*)'s%xz',sum(s%xz)
    write(*,*)'s%yz',sum(s%yz)
    write(*,*)'s%zs',sum(s%zs)
    write(*,*)'s%u',sum(s%u)
    write(*,*)'s%v',sum(s%v)
    write(*,*)'s%ue',sum(s%ue)
    write(*,*)'s%ve',sum(s%ve)
    write(*,*)'s%H',sum(s%H)
    write(*,*)'s%urms',sum(s%urms)
    write(*,*)'s%zb',sum(s%zb)
    write(*,*)'s%hh',sum(s%hh)
    write(*,*)'s%Fx',sum(s%Fx)
    write(*,*)'s%Fy',sum(s%Fy)
    write(*,*)'s%E',sum(s%E)
    write(*,*)'s%R',sum(s%R)
    write(*,*)'s%D',sum(s%D)
  end subroutine printssumso
#endif
end module debugging_module
