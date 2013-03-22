!
! NAME
!    linear_interp
! SYNOPSIS
!    Interpolate linearly into an array Y given the value of X.
!    Return both the interpolated value and the position in the
!    array.
!
! ARGUMENTS
!    - X: independent array (sorted in ascending order)
!    - Y: array whose values are a function of X
!    - XX: specified value, interpolating in first array
!    - YY: interpolation result
!    - INDINT: (optional) index in array (x(indint) <= xx <= x(indint+1))
!
! SOURCE
!
module interp
contains
  PURE SUBROUTINE linear_interp_2d(X,nx,Y,ny,Z,xx,yy,zz,method,exception)

    implicit none
    ! input/output
    integer,intent(in)                    :: nx,ny
    real*8,dimension(nx),intent(in)       :: X
    real*8,dimension(ny),intent(in)       :: Y
    real*8,dimension(nx,ny),intent(in)    :: Z
    real*8,intent(in)                     :: xx,yy
    real*8,intent(out)                    :: zz
    character(len=*),intent(in)           :: method
    real*8,intent(in)                     :: exception
    ! internal
    integer                               :: i,j,dum
    integer,dimension(4)                  :: ind
    real*8,dimension(2)                   :: yint
    real*8                                :: modx,mody,disx,disy
    logical                               :: interpX,interpY

    ! does the interpolation point fall within the data?
    if (xx>=minval(X) .and. xx<=maxval(X)) then
       interpX = .true.
    else
       interpX = .false.
    endif
    if (yy>=minval(Y) .and. yy<=maxval(Y)) then
       interpY = .true.
    else
       interpY = .false.
    endif

    if (interpX .and. interpY) then
       ! find rank position xx in X direction
       ind(1) = minloc(X,1,X>=xx)
       ind(2) = maxloc(X,1,X<=xx)
       ! find rank position yy in Y direction
       ind(3) = minloc(Y,1,Y>=yy)
       ind(4) = maxloc(Y,1,Y<=yy)
       ! distance between X points and Y points
       disx = X(ind(1))-X(ind(2))
       disy = Y(ind(3))-Y(ind(4))
       ! relative position of (xx,yy) on disx,disy
       if (disx>0.d0) then
          modx = (xx-X(ind(2)))/disx
       else
          modx = 0.d0   ! xx corresponds exactly to X point
       endif
       if (disy>0.d0) then 
          mody = (yy-Y(ind(4)))/disy
       else
          mody = 0.d0   ! yy corresponds exactly to Y point
       endif
       ! interpolate the correct Y value, based on two nearest X intersects
       ! if disy==0 then only single interpolation needed. This could also be done for X, 
       ! but since this is used by waveparams and the interp angles (Y) are more often equal
       ! this is probably faster.
       if (disy>0.d0) then
          yint(1) = (1.d0-modx)*Z(ind(2),ind(3))+modx*Z(ind(1),ind(3))
          yint(2) = (1.d0-modx)*Z(ind(2),ind(4))+modx*Z(ind(1),ind(4))
          zz = (1.d0-mody)*yint(2)+mody*yint(1)
       else
          zz = (1.d0-modx)*Z(ind(2),ind(3))+modx*Z(ind(1),ind(3))        
       endif
    else
       select case (method)
       case ('interp')
          zz = exception
       case ('extendclosest')
          ! find closest X point
          ind(1) = minloc(abs(X-xx),1)
          ! find closest Y point
          ind(3) = minloc(abs(Y-yy),1)
          ! external value given that of closest point
          zz = Z(ind(1),ind(3))
       case default
          zz = exception
       endselect
    endif

  END SUBROUTINE linear_interp_2d


  PURE SUBROUTINE LINEAR_INTERP(X, Y, N, XX, YY, INDINT)
    integer,              intent(in) :: n
    real*8, dimension(n), intent(in) :: x
    real*8, dimension(n), intent(in) :: y
    real*8, intent(in)               :: xx
    real*8, intent(out)              :: yy
    integer, intent(out), optional   :: indint
    !****
    !
    ! CODE: linear interpolation
    !
    !
    !
    real*8           :: a,  b, dyy
    integer          :: j

    yy = 0.0d0
    IF ( present(indint) ) then
       indint = 0
    ENDIF

    IF (N.LE.0) return
    !
    ! *** N GREATER THAN 0
    !
    IF (N.EQ.1) THEN
       YY = Y(1)
       return
    ENDIF

    CALL binary_search( x, N, xx, j )

    IF ( j .LE. 0 ) THEN
       yy = y(1)
    ELSEIF ( j .GE. n ) THEN
       yy = y(n)
    ELSE
       a = x (j+1)
       b = x (j)
       IF ( a .eq. b ) THEN
          dyy = 0.0d0
       ELSE
          dyy = (y(j+1) - y(j)) / (a - b)
       ENDIF
       yy = y(j) + (xx - x(j)) * dyy
    ENDIF

    IF ( present(indint) ) then
       indint = j
    ENDIF

    RETURN

  END SUBROUTINE LINEAR_INTERP

  !****f* Interpolation/binary_search
  !
  ! NAME
  !    binary_search
  ! SYNOPSIS
  !    Perform a binary search in an ordered real array
  !    to find the largest entry equal or lower than a given value:
  !    Given an array XX of length N, given value X, return a value J
  !    such that X is between XX(J) en XX (J+1)
  !
  !    XX must be monotonic, either decreasing or increasing
  !    J=0 or J=N indicates X is out of range.
  !
  ! ARGUMENTS
  !    - XX: ordered array of values
  !    - X: value to be found
  !    - J: index such that X is between XX(J) and XX(J+1)
  !
  ! SOURCE
  !
  PURE SUBROUTINE BINARY_SEARCH(XX, N, X, J) 
    integer,              intent(in) :: N
    real*8, dimension(N), intent(in) :: xx
    real*8, intent(in)               :: x
    integer, intent(out)             :: j
    !****
    !
    ! CODE: binary search in (real) arrays
    !
    ! Requirement:
    !    Parameter wp set to the proper kind
    !
    ! Subroutine from 'Numerical recipes' Fortran  edition.
    ! Given an array XX of length N, given value X, return a value J
    ! such that X is between XX(J) en XX (J+1)
    ! XX must be monotonic, either decreasing or increasin
    ! J=0 or J=N indicates X is out of range.


    !
    ! Local variables
    !
    Integer   jl, ju, jm
    LOGICAL   L1, L2

    JL = 0
    JU = N+1
10  IF (JU-JL .GT. 1) THEN
       JM = (JU+JL)/2
       L1 = XX(N) .GT. XX(1)
       L2 = X .GT. XX(JM)
       IF ( (L1.AND.L2) .OR. (.NOT. (L1 .OR. L2)) ) THEN
          JL = JM
       ELSE
          JU = JM
       ENDIF
       GOTO 10
    ENDIF

    J = JL

    RETURN

  END SUBROUTINE BINARY_SEARCH
  
  subroutine mkmap(code      ,x1        ,y1        ,m1        ,n1        , &
                   & x2        ,y2        ,n2        ,xs        ,ys        , &
                   & nrx       ,nry       ,iflag     ,nrin      ,w         , &
                   & iref      ,iprint    ,covered   ,xymiss)
    !!--description-----------------------------------------------------------------
    !
    !
    !     SUBROUTINE MKMAP
    !     Interpolation of curvilinear, numerically ordered grid (grid1) 
    !     to random points, with weighting points (grid2). 
    !
    !     J.A. Roelvink
    !     Deltares
    !     24-2-1992 (MKMAP
    !
    !     Given: numerically ordered grid M1*N1
    !     with coordinates X1 (1:M1,1:N1)
    !                  and Y1 (1:M1,1:N1)
    !
    !     Also given: random points X2(1:N2)
    !                           and Y2(1:N2)
    !
    !     To be determined:  weighting factors and pointers for bilinear interpolation
    !     Weighting factors and locations of the requested (random) points in the
    !     ordered grid are saved in resp.
    !     W(1:4,1:N2) and Iref(1:4,1:N2)
    !
    !!--pseudo code and references--------------------------------------------------
    ! NONE
    !!--declarations----------------------------------------------------------------
        !
        implicit none
    !
    ! Global variables
    !
        integer                    , intent(in)  :: iprint
        integer                    , intent(in)  :: m1
        integer                    , intent(in)  :: n1
        integer                                  :: n2
        integer , dimension(4 , n2), intent(out) :: iref
        integer , dimension(m1, n1), intent(in)  :: code
        integer , dimension(n2)                  :: covered !  0: target point is not covered by source grid (default)
                                                            !  1: target point is covered by   valid points of source grid
                                                            ! -1: target point is covered by invalid points of source grid
        integer , dimension(n2)                  :: iflag
        integer , dimension(n2)                  :: nrin
        integer , dimension(n2)                  :: nrx
        integer , dimension(n2)                  :: nry
        real                       , intent(in)  :: xymiss  ! missing value in grid 1
        real*8  , dimension(4 , n2)              :: w       ! Remains single precision!
        real*8  , dimension(m1, n1), intent(in)  :: x1
        real*8  , dimension(m1, n1), intent(in)  :: y1
        real*8  , dimension(n2)                  :: x2
        real*8  , dimension(n2)                  :: xs
        real*8  , dimension(n2)                  :: y2
        real*8  , dimension(n2)                  :: ys
    !
    ! Local variables
    !
        integer                           :: i
        integer                           :: i1
        integer                           :: i2
        integer                           :: ier
        integer                           :: iin
        integer                           :: inout
        integer                           :: ip
        integer                           :: ipt
        integer                           :: j1
        integer                           :: lomaxx
        integer                           :: lominx
        integer                           :: lomnx
        integer                           :: lomny
        integer                           :: m1max
        integer                           :: m1min
        integer                           :: n1max
        integer                           :: n1min
        integer                           :: nin
        real*8                           :: eps
        real*8                            :: xpmax
        real*8                            :: xpmean
        real*8                            :: xpmin
        real*8                            :: ypmax
        real*8                            :: ypmean
        real*8                            :: ypmin
        real*8  , dimension(5)            :: xp
        real*8  , dimension(5)            :: yp
        real*8  , dimension(4)            :: hbuff
    !
    !! executable statements -------------------------------------------------------
    !
        eps = 0.00001
        !
        !     Initialise tables
        !
        if (iprint==1) write (*, *) 'in mkmap', m1, n1, n2
        lomnx = 1
        lomny = 1
        m1min = 1
        m1max = m1
        n1min = 1
        n1max = n1
        nrx=0
        nry=0
        iflag=0
        nrin=0
        xs=0.d0
        ys=0.d0
        iref=0
        w=0.d0
        !
        ! Sort X2 en Y2
        !
        call sort(n2        ,x2        ,xs        ,nrx       )
        call sort(n2        ,y2        ,ys        ,nry       )
        !
        ! Loop over all cels of grid1
        !
        !  
        do j1 = n1min, n1max - 1
           do i1 = m1min, m1max - 1
              !
              ! Cell definition
              !
              xp(1) = x1(i1, j1)
              xp(2) = x1(i1 + 1, j1)
              xp(3) = x1(i1 + 1, j1 + 1)
              xp(4) = x1(i1, j1 + 1)
              yp(1) = y1(i1, j1)
              yp(2) = y1(i1 + 1, j1)
              yp(3) = y1(i1 + 1, j1 + 1)
              yp(4) = y1(i1, j1 + 1)
              if (     (xp(1) > xymiss-eps .and. xp(1) < xymiss+eps) &
                & .or. (xp(2) > xymiss-eps .and. xp(2) < xymiss+eps) &
                & .or. (xp(3) > xymiss-eps .and. xp(3) < xymiss+eps) &
                & .or. (xp(4) > xymiss-eps .and. xp(4) < xymiss+eps) &
                & .or. (yp(1) > xymiss-eps .and. yp(1) < xymiss+eps) &
                & .or. (yp(2) > xymiss-eps .and. yp(2) < xymiss+eps) &
                & .or. (yp(3) > xymiss-eps .and. yp(3) < xymiss+eps) &
                & .or. (yp(4) > xymiss-eps .and. yp(4) < xymiss+eps) ) cycle
              !
              ! Determine minimum and maximum X and Y of the cell
              !
              xpmin =  1.e10
              xpmax = -1.e10
              ypmin =  1.e10
              ypmax = -1.e10
              do ip = 1, 4
                 xpmin = min(xp(ip), xpmin)
                 xpmax = max(xp(ip), xpmax)
                 ypmin = min(yp(ip), ypmin)
                 ypmax = max(yp(ip), ypmax)
              enddo
              xpmean = .5*(xpmin + xpmax)
              ypmean = .5*(ypmin + ypmax)
              !
              ! First selection of points of grid2 that can be located in the cell
              !
              ! Find centre of the cell in tables Xs and Ys
              !
              call hunt(xs        ,n2        ,xpmean    ,lomnx     )
              call hunt(ys        ,n2        ,ypmean    ,lomny     )
              !
              ! For points with X-values between Xpmin and Xpmax set: iFlag(i)=1
              !
              lominx = lomnx
              lomaxx = lomnx
              do i = lomnx, 1, -1
                 if (xs(i)>=xpmin) then
                    lominx        = i
                    iflag(nrx(i)) = 1
                 else
                    exit
                 endif
              enddo
              do i = lomnx + 1, n2
                 if (xs(i)<=xpmax) then
                    lomaxx        = i
                    iflag(nrx(i)) = 1
                 else
                    exit
                 endif
              enddo
              !
              ! For the points with Y-values between Ypmin and Ypmax,
              ! that also lie between Xpmin and Xpmax: Save them in NrIn
              !
              iin = 1
              do i = lomny, 1, -1
                 if (ys(i)>=ypmin) then
                    nrin(iin) = nry(i)*iflag(nry(i))
                    iin       = iin + iflag(nry(i))
                 else
                    exit
                 endif
              enddo
              do i = lomny + 1, n2
                 if (ys(i)<=ypmax) then
                    nrin(iin) = nry(i)*iflag(nry(i))
                    iin       = iin + iflag(nry(i))
                 else
                    exit
                 endif
              enddo
              nin = iin - 1
              !
              ! Put iFlag back to 0
              !
              do i = lominx, lomaxx
                 if (i/=0) iflag(nrx(i)) = 0
              enddo
              !
              ! Check whether selected points of grid2 lie within the cell
              ! using subroutine IPON; if so, determine weights W of the surrounding
              ! values in grid1 using subroutine INTRP4. Save the weights in Wtab
              ! The reference to grid1 is saved in arrays Iref and Jref.
              !
              do iin = 1, nin
                 i2 = nrin(iin)
                 inout = -1
                 call ipon(xp        ,yp        ,4         ,x2(i2)    , &
                &          y2(i2)    ,inout     )
                 if (inout>=0) then
                    !
                    ! Check on point with nonsense information on source grid 1;
                    ! This depends on agreements
                    !
                    if (      (code(i1    , j1    )==1 .or. code(i1    , j1    )==3) &
                      & .and. (code(i1 + 1, j1    )==1 .or. code(i1 + 1, j1    )==3) &
                      & .and. (code(i1 + 1, j1 + 1)==1 .or. code(i1 + 1, j1 + 1)==3) &
                      & .and. (code(i1    , j1 + 1)==1 .or. code(i1    , j1 + 1)==3) ) then
                       !
                       ! grid2 point is covered by 4 valid grid1 points
                       !
                       covered(i2) = 1
                       call bilin5(xp, yp, x2(i2), y2(i2), hbuff, ier)
                       w(:, i2) = real(hbuff(:),4)
                       !
                       iref(1, i2) = i1 + (j1 - 1)*m1
                       iref(2, i2) = i1 + 1 + (j1 - 1)*m1
                       iref(3, i2) = i1 + 1 + j1*m1
                       iref(4, i2) = i1 + j1*m1
                    elseif (code(i1, j1)>0 .and. code(i1 + 1, j1)>0 .and. &
                      & code(i1 + 1, j1 + 1)>0 .and. code(i1, j1 + 1)>0) then
                       !
                       ! Grid2 point is covered by 4 valid grid1 points, containing
                       ! one or more boundary points:
                       ! - do not use grid1 values
                       ! - use extrapolation values
                       !
                       covered(i2) = 0
                    else
                       !
                       ! Grid2 point is not covered by 4 valid grid1 points:
                       ! - do not use grid1 values
                       ! - do not use extrapolation values
                       !
                       covered(i2) = -1
                    endif
                 endif
              enddo
           enddo
        enddo
    end subroutine mkmap
  
  subroutine grmap(f1        ,n1        ,f2        ,n2        ,iref      , &
                 & w         ,np        ,iprint    )
  !!--description-----------------------------------------------------------------
  !
  ! compute interpolated values for all points on grid 2
  !
  ! special treatment of points on grid 2 that are outside
  ! grid 1; in that case iref(1,i2)=0 AND w(ip,i2)=0 for all ip
  !
  ! Iref(1,i2)   i1    ifac   F2(i2)*ifac     Result
  !
  !      0        1      1      F2(i2)        Old value is kept
  !    j,j>0      j      0       0.           F2 is initialized
  !
  !!--pseudo code and references--------------------------------------------------
  ! NONE
  !!--declarations----------------------------------------------------------------
      implicit none
  !
  ! Global variables
  !
      integer                   , intent(in)  :: iprint
      integer                   , intent(in)  :: n1
      integer                   , intent(in)  :: n2
      integer                   , intent(in)  :: np
      integer, dimension(np, n2), intent(in)  :: iref
      real*8 , dimension(n1)    , intent(in)  :: f1
      real*8 , dimension(n2)                  :: f2
      real*8 , dimension(np, n2), intent(in)  :: w
  !
  ! Local variables
  !
      integer :: i
      integer :: i1
      integer :: i2
      integer :: ifac
      integer :: ip
  !
  !! executable statements -------------------------------------------------------
  !
      if (iprint==1) write (*, *) 'in grmap n1 n2', n1, n2
      do i2 = 1, n2
         i = iref(1, i2)
         i1 = max(i, 1)
         ifac = 1 - i/i1
         f2(i2) = f2(i2)*ifac
         !
         ! Function values at grid 2 are expressed as weighted average
         ! of function values in Np surrounding points of grid 1
         !
         if (iprint==1 .and. i2<=n2) &
          & write (*, '(1X,A,I6,4(1X,E10.4))') ' i2 w ', i2, (w(ip, i2) , ip = 1,  &
          & np)
         do ip = 1, np
            i = iref(ip, i2)
            i1 = max(i, 1)
            if (iprint==1 .and. i2<=n2) write (*, *) ' i1,f1(i1) ', i1, f1(i1)
            f2(i2) = f2(i2) + w(ip, i2)*f1(i1)
         enddo
      enddo
  end subroutine grmap
 
 subroutine grmap2(f1, cellsz1i     , n1        ,f2  ,cellsz2,  n2        ,iref      , &
                 & w         ,np     )
  !!--description-----------------------------------------------------------------
  !
  ! compute interpolated values for all points on grid 1 given reference table
  ! for grid 2; this works the other way round from GRMAP. Assumption is that
  ! grid 2 is much finer than grid 1. For each point in grid 2 we know the
  ! surrounding points in grid 1 and the related weights. Instead of using this to 
  ! interpolate from 1 to 2 we now integrate from 2 to 1 using the same weights.
  !
   !
  !!--pseudo code and references--------------------------------------------------
  ! NONE
  !!--declarations----------------------------------------------------------------
      implicit none
  !
  ! Global variables
  !
      integer                   , intent(in)  :: n1
      integer                   , intent(in)  :: n2
      integer                   , intent(in)  :: np
      integer, dimension(np, n2), intent(in)  :: iref
      real*8 , dimension(n1)                  :: f1
      real*8 , dimension(n1)    , intent(in)  :: cellsz1i   !array with 1/cell size
      real*8                    , intent(in)  :: cellsz2
      real*8 , dimension(n2)    , intent(in)  :: f2
      real*8 , dimension(np, n2), intent(in)  :: w
  !
  ! Local variables
  !
      integer :: i
      integer :: i1
      integer :: i2
      integer :: ifac
      integer :: ip
  !
  !! executable statements -------------------------------------------------------
  !
         do ip = 1, np
            do i2=1,n2
               i1 = iref(ip, i2)
               if (i1>0) then
                  !f2(i2) = f2(i2) + w(ip, i2)*f1(i1)
                   f1(i1) = f1(i1) + w(ip, i2)*f2(i2)*cellsz2*cellsz1i(i1)
               endif
            enddo
         enddo
  end subroutine grmap2
  
  subroutine ipon(xq     ,yq     ,n      ,xp     ,yp     ,inout     )
  !--description----------------------------------------------------------------
  !
  ! Deltares                                                               *
  ! AUTHOR : J.A.ROELVINK                                                  *
  ! DATE   : 22-12-1988                                                    *
  ! DETERMINE WHETHER POINT (xp,yp) LIES IN POLYGON (x,y) OF n POINTS      *
  ! POINT n+1 IS SET EQUAL TO POINT 1                                      *
  ! (ARRAY MUST HAVE DIMENSION n+1 IN MAIN PROGRAMME                       *
  ! inpout = -1 :  OUTSIDE POLYGON                                         *
  ! inpout =  0 :  ON EDGE OF POLYGON                                      *
  ! inpout =  1 :  INSIDE POLYGON                                          *
  ! USED METHOD :         - DRAW A VERTICAL LINE THROUGH (xp,yp)           *
  !                       - DETERMINE NUMBER OF INTERSECTIONS WITH POLYGON *
  !                         UNDER yp : nunder                              *
  !                       - IF nunder IS EVEN, THEN THE POINT LIES OUTSIDE *
  !                         THE POLYGON, OTHERWISE IT LIES INSIDE          *
  !                       - THE EDGE IS TREATED SEPARATELY                 *
  !
  !--pseudo code and references-------------------------------------------------
  ! NONE
  !--declarations---------------------------------------------------------------
      !
      implicit none
  !
  ! Global variables
  !
      integer               , intent(out) :: inout
      integer               , intent(in)  :: n
      real*8                , intent(in)  :: xp
      real*8                , intent(in)  :: yp
      real*8  , dimension(*)              :: xq
      real*8  , dimension(*)              :: yq
  !
  ! Local variables
  !
      integer                             :: i
      integer                             :: ierr
      integer                             :: nunder
      real*4                              :: ysn
      real*4  , dimension(:), allocatable :: x
      real*4  , dimension(:), allocatable :: y
  !
  ! executable statements ------------------------------------------------------
  !
      allocate(x(n+1))
      allocate(y(n+1))
      do i = 1, n
         x(i) = real( xq(i)-xp , 4)
         y(i) = real( yq(i)-yp , 4)
      enddo
      x(n + 1) = x(1)
      y(n + 1) = y(1)
      nunder   = 0
      do i = 1, n
         if ((x(i)<0. .and. x(i + 1)>=0.).or.(x(i + 1)<0. .and. x(i)>=0.)) then
            if (y(i)<0. .and. y(i + 1)<0.) then
               nunder = nunder + 1
            elseif ((y(i)<=0. .and. y(i + 1)>=0.) .or.                       &
                  & (y(i + 1)<=0. .and. y(i)>=0.)) then
               ysn = (y(i)*x(i + 1) - x(i)*y(i + 1))/(x(i + 1) - x(i))
               if (ysn<0.) then
                  nunder = nunder + 1
               elseif (ysn<=0.) then
                  !
                  ! Edge
                  !
                  inout = 0
                  goto 100
               else
               endif
            else
            endif
         elseif (abs(x(i))<1.0E-8 .and. abs(x(i + 1))<1.0E-8) then
            if ((y(i)<=0. .and. y(i + 1)>=0.).or.(y(i + 1)<=0..and.y(i)>=0.)) &
              & then
               !
               ! Edge
               !
               inout = 0
               goto 100
            endif
         else
         endif
      enddo
      if (mod(nunder, 2)==0) then
         !
         ! Outside
         !
         inout = -1
      else
         !
         ! Inside
         !
         inout = 1
      endif
    100 continue
    deallocate(x, stat=ierr)
    deallocate(y, stat=ierr)
  end subroutine ipon
  
  subroutine hunt(xx        ,n         ,x         ,jlo       )
  !!--description-----------------------------------------------------------------
  ! NONE
  !!--pseudo code and references--------------------------------------------------
  ! NONE
  !!--declarations----------------------------------------------------------------
      !
      implicit none
  !
  ! Global variables
  !
      integer                                  :: jlo
      integer                    , intent(in)  :: n
      real*8              , intent(in)  :: x
      real*8, dimension(n), intent(in)  :: xx
  !
  ! Local variables
  !
      integer :: inc
      integer :: jhi
      integer :: jm
      logical :: ascnd
  !
  !! executable statements -------------------------------------------------------
  !
      ascnd = xx(n)>=xx(1)
      if (jlo<=0 .or. jlo>n) then
         jlo = 0
         jhi = n + 1
         goto 3
      endif
      inc = 1
      if (x>=xx(jlo) .eqv. ascnd) then
      1  continue
         jhi = jlo + inc
         if (jhi>n) then
            jhi = n + 1
         elseif (x>=xx(jhi) .eqv. ascnd) then
            jlo = jhi
            inc = inc + inc
            goto 1
         else
         endif
      else
         jhi = jlo
      2  continue
         jlo = jhi - inc
         if (jlo<1) then
            jlo = 0
         elseif (x<xx(jlo) .eqv. ascnd) then
            jhi = jlo
            inc = inc + inc
            goto 2
         else
         endif
      endif
      3 continue
      if (jhi - jlo==1) then
         return
      endif
      jm = (jhi + jlo)/2
      if (x>xx(jm) .eqv. ascnd) then
         jlo = jm
      else
         jhi = jm
      endif
      goto 3
  end subroutine hunt
  
  subroutine indexx(n         ,arrin     ,indx      )
  !----- GPL ---------------------------------------------------------------------
  !                                                                               
  !  Copyright (C)  Stichting Deltares, 2011.                                     
  !                                                                               
  !  This program is free software: you can redistribute it and/or modify         
  !  it under the terms of the GNU General Public License as published by         
  !  the Free Software Foundation version 3.                                      
  !                                                                               
  !  This program is distributed in the hope that it will be useful,              
  !  but WITHOUT ANY WARRANTY; without even the implied warranty of               
  !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
  !  GNU General Public License for more details.                                 
  !                                                                               
  !  You should have received a copy of the GNU General Public License            
  !  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
  !                                                                               
  !  contact: delft3d.support@deltares.nl                                         
  !  Stichting Deltares                                                           
  !  P.O. Box 177                                                                 
  !  2600 MH Delft, The Netherlands                                               
  !                                                                               
  !  All indications and logos of, and references to, "Delft3D" and "Deltares"    
  !  are registered trademarks of Stichting Deltares, and remain the property of  
  !  Stichting Deltares. All rights reserved.                                     
  !                                                                               
  !-------------------------------------------------------------------------------
  !  $Id$
  !  $HeadURL$
  !!--description-----------------------------------------------------------------
  ! NONE
  !!--pseudo code and references--------------------------------------------------
  ! NONE
  !!--declarations----------------------------------------------------------------
      !
      implicit none
  !
  ! Global variables
  !
      integer                       , intent(in)  :: n
      integer, dimension(n)                       :: indx
      real*8   , dimension(n), intent(in)  :: arrin
  !
  ! Local variables
  !
      integer :: i
      integer :: indxt
      integer :: ir
      integer :: j
      integer :: l
      real*8 :: q
  !
  !! executable statements -------------------------------------------------------
  !
      do j = 1, n
         indx(j) = j
      enddo
      l = n/2 + 1
      ir = n
     10 continue
      if (l>1) then
         l = l - 1
         indxt = indx(l)
         q = arrin(indxt)
      else
         indxt = indx(ir)
         q = arrin(indxt)
         indx(ir) = indx(1)
         ir = ir - 1
         if (ir==1) then
            indx(1) = indxt
            return
         endif
      endif
      i = l
      j = l + l
     20 continue
      if (j<=ir) then
         if (j<ir) then
            if (arrin(indx(j))<arrin(indx(j + 1))) j = j + 1
         endif
         if (q<arrin(indx(j))) then
            indx(i) = indx(j)
            i = j
            j = j + j
         else
            j = ir + 1
         endif
         goto 20
      endif
      indx(i) = indxt
      goto 10
  end subroutine indexx
  
  subroutine bilin5(xa, ya, x0, y0, w, ier)
  !!--description-----------------------------------------------------------------
  ! NONE
  !!--pseudo code and references--------------------------------------------------
  !
  ! Author: H. Petit
  !
  !!--declarations----------------------------------------------------------------
      !
      implicit none
  !
  ! Global variables
  !
      integer               , intent(out) :: ier
      real*8                , intent(in)  :: x0
      real*8                , intent(in)  :: y0
      real*8  , dimension(4), intent(out) :: w
      real*8  , dimension(4), intent(in)  :: xa
      real*8  , dimension(4), intent(in)  :: ya
  !
  ! Local variables
  !
      real*8   :: a
      real*8   :: a21
      real*8   :: a22
      real*8   :: a31
      real*8   :: a32
      real*8   :: a41
      real*8   :: a42
      real*8   :: b
      real*8   :: c
      real*8   :: det
      real*8   :: discr
      real*8   :: eta
      real*8   :: x
      real*8   :: x1
      real*8   :: x2
      real*8   :: x3
      real*8   :: x3t
      real*8   :: x4
      real*8   :: xi
      real*8   :: xt
      real*8   :: y
      real*8   :: y1
      real*8   :: y2
      real*8   :: y3
      real*8   :: y3t
      real*8   :: y4
      real*8   :: yt
  !
  !! executable statements -------------------------------------------------------
  !
      ! read(12,*)x1,y1,f1
      x1 = xa(1)
      y1 = ya(1)
      ! read(12,*)x2,y2,f2
      x2 = xa(2)
      y2 = ya(2)
      ! read(12,*)x3,y3,f3
      x3 = xa(3)
      y3 = ya(3)
      ! read(12,*)x4,y4,f4
      x4 = xa(4)
      y4 = ya(4)
      x  = x0
      y  = y0
      !
      ! The bilinear interpolation problem is first transformed
      ! to the quadrangle with nodes (0,0),(1,0),(x3t,y3t),(0,1)
      ! and required location (xt,yt)
      !
      a21 = x2 - x1
      a22 = y2 - y1
      a31 = x3 - x1
      a32 = y3 - y1
      a41 = x4 - x1
      a42 = y4 - y1
      det = a21*a42 - a22*a41
      if (abs(det) < 1.0e-20) then
         ier = 1
         goto 99999
      endif
      x3t = (  a42*a31      - a41*a32     ) / det
      y3t = ( -a22*a31      + a21*a32     ) / det
      xt  = (  a42*(x - x1) - a41*(y - y1)) / det
      yt  = ( -a22*(x - x1) + a21*(y - y1)) / det
      if ((x3t < 0.0) .or. (y3t < 0.0)) then
         ! write (*, *) 'distorted quadrangle'
         ier = 1
         goto 99999
      endif
      if (abs(x3t - 1.0d0) < 1.0e-7) then
         xi = xt
         if (abs(y3t - 1.0d0) < 1.0e-7) then
            eta = yt
         elseif (abs(1.0d0 + (y3t - 1.0d0)*xt) < 1.0e-6) then
            ! write (*, *) 'extrapolation over too large a distance'
            ier = 1
            goto 99999
         else
            eta = yt / (1.0d0 + (y3t - 1.0d0)*xt)
         endif
      elseif (abs(y3t - 1.0d0) < 1.0e-6) then
         eta = yt
         if (abs(1.0d0 + (x3t - 1.0d0)*yt) < 1.0e-6) then
            ! write (*, *) 'extrapolation over too large a distance'
            ier = 1
            goto 99999
         else
            xi = xt / (1.0d0 + (x3t - 1.0d0)*yt)
         endif
      else
         a     = y3t - 1.0d0
         b     = 1.0d0 + (x3t - 1.0d0)*yt - (y3t - 1.0d0)*xt
         c     = -xt
         discr = b*b - 4.0d0*a*c
         if (discr < 1.0e-6) then
            ! write (*, *) 'extrapolation over too large a distance'
            ier = 1
            goto 99999
         endif
         xi  = (-b + sqrt(discr)) / (2.0d0*a)
         eta = ((y3t - 1.0d0)*(xi - xt) + (x3t - 1.0d0)*yt) / (x3t - 1.0d0)
      endif
      w(1) = (1.0d0-xi) * (1.0d0-eta)
      w(2) =         xi  * (1.0d0-eta)
      w(3) =         xi  *         eta 
      w(4) =        eta  * (1.0d0-xi ) 
      return
  99999 continue
  end subroutine bilin5
  
  subroutine sort(n         ,ra        ,wksp      ,iwksp     )
  !!--description-----------------------------------------------------------------
  ! Sorts an array, routine from Numerical Recipes
  !!--pseudo code and references--------------------------------------------------
  ! NONE
  !!--declarations----------------------------------------------------------------
      !
      implicit none
  !
  ! Global variables
  !
      integer                                     :: n
      integer, dimension(n)                       :: iwksp
      real*8   , dimension(n)              :: ra
      real*8   , dimension(n), intent(out) :: wksp
  !
  ! Local variables
  !
      integer :: j
  !
  !! executable statements -------------------------------------------------------
  !
      call indexx(n         ,ra        ,iwksp     )
      do j = 1, n
         wksp(j) = ra(iwksp(j))
      enddo
  end subroutine sort
  

end module interp
