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
     SUBROUTINE linear_interp_2d(X,nx,Y,ny,Z,xx,yy,zz,method,exception)
     
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


     SUBROUTINE LINEAR_INTERP(X, Y, N, XX, YY, INDINT)
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
      SUBROUTINE BINARY_SEARCH(XX, N, X, J) 
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
10      IF (JU-JL .GT. 1) THEN
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

end module interp
