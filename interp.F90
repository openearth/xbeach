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
