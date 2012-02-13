module interp_tests

    use ftnunit
    use interp

    implicit none

contains

!!! Run all interp tests 
subroutine allInterpTests

    call test( interpSimpleValue,     "Interpolate simple value" )
    
end subroutine allInterpTests


!!! tests section

! interpolate linear (1d)
subroutine interpSimpleValue
    
    real*8, dimension(2)             :: x, y
    real*8                           :: yy, xx

    xx = 5.0d0
    x(1) = 0
    x(2) = 10
    y(1) = 0
    y(2) = 1

    call linear_interp(x, y, 2, xx, yy)
     
    call assert_true( yy == 0.5, "Interpolated value is 0.5" )

end subroutine interpSimpleValue


end module interp_tests