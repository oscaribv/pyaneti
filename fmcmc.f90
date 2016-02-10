!Template to write a fortran90 program
!written by Oscar Barragán, Jan 13 2014

!-------------------------------------------------------------------------
!Name of the program
!Here write a simple description of the program
!Date -->
!-------------------------------------------------------------------------

subroutine find_anomaly(man,ta,ec,delta,imax,dman)
implicit none

integer, intent(in) :: dman
double precision, intent(in) , dimension(0:dman-1) :: man
double precision, intent(out), dimension(0:dman-1) :: ta
double precision, intent(in) :: ec, delta
integer, intent(in) :: imax

integer :: i,j,n
double precision, dimension(0:dman-1) :: f, df

  ta(:)  = 0.0
  f(:)   = ta(:) - ec * sin(ta(:)) - man(:)
  df(:)  =   1.0 - ec * cos(ta(:))
  n = 0

  do i = 0, dman-1
    do while ( abs(f(i)) >= delta .and. n <= imax )
      ta(i)  = ta(i) - f(i) / df(i)
      f(i)   = ta(i) - ec * sin(ta(i)) - man(i)
      df(i)  =   1.0 - ec * cos(ta(i))
      n = n + 1
    end do
  end do

  if ( n > imax ) then
    print *, 'I am tired, so much N-R'
    stop
  end if 

  ta(:) = sqrt( (1. + ec) / (1. -ec) ) * tan(ta(:)*0.5)
  ta(:) = 2. * atan(ta(:))

end subroutine

