!------------------------------------------------------------
!                         todo.f90
! This file has a some generic functions to be called for more
!    general subroutines or programs.
! The subroutines can be called from python by using f2py
! They also can be used in a fortran program
!              Date --> Feb  2016, Oscar Barragán
!------------------------------------------------------------

!----------------------------------------------------------
! This subroutine generates pseudo-random seeds
! Taken from
! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
!----------------------------------------------------------
subroutine init_random_seed()
use iso_fortran_env, only: int64
implicit none
integer, allocatable :: seed(:)
integer :: i, n, un, istat, dt(8), pid
integer(int64) :: t

call random_seed(size = n)
allocate(seed(n))
! First try if the OS provides a random number generator
open(newunit=un, file="/dev/urandom", access="stream", &
     form="unformatted", action="read", status="old", iostat=istat)
if (istat == 0) then
   read(un) seed
   close(un)
else
   ! Fallback to XOR:ing the current time and pid. The PID is
   ! useful in case one launches multiple instances of the same
   ! program in parallel.
   call system_clock(t)
   if (t == 0) then
      call date_and_time(values=dt)
      t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_int64 * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
   end if
   pid = getpid()
   t = ieor(t, int(pid, kind(t)))
   do i = 1, n
      seed(i) = lcg(t)
   end do
end if
call random_seed(put=seed)

contains
! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
function lcg(s)
  integer :: lcg
  integer(int64) :: s
  if (s == 0) then
     s = 104729
  else
     s = mod(s, 4294967296_int64)
  end if
  s = mod(s * 279470273_int64, 4294967291_int64)
  lcg = int(mod(s, int(huge(0), int64)), kind(0))
end function lcg
end subroutine init_random_seed

!------------------------------------------------------------
!This subroutine finds the true anomaly of an eccentric orbit
!by using the Newton-Raphson (NR)  algorithm 
!The input parameters are:
! man -> mean anomaly, ec -> eccentricity, delta -> NR limit
! imax -> iteration limit for NR, dman -> man dimension
!The output parameters are:
! ta -> True anomaly (vector with the same dimension that man)
!------------------------------------------------------------
subroutine find_anomaly(man,ta,ec,delta,imax,dman)
implicit none
!In/Out variables
  integer, intent(in) :: dman
  double precision, intent(in) , dimension(0:dman-1) :: man
  double precision, intent(out), dimension(0:dman-1) :: ta
  double precision, intent(in) :: ec, delta
  integer, intent(in) :: imax
!Local variables
  integer :: i,n
  double precision, dimension(0:dman-1) :: f, df
!
  ta(:)  = 0.0
  f(:)   = ta(:) - ec * sin(ta(:)) - man(:)
  df(:)  =   1.0 - ec * cos(ta(:))
  n = 0

  do i = 0, dman-1
    do while ( abs(f(i)) >= delta .and. n <= imax )
      !This can be improved, do it!
      ta(i)  = ta(i) - f(i) / df(i)
      f(i)   = ta(i) - ec * sin(ta(i)) - man(i)
      df(i)  =   1.0 - ec * cos(ta(i))
      n = n + 1
    end do
  end do

  if ( n > imax ) then
    print *, 'I am tired, too much Newton-Raphson for me!'
    stop
  end if 

  ta(:) = sqrt( (1. + ec) / (1. - ec) ) * tan(ta(:)*0.5)
  ta(:) = 2. * atan(ta(:))

end subroutine

!Subroutine to check if the chi2 minimization is working
subroutine fit_a_line(x_vec,chi2_vec,a,b,svec)
implicit none

!In/Out variables
  integer, intent(in) :: svec
  double precision, intent(in), dimension(0:svec-1) :: chi2_vec, x_vec
  double precision, intent(out) :: a, b

!Local variables
  double precision :: meanx, meany
  double precision, dimension(0:svec-1) :: resx, resy

!Compute the linear fit by the least square method
  meanx = sum(x_vec) / float(svec)
  meany = sum(chi2_vec) / float(svec)

  resx = x_vec - meanx
  resy = chi2_vec - meany

  b = sum( resx * resy )
  b = b / sum ( resx * resx )

  a = meany - b * meanx

end subroutine
