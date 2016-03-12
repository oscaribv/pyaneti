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
!This subrouotine finds the time of periastron passage
!by knowing the transit time
!------------------------------------------------------------
subroutine find_tp(t0, e, w, P, tp)
implicit none
!In/Out variables
  double precision, intent(in) :: t0, e, w, P
  double precision, intent(out) :: tp
!Local variables
  double precision :: theta_p
  double precision, parameter :: pi = 3.1415926535897932384626

  theta_p = atan2( sqrt(1.-e*e) * sin( pi/2. - w ), e + cos( pi/2. - w) )
  theta_p = theta_p - e * sin(theta_p)

  tp = t0 - theta_p * p / 2. / pi

end subroutine

!------------------------------------------------------------
!This subroutine finds the true anomaly of an eccentric orbit
!by using the Newton-Raphson (NR)  algorithm 
!The input parameters are:
! man -> mean anomaly, ec -> eccentricity, delta -> NR limit
! imax -> iteration limit for NR, dman -> man dimension
!The output parameters are:
! ta -> True anomaly (vector with the same dimension that man)
!------------------------------------------------------------
subroutine find_anomaly(t,t0,e,w,P,ta,delta,imax,dt)
implicit none
!In/Out variables
  integer, intent(in) :: dt
  double precision, intent(in) , dimension(0:dt-1) :: t
  double precision, intent(out), dimension(0:dt-1) :: ta
  double precision, intent(in) :: t0, e, w, P, delta
  integer, intent(in) :: imax
!Local variables
  integer :: i,n
  double precision, dimension(0:dt-1) :: ma, f, df, esin, ecos
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision :: uno, tp
!
  uno = dble(1.)

  call find_tp(t0,e,w,P,tp)

  !Calculate the mean anomaly
  ma = 2. * pi * ( t - tp ) / P

  !calculate the eccentric anomaly
  ta(:) = ma(:)
  f(:) = delta * 10
  n = 0

  do i = 0, dt-1
    do while ( abs(f(i)) >= delta .and. n <= imax )
      f(i)   = ta(i) - e * sin(ta(i)) - ma(i)
      df(i)  =   uno - e * cos(ta(i))
      ta(i)  = ta(i) - f(i) / df(i)
      n = n + 1
    end do
  end do

  if ( n > imax ) then
    print *, 'I am tired, too much Newton-Raphson for me!'
    stop
  end if 

  !ta(:) = sqrt(1. + e )  * sin(ta(:)*0.5) / sqrt(1. - e )
  !ta(:) = 2. * atan2( ta(:), cos(ta(:)*0.5) )

  !calculate the true anomaly
  esin = ( sqrt(uno-e*e) * sin(ta) ) / (uno-e*cos(ta))
  ecos = ( cos (ta) - e ) / (uno-e*cos(ta))
  ta = atan2(esin,ecos)
  !do i = 0, dma-1
  !  if ( ta(i) < dble(0.0) ) ta(i) = ta(i) + 2. * pi 
  !end do
  !ta(:) = sqrt( ( uno + e ) / ( uno - e ) ) * tan(ta(:)*0.5) 
  !ta(:) = 2.0 * asin( ta(:) / (sqrt( uno + ta(:)*ta(:) ) ) )

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


!Subroutine to get g(z), equation(10)
!emcee paper -> http://arxiv.org/abs/1202.3665
subroutine find_gz(z,a)
implicit none

!In/Out variables
  double precision, intent(inout) :: z
  !f2py itent(in,out) :: z
  double precision, intent(in) :: a

  if ( z >= 1./ a .and. z <= a ) then
    z = 1. / sqrt(z)
  else
    z = 0.0
  end if


end subroutine


!Subroutine to create random integers between 0 and n
subroutine random_int(r_int,n)
implicit none

  !In/Out variables
  integer, intent(in) :: n
  integer, intent(out), dimension(0:n-1) :: r_int
  !Local variables
  real :: r_real
  integer :: i

  do i = 0, n - 1
    r_int(i) = i
    do while ( r_int(i) == i )
    call random_number(r_real)
    r_int = int( r_real * ( n + 1 ) ) - 1
    end do
  end do

end subroutine
