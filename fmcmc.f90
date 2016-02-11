!Template to write a fortran90 program
!written by Oscar Barragán, Jan 13 2014

!-------------------------------------------------------------------------
!Name of the program
!Here write a simple description of the program
!Date -->
!-------------------------------------------------------------------------

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

!-----------------------------------------------------------

subroutine find_anomaly(man,ta,ec,delta,imax,dman)
implicit none

integer, intent(in) :: dman
double precision, intent(in) , dimension(0:dman-1) :: man
double precision, intent(out), dimension(0:dman-1) :: ta
double precision, intent(in) :: ec, delta
integer, intent(in) :: imax

integer :: i,n
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

!-----------------------------------------------------------

subroutine rv_circular(t,rv0,t0,k,P,rv,ts)
implicit none

integer, intent(in) :: ts
double precision, intent(in), dimension(0:ts-1)  :: t
double precision, intent(out), dimension(0:ts-1) :: rv
double precision, intent(in) :: k, rv0, t0, P

double precision, parameter :: pi = 3.1415926535897932384626

  rv(:) = rv0 - k * sin( 2.*pi*( t(:) - t0) / P )

end subroutine

!-----------------------------------------------------------

subroutine rv_curve(t,rv0,t0,k,P,ec,w,rv,ts)
implicit none

integer, intent(in) :: ts
double precision, intent(in), dimension(0:ts-1)  :: t
double precision, intent(out), dimension(0:ts-1) :: rv
double precision, intent(in) :: k, rv0, t0, P, ec, w

double precision, parameter :: pi = 3.1415926535897932384626
double precision, dimension(0:ts-1) :: ma, ta
double precision :: delta = 1e-4
integer :: imax

external :: find_anomaly

  imax = int(1e5)

  ma(:) = 2.* pi * ( t(:) - t0 ) / P
  call find_anomaly(ma,ta,ec,delta,imax,ts)

  rv(:) = rv0 + k * ( cos(ta(:) + w ) + ec * cos(w) )
  
end subroutine

!-----------------------------------------------------------

subroutine find_chi2(xd,yd,errs,rv0,k,ec,w,t0,P,chi2,datas)
implicit none

integer, intent(in) :: datas
double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
double precision, intent(in)  :: k, rv0, t0, P, ec, w
double precision, intent(out) :: chi2

double precision, parameter :: pi = 3.1415926535897932384626
double precision, dimension(0:datas-1) :: model, res

external :: rv_circular, rv_curve

  chi2 = 0.0

  !Here we could control what kind of fit we want to do

  !call rv_circular(xd,rv0,t0,k,P,model,datas)
  call rv_curve(xd,rv0,t0,k,P,ec,w,model,datas)

  res(:) = ( yd(:) - model(:) ) / errs(:)
  res(:) = res(:) * res(:)

  chi2 = sum(res)

end subroutine

!-----------------------------------------------------------

subroutine mcmc_rv(xd,yd,errs,rv0,k,ec,w,t0,P,rv0mc,kmc,ecmc,wmc,t0mc,Pmc,prec,imcmc,datas)
implicit none

integer, intent(in) :: imcmc, datas
double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
double precision, intent(in)  :: k, rv0, t0, P, ec, w,prec
double precision, intent(out), dimension(0:imcmc-1)  :: kmc, rv0mc, t0mc, Pmc, ecmc, wmc

double precision, parameter :: pi = 3.1415926535897932384626
double precision :: chi2_old, chi2_new
double precision  :: sk, srv0, st0, sP, sec, sw
double precision  :: q
integer :: i
real, dimension(0:4) :: r

external :: init_random_seed, find_chi2

  kmc(0)    = k
  rv0mc(0)  = rv0
  t0mc(0)   = t0
  Pmc(0)    = P
  ecmc(0)   = ec
  wmc(0)    = w

  sk = k * prec
  srv0 = k * prec
  st0 = t0 * prec
  sP = P * prec
  sec = 1. * prec
  sw = 2.*pi * prec

  call find_chi2(xd,yd,errs,rv0mc(0),kmc(0),ecmc(0),wmc(0),t0mc(0),Pmc(0),chi2_old,datas)

  call init_random_seed()
  !Let us create an array of random numbers to save time
  !Check the ram consumption

  do i = 1, imcmc - 1
    call random_number(r)
    r(0:3) = ( r(0:3) - 0.5) * 2.
    rv0mc(i) =   rv0mc(i-1) + r(0) * srv0
    kmc(i)   =   kmc(i-1)   + r(1) * sk
    ecmc(i)  =   ecmc(i-1)  + r(2) * sec
    wmc(i)   =   wmc(i-1)   + r(3) * sw
    call find_chi2(xd,yd,errs,rv0mc(i),kmc(i),ecmc(i),wmc(i),t0mc(i),Pmc(i),chi2_new,datas)
    q = exp( ( chi2_old - chi2_new ) * 0.5  )
    if ( q < r(4) ) then
      rv0mc(i) = rv0mc(i-1)
        kmc(i) = kmc(i-1)
       ecmc(i) = ecmc(i-1)
        wmc(i) = wmc(i-1)
    else
      chi2_old = chi2_new
    end if
  end do

end subroutine
