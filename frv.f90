!------------------------------------------------------------
!                         frv.f90
! This file contains subroutines to calculate Marcov Chain 
! Monte Carlo simulations in order to obtain planet parameters
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

  ta(:) = sqrt( (1. + ec) / (1. -ec) ) * tan(ta(:)*0.5)
  ta(:) = 2. * atan(ta(:))

end subroutine

!-----------------------------------------------------------
! This subroutine computes the circular radial velocity 
! curve for a set of values t. The subroutine returns 
! a vector (rv of the same size that t) by solving:
!  $ rv = rv0 - k [ sin ( 2*\pi ( t - t_0) / P ) ] $
! Where the parameters are the typical for a RV curve
!------------------------------------------------------------
subroutine rv_circular(t,rv0,t0,k,P,rv,ts)
implicit none

!In/Out variables
  integer, intent(in) :: ts
  double precision, intent(in), dimension(0:ts-1)  :: t
  double precision, intent(out), dimension(0:ts-1) :: rv
  double precision, intent(in) :: k, rv0, t0, P
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
!
  rv(:) = rv0 - k * sin( 2.*pi*( t(:) - t0) / P )

end subroutine

!-----------------------------------------------------------
! This subroutine computes the radial velocity 
! curve for an eccentric orbit, given a set of values t. 
!The subroutine returns  a vector (rv of the same size that t)
! by solving:
!  $ rv = rv0 + k [ cos ( \theta + \omega ) 
!             + e * cos ( \omega ) ] $
!  Where the parameters are the typical for a RV curve
!------------------------------------------------------------
subroutine rv_curve(t,rv0,t0,k,P,ec,w,rv,ts)
implicit none

!In/Out variables
  integer, intent(in) :: ts
  double precision, intent(in), dimension(0:ts-1)  :: t
  double precision, intent(out), dimension(0:ts-1) :: rv
  double precision, intent(in) :: k, rv0, t0, P, ec, w
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:ts-1) :: ma, ta
  double precision :: delta = 1e-4
  integer :: imax
!External function
  external :: find_anomaly
!

  imax = int(1e5)
  !Calculate the mean anomaly from the input values
  ma(:) = 2.* pi * ( t(:) - t0 ) / P
  !Obtain the eccentric anomaly by using find_anomaly
  call find_anomaly(ma,ta,ec,delta,imax,ts)

  rv(:) = rv0 + k * ( cos(ta(:) + w ) + ec * cos(w) )
  
end subroutine

!-----------------------------------------------------------
! This routine calculates the chi square for a RV curve
! given a set of xd-yd data points
! It takes into acount the possible difference in systematic
! velocities for different telescopes.
!Input parameters are:
! xd, yd, errs -> set of data to fit (array(datas))
! tlab -> Telescope labels (array of integers number)
! rv0  -> array for the different systemic velocities,
!         its size is the number of telescopes
! k, ec, w, t0, P -> typical planet parameters
! ics -> is circular flag, True or False.
! datas, nt -> sizes of xd,yd, errs (datas) and rv0(nt)  
!Output parameter:
! chi2 -> a double precision value with the chi2 value
!-----------------------------------------------------------
subroutine find_chi2(xd,yd,errs,tlab,rv0,k,ec,w,t0,P,chi2,isc,datas,nt)
implicit none

!In/Out variables
  integer, intent(in) :: datas, nt
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1)  :: tlab
  double precision, intent(in), dimension(0:nt-1)  :: rv0
  double precision, intent(in)  :: k, t0, P, ec, w
  logical, intent(in)  :: isc
  double precision, intent(out) :: chi2
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:datas-1) :: model, res
  integer :: i, tel
!External function
  external :: rv_circular, rv_curve

  chi2 = 0.0
  tel = 0
  i = 0
  !If we want a circular fit, let us do it!
  if ( isc ) then
    do i = 0, datas-1
      if ( tel .ne. tlab(i) ) tel = tel + 1 
      call rv_circular(xd(i),rv0(tel),t0,k,P,model(i),1)
    end do
  !If we want an eccentric fit, let us do it!
  else
    do i = 0, datas-1
      if ( tel .ne. tlab(i)  ) tel = tel + 1 
      call rv_curve(xd(i),rv0(tel),t0,k,P,ec,w,model(i),1)
    end do
  end if

  !Let us calculate the residuals
  ! chi^2 = \Sum_i ( M - O )^2 / \sigma^2
  res(:) = model(:) - yd(:)  
  res(:) = res(:) * res(:) / ( errs(:) * errs(:) )
  chi2 = sum(res)

end subroutine

!-----------------------------------------------------------
! This routine calculates a MCMC chain by using metropolis-
! hasting algorithm
! given a set of xd-yd data points
!Input parameters are:
! xd, yd, errs -> set of data to fit (array(datas))
! tlab -> Telescope labels (array of integers number)
! rv0  -> array for the different systemic velocities,
!         its size is the number of telescopes
! k, ec, w, t0, P -> typical planet parameters
! prec -> precision of the step size
! maxi -> maximum number of iterations
! thin factor -> number of steps between each output
! chi2_toler -> tolerance of the reduced chi^2 (1 + chi2_toler)
! ics -> is circular flag, True or False.
! datas, nt -> sizes of xd,yd, errs (datas) and rv0(nt)  
!Output parameter:
! This functions outputs a file called mh_rvfit.dat
!-----------------------------------------------------------
subroutine metropolis_hastings_rv(xd,yd,errs,tlab,rv0mc,kmc,ecmc,wmc,t0mc,Pmc,prec,maxi,thin_factor,chi2_toler,ics,datas,nt)
implicit none

!In/Out variables
  integer, intent(in) :: maxi, thin_factor, datas, nt
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1)  :: tlab
  double precision, intent(inout), dimension(0:nt-1)  :: rv0mc
  double precision, intent(inout)  :: kmc,t0mc, Pmc, ecmc, wmc, prec, chi2_toler
  !f2py intent(in,out)  ::rv0mc, kmc,t0mc, Pmc, ecmc, wmc, prec, chi2_toler
  logical, intent(in) :: ics
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision :: chi2_old, chi2_new, chi2_red
  double precision, dimension(0:nt-1)  :: rv0mcn
  double precision :: kmcn, t0mcn, Pmcn, ecmcn, wmcn
  double precision  :: sk, srv0, st0, sP, sec, sw
  double precision  :: q
  integer :: i, j, n, nu
  real, dimension(0:5+nt) :: r
!external calls
  external :: init_random_seed, find_chi2

  !Calculate the step size based in the actual value of the
  !parameters and the prec variable
  sk   = kmc     * prec
  srv0 = kmc     * prec
  st0  = t0mc    * prec
  sP   = Pmc     * prec
  sec  = 1.      * prec
  sw   = 2.*pi  * prec

  !Let us estimate our fist chi_2 value
  call find_chi2(xd,yd,errs,tlab,rv0mc,kmc,ecmc,wmc,t0mc,Pmc,chi2_old,ics,datas,nt)
  !Calculate the degrees of freedom
  nu = datas - 5 - nt
  !If we are fixing a circular orbit, ec and w are not used 
  if ( isc ) nu = nu + 2 
  !Print the initial cofiguration
  print *, ''
  print *, 'Starting MCMC calculation'
  print *, 'Initial Chi_2: ', chi2_old,'nu =', nu
  chi2_red = chi2_old / nu
  !Call a random seed 
  call init_random_seed()

  !Let us start the otput file
  open(unit=101,file='mh_rvfit.dat',status='unknown')
  write(101,*)'# i chi2 chi2_red k ec w t0 P rv0mc(vector)'
  write(101,*) 0,chi2_old,chi2_red,kmc, ecmc, wmc, t0mc, Pmc,rv0mc
  !Initialize the values
  n = 1
  i = 1

  !The infinite cycle starts!
  do while ( chi2_red >= 1 + chi2_toler .and. i <= maxi - 1 )
    !r will contain random numbers for each variable
    !Let us add a random shift to each parameter
    call random_number(r)
    r(1:5+nt) = ( r(1:5+nt) - 0.5) * 2.
    kmcn   =   kmc   + r(1) * sk
    ecmcn  =   ecmc  + r(2) * sec
    wmcn   =   wmc   + r(3) * sw
    t0mcn  =   t0mc  + r(4) * st0
    Pmcn   =   Pmc   + r(5) * sP
    do j = 0, nt-1
      rv0mcn(j) =   rv0mc(j) + r(6+j) * srv0
    end do
    !Let us calculate our new chi2
    call find_chi2(xd,yd,errs,tlab,rv0mcn,kmcn,ecmcn,wmcn,t0mcn,Pmcn,chi2_new,ics,datas,nt)
    !Ratio between the models
    q = exp( ( chi2_old - chi2_new ) * 0.5  )
    !If the new model is better, let us save it
    if ( q > r(0) ) then
      chi2_old = chi2_new
       rv0mc = rv0mcn
         kmc = kmcn
        ecmc = ecmcn
         wmc = wmcn
        t0mc = t0mcn
         Pmc = Pmcn
    end if
    !Calculate the reduced chi square
    chi2_red = chi2_old / nu
    !Save the data each thin_factor iteration
    if ( mod(i,thin_factor) == 0 ) then
      print *, 'iter ',i,'  of ',maxi
      print *, 'chi2 = ',chi2_old,'reduced chi2 =', chi2_red
      write(101,*) i,chi2_old,chi2_red,kmc, ecmc, wmc, t0mc, Pmc, rv0mc
      n = n + 1
    end if
    i = i + 1
  end do

  print *, 'Final chi2 = ',chi2_old,'. Final reduced chi2 =', chi2_red

  close(101)

end subroutine

!-------------------------------------------------------

