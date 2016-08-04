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
!Modified to work with f90
!use iso_fortran_env, only: int64
implicit none
integer,parameter :: int32 = selected_int_kind(32)
integer, allocatable :: seed(:)
integer :: i, n, un, istat, dt(8), pid
integer(int32) :: t

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
      t = (dt(1) - 1970) * 365_int32 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_int32 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_int32 * 60 * 60 * 1000 &
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
  integer(int32) :: s
  if (s == 0) then
     s = 104729
  else
     s = mod(s, 4294967296_int32)
  end if
  s = mod(s * 279470273_int32, 4294967291_int32)
  lcg = int(mod(s, int(huge(0), int32)), kind(0))
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
  double precision :: ereal, eimag
  double precision :: pi = 3.1415926535897d0

  ereal = e + cos( pi / 2.d0  - w)
  eimag = sqrt( 1.d0 - e * e ) * sin( pi/ 2.d0  - w )
  theta_p = atan2(eimag, ereal )
  theta_p = theta_p - e * sin( theta_p )

  tp = t0 - theta_p * p / 2.d0 / pi

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
subroutine find_anomaly(t,t0,e,w,P,ta,dt)
implicit none
!In/Out variables
  integer, intent(in) :: dt
  double precision, intent(in) , dimension(0:dt-1) :: t
  double precision, intent(out), dimension(0:dt-1) :: ta
  double precision, intent(in) :: t0, e, w, P
!Local variables
  integer :: i,n
  double precision, dimension(0:dt-1) :: ma, f, df, eimag, ereal
  double precision :: two_pi = 2.d0*3.1415926535897932384626d0
  double precision :: uno, tp
  double precision :: fmin=1.d-6
  integer :: imax = int(1e8)
!
  uno = 1.0d0

  call find_tp(t0,e,w,P,tp)

  !Calculate the mean anomaly
  ma = two_pi * ( t - tp ) / P

  !calculate the eccentric anomaly
  !Using Newthon-Raphson algorithm
  ta(:) = ma(:)
  f(:) = fmin * 1.0d1
  n = 0

  do i = 0, dt-1
    do while ( abs(f(i)) > fmin .and. n < imax )
      f(i)   = ta(i) - e * sin(ta(i)) - ma(i)
      df(i)  = uno - e * cos(ta(i))
      ta(i)  = ta(i) - f(i) / df(i)
      n = n + 1
    end do
  end do

  if ( n > imax ) then !This should never happen!
    print *, 'I am tired, too much Newton-Raphson for me!'
    print *, e, f
    stop
  end if 


  !calculate the true anomaly
  !Relation between true anomaly(ta) and eccentric anomaly(ea) is
  !tan(ta) = sqrt(1-e^2) sin (ea) / ( cos(ea) - e ) https://en.wikipedia.org/wiki/True_anomaly
  !In a complex plane, this is =  (cos(ea) - e) + i (sqrt(1-e^2) *sin(ea) ) 
  !with modulus = 1 - e cos(ea) 
  eimag = ( sqrt(uno-e*e) * sin(ta) ) !/ (uno-e*cos(ta))
  ereal = ( cos (ta) - e ) !/ (uno-e*cos(ta))
  !Therefore, the tue anomaly is 
  ta = atan2(eimag,ereal)

end subroutine

!Get semi-major axis assuming we know the stellar parameters

subroutine get_a_scaled(mstar,rstar,P,a,lenvec)
implicit none

!In/Out variables
  integer, intent(in) :: lenvec
  double precision, intent(in), dimension(0:lenvec-1) :: mstar, rstar, P
  double precision, intent(out), dimension(0:lenvec-1) :: a
!Local variables
  double precision :: pi = 3.1415926535897d0
  double precision :: S_radius_SI = 6.957d8 !R_sun
  double precision :: S_GM_SI = 1.3271244d20 ! G M_sun
  double precision, dimension(0:lenvec-1) :: R_SI, GM_SI

  R_SI  = rstar(:) * S_radius_SI
  GM_SI = mstar(:) * S_GM_SI

  !Get scaled semi-major axis from 3rd Kepler law
  a(:) = 0.25d0 * GM_SI(:) * ( P * 24 * 3600 ) * ( P * 24 * 3600 )
  a(:) = a(:) / R_SI(:) / R_SI(:) / R_SI(:) / pi / pi 
  a(:) = a(:)**(1.d0/3.d0)

end subroutine

!Gelman and Rubin statistics
subroutine gr_test(par_chains,nchains,nconv,is_cvg)
implicit none
  
!In/Out variables
  integer, intent(in) :: nchains, nconv
  double precision, intent(in), dimension(0:nconv-1,0:nchains-1) :: par_chains
  logical, intent(out) :: is_cvg
!Local variables
  double precision :: W, B, V, R
  double precision :: thetajj, delta
  double precision, dimension(0:nchains-1):: sj2, thetaj
  integer :: i

  is_cvg = .false.
  delta = 5.0d-2

  !Let us calculate the mean of each chain
  sj2(:) = 0.0d0
  do i = 0, nchains - 1
    thetaj(i) = sum(par_chains(:,i)) / nconv
      sj2(i) = dot_product(( par_chains(:,i) - thetaj(i) ), &
                           ( par_chains(:,i) - thetaj(i) ) )
    sj2(i) = sj2(i) / ( dble(nconv) - 1.d0 )
  end do

  !Whithin chain variance
  W = sum(sj2(:)) / nchains

  thetajj = sum(thetaj(:)) / nchains

  !Between chain variance
  B = dot_product( ( thetaj(:) - thetajj),( thetaj(:) - thetajj) )
  B = nconv * B / ( nchains - 1 )

  !Estimated variance
  V = W - (W - B) / nconv

  !Potential scale reduction factor
  R = sqrt ( V / W )


  if ( R < 1.d0 + delta ) then
      is_cvg = .true. 
  else
      print *, 'R     = ', R
      print *, 'delta = ', 1.d0 + delta
  end if

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

  if ( z >= 1.d0 / a .and. z <= a ) then
    z = 1.d0 / sqrt(z) 
  else
    z = 0.0d0
  end if

end subroutine

subroutine check_e(es,ec,is_good)
implicit none

  double precision, intent(in) :: es, ec
  logical, intent(out) :: is_good

  is_good = .true.

  if ( .not. (es*es + ec*ec) < 1.d0  ) then
    is_good = .false.
  end if

end subroutine

subroutine check_us(u1,u2,is_good)
implicit none

  double precision, intent(in) :: u1, u2
  logical, intent(out) :: is_good

  is_good = .true.

  if ( u1 + u2 > 1.d0 ) is_good = .false.

end subroutine



!Subroutine to create random integers between 0 and n
subroutine random_int(r_int,n)
implicit none

  !In/Out variables
  integer, intent(in) :: n
  integer, intent(out), dimension(0:n-1) :: r_int
  !Local variables
  double precision :: r_real
  integer :: i

  do i = 0, n - 1
    r_int(i) = i
    do while ( r_int(i) == i )
    call random_number(r_real)
    r_int(i) = int( r_real *  n ) 
    end do
  end do


end subroutine


!Subroutine to check the limits of the parameters
subroutine check_limits(params,limits,is_true,npar)
implicit none

  !In/Out variables
  integer, intent(in) :: npar
  double precision, intent(in) :: params(0:npar-1), &
                                  limits(0:2*npar-1)
  logical, intent(out) :: is_true
  !Local variables
  integer :: i, j

  is_true = .true.
  j = 0
  do i = 0, npar - 1
    if ( params(i) <= limits(j) .or. &
         params(i) >= limits(j+1) ) then
      is_true = .false.
      exit
    end if
    j = j + 2
  end do

end subroutine

!Create a normal distribution based on Box-Muller
subroutine gauss_random_bm(mu,sigma,valor,n)
implicit none

  !In/Out variables
  integer, intent(in) :: n
  double precision, intent(in) :: mu, sigma
  double precision, intent(out), dimension(0:n-1) :: valor
  !Local variables
  double precision, dimension(0:2*n-1) :: r_real
  double precision  :: two_pi = 2*3.1415926535897932384626d0

  call random_number(r_real)

  valor(:) = sqrt( - 2.d0 * log(r_real(0:n-1)) ) * &
             cos( two_pi * r_real(n:2*n-1))

  valor(:) = sigma * valor(:) + mu

end subroutine

subroutine print_chain_data(chi2,n)
implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: chi2(0:n-1)
  character(LEN=20) :: fto = "(A,F10.2)"

  write(*,*), '==========================='
  write(*,*), '     Chain statistics      '
  write(*,*), '==========================='
  write(*,*), 'chain |  reduced chi^2 '
  write(*,fto), ' best  : ',minval(chi2)
  write(*,fto), ' worst : ',maxval(chi2)
  write(*,fto), ' mean  : ', sum(chi2) / n
  write(*,*), '==========================='


end subroutine
