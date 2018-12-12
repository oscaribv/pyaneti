!-----------------------------------------------------------
!                         todo.f90
! This file has a some generic functions to be called for more
!    general subroutines or programs.
! The subroutines can be called from python by using f2py
! They also can be used in a fortran program
!              Date --> Feb  2016, Oscar Barragán
!------------------------------------------------------------

!http://stackoverflow.com/questions/18754438/generating-random-numbers-in-a-fortran-module
subroutine init_random_seed()

      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)

end subroutine

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

  !We know that the relation between theta_t0 = pi/2 - w
  !We have to calculate the eccentric anomaly by knowing this
  ereal = e + cos( pi / 2.d0  - w)
  eimag = sqrt( 1.d0 - e * e ) * sin( pi/ 2.d0  - w )
  theta_p = atan2(eimag, ereal )
  !Now we  have the eccentric anomaly, let us calcualte the mean anomaly
  theta_p = theta_p - e * sin( theta_p )
  !Time to talculate Tp
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
  double precision :: tp

  call find_tp(t0,e,w,P,tp)
  call find_anomaly_tp(t,tp,e,P,ta,dt)

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
subroutine find_anomaly_tp(t,tp,e,P,ta,dt)
implicit none
!In/Out variables
  integer, intent(in) :: dt
  double precision, intent(in) , dimension(0:dt-1) :: t
  double precision, intent(out), dimension(0:dt-1) :: ta
  double precision, intent(in) :: tp, e, P
!Local variables
  integer :: i,n
  double precision, dimension(0:dt-1) :: ma, f, df, eimag, ereal, sinma
  double precision :: two_pi = 2.d0*3.1415926535897932384626d0
  double precision :: uno
  double precision :: fmin=1.d-8, small = 1.d-5
  integer :: imax = int(1e8)
!
  uno = 1.0d0

  !Calculate the mean anomaly
  ma = two_pi * ( t - tp ) / P

  if ( e > small ) then !You have to calcuate your true anomaly, your orbit is not circular!

    sinma = sin(ma(:))
    !make guesses based on the expantion
    !Serie expantion, Murray and Dermott 1999, p.35
    ta(:) = ma(:) + e  * ( sinma(:) + &
            e * ( 0.5d0 * sin(2.d0*ma(:)) +  &
            e * 0.125d0 * ( 3.d0 * sin( 3.d0 * ma(:) ) - sinma(:)  ) &
            ) )

    !calculate the eccentric anomaly
    !Using Newthon-Raphson algorithm
    f(:) = ta(:) - e * sin(ta(:)) - ma(:)
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

  else

    !If your orbit is cirular, your true anomaly is the mean anomaly ;)
    ta(:) = ma(:)

  end if


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
  a(:) = 0.25d0 * GM_SI(:) * ( P * 24.d0 * 3600.d0 ) * ( P * 24.d0 * 3600.d0 )
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
  delta = 2.0d-2

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
  !equation 11.4  Gelman et al., 2004, Bayesian Data Analysis, 3rd edition
  V = W - (W - B) / nconv

  !Potential scale reduction factor
  R = sqrt ( V / W )

  if ( R < 1.d0 + delta ) is_cvg = .true.

end subroutine

!Subroutine to get Z <- g(z)
!Goodman & Weare, 2010 paper
subroutine find_gz(a,z)
implicit none

!In/Out variables
  double precision, intent(out) :: z
  double precision, intent(in) :: a
!Internal variables
  double precision :: x

  !Thesis of Kaiser, Alexander D
  !Computational Experiments in Markov Chain Monte Carlo
  call random_number(x)
  z = ( a - 2.d0 + 1.d0/a ) * x*x + 2.d0 * (1.d0 - 1.d0/a ) * x + 1.d0/a

end subroutine

subroutine check_e(es,ec,is_good)
implicit none

  double precision, intent(in) :: es, ec
  logical, intent(out) :: is_good

  is_good = .true.

  if ( .not. (es*es + ec*ec) < 1.d0  )  is_good = .false.

end subroutine

subroutine check_us(u1,u2,is_good)
implicit none

  double precision, intent(in) :: u1, u2
  logical, intent(out) :: is_good

  is_good = .true.

  if ( u1 + u2 > 1.d0 ) then
    is_good = .false.
  else if ( u1 > 1.d0 .or. u1 < 0.0d0 ) then
    is_good = .false.
  else if ( u2 > 1.d0 .or. u2 < -1.0d0 ) then
    is_good = .false.
   end if

end subroutine

!Subroutine to create random integers between 0 and n
subroutine random_int(r_int,mnv,mxv)
implicit none

  !In/Out variables
  integer, intent(in) :: mnv,mxv
  integer, intent(out), dimension(0:mxv-mnv-1) :: r_int
  !Local variables
  real :: r_real
  integer :: i, j

  do i = mnv, mxv
    j = i-mnv
    r_int(j) = i
    do while ( r_int(j) == i )
      call random_number(r_real)
      r_int(j) =  mnv + int ( r_real * ( 1 + mxv - mnv ) )
    end do
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
  double precision  :: two_pi = 2.d0*3.1415926535897932384626d0

  call random_number(r_real)

  valor(:) = sqrt( - 2.d0 * log(r_real(0:n-1)) ) * &
             cos( two_pi * r_real(n:2*n-1))

  valor(:) = sigma * valor(:) + mu

end subroutine



subroutine get_a_err(mstar_mean,mstar_sigma,rstar_mean,rstar_sigma,P,amean,aerr)
implicit none

!In/out variables
  double precision, intent(in) :: mstar_mean, mstar_sigma, rstar_mean, rstar_sigma, P
  double precision, intent(out) :: amean,aerr
!Local variables
  double precision :: dadm, dadr, per
  double precision :: S_radius_SI = 6.957d8 !R_sun
  double precision :: S_GM_SI = 1.3271244d20 ! G M_sun
  double precision :: R_SI, M_SI, G_SI
  double precision :: R_sigma_SI, M_sigma_SI
  double precision :: tercio, cpi2
  double precision  :: pi = 3.1415926535897932384626d0

  G_SI = 6.67408e-11
  R_SI = rstar_mean * S_radius_SI
  R_sigma_SI = rstar_sigma * S_radius_SI
  M_SI = mstar_mean * S_GM_SI / G_SI
  M_sigma_SI = mstar_sigma * S_GM_SI / G_SI
  per  = P * 3600.d0 * 24.d0
  tercio = 1.d0/3.d0
  cpi2 = 4.d0 * pi**2

  amean = per**2 * G_SI * M_SI / ( cpi2 * R_SI**3 )
  amean = amean**(tercio)

  dadm = ( G_SI * per**2 / ( cpi2 * M_SI**2) )**(tercio)
  dadm = dadm * tercio / R_SI

  dadr = - ( G_SI * M_SI*per**2 / cpi2 )**(tercio) / R_SI**2

  aerr = dadm**2 * M_sigma_SI**2 + dadr**2 * R_sigma_SI**2
  aerr = sqrt(aerr)

end subroutine

subroutine print_chain_data(chi2,n)
implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: chi2(0:n-1)
  character(LEN=20) :: fto = "(A,F10.2)"

  write(*,*) '=================================='
  write(*,*) '     Chain statistics      '
  write(*,*) '=================================='
  write(*,*) 'chain |  reduced chi^2 '
  write(*,fto) ' best  : ',minval(chi2)
  write(*,fto) ' worst : ',maxval(chi2)
  write(*,fto) ' mean  : ', sum(chi2) / n
  write(*,*) '=================================='

end subroutine

subroutine uniform_chains(pars,npars,wtf,lims,pars_out)
implicit none

  integer, intent(in) :: npars
  integer, intent(in), dimension(0:npars-1) :: wtf
  double precision, intent(in), dimension(0:2*npars-1) :: lims
  double precision, intent(in), dimension(0:npars-1) :: pars
  double precision, intent(out), dimension(0:npars-1) :: pars_out
!Local
  integer :: n, j
  double precision :: r_real

  j = 0
  do n = 0,  npars - 1
    if ( wtf(n) == 0 ) then
      pars_out(n) = pars(n)
    else
      call random_number(r_real)
      pars_out(n) = lims(j+1) - lims(j)
      pars_out(n) = lims(j) + r_real * pars_out(n)
    end if
    j = j + 2
  end do

end subroutine

subroutine create_chains(fit_pars,lims,pars_out,npars)
implicit none

  integer, intent(in) :: npars
  double precision, intent(in), dimension(0:2*npars-1) :: lims
  double precision, intent(out), dimension(0:npars-1) :: pars_out
  character, intent(in), dimension(0:npars-1) :: fit_pars
!Local
  integer :: j
  double precision :: r_real

  do j = 0, npars - 1
    if ( fit_pars(j) == 'f' ) then
       pars_out(j) = lims(j*2)
    else if ( fit_pars(j) == 'u' ) then
      call random_number(r_real)
      pars_out(j) = lims(2*j+1) - lims(2*j)
      pars_out(j) = lims(2*j) + r_real * pars_out(j)
    else if ( fit_pars(j) == 'g' ) then
      call gauss_random_bm(lims(2*j),lims(2*j+1),pars_out(j),1)
    end if
  end do

end subroutine

subroutine rhotoa(rho,P,a,n)
implicit none

  !In/Out variables
  integer, intent(in) :: n
  double precision, intent(in) :: rho, P(0:n-1)
  double precision, dimension(0:n-1), intent(out) :: a
  !Local variables
  double precision :: pi = 3.1415926535897932384626d0
  double precision :: G = 6.67508d-11*1.d3 !Gravitational constant

  a(:) = rho*(G*P(:)*P(:)*7464960000./3.0/pi)**(1./3.)

end subroutine

subroutine ewto(ew1,ew2,e,w,n)
implicit none

  !In/Out variables
  integer, intent(in) :: n
  double precision, intent(in), dimension(0:n-1) :: ew1, ew2
  double precision, intent(out), dimension(0:n-1) :: e, w

    e(:) = ew1(:) * ew1(:) + ew2(:) * ew2(:)
    w(:) = atan2(ew1(:),ew2(:))

end subroutine

subroutine btoi(b,a,e,w,i,n)
implicit none

  !In/Out variables
  integer, intent(in) :: n
  double precision, intent(in), dimension(0:n-1) :: b, a, e, w
  double precision, intent(out), dimension(0:n-1) :: i

  i(:) = acos( b(:) / a(:) * ( 1.d0 + e(:) * sin(w(:)) ) / ( 1.d0 - e(:)*e(:) ) )

end subroutine
