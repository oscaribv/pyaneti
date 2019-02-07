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
use constants

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
use constants
implicit none
!In/Out variables
  real(kind=mireal), intent(in) :: t0, e, w, P
  real(kind=mireal), intent(out) :: tp
!Local variables
  real(kind=mireal) :: theta_p
  real(kind=mireal) :: ereal, eimag

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
use constants
implicit none
!In/Out variables
  integer, intent(in) :: dt
  real(kind=mireal), intent(in) , dimension(0:dt-1) :: t
  real(kind=mireal), intent(out), dimension(0:dt-1) :: ta
  real(kind=mireal), intent(in) :: t0, e, w, P
!Local variables
  real(kind=mireal) :: tp

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
use constants
implicit none
!In/Out variables
  integer, intent(in) :: dt
  real(kind=mireal), intent(in) , dimension(0:dt-1) :: t
  real(kind=mireal), intent(out), dimension(0:dt-1) :: ta
  real(kind=mireal), intent(in) :: tp, e, P
!Local variables
  integer :: i,n
  real(kind=mireal), dimension(0:dt-1) :: ma, f, df, eimag, ereal, sinma

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
use constants
implicit none

!In/Out variables
  integer, intent(in) :: lenvec
  real(kind=mireal), intent(in), dimension(0:lenvec-1) :: mstar, rstar, P
  real(kind=mireal), intent(out), dimension(0:lenvec-1) :: a
!Local variables
  real(kind=mireal), dimension(0:lenvec-1) :: R_SI, GM_SI

  R_SI  = rstar(:) * S_radius_SI
  GM_SI = mstar(:) * S_GM_SI

  !Get scaled semi-major axis from 3rd Kepler law
  a(:) = 0.25d0 * GM_SI(:) * ( P * 24.d0 * 3600.d0 ) * ( P * 24.d0 * 3600.d0 )
  a(:) = a(:) / R_SI(:) / R_SI(:) / R_SI(:) / pi / pi 
  a(:) = a(:)**(1.d0/3.d0)

end subroutine

!Gelman and Rubin statistics
subroutine gr_test(par_chains,nchains,nconv,is_cvg)
use constants
implicit none
  
!In/Out variables
  integer, intent(in) :: nchains, nconv
  real(kind=mireal), intent(in), dimension(0:nconv-1,0:nchains-1) :: par_chains
  logical, intent(out) :: is_cvg
!Local variables
  real(kind=mireal) :: W, B, V, R
  real(kind=mireal) :: thetajj, delta
  real(kind=mireal), dimension(0:nchains-1):: sj2, thetaj
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
use constants
implicit none

!In/Out variables
  real(kind=mireal), intent(out) :: z
  real(kind=mireal), intent(in) :: a
!Internal variables
  real(kind=mireal) :: x

  !Thesis of Kaiser, Alexander D
  !Computational Experiments in Markov Chain Monte Carlo
  call random_number(x)
  z = ( a - 2.d0 + 1.d0/a ) * x*x + 2.d0 * (1.d0 - 1.d0/a ) * x + 1.d0/a

end subroutine

subroutine check_e(es,ec,is_good)
use constants
implicit none

  real(kind=mireal), intent(in) :: es, ec
  logical, intent(out) :: is_good

  is_good = .true.

  if ( .not. (es*es + ec*ec) < 1.d0  )  is_good = .false.

end subroutine

subroutine check_us(u1,u2,is_good)
use constants
implicit none

  real(kind=mireal), intent(in) :: u1, u2
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
use constants
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
use constants
implicit none

  !In/Out variables
  integer, intent(in) :: n
  real(kind=mireal), intent(in) :: mu, sigma
  real(kind=mireal), intent(out), dimension(0:n-1) :: valor
  !Local variables
  real(kind=mireal), dimension(0:2*n-1) :: r_real

  call random_number(r_real)

  valor(:) = sqrt( - 2.d0 * log(r_real(0:n-1)) ) * &
             cos( two_pi * r_real(n:2*n-1))

  valor(:) = sigma * valor(:) + mu

end subroutine



subroutine get_a_err(mstar_mean,mstar_sigma,rstar_mean,rstar_sigma,P,amean,aerr)
use constants
implicit none

!In/out variables
  real(kind=mireal), intent(in) :: mstar_mean, mstar_sigma, rstar_mean, rstar_sigma, P
  real(kind=mireal), intent(out) :: amean,aerr
!Local variables
  real(kind=mireal) :: dadm, dadr, per
  real(kind=mireal) :: R_SI, M_SI, G_SI
  real(kind=mireal) :: R_sigma_SI, M_sigma_SI
  real(kind=mireal) :: tercio, cpi2

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
use constants
implicit none
  integer, intent(in) :: n
  real(kind=mireal), intent(in) :: chi2(0:n-1)
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
use constants
implicit none

  integer, intent(in) :: npars
  integer, intent(in), dimension(0:npars-1) :: wtf
  real(kind=mireal), intent(in), dimension(0:2*npars-1) :: lims
  real(kind=mireal), intent(in), dimension(0:npars-1) :: pars
  real(kind=mireal), intent(out), dimension(0:npars-1) :: pars_out
!Local
  integer :: n, j
  real(kind=mireal) :: r_real

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
use constants
implicit none

  integer, intent(in) :: npars
  real(kind=mireal), intent(in), dimension(0:2*npars-1) :: lims
  real(kind=mireal), intent(out), dimension(0:npars-1) :: pars_out
  character, intent(in), dimension(0:npars-1) :: fit_pars
!Local
  integer :: j
  real(kind=mireal) :: r_real

  do j = 0, npars - 1
    if ( fit_pars(j) == 'f' ) then
       pars_out(j) = lims(j*2)
    else if ( fit_pars(j) == 'u' ) then
      call random_number(r_real)
      pars_out(j) = lims(2*j+1) - lims(2*j)
      pars_out(j) = lims(2*j) + r_real * pars_out(j)
    else if ( fit_pars(j) == 'g' ) then
      call gauss_random_bm(lims(2*j),lims(2*j+1),pars_out(j),1)
    else if ( fit_pars(j) == 'j' ) then
      call random_number(r_real)
      pars_out(j) = lims(2*j+1) - lims(2*j)
      pars_out(j) = lims(2*j) + r_real * pars_out(j)
    else if ( fit_pars(j) == 'm' ) then
      call random_number(r_real)
      pars_out(j) = r_real * lims(2*j+1)
    end if
  end do

end subroutine

subroutine rhotoa(rho,P,a,n)
use constants
implicit none

  !In/Out variables
  integer, intent(in) :: n
  real(kind=mireal), intent(in) :: rho, P(0:n-1)
  real(kind=mireal), dimension(0:n-1), intent(out) :: a

  !rho*^(1/3) parametrization
  a(:) = rho*(G_cgs*P(:)*P(:)*sind2/3.0d0/pi)**(1.d0/3.d0)
  !rho* parametrization
  !a(:) = (rho*G_cgs*P(:)*P(:)*sind2/3.0d0/pi)**(1.d0/3.d0)

end subroutine

subroutine ewto(ew1,ew2,e,w,n)
use constants
implicit none

  !In/Out variables
  integer, intent(in) :: n
  real(kind=mireal), intent(in), dimension(0:n-1) :: ew1, ew2
  real(kind=mireal), intent(out), dimension(0:n-1) :: e, w

    e(:) = ew1(:) * ew1(:) + ew2(:) * ew2(:)
    w(:) = atan2(ew1(:),ew2(:))

end subroutine

subroutine btoi(b,a,e,w,i,n)
use constants
implicit none

  !In/Out variables
  integer, intent(in) :: n
  real(kind=mireal), intent(in), dimension(0:n-1) :: b, a, e, w
  real(kind=mireal), intent(out), dimension(0:n-1) :: i

  i(:) = acos( b(:) / a(:) * ( 1.d0 + e(:) * sin(w(:)) ) / ( 1.d0 - e(:)*e(:) ) )

end subroutine


  subroutine fcdist(a,b,c,n,m)
  use constants
  implicit none
  !
  integer :: n,m
  real(kind=mireal), dimension(0:n-1) :: a
  real(kind=mireal), dimension(0:m-1) :: b
  real(kind=mireal), dimension(0:n-1,0:m-1) :: c
  !
  integer :: o

  do o = 0, n - 1
    c(o,:) = b(:) - a(o)
  end do

  end subroutine

  !https://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
  subroutine inverse(m,cc,n)
  !============================================================
  ! Inverse matrix
  ! Method: Based on Doolittle LU factorization for Ax=b
  ! Alex G. December 2009
  !-----------------------------------------------------------
  ! input ...
  ! a(n,n) - array of coefficients for matrix A
  ! n      - dimension
  ! output ...
  ! c(n,n) - inverse matrix of A
  ! comments ...
  ! the original matrix a(n,n) will be destroyed
  ! during the calculation
!===========================================================
  use constants
  implicit none
    integer, intent(in) :: n
    real(kind=mireal), intent(in) :: m(0:n-1,0:n-1)
    real(kind=mireal), intent(out) :: cc(0:n-1,0:n-1)
    !
    real(kind=mireal) :: a(n,n), c(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
    real(kind=mireal) :: coeff
    integer :: i, j, k

    a = m

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    !  step 1: forward elimination
    do k=1, n-1
      do i=k+1,n
        coeff=a(i,k)/a(k,k)
        L(i,k) = coeff
        do j=k+1,n
           a(i,j) = a(i,j)-coeff*a(k,j)
        end do
      end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
         U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do

    cc = c

  end subroutine inverse

!Subroutine to find the determinant of a square matrix
!Taken from https://github.com/ashwith/workspace/blob/master/fortran/determinant.f90
!Modified as a subroutine
!----------------------------------------------------------------------------------------
!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
  subroutine findlogddet(a,det,n)
  use constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  real(kind=mireal), DIMENSION(n,n), intent(in) :: a
  real(kind=mireal), intent (out) :: det
  real(kind=mireal) :: m, temp
  real(kind=mireal), DIMENSION(n,n) :: matrix
  INTEGER :: i, j, k, l
  LOGICAL :: DetExists = .TRUE.
  matrix = a
  l = 1
  !Convert to upper triangular form
  DO k = 1, n-1
    IF (matrix(k,k) == 0) THEN
      DetExists = .FALSE.
      DO i = k+1, n
        IF (matrix(i,k) /= 0) THEN
          DO j = 1, n
            temp = matrix(i,j)
            matrix(i,j)= matrix(k,j)
            matrix(k,j) = temp
          END DO
          DetExists = .TRUE.
          l=-l
          EXIT
        ENDIF
      END DO
      IF (DetExists .EQV. .FALSE.) THEN
        det = 0
        return
      END IF
    ENDIF
    DO j = k+1, n
      m = matrix(j,k)/matrix(k,k)
      DO i = k+1, n
        matrix(j,i) = matrix(j,i) - m*matrix(k,i)
      END DO
    END DO
  END DO

  !Calculate determinant by finding product of diagonal elements
  det = log(abs(real(l)))
  DO i = 1, n
    det = det + log(abs(matrix(i,i)))
  END DO

  end subroutine findlogddet

  subroutine fill_diag(v,M,n)
  use constants
  implicit none
  !
  integer, intent(in) :: n
  real(kind=mireal), intent(in) :: v(0:n-1)
  real(kind=mireal), intent(out) :: M(0:n-1,0:n-1)
  !
  integer :: i

  M = 0.d0
  do i = 0, n - 1
    M(i,i) = v(i)
  end do

  end subroutine fill_diag
