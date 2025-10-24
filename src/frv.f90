!-------------------------------------------------------------------------------
!                                    frv.f90
!       This file contains subroutines to calculate Marcov Chain
!     Monte Carlo simulations in order to obtain planet parameters
!      The subroutines can be called from python by using f2py
!             They also can be used in a fortran program
!              Date --> Feb  2016, Oscar Barragan
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! This subroutine computes the radial velocity for a multiplanet system
! for a given set of values t.
!The subroutine returns  a vector (rv of the same size that t)
! by solving:
!  rv = rv0 + k [ cos ( theta + w ) + e * cos ( w ) ]
!  Where the parameters are the typical for a RV curve
! Input parameters:
! t(:)   -> time vector
! rv0    -> systemic velocity of a given instrument
! t0(:)  -> Transit epoch (for not transit planets is when the planet
!           is at 90 deg with respect to the sky)
! k(:)   -> RV semi-amplitude for each planet
! P(:)   -> Period of each planet
! e(:)   -> eccentricity of each planet
! w(:)   -> periastron passage of each planet
!alpha   -> linear trend
!beta    -> quadratic trend
!rv(:)   -> set of RV meassurements
!ts      -> size of t
!npl     -> numper of planets
!-------------------------------------------------------------------------------
subroutine rv_curve_mp(t,rv0,t0,P,e,w,k,alpha,beta,rv,ts,npl)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: ts, npl
  real(kind=mireal), intent(in), dimension(0:ts-1)  :: t
  real(kind=mireal), intent(out), dimension(0:ts-1) :: rv
  real(kind=mireal), intent(in), dimension(0:npl-1) :: t0, P, e, w, k
  real(kind=mireal), intent(in) :: rv0, alpha, beta
!Local variables
  real(kind=mireal), dimension(0:ts-1) :: ta
  integer :: i
!External function
  external :: find_anomaly
!

  !Added systemic velocity and linear and quadratic trends
  rv(:) = rv0 + (t(:)-t0(0))*alpha + (t(:)-t0(0))**2*beta

  !Now add the planet influence on the star
  !each i is for a different planet
  do i = 0, npl-1
   !Obtain the true anomaly by using find_anomaly
   call find_anomaly(t,t0(i),e(i),w(i),P(i),ta,ts)
   !Now compute the RV for the planet i
   rv(:) = rv(:) + k(i) * ( cos(ta(:) + w(i) ) + e(i) * cos(w(i)) )
  end do
  !The final RV induced by all planets

end subroutine


!-------------------------------------------------------------------------------
! This routine calculates the chi square for a RV curve
! given a set of xd-yd data points
! It takes into acount the possible difference in systematic
! velocities for different telescopes.
!Input parameters are:
! xd, yd, errs -> set of data to fit (array(n_data))
! tlab -> Telescope labels (array of integers number)
! rv0  -> array for the different systemic velocities,
!         its size is the number of telescopes
! k, ec, w, t0, P -> typical planet parameters
! n_data, nt -> sizes of xd,yd, errs (n_data) and rv0(nt)
!Output parameter:
! chi2 -> a real(kind=mireal) value with the chi2 value
!-----------------------------------------------------------
subroutine find_res_rv(xd,yd,tlab,params,flag,res,n_data,nt,npl)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: n_data, nt, npl
  real(kind=mireal), intent(in), dimension(0:n_data-1)  :: xd, yd
  integer, intent(in), dimension(0:n_data-1) :: tlab
  real(kind=mireal), intent(in), dimension(0:6+nt,0:npl-1) :: params
  logical, intent(in)  :: flag(0:3)
  real(kind=mireal), intent(out) :: res(0:n_data-1)
!Local variables
  real(kind=mireal), dimension(0:npl-1) :: t0, P, e, w, k
  real(kind=mireal), dimension(0:nt-1)  :: rv0
  real(kind=mireal)  :: alpha, beta
  real(kind=mireal), dimension(0:n_data-1) :: model
  logical :: is_limit_good
!External function
  external :: rv_curve_mp

  t0(:)  = params(0,:)
  P(:)   = params(1,:)
  e(:)   = params(2,:)
  w(:)   = params(3,:)
  k(:)   = params(4,:)
  alpha  = params(5,0)
  beta   = params(6,0)
  rv0(:) = params(7:6+nt,0)

  if ( flag(0) ) P(:) = 1.d1**params(1,:)
  if ( flag(1) ) call ewto(e,w,e,w,npl)
  if ( flag(2) ) k(:) = 1.d1**params(4,:)
  if ( flag(3) ) rv0(:) = 1.d1**params(7:6+nt,0)

  is_limit_good = .true.

  !Avoid solutions with eccentricities larger than 1 and smaller than zero
  if ( any( e > 1.d0 ) .or. any(e < 0.d0 ) ) is_limit_good = .false.

  !Avoid solutions with Doppler semi-amplitudes smaller than zero
  if (  any( k < 0.d0 ) ) is_limit_good = .false.

  if ( is_limit_good ) then

      call rv_curve_mp(xd,0.d0,t0,P,e,w,k,alpha,beta,model,n_data,npl)
      model(:) = model(:) + rv0(tlab(:))
      res(:)   =  yd(:) - model(:)
  else

      res(:) = huge(0.d0)

  end if

end subroutine



!-------------------------------------------------------------------------------
! This routine calculates the chi square for a RV curve
! given a set of xd-yd data points
! same as find_res_rv_py but now the params vector is an array
! this way is more compatible with python
!-----------------------------------------------------------
subroutine find_res_rv_py(xd,yd,tlab,params,rv0,alpha,beta,flag,res,n_data,nt,npl)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: n_data, nt, npl
  real(kind=mireal), intent(in), dimension(0:n_data-1)  :: xd, yd
  integer, intent(in), dimension(0:n_data-1) :: tlab
  real(kind=mireal), intent(in), dimension(0:5*npl-1) :: params
  real(kind=mireal), intent(in), dimension(0:nt-1)  :: rv0
  real(kind=mireal), intent(in)  :: alpha, beta
  logical, intent(in)  :: flag(0:2)
  real(kind=mireal), intent(out) :: res(0:n_data-1)
!Local variables
  integer :: i
  real(kind=mireal), dimension(0:npl-1) :: t0, P, e, w, k
  real(kind=mireal), dimension(0:n_data-1) :: model
  logical :: is_limit_good
!External function
  external :: rv_curve_mp

  do i = 1, npl
    t0(i-1)  = params(1*i-1)
    P(i-1)   = params(2*i-1)
    e(i-1)   = params(3*i-1)
    w(i-1)   = params(4*i-1)
    k(i-1)   = params(5*i-1)
  end do

  if ( flag(0) ) P(:) = 1.d1**P(:)
  if ( flag(1) ) call ewto(e,w,e,w,npl)
  if ( flag(2) ) k(:) = 1.d1**k(:)

  is_limit_good = .true.

  !Avoid solutions with eccentricities larger than 1 and smaller than zero
  if ( any( e > 1.d0 ) .or. any(e < 0.d0 ) ) is_limit_good = .false.

  !Avoid solutions with Doppler semi-amplitudes smaller than zero
  if (  any( k < 0.d0 ) ) is_limit_good = .false.

  if ( is_limit_good ) then

      call rv_curve_mp(xd,0.d0,t0,P,e,w,k,alpha,beta,model,n_data,npl)
      model(:) = model(:) + rv0(tlab(:))
      res(:)   =  yd(:) - model(:)
  else

      res(:) = huge(0.d0)

  end if

end subroutine

!-------------------------------------------------------------------------------
! This routine calculates the chi square for a RV curve
! given a set of xd-yd data points
! It takes into acount the possible difference in systematic
! velocities for different telescopes.
!Input parameters are:
! xd, yd, errs -> set of data to fit (array(n_data))
! tlab -> Telescope labels (array of integers number)
! rv0  -> array for the different systemic velocities,
!         its size is the number of telescopes
! k, ec, w, t0, P -> typical planet parameters
! n_data, nt -> sizes of xd,yd, errs (n_data) and rv0(nt)
!Output parameter:
! chi2 -> a real(kind=mireal) value with the chi2 value
!-----------------------------------------------------------
subroutine find_chi2_rv(xd,yd,errs,tlab,jrvlab,params,jitter,flag,chi2,n_data,nt,nj,npl)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: n_data, nt, npl, nj
  real(kind=mireal), intent(in), dimension(0:n_data-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:n_data-1) :: tlab, jrvlab
  real(kind=mireal), intent(in), dimension(0:6+nt,0:npl-1) :: params
  real(kind=mireal), dimension(0:nj-1), intent(in) :: jitter
  logical, intent(in)  :: flag(0:3)
  real(kind=mireal), intent(out) :: chi2
!Local variables
  real(kind=mireal), dimension(0:n_data-1) :: res

  call find_res_rv(xd,yd,tlab,params,flag,res,n_data,nt,npl)
  res(:)   = res(:) / sqrt( errs(:)**2 + jitter(jrvlab(:))**2 )
  chi2 = dot_product(res,res)

end subroutine
