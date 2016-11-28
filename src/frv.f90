!-------------------------------------------------------------------------------
!                                    frv.f90
!       This file contains subroutines to calculate Marcov Chain
!     Monte Carlo simulations in order to obtain planet parameters
!      The subroutines can be called from python by using f2py
!             They also can be used in a fortran program
!              Date --> Feb  2016, Oscar Barragán
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
subroutine rv_curve_mp(t,rv0,t0,k,P,e,w,alpha,beta,rv,ts,npl)
implicit none

!In/Out variables
  integer, intent(in) :: ts, npl
  double precision, intent(in), dimension(0:ts-1)  :: t
  double precision, intent(out), dimension(0:ts-1) :: rv
  double precision, intent(in), dimension(0:npl-1) :: k, t0, P, e, w
  double precision, intent(in) :: rv0, alpha, beta
!Local variables
  double precision, dimension(0:ts-1) :: ta
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
! xd, yd, errs -> set of data to fit (array(datas))
! tlab -> Telescope labels (array of integers number)
! rv0  -> array for the different systemic velocities,
!         its size is the number of telescopes
! k, ec, w, t0, P -> typical planet parameters
! datas, nt -> sizes of xd,yd, errs (datas) and rv0(nt)  
!Output parameter:
! chi2 -> a double precision value with the chi2 value
!-----------------------------------------------------------
subroutine find_chi2_rv(xd,yd,errs,tlab,params,jitter,flag,chi2,datas,nt,npl)
implicit none

!In/Out variables
  integer, intent(in) :: datas, nt, npl
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1) :: tlab
  double precision, intent(in), dimension(0:6+nt,0:npl-1) :: params
  double precision, intent(in) :: jitter
  logical, intent(in)  :: flag(0:3)
  double precision, intent(out) :: chi2
!Local variables
  double precision, dimension(0:npl-1) :: t0, P, e, w, k
  double precision, dimension(0:nt-1)  :: rv0
  double precision  :: alpha, beta
  double precision, dimension(0:datas-1) :: model, res
  integer :: i, tel
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
  if ( flag(1) ) then
    e(:) = params(2,:) * params(2,:) + params(3,:) * params(3,:)
    w(:) = atan2( params(2,:),params(3,:) ) 
  end if
  if ( flag(2) ) k(:) = 1.d1**params(4,:)
  if ( flag(3) ) rv0(:) = 1.d1**params(7:6+nt,0)

  is_limit_good = .true.

  do i = 0, npl - 1
    if ( e(i) > 1.0d0 ) then
      is_limit_good = .false.
    end if
  end do

  if ( is_limit_good ) then

    !Telescope label counter
    tel = 0

    do i = 0, datas-1
     if ( tel .ne. tlab(i)  ) tel = tel + 1
      call rv_curve_mp(xd(i),rv0(tel),t0,k,P,e,w,alpha,beta,model(i),1,npl)
    end do

    !Let us obtain chi^2
    res(:) = ( model(:) - yd(:) ) / sqrt( errs(:)**2 + jitter**2 )
    chi2 = dot_product(res,res)

  else

    chi2 = huge(0.e0)

  end if

end subroutine
