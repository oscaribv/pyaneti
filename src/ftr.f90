!------------------------------------------------------------
!                         ftr.f90
! This file contains subroutines to calculate Marcov Chain 
! Monte Carlo simulations in order to obtain planet parameters
! from light curve fitting of transit planets
! The subroutines can be called from python by using f2py
! They also can be used in a fortran program
!              Date --> Feb  2016, Oscar Barragán
!------------------------------------------------------------

!-----------------------------------------------------------
!                     find_z
!  This suborutine finds the projected distance between
!  the star and planet centers. Eq. (5), ( z = r_sky) from
!  Winn, 2010, Transit and Occultations.
!------------------------------------------------------------
subroutine find_z(t,pars,z,ts)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: ts
  real(kind=mireal), intent(in), dimension(0:ts-1) :: t
  real(kind=mireal), intent(in), dimension(0:5) :: pars
  real(kind=mireal), intent(out), dimension(0:ts-1) :: z
!Local variables
  real(kind=mireal), dimension(0:ts-1) :: ta, swt
  real(kind=mireal) :: tp, P, e, w, i, a
  real(kind=mireal) :: si
!

  tp  = pars(0)
  P   = pars(1)
  e   = pars(2)
  w   = pars(3)
  i   = pars(4)
  a   = pars(5)

  !Obtain the true anomaly by using find_anomaly_tp
  call find_anomaly_tp(t,tp,e,P,ta,ts)

  swt = sin(w+ta)
  si = sin(i)

  where (swt > 0.0 ) !We have the planet in front of the star -> transit
  !z has been calculated
    z = a * ( 1.d0 - e * e ) * sqrt( 1.d0 - swt * swt * si * si ) &
        / ( 1.d0 + e * cos(ta) )
  elsewhere !We have the planet behind the star -> occulation
  !z has avalue which gives flux = 1 in Mandel & Agol module
      z = 1.d1
  end where

end subroutine

subroutine flux_tr(xd,trlab,pars,rps,ldc,&
           n_cad,t_cad,nbands,nradius,n_data,npl,flux_out)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: n_data, npl, nbands, nradius
  integer, intent(in) :: n_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension(0:n_data-1)  :: xd
  integer, intent(in), dimension(0:n_data-1)  :: trlab !this indicates the instrument label
  real(kind=mireal), intent(in), dimension(0:5,0:npl-1) :: pars
  real(kind=mireal), intent(in), dimension(0:nradius*npl-1) :: rps
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  real(kind=mireal), intent(in) :: t_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension (0:2*nbands-1) :: ldc
  real(kind=mireal), intent(out), dimension(0:n_data-1) :: flux_out !output flux model

  !If any band requieres binning, let us call the routines with binning
  if ( any(n_cad > 1) ) then
    !If we are fitting only one band, let's call the routine that computes it faster
    if ( nbands == 1 ) then
      call flux_tr_singleband(xd,pars,rps,ldc,n_cad(0),t_cad(0),n_data,npl,flux_out)
    !Let's  call the complicated multiband function
    else
      call flux_tr_multiband(xd,trlab,pars,rps,ldc,n_cad,t_cad,nbands,nradius,n_data,npl,flux_out)
    end if
  else !We no need binning and things can be faster!
      !If we are fitting only one band, let's call the routine that computes it faster
    if ( nbands == 1 ) then
      call flux_tr_singleband_nobin(xd,pars,rps,ldc,n_data,npl,flux_out)
    !Let's  call the complicated multiband function
    else
      call flux_tr_multiband_nobin(xd,trlab,pars,rps,ldc,nbands,nradius,n_data,npl,flux_out)
    end if
  end if

end subroutine

subroutine flux_tr_multiband(xd,trlab,pars,rps,ldc,&
           n_cad,t_cad,nbands,nradius,n_data,npl,flux_out)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: n_data, npl, nbands, nradius
  integer, intent(in) :: n_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension(0:n_data-1)  :: xd
  integer, intent(in), dimension(0:n_data-1)  :: trlab !this indicates the instrument label
  real(kind=mireal), intent(in), dimension(0:5,0:npl-1) :: pars
  real(kind=mireal), intent(in), dimension(0:nradius*npl-1) :: rps
!  real(kind=mireal), intent(in), dimension(0:nbands-1,0:npl-1) :: rps
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  real(kind=mireal), intent(in) :: t_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension (0:2*nbands-1) :: ldc
  real(kind=mireal), intent(out), dimension(0:n_data-1) :: flux_out !output flux model
!Local variables
  real(kind=mireal), dimension(0:n_data-1) :: mu
  real(kind=mireal) :: npl_dbl, u1(0:nbands-1), u2(0:nbands-1)
  real(kind=mireal), allocatable, dimension(:)  :: flux_unbinned
  real(kind=mireal)   :: flux_binned
  real(kind=mireal), allocatable, dimension(:)  :: t_unbinned, z
  integer :: i, n, j, control, lj, n_cadj
  integer, allocatable :: k(:)

  !This flag controls the multi-radius fits
  control = 1
  if (nradius == 1) control = 0

  npl_dbl = dble(npl)


  do n = 0, nbands - 1
    u1(n) = ldc(2*n)
    u2(n) = ldc(2*n+1)
  end do

  flux_out = 1.d0

  do n = 0, npl - 1

    do j = 0, n_data - 1

      !variable with the current telescope label
      lj = trlab(j)
      n_cadj = n_cad(lj)

      !Let us allocate memory in a smart way
      !If the previous point was the same telescope label, then the arrays have
      !the same dimensions, and we no need to allocate new memory
      if ( j < 1 .or. (j > 0 .and. (trlab(j-1) .ne. lj ) ) ) then !we need to allocate

        allocate(flux_unbinned(0:n_cadj-1),t_unbinned(0:n_cadj-1),z(0:n_cadj-1),k(0:n_cadj-1))

      end if

      k(:) = (/(i, i=0,n_cadj-1, 1)/)
      !Calculate the time-stamps for the binned model
      t_unbinned(:) = xd(j) + t_cad(lj)*((k(:)+1.d0)-0.5d0*(n_cadj+1.d0))/n_cadj

      !Each z is independent for each planet
      call find_z(t_unbinned,pars(0:5,n),z,n_cadj)

      if ( ALL( z > 1.d0 + rps(n*nradius+lj*control) ) .or. rps(n*nradius+lj*control) < small ) then

        flux_out(j) = flux_out(j) * 1.d0

      else

        !Now we have z, let us use Agol's routines
        call occultquad(z,u1(lj),u2(lj),rps(n*nradius+lj*control),flux_unbinned,mu,n_cadj)
        !!!!!call qpower2(z,rp(n),u1,u2,flux_ub(:,n),n_cad)

        !Bin the model if needed
        flux_binned = SUM(flux_unbinned)/n_cadj

	!Compute the final flux
        flux_out(j) = flux_out(j) * flux_binned

      end if

      !The memory is deallocated only if the upcoming point is a new telescope label
      !or if we have reached the number of iterations
      if ( lj .ne. trlab(j+1) .or. j + 1 == n_data ) deallocate(flux_unbinned,t_unbinned,z,k)

    end do

  end do

end subroutine


!This subroutine computes the flux of a star with transiting planets for different bands
!assuming that not a single band requieres binning
subroutine flux_tr_multiband_nobin(xd,trlab,pars,rps,ldc,&
           nbands,nradius,n_data,npl,flux_out)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: n_data, npl, nbands, nradius
  real(kind=mireal), intent(in), dimension(0:n_data-1)  :: xd
  integer, intent(in), dimension(0:n_data-1)  :: trlab !this indicates the instrument label
  real(kind=mireal), intent(in), dimension(0:5,0:npl-1) :: pars
  real(kind=mireal), intent(in), dimension(0:nradius*npl-1) :: rps
!  real(kind=mireal), intent(in), dimension(0:nbands-1,0:npl-1) :: rps
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  real(kind=mireal), intent(in), dimension (0:2*nbands-1) :: ldc
  real(kind=mireal), intent(out), dimension(0:n_data-1) :: flux_out !output flux model
!Local variables
  real(kind=mireal) :: npl_dbl, u1(0:nbands-1), u2(0:nbands-1)
  real(kind=mireal) :: flux_model
  real(kind=mireal) :: mu, z(1)
  integer :: n, j, control, lj

  !This flag controls the multi-radius fits
  control = 1
  if (nradius == 1) control = 0

  npl_dbl = dble(npl)


  do n = 0, nbands - 1
    u1(n) = ldc(2*n)
    u2(n) = ldc(2*n+1)
  end do

  flux_out = 1.d0

  do n = 0, npl - 1

    do j = 0, n_data - 1

      !variable with the current telescope label
      lj = trlab(j)

      !Each z is independent for each planet
      call find_z(xd(j),pars(0:5,n),z,1)

      !Now we have z, let us use Agol's routines
      call occultquad(z,u1(lj),u2(lj),rps(n*nradius+lj*control),flux_model,mu,1)
      !!!!!call qpower2(z,rp(n),u1,u2,flux_ub(:,n),n_cad)

      !Compute the final flux
      flux_out(j) = flux_out(j) * flux_model

    end do

  end do

end subroutine

!Flux for a single band fit
subroutine flux_tr_singleband(xd,pars,rps,ldc,&
           n_cad,t_cad,n_data,npl,flux_out)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: n_data, n_cad, npl
  real(kind=mireal), intent(in), dimension(0:n_data-1)  :: xd
  real(kind=mireal), intent(in), dimension(0:5,0:npl-1) :: pars
  real(kind=mireal), intent(in), dimension(0:npl-1) :: rps
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  real(kind=mireal), intent(in) :: t_cad
  real(kind=mireal), intent(in), dimension (0:1) :: ldc
  real(kind=mireal), intent(out), dimension(0:n_data-1) :: flux_out
!Local variables
  real(kind=mireal), dimension(0:n_cad-1)  :: flux_unbinned, t_unbinned, z, mu
  real(kind=mireal) :: npl_dbl, u1, u2
  real(kind=mireal)   :: flux_binned
  integer :: n, i, j, k(0:n_cad-1)
!External function
  external :: occultquad, find_z

  npl_dbl = dble(npl)

  u1 = ldc(0)
  u2 = ldc(1)

  !Create the jumps
  k(:) = (/(i, i=0,n_cad-1, 1)/)

  !This variable will contain the flux with all the planetary transits
  flux_out = 1.d0

  !Let's do the planets one by one
  do n = 0, npl - 1

    !Let's compute the planet model for the planet n
    do j = 0, n_data - 1

      !Calculate the time-stamps for the binned model
      t_unbinned(:) = xd(j) + t_cad*((k(:)+1.d0)-0.5d0*(n_cad+1.d0))/n_cad
      !Each z is independent for each planet
      call find_z(t_unbinned,pars(0:5,n),z,n_cad)

      !If there is no eclipse then let's avoid call the transit model routines
      if ( ALL( z > 1.d0 + rps(n) ) .or. rps(n) < small ) then

	!There is no transit, so no modification to the flux
        flux_out(j) = flux_out(j) * 1.d0

      else

        !Now we have z, let us use Agol's routines
        call occultquad(z,u1,u2,rps(n),flux_unbinned,mu,n_cad)
        !call qpower2(z,rps(n),u1,u2,flux_ub(:,n),n_cad)

        !Bin the model if needed
        flux_binned = SUM(flux_unbinned)/n_cad

	!Compute the final flux
        flux_out(j) = flux_out(j) * flux_binned

      end if

    end do
    !Now all the model of planet n has been computed

  end do

end subroutine


!Flux for a single band fit, when no binning in needed
subroutine flux_tr_singleband_nobin(xd,pars,rps,ldc,n_data,npl,flux_out)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: n_data, npl
  real(kind=mireal), intent(in), dimension(0:n_data-1)  :: xd
  real(kind=mireal), intent(in), dimension(0:5,0:npl-1) :: pars
  real(kind=mireal), intent(in), dimension(0:npl-1) :: rps
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  real(kind=mireal), intent(in), dimension (0:1) :: ldc
  real(kind=mireal), intent(out), dimension(0:n_data-1) :: flux_out
!Local variables
  real(kind=mireal), dimension(0:n_data-1) :: flux_model, z, mu
  real(kind=mireal) :: npl_dbl, u1, u2
  integer :: n
!External function
  external :: occultquad, find_z

  npl_dbl = dble(npl)

  u1 = ldc(0)
  u2 = ldc(1)

  !This variable will contain the flux with all the planetary transits
  flux_out = 1.d0

  !Let's do the planets one by one
  do n = 0, npl - 1

      !Each z is independent for each planet
      call find_z(xd,pars(0:5,n),z,n_data)

      !Now we have z, let us use Agol's routines
      call occultquad(z,u1,u2,rps(n),flux_model,mu,n_data)
      !call qpower2(z,rps(n),u1,u2,flux_ub(:,n),n_cad)

      !Compute the final flux
      flux_out = flux_out * flux_model

      !Now all the model of planet n has been computed

  end do

end subroutine

subroutine find_chi2_tr(xd,yd,errs,trlab,jtrlab,pars,rps,ldc,jtr,flag, &
           n_cad,t_cad,chi2,n_data,nbands,nradius,njtr,npl)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: n_data, npl, nbands,nradius, njtr
  integer, intent(in) :: n_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension(0:n_data-1)  :: xd, yd, errs
  real(kind=mireal), intent(in), dimension(0:5,0:npl-1) :: pars
  real(kind=mireal), intent(in), dimension(0:nradius*npl-1) :: rps
!  real(kind=mireal), intent(in), dimension(0:nbands-1,0:npl-1) :: rps
  integer, intent(in), dimension(0:n_data-1)  :: trlab, jtrlab
  !pars = T0, P, e, w, b, a/R*
  real(kind=mireal), intent(in) :: t_cad(0:nbands-1)
  real(kind=mireal), dimension(0:njtr-1), intent(in) :: jtr
  logical, intent(in), dimension(0:3) :: flag
  real(kind=mireal), intent(in) :: ldc(0:2*nbands-1)
  real(kind=mireal), intent(out) :: chi2
!Local variables
  real(kind=mireal), dimension(0:n_data-1) :: res

    call find_res_tr(xd,yd,trlab,pars,rps,ldc,flag, &
           n_cad,t_cad,res,n_data,nbands,nradius,npl)
    res(:) = res(:) / sqrt( errs(:)**2 + jtr(jtrlab(:))**2 )
    chi2 = dot_product(res,res)


end subroutine

subroutine find_res_tr(xd,yd,trlab,pars,rps,ldc,flag, &
           n_cad,t_cad,res,n_data,nbands,nradius,npl)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: n_data, npl, nbands, nradius
  integer, intent(in) :: n_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension(0:n_data-1)  :: xd, yd
  real(kind=mireal), intent(in), dimension(0:5,0:npl-1) :: pars
  real(kind=mireal), intent(in), dimension(0:nradius*npl-1) :: rps
!  real(kind=mireal), intent(in), dimension(0:nbands-1,0:npl-1) :: rps
  integer, intent(in), dimension(0:n_data-1)  :: trlab
  !pars = T0, P, e, w, b, a/R*
  real(kind=mireal), intent(in) :: t_cad(0:nbands-1)
  logical, intent(in), dimension(0:3) :: flag
  real(kind=mireal), intent(in) :: ldc(0:2*nbands-1)
  real(kind=mireal), intent(out) :: res(0:n_data-1)
!Local variables
  real(kind=mireal), dimension(0:n_data-1) :: flux_out
  real(kind=mireal), dimension(0:5,0:npl-1) :: up_pars !updated parameters
  real(kind=mireal), dimension(0:npl-1) :: t0, P, e, w, i, a, tp
  real(kind=mireal), dimension (0:2*nbands-1) :: up_ldc
  logical :: is_good
  integer :: n


  t0  = pars(0,:)
  P   = pars(1,:)
  e   = pars(2,:)
  w   = pars(3,:)
  i   = pars(4,:)
  a   = pars(5,:)

  if (flag(0)) P = 1.d0**pars(1,:)
  if (flag(1)) call ewto(e,w,e,w,npl)
  if (flag(3)) call rhotoa(a(0),P(:),a(:),npl)
  if (flag(2)) call btoi(i,a,e,w,i,npl)


  !Update limb darkening coefficients, pass from q's to u's
  is_good = .true.

  do n = 0, nbands - 1
    call get_us(ldc(2*n),ldc(2*n+1),up_ldc(2*n),up_ldc(2*n+1),1)
    call check_us(up_ldc(2*n),up_ldc(2*n+1),is_good)
    if ( .not. is_good ) exit
  end do

  !Avoid solutions with eccentricities larger than 1 and smaller than 0
  if ( any( e > 1.d0 ) .or. any(e < 0.d0 ) ) is_good = .false.

  !Avoid solutions when the planet orbits falls inside the star
  if ( any(a < 1.1) ) is_good = .false.

  if ( is_good ) then

    do n = 0, npl - 1
      call find_tp(t0(n),e(n),w(n),P(n),tp(n))
    end do

    !At this point the parameters to fit are tp,P,e,w,i,a without parametrization
    up_pars(0,:) = tp
    up_pars(1,:) = P
    up_pars(2,:) = e
    up_pars(3,:) = w
    up_pars(4,:) = i
    up_pars(5,:) = a

    !Here we have a vector for the radius called rps

    call flux_tr(xd,trlab,up_pars,rps,up_ldc,&
           n_cad,t_cad,nbands,nradius,n_data,npl,flux_out)
    res(:) =  flux_out(:) - yd(:)

  else

    res(:) = huge(0.d0)

  end if

end subroutine
