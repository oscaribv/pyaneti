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
implicit none

!In/Out variables
  integer, intent(in) :: ts
  double precision, intent(in), dimension(0:ts-1) :: t
  double precision, intent(in), dimension(0:5) :: pars
  double precision, intent(out), dimension(0:ts-1) :: z
!Local variables
  double precision, dimension(0:ts-1) :: ta, swt
  double precision :: tp, P, e, w, i, a
  double precision :: si
!External function
  external :: find_anomaly
!

  tp  = pars(0)
  P   = pars(1)
  e   = pars(2)
  w   = pars(3)
  i   = pars(4)
  a   = pars(5)

  !Obtain the true anomaly by using find_anomaly
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

subroutine flux_tr(xd,pars,ldc,&
           n_cad,t_cad,datas,npl,muld)
implicit none

!In/Out variables
  integer, intent(in) :: datas, n_cad, npl
  double precision, intent(in), dimension(0:datas-1)  :: xd
  double precision, intent(in), dimension(0:6,0:npl-1) :: pars
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  double precision, intent(in) :: t_cad
  double precision, intent(in), dimension (0:1) :: ldc
  double precision, intent(out), dimension(0:datas-1) :: muld
!Local variables
  double precision, dimension(0:datas-1) :: muld_npl
  double precision, dimension(0:datas-1) :: mu
  double precision :: npl_dbl, small, u1, u2, rp(0:npl-1)
  double precision, dimension(0:n_cad-1,0:npl-1)  :: flux_ub
  double precision, dimension(0:n_cad-1)  :: xd_ub, z, fmultip
  integer :: n, j, k(0:n_cad-1)
!External function
  external :: occultquad, find_z

  small = 1.d-5
  npl_dbl = dble(npl)

  u1 = ldc(0)
  u2 = ldc(1)

  !Get planet radius
  rp(:) = pars(6,:)

  do j = 0, n_cad - 1
    k(j) = j
  end do

  muld_npl(:) = 0.d0
  flux_ub(:,:) = 0.d0
  do j = 0, datas - 1

    !Calculate the time-stamps for the binned model
    xd_ub(:) = xd(j) + t_cad*((k(:)+1.d0)-0.5d0*(n_cad+1.d0))/n_cad

    !control the label of the planet
    do n = 0, npl - 1

      !Each z is independent for each planet
      call find_z(xd_ub,pars(0:5,n),z,n_cad)

      if ( ALL( z > 1.d0 + rp(n) ) .or. rp(n) < small ) then

        muld_npl(j) = muld_npl(j) + 1.d0 !This is not eclipse
        flux_ub(:,n) = 0.d0

      else

        !Now we have z, let us use Agol's routines
        call occultquad(z,u1,u2,rp(n),flux_ub(:,n),mu,n_cad)
        !call qpower2(z,rp(n),u1,u2,flux_ub(:,n),n_cad)

      end if

    end do !planets

    fmultip(:) = 0.0
    !Sum the flux of all each sub-division of the model due to each planet
    do n = 0, n_cad - 1
      fmultip(n) =  SUM(flux_ub(n,:))
    end do

    !Re-bin the model
    muld_npl(j) = muld_npl(j) +  sum(fmultip) / n_cad

    !Calcualte the flux received taking into account the transit of all planets
    muld(j) =  1.0d0 + muld_npl(j) - npl_dbl

    !Restart flux_ub
    flux_ub(:,:) = 0.0

  end do !datas

end subroutine


subroutine find_chi2_tr(xd,yd,errs,pars,jitter,flag,ldc,&
           n_cad,t_cad,chi2,datas,npl)
implicit none

!In/Out variables
  integer, intent(in) :: datas, n_cad, npl
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  double precision, intent(in), dimension(0:6,0:npl-1) :: pars
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  double precision, intent(in) :: t_cad
  double precision, intent(in) :: jitter
  logical, intent(in), dimension(0:3) :: flag
  double precision, intent(in), dimension (0:1) :: ldc
  double precision, intent(out) :: chi2
!Local variables
  double precision, dimension(0:datas-1) :: res, muld
  double precision, dimension(0:6,0:npl-1) :: up_pars !updated parameters
  double precision, dimension(0:npl-1) :: t0, P, e, w, i, a, rp, tp
  double precision :: u1, u2, q1k, q2k
  double precision, dimension (0:1) :: up_ldc
  logical :: is_good
  integer :: n
!External function
  external :: flux_tr

  t0  = pars(0,:)
  P   = pars(1,:)
  e   = pars(2,:)
  w   = pars(3,:)
  i   = pars(4,:)
  a   = pars(5,:)
  rp  = pars(6,:)

  if (flag(0)) P = 1.d0**pars(1,:)
  if (flag(1)) call ewto(e,w,e,w,npl)
  if (flag(3)) call rhotoa(a(0),P(:),a(:),npl)
  if (flag(2)) call btoi(i,a,e,w,i,npl)

  !Update limb darkening coefficients, pass from q's to u's
  q1k = ldc(0)
  q2k = ldc(1)
  !re-transform the parameters to u1 and u2
  u1 = sqrt(q1k)
  u2 = u1*( 1.d0 - 2.d0*q2k)
  u1 = 2.d0*u1*q2k
  up_ldc = (/ u1 , u2 /)
 !are the u1 and u2 within a physical solution
  call check_us(u1,u2,is_good)

  if ( any( e > 1.d0 ) ) is_good = .false.

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
    up_pars(6,:) = rp

    call flux_tr(xd,up_pars,up_ldc,n_cad,t_cad,datas,npl,muld)
    res(:) = ( muld(:) - yd(:) ) / sqrt( errs(:)**2 + jitter**2 )
    chi2 = dot_product(res,res)

  else

    chi2 = huge(0.e0)

  end if

end subroutine