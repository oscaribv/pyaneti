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
subroutine find_z(t,pars,flag,z,ts)
implicit none

!In/Out variables
  integer, intent(in) :: ts
  double precision, intent(in), dimension(0:ts-1) :: t
  double precision, intent(in), dimension(0:5) :: pars
  double precision, intent(out), dimension(0:ts-1) :: z
  logical, intent(in), dimension(0:3) :: flag 
!Local variables
  double precision, dimension(0:ts-1) :: ta, swt
  double precision :: t0, P, e, w, i, a
  double precision :: wp, ws, si
  double precision :: pi = 3.1415926535897932384626d0
  integer :: n, j
!External function
  external :: find_anomaly
!

  t0  = pars(0)
  P   = pars(1)
  e   = pars(2)
  w   = pars(3)
  i   = pars(4)
  a   = pars(5)

  if ( flag(0) ) P = 1.d0**pars(1)
  if ( flag(1) ) then
    e = pars(2) * pars(2) + pars(3) * pars(3)
    w = atan2(pars(2),pars(3))
  end if
  ws  = w       !star periastron
  wp  = ws + pi !planet periastron
  if (flag(3)) a = 10.d0**a
  if (flag(2)) i = acos( i / a * ( 1.d0 + e * sin(wp) ) / ( 1.d0 - e*e ) )

  !Obtain the true anomaly by using find_anomaly
  call find_anomaly(t,t0,e,ws,P,ta,ts)

  swt = sin(wp+ta)
  si = sin(i)

  where (swt < 0.0 ) !We have the planet in front of the star -> transit
  !z has been calculated
    z = a * ( 1.d0 - e * e ) * sqrt( 1.d0 - swt * swt * si * si ) &
        / ( 1.d0 + e * cos(ta) )
  elsewhere !We have the planet behind the star -> occulation
  !z has avalue which gives flux = 1 in Mandel & Agol module
      z = 1.d1
  end where

end subroutine

subroutine flux_tr(xd,pars,jitter,flag,ldc,&
           n_cad,t_cad,datas,npl,muld)
implicit none

!In/Out variables
  integer, intent(in) :: datas, n_cad, npl
  double precision, intent(in), dimension(0:datas-1)  :: xd
  double precision, intent(in), dimension(0:6,0:npl-1) :: pars
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  double precision, intent(in) :: t_cad
  double precision, intent(in) :: jitter
  logical, intent(in), dimension(0:3) :: flag
  double precision, intent(in), dimension (0:1) :: ldc
  double precision, intent(out), dimension(0:datas-1) :: muld
!Local variables
  double precision, dimension(0:datas-1) :: muld_npl
  double precision, dimension(0:datas-1) :: res, mu
  double precision :: npl_dbl, small, q1k, q2k, u1, u2, pz(0:npl-1)
  double precision, dimension(0:n_cad-1,0:npl-1)  :: flux_ub
  double precision, dimension(0:n_cad-1)  :: xd_ub, z, fmultip
  integer :: n, j, k
  logical :: is_good
!External function
  external :: occultquad, find_z

  small = 1.d-5
  npl_dbl = dble(npl)

  q1k = ldc(0)
  q2k = ldc(1)
  !re-transform the parameters to u1 and u2
  u1 = sqrt(q1k)
  u2 = u1*( 1.d0 - 2.d0*q2k)
  u1 = 2.d0*u1*q2k

  !Get planet radius
  pz(:) = pars(6,:)

  !Selective re-sampling
  n = 0
  muld_npl(:) = 0.d0
  do j = 0, datas - 1

    !Calculate the time-stamps for the binned model
    do k = 0, n_cad - 1
      xd_ub(k) = xd(j) + t_cad*((k+1.d0)-0.5d0*(n_cad+1.d0))/n_cad
    end do

    !control the label of the planet
    do n = 0, npl - 1

      !Each z is independent for each planet
      call find_z(xd_ub,pars(0:5,n),flag,z,n_cad)

      if ( ALL( z > 1.d0 + pz(n) ) .or. pz(n) < small ) then

        muld_npl(j) = muld_npl(j) + 1.d0 !This is not eclipse

      else

        !Now we have z, let us use Agol's routines
        call occultquad(z,u1,u2,pz(n),flux_ub(:,n),mu,n_cad)

      end if

    end do !planets

    fmultip(:) = 0.0
    !Sum the flux of all each sub-division of the model due to each planet
    do k = 0, n_cad - 1
      fmultip(k) =  SUM(flux_ub(k,:))
    end do

    !Re-bin the model
    muld_npl(j) = muld_npl(j) +  sum(fmultip) / n_cad

    !Calcualte the flux received taking into account the transit of all planets
    muld(j) =  1.0d0 + muld_npl(j) - npl_dbl

    !Restart flux_ub
    flux_ub(:,:) = 0.0

  end do !datas

end subroutine

subroutine find_chi2_tr(xd,yd,errs,plab_tr,pars,jitter,flag,ldc,&
           n_cad,t_cad,chi2,datas,npl)
implicit none

!In/Out variables
  integer, intent(in) :: datas, n_cad, npl
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1)  :: plab_tr
  double precision, intent(in), dimension(0:6,0:npl-1) :: pars
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  double precision, intent(in) :: t_cad
  double precision, intent(in) :: jitter
  logical, intent(in), dimension(0:3) :: flag
  double precision, intent(in), dimension (0:1) :: ldc
  double precision, intent(out) :: chi2
!Local variables
  double precision, dimension(0:datas-1) :: res, muld
  double precision :: npl_dbl, small, u1, u2, pz(0:npl-1), q1k, q2k
  logical :: is_good
  integer :: n
!External function
  external :: occultquad, find_z, flux_tr

  small = 1.d-5
  npl_dbl = dble(npl)

  q1k = ldc(0)
  q2k = ldc(1)
  !re-transform the parameters to u1 and u2
  u1 = sqrt(q1k)
  u2 = u1*( 1.d0 - 2.d0*q2k)
  u1 = 2.d0*u1*q2k

!are the u1 and u2 within a physical solution
  call check_us(u1,u2,is_good)
  if ( flag(1) ) then
    do n = 0, npl - 1
     call check_e(pars(2,n),pars(3,n),is_good)
     if ( .not. is_good ) exit
    end do
  end if

  if ( is_good ) then

  call flux_tr(xd,pars,jitter,flag,ldc,n_cad,t_cad,datas,npl,muld)

  !Let us calculate the residuals
  ! chi^2 = \Sum_i ( M - O )^2 / \sigma^2
  !Here I am assuming that we want limb darkening
  !If this is not true, use mu
  res(:) = ( muld(:) - yd(:) ) / sqrt( errs(:)**2 + jitter**2 )
  chi2 = dot_product(res,res)

  else

    chi2 = huge(0.e0)

  end if

end subroutine
