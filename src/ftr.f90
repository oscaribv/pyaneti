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

subroutine find_chi2_tr(xd,yd,errs,trlab,jtrlab,pars,rps,ldc,jtr,flag, &
           n_cad,t_cad,chi2,datas,nbands,njtr,npl)
implicit none

!In/Out variables
  integer, intent(in) :: datas, n_cad, npl, nbands, njtr
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  double precision, intent(in), dimension(0:5,0:npl-1) :: pars
  double precision, intent(in), dimension(0:nbands*npl-1) :: rps
!  double precision, intent(in), dimension(0:nbands-1,0:npl-1) :: rps
  integer, intent(in), dimension(0:datas-1)  :: trlab, jtrlab
  !pars = T0, P, e, w, b, a/R*
  double precision, intent(in) :: t_cad
  double precision, dimension(0:njtr-1), intent(in) :: jtr
  logical, intent(in), dimension(0:3) :: flag
  double precision, intent(in) :: ldc(0:2*nbands-1)
  double precision, intent(out) :: chi2
!Local variables
  double precision, dimension(0:datas-1) :: res, muld
  double precision, dimension(0:5,0:npl-1) :: up_pars !updated parameters
  double precision, dimension(0:npl-1) :: t0, P, e, w, i, a, tp
  double precision, dimension(0:nbands-1) :: u1, u2, q1k, q2k
  double precision, dimension (0:2*nbands-1) :: up_ldc
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
  do n = 0, nbands - 1
    q1k(n) = ldc(2*n)
    q2k(n) = ldc(2*n+1)
    u1(n) = 2.d0*q1k(n)*sqrt(q2k(n))
    u2(n) = sqrt(q1k(n))*(1.d0 - 2.d0*q2k(n))
    up_ldc(2*n)   = u1(n)
    up_ldc(2*n+1) = u2(n)
  end do


  is_good = .true.
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

    !Here we have a vector for the radius called rps

    call flux_tr(xd,trlab,up_pars,rps,up_ldc,&
           n_cad,t_cad,nbands,datas,npl,muld)
    res(:) = ( muld(:) - yd(:) ) / sqrt( errs(:)**2 + jtr(jtrlab(:))**2 )
    chi2 = dot_product(res,res)

  else

    chi2 = huge(0.d0)

  end if

end subroutine

subroutine flux_tr(xd,trlab,pars,rps,ldc,&
           n_cad,t_cad,nbands,datas,npl,muld)
implicit none

!In/Out variables
  integer, intent(in) :: datas, n_cad, npl, nbands
  double precision, intent(in), dimension(0:datas-1)  :: xd
  integer, intent(in), dimension(0:datas-1)  :: trlab !this indicates the instrument label
  double precision, intent(in), dimension(0:5,0:npl-1) :: pars
  double precision, intent(in), dimension(0:nbands*npl-1) :: rps
!  double precision, intent(in), dimension(0:nbands-1,0:npl-1) :: rps
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  double precision, intent(in) :: t_cad
  double precision, intent(in), dimension (0:2*nbands-1) :: ldc
  double precision, intent(out), dimension(0:datas-1) :: muld !output flux model
!Local variables
  double precision, dimension(0:datas-1) :: muld_npl
  double precision, dimension(0:datas-1) :: mu
  double precision :: npl_dbl, small, u1(0:nbands-1), u2(0:nbands-1)
  double precision, dimension(0:n_cad-1,0:npl-1)  :: flux_ub
  double precision, dimension(0:n_cad-1)  :: xd_ub, z, fmultip
  integer :: n, j, k(0:n_cad-1)
!External function
  external :: occultquad, find_z

  small = 1.d-5
  npl_dbl = dble(npl)

  do n = 0, nbands - 1
    u1(n) = ldc(2*n)
    u2(n) = ldc(2*n+1)
  end do

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


      if ( ALL( z > 1.d0 + rps(n*nbands+trlab(j)) ) .or. rps(n*nbands+trlab(j)) < small ) then


        muld_npl(j) = muld_npl(j) + 1.d0 !This is not eclipse
        flux_ub(:,n) = 0.d0

      else

        !Now we have z, let us use Agol's routines
        call occultquad(z,u1(trlab(j)),u2(trlab(j)),rps(n*nbands+trlab(j)),flux_ub(:,n),mu,n_cad)
        !!!!!call qpower2(z,rp(n),u1,u2,flux_ub(:,n),n_cad)

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