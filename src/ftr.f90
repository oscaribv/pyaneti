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

subroutine flux_tr(xd,trlab,pars,rps,ldc,&
           n_cad,t_cad,nbands,datas,npl,muld)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: datas, npl, nbands
  integer, intent(in) :: n_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension(0:datas-1)  :: xd
  integer, intent(in), dimension(0:datas-1)  :: trlab !this indicates the instrument label
  real(kind=mireal), intent(in), dimension(0:5,0:npl-1) :: pars
  real(kind=mireal), intent(in), dimension(0:nbands*npl-1) :: rps
!  real(kind=mireal), intent(in), dimension(0:nbands-1,0:npl-1) :: rps
  !pars = T0, P, e, w, b, a/R*, Rp/R*
  real(kind=mireal), intent(in) :: t_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension (0:2*nbands-1) :: ldc
  real(kind=mireal), intent(out), dimension(0:datas-1) :: muld !output flux model
!Local variables
  real(kind=mireal), dimension(0:datas-1) :: muld_npl
  real(kind=mireal), dimension(0:datas-1) :: mu
  real(kind=mireal) :: npl_dbl, u1(0:nbands-1), u2(0:nbands-1)
  real(kind=mireal), allocatable, dimension(:,:)  :: flux_ub
  real(kind=mireal), allocatable, dimension(:)  :: xd_ub, z, fmultip
  integer :: n, j
  integer, allocatable :: k(:)

  npl_dbl = dble(npl)

  do n = 0, nbands - 1
    u1(n) = ldc(2*n)
    u2(n) = ldc(2*n+1)
  end do


  muld_npl(:) = 0.d0
  do j = 0, datas - 1

   allocate (flux_ub(0:n_cad(trlab(j))-1,0:npl-1))
   allocate (xd_ub(0:n_cad(trlab(j))-1),z(0:n_cad(trlab(j))-1),fmultip(0:n_cad(trlab(j))-1))
   allocate (k(0:n_cad(trlab(j))-1))



    do n = 0, n_cad(trlab(j)) - 1
      k(n) = n
    end do


    !Calculate the time-stamps for the binned model
    xd_ub(:) = xd(j) + t_cad(trlab(j))*((k(:)+1.d0)-0.5d0*(n_cad(trlab(j))+1.d0))/n_cad(trlab(j))

    !control the label of the planet
    do n = 0, npl - 1

      !Each z is independent for each planet
      call find_z(xd_ub,pars(0:5,n),z,n_cad(trlab(j)))


      if ( ALL( z > 1.d0 + rps(n*nbands+trlab(j)) ) .or. rps(n*nbands+trlab(j)) < small ) then


        muld_npl(j) = muld_npl(j) + 1.d0 !This is not eclipse
        flux_ub(:,n) = 0.d0

      else

        !Now we have z, let us use Agol's routines
        call occultquad(z,u1(trlab(j)),u2(trlab(j)),rps(n*nbands+trlab(j)),flux_ub(:,n),mu,n_cad(trlab(j)))
        !!!!!call qpower2(z,rp(n),u1,u2,flux_ub(:,n),n_cad)

      end if

    end do !planets

    fmultip(:) = 0.0
    !Sum the flux of all each sub-division of the model due to each planet
    do n = 0, n_cad(trlab(j)) - 1
      fmultip(n) =  SUM(flux_ub(n,:))
    end do

    !Re-bin the model
    muld_npl(j) = muld_npl(j) +  sum(fmultip) / n_cad(trlab(j))

    !Calcualte the flux received taking into account the transit of all planets
    muld(j) =  1.0d0 + muld_npl(j) - npl_dbl

    !Restart flux_ub
    flux_ub(:,:) = 0.0

    deallocate(flux_ub,xd_ub,z,fmultip,k)

  end do !datas

end subroutine

subroutine find_chi2_tr(xd,yd,errs,trlab,jtrlab,pars,rps,ldc,jtr,flag, &
           n_cad,t_cad,chi2,datas,nbands,njtr,npl)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: datas, npl, nbands, njtr
  integer, intent(in) :: n_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension(0:datas-1)  :: xd, yd, errs
  real(kind=mireal), intent(in), dimension(0:5,0:npl-1) :: pars
  real(kind=mireal), intent(in), dimension(0:nbands*npl-1) :: rps
!  real(kind=mireal), intent(in), dimension(0:nbands-1,0:npl-1) :: rps
  integer, intent(in), dimension(0:datas-1)  :: trlab, jtrlab
  !pars = T0, P, e, w, b, a/R*
  real(kind=mireal), intent(in) :: t_cad(0:nbands-1)
  real(kind=mireal), dimension(0:njtr-1), intent(in) :: jtr
  logical, intent(in), dimension(0:3) :: flag
  real(kind=mireal), intent(in) :: ldc(0:2*nbands-1)
  real(kind=mireal), intent(out) :: chi2
!Local variables
  real(kind=mireal), dimension(0:datas-1) :: res

    call find_res_tr(xd,yd,trlab,pars,rps,ldc,flag, &
           n_cad,t_cad,res,datas,nbands,npl)
    res(:) = res(:) / sqrt( errs(:)**2 + jtr(jtrlab(:))**2 )
    chi2 = dot_product(res,res)


end subroutine

subroutine find_res_tr(xd,yd,trlab,pars,rps,ldc,flag, &
           n_cad,t_cad,res,datas,nbands,npl)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: datas, npl, nbands
  integer, intent(in) :: n_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension(0:datas-1)  :: xd, yd
  real(kind=mireal), intent(in), dimension(0:5,0:npl-1) :: pars
  real(kind=mireal), intent(in), dimension(0:nbands*npl-1) :: rps
!  real(kind=mireal), intent(in), dimension(0:nbands-1,0:npl-1) :: rps
  integer, intent(in), dimension(0:datas-1)  :: trlab
  !pars = T0, P, e, w, b, a/R*
  real(kind=mireal), intent(in) :: t_cad(0:nbands-1)
  logical, intent(in), dimension(0:3) :: flag
  real(kind=mireal), intent(in) :: ldc(0:2*nbands-1)
  real(kind=mireal), intent(out) :: res(0:datas-1)
!Local variables
  real(kind=mireal), dimension(0:datas-1) :: muld
  real(kind=mireal), dimension(0:5,0:npl-1) :: up_pars !updated parameters
  real(kind=mireal), dimension(0:npl-1) :: t0, P, e, w, i, a, tp
  real(kind=mireal), dimension(0:nbands-1) :: u1, u2, q1k, q2k
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
    q1k(n) = ldc(2*n)
    q2k(n) = ldc(2*n+1)
    u1(n) = 2.d0*q1k(n)*sqrt(q2k(n))
    u2(n) = sqrt(q1k(n))*(1.d0 - 2.d0*q2k(n))
    up_ldc(2*n)   = u1(n)
    up_ldc(2*n+1) = u2(n)
    call check_us(u1(n),u2(n),is_good)
    if ( .not. is_good ) exit
  end do

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
    res(:) =  muld(:) - yd(:)

  else

    res(:) = huge(0.d0)

  end if

end subroutine