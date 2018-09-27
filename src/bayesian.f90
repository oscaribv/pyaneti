subroutine get_loglike(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           tlab,jrvlab,tff,flags,&
           t_cad,n_cad,pars,rvs,ldc,trends,jrv,jtr, &
           loglike,chi2_rv,chi2_tr,npl,n_tel,n_jrv,size_rv,size_tr)
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npl, n_tel,n_jrv,n_cad !size of RV and LC data
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: tlab, jrvlab
!pars = T0, P, e, w, b, a/R*, Rp/R*, K -> for each planet
  double precision, intent(in) :: pars(0:8*npl-1), rvs(0:n_tel-1), ldc(0:1)
  double precision, intent(in) :: t_cad
  double precision, intent(in) :: trends(0:1)
  double precision, dimension(0:n_jrv-1), intent(in) :: jrv
  double precision, intent(in) :: jtr
  logical, intent(in) :: flags(0:5)
  logical, intent(in) :: tff(0:1) !total_fit_flag
  double precision, intent(out) :: loglike, chi2_rv, chi2_tr
!Local variables
  double precision :: chi2_total
  double precision :: log_errs
  double precision :: two_pi = 2.d0*3.1415926535897932384626d0
  integer :: m
  external:: get_total_chi2

  !Calcualte the chi2
  call get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           tlab,jrvlab,tff,flags,&
           t_cad,n_cad,pars,rvs,ldc,trends,jrv,jtr, &
           chi2_rv,chi2_tr,npl,n_tel,n_jrv,size_rv,size_tr)


  !Calculate the normalization term
    log_errs = 0.0
    if ( tff(0)  .and. size_rv > 1) then
      do m = 0, size_rv - 1
        log_errs = log_errs + &
        log( 1.0d0/sqrt( two_pi * ( e_rv(m)**2 + jrv(jrvlab(m))**2 ) ) )
      end do
    else
      chi2_rv = 0.d0
    end if

    if ( tff(1) .and. size_tr > 1 ) then
      log_errs = log_errs + &
      sum(log( 1.0d0/sqrt( two_pi * ( e_tr(:)**2 + jtr**2 ) ) ) )
    else
      chi2_tr = 0.d0
    end if

    chi2_total = chi2_rv + chi2_tr

    loglike = log_errs - 0.5d0 * chi2_total

end subroutine


subroutine get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           tlab,jrvlab,tff,flags,&
           t_cad,n_cad,pars,rvs,ldc,trends,jrv,jtr, &
           chi2_rv,chi2_tr,npl,n_tel,n_jrv,size_rv,size_tr)
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npl, n_tel,n_jrv,n_cad !size of RV and LC data
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: tlab, jrvlab
!pars = T0, P, e, w, b, a/R*, Rp/R*, K -> for each planet
  double precision, intent(in) :: pars(0:8*npl-1), rvs(0:n_tel-1), ldc(0:1)
  double precision, intent(in) :: t_cad
  double precision, intent(in) :: trends(0:1)
  double precision, dimension(0:n_jrv-1), intent(in) :: jrv
  double precision, intent(in) :: jtr
  logical, intent(in) :: flags(0:5)
  logical, intent(in) :: tff(0:1) !total_fit_flag
  double precision, intent(out) :: chi2_rv, chi2_tr
!Local variables
  double precision :: pars_rv(0:7+n_tel-1,0:npl-1)
  double precision :: pars_tr(0:6,0:npl-1)
  logical :: flag_rv(0:3), flag_tr(0:3)
  integer :: i, j

  !Create the parameter variables for rv and tr
  do i = 0, npl - 1
    j = 8*i !the new data start point
    pars_tr(:,i)   = pars(j:j+6)
    pars_rv(0:3,i) = pars(j:j+3)
    pars_rv(4,i) = pars(j+7)
    pars_rv(5:6,i) = trends(:)
    pars_rv(7:7+n_tel-1,i) = rvs(:)
  end do

  !Put the correct flags
  flag_tr(:)   = flags(0:3)
  flag_rv(0:1) = flags(0:1)
  flag_rv(2:3) = flags(4:5)

  !Let us calculate chi2
  chi2_rv = 0.d0
  chi2_tr = 0.d0

  if (tff(1) ) &
  call find_chi2_tr(x_tr,y_tr,e_tr,pars_tr,jtr,flag_tr,&
                        ldc,n_cad,t_cad,chi2_tr,size_tr,npl)
  if (tff(0) ) &
  call find_chi2_rv(x_rv,y_rv,e_rv,tlab,jrvlab,pars_rv,jrv,&
                    flag_rv,chi2_rv,size_rv,n_tel,n_jrv,npl)


end subroutine

subroutine gauss_prior(mu,sigma,x,prob)
implicit none

  !In/Out variables
  double precision, intent(in) :: mu, sigma, x
  double precision, intent(out)  :: prob
  !Local variables
  double precision  :: sigma2
  double precision  :: two_pi = 2.d0*3.1415926535897932384626d0

  sigma2 = sigma*sigma
  prob = sqrt(two_pi*sigma2)
  prob = exp(- 0.5d0 * (x - mu)**2 / sigma2) / prob

end subroutine


subroutine uniform_prior(lx,rx,x,prob)
implicit none

  !In/Out variables
  double precision, intent(in) :: lx, rx, x
  double precision, intent(out)  :: prob

  prob = 1.d0/abs(rx - lx)
  if ( x < lx .or. x > rx ) prob = 0.d0

end subroutine

subroutine get_priors(fit_pars,lims,pars_in,priors_out,npars)
implicit none

  integer, intent(in) :: npars
  double precision, intent(in), dimension(0:2*npars-1) :: lims
  double precision, intent(in), dimension(0:npars-1) :: pars_in
  double precision, intent(out), dimension(0:npars-1) :: priors_out
  character, intent(in), dimension(0:npars-1) :: fit_pars
!Local
  integer :: j

  priors_out = 1.d0
  do j = 0, npars - 1
    if ( fit_pars(j) == 'u' ) then
      call uniform_prior(lims(2*j),lims(2*j+1),pars_in(j),priors_out(j))
      !end if
    else if ( fit_pars(j) == 'g' ) then
      call gauss_prior(lims(2*j),lims(2*j+1),pars_in(j),priors_out(j))
    end if
  end do

end subroutine
