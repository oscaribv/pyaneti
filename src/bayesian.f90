subroutine get_loglike(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           rvlab,jrvlab,trlab,jtrlab,tff,flags,&
           pars,model_int,model_double,&
           npars,loglike,chi2_rv,chi2_tr,size_rv,size_tr)
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npars !size of RV and LC data
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: rvlab, jrvlab
  integer, intent(in), dimension(0:size_tr-1) :: trlab, jtrlab
!npars = 7*npl + (npl + LDC)*nbands + noffsets + njitter + ntrends
  double precision, intent(in) :: pars(0:npars-1)
  double precision, intent(in) :: model_double(0:0)
  integer, intent(in) :: model_int(0:6)
  logical, intent(in) :: flags(0:5)
  logical, intent(in) :: tff(0:1) !total_fit_flag
  double precision, intent(out) :: loglike, chi2_rv, chi2_tr
!Local variables
  !Model int variables
  integer :: npl, nbands, ntels, nldc, njrv, njtr, n_cad
  !Model double variables
  double precision :: t_cad
  !
  double precision :: chi2_total
  double precision :: log_errs
  external:: get_total_chi2

  !Integer variables
  npl    = model_int(0)
  ntels  = model_int(1)
  nbands = model_int(2)
  nldc   = model_int(3)
  njrv   = model_int(4)
  njtr   = model_int(5)
  n_cad  = model_int(6)
  !Double variables
  t_cad  = model_double(0)

  !Calculate the chi2
   call get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
                           rvlab,jrvlab,trlab,jtrlab,tff,flags,&
                           t_cad,n_cad,pars,chi2_rv,chi2_tr,log_errs,npars,&
                           npl,ntels,nbands,nldc,njrv,njtr,size_rv,size_tr)

    chi2_total = chi2_rv + chi2_tr

    loglike = log_errs - 0.5d0 * chi2_total

end subroutine


subroutine get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           rvlab,jrvlab,trlab,jtrlab,tff,flags,&
           t_cad,n_cad,pars, &
           chi2_rv,chi2_tr,log_errs,npars,npl,ntels,nbands,nldc,njrv,njtr,size_rv,size_tr)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr,npars, npl, ntels,nbands,nldc,njrv,njtr,n_cad !size of RV and LC data
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: rvlab, jrvlab
  integer, intent(in), dimension(0:size_tr-1) :: trlab, jtrlab
!pars = T0, P, e, w, b, a/R*, Rp/R*, K -> for each planet
  double precision, intent(in) :: t_cad
  double precision, intent(in) :: pars(0:npars-1)
  logical, intent(in) :: flags(0:5)
  logical, intent(in) :: tff(0:1) !total_fit_flag
  double precision, intent(out) :: chi2_rv, chi2_tr, log_errs
!Local variables
  double precision :: ldc(0:nldc*nbands-1)
  double precision, dimension(0:njrv-1) :: jrv
  double precision, dimension(0:njtr-1) :: jtr
  double precision :: pars_rv(0:7+ntels-1,0:npl-1)
  double precision :: pars_tr(0:5,0:npl-1), rps(0:nbands*npl-1)
  logical :: flag_rv(0:3), flag_tr(0:3)
  integer :: i, j
  integer :: srp, sldc, srv, sjitrv, sjittr, strends

  srp = 7*npl                !start of planet radius
  sldc = srp + npl*nbands    !start of LDC
  srv = sldc + nldc*nbands   !start of RV offsets
  sjitrv = srv + ntels       !start of RV jitter
  sjittr = sjitrv + njrv     !start of TR jitter
  strends = sjittr + njtr    !start of RV trends

  !Create the parameter variables for rv and tr
  do i = 0, npl - 1
    j = 7*i !the new data start point
    pars_tr(0:5,i)   = pars(j:j+5)
    rps(i*nbands:(i+1)*nbands-1) = pars(srp+nbands*i:srp+nbands*(i+1)-1)
    pars_rv(0:3,i)   = pars(j:j+3)
    pars_rv(4,i)     = pars(j+6)
    pars_rv(5:6,i)   = pars(strends:strends+1)
    pars_rv(7:7+ntels-1,i) = pars(srv:srv+ntels-1)
  end do

   ldc(:) = pars(sldc:sldc+nldc*nbands-1)
   jrv(:) = pars(sjitrv:sjitrv+njrv-1)
   jtr(:) = pars(sjittr:sjittr+njtr-1)

  !Put the correct flags
  flag_tr(:)   = flags(0:3)
  flag_rv(0:1) = flags(0:1)
  flag_rv(2:3) = flags(4:5)

  !Let us calculate chi2
  chi2_rv = 0.d0
  chi2_tr = 0.d0

  if (tff(1) ) &
  call find_chi2_tr(x_tr,y_tr,e_tr,trlab,jtrlab,pars_tr,rps,ldc,jtr,flag_tr, &
           n_cad,t_cad,chi2_tr,size_tr,nbands,njtr,npl)
  if (tff(0) ) &
  call find_chi2_rv(x_rv,y_rv,e_rv,rvlab,jrvlab,pars_rv,jrv,&
                    flag_rv,chi2_rv,size_rv,ntels,njrv,npl)


  !Calculate the normalization term
    log_errs = 0.0
    if ( tff(0)  .and. size_rv > 1) then
      log_errs = log_errs + sum( log( 1.0d0/sqrt( two_pi * ( e_rv(:)**2 + jrv(jrvlab(:))**2 ) ) ) )
    else
      chi2_rv = 0.d0
    end if

    if ( tff(1) .and. size_tr > 1 ) then
      log_errs = log_errs + sum( log( 1.0d0/sqrt( two_pi * ( e_tr(:)**2 + jtr(jtrlab(:))**2 ) ) ) )
    else
      chi2_tr = 0.d0
    end if


end subroutine

subroutine gauss_prior(mu,sigma,x,prob)
use constants
implicit none

  !In/Out variables
  double precision, intent(in) :: mu, sigma, x
  double precision, intent(out)  :: prob
  !Local variables
  double precision  :: sigma2

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
