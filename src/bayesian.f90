subroutine get_loglike(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           rvlab,jrvlab,trlab,jtrlab,tff,flags,kernels,&
           pars,model_int,model_double,nmodel_int,nmodel_double,&
           npars,loglike,chi2_rv,chi2_tr,size_rv,size_tr)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npars, nmodel_int, nmodel_double !size of RV and LC data
  real(kind=mireal), intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  real(kind=mireal), intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: rvlab, jrvlab
  integer, intent(in), dimension(0:size_tr-1) :: trlab, jtrlab
  character(len=6), intent(in) :: kernels
  !npars = 7*npl + (npl + LDC)*nbands + noffsets + njitter + ntrends
  real(kind=mireal), intent(in) :: pars(0:npars-1)
  real(kind=mireal), intent(in) :: model_double(0:nmodel_double-1)
  integer, intent(in) :: model_int(0:nmodel_int-1)
  logical, intent(in) :: flags(0:5)
  logical, intent(in) :: tff(0:1) !total_fit_flag
  real(kind=mireal), intent(out) :: loglike, chi2_rv, chi2_tr
  !Local variables
  !Model int variables
  integer :: npl, nbands,nradius, ntels, nldc, njrv, njtr, n_cad(0:nmodel_int-8-1), np_rv, np_tr
  !Model double variables
  real(kind=mireal) :: t_cad(0:nmodel_double-1)
  !
  real(kind=mireal) :: chi2_total
  real(kind=mireal) :: log_errs
  character(len=3) :: kernel_rv, kernel_tr

  !Integer variables
  npl    = model_int(0)
  ntels  = model_int(1)
  nbands = model_int(2)
  nradius= model_int(3)
  nldc   = model_int(4)
  njrv   = model_int(5)
  njtr   = model_int(6)
  np_rv  = model_int(7) !number of kernel parameters for rv fit
  np_tr  = model_int(8) !number of kernel parameters for tr fit
  n_cad(:)  = model_int(9:9+nbands-1)
  !Double variables
  t_cad(:)  = model_double(:)
  !
  kernel_rv = kernels(1:3)
  kernel_tr = kernels(4:6)

  if ( kernel_rv == 'Non' .and. kernel_tr == 'Non' ) then

  !We are only working with white noise
  !Let us compute the chi2
  call get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
                           rvlab,jrvlab,trlab,jtrlab,tff,flags,&
                           t_cad,n_cad,pars,chi2_rv,chi2_tr,log_errs,npars,&
                           npl,ntels,nbands,nradius,nldc,njrv,njtr,size_rv,size_tr)

    chi2_total = chi2_rv + chi2_tr

    loglike = log_errs - 0.5d0 * chi2_total

 else

  !We are working with correlated noise, let's do GP!
   call get_total_GP_likelihood(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
        rvlab,jrvlab,trlab,jtrlab,tff,flags,t_cad,n_cad,pars,kernels, &
        chi2_rv,chi2_tr,loglike,npars,npl,ntels,nbands,nradius,nldc,njrv,njtr,np_rv,&
        np_tr,size_rv,size_tr)

  end if


end subroutine

!--------------------------------------------------------------------------------------
!This subroutine calculates the Likelihood for correlated noise in RV and/or TR data
!--------------------------------------------------------------------------------------
subroutine get_total_GP_likelihood(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           rvlab,jrvlab,trlab,jtrlab,tff,flags,t_cad,n_cad,pars,kernels, &
       chi2_rv,chi2_tr,log_like_total,npars,npl,ntels,nbands,nradius,nldc,njrv,njtr,&
           np_rv,np_tr,size_rv,size_tr)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr,npars, npl, ntels,nbands,nradius,nldc,njrv,njtr, np_rv, np_tr !size of RV and LC data
  integer, intent(in) :: n_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  real(kind=mireal), intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: rvlab, jrvlab
  integer, intent(in), dimension(0:size_tr-1) :: trlab, jtrlab
  character(len=6), intent(in) :: kernels
!pars = T0, P, e, w, b, a/R*, Rp/R*, K -> for each planet
  real(kind=mireal), intent(in) :: t_cad(0:nbands-1)
  real(kind=mireal), intent(in) :: pars(0:npars-1)
  logical, intent(in) :: flags(0:5)
  logical, intent(in) :: tff(0:1) !total_fit_flag
  real(kind=mireal), intent(out) :: chi2_rv, chi2_tr, log_like_total
!Local variables
  real(kind=mireal) :: ldc(0:nldc*nbands-1)
  real(kind=mireal), dimension(0:njrv-1) :: jrv
  real(kind=mireal), dimension(0:njtr-1) :: jtr
  real(kind=mireal) :: pars_rv(0:7+ntels-1,0:npl-1)
  real(kind=mireal) :: pars_tr(0:5,0:npl-1), rps(0:nradius*npl-1)
  real(kind=mireal), dimension(0:size_rv-1) :: res_rv
  real(kind=mireal), dimension(0:size_tr-1) :: res_tr
  real(kind=mireal), dimension(0:np_rv-1) :: pk_rv
  real(kind=mireal), dimension(0:np_tr-1) :: pk_tr
  real(kind=mireal) :: log_errs, nll_rv, nll_tr
  logical :: flag_rv(0:3), flag_tr(0:3)
  integer :: i, j, mk
  character(len=3) :: kernel_rv, kernel_tr
  integer :: srp, sldc, srv, sjitrv, sjittr, strends, skrv, sktr

  srp = 7*npl                !start of planet radius
  sldc = srp + npl*nradius    !start of LDC
  srv = sldc + nldc*nbands   !start of RV offsets
  sjitrv = srv + ntels       !start of RV jitter
  sjittr = sjitrv + njrv     !start of TR jitter
  strends = sjittr + njtr    !start of RV trends
  skrv = strends + 2         !start of RV kernel parameters
  sktr = skrv + np_rv        !start of TR kernel parameters

  !Create the parameter variables for rv and tr
  do i = 0, npl - 1
    j = 7*i !the new data start point
    pars_tr(0:5,i)   = pars(j:j+5)
    rps(i*nradius:(i+1)*nradius-1) = pars(srp+nradius*i:srp+nradius*(i+1)-1)
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

  kernel_rv = kernels(1:3)
  kernel_tr = kernels(4:6)

  pk_rv = pars(skrv:skrv+np_rv)
  pk_tr = pars(sktr:sktr+np_tr)

  !Let us calculate chi2
  chi2_rv = 0.d0
  chi2_tr = 0.d0

  log_errs = 0.0
  if (tff(0) ) then

    if (kernel_rv == 'Non' ) then
      call find_chi2_rv(x_rv,y_rv,e_rv,rvlab,jrvlab,pars_rv,jrv,&
                        flag_rv,chi2_rv,size_rv,ntels,njrv,npl)
      !Calculate the normalization term
      log_errs = log_errs + sum( log( 1.0d0/sqrt( two_pi * ( e_rv(:)**2 + jrv(jrvlab(:))**2 ) ) ) )
      !chi2_rv = 0.d0
      nll_rv = log_errs - 0.5 * chi2_rv

    else if ( kernel_rv(1:2) == 'MQ' .or. kernel_rv(1:2) == 'EX' .or. kernel_rv == 'MF' ) then

      !Check how many time-series we have
      !Solution seen at https://stackoverflow.com/questions/9900417/character-to-integer-conversion-in-fortran
      read(kernel_rv(3:3),'(i1)') mk

      !This is the V. Rajpaul Framework
      !Find the residuals of the RV time-series, this takes into account all the modelled planets
      !And all the RV offsets
      res_rv = 0.0
      call find_res_rv(x_rv(0:size_rv/mk-1),y_rv(0:size_rv/mk-1),rvlab(0:size_rv/mk-1),&
                     pars_rv,flag_rv,res_rv(0:size_rv/mk-1),size_rv/mk,ntels,npl)
      !Find the residuals of the Activity indicators time-series, this takes into account aonly offsets
      do i = 1, mk - 1
          res_rv(i*size_rv/mk:(i+1)*size_rv/mk-1) = y_rv(i*size_rv/mk:(i+1)*size_rv/mk-1) &
                                                          - pars(srv+rvlab(i*size_rv/mk:(i+1)*size_rv/mk-1))

      end do

      !Compute the log likelihood using the multi-dimensional GP Regression
      call NLL_GP(pk_rv,kernel_rv,x_rv,res_rv,e_rv,jrv,jrvlab,nll_rv,chi2_rv,np_rv,size_rv,njrv)

    else

      !RV GP
      call find_res_rv(x_rv,y_rv,rvlab,pars_rv,flag_rv,res_rv,size_rv,ntels,npl)
      call NLL_GP(pk_rv,kernel_rv,x_rv,res_rv,e_rv,jrv,jrvlab,nll_rv,chi2_rv,np_rv,size_rv,njrv)

    end if

  else

  !There is no RV fit
  nll_rv = 0.0

  end if

  if (tff(1) ) then

    if (kernel_tr == 'Non' ) then

      call find_chi2_tr(x_tr,y_tr,e_tr,trlab,jtrlab,pars_tr,rps,ldc,jtr,flag_tr, &
           n_cad,t_cad,chi2_tr,size_tr,nbands,nradius,njtr,npl)
      log_errs =  sum( log( 1.0d0/sqrt( two_pi * ( e_tr(:)**2 + jtr(jtrlab(:))**2 ) ) ) )
      nll_tr = log_errs - 0.5 * chi2_tr

    else

      !TR GP
      call find_res_tr(x_tr,y_tr,trlab,pars_tr,rps,ldc,flag_tr,n_cad,t_cad,res_tr,size_tr,nbands,nradius,npl)
      call NLL_GP(pk_tr,kernel_tr,x_tr,res_tr,e_tr,jtr,jtrlab,nll_tr,chi2_tr,np_tr,size_tr,njtr)

    end if

  else

    !There is no TR fit
    nll_tr = 0.0

  end if

  log_like_total = nll_rv + nll_tr

end subroutine


subroutine get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           rvlab,jrvlab,trlab,jtrlab,tff,flags,&
           t_cad,n_cad,pars, &
           chi2_rv,chi2_tr,log_errs,npars,npl,ntels,nbands,nradius,nldc,njrv,njtr,size_rv,size_tr)
use constants
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr,npars, npl, ntels,nbands,nradius,nldc,njrv,njtr !size of RV and LC data
  integer, intent(in) :: n_cad(0:nbands-1)
  real(kind=mireal), intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  real(kind=mireal), intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: rvlab, jrvlab
  integer, intent(in), dimension(0:size_tr-1) :: trlab, jtrlab
!pars = T0, P, e, w, b, a/R*, Rp/R*, K -> for each planet
  real(kind=mireal), intent(in) :: t_cad(0:nbands-1)
  real(kind=mireal), intent(in) :: pars(0:npars-1)
  logical, intent(in) :: flags(0:5)
  logical, intent(in) :: tff(0:1) !total_fit_flag
  real(kind=mireal), intent(out) :: chi2_rv, chi2_tr, log_errs
!Local variables
  real(kind=mireal) :: ldc(0:nldc*nbands-1)
  real(kind=mireal), dimension(0:njrv-1) :: jrv
  real(kind=mireal), dimension(0:njtr-1) :: jtr
  real(kind=mireal) :: pars_rv(0:7+ntels-1,0:npl-1)
  real(kind=mireal) :: pars_tr(0:5,0:npl-1), rps(0:nradius*npl-1)
  logical :: flag_rv(0:3), flag_tr(0:3)
  integer :: i, j
  integer :: srp, sldc, srv, sjitrv, sjittr, strends

  srp = 7*npl                !start of planet radius
  sldc = srp + npl*nradius   !start of LDC
  srv = sldc + nldc*nbands   !start of RV offsets
  sjitrv = srv + ntels       !start of RV jitter
  sjittr = sjitrv + njrv     !start of TR jitter
  strends = sjittr + njtr    !start of RV trends

  !Create the parameter variables for rv and tr
  do i = 0, npl - 1
    j = 7*i !the new data start point
    pars_tr(0:5,i)   = pars(j:j+5)
    rps(i*nradius:(i+1)*nradius-1) = pars(srp+nradius*i:srp+nradius*(i+1)-1)
    pars_rv(0:3,i)   = pars(j:j+3)
    pars_rv(4,i)     = pars(j+6)
    pars_rv(5:6,i)   = pars(strends:strends+1)
    pars_rv(7:7+ntels-1,i) = pars(srv:srv+ntels-1)
  end do

   !print *, nradius
   !print *, rps(0)

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
           n_cad,t_cad,chi2_tr,size_tr,nbands,nradius,njtr,npl)
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
  real(kind=mireal), intent(in) :: mu, sigma, x
  real(kind=mireal), intent(out)  :: prob
  !Local variables
  real(kind=mireal)  :: sigma2

  sigma2 = sigma*sigma
  prob = sqrt(two_pi*sigma2)
  prob = exp(- 0.5d0 * (x - mu)**2 / sigma2) / prob

end subroutine


subroutine uniform_prior(lx,rx,x,prob)
use constants
implicit none

  !In/Out variables
  real(kind=mireal), intent(in) :: lx, rx, x
  real(kind=mireal), intent(out)  :: prob

  prob = 1.d0/abs(rx - lx)
  if ( x < lx .or. x > rx ) prob = 0.d0

end subroutine


subroutine uniform_prior_constant(lx,rx,x,prob)
use constants
implicit none

  !In/Out variables
  real(kind=mireal), intent(in) :: lx, rx, x
  real(kind=mireal), intent(out)  :: prob

  prob = 1.d0
  if ( x < lx .or. x > rx ) prob = 0.d0

end subroutine

!jeffreys_prior between [lx,rx]
subroutine jeffreys_prior(lx,rx,x,prob)
use constants
implicit none

  !In/Out variables
  real(kind=mireal), intent(in) :: lx, rx, x
  real(kind=mireal), intent(out)  :: prob

  prob = 1./x/log(rx/lx)
  if (x < 0) prob = 0.0

end subroutine

!Modifed version of Jeffreys priors as given by Gregory 2005
!https://iopscience.iop.org/article/10.1086/432594/pdf
subroutine modjeffreys_prior(lx,rx,x,prob)
use constants
implicit none

  !In/Out variables
  real(kind=mireal), intent(in) :: lx, rx, x
  real(kind=mireal), intent(out)  :: prob

  prob = 1.d0/(x+lx)/log((lx+rx)/lx)
  if (x < 0 .or. x > rx) prob = 0.0


end subroutine

!Modifed version of Jeffreys priors as given by Gregory 2005
!This subroutine allows for negative values
!https://iopscience.iop.org/article/10.1086/432594/pdf
subroutine negmodjeffreys_prior(lx,rx,x,prob)
use constants
implicit none

  !In/Out variables
  real(kind=mireal), intent(in) :: lx, rx, x
  real(kind=mireal), intent(out)  :: prob
  real(kind=mireal) ::xabs

  xabs = abs(x)
  prob = 1.d0/(xabs+lx)/log((lx+rx)/lx)
  if (x < -rx .or. x > rx) prob = 0.0


end subroutine

!Beta function prior as defined by
!https://academic.oup.com/mnrasl/article/434/1/L51/1166128
subroutine beta_prior(a,b,x,prob)
use constants
implicit none

  !In/Out variables
  real(kind=mireal), intent(in) :: a, b, x
  real(kind=mireal), intent(out)  :: prob

  prob = gamma(a+b)/(gamma(a)+gamma(b))
  prob = prob*x**(a-1.)*(1.-x)**(b-1.)

end subroutine

!Prior to exclude multiple solutions for a given period
subroutine t0_prior(a,b,x,prob)
use constants
implicit none

  !In/Out variables
  real(kind=mireal), intent(in) :: a, b, x
  real(kind=mireal), intent(out)  :: prob

  prob = gamma(a+b)/(gamma(a)+gamma(b))
  prob = prob*x**(a-1.)*(1.-x)**(b-1.)

end subroutine

subroutine get_priors(fit_pars,lims,pars_in,priors_out,npars)
use constants
implicit none

  integer, intent(in) :: npars
  real(kind=mireal), intent(in), dimension(0:2*npars-1) :: lims
  real(kind=mireal), intent(in), dimension(0:npars-1) :: pars_in
  real(kind=mireal), intent(out), dimension(0:npars-1) :: priors_out
  character, intent(in), dimension(0:npars-1) :: fit_pars
!Local
  integer :: j

  priors_out = 1.d0
  do j = 0, npars - 1
    if ( fit_pars(j) == 'u' ) then
      call uniform_prior(lims(2*j),lims(2*j+1),pars_in(j),priors_out(j))
    else if ( fit_pars(j) == 'f' ) then
      !nothing
    else if ( fit_pars(j) == 'g' ) then
      call gauss_prior(lims(2*j),lims(2*j+1),pars_in(j),priors_out(j))
    else if ( fit_pars(j) == 'j' ) then
      call jeffreys_prior(lims(2*j),lims(2*j+1),pars_in(j),priors_out(j))
    else if ( fit_pars(j) == 'm' ) then
      call modjeffreys_prior(lims(2*j),lims(2*j+1),pars_in(j),priors_out(j))
    else if ( fit_pars(j) == 'b' ) then
      call beta_prior(lims(2*j),lims(2*j+1),pars_in(j),priors_out(j))
    else
      print*, 'The prior ', fit_pars(j), ' has not been defined in get_priors'
      stop
    end if
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
      if (lims(2*j+1) < lims(2*j)) then
        print *, "One of your uniform U[a,b] priors is wrong!"
        print *, "Such that a > b"
        stop
      end if
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
      pars_out(j) = lims(2*j+1) * r_real
    else if ( fit_pars(j) == 'n' ) then
      call random_number(r_real)
      pars_out(j) = lims(2*j+1) * (r_real - 0.5)
    else if ( fit_pars(j) == 'r' ) then
      call random_number(r_real)
      pars_out(j) = r_real
    else if ( fit_pars(j) == 'o' ) then
      call random_number(r_real)
      pars_out(j) = lims(2*j+1) - lims(2*j)
      pars_out(j) = lims(2*j) + r_real * pars_out(j)
    else if ( fit_pars(j) == 'b' ) then
      call random_number(r_real)
      pars_out(j) = r_real
    else if ( fit_pars(j) == 't' ) then
      call random_number(r_real)
      pars_out(j) = lims(2*j+1) - lims(2*j)
      pars_out(j) = lims(2*j) + r_real * pars_out(j)
    else
      print*, 'The prior ', fit_pars(j), ' has not been defined in create_chains'
      stop
    end if

  end do

end subroutine
