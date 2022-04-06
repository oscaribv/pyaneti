!--------------------------------------------------------------------------------
subroutine mcmc_stretch_move(                &
           x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,    &
           rvlab,jrvlab,trlab,jtrlab,        &
           flags, total_fit_flag,            &
           prior_flags, prior_vals,          &
           kernels,                          &
           model_int,                        &
           model_double,npars,nmodel_int,nmodel_double, &
           nwalks, maxi, thin_factor, nconv, &
           size_rv, size_tr)
use constants
implicit none

!npars = 7*npl + (npl + LDC)*nbands + noffsets + njitter + ntrends + GP_hyper_parameters

!In/Out variables
  integer, intent(in) :: size_rv, size_tr
  integer, intent(in) :: npars, nmodel_int, nmodel_double
  integer, intent(in) :: model_int(0:nmodel_int-1)
  !mcmc_int parameters
  integer :: nwalks, maxi, thin_factor, nconv
  real(kind=mireal), intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  real(kind=mireal), intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: rvlab, jrvlab
  integer, intent(in), dimension(0:size_tr-1) :: trlab, jtrlab
  real(kind=mireal), intent(in), dimension(0:2*npars - 1):: prior_vals
  real(kind=mireal), intent(in) ::  model_double(0:nmodel_double-1)
  character, intent(in) :: prior_flags(0:npars-1)
  character(len=6), intent(in) :: kernels
  logical, intent(in) :: flags(0:5), total_fit_flag(0:1)
!Local variables
  real(kind=mireal), dimension(0:nwalks-1,0:npars-1) :: pars_old, pars_new
  real(kind=mireal), dimension(0:nwalks-1,0:npars-1) :: priors_old, priors_new
  real(kind=mireal), dimension(0:nwalks-1) :: r_rand, z_rand
  real(kind=mireal), dimension(0:nwalks-1) :: chi2_old_total, chi2_new_total, chi2_red
  real(kind=mireal), dimension(0:nwalks-1) :: chi2_old_rv, chi2_old_tr
  real(kind=mireal), dimension(0:nwalks-1) :: chi2_new_rv, chi2_new_tr
  real(kind=mireal), dimension(0:nwalks-1,0:npars-1,0:nconv-1) :: pars_chains
  real(kind=mireal), dimension(0:nwalks-1,0:nconv-1) :: chi2_rv_chains, chi2_tr_chains, loglike_chains
  real(kind=mireal), dimension(0:nwalks-1) :: log_prior_old, log_prior_new
  real(kind=mireal), dimension(0:nwalks-1) :: log_likelihood_old, log_likelihood_new
  integer, dimension(0:nwalks/2-1) :: r_int
!  real(kind=mireal)  :: a_factor, dof, qq
  real(kind=mireal)  :: dof, qq
  real(kind=mireal) :: limit_prior
  logical :: continua, is_limit_good, is_cvg
  integer :: nk, j, n, o, n_burn, spar, spar1, iensemble, inverted(0:1)
  integer :: nks, nke
!external calls
  external :: init_random_seed


  !call the random seed
  print *, 'CREATING RANDOM SEED'
  call init_random_seed()

  !spar -> size of parameters, dof -> degrees of freedom
  spar = 0
  do o = 0, SIZE(prior_flags) - 1
    if (prior_flags(o) .ne. 'f') spar = spar + 1
  end do

  spar1 = spar - 1

  dof = 0
  if ( total_fit_flag(0) ) dof = dof + size_rv
  if ( total_fit_flag(1) ) dof = dof + size_tr
  dof  = dof - spar

  print *, 'CREATING CHAINS'

  priors_old(:,:) = 1.d0
  priors_new(:,:) = 1.d0

  is_limit_good = .false.
  !!$OMP PARALLEL &
  !!$OMP PRIVATE(is_limit_good,limit_prior)
  !!$OMP DO SCHEDULE(DYNAMIC)
  !do nk = 0, nwalks - 1
  nk = 0
  do while (nk < nwalks)

      call create_chains(prior_flags,prior_vals,pars_old(nk,:),npars)

      call get_priors(prior_flags,prior_vals,pars_old(nk,:),priors_old(nk,:),npars)

      log_prior_old(nk) = sum( log(priors_old(nk,:) ) )

      call get_loglike(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           rvlab,jrvlab,trlab,jtrlab,total_fit_flag,flags,kernels,&
           pars_old(nk,:),model_int,model_double,nmodel_int,nmodel_double,&
           npars,log_likelihood_old(nk),chi2_old_rv(nk),chi2_old_tr(nk),size_rv,size_tr)

      chi2_old_total(nk) = chi2_old_rv(nk) + chi2_old_tr(nk)
      log_likelihood_old(nk) = log_prior_old(nk) + log_likelihood_old(nk)

      if ( .not. ISNAN(chi2_old_total(nk))) nk = nk + 1

  end do
  !!$OMP END PARALLEL

  chi2_red(:) = chi2_old_total(:) / dof

  !Print the initial cofiguration
  print *, ''
  print *, 'STARTING MCMC CALCULATION'
  if ( total_fit_flag(0) ) &
  print *, 'RV datapoints  = ', size_rv
  if ( total_fit_flag(1) ) &
  print *, 'TR datapoints  = ', size_tr
  print *, 'No. parameters = ', int(spar)
  print *, 'dof            = ', int(dof)
  call print_chain_data(chi2_red,nwalks)

  !Initialize the values
  j = 1
  n = 0
  continua = .true.
  !a_factor = 2.d0
  n_burn = 1
  inverted = (/ 1 , 0 /)

  !The infinite cycle starts!
  print *, 'STARTING MCMC SAMPLING!'
  do while ( continua )

    !Creating the random variables
    !Do the work for all the walkers
    call random_number(r_rand)

    !Perform the paralelization following Foreman-Mackey, 2013
    do iensemble = 0, 1

      nks = iensemble * nwalks/2
      nke = (iensemble + 1) * nwalks/2 - 1

      !Create random integers to be the index of the walks
      !Note that r_ink(i) != i (avoid copy the same walker)
      call random_int(r_int,inverted(iensemble)*nwalks/2,(inverted(iensemble)+1)*nwalks/2 - 1)


      !Pick a random walker from the complementary ensemble
      do nk = nks, nke
        pars_new(nk,:) = pars_old(r_int(nk-nks),:)
      end do

    !Paralellization calls
    !$OMP PARALLEL &
    !$OMP PRIVATE(is_limit_good,qq,limit_prior)
    !$OMP DO SCHEDULE(DYNAMIC)
    do nk = nks, nke

      !Generate the random step to perform the stretch move
      call find_gz(a_factor,z_rand(nk))

      !Perform the stretch move
      !Eq. (7), Goodman & Weare (2010)
      pars_new(nk,:)    = pars_new(nk,:) +  z_rand(nk) * ( pars_old(nk,:) - pars_new(nk,:) )

      call get_priors(prior_flags,prior_vals,pars_new(nk,:),priors_new(nk,:),npars)

      limit_prior = PRODUCT(priors_new(nk,:))

      !Let us check if the new parameters are inside the limits
      is_limit_good = .true.
      if ( limit_prior < 1.d-100 ) is_limit_good = .false.

      chi2_new_total(nk) = huge(0.0d0) !A really big number!
      log_likelihood_new(nk) = -huge(0.e0)

      if ( is_limit_good ) then !If we are inside the limits, let us calculate chi^2

      call get_loglike(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,    &
           rvlab,jrvlab,trlab,jtrlab,total_fit_flag,flags,kernels,   &
           pars_new(nk,:),model_int,model_double,nmodel_int,nmodel_double,            &
           npars,log_likelihood_new(nk),chi2_new_rv(nk),chi2_new_tr(nk),size_rv,size_tr)

      chi2_new_total(nk) = chi2_new_rv(nk) + chi2_new_tr(nk)

      end if

      log_prior_new(nk) = sum( log(priors_new(nk,:) ) )

      log_likelihood_new(nk) = log_prior_new(nk) + log_likelihood_new(nk)

      qq = log_likelihood_new(nk) - log_likelihood_old(nk)
      !z^(pars-1) normalization factor needed to perform the stretch move
      !Goodman & Weare (2010)
      qq = z_rand(nk)**spar1 * exp(qq)

      !Check the acceptance of the new model
      if ( qq > r_rand(nk) ) then
        !If yes, let us save it as the old vectors
        log_likelihood_old(nk) = log_likelihood_new(nk)
        chi2_old_total(nk)     = chi2_new_total(nk)
        chi2_old_rv(nk)        = chi2_new_rv(nk)
        chi2_old_tr(nk)        = chi2_new_tr(nk)
        pars_old(nk,:)         = pars_new(nk,:)
        priors_old(nk,:)       = priors_new(nk,:)
      end if

    end do !walkers
    !$OMP END PARALLEL

    end do !iensemble

    !Compute the reduced chi square
    chi2_red(:) = chi2_old_total(:) / dof

    if ( mod(j,thin_factor) == 0 ) then
      !If the chains have not converged, let us check convergence
      !Let us save a 3D array with the informations of the parameters,
      !the nk and the iteration. This array is used to perform GR test
      pars_chains(:,:,n)      = pars_old(:,:)
      chi2_rv_chains(:,n)     = chi2_old_rv(:)
      chi2_tr_chains(:,n)     = chi2_old_tr(:)
      loglike_chains(:,n)     = log_likelihood_old(:)
      n = n + 1
      !Is it time to check covergence=
      if ( n == nconv ) then
        !Perform G-R test
        n = 0 !reinitilize n
        call print_chain_data(chi2_red,nwalks)
        print *, '=================================='
        print *, '     PERFOMING GELMAN-RUBIN'
        print *, '      TEST FOR CONVERGENCE'
        print *, '=================================='
        !Check convergence for all the parameters
        is_cvg = .true.
        do o = 0, npars-1
          if ( prior_flags(o) .ne. 'f' ) then !perform G-R test
            call gr_test(pars_chains(:,o,:),nwalks,nconv,is_cvg)
          end if
          if ( .not. is_cvg ) exit
        end do

        if ( j < thin_factor*nconv + 1 ) is_cvg = .False.

        if ( .not. is_cvg  ) then
          print *, '=================================='
          print *, 'CHAINS HAVE NOT CONVERGED YET!'
          print *,  nconv*thin_factor,' ITERATIONS MORE!'
          print *, '=================================='
        else
          print *, '=================================='
          print *, '      CHAINS HAVE CONVERGED'
          print *, '=================================='
          print *, '   CREATING OUTPUT DATA FILES'
          print *, '=================================='
          !Let us start the otput file
          continua = .false.
        end if ! is_cvg
      end if !nconv
    end if !j/thin_factor

    !check if we exceed the maximum number of iterations
    if ( j > maxi ) then
      print *, 'Maximum number of iteration reached!'
      continua = .FALSE.
    end if

  j = j + 1

  end do !infinite loop
  !the MCMC part has ended

  !Let us create the output file
  open(unit=101,file='all_data.dat',status='unknown')

  do n = 0, nconv - 1
    do nk = 0, nwalks - 1
      write(101,*) n, loglike_chains(nk,n), chi2_rv_chains(nk,n),chi2_tr_chains(nk,n), pars_chains(nk,:,n)
    end do
  end do

  !Close file
  close(101)

end subroutine