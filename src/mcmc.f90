subroutine mcmc_stretch_move( &
           x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab,jrvlab, &  !Data vars
           flags, total_fit_flag,is_jit, &              !flags
           fit_all, fit_rvs, fit_ldc,fit_trends, &      !fitting controls
           nwalks, maxi, thin_factor, nconv, &          !mcmc evolution controls
           lims, lims_rvs, lims_ldc, &                  !prior limits
           n_cad, t_cad, &                              !cadence cotrols
           npl, n_tel, n_jrv, &                          !planets and telescopes
           size_rv, size_tr &                           !data sizes
           )
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npl, n_tel, n_jrv !size of RV and LC data
  integer, intent(in) :: nwalks, maxi, thin_factor, nconv, n_cad
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: tlab, jrvlab
  double precision, intent(in), dimension(0:2*8*npl - 1):: lims !, lims_p
  double precision, intent(in), dimension(0:2*n_tel - 1) :: lims_rvs !, lims_p_rvs
  double precision, intent(in), dimension(0:3) :: lims_ldc !, lims_p_ldc
  double precision, intent(in) ::  t_cad
  character, intent(in) :: fit_trends(0:1)
  character, intent(in) :: fit_all(0:8*npl-1), fit_rvs(0:n_tel-1), fit_ldc(0:1)
  logical, intent(in) :: flags(0:5), total_fit_flag(0:1) !CHECK THE SIZE
  logical, intent(in) :: is_jit(0:1)
!Local variables
  double precision, dimension(0:nwalks-1,0:8*npl-1) :: pars_old, pars_new
  double precision, dimension(0:nwalks-1,0:8*npl-1) :: priors_old, priors_new
  double precision, dimension(0:nwalks-1,0:1) :: priors_ldc_old, priors_ldc_new
  double precision, dimension(0:nwalks-1,0:n_tel-1) :: rvs_old, rvs_new
  double precision, dimension(0:nwalks-1,0:1) :: ldc_old, ldc_new
  double precision, dimension(0:nwalks-1) :: r_rand, z_rand
  double precision, dimension(0:nwalks-1) :: chi2_old_total, chi2_new_total, chi2_red
  double precision, dimension(0:nwalks-1) :: chi2_old_rv, chi2_old_tr
  double precision, dimension(0:nwalks-1) :: chi2_new_rv, chi2_new_tr
  double precision, dimension(0:nwalks-1) :: jitter_tr_old, jitter_tr_new
  double precision, dimension(0:nwalks-1,0:n_jrv-1) :: jitter_rv_old, jitter_rv_new
  double precision, dimension(0:nwalks-1,0:1) :: tds_old, tds_new !linear and quadratic terms
  double precision, dimension(0:nwalks-1,0:8*npl-1,0:nconv-1) :: pars_chains
  double precision, dimension(0:nwalks-1,0:nconv-1) :: chi2_rv_chains, chi2_tr_chains, loglike_chains
  double precision, dimension(0:nwalks-1,0:nconv-1) :: jitter_tr_chains
  double precision, dimension(0:nwalks-1,0:n_tel-1,0:nconv-1) :: jitter_rv_chains
  double precision, dimension(0:nwalks-1,0:1,0:nconv-1) :: tds_chains, ldc_chains
  double precision, dimension(0:nwalks-1,0:n_tel-1,0:nconv-1) :: rvs_chains
  double precision, dimension(0:3) :: lims_trends
  double precision, dimension(0:nwalks-1) :: log_prior_old, log_prior_new
  double precision, dimension(0:nwalks-1) :: log_likelihood_old, log_likelihood_new
  integer, dimension(0:nwalks/2-1) :: r_int
  double precision  :: a_factor, dof, tds, qq
  double precision  :: lims_e_dynamic(0:1,0:npl-1)
  double precision  :: a_mean(0:npl-1), a_sigma(0:npl-1)
  double precision :: limit_prior
  logical :: continua, is_limit_good, is_cvg
  integer :: nk, j, n, m, o, n_burn, spar, spar1, iensemble, inverted(0:1)
  integer :: nks, nke
  integer :: wtf_trends(0:1)
  integer :: wtf_all(0:8*npl-1), wtf_rvs(0:n_tel-1), wtf_ldc(0:1)
!external calls
  external :: init_random_seed, find_chi2_tr, find_chi2_rv

  !call the random seed
  print *, 'CREATING RANDOM SEED'
  call init_random_seed()

  !Let us fill the what-to-fit vectors
  wtf_all = 0
  do o = 0, 8*npl-1
    if ( fit_all(o) .ne. 'f' ) wtf_all(o) = 1
  end do
  wtf_ldc = 0
  do o = 0, 1
    if ( fit_ldc(o) .ne. 'f' ) wtf_ldc(o) = 1
  end do
  wtf_rvs = 0
  do o = 0, n_tel-1
    if ( fit_rvs(o) .ne. 'f' ) wtf_rvs(o) = 1
  end do
  wtf_trends = 0
  do o = 0, 1
    if ( fit_trends(o) .ne. 'f' ) wtf_trends(o) = 1
  end do

  spar = sum(wtf_all) + sum(wtf_ldc) + sum(wtf_rvs) + sum(wtf_trends)
  !spar -> size of parameters, dof -> degrees of freedom
  if ( is_jit(0) ) spar = spar + 1*n_jrv
  if ( is_jit(1) ) spar = spar + 1

  spar1 = spar - 1

  dof = 0
  if ( total_fit_flag(0) ) dof = dof + size_rv
  if ( total_fit_flag(1) ) dof = dof + size_tr
  dof  = dof - spar

  !Jitter vars
  jitter_rv_old(:,:) = 0.0d0
  jitter_tr_old(:) = 0.0d0
  jitter_rv_new(:,:) = 0.0d0
  jitter_tr_new(:) = 0.0d0

  if ( is_jit(0) ) then
    do m = 0, n_jrv - 1
      do nk = 0, nwalks - 1
        call create_chains('u',(/0.d0,e_rv(0)/),jitter_rv_old(nk,m),1)
      end do
    end do
  end if
  if ( is_jit(1) ) then
    do nk = 0, nwalks - 1
      call create_chains('u',(/0.d0,e_tr(0)/),jitter_tr_old(nk),1)
    end do
  end if

  !Linear and quadratic terms
  lims_trends(:) = 0.0d0
  tds = 0.0d0
  if ( wtf_trends(0) == 1 ) then !linear trend
    lims_trends(0) = -1.0d-1
    lims_trends(1) =  1.0d-1
  end if
  if ( wtf_trends(1) == 1) then !quadratic trend
    lims_trends(2) = -1.0d-1
    lims_trends(3) =  1.0d-1
  end if


  print *, 'CREATING CHAINS'

  priors_old(:,:) = 1.d0
  priors_new(:,:) = 1.d0
  priors_ldc_old(:,:) = 1.d0
  priors_ldc_new(:,:) = 1.d0

  !Let us create uniformative random priors
  is_limit_good = .false.
  !$OMP PARALLEL &
  !$OMP PRIVATE(is_limit_good,m,limit_prior,a_mean,a_sigma)
  !$OMP DO SCHEDULE(DYNAMIC)
  do nk = 0, nwalks - 1

      call create_chains(fit_all,lims,pars_old(nk,:),8*npl)

      !If we are using e and w parameterization, let us be sure we do not have e > 1
      if ( flags(1) ) then
        do m = 0, npl-1
          lims_e_dynamic(:,m) = sqrt( 1.d0 - pars_old(nk,2+8*m)**2 )
          lims_e_dynamic(0,m) = - lims_e_dynamic(0,m)
          if ( fit_all(3+8*m) == 'f' ) lims_e_dynamic(0,m) = lims(3*2+8*m*2)
          call create_chains(fit_all(3+8*m),lims_e_dynamic(:,m),pars_old(nk,3+8*m),1)
        end do
      end if

      call get_priors(fit_all,lims,pars_old(nk,:),priors_old(nk,:),8*npl)

      call create_chains(fit_trends,lims_trends,tds_old(nk,:),2)
      !call get_priors(fit_trends,lims_trends,tds_old(nk,:),priors_trends(nk,:),2)

      call create_chains(fit_rvs,lims_rvs,rvs_old(nk,:),n_tel)
      !call get_priors(fit_rvs,lims_rvs,rvs_old(nk,:),priors_rvs(nk,:),n_tel)

      call create_chains(fit_ldc,lims_ldc,ldc_old(nk,:),2)
      call get_priors(fit_ldc,lims_ldc,ldc_old(nk,:),priors_ldc_old(nk,:),2)

      log_prior_old(nk) = sum( log(priors_old(nk,:) ) + sum( log(priors_ldc_old(nk,:) ) ) )

      call get_loglike(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab,jrvlab, &
           total_fit_flag,flags,t_cad,n_cad,pars_old(nk,:),rvs_old(nk,:), &
           ldc_old(nk,:),tds_old(nk,:),jitter_rv_old(nk,:),jitter_tr_old(nk),&
           log_likelihood_old(nk),chi2_old_rv(nk),chi2_old_tr(nk),npl,n_tel,n_jrv,size_rv,size_tr)

      chi2_old_total(nk) = chi2_old_rv(nk) + chi2_old_tr(nk)
      log_likelihood_old(nk) = log_prior_old(nk) + log_likelihood_old(nk)

  end do
  !$OMP END PARALLEL

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
  a_factor = 2.d0
  n_burn = 1
  inverted = (/ 1 , 0 /)

  !The infinite cycle starts!
  print *, 'STARTING INFINITE LOOP!'
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
        pars_new(nk,:)      = pars_old(r_int(nk-nks),:)
        rvs_new(nk,:)       = rvs_old(r_int(nk-nks),:)
        ldc_new(nk,:)       = ldc_old(r_int(nk-nks),:)
        tds_new(nk,:)       = tds_old(r_int(nk-nks),:)
        jitter_rv_new(nk,:) = jitter_rv_old(r_int(nk-nks),:)
        jitter_tr_new(nk)   = jitter_tr_old(r_int(nk-nks))
      end do

    !Paralellization calls
    !$OMP PARALLEL &
    !$OMP PRIVATE(is_limit_good,qq,m,limit_prior,a_mean,a_sigma)
    !$OMP DO SCHEDULE(DYNAMIC)
    do nk = nks, nke

      !Generate the random step to perform the stretch move
      call find_gz(a_factor,z_rand(nk))

      !Perform the stretch move
      !Eq. (7), Goodman & Weare (2010)
      pars_new(nk,:)    = pars_new(nk,:) +  z_rand(nk) *   &
                        ( pars_old(nk,:) - pars_new(nk,:) )
      rvs_new(nk,:)     = rvs_new(nk,:) +  z_rand(nk) *    &
                        ( rvs_old(nk,:) - rvs_new(nk,:) )
      ldc_new(nk,:)     = ldc_new(nk,:) +  z_rand(nk) *    &
                        ( ldc_old(nk,:) - ldc_new(nk,:) )
      tds_new(nk,:)     = tds_new(nk,:) +  z_rand(nk) * &
                        ( tds_old(nk,:) - tds_new(nk,:) )
      jitter_rv_new(nk,:) = jitter_rv_new(nk,:) + z_rand(nk) *             &
                         ( jitter_rv_old(nk,:) - jitter_rv_new(nk,:) )
      jitter_tr_new(nk) = jitter_tr_new(nk) + z_rand(nk) *             &
                         ( jitter_tr_old(nk) - jitter_tr_new(nk) )

      call get_priors(fit_all,lims,pars_new(nk,:),priors_new(nk,:),8*npl)
      call get_priors(fit_ldc,lims_ldc,ldc_new(nk,:),priors_ldc_new(nk,:),2)


      limit_prior = PRODUCT(priors_new(nk,:)) * PRODUCT(priors_ldc_new(nk,:) )

      !Let us check if the new parameters are inside the limits
      is_limit_good = .true.
      if ( limit_prior < 1.d-100 ) is_limit_good = .false.
      if ( is_limit_good ) then
        if ( ANY( jitter_rv_new(nk,:) < 0.0d0 ) .or. jitter_tr_new(nk) < 0.0d0 ) is_limit_good = .false.
      end if

      chi2_new_total(nk) = huge(0.0d0) !A really big number!
      log_likelihood_new(nk) = -huge(0.e0)

      if ( is_limit_good ) then !If we are inside the limits, let us calculate chi^2

      call get_loglike(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab,jrvlab, &
           total_fit_flag,flags,t_cad,n_cad,pars_new(nk,:),rvs_new(nk,:), &
           ldc_new(nk,:),tds_new(nk,:),jitter_rv_new(nk,:),jitter_tr_new(nk),&
           log_likelihood_new(nk),chi2_new_rv(nk),chi2_new_tr(nk),npl,n_tel,&
           n_jrv,size_rv,size_tr)

      chi2_new_total(nk) = chi2_new_rv(nk) + chi2_new_tr(nk)

      end if

      log_prior_new(nk) = sum( log(priors_new(nk,:) ) ) + &
                          sum( log(priors_ldc_new(nk,:) ) )

      log_likelihood_new(nk) = log_prior_new(nk) + log_likelihood_new(nk)

      qq = log_likelihood_new(nk) - log_likelihood_old(nk)
      !z^(pars-1) normalization factor needed to perform the stretch move
      !Goodman & Weare (2010)
      qq = z_rand(nk)**spar1 * exp(qq)

      !Check if the new likelihood is better
      if ( qq > r_rand(nk) ) then
        !If yes, let us save it as the old vectors
        log_likelihood_old(nk) = log_likelihood_new(nk)
        chi2_old_total(nk)     = chi2_new_total(nk)
        chi2_old_rv(nk)        = chi2_new_rv(nk)
        chi2_old_tr(nk)        = chi2_new_tr(nk)
        pars_old(nk,:)         = pars_new(nk,:)
        rvs_old(nk,:)          = rvs_new(nk,:)
        ldc_old(nk,:)          = ldc_new(nk,:)
        tds_old(nk,:)          = tds_new(nk,:)
        jitter_rv_old(nk,:)    = jitter_rv_new(nk,:)
        jitter_tr_old(nk)      = jitter_tr_new(nk)
        priors_old(nk,:)       = priors_new(nk,:)
        priors_ldc_old(nk,:)   = priors_ldc_new(nk,:)
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
      ldc_chains(:,:,n)       = ldc_old(:,:)
      rvs_chains(:,:,n)       = rvs_old(:,:)
      tds_chains(:,:,n)       = tds_old(:,:)
      jitter_rv_chains(:,:,n) = jitter_rv_old(:,:)
      jitter_tr_chains(:,n)   = jitter_tr_old(:)
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
        do o = 0, 8*npl-1
          if (wtf_all(o) == 1 ) then !perform G-R test
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
  open(unit=201,file='jitter_data.dat',status='unknown')
  open(unit=301,file='trends_data.dat',status='unknown')

  do n = 0, nconv - 1
    do nk = 0, nwalks - 1
      write(101,*) n, loglike_chains(nk,n), chi2_rv_chains(nk,n),chi2_tr_chains(nk,n), pars_chains(nk,:,n), &
                   ldc_chains(nk,:,n), rvs_chains(nk,:,n)
      if ( is_jit(0) .or.  is_jit(1) ) &
          write(201,*) jitter_rv_chains(nk,:,n), jitter_tr_chains(nk,n)
      if ( sum(wtf_trends)  > 0 ) &
         write(301,*) tds_chains(nk,:,n)
    end do
  end do

  !Close file
  close(101)
  close(201)
  close(301)

end subroutine

!#--------------------------------------------------------------------------------

subroutine mcmc_stretch_move_new(            &
           x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,    &
           rvlab,jrvlab,trlab,jtrlab,        &
           flags, total_fit_flag,            &
           prior_flags, prior_vals,          &
           npars, model_int,                 &
           model_double,                     &
           nwalks, maxi, thin_factor, nconv, &
           size_rv, size_tr)
implicit none

!npars = 7*npl + (npl + LDC)*nbands + noffsets + njitter + ntrends

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npars, model_int(0:6)
  !mcmc_int parameters
  integer :: nwalks, maxi, thin_factor, nconv
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: rvlab, jrvlab
  integer, intent(in), dimension(0:size_tr-1) :: trlab, jtrlab
  double precision, intent(in), dimension(0:2*npars - 1):: prior_vals
  double precision, intent(in) ::  model_double(0)
  character, intent(in) :: prior_flags(0:npars-1)
  logical, intent(in) :: flags(0:5), total_fit_flag(0:1)
!Local variables
  double precision, dimension(0:nwalks-1,0:npars-1) :: pars_old, pars_new
  double precision, dimension(0:nwalks-1,0:npars-1) :: priors_old, priors_new
  double precision, dimension(0:nwalks-1) :: r_rand, z_rand
  double precision, dimension(0:nwalks-1) :: chi2_old_total, chi2_new_total, chi2_red
  double precision, dimension(0:nwalks-1) :: chi2_old_rv, chi2_old_tr
  double precision, dimension(0:nwalks-1) :: chi2_new_rv, chi2_new_tr
  double precision, dimension(0:nwalks-1,0:npars-1,0:nconv-1) :: pars_chains
  double precision, dimension(0:nwalks-1,0:nconv-1) :: chi2_rv_chains, chi2_tr_chains, loglike_chains
  double precision, dimension(0:nwalks-1) :: log_prior_old, log_prior_new
  double precision, dimension(0:nwalks-1) :: log_likelihood_old, log_likelihood_new
  integer, dimension(0:nwalks/2-1) :: r_int
  double precision  :: a_factor, dof, tds, qq
  double precision :: limit_prior
  logical :: continua, is_limit_good, is_cvg
  integer :: nk, j, n, m, o, n_burn, spar, spar1, iensemble, inverted(0:1)
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
  !$OMP PARALLEL &
  !$OMP PRIVATE(is_limit_good,m,limit_prior)
  !$OMP DO SCHEDULE(DYNAMIC)
  do nk = 0, nwalks - 1

      call create_chains(prior_flags,prior_vals,pars_old(nk,:),npars)

      call get_priors(prior_flags,prior_vals,pars_old(nk,:),priors_old(nk,:),npars)

      log_prior_old(nk) = sum( log(priors_old(nk,:) ) )

      call get_loglike_new(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           rvlab,jrvlab,trlab,jtrlab,total_fit_flag,flags,&
           pars_old(nk,:),model_int,model_double,&
           npars,log_likelihood_old(nk),chi2_old_rv(nk),chi2_old_tr(nk),size_rv,size_tr)

      chi2_old_total(nk) = chi2_old_rv(nk) + chi2_old_tr(nk)
      log_likelihood_old(nk) = log_prior_old(nk) + log_likelihood_old(nk)

  end do
  !$OMP END PARALLEL

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
  a_factor = 2.d0
  n_burn = 1
  inverted = (/ 1 , 0 /)

  !The infinite cycle starts!
  print *, 'STARTING INFINITE LOOP!'
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
    !$OMP PRIVATE(is_limit_good,qq,m,limit_prior)
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

      call get_loglike_new(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,    &
           rvlab,jrvlab,trlab,jtrlab,total_fit_flag,flags,   &
           pars_new(nk,:),model_int,model_double,            &
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


subroutine get_loglike_new(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
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
  integer :: m
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
   call get_total_chi2_new(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
                           rvlab,jrvlab,trlab,jtrlab,tff,flags,&
                           t_cad,n_cad,pars,chi2_rv,chi2_tr,log_errs,npars,&
                           npl,ntels,nbands,nldc,njrv,njtr,size_rv,size_tr)

    chi2_total = chi2_rv + chi2_tr

    loglike = log_errs - 0.5d0 * chi2_total

end subroutine


subroutine get_total_chi2_new(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr, &
           rvlab,jrvlab,trlab,jtrlab,tff,flags,&
           t_cad,n_cad,pars, &
           chi2_rv,chi2_tr,log_errs,npars,npl,ntels,nbands,nldc,njrv,njtr,size_rv,size_tr)
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
  double precision :: ppars(0:8*npl-1), rvs(0:ntels-1)
  double precision :: trends(0:1)
  double precision :: ldc(0:nldc*nbands-1)
  double precision, dimension(0:njrv-1) :: jrv
  double precision, dimension(0:njtr-1) :: jtr
  double precision :: pars_rv(0:7+ntels-1,0:npl-1)
  double precision :: pars_tr(0:5,0:npl-1), rps(0:nbands-1,0:npl-1)
  double precision :: two_pi = 2.d0*3.1415926535897932384626d0
  logical :: flag_rv(0:3), flag_tr(0:3)
  integer :: i, j, m
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
    rps(:,i)         = pars(srp+nbands*i:srp+nbands*(i+1)-1)
    pars_rv(0:3,i)   = pars(j:j+3)
    pars_rv(4,i)     = pars(j+6)
    pars_rv(5:6,i)   = pars(strends:strends+1)
    pars_rv(7:7+ntels-1,i) = pars(srv:srv+ntels-1)
!  print *, srp+nbands*i,srp+nbands*(i+1)-1
!  print *, pars(srp+nbands*i:srp+nbands*(i+1)-1)
  end do
!  print *, pars_tr
!  print *, pars_rv
!  print *, rps
!  stop

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

!  print *, trlab
!  print *, jtrlab
!  print *, pars_tr
!  print *, rps
!  print *, ldc
!  print *, jtr

  if (tff(1) ) &
  call find_chi2_tr_nuevo(x_tr,y_tr,e_tr,trlab,jtrlab,pars_tr,rps,ldc,jtr,flag_tr, &
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