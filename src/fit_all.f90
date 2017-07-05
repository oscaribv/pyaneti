subroutine get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab,tff,flags,&
           t_cad,n_cad,pars,rvs,ldc,trends,jrv,jtr,chi2_rv,chi2_tr,npl,n_tel,size_rv,size_tr)
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npl, n_tel, n_cad !size of RV and LC data
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: tlab
!pars = T0, P, e, w, b, a/R*, Rp/R*, K -> for each planet
  double precision, intent(in) :: pars(0:8*npl-1), rvs(0:n_tel-1), ldc(0:1)
  double precision, intent(in) :: t_cad
  double precision, intent(in) :: trends(0:1)
  double precision, dimension(0:n_tel-1), intent(in) :: jrv
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
  call find_chi2_rv(x_rv,y_rv,e_rv,tlab,pars_rv,jrv,&
                    flag_rv,chi2_rv,size_rv,n_tel,npl)


end subroutine

subroutine multi_all_stretch_move( &
           x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab, &  !Data vars
           stellar_pars,afk,&                     !Stellar parameters and flag|
           flags, total_fit_flag,is_jit, &        !flags
           fit_all, fit_rvs, fit_ldc,fit_trends, &!fitting controls
           nwalks, maxi, thin_factor, nconv, &    !mcmc evolution controls
           lims, lims_rvs, lims_ldc, &            !prior limits
           n_cad, t_cad, &                        !cadence cotrols
           npl, n_tel, &                          !planets and telescopes
           size_rv, size_tr &                     !data sizes
           )
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npl, n_tel !size of RV and LC data
  integer, intent(in) :: nwalks, maxi, thin_factor, nconv, n_cad
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: tlab
  double precision, intent(in), dimension(0:3) :: stellar_pars
  double precision, intent(in), dimension(0:2*8*npl - 1):: lims !, lims_p
  double precision, intent(in), dimension(0:2*n_tel - 1) :: lims_rvs !, lims_p_rvs
  double precision, intent(in), dimension(0:3) :: lims_ldc !, lims_p_ldc
  double precision, intent(in) ::  t_cad
  character, intent(in) :: fit_trends(0:1)
  character, intent(in) :: fit_all(0:8*npl-1), fit_rvs(0:n_tel-1), fit_ldc(0:1)
  logical, intent(in) :: flags(0:5), total_fit_flag(0:1) !CHECK THE SIZE
  logical, intent(in) :: afk(0:npl-1), is_jit(0:1)
!Local variables
  double precision, dimension(0:nwalks-1,0:8*npl-1) :: pars_old, pars_new
  double precision, dimension(0:nwalks-1,0:8*npl-1) :: priors_old, priors_new
  double precision, dimension(0:nwalks-1,0:1) :: priors_ldc_old, priors_ldc_new
  double precision, dimension(0:nwalks-1,0:n_tel-1) :: rvs_old, rvs_new
  double precision, dimension(0:nwalks-1,0:1) :: ldc_old, ldc_new
  double precision, dimension(0:nwalks-1) :: r_rand, z_rand, mstar, rstar
  double precision, dimension(0:nwalks-1) :: chi2_old_total, chi2_new_total, chi2_red
  double precision, dimension(0:nwalks-1) :: chi2_old_rv, chi2_old_tr
  double precision, dimension(0:nwalks-1) :: chi2_new_rv, chi2_new_tr
  double precision, dimension(0:nwalks-1) :: jitter_tr_old, jitter_tr_new
  double precision, dimension(0:nwalks-1,0:n_tel-1) :: jitter_rv_old, jitter_rv_new
  double precision, dimension(0:nwalks-1,0:1) :: tds_old, tds_new !linear and quadratic terms
  double precision, dimension(0:nwalks-1,0:8*npl-1,0:nconv-1) :: pars_chains
  double precision, dimension(0:nwalks-1,0:nconv-1) :: chi2_rv_chains, chi2_tr_chains
  double precision, dimension(0:nwalks-1,0:nconv-1) :: jitter_tr_chains
  double precision, dimension(0:nwalks-1,0:n_tel-1,0:nconv-1) :: jitter_rv_chains
  double precision, dimension(0:nwalks-1,0:1,0:nconv-1) :: tds_chains, ldc_chains
  double precision, dimension(0:nwalks-1,0:n_tel-1,0:nconv-1) :: rvs_chains
  double precision, dimension(0:3) :: lims_trends
  double precision, dimension(0:nwalks-1) :: log_prior_old, log_prior_new
  double precision, dimension(0:nwalks-1) :: log_likelihood_old, log_likelihood_new
  double precision, dimension(0:nwalks-1) :: log_errs_old, log_errs_new
  integer, dimension(0:nwalks/2-1) :: r_int
  double precision  :: a_factor, dof, tds, qq
  double precision  :: lims_e_dynamic(0:1,0:npl-1)
  double precision  :: mstar_mean, mstar_sigma, rstar_mean, rstar_sigma
  double precision  :: a_mean(0:npl-1), a_sigma(0:npl-1)
  double precision :: two_pi = 2.d0*3.1415926535897932384626d0
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

  !Get the stellar parameters
  mstar_mean  = stellar_pars(0)
  mstar_sigma = stellar_pars(1)
  rstar_mean  = stellar_pars(2)
  rstar_sigma = stellar_pars(3)

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
  if ( is_jit(0) ) spar = spar + 1*n_tel
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
    do m = 0, n_tel - 1
      call gauss_random_bm(e_rv(0)*1e-1,e_rv(0)*1e-2,jitter_rv_old(:,m),nwalks)
    end do
  end if
  if ( is_jit(1) ) &
    call gauss_random_bm(e_tr(0)*1e-1,e_tr(0)*1e-2,jitter_tr_old,nwalks)
  log_errs_new = 0.0d0
  log_errs_old = 0.0d0
  do nk = 0, nwalks - 1
    if ( total_fit_flag(0) ) then
      do m = 0, size_rv - 1
        log_errs_old(nk) = log_errs_old(nk) + &
        log( 1.0d0/sqrt( two_pi * ( e_rv(m)**2 + jitter_rv_old(nk,tlab(m))**2 ) ) )
      end do
    end if
    if ( total_fit_flag(1) ) &
    log_errs_old(nk) = log_errs_old(nk) + &
    sum(log( 1.0d0/sqrt( two_pi * ( e_tr(:)**2 + jitter_tr_old(nk)**2 ) ) ) )
  end do

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

  call gauss_random_bm(mstar_mean,mstar_sigma,mstar,nwalks)
  call gauss_random_bm(rstar_mean,rstar_sigma,rstar,nwalks)

  priors_old(:,:) = 1.d0
  priors_new(:,:) = 1.d0
  priors_ldc_old(:,:) = 1.d0
  priors_ldc_new(:,:) = 1.d0

  !Let us create uniformative random priors
  is_limit_good = .false.
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

      !Will we use spectroscopic priors for some planets?
      do m = 0, npl - 1
        if ( afk(m) ) then
          !The parameter comes from 3rd Kepler law !pars_old(1) is the period
          call get_a_scaled(mstar(nk),rstar(nk),pars_old(nk,1+8*m),pars_old(nk,5+8*m),1)
          call get_a_err(mstar_mean,mstar_sigma,rstar_mean,rstar_sigma,&
               pars_old(nk,1+8*m),a_mean(m),a_sigma(m))
          call gauss_prior(a_mean(m),a_sigma(m),pars_old(nk,5+8*m),priors_old(nk,5+8*m))
        end if
      end do

      log_prior_old(nk) = sum( log(priors_old(nk,:) ) + sum( log(priors_ldc_old(nk,:) ) ) )

      call get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab, &
           total_fit_flag,flags,t_cad,n_cad,pars_old(nk,:),rvs_old(nk,:), &
           ldc_old(nk,:),tds_old(nk,:),jitter_rv_old(nk,:),jitter_tr_old(nk),&
           chi2_old_rv(nk),chi2_old_tr(nk),npl,n_tel,size_rv,size_tr)

      chi2_old_total(nk) = chi2_old_rv(nk) + chi2_old_tr(nk)

      log_likelihood_old(nk) = log_prior_old(nk) + log_errs_old(nk) - 0.5d0 * chi2_old_total(nk)

  end do

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
      pars_new(nk,:)    = pars_new(nk,:) + wtf_all(:) * z_rand(nk) *   &
                        ( pars_old(nk,:) - pars_new(nk,:) )
      rvs_new(nk,:)     = rvs_new(nk,:) + wtf_rvs(:) * z_rand(nk) *    &
                        ( rvs_old(nk,:) - rvs_new(nk,:) )
      ldc_new(nk,:)     = ldc_new(nk,:) + wtf_ldc(:) * z_rand(nk) *    &
                        ( ldc_old(nk,:) - ldc_new(nk,:) )
      tds_new(nk,:)     = tds_new(nk,:) + wtf_trends(:) * z_rand(nk) * &
                        ( tds_old(nk,:) - tds_new(nk,:) )
      jitter_rv_new(nk,:) = jitter_rv_new(nk,:) + z_rand(nk) *             &
                         ( jitter_rv_old(nk,:) - jitter_rv_new(nk,:) )
      jitter_tr_new(nk) = jitter_tr_new(nk) + z_rand(nk) *             &
                         ( jitter_tr_old(nk) - jitter_tr_new(nk) )

      call get_priors(fit_all,lims,pars_new(nk,:),priors_new(nk,:),8*npl)
      call get_priors(fit_ldc,lims_ldc,ldc_new(nk,:),priors_ldc_new(nk,:),2)

      do m = 0, npl - 1
        if ( afk(m) ) then
          !The parameter comes from 3rd Kepler law
          call get_a_err(mstar_mean,mstar_sigma,rstar_mean,rstar_sigma,&
               pars_new(nk,1+8*m),a_mean(m),a_sigma(m))
          call gauss_prior(a_mean(m),a_sigma(m),pars_new(nk,5+8*m),priors_new(nk,5+8*m))
        end if
      end do

      limit_prior = PRODUCT(priors_new(nk,:)) * PRODUCT(priors_ldc_new(nk,:) )

      !Let us check if the new parameters are inside the limits
      is_limit_good = .true.
      if ( limit_prior < 1.d-20 ) is_limit_good = .false.
      if ( is_limit_good ) then
        if ( ANY( jitter_rv_new(nk,:) < 0.0d0 ) .or. jitter_tr_new(nk) < 0.0d0 ) is_limit_good = .false.
      end if

      chi2_new_total(nk) = huge(0.0d0) !A really big number!

      if ( is_limit_good ) then !If we are inside the limits, let us calculate chi^2

        call get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab,       &
             total_fit_flag,flags, t_cad,n_cad,pars_new(nk,:),rvs_new(nk,:),  &
             ldc_new(nk,:),tds_new(nk,:),jitter_rv_new(nk,:),jitter_tr_new(nk), &
             chi2_new_rv(nk),chi2_new_tr(nk),npl,n_tel,size_rv,size_tr)

        chi2_new_total(nk) = chi2_new_rv(nk) + chi2_new_tr(nk)

      end if

      !ADD JITTER TERMS
      log_errs_new(nk) = 0.0d0
      if ( total_fit_flag(0) ) then
        do m = 0, size_rv - 1
          log_errs_new(nk) = log_errs_new(nk) + &
          log( 1.0d0/sqrt( two_pi * ( e_rv(m)**2 + jitter_rv_new(nk,tlab(m))**2 ) ) )
        end do
      end if
      if ( total_fit_flag(1) ) &
         log_errs_new(nk) = log_errs_new(nk) + &
         sum(log( 1.0d0/sqrt( two_pi * ( e_tr(:)**2 + jitter_tr_new(nk)**2 ) ) ) )

      log_prior_new(nk) = sum( log(priors_new(nk,:) ) ) + &
                          sum( log(priors_ldc_new(nk,:) ) )

      log_likelihood_new(nk) = log_prior_new(nk) + log_errs_new(nk) - 0.5d0 * chi2_new_total(nk)

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
      write(101,*) n, nk, chi2_rv_chains(nk,n),chi2_tr_chains(nk,n), pars_chains(nk,:,n), &
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
