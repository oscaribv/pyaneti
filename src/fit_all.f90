subroutine get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab,plab_tr,tff,flags,&
           t_cad,n_cad,pars,rvs,ldc,trends,jrv,jtr,chi2_rv,chi2_tr,npl,n_tel,size_rv,size_tr)
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npl, n_tel, n_cad !size of RV and LC data
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_tr-1) :: plab_tr
  integer, intent(in), dimension(0:size_rv-1) :: tlab
  !pars = T0, P, e, w, b, a/R*, Rp/R*, K -> for each planet
  double precision, intent(in) :: pars(0:8*npl-1), rvs(0:n_tel-1), ldc(0:1)
  double precision, intent(in) :: t_cad
  double precision, intent(in) :: trends(0:1)
  double precision, intent(in) :: jrv, jtr
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
  call find_chi2_tr(x_tr,y_tr,e_tr,plab_tr,pars_tr,jtr,flag_tr,&
                        ldc,n_cad,t_cad,chi2_tr,size_tr,npl)
  if (tff(0) ) &
  call find_chi2_rv(x_rv,y_rv,e_rv,tlab,pars_rv,jrv,&
                    flag_rv,chi2_rv,size_rv,n_tel,npl)


end subroutine

subroutine multi_all_stretch_move( &
           x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab, &  !Data vars
           plab_tr,pars, rvs, ldc, &              !parameters vars
           stellar_pars,afk,&
           flags, total_fit_flag,is_jit, &               !flags
           wtf_all, wtf_rvs, wtf_ldc,wtf_trends, &           !fitting controls
           nwalks, maxi, thin_factor, nconv, &    !
           lims, lims_rvs, lims_ldc, &            !prior limits
           lims_p, lims_p_rvs, lims_p_ldc, &      !physical limits
           n_cad, t_cad, &                       !cadence cotrols
           npl, n_tel, & !integers                !planets and telescopes
           size_rv, size_tr &                     !data sizes
           )
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npl, n_tel !size of RV and LC data
  integer, intent(in) :: nwalks, maxi, thin_factor, nconv, n_cad
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_tr-1) :: plab_tr
  integer, intent(in), dimension(0:size_rv-1) :: tlab
  double precision, intent(in), dimension(0:3) :: stellar_pars
  double precision, intent(in), dimension(0:8*npl - 1) :: pars
  double precision, intent(in), dimension(0:2*8*npl - 1):: lims, lims_p
  double precision, intent(in), dimension(0:n_tel - 1) :: rvs
  double precision, intent(in), dimension(0:1) :: ldc
  double precision, intent(in), dimension(0:2*n_tel - 1) :: lims_rvs, lims_p_rvs
  double precision, intent(in), dimension(0:3) :: lims_ldc, lims_p_ldc
  double precision, intent(in) ::  t_cad
  integer, intent(in) :: wtf_trends(0:1)
  integer, intent(in) :: wtf_all(0:8*npl-1), wtf_rvs(0:n_tel-1), wtf_ldc(0:1)
  logical, intent(in) :: flags(0:5), total_fit_flag(0:1) !CHECK THE SIZE
  logical, intent(in) :: afk(0:npl-1), is_jit(0:1)
!Local variables
  double precision, dimension(0:nwalks-1,0:8*npl-1) :: pars_old, pars_new
  double precision, dimension(0:nwalks-1,0:n_tel-1) :: rvs_old, rvs_new
  double precision, dimension(0:nwalks-1,0:1) :: ldc_old, ldc_new
  double precision, dimension(0:nwalks-1) :: r_rand, z_rand, mstar, rstar
  double precision, dimension(0:nwalks-1) :: chi2_old_total, chi2_new_total, chi2_red
  double precision, dimension(0:nwalks-1) :: chi2_old_rv, chi2_old_tr
  double precision, dimension(0:nwalks-1) :: chi2_new_rv, chi2_new_tr
  double precision, dimension(0:nwalks-1,0:8*npl-1,0:nconv-1) :: pars_chains
  double precision, dimension(0:nwalks-1) :: jitter_rv_old, jitter_rv_new, jitter_tr_old, jitter_tr_new
  double precision, dimension(0:nwalks-1,0:1) :: tds_old, tds_new !linear and quadratic terms
  double precision, dimension(0:3) :: lims_trends
  double precision, dimension(0:nwalks-1,0:size_rv+size_tr-1) :: mult_old, mult_new, mult_total
  integer, dimension(0:nwalks-1) :: r_int
  double precision  :: q, a_factor, dof, tds
  double precision  :: lims_ldc_dynamic(0:1), lims_e_dynamic(0:1,0:npl-1)
  double precision  :: mstar_mean, mstar_sigma, rstar_mean, rstar_sigma
  double precision  :: jrrv, jrtr
  logical :: continua, is_burn, is_limit_good, is_cvg
  logical :: is_kepler
  integer :: nk, j, n, m, o, n_burn, spar, spar1
!external calls
  external :: init_random_seed, find_chi2_tr, find_chi2_rv

  !Get the stellar parameters
  mstar_mean  = stellar_pars(0)
  mstar_sigma = stellar_pars(1)
  rstar_mean  = stellar_pars(2)
  rstar_sigma = stellar_pars(3)

  spar = sum(wtf_all) + sum(wtf_ldc) + sum(wtf_rvs) + sum(wtf_trends)
  !spar -> size of parameters, dof -> degrees of freedom
  if ( is_jit(0) ) spar = spar + 1
  if ( is_jit(1) ) spar = spar + 1
  spar1 = spar - 1

  dof = 0
  if ( total_fit_flag(0) ) dof = dof + size_rv
  if ( total_fit_flag(1) ) dof = dof + size_tr
  dof  = dof - spar

  !call the random seed
  print *, 'CREATING RANDOM SEED'
  call init_random_seed()

  !Jitter vars
  jitter_rv_old(:) = 0.0d0
  jitter_tr_old(:) = 0.0d0
  jitter_rv_new(:) = 0.0d0
  jitter_tr_new(:) = 0.0d0
  if ( is_jit(0) ) &
     call gauss_random_bm(e_rv(0)*1e-1,e_rv(0)*1e-2,jitter_rv_old,nwalks)
  if ( is_jit(1) ) &
    call gauss_random_bm(e_tr(0)*1e-1,e_tr(0)*1e-2,jitter_tr_old,nwalks)
  mult_old(:,:) = 1.0d0
  mult_new(:,:) = 1.0d0
  if ( is_jit(0) .or. is_jit(1) ) then
    do nk = 0, nwalks - 1
      do j = 0, size_rv-1
        mult_old(nk,j) = 1.0d0/sqrt( e_rv(j)**2 + jitter_rv_old(nk)**2  )
        if ( .not. is_jit(0) ) mult_old(nk,j) = 1.0d0
      end do
      do j = size_rv, size_rv+size_tr-1
        mult_old(nk,j) = 1.0d0/sqrt( e_tr(j-size_rv)**2 + jitter_tr_old(nk)**2)
        if ( .not. is_jit(1) ) mult_old(nk,j) = 1.0d0
      end do
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


  print *, 'CREATING PRIORS'
  call gauss_random_bm(mstar_mean,mstar_sigma,mstar,nwalks)
  call gauss_random_bm(rstar_mean,rstar_sigma,rstar,nwalks)
  !Let us create uniformative random priors
  is_limit_good = .false.
  do nk = 0, nwalks - 1
    !random parameters
      call uniform_priors(pars,8*npl,wtf_all,lims,pars_old(nk,:))
      !If we are using the parametrization for e and omega
      if ( flags(1) ) then
        do m = 0, npl-1
          lims_e_dynamic(:,m) = sqrt( 1.d0 - pars_old(nk,2+8*m)**2 )
          lims_e_dynamic(0,m) = - lims_e_dynamic(0,m)
          call uniform_priors(pars(3),1,wtf_all(3),lims_e_dynamic(:,m),pars_old(nk,3+8*m))
        end do
      end if

      call uniform_priors(tds,2,wtf_trends,lims_trends,tds_old(nk,:))

      call uniform_priors(rvs,n_tel,wtf_rvs,lims_rvs,rvs_old(nk,:))

      call uniform_priors(ldc(0),1,wtf_ldc(0),lims_ldc(0:1),ldc_old(nk,0))
      lims_ldc_dynamic(:) = (/ lims_ldc(2), 1.0d0 - ldc_old(nk,0) /)

      call uniform_priors(ldc(1),1,wtf_ldc(1),lims_ldc_dynamic,ldc_old(nk,1))

      !Will we use spectroscopic priors for some planets?
      do m = 0, npl - 1
        if ( afk(m) ) then
          !The parameter comes from 3rd Kepler law !pars_old(1) is the period
          call get_a_scaled(mstar(nk),rstar(nk),pars_old(nk,1+8*m),pars_old(nk,5+8*m),1)
        end if
      end do

      call get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab,plab_tr, &
           total_fit_flag,flags,t_cad,n_cad,pars_old(nk,:),rvs_old(nk,:), &
           ldc_old(nk,:),tds_old(nk,:),jitter_rv_old(nk),jitter_tr_old(nk),&
           chi2_old_rv(nk),chi2_old_tr(nk),npl,n_tel,size_rv,size_tr)

           chi2_old_total(nk) = chi2_old_rv(nk) + chi2_old_tr(nk)


  end do

 chi2_red(:) = chi2_old_total(:) / dof

  !Print the initial cofiguration
  print *, ''
  print *, 'STARTING MCMC CALCULATION'
  print *, 'RV datapoints  = ', size_rv
  print *, 'TR datapoints  = ', size_tr
  print *, 'No. parameters = ', int(spar)
  print *, 'dof            = ', int(dof)
  call print_chain_data(chi2_red,nwalks)

  !Initialize the values
  j = 1
  n = 0
  continua = .true.
  is_burn = .false.
  a_factor = 2.d0
  n_burn = 1
  is_kepler = .false.

  !The infinite cycle starts!
  print *, 'STARTING INFINITE LOOP!'
  do while ( continua )

    !Creating the random variables
    !Do the work for all the walkers
    call random_number(z_rand)
    call random_number(r_rand)
    !Create random integers to be the index of the walks
    !Note that r_ink(i) != i (avoid copy the same walker)
    call random_int(r_int,nwalks)
    call gauss_random_bm(mstar_mean,mstar_sigma,mstar,nwalks)
    call gauss_random_bm(rstar_mean,rstar_sigma,rstar,nwalks)

    !Pick a random walker
    do nk = 0, nwalks - 1 !walkers
      pars_new(nk,:) = pars_old(r_int(nk),:)
      rvs_new(nk,:) = rvs_old(r_int(nk),:)
      ldc_new(nk,:) = ldc_old(r_int(nk),:)
      tds_new(nk,:) = tds_old(r_int(nk),:)
      jitter_rv_new(nk) = jitter_rv_old(r_int(nk))
      jitter_tr_new(nk) = jitter_tr_old(r_int(nk))
    end do

    !Paralellization calls
    !$OMP PARALLEL &
    !$OMP PRIVATE(is_limit_good,q,m,jrrv,jrtr)
    !$OMP DO SCHEDULE(DYNAMIC)
    do nk = 0, nwalks - 1

      !Generate the random step to perform the stretch move
      z_rand(nk) = z_rand(nk) * a_factor
      call find_gz(z_rand(nk),a_factor)

      !Perform the stretch move
      !Eq. (7), Goodman & Weare (2010)
      pars_new(nk,:) = pars_new(nk,:) + wtf_all(:) * z_rand(nk) * &
                     ( pars_old(nk,:) - pars_new(nk,:) )
      rvs_new(nk,:)  = rvs_new(nk,:) + wtf_rvs(:) * z_rand(nk) * &
                     ( rvs_old(nk,:) - rvs_new(nk,:) )
      ldc_new(nk,:)  = ldc_new(nk,:) + wtf_ldc(:) * z_rand(nk) * &
                     ( ldc_old(nk,:) - ldc_new(nk,:) )
      tds_new(nk,:)  = tds_new(nk,:) + wtf_trends(:) * z_rand(nk) * &
                     ( tds_old(nk,:) - tds_new(nk,:) )
      jitter_rv_new(nk) = jitter_rv_new(nk) + z_rand(nk) * &
                         ( jitter_rv_old(nk) - jitter_rv_new(nk) )
      jitter_tr_new(nk) = jitter_tr_new(nk) + z_rand(nk) * &
                         ( jitter_tr_old(nk) - jitter_tr_new(nk) )

      !For the fixed parameters, we assume a gaussian distribution with error bars
      do m = 0, 8*npl-1
       if ( wtf_all(m) == 0 ) then
          call gauss_random_bm(pars(m),pars(m)-lims(2*m),pars_new(nk,m),1)
       end if
      end do
      do m = 0, 1
        if ( wtf_ldc(m) == 0 ) call gauss_random_bm(ldc(m),ldc(m)-lims_ldc(m+2),ldc_new(nk,m),1)
      end do

      do m = 0, npl - 1
        if ( afk(m) ) then
          !The parameter comes from 3rd Kepler law
         !pars_old(1) is the period
        call get_a_scaled(mstar(nk),rstar(nk),pars_new(nk,1+8*m),pars_new(nk,5+8*m),1)
        end if
      end do

      !Let us check if the new parameters are inside the limits
      is_limit_good = .true.
      call check_limits(pars_new(nk,:),lims_p,is_limit_good,8*npl)
      if (is_limit_good ) call check_limits(rvs_new(nk,:),lims_p_rvs,is_limit_good,n_tel)
      if (is_limit_good ) call check_limits(ldc_new(nk,:),lims_p_ldc,is_limit_good,2)
      if ( is_limit_good ) then
        if ( jitter_rv_new(nk) < 0.0d0 .or. jitter_tr_new(nk) < 0.0d0 ) is_limit_good = .false.
      end if

      chi2_new_total(nk) = huge(0.0d0) !A really big number!
      if ( is_limit_good ) then !If we are inside the limits, let us calculate chi^2
        call get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab,plab_tr, &
             total_fit_flag,flags, t_cad,n_cad,pars_new(nk,:),rvs_new(nk,:), &
             ldc_new(nk,:),tds_new(nk,:),jitter_rv_new(nk),jitter_tr_new(nk), &
             chi2_new_rv(nk),chi2_new_tr(nk),npl,n_tel,size_rv,size_tr)

        chi2_new_total(nk) = chi2_new_rv(nk) + chi2_new_tr(nk)
     end if

      mult_new(nk,:) = 1.0d0
      q = 1.0d0
      if ( is_jit(0) .or. is_jit(1) ) then

        jrrv = jitter_rv_new(nk)
        jrrv = jrrv*jrrv
        jrtr = jitter_tr_new(nk)
        jrtr = jrtr*jrtr
        do m = 0, size_rv-1
          mult_new(nk,m) = 1.0d0/sqrt( e_rv(m)**2 + jrrv  )
          if ( .not. is_jit(0) ) mult_new(nk,m) = 1.0d0
        end do
        do m = size_rv, size_rv+size_tr-1
          mult_new(nk,m) = 1.0d0/sqrt( e_tr(m-size_rv)**2 + jrtr  )
          if ( .not. is_jit(1) ) mult_new(nk,m) = 1.0d0
        end do

        mult_total(nk,:) = mult_new(nk,:) / mult_old(nk,:)

        !Add the jitter terms
        do m = 0, size_rv+size_tr-1
          q = q * mult_total(nk,m)
        end do

      end if

      !Let us compare our models
      !Compute the likelihood
      q = q * z_rand(nk)**spar1 * &
              exp( ( chi2_old_total(nk) - chi2_new_total(nk) ) * 0.5d0  )

      !Check if the new likelihood is better
      if ( q > r_rand(nk) ) then
        !If yes, let us save it as the old vectors
        chi2_old_total(nk) = chi2_new_total(nk)
        chi2_old_rv(nk) = chi2_new_rv(nk)
        chi2_old_tr(nk) = chi2_new_tr(nk)
        pars_old(nk,:) = pars_new(nk,:)
        rvs_old(nk,:) = rvs_new(nk,:)
        ldc_old(nk,:) = ldc_new(nk,:)
        tds_old(nk,:) = tds_new(nk,:)
        jitter_rv_old(nk) = jitter_rv_new(nk)
        jitter_tr_old(nk) = jitter_tr_new(nk)
        mult_old(nk,:) = mult_new(nk,:)
      end if

      !Compute the reduced chi square
      chi2_red(nk) = chi2_old_total(nk) / dof

      !Start to burn-in
      if ( is_burn ) then
        if ( mod(j,thin_factor) == 0 ) then
         ! !$OMP CRITICAL
         write(101,*) j, nk, chi2_old_rv(nk),chi2_old_tr(nk), pars_old(nk,:), &
                      ldc_old(nk,:), rvs_old(nk,:)
         if ( is_jit(0) .or.  is_jit(1) ) &
             write(201,*) jitter_rv_old(nk), jitter_tr_old(nk)
         if ( sum(wtf_trends)  > 0 ) &
            write(301,*) tds_old(nk,:)
         ! !$OMP END CRITICAL
        end if
      end if
      !End burn-in

    end do !walkers
    !$OMP END PARALLEL

    !Evolve burning-in controls
    if ( is_burn ) then
      !HERE EVOLVE BURNING-IN CONTROLS
      !if ( mod(j,thin_factor) == 0 ) print *, n_burn, 'of', nconv
      if ( mod(j,thin_factor) == 0 ) n_burn = n_burn + 1
      if ( n_burn > nconv ) continua = .false.
    else
      if ( mod(j,thin_factor) == 0 ) then
        !If the chains have not converged, let us check convergence
        !Let us save a 3D array with the informations of the parameters,
        !the nk and the iteration. This array is used to perform GR test
        pars_chains(:,:,n) = pars_old(:,:)
        n = n + 1
        !Is it time to check covergence=
        if ( n == nconv ) then
          !Perform G-R test
          n = 0 !reinitilize n
          call print_chain_data(chi2_red,nwalks)
          print *, '==========================='
          print *, '  PERFOMING GELMAN-RUBIN'
          print *, '   TEST FOR CONVERGENCE'
          print *, '==========================='
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
            print *, '==========================='
            print *, '  CHAINS HAVE CONVERGED'
            print *, '==========================='
            print *, 'STARTING BURN-IN PHASE'
            print *, '==========================='
            is_burn = .True.
            print *, 'CREATING BURN-IN DATA FILE'
            !Let us start the otput file
            open(unit=101,file='all_data.dat',status='unknown')
            open(unit=201,file='jitter_data.dat',status='unknown')
            open(unit=301,file='trends_data.dat',status='unknown')
          end if ! is_cvg
        end if !nconv
      end if !j/thin_factor
    end if !is_burn

    !check if we exceed the maximum number of iterations
    if ( j > maxi ) then
      print *, 'Maximum number of iteration reached!'
      continua = .FALSE.
    end if

  j = j + 1

  end do !infinite loop
  !the MCMC part has ended

  !Close file
  close(101)
  close(201)
  close(301)

end subroutine
