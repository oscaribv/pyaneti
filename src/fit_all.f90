subroutine multi_all_stretch_move( &
           x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab, & !Data vars
           pars, rvs, ldc, &                      !parameters vars
           flags, &                               !flags
           wtf_all, wtf_rvs, wtf_ldc, &           !fitting controls
           nwalks, maxi, thin_factor, nconv, &    !
           lims, lims_rvs, lims_ldc, &            !prior limits
           lims_p, lims_p_rvs, lims_p_ldc, &      !physical limits
           n_cad, t_cad, &                        !cadence cotrols
           npl, n_tel, & !integers                !planets and telescopes
           size_rv, size_tr &                     !data sizes
           )

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npl, n_tel !size of RV and LC data
  integer, intent(in) :: nwalks, maxi, thin_factor, n_cad
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  double precision, intent(in), dimension(0:8*npl - 1) :: pars, lims, lims_p
  double precision, intent(in), dimension(0:size_rv-1) :: tlab
  double precision, intent(in), dimension(0:n_tel - 1) :: rvs, lims_rvs, lims_p_rvs
  double precision, intent(in), dimension(0:1) :: ldc, lims_ldc, lims_p_ldc
  double precision, intent(in) ::  t_cad
  integer, intent(in) :: wtf_all(0:8*npl-1), wtf_rvs(0:n_tel-1), wtf_ldc(0:1)
  logical, intent(in) :: flags(0:3) !CHECK THE SIZE
!Local variables
  double precision, dimension(0:nwalks-1,0:8*npl-1) :: pars_old, pars_new
  double precision, dimension(0:nwalks-1,0:n_tel - 1) :: rvs_old, rvs_new
  double precision, dimension(0:nwalks-1,0:1) :: ldc_old, ldc_new
  double precision, dimension(0:nwalks-1) :: r_rand, z_rand
  double precision, dimension(0:nwalks-1) :: chi2_old_total, chi2_new_total, chi2_red
  double precision, dimension(0:nwalks-1) :: chi2_old_rv, chi2_new_rv
  double precision, dimension(0:nwalks-1) :: chi2_old_tr, chi2_new_tr
  double precision, dimension(0:nwalks-1,0:8*npl-1,0:nconv-1) :: pars_chains
  integer, dimension(0:nwalks-1) :: r_int
  double precision  :: q, spar, a_factor, dof
  logical :: continua, is_burn, is_limit_good, is_cvg
  integer :: nk, j, n, m, o, n_burn, new_thin_factor
!external calls
  external :: init_random_seed, find_chi2_tr, find_chi2_rv


  !call the random seed
  print *, 'CREATING RANDOM SEED'
  call init_random_seed()

  print *, 'CREATING UNIFORM UNIFORMATIVE PRIORS'
  !Let us create uniformative random priors
  do nk = 0, nwalks - 1
    !random parameters
    call uniform_priors(pars,8*npl,wtf_all,lims,pars_old(nk,:))
    call uniform_priors(rvs,n_tel,wtf_rvs,lims_rvs,rvs_old(nk,:))
    call uniform_priors(ldc,2,wtf_ldc,lims_ldc,ldc_old(nk,:))
    !Let us calculate the initla chi2 values for each chain
    !call find_chi2_rv(x_rv,y_rv,e_rv,tlab,params_rv,0.0d0, &
    !     flag_rv,chi2_old_rv(nk),size_rv,n_tel,npl)
    !call find_chi2_tr(x_tr,y_tr,e_tr,params_tr(0:5),0.0d0, &
    !       flag_tr, params_tr(6:8),n_cad,t_cad,chi2_old_tr(nk),size_tr)
    !chi2_old_total(nk) = chi2_old_tr(nk) + chi2_old_rv(nk)
  end do

 chi2_red(:) = chi2_old_total(:) / nu !nu = dof

  !Print the initial cofiguration
  print *, ''
  print *, 'STARTING MCMC'
  print *, 'dof = ', nu
  call print_chain_data(chi2_red,nwalks)

  !Initialize the values

  j = 1
  n = 0
  continua = .TRUE.
  is_burn = .FALSE.
  is_jitter = jit
  aa = a_factor
  n_burn = 1

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

    !Pick a random walker
    do nk = 0, nwalks - 1 !walkers
      pars_new(nk,:) = pars_old(r_int(nk),:)
      rvs_new(nk,:) = rvs_old(r_int(nk),:)
      ldc_new(nk,:) = ldc_old(r_int(nk),:)
    end do

    !Paralellization calls
    ! !$OMP PARALLEL &
    ! !$OMP PRIVATE(is_limit_good,q,m,params_tr,params_rv)
    ! !$OMP DO SCHEDULE(DYNAMIC)
    do nk = 0, nwalks - 1

      !Generate the random step to perform the stretch move
      z_rand(nk) = z_rand(nk) * aa
      call find_gz(z_rand(nk),aa)

      !Perform the stretch move
         !This evolves all the parameters at the same time
      pars_new(nk,:) = pars_new(nk,:) + wtf_all(:) * z_rand(nk) * &
                     ( pars_old(nk,:) - pars_new(nk,:) )
      rvs_new(nk,:) = rvs_new(nk,:) + wtf_all(:) * z_rand(nk) * &
                    ( rvs_old(nk,:) - rvs_new(nk,:) )
      ldc_new(nk,:) = ldc_new(nk,:) + wtf_all(:) * z_rand(nk) * &
                     ( ldc_old(nk,:) - ldc_new(nk,:) )

      !Let us check if the new parameters are inside the limits
      call check_limits(pars_new,lims_p,is_limit_good,8*npl)
      !CHECK ALSO THE LIMITS FOR RV AND LDC

      if ( is_limit_good ) then
        !CALCULATE THE CHI2
      else
         chi2_new_total(nk) = huge(0.0d0) !A really big number!
      end if

    !Let us compare our models
    !Compute the likelihood
    q = 1.0d0 !add jitter later
    q = q * z_rand(nk)**( int(spar - 1) ) * &
          exp( ( chi2_old_total(nk) - chi2_new_total(nk) ) * 0.5  )

    !Check if the new likelihood is better
    if ( q > r_rand(nk) ) then
      !If yes, let us save it as the old vectors
      chi2_old_total(nk) = chi2_new_total(nk)
      pars_old(nk,:) = pars_new(nk,:)
      rvs_old(nk,:) = rvs_new(nk,:)
      ldc_old(nk,:) = ldc_new(nk,:)
      !WE CAN ADD JITTER LATER
    end if

    !Compute the reduced chi square

    chi2_red(nk) = chi2_old_total(nk) / nu

    !start to burn-in
    !ONCE THE CHAINS HAVE CONVERGED, WRITE DATA HERE

    end do !walkers
    ! !$OMP END PARALLEL

    !Evolve burning-in controls
    if ( is_burn ) then
      !HERE EVOLVE BURNING-IN CONTROLS
      if ( mod(j,new_thin_factor) == 0 ) n_burn = n_burn + 1
      if ( n_burn > nconv ) continua = .false.
    else
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
          print *, 'STARTING BURNING-IN PHASE'
          print *, '==========================='
          is_burn = .True.
          new_thin_factor = thin_factor
!          print *, 'Creating ', output_files
!          !Let us start the otput file
!          open(unit=101,file='mh_fit.dat',status='unknown')
        end if ! is_cvg
      end if !nconv
    end if !is_burn

    !check if we exceed the maximum number of iterations
    if ( j > maxi ) then
      print *, 'Maximum number of iteration reached!'
      continua = .FALSE.
    end if

  end do !infinite loop
  !the MCMC part has ended


end subroutine
