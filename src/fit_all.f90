subroutine get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab, &
           t_cad,n_cad,pars,rvs,ldc,chi2,npl,n_tel,size_rv,size_tr)

implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npl, n_tel, n_cad !size of RV and LC data
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: tlab
  double precision, intent(in) :: pars(0:8*npl-1), rvs(0:n_tel-1), ldc(0:1)
  double precision, intent(in) :: t_cad
  double precision, intent(out) :: chi2
  !pars = T0, P, e, w, b, a/R*, Rp/R*, K -> for each parameter
  !THIS DOES NOT ADD alpha and beta
!Local variables
  double precision :: pars_rv(0:7+n_tel-1,0:npl-1)
  double precision :: pars_tr(0:6,0:npl-1)
  double precision :: chi2_rv, chi2_tr
  logical :: flag_rv(0:3), flag_tr(0:3)
  integer :: i, j


  !Create the parameter variables for rv and tr
  do i = 0, npl - 1
    j = 8*i !the new data start point
    pars_tr(:,i)   = pars(j:j+6)
    pars_rv(0:3,i) = pars(j:j+3)
    pars_rv(4,i) = pars(j+7)
    pars_rv(5:6,i) = 0.d0
    !ADD alpha and beta
    pars_rv(7:7+n_tel-1,i) = rvs(:)
  end do

  flag_rv = (/ .false., .false., .false., .false. /)
  flag_tr = (/ .false., .false., .true. , .false. /)

  !Let us calculate chi2
  call find_chi2_tr_new(x_tr,y_tr,e_tr,pars_tr,0.d0,flag_tr,&
                        ldc,n_cad,t_cad,chi2_tr,size_tr,npl)
  call find_chi2_rv(x_rv,y_rv,e_rv,tlab,pars_rv,0.d0,&
                    flag_rv,chi2_rv,size_rv,n_tel,npl)

  chi2 = chi2_rv + chi2_tr

end subroutine

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
implicit none

!In/Out variables
  integer, intent(in) :: size_rv, size_tr, npl, n_tel !size of RV and LC data
  integer, intent(in) :: nwalks, maxi, thin_factor, nconv, n_cad
  double precision, intent(in), dimension(0:size_rv-1) :: x_rv, y_rv, e_rv
  double precision, intent(in), dimension(0:size_tr-1) :: x_tr, y_tr, e_tr
  integer, intent(in), dimension(0:size_rv-1) :: tlab
  double precision, intent(in), dimension(0:8*npl - 1) :: pars
  double precision, intent(in), dimension(0:2*8*npl - 1):: lims, lims_p
  double precision, intent(in), dimension(0:n_tel - 1) :: rvs
  double precision, intent(in), dimension(0:1) :: ldc
  double precision, intent(in), dimension(0:2*n_tel - 1) :: lims_rvs, lims_p_rvs
  double precision, intent(in), dimension(0:3) :: lims_ldc, lims_p_ldc
  double precision, intent(in) ::  t_cad
  integer, intent(in) :: wtf_all(0:8*npl-1), wtf_rvs(0:n_tel-1), wtf_ldc(0:1)
  logical, intent(in) :: flags(0:5) !CHECK THE SIZE
!Local variables
  double precision, dimension(0:nwalks-1,0:8*npl-1) :: pars_old, pars_new
  double precision, dimension(0:nwalks-1,0:n_tel-1) :: rvs_old, rvs_new
  double precision, dimension(0:nwalks-1,0:1) :: ldc_old, ldc_new
  double precision, dimension(0:nwalks-1) :: r_rand, z_rand
  double precision, dimension(0:nwalks-1) :: chi2_old_total, chi2_new_total, chi2_red
  double precision, dimension(0:nwalks-1,0:8*npl-1,0:nconv-1) :: pars_chains
  integer, dimension(0:nwalks-1) :: r_int
  double precision  :: q, spar, a_factor, dof
  double precision  :: mean_chi2
  logical :: continua, is_burn, is_limit_good, is_cvg
  integer :: nk, j, n, m, o, n_burn, new_thin_factor
!external calls
  external :: init_random_seed, find_chi2_tr, find_chi2_rv

  !spar -> size of parameters, dof -> degrees of freedom
  spar = sum(wtf_all) + sum(wtf_ldc) + sum(wtf_rvs)
  dof  = size_rv + size_tr - spar

  !call the random seed
  print *, 'CREATING RANDOM SEED'
  call init_random_seed()

  print *, 'CREATING UNIFORM UNIFORMATIVE PRIORS'
  !Let us create uniformative random priors
  is_limit_good = .false.
  do nk = 0, nwalks - 1
    !random parameters
      call uniform_priors(pars,8*npl,wtf_all,lims,pars_old(nk,:))
      call uniform_priors(rvs,n_tel,wtf_rvs,lims_rvs,rvs_old(nk,:))
      call uniform_priors(ldc,2,wtf_ldc,lims_ldc,ldc_old(nk,:))
      call get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab, &
         t_cad,n_cad,pars_old(nk,:),rvs_old(nk,:),ldc_old(nk,:), &
         chi2_old_total(nk),npl,n_tel,size_rv,size_tr)
  end do

 chi2_red(:) = chi2_old_total(:) / dof

  !Print the initial cofiguration
  print *, ''
  print *, 'STARTING MCMC CALCULATION'
  print *, 'dof = ', dof
  call print_chain_data(chi2_red,nwalks)

  !Initialize the values
  j = 1
  n = 0
  continua = .TRUE.
  is_burn = .FALSE.
  !is_jitter = jit
  a_factor = 2.d0
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
    !$OMP PARALLEL &
    !$OMP PRIVATE(is_limit_good,q)
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

      !Let us check if the new parameters are inside the limits
      is_limit_good = .true.
      call check_limits(pars_new(nk,:),lims_p,is_limit_good,8*npl)
      if (is_limit_good ) call check_limits(rvs_new(nk,:),lims_p_rvs,is_limit_good,n_tel)
      if (is_limit_good ) call check_limits(ldc_new(nk,:),lims_p_ldc,is_limit_good,2)
      !CHECK ALSO THE LIMITS FOR RV AND LDC

      if ( is_limit_good ) then
        call get_total_chi2(x_rv,y_rv,x_tr,y_tr,e_rv,e_tr,tlab, &
             t_cad,n_cad,pars_new(nk,:),rvs_new(nk,:),ldc_new(nk,:), &
             chi2_new_total(nk),npl,n_tel,size_rv,size_tr)
      else
         chi2_new_total(nk) = huge(0.0d0) !A really big number!
      end if

    !Let us compare our models
    !Compute the likelihood
    q = 1.0d0 !add jitter later
    q = q * z_rand(nk)**( int(spar - 1) ) * &
          exp( ( chi2_old_total(nk) - chi2_new_total(nk) ) * 0.5d0  )

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
    chi2_red(nk) = chi2_old_total(nk) / dof

    !Start to burn-in
    if ( is_burn ) then
      if ( mod(j,new_thin_factor) == 0 ) then
       ! !$OMP CRITICAL
       !write(*,*)   j, nk, chi2_old_total(nk), pars_old(nk,:), &
       !             ldc_old(nk,:), rvs_old(nk,:)
       write(101,*) j, nk, chi2_old_total(nk), pars_old(nk,:), &
                    ldc_old(nk,:), rvs_old(nk,:)
       ! !$OMP END CRITICAL
      end if
    end if
    !End burn-i


    end do !walkers
    !$OMP END PARALLEL

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
          print *, 'CREATING BURNING-IN DATA FILE'
          !Let us start the otput file
          open(unit=101,file='all_data.dat',status='unknown')
        end if ! is_cvg
      end if !nconv
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

end subroutine
