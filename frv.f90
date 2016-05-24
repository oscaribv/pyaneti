!------------------------------------------------------------
!                         frv.f90
! This file contains subroutines to calculate Marcov Chain 
! Monte Carlo simulations in order to obtain planet parameters
! The subroutines can be called from python by using f2py
! They also can be used in a fortran program
!              Date --> Feb  2016, Oscar Barragán
!------------------------------------------------------------

!-----------------------------------------------------------
! This subroutine computes the circular radial velocity 
! curve for a set of values t. The subroutine returns 
! a vector (rv of the same size that t) by solving:
!  $ rv = rv0 - k [ sin ( 2* \pi ( t - t_0) / P ) ] $
! Where the parameters are the typical for a RV curve
!------------------------------------------------------------
subroutine rv_circular(t,rv0,t0,k,P,rv,ts)
implicit none

!In/Out variables
  integer, intent(in) :: ts
  double precision, intent(in), dimension(0:ts-1)  :: t
  double precision, intent(out), dimension(0:ts-1) :: rv
  double precision, intent(in) :: k, rv0, t0, P
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626d0
!
  rv(:) = rv0 - k * sin( 2.d0*pi*( t(:) - t0) / P )

end subroutine

!-----------------------------------------------------------
! This subroutine computes the radial velocity 
! curve for an eccentric orbit, given a set of values t. 
!The subroutine returns  a vector (rv of the same size that t)
! by solving:
!  $ rv = rv0 + k [ cos ( \theta + \omega ) 
!             + e * cos ( \omega ) ] $
!  Where the parameters are the typical for a RV curve
!------------------------------------------------------------
subroutine rv_curve(t,rv0,t0,k,P,e,w,rv,ts)
implicit none

!In/Out variables
  integer, intent(in) :: ts
  double precision, intent(in), dimension(0:ts-1)  :: t
  double precision, intent(in) :: rv0, t0, k, P, e, w
  double precision, intent(out), dimension(0:ts-1) :: rv
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626d0
  double precision, dimension(0:ts-1) :: ta
!External function
  external :: find_anomaly
!
  !Obtain the eccentric anomaly by using find_anomaly
  call find_anomaly(t,t0,e,w,P,ta,ts)

  rv(:) = rv0 + k * ( cos( ta(:) + w ) + e * cos(w) )
  
end subroutine

!-----------------------------------------------------------
!CHANGE THIS HEADER
! This subroutine computes the radial velocity for a multiplanet system
! curve for an eccentric orbit, given a set of values t. 
!The subroutine returns  a vector (rv of the same size that t)
! by solving:
!  $ rv = rv0 + k [ cos ( \theta + \omega ) 
!             + e * cos ( \omega ) ] $
!  Where the parameters are the typical for a RV curve
!------------------------------------------------------------
subroutine rv_curve_mp(t,rv0,t0,k,P,e,w,rv,ts,npl)
implicit none

!In/Out variables
  integer, intent(in) :: ts, npl
  double precision, intent(in), dimension(0:ts-1)  :: t
  double precision, intent(out), dimension(0:ts-1) :: rv
  double precision, intent(in), dimension(0:npl-1) :: k, t0, P, e, w
  !Here rv0 depends on the telescope systemic not the planets
  double precision, intent(in) :: rv0
!Local variables
  double precision, dimension(0:ts-1) :: ta
  integer :: i
!External function
  external :: find_anomaly
!

  rv(:) = rv0
  do i = 0, npl-1
   !Obtain the true anomaly by using find_anomaly
   !each i is for a different planet
   call find_anomaly(t,t0(i),e(i),w(i),P(i),ta,ts)
   rv(:) = rv(:) + k(i) * ( cos(ta(:) + w(i) ) + e(i) * cos(w(i)) )
   ta(:) = 0.0d0
  end do
 
end subroutine

!-----------------------------------------------------------
! This routine calculates the chi square for a RV curve
! given a set of xd-yd data points
! It takes into acount the possible difference in systematic
! velocities for different telescopes.
!Input parameters are:
! xd, yd, errs -> set of data to fit (array(datas))
! tlab -> Telescope labels (array of integers number)
! rv0  -> array for the different systemic velocities,
!         its size is the number of telescopes
! k, ec, w, t0, P -> typical planet parameters
! datas, nt -> sizes of xd,yd, errs (datas) and rv0(nt)  
!Output parameter:
! chi2 -> a double precision value with the chi2 value
!-----------------------------------------------------------
subroutine find_chi2_rv(xd,yd,errs,tlab,params,flag,chi2,datas,nt,npl)
implicit none

!In/Out variables
  integer, intent(in) :: datas, nt, npl
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1) :: tlab
  double precision, intent(in), dimension(0:4+nt,0:npl-1) :: params
  logical, intent(in)  :: flag(0:3)
  double precision, intent(out) :: chi2
!Local variables
  double precision, dimension(0:npl-1) :: t0, P, e, w, k
  double precision, dimension(0:nt-1)  :: rv0
  double precision, parameter :: pi = 3.1415926535897932384626d0
  double precision, dimension(0:datas-1) :: model, res
  integer :: i, tel
!External function
  external :: rv_circular, rv_curve_mp

  t0(:)  = params(0,:)
  P(:)   = params(1,:)
  e(:)   = params(2,:)
  w(:)   = params(3,:)
  k(:)   = params(4,:)
  rv0(:) = params(5:4+nt,0) 

  if ( flag(0) ) P(:) = 1.d1**params(1,:)
  if ( flag(1) ) then
    e(:) = params(2,:) * params(2,:) + params(3,:) * params(3,:)
    w(:) = atan2( params(2,:),params(3,:) ) 
  end if
  if ( flag(2) ) k(:) = 1.d1**params(4,:)
  if ( flag(3) ) rv0(:) = 1.d1**params(5:4+nt,0)

  !Telescope label counter
  tel = 0

  do i = 0, datas-1
   if ( tel .ne. tlab(i)  ) tel = tel + 1 
    call rv_curve_mp(xd(i),rv0(tel),t0,k,P,e,w,model(i),1,npl)
  end do

  !Let us obtain chi^2
  res(:) = ( model(:) - yd(:) ) / errs(:) 
  chi2 = dot_product(res,res)

end subroutine

!-----------------------------------------------------------
! This routine calculates a MCMC chain by using metropolis-
! hasting algorithm
! given a set of xd-yd data points
!Input parameters are:
! xd, yd, errs -> set of data to fit (array(datas))
! tlab -> Telescope labels (array of integers number)
! rv0  -> array for the different systemic velocities,
!         its size is the number of telescopes
! k, ec, w, t0, P -> typical planet parameters
! prec -> precision of the step size
! maxi -> maximum number of iterations
! thin factor -> number of steps between each output
! chi2_toler -> tolerance of the reduced chi^2 (1 + chi2_toler)
! datas, nt -> sizes of xd,yd, errs (datas) and rv0(nt)  
!Output parameter:
! This functions outputs a file called mh_rvfit.dat
!-----------------------------------------------------------
subroutine metropolis_hastings_rv(xd,yd,errs,tlab,params,prec,maxi, &
thin_factor,wtf,flag,nconv,datas,nt)
implicit none

!In/Out variables
  integer, intent(in) :: maxi, thin_factor, datas, nt, nconv
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1)  :: tlab
  double precision, intent(inout), dimension(0:4+nt) :: params
  !f2py intent(in,out)  :: params
  double precision, intent(in), dimension(0:5) :: wtf
  double precision, intent(in)  :: prec
  logical, intent(in) :: flag(0:3)
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision :: chi2_old, chi2_new, chi2_red
  double precision, dimension(0:4+nt) :: params_new
  double precision, dimension(0:nconv-1) :: chi2_vec, x_vec
  double precision  :: q, chi2_y, chi2_slope, toler_slope
  double precision  :: prec_init, esin, ecos
  integer :: j, nu, n
  logical :: get_out
  !double precision, dimension(0:5+nt) :: r
  !Let us add a plus random generator
  double precision, dimension(0:6+nt) :: r
!external calls
  external :: init_random_seed, find_chi2_rv


  if ( flag(0) ) params(1) = log10(params(1))
  if ( flag(1) ) then
    esin = sqrt(params(2)) * sin(params(3))
    ecos = sqrt(params(2)) * cos(params(3))
    params(2) = esin
    params(3) = ecos
  end if
  if ( flag(2) ) params(4) = log10(params(4))
  if ( flag(3) ) params(5:4+nt) = log10(params(5:4+nt))

  !Let us estimate our fist chi_2 value
  call find_chi2_rv(xd,yd,errs,tlab,params,flag,chi2_old,datas,nt,1)
  !Calculate the degrees of freedom
  nu = datas - size(params)
  !If we are fixing a circular orbit, ec and w are not used 
  chi2_red = chi2_old / nu
  !Print the initial cofiguration
  print *, ''
  print *, 'Starting MCMC calculation'
  print *, 'Initial Chi2_red= ', chi2_red,'nu =', nu

    !Let us start the otput file
  open(unit=101,file='mh_rvfit.dat',status='unknown')
  !Initialize the values

  toler_slope = prec
  j = 1
  n = 0
  get_out = .TRUE.

   !Call a random seed 
  print *, 'CREATING RANDOM SEED'
  call init_random_seed()


  !The infinite cycle starts!
  print *, 'STARTING INFINITE LOOP!'
  do while ( get_out )
    !r will contain random numbers for each variable
    !Let us add a random shift to each parameter
    call random_number(r)
    !prec_init = prec * sqrt(chi2_red) 
    !if ( chi2_red < 10.0 ) prec_init = prec
    !prec_init = prec * chi2_red * r(6+nt)
    prec_init = prec
    r(1:5+nt) = ( r(1:5+nt) - 0.5 ) * 2.0
    params_new(0:4)    = params(0:4)    + r(1:5)    * prec_init * wtf(0:4)
    params_new(5:4+nt) = params(5:4+nt) + r(6:5+nt) * prec_init * wtf(5)
    !Let us calculate our new chi2
    call find_chi2_rv(xd,yd,errs,tlab,params_new,flag,chi2_new,datas,nt,1)
    !Ratio between the models
    q = exp( ( chi2_old - chi2_new ) * 0.5  )
    !If the new model is better, let us save it
    if ( q > r(0) ) then
      chi2_old = chi2_new
      params = params_new
    end if
   !Calculate the reduced chi square
   chi2_red = chi2_old / nu
   !Save the data each thin_factor iteration
   if ( mod(j,thin_factor) == 0 ) then
    print *, 'iter ',j,', Chi2_red =', chi2_red
    write(101,*) j, chi2_old, chi2_red, params
    !Check convergence here
    chi2_vec(n) = chi2_red
    x_vec(n) = n
    n = n + 1

    if ( n == size(chi2_vec) ) then
      call fit_a_line(x_vec,chi2_vec,chi2_y,chi2_slope,nconv)
      n = 0
      !If chi2_red has not changed the last nconv iterations
      print *, abs(chi2_slope), toler_slope / (chi2_y)
      if ( abs(chi2_slope) < ( toler_slope / (chi2_y) ) ) then
        print *, 'I did my best to converge, chi2_red =', &
                  chi2_y
        get_out = .FALSE.
      end if

      if ( j > maxi ) then
        print *, 'Maximum number of iteration reached!'
        get_out = .FALSE.
      end if

    end if
      !I checked covergence
  end if
  j = j + 1
end do

  close(101)

end subroutine

!-------------------------------------------------------

!-----------------------------------------------------------
subroutine stretch_move_rv(xd,yd,errs,tlab,pars,lims,nwalks,maxi, &
thin_factor,wtf,flag,nconv,datas,nt,npl)
implicit none

!In/Out variables
  integer, intent(in) :: nwalks, maxi, thin_factor, datas, nt, nconv, npl
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1)  :: tlab
  double precision, intent(in), dimension(0:npl*(5+nt)-1) :: pars
  double precision, intent(in), dimension(0:2*npl*(5+nt)-1) :: lims
  integer, intent(in), dimension(0:6*npl-1) :: wtf
  logical, intent(in) ::  flag(0:3)
!Local variables
  double precision, dimension(0:5+nt-1,0:npl-1) :: params
  double precision, dimension(0:2*(5+nt)-1,0:npl-1) :: limits
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:nwalks-1) :: chi2_old, chi2_new, chi2_red
  double precision, dimension(0:4+nt,0:npl-1,0:nwalks-1) :: params_old, params_new
  double precision, dimension(0:4+nt,0:npl-1,0:nwalks-1,0:nconv-1) :: params_chains
  double precision  :: q
  double precision  :: esin(0:npl-1), ecos(0:npl-1), aa, chi2_red_min
  integer :: l,o,m, j, n, nu, nk, n_burn, spar, new_thin_factor
  logical :: get_out, is_burn, is_limit_good, is_cvg
  double precision :: r_rand(0:nwalks-1), z_rand(0:nwalks-1)
  integer, dimension(0:nwalks-1) :: r_int
  integer, dimension(0:5+nt-1,0:npl-1) :: wtf_all 
  double precision :: r_real
  character(len=11), dimension(0:npl-1) :: output_files
!external calls
  external :: init_random_seed, find_chi2_rv

!I have to do this automatically
  output_files(0) = 'planet1.dat'
  output_files(1) = 'planet2.dat'

  !Let us convert all the big arrays to 
  !matrix form 
  do m = 0, npl - 1
    params(:,m) = pars(m*(5+nt):(m+1)*(5+nt)-1)
    limits(:,m) = lims(m*2*(5+nt):(m+1)*2*(5+nt)-1)
    wtf_all(0:4,m) = wtf(m*6:(m+1)*4-1)
    wtf_all(5:5+nt-1,m) = wtf(5*(m+1))
  end do

  !size of parameters (only parameters to fit!)
  !The telescopes are the same for all the planets
  spar = nt
  do m = 0, npl - 1
    !what parameters are we fitting?
    spar = spar + sum(wtf_all(0:4,m))
  end do

  !Check if there are flags to evolve modified parameters

  !Period
  if ( flag(0) )  then
    params(1,:)   = log10(params(1,:))
    limits(2:3,:) = log10(limits(2:3,:))
  end if
  ! e and w
  if ( flag(1) ) then
    esin(:) = sqrt(params(2,:)) * sin(params(3,:))
    ecos(:) = sqrt(params(2,:)) * cos(params(3,:))
    params(2,:) = esin
    params(3,:) = ecos
    limits(4,:) = - sqrt(limits(5,:))
    limits(5,:) =   sqrt(limits(5,:))
    limits(6,:) = limits(4,:)
    limits(7,:) = limits(5,:)
  end if
  !k
  if ( flag(2) ) then
    params(4,:)  = log10(params(4,:))
    limits(8:9,:) = log10(limits(8:9,:))
  end if
  !rv0's
  if ( flag(3) ) then
    params(5:4+nt,:) = log10(params(5:4+nt,:))
    limits(10:2*(5+nt)-1,:) = log10(limits(10:2*(5+nt)-1,:))
  end if

  print *, 'CREATING RANDOM SEED'
  call init_random_seed()

  print *, 'CREATING RANDOM (UNIFORM) UNIFORMATIVE PRIORS'
  !Let us create uniformative random priors
  do nk = 0, nwalks - 1

    do m = 0, npl - 1

      !Counter for the limits
      j = 0

      do n = 0, 4 + nt
        if ( wtf_all(n,m) == 0 ) then
          !this parameter does not change for the planet m
          params_old(n,m,nk) = params(n,m)
        else
          call random_number(r_real)
          params_old(n,m,nk) = limits(j+1,m) - limits(j,m)
          params_old(n,m,nk) = limits(j,m) + r_real*params_old(n,m,nk) 
          !print *,n, params_old(n,m,nk), limits(j+1,m) - limits(j,m)
          !print *, params_old(2,m,nk)
        end if
          !print *, params_old(2,m,nk), limits(5,m) - limits(4,m)
        j = j + 2 !two limits for each parameter

      end do

      !Check that e < 1 for ew
      if ( flag(1) ) then
          is_limit_good = .false.
          do while ( .not. is_limit_good )
            !print *, params_old(2,m,nk), params_old(3,m,nk)
            call check_e(params_old(2,m,nk),params_old(3,m,nk),limits(5,m)**2,is_limit_good)
            if ( .not. is_limit_good  ) then
              params_old(2,m,nk) = params_old(2,m,nk) * params_old(2,m,nk)
              params_old(3,m,nk) = params_old(3,m,nk) * params_old(3,m,nk)
            end if
        end do
      end if

    end do
        
    !Each walker is a point in a parameter space
    !Each point contains the information of all the planets
    !Let us estimate our first chi_2 value for each walker
    call find_chi2_rv(xd,yd,errs,tlab,params_old(:,:,nk), &
                      flag,chi2_old(nk),datas,nt,npl)
  end do

  !Calculate the degrees of freedom
  nu = datas - spar
  if ( nu <= 0 ) then
    print *, 'Your number of parameters is larger than your datapoints!'
    print *, 'I cannot fit that!'
    stop
  end if
  !If we are fixing a circular orbit, ec and w are not used 
  chi2_red = chi2_old / nu
  !Print the initial cofiguration
  print *, ''
  print *, 'Starting stretch move MCMC calculation'
  print *, 'Initial Chi2_red= ', sum(chi2_red) / nwalks ,'DOF =', nu

  !Let us start the otput files
  do m = 0, npl - 1 
    open(unit=m,file=output_files(m),status='unknown')
  end do

  !Initialize the values
  j = 1
  n = 0
  get_out = .TRUE.
  is_burn = .FALSE.
  aa = 2.0d0 !this is suggested by the original paper
  n_burn = 1

  !The infinite cycle starts!
  print *, 'STARTING INFINITE LOOP!'
  do while ( get_out )

    !Do the work for all the walkers
    call random_number(z_rand)
    call random_number(r_rand)
    !Create random integers to be the index of the walks
    !Note that r_ink(i) != i (avoid copy the same walker)
    call random_int(r_int,nwalks)

    call random_number(aa)
    aa = 1.d0 + 9.9d1 * aa

    do nk = 0, nwalks - 1 !walkers
      do m = 0, npl - 1
        params_new(:,m,nk) = params_old(:,m,r_int(nk))
      end do
    end do

    do nk = 0, nwalks - 1 !walkers

    !Draw the random walker nk, from the complemetary walkers
      z_rand(nk) = z_rand(nk) * aa 

      !The gz function to mantain the affine variance codition in the walks
      call find_gz(z_rand(nk),aa) 

      !Evolve for all the planets for all the parameters
      do m = 0, npl - 1
!        params_new(:,m,nk) = params_old(:,m,r_int(nk))
        params_new(:,m,nk) = params_new(:,m,nk) + wtf_all(:,m) * z_rand(nk) * &
                           ( params_old(:,m,nk) - params_new(:,m,nk) )

        !Let us check the limits
        call check_limits(params_new(:,m,nk),limits(:,m), &
                          is_limit_good,5+nt)

        !Check that e < 1 for ew
        if ( flag(1) ) then
          call check_e(params_new(2,m,nk),params_new(3,m,nk),dble(1.0), is_limit_good )
        end if          

        !If we are out of limits this chain is bad
        if ( .not. is_limit_good ) exit

      end do !number of planets loop

      if ( is_limit_good ) then !evaluate chi2
        call find_chi2_rv(xd,yd,errs,tlab,params_new(:,:,nk),flag,chi2_new(nk),datas,nt,npl)
      else !we do not have a good model
        chi2_new(nk) = huge(dble(0.0)) !a really big number
      end if

      !Compare the models 
      q = z_rand(nk)**(spar - 1) * &
          exp( ( chi2_old(nk) - chi2_new(nk) ) * 0.5  )

      if ( q >= r_rand(nk) ) then !is the new model better?
          chi2_old(nk) = chi2_new(nk)
          params_old(:,:,nk) = params_new(:,:,nk)
      end if

      chi2_red(nk) = chi2_old(nk) / nu

      !Start to burn-in 
      if ( is_burn ) then
        if ( mod(j,new_thin_factor) == 0 ) then
          do m = 0, npl - 1 !Print a file with data of each planet 
            write(m,*) j, chi2_old(nk), chi2_red(nk), params_old(:,m,nk)
          end do
        end if
      end if
      !End burn-in

    end do !walkers

     if ( is_burn ) then

        if ( mod(j,new_thin_factor) == 0 ) then
            print *, 'Iter ',j
            n_burn = n_burn + 1
        end if

        if ( n_burn > nconv ) get_out = .false.

     else

      if ( mod(j,thin_factor) == 0 ) then

        !Obtain the chi2 mean of all the variables
        chi2_red_min = sum(chi2_red) / nwalks

        print *, 'Iter ',j,', Chi^2_red =', chi2_red_min

        !Create the 4D array to use the Gelman-Rubin test
        !The first two elemets are the parameters for mp fit
        !third is the information of all chains
        !fourth is the chains each iteration
        params_chains(:,:,:,n) = params_old(:,:,:)
        
        n = n + 1

        if ( n == nconv ) then !Perform GR test

          n = 0

          print *, '==========================='
          print *, '  PERFOMING GELMAN-RUBIN'
          print *, '   TEST FOR CONVERGENCE'
          print *, '==========================='

          !Let us check convergence for all the parameters
          is_cvg = .true.
          do o = 0, 4 + nt !For all parameters 
            do l = 0, npl - 1 !for all planets
            !Do the test to the parameters that we are fitting
              if ( wtf_all(o,l) == 1 ) then
                !do the Gelman and Rubin statistics
                call gr_test(params_chains(o,l,:,:),nwalks,nconv,is_cvg)
                ! If only a chain for a given parameter does
                ! not converge is enoug to keep iterating
              end if
              if ( .not. is_cvg ) exit
            end do
            if ( .not. is_cvg ) exit
          end do

          if ( .not. is_cvg  ) then
            print *, '=================================='
            print *, 'CHAINS HAVE NOT CONVERGED YET!'
            print *,  nconv,' ITERATIONS MORE!'
            print *, '=================================='
          else
            print *, '==========================='
            print *, '  CHAINS HAVE CONVERGED'
            print *, '==========================='
            print *, 'STARTING BURNING-IN PHASE'
            print *, '==========================='
            is_burn = .True.
            new_thin_factor = 50
            do m = 0, npl - 1
              print *, 'Creating ', output_files(m)
            end do
          end if

          if ( j > maxi ) then
            print *, 'Maximum number of iteration reached!'
            get_out = .FALSE.
          end if

        end if
      !I checked covergence

      end if

    end if

  j = j + 1

  end do

  !Close all the files
  do m = 0, npl - 1 
    close(m)
  end do

end subroutine


!-----------------------------------------------------------
subroutine stretch_move(xd_rv,yd_rv,errs_rv,tlab,xd_tr,yd_tr,errs_tr,params, &
limits,limits_physical,nwalks,a_factor,maxi,thin_factor,n_cad,t_cad,wtf,flag,nconv,drv,dtr,nt)
implicit none

!In/Out variables
  integer, intent(in) :: nwalks, maxi, thin_factor,n_cad, drv, dtr, nt, nconv
  double precision, intent(in), dimension(0:drv-1)  :: xd_rv, yd_rv, errs_rv
  double precision, intent(in), dimension(0:dtr-1)  :: xd_tr, yd_tr, errs_tr
  integer, intent(in), dimension(0:drv-1)  :: tlab
  double precision, intent(inout), dimension(0:10+nt-1) :: params
  !f2py intent(in,out)  :: params
  double precision, intent(inout), dimension(0:2*(10+nt)-1) :: limits
  !f2py intent(in,out)  :: limits
  double precision, intent(inout), dimension(0:2*(10+nt)-1) :: limits_physical
  !f2py intent(in,out)  :: limits_physical
  integer, intent(in), dimension(0:10) :: wtf
  double precision, intent(in)  :: a_factor, t_cad
  logical, intent(in) :: flag(0:5)
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626d0
  double precision, dimension(0:nwalks-1) :: chi2_old_rv, &
  chi2_new_rv, chi2_old_tr, chi2_new_tr, &
  chi2_old_total, chi2_new_total, chi2_red
  double precision, dimension(0:dtr-1)  :: zr
  double precision, dimension(0:4+nt) :: params_rv
  double precision, dimension(0:8) :: params_tr
  double precision, dimension(0:10+nt-1,0:nwalks-1,0:nconv-1) :: params_chains
  double precision, dimension(0:10+nt-1,0:nwalks-1) :: params_old, params_new
  double precision  :: q, q1k, q2k
  double precision  :: esin, ecos, aa, chi2_red_min
  integer :: o,j, n, nu, nk, n_burn, spar, new_thin_factor
  logical :: get_out, is_burn, is_limit_good, flag_rv(0:3), flag_tr(0:3), is_cvg
  !Let us add a plus random generator
  double precision, dimension(0:nwalks-1) :: r_rand, z_rand
  integer, dimension(0:nwalks-1) :: r_int
  integer, dimension(0:10+nt-1) :: wtf_all 
  real :: r_real
  character (len=15) :: output_files
!external calls
  external :: init_random_seed, find_chi2_rv

  output_files = 'mh_fit.dat'

  !What are we going to fit?
  wtf_all(0:9) = wtf(0:9)
  wtf_all(10:10+nt-1) = wtf(10)

  flag_tr(:)   = flag(0:3)
  flag_rv(0:1) = flag(0:1)
  flag_rv(2:3) = flag(4:5)

  spar = sum(wtf_all)

  !Period
  if ( flag(0) )  then
    params(1) = log10(params(1))
    limits(2:3) = log10(limits(2:3))
    limits_physical(2:3) = log10(limits_physical(2:3))
  end if
  ! e and w
  if ( flag(1) ) then
    esin = sqrt(params(2)) * sin(params(3))
    ecos = sqrt(params(2)) * cos(params(3))
    params(2) = esin
    params(3) = ecos
    limits(4) = -sqrt(limits(5))
    limits(5) =  sqrt(limits(5))
    limits(6) = limits(4)
    limits(7) = limits(5)
    limits_physical(4) = -sqrt(limits_physical(5))
    limits_physical(5) =  sqrt(limits_physical(5))
    limits_physical(6) = limits_physical(4)
    limits_physical(7) = limits_physical(5)
  end if
 !i
  if ( flag(2) ) then
    params(4) = sin(params(4))
    limits(8:9) = sin(limits(8:9))
    limits_physical(8:9) = sin(limits_physical(8:9))
  end if
  !a = rp/r*
  if ( flag(3) ) then
    params(5) = log10(params(5))
    limits(10:11) = log10(limits(10:11))
    limits_physical(10:11) = log10(limits_physical(10:11))
  end if
  !u1 and u2
  !Change the  u1 and u2 limb darkening coefficints following Kipping 2013
  q1k = ( params(6) + params(7) )
  q2k = 0.5 * params(6) / q1k
  q1k = q1k*q1k
  !the limits are by default [0,1] take care with this
  params(6) = q1k
  params(7) = q2k
  limits(12) = 0.0
  limits_physical(12) = 0.0
  limits(13) = 1.0
  limits_physical(13) = 1.0
  limits(14) = 0.0
  limits_physical(14) = 0.0
  limits(15) = 1.0
  limits_physical(15) = 1.0
  !k
  if ( flag(4) ) then
    params(9) = log10(params(9))
    limits(18:19) = log10(limits(18:19))
    limits_physical(18:19) = log10(limits_physical(18:19))
  end if
  !rv0's
  if ( flag(5) ) then
    params(10:10+nt-1) = log10(params(10:10+nt-1))
    limits(20:2*(10+nt)-1) = log10(limits(20:2*(10+nt)-1))
    limits_physical(20:2*(10+nt)-1) = log10(limits_physical(20:2*(10+nt)-1))
  end if

  !Call a random seed 
  print *, 'CREATING RANDOM SEED'
  call init_random_seed()


  print *, 'CREATING RANDOM UNIFORMATIVE PRIORS'
  !Let us create uniformative random priors
  chi2_old_total(:) = huge(1.d4)
  aa = huge(1.e4)

  do o = 1, 1

  do nk = 0, nwalks - 1

    print *, 'creating walker ', nk + 1

    do while ( chi2_old_total(nk) > aa )

    j = 0

    do n = 0, 10 + nt - 1


      if ( wtf_all(n) == 0 ) then
        params_old(n,nk) = params(n)
      else
        call random_number(r_real)
        params_old(n,nk) = limits(j+1) - limits(j)
        params_old(n,nk) = limits(j) + r_real*params_old(n,nk) 
      end if
      !print *, params_old(n,nk), limits(j), limits(j+1)
      j = j + 2
    end do

   !Check that u1 + u2 < 1 for ew
!   is_limit_good = .false.
!   do while ( .not. is_limit_good )
!     call check_us(params_old(6,nk),params_old(7,nk),is_limit_good)
!     if ( .not. is_limit_good  ) then
!      params_old(6,nk) = params_old(6,nk) * params_old(6,nk)
!      params_old(7,nk) = params_old(7,nk) * params_old(7,nk)
!     end if
!   end do

   !Check that e < 1 for ew
   if ( flag(1) ) then
     is_limit_good = .false.
     do while ( .not. is_limit_good )
      call check_e(params_old(2,nk),params_old(3,nk),limits(5)**2,is_limit_good)
      if ( .not. is_limit_good  ) then
        params_old(2,nk) = params_old(2,nk) * params_old(2,nk)
        params_old(3,nk) = params_old(3,nk) * params_old(3,nk)
       end if
     end do
   end if

    params_tr = params_old(0:8,nk)
    params_rv(0:3) = params_old(0:3,nk)
    params_rv(4:4+nt) = params_old(9:9+nt,nk)
    !Find the chi2 for each case
    !call find_z(xd_tr,params_tr(0:5),flag_tr,zr,dtr)
    call find_chi2_tr(xd_tr,yd_tr,errs_tr,params_tr(0:5),flag_tr,params_tr(6:8),n_cad,t_cad,chi2_old_tr(nk),dtr)
    call find_chi2_rv(xd_rv,yd_rv,errs_rv,tlab,params_rv,flag_rv,chi2_old_rv(nk),drv,nt,1)
    chi2_old_total(nk) = chi2_old_tr(nk) + chi2_old_rv(nk)
    !print *, chi2_old_total(nk)

  end do

  end do
  
  aa = sum(chi2_old_total)/nwalks
  print *, 'mean chi2 =', aa

  end do

  !stop

  !Calculate the degrees of freedom
  nu = drv + dtr - spar
  if ( nu <= 0 ) then
    print *, 'Your number of parameters is larger than your datapoints!'
    print *, 'I cannot fit that!'
    stop
  end if

  chi2_red(:) = chi2_old_total(:) / nu

  !Print the initial cofiguration
  print *, ''
  print *, 'Starting stretch move MCMC calculation'
  print *, 'Initial Chi2_red= ', sum(chi2_red) / nwalks,'DOF =', nu

  !Let us start the otput file
  open(unit=101,file='mh_fit.dat',status='unknown')
  !Initialize the values

  j = 1
  n = 0
  get_out = .TRUE.
  is_burn = .FALSE.
  aa = a_factor
  n_burn = 1

  !The infinite cycle starts!
  print *, 'STARTING INFINITE LOOP!'
  do while ( get_out )

    !Do the work for all the walkers
    call random_number(z_rand)
    call random_number(r_rand)
    !Create random integers to be the index of the walks
    !Note that r_ink(i) != i (avoid copy the same walker)
    call random_int(r_int,nwalks)

    do nk = 0, nwalks - 1 !walkers
        params_new(:,nk) = params_old(:,r_int(nk))
    end do

    !Let us vary aa randomlly
    if ( a_factor < 1.0d0 ) then
      call random_number(aa)
      aa = 1.0d0 + thin_factor * aa 
    end if

    do nk = 0, nwalks - 1

      !Let us generate a random step
      z_rand(nk) = z_rand(nk) * aa 
      call find_gz(z_rand(nk),aa) 
    
      !Now we can have the evolved walker
      do o = 0, 9 + nt
      !params_new(o,nk) = params_new(o,nk) + wtf_all(o) * z_rand(nk) * &
      !                 ( params_old(o,nk) - params_new(o,nk) )
      params_new(o,nk) = params_new(o,nk) + z_rand(nk) * &
                       ( params_old(o,nk) - params_new(o,nk) )
      !params_new(:,nk) = params_new(:,nk) + wtf_all(:) * z_rand(nk) * &
      !                 ( params_old(:,nk) - params_new(:,nk) )
      end do

      !Let us check the limits
      !call check_limits(params_new(:,nk),limits, &
      !is_limit_good,4+nt)

        !Let us check the limits
        call check_limits(params_new(:,nk),limits_physical(:), &
                          is_limit_good,10+nt)
        if ( is_limit_good ) then
          ! u1 + u2 < 1 limit
!          call check_us(params_new(6,nk),params_new(7,nk),is_limit_good)
!          if (is_limit_good ) then
            !Check that e < 1 for ew
            if ( flag(1) ) then
              call check_e(params_new(2,nk),params_new(3,nk),dble(1.0), is_limit_good )
            end if          
!          end if
        end if

      if ( is_limit_good ) then !evaluate chi2
        !Obtain the new chi square 
        params_tr = params_new(0:8,nk)
        params_rv(0:3) = params_new(0:3,nk)
        params_rv(4:4+nt) = params_new(9:9+nt,nk)
        !Find the chi2 for each case
        !call find_z(xd_tr,params_tr(0:5),flag_tr,zr,dtr)
        call find_chi2_tr(xd_tr,yd_tr,errs_tr,params_tr(0:5),flag_tr,params_tr(6:8),n_cad,t_cad,chi2_new_tr(nk),dtr)
        call find_chi2_rv(xd_rv,yd_rv,errs_rv,tlab,params_rv,flag_rv,chi2_new_rv(nk),drv,nt,1)
        chi2_new_total(nk) = chi2_new_tr(nk) + chi2_new_rv(nk)
      else !we do not have a good model
        chi2_new_total(nk) = huge(dble(0.0))
        !print *, 'I almost collapse!'
      end if

      !Is the new model better? 
      q = z_rand(nk)**( spar - 1 ) * &
          dexp( ( chi2_old_total(nk) - chi2_new_total(nk) ) * 0.5  )

      if ( q >= r_rand(nk) ) then !is the new model better?
        chi2_old_total(nk) = chi2_new_total(nk)
        params_old(:,nk) = params_new(:,nk)
      end if

      chi2_red(nk) = chi2_old_total(nk) / nu

      !Start to burn-in 
      if ( is_burn ) then
        if ( mod(j,new_thin_factor) == 0 ) then
         write(101,*) j, chi2_old_total(nk), chi2_red(nk), params_old(:,nk)
        end if
      end if
      !End burn-in

    end do !walkers


     if ( is_burn ) then

        if ( mod(j,new_thin_factor) == 0 ) then
          print *, 'Iter ', j
           n_burn = n_burn + 1
        end if
        if ( n_burn > nconv ) get_out = .false.

     else

      if ( mod(j,thin_factor) == 0 ) then

        !Obtain the chi2 mean of all the variables
        chi2_red_min = sum(chi2_red) / nwalks

        print *, 'Iter ',j,', Chi^2_red =', chi2_red_min

        !Create the 3D array to use the Gelman-Rubin test
        !The first elemets are the parameters for transit fit
        !second is the information of all chains
        !third is the chains each iteration
        params_chains(:,:,n) = params_old(:,:)
        
        n = n + 1

        if ( n == nconv ) then !Perform GR test

          n = 0

          print *, '==========================='
          print *, '  PERFOMING GELMAN-RUBIN'
          print *, '   TEST FOR CONVERGENCE'
          print *, '==========================='

          !Let us check convergence for all the parameters
          is_cvg = .true.
          do o = 0, 10 + nt - 1 !For all parameters 
              !Do the test to the parameters that we are fitting
              if ( wtf_all(o) == 1 ) then
                !do the Gelman and Rubin statistics
                call gr_test(params_chains(o,:,:),nwalks,nconv,is_cvg)
                ! If only a chain for a given parameter does
                ! not converge is enoug to keep iterating
              end if
            if ( .not. is_cvg ) exit
          end do

          if ( .not. is_cvg  ) then
            print *, '=================================='
            print *, 'CHAINS HAVE NOT CONVERGED YET!'
            print *,  nconv,' ITERATIONS MORE!'
            print *, '=================================='
          else
            print *, '==========================='
            print *, '  CHAINS HAVE CONVERGED'
            print *, '==========================='
            print *, 'STARTING BURNING-IN PHASE'
            print *, '==========================='
            is_burn = .True.
            new_thin_factor = 50
            print *, 'Creating ', output_files
          end if

          if ( j > maxi ) then
            print *, 'Maximum number of iteration reached!'
            get_out = .FALSE.
          end if

        end if
      !I checked covergence

      end if

    end if

  j = j + 1

  end do

  close(101)

end subroutine

