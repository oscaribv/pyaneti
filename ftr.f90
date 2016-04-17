!------------------------------------------------------------
!                         ftr.f90
! This file contains subroutines to calculate Marcov Chain 
! Monte Carlo simulations in order to obtain planet parameters
! from light curve fitting of transit planets
! The subroutines can be called from python by using f2py
! They also can be used in a fortran program
!              Date --> Feb  2016, Oscar Barragán
!------------------------------------------------------------

!-----------------------------------------------------------
!  Find z suborutine
!------------------------------------------------------------
subroutine find_z(t,pars,flag,z,ts)
implicit none

!In/Out variables
  integer, intent(in) :: ts
  double precision, intent(in), dimension(0:ts-1) :: t
  double precision, intent(in), dimension(0:5) :: pars
  double precision, intent(out), dimension(0:ts-1) :: z
  logical, intent(in), dimension(0:3) :: flag 
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:ts-1) :: ta, swt
  double precision :: t0, P, e, w, i, a
  double precision :: si, delta = 1.e-7
  integer :: imax = int(1e7)
!External function
  external :: find_anomaly
!

  t0  = pars(0)
  P   = pars(1)
  e   = pars(2)
  w   = pars(3)
  i   = pars(4)
  a   = pars(5)

  if ( flag(0) ) P = 10**pars(1)
  if ( flag(1) ) then
    e = pars(2) * pars(2) + pars(3) * pars(3)
    w = atan2(pars(2),pars(3))
  end if
  if (flag(2)) i = asin(pars(4))
  if (flag(3)) a = 10**pars(5)

  !Obtain the eccentric anomaly by using find_anomaly
  call find_anomaly(t,t0,e,w,P,ta,delta,imax,ts)
  swt = sin(w+ta)

  si = sin(i)
  z = a * ( 1. - e * e ) * sqrt( 1. - swt * swt * si * si ) &
      / ( 1.d0 + e * cos(ta) ) 
  !z has been calculated
  
end subroutine

subroutine check_eclipse(z,pz,is_good,sizez)
implicit none

!In/Out variables
  integer, intent(in) :: sizez
  double precision, intent(in), dimension(0:sizez) :: z
  double precision, intent(in) :: pz
  logical, intent(out) :: is_good
!Local variables
  integer :: i

  is_good = .false.
  !At least we have to have one eclipse condition
  do i = 0, sizez - 1
    if ( z(i) < 1.d0 + pz ) then
      is_good = .true.
      exit
    end if
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
! ics -> is circular flag, True or False.
! datas, nt -> sizes of xd,yd, errs (datas) and rv0(nt)  
!Output parameter:
! chi2 -> a double precision value with the chi2 value
!-----------------------------------------------------------
!subroutine find_chi2_tr(xd,yd,errs,t0,P,e,w,i,a,u1,u2,pz,chi2,isc,datas)
subroutine find_chi2_tr(xd,yd,errs,z,params,chi2,datas)
implicit none

!In/Out variables
  integer, intent(in) :: datas
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs, z
  double precision, intent(in), dimension (0:2) :: params
  double precision, intent(out) :: chi2
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:datas-1) :: res, muld, mu
  double precision :: u1, u2, pz
!  integer :: n
!External function
  external :: occultquad

  u1  = params(0)
  u2  = params(1)
  pz  = params(2) 

  !Now we have z, let us use Agol's routines
  call occultquad(z,u1,u2,pz,muld,mu,datas)

  !print *, muld
  !Let us calculate the residuals
  ! chi^2 = \Sum_i ( M - O )^2 / \sigma^2
  !Here I am assuming that we want limb darkening
  !If this is not true, use mu 
  res(:) = ( muld(:) - yd(:) ) / errs(:) 
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
! ics -> is circular flag, True or False.
! datas, nt -> sizes of xd,yd, errs (datas) and rv0(nt)  
!Output parameter:
! This functions outputs a file called mh_rvfit.dat
!-----------------------------------------------------------
!subroutine metropolis_hastings_tr(xd,yd,errs,params,prec,maxi,thin_factor,wtf,flag,nconv,datas)
!implicit none
!
!!In/Out variables
!  integer, intent(in) :: maxi, thin_factor, nconv, datas
!  integer, intent(in),dimension(0:8) :: wtf
!  double precision, intent(in), dimension(0:datas-1) :: xd, yd, errs
!  double precision, intent(in)  :: prec
!  !Params vector contains all the parameters
!  !t0,P,e,w,i,a,u1,u2,pz
!  double precision, intent(inout), dimension(0:8) :: params
!  !f2py intent(in,out)  :: params
!  logical, intent(in) :: flag(0:3)
!!Local variables
!  double precision, parameter :: pi = 3.1415926535897932384626
!  double precision :: chi2_old, chi2_new, chi2_red
!  double precision, dimension(0:8) :: params_new
!  double precision, dimension(0:nconv-1) :: chi2_vec, x_vec
!  double precision  :: q, chi2_y, chi2_slope, toler_slope
!  double precision :: ecos, esin, prec_init
!  integer :: j, nu, n
!  double precision, dimension(0:9) :: r
!  logical :: get_out
!!external calls
!  external :: init_random_seed, find_chi2_tr
!
!  !Period
!  if ( flag(0) ) params(1) = log10(params(1))
!  !eccentricity and periastron
!  if ( flag(1) ) then
!    esin = sqrt(params(2)) * sin(params(3))
!    ecos = sqrt(params(2)) * cos(params(3))
!    params(2) = esin
!    params(3) = ecos
!  end if
!  if (flag(2) ) params(4) = sin(params(4))
!  if (flag(3) ) params(5) = log10(params(5))
!
!  !Let us estimate our fist chi_2 value
!  !call find_chi2_tr(xd,yd,errs,t0,P,e,w,i,a,u1,u2,pz,chi2_old,ics,datas)
!  call find_chi2_tr(xd,yd,errs,params,flag,chi2_old,datas)
!  !Calculate the degrees of freedom
!  nu = datas - size(params)
!  chi2_red = chi2_old / nu
!  !Print the initial cofiguration
!  print *, ''
!  print *, 'Starting MCMC calculation'
!  print *, 'Initial Chi2_red= ',chi2_red,'nu =',nu
!
!
!  !Let us start the otput file
!  open(unit=101,file='mh_trfit.dat',status='unknown')
!
!  !Initialize the values
!  toler_slope = prec 
!  j = 1
!  n = 0
!  get_out = .TRUE.
!
!  call init_random_seed()
!
!  !The infinite cycle starts!
!  do while ( get_out )
!    !r will contain random numbers for each variable
!    !Let us add a random shift to each parameter
!    !Call a random seed 
!    call random_number(r)
!    !prec_init = prec * sqrt(chi2_red)
!    prec_init = prec
!    r(1:9) = ( r(1:9) - 0.5) * 2.
!    params_new = params + r(1:9) * prec_init * wtf
!    !print *, ''
!    !print *, params
!    !print *, params_new
!    !Let us calculate our new chi2
!    call find_chi2_tr(xd,yd,errs,params_new,flag,chi2_new,datas)
!    !print *, prec, wtf
!    !print *, chi2_old, chi2_new
!    !Ratio between the models
!    q = exp( ( chi2_old - chi2_new ) * 0.5  )
!    !If the new model is better, let us save it
!    if ( q > r(0) ) then
!      chi2_old = chi2_new
!      params = params_new
!    end if
!    !Calculate the reduced chi square
!    chi2_red = chi2_old / nu
!    !stop
!    !Save the data each thin_factor iteration
!    if ( mod(j,thin_factor) == 0 ) then
!      print *, 'iter ',j,', Chi2_red =', chi2_red
!      write(101,*) j, chi2_old, chi2_red, params
!      !Check convergence here
!      chi2_vec(n) = chi2_red
!      x_vec(n) = n
!      n = n + 1
!
!      if ( n == size(chi2_vec) ) then
!        call fit_a_line(x_vec,chi2_vec,chi2_y,chi2_slope,nconv)        
!        n = 0
!        !If chi2_red has not changed the last nconv iterations
!        print *, abs(chi2_slope), toler_slope
!        if ( abs(chi2_slope) < toler_slope ) then
!          print *, 'I did my best to converge, chi2_red =', &
!                    chi2_y
!          get_out = .FALSE.
!        end if
!      end if
!      if ( j > maxi ) then
!        print *, 'Maximum number of iteration reached!'
!        get_out = .FALSE.
!      end if
!    end if
!    !I checked covergence
!
!    j = j + 1
!  end do
!
!  close(101)
!
!end subroutine
!!
!
!subroutine metropolis_hastings(xd_rv,yd_rv,errs_rv,tlab,xd_tr,yd_tr,errs_tr,params,prec,maxi,thin_factor, &
!           wtf,flag,nconv,drv,dtr,ntl)
!implicit none
!
!!In/Out variables
!  integer, intent(in) :: maxi, thin_factor, nconv, drv, dtr, ntl
!  integer, intent(in),dimension(0:10) :: wtf
!  double precision, intent(in), dimension(0:drv-1)::xd_rv,yd_rv,errs_rv
!  integer, intent(in), dimension(0:drv-1):: tlab
!  double precision, intent(in), dimension(0:dtr-1) :: xd_tr, yd_tr, errs_tr
!  double precision, intent(in)  :: prec
!  !Params vector contains all the parameters
!  !t0,P,e,w,i,a,u1,u2,pz
!  double precision, intent(inout), dimension(0:9+ntl) :: params
!  !f2py intent(in,out)  :: params
!  logical, intent(in) :: flag(0:5)
!!Local variables
!  double precision, parameter :: pi = 3.1415926535897932384626
!  double precision :: chi2_rv_old, chi2_tr_old, chi2_rv_new, chi2_tr_new, chi2_red
!  double precision :: chi2_old_total, chi2_new_total
!  double precision, dimension(0:4+ntl) :: params_rv, params_rv_new
!  double precision, dimension(0:8) :: params_tr, params_tr_new
!  double precision, dimension(0:9+ntl) :: params_new
!  double precision, dimension(0:nconv-1) :: chi2_vec, x_vec
!  double precision  :: q, chi2_y, chi2_slope, toler_slope
!  double precision  :: esin, ecos, prec_init
!  integer :: j, nu, n
!  double precision, dimension(0:10+ntl) :: r
!  logical :: get_out, flag_rv(0:3), flag_tr(0:3)
!!external calls
!  external :: init_random_seed, find_chi2_tr
!
!
!  if ( flag(0) ) params(1) = log10(params(1))
!  if ( flag(1) ) then
!    esin = sqrt(params(2)) * sin(params(3))
!    ecos = sqrt(params(2)) * cos(params(3))
!    params(2) = esin
!    params(3) = ecos
!  end if
!  if ( flag(2) ) params(4) = sin(params(4))
!  if ( flag(3) ) params(5) = log10(params(5))
!  if ( flag(4) ) params(9) = log10(params(9))
!  if ( flag(5) ) params(10:9+ntl) = log10(params(10:9+ntl))
!
!  flag_tr = flag(0:3)
!  flag_rv(0:1) = flag(0:1)
!  flag_rv(2:3) = flag(4:5)
!
!  params_tr = params(0:8)
!  params_rv(0:3) = params(0:3)
!  params_rv(4:4+ntl) = params(9:9+ntl)
!
!  !Let us estimate our fist chi_2 value
!  call find_chi2_tr(xd_tr,yd_tr,errs_tr,params_tr,flag_tr,chi2_tr_old,dtr)
!  call find_chi2_rv(xd_rv,yd_rv,errs_rv,tlab,params_rv,flag_rv,chi2_rv_old,drv,ntl,1)
!  print *, chi2_rv_old, chi2_tr_old
!  chi2_old_total = chi2_tr_old + chi2_rv_old
!  !Calculate the degrees of freedom
!  nu = dtr + drv - size(params)
!  chi2_red = ( chi2_old_total ) / nu
!  !Print the initial cofiguration
!  print *, ''
!  print *, 'Starting MCMC calculation'
!  print *, 'Initial Chi2_red= ',chi2_red,'nu =',nu
!
!
!  !Let us start the otput file
!  open(unit=101,file='mh_fit.dat',status='unknown')
!
!  !Initialize the values
!  toler_slope = prec 
!  j = 1
!  n = 0
!  get_out = .TRUE.
!
!  !Call a random seed 
!  call init_random_seed()
!
!  !The infinite cycle starts!
!  do while ( get_out )
!    !r will contain random numbers for each variable
!    !Let us add a random shift to each parameter
!    call random_number(r)
!    prec_init = prec * sqrt(chi2_old_total) * r(10)
!    !prec_init = prec 
!    r(1:10+ntl) = ( r(1:10+ntl) - 0.5) * 2.
!    params_new(0:9)      = params(0:9)      + r(1:10)     * prec_init * wtf(0:9)
!    params_new(10:9+ntl) = params(10:9+ntl) + r(11:10+ntl) * prec_init * wtf(10)
!    !Let us calculate our new chi2
!    params_tr_new = params_new(0:8)
!    params_rv_new(0:3) = params_new(0:3)
!    params_rv_new(4:4+ntl) = params_new(9:9+ntl)
!    call find_chi2_tr(xd_tr,yd_tr,errs_tr,params_tr_new,flag_tr,chi2_tr_new,dtr)
!    call find_chi2_rv(xd_rv,yd_rv,errs_rv,tlab,params_rv_new,flag_rv,chi2_rv_new,drv,ntl,1)
!    chi2_new_total = chi2_tr_new + chi2_rv_new
!    !Ratio between the models
!    q = exp( ( chi2_old_total - chi2_new_total ) * 0.5  )
!    !If the new model is better, let us save it
!    if ( q > r(0) ) then
!      chi2_old_total = chi2_new_total
!      params = params_new
!    end if
!    !Calculate the reduced chi square
!    chi2_red = chi2_old_total / nu
!    !Save the data each thin_factor iteration
!    if ( mod(j,thin_factor) == 0 ) then
!      print *, 'c2 tr =', chi2_tr_new, 'c2 rv =', chi2_rv_new
!      print *, 'iter ',j,', Chi2_red =', chi2_red
!      write(101,*) j, chi2_old_total, chi2_red, params
!      !Check convergence here
!      chi2_vec(n) = chi2_red
!      x_vec(n) = n
!      n = n + 1
!      if ( n == size(chi2_vec) ) then
!        call fit_a_line(x_vec,chi2_vec,chi2_y,chi2_slope,nconv)        
!        n = 0
!        !If chi2_red has not changed the last nconv iterations
!        print *, abs(chi2_slope), toler_slope
!        if ( abs(chi2_slope) < toler_slope ) then
!          print *, 'I did my best to converge, chi2_red =', &
!                    chi2_y
!          get_out = .FALSE.
!        end if
!
!      if ( j > maxi ) then
!        print *, 'Maximum number of iteration reached!'
!        get_out = .FALSE.
!      end if
!
!      end if
!      !I checked covergence
!    end if
!    j = j + 1
!  end do
!
!  close(101)
!
!end subroutine
!
!-----------------------------------------------------------
subroutine stretch_move_tr(xd,yd,errs,pars,lims,limits_physical, &
nwalks,maxi,thin_factor,wtf,flag,nconv,datas)
implicit none

!In/Out variables
  integer, intent(in) :: nwalks, maxi, thin_factor, datas, nconv
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  double precision, intent(in), dimension(0:8) :: pars
  double precision, intent(in), dimension(0:(2*9)-1) :: lims
  integer, intent(in), dimension(0:8) :: wtf
  logical, intent(in) :: flag(0:3)
!Local variables
  double precision, dimension(0:8) :: params
  double precision, dimension(0:datas-1) :: zr
  double precision, dimension(0:2*(9)-1) :: limits, limits_physical
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:nwalks-1) :: chi2_old, chi2_new, chi2_red
  double precision, dimension(0:8,0:nwalks-1) :: params_old, params_new
  double precision, dimension(0:8,0:nwalks-1,0:nconv-1) :: params_chains
  double precision  :: q
  double precision  :: esin, ecos, aa, chi2_red_min
  integer :: o, j, n, nu, nk, n_burn, spar, new_thin_factor
  logical :: get_out, is_burn, is_limit_good, is_cvg
  double precision :: r_rand(0:nwalks-1), z_rand(0:nwalks-1)
  integer, dimension(0:nwalks-1) :: r_int
  integer, dimension(0:8) :: wtf_all 
  double precision :: r_real
  character(len=15) :: output_files
  integer :: bc, wc
!external calls
  external :: init_random_seed, find_chi2_tr

  output_files = 'mh_trfit.dat'

  params(:) = pars(:)
  limits(:) = lims(:)
  wtf_all(:) = wtf(:)

  !size of parameters (only parameters to fit!)
  !what parameters are we fitting?
  spar = sum(wtf_all(:))

  !Check if there are flags to evolve modified parameters

  !Period
  if ( flag(0) )  then
    params(1)   = log10(params(1))
    limits(2:3) = log10(limits(2:3))
    limits_physical(2:3) = log10(limits_physical(2:3))
  end if
  ! e and w
  if ( flag(1) ) then
    esin = sqrt(params(2)) * sin(params(3))
    ecos = sqrt(params(2)) * cos(params(3))
    params(2) = esin
    params(3) = ecos
    limits(4) = - sqrt(limits(5))
    limits(5) =   sqrt(limits(5))
    limits(6) = limits(4)
    limits(7) = limits(5)
    limits_physical(4) = - sqrt(limits_physical(5))
    limits_physical(5) =   sqrt(limits_physical(5))
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

  print *, 'CREATING RANDOM SEED'
  call init_random_seed()

  print *, 'CREATING RANDOM (UNIFORM) UNIFORMATIVE PRIORS'
  !Let us create uniformative uniform random priors

  do nk = 0, nwalks - 1

      !Counter for the limits
      j = 0

      do n = 0, 8
        if ( wtf_all(n) == 0 ) then
          !this parameter does not change
          params_old(n,nk) = params(n)
        else
          call random_number(r_real)
          params_old(n,nk) = limits(j+1) - limits(j)
          params_old(n,nk) = limits(j) + r_real*params_old(n,nk) 
          !print *,n, params_old(n,m,nk), limits(j+1,m) - limits(j,m)
          !print *, params_old(2,m,nk)
        end if
          !print *, params_old(2,m,nk), limits(5,m) - limits(4,m)
        j = j + 2 !two limits for each parameter

      end do

      !Check that u1 + u2 < 1 for ew
          is_limit_good = .false.
          do while ( .not. is_limit_good )
            call check_us(params_old(6,nk),params_old(7,nk),is_limit_good)
            if ( .not. is_limit_good  ) then
              params_old(6,nk) = params_old(6,nk) * params_old(6,nk)
              params_old(7,nk) = params_old(7,nk) * params_old(7,nk)
            end if
        end do

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
       

    !Let us find z
    call find_z(xd,params_old(0:5,nk),flag,zr,datas) 
    !Each walker is a point in a parameter space
    !Let us estimate our first chi_2 value for each walker
    call find_chi2_tr(xd,yd,errs,zr,params_old(6:8,nk), &
                      chi2_old(nk),datas)
  end do

  !Calculate the degrees of freedom
  nu = datas - spar
  if ( nu <= 0 ) then
    print *, 'Your number of parameters is larger than your datapoints!'
    print *, 'I cannot fit that!'
    stop
  end if
  chi2_red = chi2_old / nu
  !Print the initial cofiguration
  print *, ''
  print *, 'Starting stretch move MCMC calculation'
  print *, 'Initial Chi2_red= ', sum(chi2_red) / nwalks ,'DOF =', nu

  !Let us start the otput files
  open(unit=123,file=output_files,status='unknown')

  !Initialize the values
  j = 1
  n = 0
  get_out = .TRUE.
  is_burn = .FALSE.
  aa = 2.0 !this is suggested by the original paper
  n_burn = 1
  bc = -1

  open(unit=23,file='worst.dat')
  open(unit=24,file='better.dat')

  !The infinite cycle starts!
  print *, 'STARTING INFINITE LOOP!'
  do while ( get_out )

    !Do the work for all the walkers
    call random_number(z_rand)
    call random_number(r_rand)
    !Create random integers to be the index of the walks
    !Note that r_ink(i) != i (avoid copy the same walker)
    call random_int(r_int,nwalks)

    !call random_number(aa)
    !aa = 1.d0 + 9.9d1 * aa

    do nk = 0, nwalks - 1 !walkers
       ! if (mod(j,thin_factor) == 0) write(23,*) j, nk, r_int(nk)
        params_new(:,nk) = params_old(:,r_int(nk))
    end do

    do nk = 0, nwalks - 1 !walkers

    !Draw the random walker nk, from the complemetary walkers
      z_rand(nk) = z_rand(nk) * aa 

      !The gz function to mantain the affine variance codition in the walks
      call find_gz(z_rand(nk),aa) 

      !Evolve for all the planets for all the parameters
        !params_new(:,nk) = params_old(:,r_int(nk))
      params_new(:,nk) = params_new(:,nk) + wtf_all(:) * z_rand(nk) * &
                           ( params_old(:,nk) - params_new(:,nk) )

        !Let us check the limits
        call check_limits(params_new(:,nk),limits_physical(:), &
                          is_limit_good,9)
        if ( is_limit_good ) then
          ! u1 + u2 < 1 limit
          call check_us(params_new(6,nk),params_new(7,nk),is_limit_good)
          if (is_limit_good ) then
            !Check that e < 1 for ew
            if ( flag(1) ) then
              call check_e(params_new(2,nk),params_new(3,nk),dble(0.9d0),is_limit_good)
            end if          
          end if
        end if


      if ( is_limit_good ) then !evaluate chi2
        !Let us find z
        call find_z(xd,params_new(0:5,nk),flag,zr,datas) 
        !HERE I HAVE TO WRITE THE CHECK ECLIPSE ROUTINE
        !find chi2
        call find_chi2_tr(xd,yd,errs,zr,params_new(6:8,nk),chi2_new(nk),datas)
      else !we do not have a good model
        chi2_new(nk) = huge(dble(0.0)) !a really big number
      end if

      if (  mod(j,thin_factor) == 0  ) then
        if (nk == bc ) then
          print *, 'here'
          write(24,*) j, bc, chi2_old(bc), chi2_new(bc)
          write(*,*) j, bc, chi2_old(bc), chi2_new(bc)
        end if
        if (nk == wc ) then
          write(23,*) j, wc, chi2_old(wc), chi2_new(wc)
        end if
      end if

      !Compare the models 
      q = z_rand(nk)**(spar - 1) * &
          exp( ( chi2_old(nk) - chi2_new(nk) ) * 0.5  )

      if ( q >= r_rand(nk) ) then !is the new model better?
          chi2_old(nk) = chi2_new(nk)
          params_old(:,nk) = params_new(:,nk)
      end if

      chi2_red(nk) = chi2_old(nk) / nu

      !Start to burn-in 
      if ( is_burn ) then
        if ( mod(j,new_thin_factor) == 0 ) then
            write(123,*) j, chi2_old(nk), chi2_red(nk), params_old(:,nk)
        end if
      end if
      !End burn-in

    end do !walkers

   wc = maxloc(chi2_old,dim=1) - 1
   bc = minloc(chi2_old,dim=1) - 1
   !What is our worst chain
   !bc = maxloc(chi2_old,dim=1) - 1
   !print *, bc


     if ( is_burn ) then

        if ( mod(j,new_thin_factor) == 0 ) then
          print *, 'Iter', j
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
          do o = 0, 8 !For all parameters 
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

      !o = minloc(chi2_old,dim=1) - 1
      !nk= maxloc(chi2_old,dim=1) - 1
      !params_old(:,nk) = params_old(:,o)
      !chi2_old(nk) = chi2_old(o)
      !chi2_red(nk) = chi2_red(o)


      end if

            

    end if

  j = j + 1

  end do

  !Close all the files
  close(123)
  close(23)
  close(24)

end subroutine

