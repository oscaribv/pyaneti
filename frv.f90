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
!  $ rv = rv0 - k [ sin ( 2*\pi ( t - t_0) / P ) ] $
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
  double precision, parameter :: pi = 3.1415926535897932384626
!
  rv(:) = rv0 - k * sin( 2.*pi*( t(:) - t0) / P )

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
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:ts-1) :: ta
  double precision :: delta = 1.e-7
  integer :: imax = int(1e5)
!External function
  external :: find_anomaly
!
  !Obtain the eccentric anomaly by using find_anomaly
  call find_anomaly(t,t0,e,w,P,ta,delta,imax,ts)

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
!subroutine rv_curve_mp(t,rv0,t0,k,P,ec,w,rv,ts,np)
!implicit none
!
!!In/Out variables
!  integer, intent(in) :: ts, np
!  double precision, intent(in), dimension(0:ts-1)  :: t
!  double precision, intent(out), dimension(0:ts-1) :: rv
!  double precision, intent(in), dimension(0:np-1) :: k, t0, P, ec, w
!  !Here rv0 depends on the telescope systemic not the planets
!  double precision :: rv0
!!Local variables
!!  double precision, parameter :: pi = 3.1415926535897932384626
!  double precision, dimension(0:ts-1) :: ma, ta
!  double precision :: delta = 1e-5
!  integer :: imax, i
!!External function
!  external :: find_anomaly
!!
!  imax = int(1e5)
!  !Calculate the mean anomaly from the input values
!  rv(:) = rv0
!  do i = 0, np-1
!    ma(:) = mod(2.*pi*(t(:)-t0(i))/P(i),2.*pi)
!!   !Obtain the eccentric anomaly by using find_anomaly
!   call find_anomaly(ma(:),ta(:),ec(i),delta,imax,ts)
!   rv(:) = rv(:) + k(i) * ( cos(ta(:) + w(i) ) + ec(i) * cos(w(i)) )
!  end do
!  
!end subroutine
!
!!-----------------------------------------------------------
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
!subroutine find_chi2_rv(xd,yd,errs,tlab,rv0,k,ec,w,t0,P,chi2,isc,datas,nt,np)
subroutine find_chi2_rv(xd,yd,errs,tlab,params,flag,chi2,isc,datas,nt)
implicit none

!In/Out variables
  integer, intent(in) :: datas, nt
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1) :: tlab
  double precision, intent(in), dimension(0:4+nt) :: params
  logical, intent(in)  :: isc, flag(0:3)
  double precision, intent(out) :: chi2
!Local variables
  double precision :: t0, P, e, w, k
  double precision, dimension(0:nt-1)  :: rv0
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:datas-1) :: model, res
  integer :: i, tel
!External function
  external :: rv_circular, rv_curve

  t0  = params(0)
  P   = params(1)
  e   = params(2)
  w   = params(3)
  k   = params(4)
  rv0 = params(5:4+nt) 

  if ( flag(0) ) P = 10**params(1)
  if ( flag(1) ) then
    e = params(2) * params(2) + params(3) * params(3)
    w = atan2( params(2),params(3) ) 
  end if
  if ( flag(2) ) k = 10**params(4)
  if ( flag(3) ) rv0(:) = 10**params(5:4+nt)

  tel = 0
  !If we want a circular fit, let us do it!
  if ( isc ) then
    do i = 0, datas-1
      if ( tel .ne. tlab(i) ) tel = tel + 1 
      call rv_circular(xd(i),rv0(tel),t0,k,P,model(i),1)
    end do
  !If we want an eccentric fit, let us do it!
  else !if (np == 1) then
    do i = 0, datas-1
      if ( tel .ne. tlab(i)  ) tel = tel + 1 
      call rv_curve(xd(i),rv0(tel),t0,k,P,e,w,model(i),1)
    end do
  !Multiplanet fit
!  else
!    do i = 0, datas-1
!      if ( tel .ne. tlab(i)  ) tel = tel + 1 
!      call rv_curve_mp(xd(i),rv0(tel),t0,k,P,ec,w,model(i),1,np)
!    end do

  end if

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
! ics -> is circular flag, True or False.
! datas, nt -> sizes of xd,yd, errs (datas) and rv0(nt)  
!Output parameter:
! This functions outputs a file called mh_rvfit.dat
!-----------------------------------------------------------
subroutine metropolis_hastings_rv(xd,yd,errs,tlab,params,prec,maxi,thin_factor,ics,wtf,flag,nconv,datas,nt)
implicit none

!In/Out variables
  integer, intent(in) :: maxi, thin_factor, datas, nt, nconv
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1)  :: tlab
  double precision, intent(inout), dimension(0:4+nt) :: params
  !f2py intent(in,out)  :: params
  double precision, intent(in), dimension(0:5) :: wtf
  double precision, intent(in)  :: prec
  logical, intent(in) :: ics, flag(0:3)
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
    esin = sqrt(params(2)) * dsin(params(3))
    ecos = sqrt(params(2)) * dcos(params(3))
    params(2) = esin
    params(3) = ecos
  end if
  if ( flag(2) ) params(4) = log10(params(4))
  if ( flag(3) ) params(5:4+nt) = log10(params(5:4+nt))

  !Let us estimate our fist chi_2 value
  call find_chi2_rv(xd,yd,errs,tlab,params,flag,chi2_old,ics,datas,nt)
  !Calculate the degrees of freedom
  nu = datas - size(params)
  !If we are fixing a circular orbit, ec and w are not used 
  if ( ics ) nu = nu + 2 
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
    call find_chi2_rv(xd,yd,errs,tlab,params_new,flag,chi2_new,ics,datas,nt)
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
subroutine stretch_move_rv(xd,yd,errs,tlab,params,nwalks,prec,maxi,thin_factor,ics,wtf,flag,nconv,datas,nt)
implicit none

!In/Out variables
  integer, intent(in) :: nwalks, maxi, thin_factor, datas, nt, nconv
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1)  :: tlab
  double precision, intent(inout), dimension(0:4+nt) :: params
  !f2py intent(in,out)  :: params
  integer, intent(in), dimension(0:5) :: wtf
  double precision, intent(in)  :: prec
  logical, intent(in) :: ics, flag(0:3)
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:nwalks-1) :: chi2_old, chi2_new, chi2_red
  double precision, dimension(0:4+nt,0:nwalks-1) :: params_old, params_new
  double precision, dimension(0:nconv-1) :: chi2_vec, x_vec
  double precision  :: q, chi2_y, chi2_slope, toler_slope
  double precision, dimension(0:4+nt) :: wtf_all
  double precision  :: prec_init, esin, ecos, nu, aa, chi2_red_min
  integer :: j, n, nk, n_burn, spar, chi2_walker, new_thin_factor,good_chain(0:nwalks-1)
  logical :: get_out, is_burn
  !Let us add a plus random generator
  double precision, dimension(0:nwalks-1) :: r_rand, z_rand, r_real
  integer, dimension(0:nwalks-1) :: r_int
!external calls
  external :: init_random_seed, find_chi2_rv

  spar = size(params)

  if ( flag(0) ) params(1) = log10(params(1))
  if ( flag(1) ) then
    esin = sqrt(params(2)) * dsin(params(3))
    ecos = sqrt(params(2)) * dcos(params(3))
    params(2) = esin
    params(3) = ecos
  end if
  if ( flag(2) ) params(4) = log10(params(4))
  if ( flag(3) ) params(5:4+nt) = log10(params(5:4+nt))

  wtf_all(0:4) = wtf(0:4)
  wtf_all(5:4+nt) = wtf(5)

  !Call a random seed 
  print *, 'CREATING RANDOM SEED'
  call init_random_seed()

  call random_number(r_real)
  !Let us initinate all the walkers with the priors
  !I am not sure it this works, it is just a test
  r_real = ( r_real - 0.5 ) * 2.0
  do j = 0, nwalks - 1
    params_old(:,j) = params * ( 1 + wtf_all*r_real(j)*0.01 ) 
    !Let us estimate our first chi_2 value
    call find_chi2_rv(xd,yd,errs,tlab,params_old(:,j),flag,chi2_old(j),ics,datas,nt)
  end do
  

  !Calculate the degrees of freedom
  nu = dble(datas - spar)
  !If we are fixing a circular orbit, ec and w are not used 
  if ( ics ) nu = nu + 2 
  chi2_red = chi2_old / nu
  !Print the initial cofiguration
  print *, ''
  print *, 'Starting stretch move MCMC calculation'
  print *, 'Initial Chi2_red= ', chi2_red(0),'nu =', nu

  !Let us start the otput file
  open(unit=101,file='mh_rvfit.dat',status='unknown')
  !Initialize the values

  toler_slope = prec
  j = 1
  n = 0
  get_out = .TRUE.
  is_burn = .FALSE.

  aa = 2.0
  n_burn = 1

  !The infinite cycle starts!
  print *, 'STARTING INFINITE LOOP!'

  do while ( get_out )

    !Do the work for all the walkers

    !--------------------------------------------
    !AVOID THAT r_int(nk) is equal to nk
    !--------------------------------------------

    call random_number(z_rand)
    call random_number(r_rand)
    call random_number(r_real)
    !We have to add one and then extract one
    !In this way r_int can be also zero
    r_int = int( r_real * (nwalks + 1 ) )
    r_int = r_int - 1

    do nk = 0, nwalks - 1
    !Draw the random walker nk, from the complemetary walkers
    !This definition does not avoid to copy the same k walker
      params_new(:,nk) = params_old(:,r_int(nk))
      !Let us generate a random step
      z_rand(nk) = z_rand(nk) * aa ! a = 2 as suggested by emcee paper
      call find_gz(z_rand(nk),aa) 
    
      !Now we can have the evolved walker
      params_new(:,nk) = params_new(:,nk) + wtf_all * z_rand(nk) * &
                       ( params_old(:,nk) - params_new(:,nk) )

      !Obtain the new chi square 
      call find_chi2_rv(xd,yd,errs,tlab,params_new(:,nk),flag,chi2_new(nk),ics,datas,nt)

      !Is the new model better? 
      q = z_rand(nk)**( spar - 1.) * &
          exp( ( chi2_old(nk) - chi2_new(nk) ) * 0.5  )

      if ( q >= r_rand(nk) ) then
        chi2_old(nk) = chi2_new(nk)
        params_old(:,nk) = params_new(:,nk)
      end if

      chi2_red(nk) = chi2_old(nk) / nu

      !Start to burn-in 
      if ( is_burn .and. mod(j,new_thin_factor) == 0 .and. nk == good_chain(0) ) then
          write(101,*) n_burn, chi2_old(nk), chi2_red(nk), params_old(:,nk)
          print *, n_burn, chi2_red(nk)
      end if
      !End burn-in

    end do

     if ( is_burn .and. mod(j,new_thin_factor) == 0 ) n_burn = n_burn + 1
     if ( n_burn > nconv ) get_out = .false.

    !chi2_red_min = minval(chi2_red)   

    !Save the data each thin_factor iteration
    if ( mod(j,thin_factor) == 0 .and. .not. is_burn ) then

      good_chain = minloc(chi2_red)
      chi2_red_min = minval(chi2_red)
      print *, 'iter ',j,', Chi2_red =', chi2_red_min, good_chain(0)
      !Check convergence here
      chi2_vec(n) = chi2_red_min
      x_vec(n) = n
      n = n + 1

      if ( n == size(chi2_vec) ) then

        call fit_a_line(x_vec,chi2_vec,chi2_y,chi2_slope,nconv)
        n = 0
        !If chi2_red has not changed the last nconv iterations
        print *, abs(chi2_slope), toler_slope / (chi2_y)
        if ( abs(chi2_slope) < ( toler_slope / (chi2_y) ) ) then
          print *, 'THE CHAIN HAS CONVERGED'
          print *, 'Starting burning-in phase'
          is_burn = .True.
          new_thin_factor = 10
          print *, 'The best chain is', good_chain(0), &
          'with chi2_red =', chi2_red(good_chain(0))
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

