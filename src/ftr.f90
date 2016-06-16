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
  double precision, dimension(0:ts-1) :: ta, swt
  double precision :: t0, P, e, w, i, a
  double precision :: si
!External function
  external :: find_anomaly
!

  t0  = pars(0)
  P   = pars(1)
  e   = pars(2)
  w   = pars(3)
  i   = pars(4)
  a   = pars(5)

  if ( flag(0) ) P = 1.d0**pars(1)
  if ( flag(1) ) then
    e = pars(2) * pars(2) + pars(3) * pars(3)
    w = atan2(pars(2),pars(3))
  end if
  if (flag(2)) i = asin(pars(4))
  if (flag(3)) a = 1.d0**pars(5)

  !Obtain the eccentric anomaly by using find_anomaly
  call find_anomaly(t,t0,e,w,P,ta,ts)
  swt = sin(w+ta)

  si = sin(i)
  z = a * ( 1.d0 - e * e ) * sqrt( 1.d0 - swt * swt * si * si ) &
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
  double precision :: limit

  !This works only for not grazing transits
  limit = 1.d0 + pz
  is_good = .false.
  !At least we have to have one eclipse condition
  do i = 1, sizez - 2
    if ( z(i) < limit ) then
      is_good = .true.
      exit
    end if
  end do

  !if ( z(0) < limit .and. z(sizez-1) < limit ) is_good = .false.

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
subroutine find_chi2_tr(xd,yd,errs,pars,flag,params,n_cad,t_cad,chi2,datas)
!subroutine find_chi2_tr(yd,errs,z,params,chi2,datas)
implicit none

!In/Out variables
  integer, intent(in) :: datas, n_cad
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  double precision, intent(in), dimension(0:5) :: pars
  double precision, intent(in) :: t_cad
  logical, intent(in), dimension(0:3) :: flag 
  double precision, intent(in), dimension (0:2) :: params
  double precision, intent(out) :: chi2
!Local variables
  double precision, dimension(0:datas-1) :: res, muld, mu
  double precision :: u1, u2, pz, q1k, q2k, zdum(0:0)
  !double precision, dimension(0:datas-1,0:n_cad-1)  :: xd_ub, z, flux_ub
  double precision, dimension(0:n_cad-1)  :: xd_ub, z, flux_ub
  integer ::  j, k
!External function
  external :: occultquad

  q1k = params(0)
  q2k = params(1)
  pz  = params(2) 
  !re-transform the parameters to u1 and u2
  u1 = sqrt(q1k)
  u2 = u1*( 1.d0 - 2.d0*q2k)
  u1 = 2.d0*u1*q2k

  !Selective re-sampling
  do j = 0, datas - 1

    !Are we generating an eclipse?
    call find_z(xd(j),pars,flag,zdum,1) 

    if ( zdum(0) > 1.0 + pz ) then

      muld(j) = 1.d0 !This is not eclipse

    else   

      do k = 0, n_cad - 1
        xd_ub(k) = xd(j) + t_cad*((k+1)-0.5*(n_cad+1))/n_cad 
      end do

      call find_z(xd_ub,pars,flag,z,n_cad) 
      !Now we have z, let us use Agol's routines
      call occultquad(z,u1,u2,pz,flux_ub,mu,n_cad)

      !Re-bin the data
      muld(j) = sum(flux_ub) / n_cad

    end if

  end do

    !Let us calculate the residuals
  ! chi^2 = \Sum_i ( M - O )^2 / \sigma^2
  !Here I am assuming that we want limb darkening
  !If this is not true, use mu 
  res(:) = ( muld(:) - yd(:) ) / errs(:) 
  chi2 = dot_product(res,res)

end subroutine
!-----------------------------------------------------------
subroutine stretch_move_tr(xd,yd,errs,pars,lims,limits_physical, &
nwalks,a_factor,maxi,thin_factor,n_cad,t_cad,wtf,flag,nconv,datas)
implicit none

!In/Out variables
  integer, intent(in) :: nwalks, maxi, thin_factor, n_cad, datas, nconv
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  double precision, intent(in), dimension(0:8) :: pars
  double precision, intent(in), dimension(0:(2*9)-1) :: lims
  double precision, intent(in) :: a_factor, t_cad
  integer, intent(in), dimension(0:8) :: wtf
  logical, intent(in) :: flag(0:3)
!Local variables
  double precision, dimension(0:8) :: params
  double precision, dimension(0:datas-1) :: zr
  double precision, dimension(0:2*(9)-1) :: limits, limits_physical
  double precision, dimension(0:nwalks-1) :: chi2_old, chi2_new, chi2_red
  double precision, dimension(0:8,0:nwalks-1) :: params_old, params_new
  double precision, dimension(0:8,0:nwalks-1,0:nconv-1) :: params_chains
  double precision  :: q, q1k, q2k
  double precision  :: esin, ecos, aa, chi2_red_min
  integer :: o, j, n, nu, nk, n_burn, spar, new_thin_factor
  logical :: get_out, is_burn, is_limit_good, is_cvg, is_eclipse
  double precision :: r_rand(0:nwalks-1), z_rand(0:nwalks-1)
  integer, dimension(0:nwalks-1) :: r_int
  integer, dimension(0:8) :: wtf_all 
  double precision :: r_real
  character(len=15) :: output_files
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
  !q1 and q2
  limits(12) = 0.0
  limits_physical(12) = 0.0
  limits(13) = 1.0
  limits_physical(13) = 1.0
  limits(14) = 0.0
  limits_physical(14) = 0.0
  limits(15) = 1.0
  limits_physical(15) = 1.0

  print *, 'CREATING RANDOM SEED'
  call init_random_seed()

  print *, 'CREATING RANDOM (UNIFORM) UNIFORMATIVE PRIORS'
  !Let us create uniformative uniform random priors

  nk = 0
  do while (nk < nwalks )

      print *, 'creating walker ', nk + 1
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
!          is_limit_good = .false.
!          do while ( .not. is_limit_good )
!            call check_us(params_old(6,nk),params_old(7,nk),is_limit_good)
!            if ( .not. is_limit_good  ) then
!              params_old(6,nk) = params_old(6,nk) * params_old(6,nk)
!              params_old(7,nk) = params_old(7,nk) * params_old(7,nk)
!            end if
!        end do

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

    !Let us check if there is an eclipse
    call check_eclipse(zr,params_old(8,nk),is_eclipse,datas)
    !if we do not have an eclipse, let us recalculate this wlaker
    if ( .not. is_eclipse ) then 
      nk = nk
    else
    !Each walker is a point in a parameter space
    !Let us estimate our first chi_2 value for each walker
    call find_chi2_tr(xd,yd,errs,params_old(0:5,nk),flag,params_old(6:8,nk), &
                      n_cad, t_cad, chi2_old(nk),datas)
      nk = nk + 1
    end if

  end do
!  stop
  

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

    if ( a_factor <= 1.d0 ) then
      call random_number(aa)
      aa = 1.d0 + abs(aa*a_factor)
    end if

    do nk = 0, nwalks - 1 !walkers
        params_new(:,nk) = params_old(:,r_int(nk))
    end do

    do nk = 0, nwalks - 1 !walkers

    !Draw the random walker nk, from the complemetary walkers
      z_rand(nk) = z_rand(nk) * aa 

      !The gz function to mantain the affine variance codition in the walks
      call find_gz(z_rand(nk),aa) 

      !Evolve for all the planets for all the parameters
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
        call check_eclipse(zr,params_new(8,nk),is_eclipse,datas)
        !find chi2
        if ( is_eclipse ) then
        call find_chi2_tr(xd,yd,errs,params_new(0:5,nk),flag,params_new(6:8,nk),n_cad,t_cad,chi2_new(nk),datas)
        else
          chi2_new(nk) = huge(0.0d0)
        end if
      else !we do not have a good model
        chi2_new(nk) = huge(dble(0.0)) !a really big number
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
            new_thin_factor = thin_factor
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

  !Close all the files
  close(123)

end subroutine

