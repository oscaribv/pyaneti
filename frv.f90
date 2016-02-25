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
subroutine rv_curve(t,rv0,t0,k,P,ec,w,rv,ts)
implicit none

!In/Out variables
  integer, intent(in) :: ts
  double precision, intent(in), dimension(0:ts-1)  :: t
  double precision, intent(out), dimension(0:ts-1) :: rv
  double precision, intent(in) :: k, rv0, t0, P, ec, w
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:ts-1) :: ma, ta
  double precision :: delta = 1e-4
  integer :: imax
!External function
  external :: find_anomaly
!

  imax = int(1e5)
  !Calculate the mean anomaly from the input values
  ma(:) = 2.* pi * ( t(:) - t0 ) / P
  !Obtain the eccentric anomaly by using find_anomaly
  call find_anomaly(ma,ta,ec,delta,imax,ts)

  rv(:) = rv0 + k * ( cos(ta(:) + w ) + ec * cos(w) )
  
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
subroutine rv_curve_mp(t,rv0,t0,k,P,ec,w,rv,ts,np)
implicit none

!In/Out variables
  integer, intent(in) :: ts, np
  double precision, intent(in), dimension(0:ts-1)  :: t
  double precision, intent(out), dimension(0:ts-1) :: rv
  double precision, intent(in), dimension(0:np-1) :: k, t0, P, ec, w
  !Here rv0 depends on the telescope systemic not the planets
  double precision :: rv0
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:ts-1) :: ma, ta
  double precision :: delta = 1e-5
  integer :: imax, i
!External function
  external :: find_anomaly
!
  imax = int(1e5)
  !Calculate the mean anomaly from the input values
  rv(:) = rv0
  do i = 0, np-1
    ma(:) = 2.* pi * ( t(:) - t0(i) ) / P(i)
   !Obtain the eccentric anomaly by using find_anomaly
   call find_anomaly(ma(:),ta(:),ec(i),delta,imax,ts)
   rv(:) = rv(:) + k(i) * ( cos(ta(:) + w(i) ) + ec(i) * cos(w(i)) )
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
subroutine find_chi2_rv(xd,yd,errs,tlab,rv0,k,ec,w,t0,P,chi2,isc,datas,nt,np)
implicit none

!In/Out variables
  integer, intent(in) :: datas, nt, np
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1)  :: tlab
  double precision, intent(in), dimension(0:nt-1)  :: rv0
  double precision, intent(in), dimension(0:np-1)  :: k, t0, P, ec, w
  logical, intent(in)  :: isc
  double precision, intent(out) :: chi2
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:datas-1) :: model, res
  integer :: i, j, tel
!External function
  external :: rv_circular, rv_curve

  chi2 = 0.0
  tel = 0
  i = 0
  !If we want a circular fit, let us do it!
  if ( isc ) then
    do i = 0, datas-1
      if ( tel .ne. tlab(i) ) tel = tel + 1 
      call rv_circular(xd(i),rv0(tel),t0(0),k(0),P(0),model(i),1)
    end do
  !If we want an eccentric fit, let us do it!
  else if (np == 1) then
    do i = 0, datas-1
      if ( tel .ne. tlab(i)  ) tel = tel + 1 
      call rv_curve(xd(i),rv0(tel),t0(0),k(0),P(0),ec(0),w(0),model(i),1)
    end do
  !Multiplanet fit
  else
    do i = 0, datas-1
      if ( tel .ne. tlab(i)  ) tel = tel + 1 
      call rv_curve_mp(xd(i),rv0(tel),t0,k,P,ec,w,model(i),1,np)
    end do

  end if

  !Let us calculate the residuals
  ! chi^2 = \Sum_i ( M - O )^2 / \sigma^2
  res(:) = model(:) - yd(:)  
  res(:) = res(:) * res(:) / ( errs(:) * errs(:) )
  chi2 = sum(res)

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
!I could deal with differences in parameters by ussing vectors insead of
!independient float parameters, PLEASE OSCAR DO THIS!
subroutine metropolis_hastings_rv(xd,yd,errs,tlab,rv0mc,kmc,ecmc,wmc,t0mc,Pmc,prec,maxi,thin_factor,chi2_toler,ics,datas,nt,np)
implicit none

!In/Out variables
  integer, intent(in) :: maxi, thin_factor, datas, nt, np
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  integer, intent(in), dimension(0:datas-1)  :: tlab
  double precision, intent(inout), dimension(0:nt-1)  :: rv0mc
  double precision, intent(in)  :: prec, chi2_toler
  double precision, intent(inout), dimension(0:np-1)  :: kmc,t0mc, Pmc, ecmc, wmc
  !f2py intent(in,out)  ::rv0mc, kmc,t0mc, Pmc, ecmc, wmc
  logical, intent(in) :: ics
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision :: chi2_old, chi2_new, chi2_red
  double precision, dimension(0:nt-1)  :: rv0mcn
  double precision, dimension(0:np-1) :: kmcn, t0mcn, Pmcn, ecmcn, wmcn
  double precision  :: sk, srv0, st0, sP, sec, sw
  double precision  :: q
  integer :: i, j, nu, yo
  real, dimension(0:5+nt) :: r
!external calls
  external :: init_random_seed, find_chi2_rv

  !Calculate the step size based in the actual value of the
  !parameters and the prec variable
  sk   = kmc(0) * prec
  srv0 = kmc(0) * prec
  st0  = t0mc(0)* prec
  sP   = Pmc(0) * prec
  sec  = 1.     * prec
  sw   = 2.*pi  * prec

  !Let us estimate our fist chi_2 value
  call find_chi2_rv(xd,yd,errs,tlab,rv0mc,kmc,ecmc,wmc,t0mc,Pmc,chi2_old,ics,datas,nt,np)
  !Calculate the degrees of freedom
  nu = datas - 5 - nt
  !If we are fixing a circular orbit, ec and w are not used 
  if ( ics ) nu = nu + 2 
  !Print the initial cofiguration
  print *, ''
  print *, 'Starting MCMC calculation'
  print *, 'Initial Chi_2: ', chi2_old,'nu =', nu
  chi2_red = chi2_old / nu
  !Call a random seed 
  call init_random_seed()

  !Let us start the otput file
  open(unit=101,file='mh_rvfit.dat',status='unknown')
  write(101,*)'# i chi2 chi2_red k ec w t0 P rv0mc(vector)'
  write(101,*) 0,chi2_old,chi2_red,kmc, ecmc, wmc, t0mc, Pmc,rv0mc
  !Initialize the values
  i = 1
  yo = 0

  !The infinite cycle starts!
  do while ( chi2_red >= 1. + chi2_toler .and. i <= maxi )
    !r will contain random numbers for each variable
    !Let us add a random shift to each parameter
    call random_number(r)
    r(1:5+nt) = ( r(1:5+nt) - 0.5) * 2.
    kmcn(yo)   =   kmc(yo)   + r(1) * sk
    ecmcn(yo)  =   ecmc(yo)  + r(2) * sec
    ecmcn(yo)  =   abs(ecmcn(yo))  
    wmcn(yo)   =   wmc(yo)   + r(3) * sw
    t0mcn(yo)  =   t0mc(yo)  + r(4) * st0
    Pmcn(yo)   =   Pmc(yo)   + r(5) * sP
    do j = 0, nt-1
      rv0mcn(j) =   rv0mc(j) + r(6+j) * srv0
    end do
    !Let us calculate our new chi2
    call find_chi2_rv(xd,yd,errs,tlab,rv0mcn,kmcn,ecmcn,wmcn,t0mcn,Pmcn,chi2_new,ics,datas,nt,np)
    !Ratio between the models
    q = exp( ( chi2_old - chi2_new ) * 0.5  )
    !If the new model is better, let us save it
    if ( q > r(0) ) then
      chi2_old = chi2_new
       rv0mc(yo) = rv0mcn(yo)
         kmc(yo) = kmcn(yo)
        ecmc(yo) = ecmcn(yo)
         wmc(yo) = wmcn(yo)
        t0mc(yo) = t0mcn(yo)
         Pmc(yo) = Pmcn(yo)
    end if
    yo = mod(i,np)
    !Calculate the reduced chi square
    chi2_red = chi2_old / nu
    !Save the data each thin_factor iteration
    if ( mod(i,thin_factor) == 0 ) then
      print *, 'iter ',i,'  of ',maxi
      print *, 'chi2 = ',chi2_old,'reduced chi2 =', chi2_red
      write(101,*) i,chi2_old,chi2_red,kmc, ecmc, wmc, t0mc, Pmc, rv0mc
    end if
    i = i + 1
  end do

  print *, 'Final chi2 = ',chi2_old,'. Final reduced chi2 =', chi2_red

  close(101)

end subroutine

!-------------------------------------------------------

