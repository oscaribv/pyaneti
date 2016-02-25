!------------------------------------------------------------
!                         frv.f90
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
subroutine find_z(t,t0,e,w,P,i,a,z,ts)
implicit none

!In/Out variables
  integer, intent(in) :: ts
  double precision, intent(in), dimension(0:ts-1)  :: t
  double precision, intent(in) :: t0, e, w, P, i, a
  double precision, intent(out), dimension(0:ts-1) :: z
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:ts-1) :: ma, ta
  double precision :: delta = 1e-4
  double precision :: si
  integer :: imax
!External function
  external :: find_anomaly
!

  si = sin(i)

  imax = int(1e7)
  !Calculate the mean anomaly from the input values
  ma(:) = 2.* pi * ( t(:) - t0 ) / P
  !Obtain the eccentric anomaly by using find_anomaly
  call find_anomaly(ma,ta,e,delta,imax,ts)

  z(:) = a * ( 1- e*e ) / (1 + e * cos(ta(:))) * &
         sqrt( 1. - sin(w+ta(:))*sin(w+ta(:))*si*si) 

  !z has been calculated
  
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
subroutine find_chi2_tr(xd,yd,errs,t0,e,w,P,i,a,u1,u2,pz,chi2,isc,datas)
implicit none

!In/Out variables
  integer, intent(in) :: datas
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  double precision, intent(in) :: t0, e, w, P, i, a, u1, u2, pz
  logical, intent(in)  :: isc
  double precision, intent(out) :: chi2
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision, dimension(0:datas-1) :: z, res, muld, mu
!External function
  external :: occultquad, find_z

  !Let us find the projected distance z
  call find_z(xd,t0,e,w,P,i,a,z,datas)
  !Now we have z, let us use Agol's routines
  !If we want a circular fit, let us do it!
  if ( isc ) then
    print *, "This option is obsolete now"
    stop
  else
    call occultquad(z,u1,u2,pz,muld,mu,datas)
  end if

  !Let us calculate the residuals
  ! chi^2 = \Sum_i ( M - O )^2 / \sigma^2
  !Here I am assuming that we want linmd darkening
  !If this is not true, use mu 
  res(:) = muld(:) - yd(:)  
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
subroutine metropolis_hastings_tr(xd,yd,errs,t0,e,w,P,i,a,u1,u2,pz,prec,maxi,thin_factor,chi2_toler,ics,datas)
implicit none

!In/Out variables
  integer, intent(in) :: maxi, thin_factor, datas
  double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
  double precision, intent(inout)  :: t0,e,w,P,i,a,u1,u2,pz,prec,chi2_toler
  !f2py intent(in,out)  :: t0,e,w,P,i,a,u1,u2,pz,prec,chi2_toler
  logical, intent(in) :: ics
!Local variables
  double precision, parameter :: pi = 3.1415926535897932384626
  double precision :: chi2_old, chi2_new, chi2_red
  double precision :: t0n,en,wn,Pn,ni,an,u1n,u2n,pzn
  double precision :: st0,se,sw,sP,si,sa,su1,su2,spz
  double precision  :: q
  integer :: j, n, nu
  real, dimension(0:9) :: r
!external calls
  external :: init_random_seed, find_chi2_tr
 
  !Calculate the step size based in the actual value of the
  !parameters and the prec variable
  st0 = t0*prec
  se  = e*prec
  sw  = w*prec
  sP  = P*prec
  si  = i*prec
  sa  = a*prec
  su1 = u1*prec
  su2 = u2*prec
  spz = pz*prec

  !Let us estimate our fist chi_2 value
  call find_chi2_tr(xd,yd,errs,t0,e,w,P,i,a,u1,u2,pz,chi2_old,ics,datas)
  !Calculate the degrees of freedom
  nu = datas - 9
  !Print the initial cofiguration
  print *, ''
  print *, 'Starting MCMC calculation'
  print *, 'Initial Chi_2: ', chi2_old,'nu =', nu
  chi2_red = chi2_old / nu
  !Call a random seed 
  call init_random_seed()

  !Let us start the otput file
  open(unit=101,file='mh_trfit.dat',status='unknown')
  write(101,*)'# i chi2 chi2_red k ec w t0 P rv0mc(vector)'
  write(101,*) 0,chi2_old,chi2_red
  !Initialize the values
  n = 1
  j = 1

  !The infinite cycle starts!
  do while ( chi2_red >= 1 + chi2_toler .and. j <= maxi - 1 )
    !r will contain random numbers for each variable
    !Let us add a random shift to each parameter
    call random_number(r)
    r(1:9) = ( r(1:9) - 0.5) * 2.
     t0n = t0 + r(1) * st0
     en  = e  + r(2) * se
     wn  = w  + r(3) * sw
     Pn  = P  + r(4) * sP
     ni  = i  + r(5) * si
     an  = a  + r(6) * sa
     u1n = u1 !+ r(7) * su1
     u2n = u2 !+ r(8) * su2
     pzn = pz + r(9) * spz
      !print *, u1n, u2n
    !Let us calculate our new chi2
    call find_chi2_tr(xd,yd,errs,t0n,en,wn,Pn,ni,an,u1n,u2n,pzn,chi2_new,ics,datas)
    !Ratio between the models
    q = exp( ( chi2_old - chi2_new ) * 0.5  )
    !If the new model is better, let us save it
    if ( q > r(0) ) then
      chi2_old = chi2_new
      t0 = t0n
      e  = en
      w  = wn
      P  = Pn
      i  = ni
      a  = an
!      u1 = u1n
!      u2 = u2n
      pz = pzn
    end if
    !Calculate the reduced chi square
    chi2_red = chi2_old / nu
    !Save the data each thin_factor iteration
    if ( mod(j,thin_factor) == 0 ) then
      print *, 'iter ',j,'  of ',maxi
      print *, 'chi2 = ',chi2_old,'reduced chi2 =', chi2_red
      write(101,*) j,chi2_old,chi2_red,t0,e,w,P,i,a,u1,u2,pz
      n = n + 1
    end if
    j = j + 1
  end do

  print *, 'Final chi2 = ',chi2_old,'. Final reduced chi2 =', chi2_red

  close(101)

end subroutine

