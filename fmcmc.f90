!Template to write a fortran90 program
!written by Oscar Barragán, Jan 13 2014

!-------------------------------------------------------------------------
!Name of the program
!Here write a simple description of the program
!Date -->
!-------------------------------------------------------------------------

subroutine find_anomaly(man,ta,ec,delta,imax,dman)
implicit none

integer, intent(in) :: dman
double precision, intent(in) , dimension(0:dman-1) :: man
double precision, intent(out), dimension(0:dman-1) :: ta
double precision, intent(in) :: ec, delta
integer, intent(in) :: imax

integer :: i,j,n
double precision, dimension(0:dman-1) :: f, df

  ta(:)  = 0.0
  f(:)   = ta(:) - ec * sin(ta(:)) - man(:)
  df(:)  =   1.0 - ec * cos(ta(:))
  n = 0

  do i = 0, dman-1
    do while ( abs(f(i)) >= delta .and. n <= imax )
      ta(i)  = ta(i) - f(i) / df(i)
      f(i)   = ta(i) - ec * sin(ta(i)) - man(i)
      df(i)  =   1.0 - ec * cos(ta(i))
      n = n + 1
    end do
  end do

  if ( n > imax ) then
    print *, 'I am tired, so much N-R'
    stop
  end if 

  ta(:) = sqrt( (1. + ec) / (1. -ec) ) * tan(ta(:)*0.5)
  ta(:) = 2. * atan(ta(:))

end subroutine

!-----------------------------------------------------------

subroutine rv_circular(t,rv0,t0,k,P,rv,ts)
implicit none

integer, intent(in) :: ts
double precision, intent(in), dimension(0:ts-1)  :: t
double precision, intent(out), dimension(0:ts-1) :: rv
double precision, intent(in) :: k, rv0, t0, P

double precision, parameter :: pi = 3.1415926535897932384626

  rv(:) = rv0 - k * sin( 2.*pi*( t(:) - t0) / P )

end subroutine

!-----------------------------------------------------------

subroutine rv_curve(t,rv0,t0,k,P,ec,w,rv,ts)
implicit none

integer, intent(in) :: ts
double precision, intent(in), dimension(0:ts-1)  :: t
double precision, intent(out), dimension(0:ts-1) :: rv
double precision, intent(in) :: k, rv0, t0, P, ec, w

double precision, parameter :: pi = 3.1415926535897932384626
double precision, dimension(0:ts-1) :: ma, ta
double precision :: delta = 1e-4
integer :: imax = 1e5

external :: find_anomaly

  ma(:) = 2.* pi * ( t(:) - t0 ) / P
  call find_anomaly(ma,ta,ec,delta,imax,ts)

  rv(:) = rv0 + k * ( cos(ta(:) + w ) + ec * cos(w) )
  
end subroutine

!-----------------------------------------------------------

subroutine find_chi2(xd,yd,errs,rv0,k,ec,w,t0,P,chi2,datas)
implicit none

integer, intent(in) :: datas
double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
double precision, intent(in)  :: k, rv0, t0, P, ec, w
double precision, intent(out) :: chi2

double precision, parameter :: pi = 3.1415926535897932384626
double precision, dimension(0:datas-1) :: model, res

external :: rv_circular, rv_curve

  chi2 = 0.0

  call rv_circular(xd,rv0,t0,k,P,model,datas)

  res(:) = ( yd(:) - model(:) ) / errs(:)
  res(:) = res(:) * res(:)

  chi2 = sum(res)

end subroutine

!-----------------------------------------------------------

subroutine mcmc_rv(xd,yd,errs,rv0,k,ec,w,t0,P,rv0mc,kmc,ecmc,wmc,t0mc,Pmc,imcmc,datas)
implicit none

integer, intent(in) :: imcmc, datas
double precision, intent(in), dimension(0:datas-1)  :: xd, yd, errs
double precision, intent(in)  :: k, rv0, t0, P, ec, w
double precision, intent(out), dimension(0:imcmc-1)  :: kmc, rv0mc, t0mc, Pmc, ecmc, wmc

double precision, parameter :: pi = 3.1415926535897932384626
double precision :: chi2_old, chi2_new, prec = 1e-4
double precision  :: sk, srv0, st0, sP, sec, sw
double precision  :: q
integer :: i
real, dimension(0:imcmc-1,0:2) :: r

external :: find_chi2

  kmc(0)    = k
  rv0mc(0)  = rv0
  t0mc(0)   = t0
  Pmc(0)    = P
  ecmc(0)   = ec
  wmc(0)    = w

  sk = k * prec
  srv0 = rv0 * prec
  st0 = st0 * prec
  sP = sP * prec
  sec = 1. * prec
  sw = 2*pi * prec

  call find_chi2(xd,yd,errs,rv0mc(0),kmc(0),ecmc(0),wmc(0),t0mc(0),Pmc(0),chi2_old,datas)

  call init_random_seed()
  !Let us create an array of random numbers to save time
  !Check the ram consumption
  call random_number(r)
  r = (r - 0.5) * 2.

  do i = 1, imcmc - 1
    rv0mc(i) = rv0mc(i-1) + r(i,0) * srv0
    kmc(i)   =   kmc(i-1) + r(i,1) * sk
    call find_chi2(xd,yd,errs,rv0mc(i),kmc(i),ecmc(0),wmc(0),t0mc(0),Pmc(0),chi2_new,datas)
    q = exp( ( chi2_old - chi2_new ) * 0.5  )
    if ( q < abs(r(i,2)) ) then
      rv0mc(i) = rv0mc(i-1)
        kmc(i) = kmc(i-1)
    else
      chi2_old = chi2_new
    end if
  end do

end subroutine
