!
subroutine get_QP_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(in) ::  pars(0:npars-1)
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: titj, tau, sintau, alpha, beta
  real(kind=mireal) :: le, lp, P, le2, lp2, P2

  le = pars(0)
  lp = pars(1)
  P  = pars(2)

  le2 = le*le
  lp2 = lp*lp
  P2  = P*P

  call fcdist(x1,x2,titj,nx1,nx2)

  tau = 2.*pi*titj/P

  sintau = sin(tau)

  alpha = pi*sintau/(2.*P*lp2)

  beta = titj/le2

  gamma_g_g = - (sin(tau/2.))**2/2./lp2 - titj*titj/2./le2
  gamma_g_g = exp(gamma_g_g)

  gamma_g_dg = - gamma_g_g * (alpha + beta)

  gamma_dg_g = - gamma_g_dg

  gamma_dg_dg = - alpha*alpha &
                + pi*pi*cos(tau)/P2/lp2         &
                - tau*sintau/(2.*lp2*le2)     &
                - beta*beta                     &
                + 1./le2
  gamma_dg_dg = gamma_g_g * gamma_dg_dg

end subroutine get_QP_gammas

subroutine get_EXP_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(in) ::  pars(0:npars-1)
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: titj
  real(kind=mireal) :: l !Lambda parameter for the square exponential kernel

  l = pars(0)

  call fcdist(x1,x2,titj,nx1,nx2)

  gamma_g_g = - 0.5*titj*titj/l/l
  gamma_g_g = exp(gamma_g_g)

  gamma_g_dg = titj/l/l
  gamma_g_dg = gamma_g_g * gamma_g_dg

  gamma_dg_g = - gamma_g_dg

  gamma_dg_dg = ( 1./l/l - titj*titj/l/l/l/l)
  gamma_dg_dg = gamma_g_g * gamma_dg_dg

end subroutine get_EXP_gammas

subroutine get_M52_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(in) ::  pars(0:npars-1)
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: titj, expt, dt, sgn
  real(kind=mireal) :: sq5
  real(kind=mireal) :: l !Lambda parameter for the square exponential kernel


  l = pars(0)

  call fcdist(x1,x2,titj,nx1,nx2)

  sq5  = sqrt(5.)
  dt   = sq5*abs(titj)/l
  expt = exp(-dt)

  !Find the signs required for the absolute value derivative
  sgn = titj/abs(titj)
  !Fill with zero all the divisions by zero
  where (sgn .ne. sgn)
    sgn = 0.0
  end where

  gamma_g_g = expt * ( 1. + dt + dt*dt/3. )

  gamma_g_dg = 1. / 3. * sgn * expt * dt * (dt + 1.)
  gamma_g_dg = sq5 / l * gamma_g_dg

  gamma_dg_g = - gamma_g_dg

  gamma_dg_dg = 1. / 3. * (dt*dt - dt - 1.) * expt
  gamma_dg_dg =  - 5. / l / l * gamma_dg_dg

end subroutine get_M52_gammas

!This subroutine allow us to select what gammas do we want
subroutine get_gammas(x1,x2,pars,kernel,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  real(kind=mireal), intent(in) ::  pars(0:npars-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out), dimension(0:nx1-1,0:nx2-1) ::  gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  !

  if (kernel == 'MQ') then !Multi-Quasi-periodic Kernel
       call get_QP_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,3)
  else if (kernel == 'ME') then !Multi-Exponential Kernel
      call get_EXP_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,1)
  else if (kernel == 'MM') then !Multi-Matern 5/2 Kernel
      call get_M52_gammas(x1,x2,pars,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1,nx2,1)
  end if

end subroutine

subroutine MultidimCov(pars,x1,x2,kernel,ndim,cov,nx1,nx2)
  use constants
  implicit none

  !
  integer, intent(in) :: nx1, nx2, ndim
  real(kind=mireal), intent(in) :: pars(0:ndim*2+3-1)
  real(kind=mireal), intent(in) :: x1(0:nx1-1)
  real(kind=mireal), intent(in) :: x2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !

  if ( ndim == 1 ) then
     call Multi1(pars,x1,x2,kernel(1:2),cov,nx1,nx2)
  else if ( ndim == 2 ) then
     call Multi2(pars,x1,x2,kernel(1:2),cov,nx1,nx2)
  else if ( ndim == 3 ) then
     call Multi3(pars,x1,x2,kernel(1:2),cov,nx1,nx2)
  else if ( ndim == 4 ) then
     call Multi4(pars,x1,x2,kernel(1:2),cov,nx1,nx2)
  else if ( ndim == 5 ) then
     call Multi5(pars,x1,x2,kernel(1:2),cov,nx1,nx2)
  else
     call Multin(pars,x1,x2,kernel(1:2),ndim,cov,nx1,nx2)
  end if

end subroutine

!This subroutine computes the multi-dimensional GP matrix for one time-serie
!With a QP kernel
subroutine Multi1(pars,x1,x2,kernel,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:4)
  real(kind=mireal), intent(in) :: x1(0:nx1-1)
  real(kind=mireal), intent(in) :: x2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  integer, parameter :: m = 1
  integer :: npars
  real(kind=mireal), dimension(0:nx1-1,0:nx2-1) :: gamma_dg_dg, gamma_g_g,gamma_dg_g,gamma_g_dg
  real(kind=mireal) :: Ac, Ar
  real(kind=mireal) :: hyperpars(0:2)

  npars = size(pars) - 2

  Ac = pars(0)
  Ar = pars(1)
  hyperpars(:) = pars(2:)

  call get_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),hyperpars,kernel,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m,npars)

  cov = Ac*Ac*gamma_g_g + Ar*Ar*gamma_dg_dg

end subroutine Multi1

!This subroutine computes the multi-dimensional GP matrixor two time-series
!With a QP kernel
subroutine Multi2(pars,x1,x2,kernel,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:6) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  integer, parameter :: m = 2
  integer :: npars
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k11, k12
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k21, k22
  real(kind=mireal) :: Vc, Vr, Lc, Lr
  real(kind=mireal) :: hyperpars(0:2)

  npars = size(pars) - 4

  Vc = pars(0); Vr = pars(1); Lc = pars(2); Lr = pars(3)
  hyperpars(:) = pars(4:)

  call get_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),hyperpars,kernel,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m,npars)

  k11 = Vc**2 * gamma_g_g + Vr**2 * gamma_dg_dg

  k22 = Lc**2 * gamma_g_g + Lr**2 * gamma_dg_dg

  k12 = Vc*Lc*gamma_g_g + Vr*Lr*gamma_dg_dg + (Vc*Lr - Vr*Lc)*gamma_g_dg

  if (nx1 == nx2) then
    k21 = transpose(k12)
  else
    k21 = Lc*Vc*gamma_g_g + Lr*Vr*gamma_dg_dg + (Lc*Vr - Lr*Vc)*gamma_g_dg
  end if

  !
  cov(0*nx1/m:1*nx1/m-1,0*nx2/m:1*nx2/m-1) = k11
  cov(0*nx1/m:1*nx1/m-1,1*nx2/m:2*nx2/m-1) = k12
  !
  cov(1*nx1/m:2*nx1/m-1,0*nx2/m:1*nx2/m-1) = k21
  cov(1*nx1/m:2*nx1/m-1,1*nx2/m:2*nx2/m-1) = k22
  !

end subroutine Multi2

!
subroutine Multi3(pars,x1,x2,kernel,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:8) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  integer, parameter :: m = 3
  integer :: npars
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k11, k12, k13
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k21, k22, k23
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k31, k32, k33
  real(kind=mireal) :: a1, b1, a2, b2, a3, b3
  real(kind=mireal) :: hyperpars(0:2)

  npars = size(pars) - 6

  a1 = pars(0); b1 = pars(1); a2 = pars(2); b2 = pars(3); a3 = pars(4); b3 = pars(5)
  hyperpars(:) = pars(6:)

  call get_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),hyperpars,kernel,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m,npars)

  k11 = a1**2 * gamma_g_g + b1**2 * gamma_dg_dg

  k22 = a2**2 * gamma_g_g + b2**2 * gamma_dg_dg

  k33 = a3**2 * gamma_g_g + b3**2 * gamma_dg_dg

  k12 = a1*a2*gamma_g_g + b1*b2*gamma_dg_dg + (a1*b2 - b1*a2)*gamma_g_dg

  k13 = a1*a3*gamma_g_g + b1*b3*gamma_dg_dg + (a1*b3 - b1*a3)*gamma_g_dg

  k23 = a2*a3*gamma_g_g + b2*b3*gamma_dg_dg + (a2*b3 - b2*a3)*gamma_g_dg

  if (nx1 == nx2) then
    k21 = transpose(k12)
    k31 = transpose(k13)
    k32 = transpose(k23)
  else
    k21 = a2*a1*gamma_g_g + b2*b1*gamma_dg_dg + (a2*b1 - a1*b2)*gamma_g_dg
    k31 = a3*a1*gamma_g_g + b3*b1*gamma_dg_dg + (a3*b1 - a1*b3)*gamma_g_dg
    k32 = a3*a2*gamma_g_g + b3*b2*gamma_dg_dg + (a3*b2 - a2*b3)*gamma_g_dg
  end if


  cov(0*nx1/m:1*nx1/m-1,0*nx2/m:1*nx2/m-1) = k11
  cov(0*nx1/m:1*nx1/m-1,1*nx2/m:2*nx2/m-1) = k12
  cov(0*nx1/m:1*nx1/m-1,2*nx2/m:3*nx2/m-1) = k13
  !
  cov(1*nx1/m:2*nx1/m-1,0*nx2/m:1*nx2/m-1) = k21
  cov(1*nx1/m:2*nx1/m-1,1*nx2/m:2*nx2/m-1) = k22
  cov(1*nx1/m:2*nx1/m-1,2*nx2/m:3*nx2/m-1) = k23
  !
  cov(2*nx1/m:3*nx1/m-1,0*nx2/m:1*nx2/m-1) = k31
  cov(2*nx1/m:3*nx1/m-1,1*nx2/m:2*nx2/m-1) = k32
  cov(2*nx1/m:3*nx1/m-1,2*nx2/m:3*nx2/m-1) = k33


end subroutine Multi3
!!
subroutine Multi4(pars,x1,x2,kernel,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:10)
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  integer, parameter :: m = 4
  integer :: npars
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k11, k12, k13, k14
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k21, k22, k23, k24
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k31, k32, k33, k34
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k41, k42, k43, k44
  real(kind=mireal) :: a1, b1, a2, b2, a3, b3, a4, b4
  real(kind=mireal) :: hyperpars(0:2)

  npars = size(pars) - 8

  a1 = pars(0); b1 = pars(1); a2 = pars(2); b2 = pars(3); a3 = pars(4); b3 = pars(5); a4 = pars(6); b4 = pars(7)
  hyperpars(:) = pars(8:)

  call get_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),hyperpars,kernel,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m,npars)

  k11 = a1**2 * gamma_g_g + b1**2 * gamma_dg_dg

  k22 = a2**2 * gamma_g_g + b2**2 * gamma_dg_dg

  k33 = a3**2 * gamma_g_g + b3**2 * gamma_dg_dg

  k44 = a4**2 * gamma_g_g + b4**2 * gamma_dg_dg

  k12 = a1*a2*gamma_g_g + b1*b2*gamma_dg_dg + (a1*b2 - a2*b1)*gamma_g_dg

  k13 = a1*a3*gamma_g_g + b1*b3*gamma_dg_dg + (a1*b3 - a3*b1)*gamma_g_dg

  k14 = a1*a4*gamma_g_g + b1*b4*gamma_dg_dg + (a1*b4 - a4*b1)*gamma_g_dg

  k23 = a2*a3*gamma_g_g + b2*b3*gamma_dg_dg + (a2*b3 - a3*b2)*gamma_g_dg

  k24 = a2*a4*gamma_g_g + b2*b4*gamma_dg_dg + (a2*b4 - a4*b2)*gamma_g_dg

  k34 = a3*a4*gamma_g_g + b3*b4*gamma_dg_dg + (a3*b4 - a4*b3)*gamma_g_dg

  if (nx1 == nx2) then
    k21 = transpose(k12)
    k31 = transpose(k13)
    k32 = transpose(k23)
    k41 = transpose(k14)
    k42 = transpose(k24)
    k43 = transpose(k34)
  else
    k21 = a2*a1*gamma_g_g + b2*b1*gamma_dg_dg + (a2*b1 - a1*b2)*gamma_g_dg
    k31 = a3*a1*gamma_g_g + b3*b1*gamma_dg_dg + (a3*b1 - a1*b3)*gamma_g_dg
    k32 = a3*a2*gamma_g_g + b3*b2*gamma_dg_dg + (a3*b2 - a2*b3)*gamma_g_dg
    k41 = a4*a1*gamma_g_g + b4*b1*gamma_dg_dg + (a4*b1 - a1*b4)*gamma_g_dg
    k42 = a4*a2*gamma_g_g + b4*b2*gamma_dg_dg + (a4*b2 - a2*b4)*gamma_g_dg
    k43 = a4*a3*gamma_g_g + b4*b3*gamma_dg_dg + (a4*b3 - a3*b4)*gamma_g_dg
  end if


  !
  !Time to fill the covariance matrix
  !
  cov(0*nx1/m:1*nx1/m-1,0*nx2/m:1*nx2/m-1) = k11
  cov(0*nx1/m:1*nx1/m-1,1*nx2/m:2*nx2/m-1) = k12
  cov(0*nx1/m:1*nx1/m-1,2*nx2/m:3*nx2/m-1) = k13
  cov(0*nx1/m:1*nx1/m-1,3*nx2/m:4*nx2/m-1) = k14
  !
  cov(1*nx1/m:2*nx1/m-1,0*nx2/m:1*nx2/m-1) = k21
  cov(1*nx1/m:2*nx1/m-1,1*nx2/m:2*nx2/m-1) = k22
  cov(1*nx1/m:2*nx1/m-1,2*nx2/m:3*nx2/m-1) = k23
  cov(1*nx1/m:2*nx1/m-1,3*nx2/m:4*nx2/m-1) = k24
  !
  cov(2*nx1/m:3*nx1/m-1,0*nx2/m:1*nx2/m-1) = k31
  cov(2*nx1/m:3*nx1/m-1,1*nx2/m:2*nx2/m-1) = k32
  cov(2*nx1/m:3*nx1/m-1,2*nx2/m:3*nx2/m-1) = k33
  cov(2*nx1/m:3*nx1/m-1,3*nx2/m:4*nx2/m-1) = k34
  !
  cov(3*nx1/m:4*nx1/m-1,0*nx2/m:1*nx2/m-1) = k41
  cov(3*nx1/m:4*nx1/m-1,1*nx2/m:2*nx2/m-1) = k42
  cov(3*nx1/m:4*nx1/m-1,2*nx2/m:3*nx2/m-1) = k43
  cov(3*nx1/m:4*nx1/m-1,3*nx2/m:4*nx2/m-1) = k44
  !


end subroutine Multi4
!
subroutine Multi5(pars,x1,x2,kernel,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:12)
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !
  integer, parameter :: m = 5
  integer :: npars
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_dg_g, gamma_g_dg
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k11, k12, k13, k14, k15
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k21, k22, k23, k24, k25
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k31, k32, k33, k34, k35
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k41, k42, k43, k44, k45
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: k51, k52, k53, k54, k55
  real(kind=mireal) :: a1, b1, a2, b2, a3, b3, a4, b4, a5, b5
  real(kind=mireal) :: hyperpars(0:2)

  npars = size(pars) - 10

  a1 = pars(0); b1 = pars(1); a2 = pars(2); b2 = pars(3); a3 = pars(4); b3 = pars(5); a4 = pars(6); b4 = pars(7)
  a5 = pars(8); b5 = pars(9)
  hyperpars(:) = pars(10:)

  call get_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),hyperpars,kernel,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m,npars)

  k11 = a1**2 * gamma_g_g + b1**2 * gamma_dg_dg

  k22 = a2**2 * gamma_g_g + b2**2 * gamma_dg_dg

  k33 = a3**2 * gamma_g_g + b3**2 * gamma_dg_dg

  k44 = a4**2 * gamma_g_g + b4**2 * gamma_dg_dg

  k55 = a5**2 * gamma_g_g + b5**2 * gamma_dg_dg

  k12 = a1*a2*gamma_g_g + b1*b2*gamma_dg_dg + (a1*b2 - a2*b1)*gamma_g_dg

  k13 = a1*a3*gamma_g_g + b1*b3*gamma_dg_dg + (a1*b3 - a3*b1)*gamma_g_dg

  k14 = a1*a4*gamma_g_g + b1*b4*gamma_dg_dg + (a1*b4 - a4*b1)*gamma_g_dg

  k15 = a1*a5*gamma_g_g + b1*b5*gamma_dg_dg + (a1*b5 - a5*b1)*gamma_g_dg

  k23 = a2*a3*gamma_g_g + b2*b3*gamma_dg_dg + (a2*b3 - a3*b2)*gamma_g_dg

  k24 = a2*a4*gamma_g_g + b2*b4*gamma_dg_dg + (a2*b4 - a4*b2)*gamma_g_dg

  k25 = a2*a5*gamma_g_g + b2*b5*gamma_dg_dg + (a2*b5 - a5*b2)*gamma_g_dg

  k34 = a3*a4*gamma_g_g + b3*b4*gamma_dg_dg + (a3*b4 - a4*b3)*gamma_g_dg

  k35 = a3*a5*gamma_g_g + b3*b5*gamma_dg_dg + (a3*b5 - a5*b3)*gamma_g_dg

  k45 = a4*a5*gamma_g_g + b4*b5*gamma_dg_dg + (a4*b5 - a5*b4)*gamma_g_dg

  if (nx1 == nx2) then
    k21 = transpose(k12)
    k31 = transpose(k13)
    k32 = transpose(k23)
    k41 = transpose(k14)
    k42 = transpose(k24)
    k43 = transpose(k34)
    k51 = transpose(k15)
    k52 = transpose(k25)
    k53 = transpose(k35)
    k54 = transpose(k45)
  else
    k21 = a2*a1*gamma_g_g + b2*b1*gamma_dg_dg + (a2*b1 - a1*b2)*gamma_g_dg
    k31 = a3*a1*gamma_g_g + b3*b1*gamma_dg_dg + (a3*b1 - a1*b3)*gamma_g_dg
    k32 = a3*a2*gamma_g_g + b3*b2*gamma_dg_dg + (a3*b2 - a2*b3)*gamma_g_dg
    k41 = a4*a1*gamma_g_g + b4*b1*gamma_dg_dg + (a4*b1 - a1*b4)*gamma_g_dg
    k42 = a4*a2*gamma_g_g + b4*b2*gamma_dg_dg + (a4*b2 - a2*b4)*gamma_g_dg
    k43 = a4*a3*gamma_g_g + b4*b3*gamma_dg_dg + (a4*b3 - a3*b4)*gamma_g_dg
    k51 = a5*a1*gamma_g_g + b5*b1*gamma_dg_dg + (a5*b1 - a1*b5)*gamma_g_dg
    k52 = a5*a2*gamma_g_g + b5*b2*gamma_dg_dg + (a5*b2 - a2*b5)*gamma_g_dg
    k53 = a5*a3*gamma_g_g + b5*b3*gamma_dg_dg + (a5*b3 - a3*b5)*gamma_g_dg
    k54 = a5*a4*gamma_g_g + b5*b4*gamma_dg_dg + (a5*b4 - a4*b5)*gamma_g_dg
  end if
  !
  !Time to fill the covariance matrix
  !
  cov(0*nx1/m:1*nx1/m-1,0*nx2/m:1*nx2/m-1) = k11
  cov(0*nx1/m:1*nx1/m-1,1*nx2/m:2*nx2/m-1) = k12
  cov(0*nx1/m:1*nx1/m-1,2*nx2/m:3*nx2/m-1) = k13
  cov(0*nx1/m:1*nx1/m-1,3*nx2/m:4*nx2/m-1) = k14
  cov(0*nx1/m:1*nx1/m-1,4*nx2/m:5*nx2/m-1) = k15
  !
  cov(1*nx1/m:2*nx1/m-1,0*nx2/m:1*nx2/m-1) = k21
  cov(1*nx1/m:2*nx1/m-1,1*nx2/m:2*nx2/m-1) = k22
  cov(1*nx1/m:2*nx1/m-1,2*nx2/m:3*nx2/m-1) = k23
  cov(1*nx1/m:2*nx1/m-1,3*nx2/m:4*nx2/m-1) = k24
  cov(1*nx1/m:2*nx1/m-1,4*nx2/m:5*nx2/m-1) = k25
  !
  cov(2*nx1/m:3*nx1/m-1,0*nx2/m:1*nx2/m-1) = k31
  cov(2*nx1/m:3*nx1/m-1,1*nx2/m:2*nx2/m-1) = k32
  cov(2*nx1/m:3*nx1/m-1,2*nx2/m:3*nx2/m-1) = k33
  cov(2*nx1/m:3*nx1/m-1,3*nx2/m:4*nx2/m-1) = k34
  cov(2*nx1/m:3*nx1/m-1,4*nx2/m:5*nx2/m-1) = k35
  !
  cov(3*nx1/m:4*nx1/m-1,0*nx2/m:1*nx2/m-1) = k41
  cov(3*nx1/m:4*nx1/m-1,1*nx2/m:2*nx2/m-1) = k42
  cov(3*nx1/m:4*nx1/m-1,2*nx2/m:3*nx2/m-1) = k43
  cov(3*nx1/m:4*nx1/m-1,3*nx2/m:4*nx2/m-1) = k44
  cov(3*nx1/m:4*nx1/m-1,4*nx2/m:5*nx2/m-1) = k45
  !
  cov(4*nx1/m:5*nx1/m-1,0*nx2/m:1*nx2/m-1) = k51
  cov(4*nx1/m:5*nx1/m-1,1*nx2/m:2*nx2/m-1) = k52
  cov(4*nx1/m:5*nx1/m-1,2*nx2/m:3*nx2/m-1) = k53
  cov(4*nx1/m:5*nx1/m-1,3*nx2/m:4*nx2/m-1) = k54
  cov(4*nx1/m:5*nx1/m-1,4*nx2/m:5*nx2/m-1) = k55
  !


end subroutine Multi5

subroutine Multin(pars,x1,x2,kernel,m,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, m
  real(kind=mireal), intent(in) :: pars(0:m*2+3-1) !m*2 amplitudes, 3 parameters for the QP kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1), x2(0:nx2-1)
  character(len=2), intent(in)  :: kernel
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
 !
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1) :: gamma_g_g, gamma_dg_dg, gamma_g_dg, gamma_dg_g
  real(kind=mireal), dimension(0:nx1/m-1,0:nx2/m-1,0:m-1,0:m-1) :: kas
  real(kind=mireal) :: Amplitudes(0:m*2-1)
  real(kind=mireal) :: hyperpars(0:2)
  integer :: i, j, npars

  !Read the parameters
  Amplitudes(:) = pars(0:m*2-1)
  npars = size(pars) - m*2
  hyperpars = pars(m*2:)

  call get_gammas(x1(0:nx1/m-1),x2(0:nx2/m-1),hyperpars,kernel,gamma_g_g,gamma_g_dg,gamma_dg_g,gamma_dg_dg,nx1/m,nx2/m,npars)

  if ( nx1 .ne. nx2  ) then !Compute the K's for not squared matrices

  do i = 0, m - 1
    do j = 0, m - 1
      kas(:,:,i,j) = amplitudes(i*2)*amplitudes(j*2)*gamma_g_g + &
                     amplitudes(i*2+1)*amplitudes(j*2+1)*gamma_dg_dg + &
                     ( amplitudes(i*2)*amplitudes(j*2+1) - amplitudes(i*2+1)*amplitudes(j*2) )* gamma_g_dg
      cov(i*nx1/m:(i+1)*nx1/m-1,j*nx2/m:(j+1)*nx2/m-1) = kas(:,:,i,j)
    end do
  end do

  else !Compute the K's for square matrices

  do i = 0, m - 1
    do j = i, m - 1
      kas(:,:,i,j) = amplitudes(i*2)*amplitudes(j*2)*gamma_g_g + &
                     amplitudes(i*2+1)*amplitudes(j*2+1)*gamma_dg_dg + &
                     ( amplitudes(i*2)*amplitudes(j*2+1) - amplitudes(i*2+1)*amplitudes(j*2) )* gamma_g_dg
      cov(i*nx1/m:(i+1)*nx1/m-1,j*nx2/m:(j+1)*nx2/m-1) = kas(:,:,i,j)
    end do
  end do

  do i = 1, m - 1
    do j = 0, i - 1
      kas(:,:,i,j) = transpose(kas(:,:,j,i))
      cov(i*nx1/m:(i+1)*nx1/m-1,j*nx2/m:(j+1)*nx2/m-1) = kas(:,:,i,j)
    end do
  end do

  end if

end subroutine
