!---------------------------------------------------------------------
!          General routines to compute bla bla
!---------------------------------------------------------------------
  !This routine is used to select the desired kernel
  subroutine covfunc(kernel,pars,x1,x2,cov,nx1,nx2,npars)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  character (len=3), intent(in) :: kernel
  real(kind=mireal), intent(in) :: pars(0:npars-1) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1)
  real(kind=mireal), intent(in) :: x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)

  if ( kernel == 'SEK' ) then
     call SEKernel(pars,x1,x2,cov,nx1,nx2)
  else if ( kernel == 'M32' ) then
     call M32Kernel(pars,x1,x2,cov,nx1,nx2)
  else if ( kernel == 'QPK' ) then
     call QPKernel(pars,x1,x2,cov,nx1,nx2)
  else if ( kernel == 'QP2' ) then
     call QPKernel2(pars,x1,x2,cov,nx1,nx2)
  else if ( kernel == 'VRF' ) then
     call QPsuperKernel(pars,x1,cov,nx1)
  else
     print *, 'ERROR: Kernel ', kernel,' is not defined!'
     stop
  end if

  end subroutine covfunc

  subroutine pred_gp(kernel,covpar,xobs,yobs,eobs,xtest,jit,ljit,mvec,cov,nobs,ntest,npar,njit)
  use constants
  implicit none
  !
  integer, intent(in) :: nobs, ntest, npar, njit
  character (len=3), intent(in) :: kernel
  real(kind=mireal), intent(in) :: covpar(0:npar-1)
  real(kind=mireal), dimension(0:nobs-1), intent(in) :: xobs, yobs, eobs
  real(kind=mireal), dimension(0:ntest-1), intent(in) :: xtest
  integer, dimension(0:nobs-1), intent(in):: ljit
  real(kind=mireal), dimension(0:njit-1), intent(in) :: jit
  real(kind=mireal), dimension(0:ntest-1), intent(out) :: mvec
  real(kind=mireal), dimension(0:ntest-1,0:ntest-1), intent(out) :: cov
  !local variables
  real(kind=mireal), dimension(0:nobs-1) :: s2
  real(kind=mireal) :: K(0:nobs-1,0:nobs-1)
  real(kind=mireal) :: dummy(0:nobs-1,0:nobs-1)
  real(kind=mireal) :: Ki(0:nobs-1,0:nobs-1)
  real(kind=mireal) :: Kss(0:ntest-1,0:ntest-1)
  real(kind=mireal) :: Ks(0:ntest-1,0:nobs-1)

  !covariance matrix for observed vector
  call covfunc(kernel,covpar,xobs,xobs,K,nobs,nobs,npar)
  s2 = eobs**2 + jit(ljit)**2
  call fill_diag(s2,dummy,nobs)
  K = K + dummy
  !covariance matrix for test vector
  call covfunc(kernel,covpar,xtest,xtest,Kss,ntest,ntest,npar)
  !covariance matrix for cross vector
  call covfunc(kernel,covpar,xtest,xobs,Ks,ntest,nobs,npar)
  !Get the inverse matrix
  dummy = K
  call inverse(dummy,Ki,nobs)

  !
  mvec = matmul(Ks,matmul(Ki,yobs))
  !
  cov = Kss - matmul(Ks,matmul(Ki,transpose(ks)))

  end subroutine pred_gp

  subroutine NLL_GP(p,kernel,x,y,e,jit,ljit,nll,chi2,np,nx,njit)
  use constants
  implicit none
  !
  integer, intent(in) :: np, nx, njit
  real(kind=mireal), dimension(0:np-1), intent(in) :: p
  real(kind=mireal), dimension(0:njit-1), intent(in) :: jit
  real(kind=mireal), dimension(0:nx-1), intent(in) :: x, y, e
  integer, dimension(0:nx-1), intent(in):: ljit
  character (len=3), intent(in) :: kernel
  real(kind=mireal), intent(out) :: nll, chi2
  !
  !local variables
  real(kind=mireal), dimension(0:nx-1) :: s2
  real(kind=mireal) :: K(0:nx-1,0:nx-1)
  real(kind=mireal) :: dummy(0:nx-1,0:nx-1)
  real(kind=mireal) :: Ki(0:nx-1,0:nx-1)
  real(kind=mireal) :: nl1, nl2
  integer :: i

  !covariance matrix for observed vector
  s2 = e**2 + jit(ljit)**2
  call covfunc(kernel,p,x,x,K,nx,nx,np)
  dummy = K
  call fill_diag(s2,dummy,nx)
  K = K + dummy
  !Get the inverse matrix
  dummy = K
  call inverse(dummy,Ki,nx)

  nl1 = dot_product(y,matmul(Ki,y))
  !
  call findlogddet(K,nl2,nx)
  !
  chi2 = nl1

  !If the code arrives here, is because the computed covariance matrix is not
  !a positive-definite matrix, let us avoid this solution here
  if (chi2 < 0. ) chi2 = huge(chi2)

  nll = 5.d-1*(chi2 + nl2 + nx * log_two_pi)


  end subroutine NLL_GP


  subroutine NLL_GP_py(p,kernel,x,y,e,jit,ljit,nll,np,nx,njit)
  use constants
  implicit none
  !
  integer, intent(in) :: np, nx, njit
  real(kind=mireal), dimension(0:np-1), intent(in) :: p
  real(kind=mireal), dimension(0:njit-1), intent(in) :: jit
  real(kind=mireal), dimension(0:nx-1), intent(in) :: x, y, e
  integer, dimension(0:nx-1), intent(in):: ljit
  character (len=3), intent(in) :: kernel
  real(kind=mireal), intent(out) :: nll
  !
  !local variables
  real(kind=mireal) :: chi2

  call NLL_GP(p,kernel,x,y,e,jit,ljit,nll,chi2,np,nx,njit)

  end subroutine NLL_GP_py

!---------------------------------------------------------------------
!                         KERNELS
!---------------------------------------------------------------------

  subroutine SEKernel(pars,x1,x2,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:1) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1)
  real(kind=mireal), intent(in) :: x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)

  !Get the x_i - x_j
  call fcdist(x1,x2,cov,nx1,nx2)
  cov = pars(0) * exp( - pars(1) * cov * cov )

  end subroutine SEKernel

  subroutine M32Kernel(pars,x1,x2,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:1) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1)
  real(kind=mireal), intent(in) :: x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !A = pars(0), Gamma = pars(1)

  !Get the x_i - x_j
  call fcdist(x1,x2,cov,nx1,nx2)
  cov = pars(1)*cov*cov
  cov = pars(0) * ( 1.d0 + sqrt(3.d0*cov)) * exp(-sqrt(3.d0*cov))

  end subroutine M32Kernel

  subroutine QPKernel(pars,x1,x2,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:3) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1)
  real(kind=mireal), intent(in) :: x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !A = pars(0), lambda_p = pars(1), lambda_q = pars(2), P = pars(3)
  !
  real(kind=mireal) :: A, lp, le, P

  A  = pars(0)
  le = pars(1)
  lp = pars(2)
  P  = pars(3)

  !Get the x_i - x_j
  call fcdist(x1,x2,cov,nx1,nx2)
  cov = - (sin(pi*cov/P))**2/2./lp**2 - cov**2/2./le**2
  cov = A * exp(cov)

  end subroutine QPKernel


subroutine QPKernel2(pars,x1,x2,cov,nx1,nx2)
  use constants
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  real(kind=mireal), intent(in) :: pars(0:3) !There are only two parameters for this kernel
  real(kind=mireal), intent(in) :: x1(0:nx1-1)
  real(kind=mireal), intent(in) :: x2(0:nx2-1)
  real(kind=mireal), intent(out) :: cov(0:nx1-1,0:nx2-1)
  !h = pars(0), Gamma_1 = pars(1), Gamma_2 = pars(2), P = pars(3)

  !Get the x_i - x_j
  call fcdist(x1,x2,cov,nx1,nx2)
  cov = pars(0)**2 * exp( &
        - pars(1) * ( sin( pi * cov / pars(3) ) )**2 &
        - cov * cov / (2.*pars(2)**2) )

  end subroutine QPKernel2