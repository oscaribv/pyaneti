module kernels
implicit none

  contains
!---------------------------------------------------------------------
!          General routines to compute bla bla
!---------------------------------------------------------------------
  !This routine is used to select the desired kernel
  subroutine covfunc(kernel,pars,x1,x2,cov,nx1,nx2,npars)
  implicit none
  !
  integer, intent(in) :: nx1, nx2, npars
  character (len=30), intent(in) :: kernel
  double precision, intent(in) :: pars(0:npars-1) !There are only two parameters for this kernel
  double precision, intent(in) :: x1(0:nx1-1)
  double precision, intent(in) :: x2(0:nx2-1)
  double precision, intent(out) :: cov(0:nx1-1,0:nx2-1)

  if ( kernel == 'SE' ) then
     call SEKernel(pars,x1,x2,cov,nx1,nx2)
  else if ( kernel == 'M32' ) then
     call M32Kernel(pars,x1,x2,cov,nx1,nx2)
  else if ( kernel == 'QP' ) then
  else
     call QPKernel(pars,x1,x2,cov,nx1,nx2)
     print *, 'Kernel ', kernel,' is not defined!'
     stop
  end if

  end subroutine covfunc

  subroutine pred_gp(kernel,covpar,xobs,yobs,eobs,xtest,mvec,cov,nobs,ntest,npar)
  implicit none
  !
  integer, intent(in) :: nobs, ntest, npar
  character (len=30), intent(in) :: kernel
  double precision, intent(in) :: covpar(0:npar-1)
  double precision, dimension(0:nobs-1), intent(in) :: xobs, yobs, eobs
  double precision, dimension(0:ntest-1), intent(in) :: xtest
  double precision, dimension(0:ntest-1), intent(out) :: mvec
  double precision, dimension(0:ntest-1,0:ntest-1), intent(out) :: cov
  !local variables
  double precision :: K(0:nobs-1,0:nobs-1)
  double precision :: dummy(0:nobs-1,0:nobs-1)
  double precision :: Ki(0:nobs-1,0:nobs-1)
  double precision :: Kss(0:ntest-1,0:ntest-1)
  double precision :: Ks(0:ntest-1,0:nobs-1)

  !covariance matrix for observed vector
  call covfunc(kernel,covpar,xobs,xobs,K,nobs,nobs,npar)
  call fill_diag(eobs*eobs,dummy,nobs)
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

  subroutine NLL_GP(p,kernel,x,y,e,nll,np,nx)
  implicit none
  !
  integer, intent(in) :: np, nx
  double precision, dimension(0:np-1), intent(in) :: p
  double precision, dimension(0:nx-1), intent(in) :: x, y, e
  character (len=30), intent(in) :: kernel
  double precision, intent(out) :: nll
  !
  !local variables
  double precision :: K(0:nx-1,0:nx-1)
  double precision :: dummy(0:nx-1,0:nx-1)
  double precision :: Ki(0:nx-1,0:nx-1)
  double precision :: nl1, nl2
  double precision, parameter :: log_two_pi = log(2.d0*acos(-1.d0))

  !covariance matrix for observed vector
  call covfunc(kernel,p,x,x,K,nx,nx,np)
  call fill_diag(e*e,dummy,nx)
  K = K + dummy
  !Get the inverse matrix
  dummy = K
  call inverse(dummy,Ki,nx)

  nl1 = dot_product(y,matmul(Ki,y))
  print *, 'nl1', nl1
  !
  call findlogddet(K,nl2,nx)
  print *, 'nl2', nl2
  !
  nll = 5.d-1*(nl1 + nl2 + nx * log_two_pi)
  print *,'nll', nll

  end subroutine NLL_GP

!---------------------------------------------------------------------
!                         KERNELS
!---------------------------------------------------------------------

  subroutine SEKernel(pars,x1,x2,cov,nx1,nx2)
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  double precision, intent(in) :: pars(0:1) !There are only two parameters for this kernel
  double precision, intent(in) :: x1(0:nx1-1)
  double precision, intent(in) :: x2(0:nx2-1)
  double precision, intent(out) :: cov(0:nx1-1,0:nx2-1)

  !Get the x_i - x_j
  call fcdist(x1,x2,cov,nx1,nx2)
  cov = pars(0) * exp( - pars(1) * cov * cov )

  end subroutine SEKernel

  subroutine M32Kernel(pars,x1,x2,cov,nx1,nx2)
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  double precision, intent(in) :: pars(0:1) !There are only two parameters for this kernel
  double precision, intent(in) :: x1(0:nx1-1)
  double precision, intent(in) :: x2(0:nx2-1)
  double precision, intent(out) :: cov(0:nx1-1,0:nx2-1)
  !A = pars(0), Gamma = pars(1)

  !Get the x_i - x_j
  call fcdist(x1,x2,cov,nx1,nx2)
  cov = pars(1)*cov*cov
  cov = pars(0) * ( 1.d0 + sqrt(3.d0*cov)) * exp(-sqrt(3.d0*cov))

  end subroutine M32Kernel

  subroutine QPKernel(pars,x1,x2,cov,nx1,nx2)
  implicit none
  !
  integer, intent(in) :: nx1, nx2
  double precision, intent(in) :: pars(0:3) !There are only two parameters for this kernel
  double precision, intent(in) :: x1(0:nx1-1)
  double precision, intent(in) :: x2(0:nx2-1)
  double precision, intent(out) :: cov(0:nx1-1,0:nx2-1)
  !A = pars(0), Gamma_1 = pars(1), Gamma_2 = pars(2), P = pars(3)
  !
  double precision, parameter :: pi = acos(-1.d0)

  !Get the x_i - x_j
  call fcdist(x1,x2,cov,nx1,nx2)
  cov = pars(0) * exp( &
        - pars(1) * ( sin( pi * cov / pars(3) ) )**2 &
        - pars(2) * cov * cov )

  end subroutine QPKernel

end module kernels