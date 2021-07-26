 !https://ww2.odu.edu/~agodunov/computing/programs/book2/Ch06/Inverse.f90
  subroutine inverse(m,cc,n)
  !============================================================
  ! Inverse matrix
  ! Method: Based on Doolittle LU factorization for Ax=b
  ! Alex G. December 2009
  !-----------------------------------------------------------
  ! input ...
  ! a(n,n) - array of coefficients for matrix A
  ! n      - dimension
  ! output ...
  ! c(n,n) - inverse matrix of A
  ! comments ...
  ! the original matrix a(n,n) will be destroyed
  ! during the calculation
!===========================================================
  use constants
  implicit none
    integer, intent(in) :: n
    real(kind=mireal), intent(in) :: m(0:n-1,0:n-1)
    real(kind=mireal), intent(out) :: cc(0:n-1,0:n-1)
    !
    real(kind=mireal) :: a(n,n), c(n,n), L(n,n), U(n,n), b(n), d(n), x(n)
    real(kind=mireal) :: coeff
    integer :: i, j, k

    a = m

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    !  step 1: forward elimination
    do k=1, n-1
      do i=k+1,n
        coeff=a(i,k)/a(k,k)
        L(i,k) = coeff
        do j=k+1,n
           a(i,j) = a(i,j)-coeff*a(k,j)
        end do
      end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
         U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do

    cc = c

  end subroutine inverse

!Subroutine to find the determinant of a square matrix
!Taken from https://github.com/ashwith/workspace/blob/master/fortran/determinant.f90
!Modified as a subroutine
!----------------------------------------------------------------------------------------
!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
  subroutine findlogddet(a,det,n)
  use constants
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  real(kind=mireal), DIMENSION(n,n), intent(in) :: a
  real(kind=mireal), intent (out) :: det
  real(kind=mireal) :: m, temp
  real(kind=mireal), DIMENSION(n,n) :: matrix
  INTEGER :: i, j, k, l
  LOGICAL :: DetExists = .TRUE.
  matrix = a
  l = 1
  !Convert to upper triangular form
  DO k = 1, n-1
    IF (matrix(k,k) == 0) THEN
      DetExists = .FALSE.
      DO i = k+1, n
        IF (matrix(i,k) /= 0) THEN
          DO j = 1, n
            temp = matrix(i,j)
            matrix(i,j)= matrix(k,j)
            matrix(k,j) = temp
          END DO
          DetExists = .TRUE.
          l=-l
          EXIT
        ENDIF
      END DO
      IF (DetExists .EQV. .FALSE.) THEN
        det = 0
        return
      END IF
    ENDIF
    DO j = k+1, n
      m = matrix(j,k)/matrix(k,k)
      DO i = k+1, n
        matrix(j,i) = matrix(j,i) - m*matrix(k,i)
      END DO
    END DO
  END DO

  !Calculate determinant by finding product of diagonal elements
  det = log(abs(real(l)))
  DO i = 1, n
    det = det + log(abs(matrix(i,i)))
  END DO

  end subroutine findlogddet

  subroutine fill_diag(v,M,n)
  use constants
  implicit none
  !
  integer, intent(in) :: n
  real(kind=mireal), intent(in) :: v(0:n-1)
  real(kind=mireal), intent(out) :: M(0:n-1,0:n-1)
  !
  integer :: i

  M = 0.d0
  do i = 0, n - 1
    M(i,i) = v(i)
  end do

  end subroutine fill_diag

  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  ! Taken and modified from http://fortranwiki.org/fortran/show/Matrix+inversion
  subroutine cholesky_inv_det(A,Ainv,detlog,n)
   use constants
    integer, intent(in) :: n
    real(kind=mireal), dimension(n,n), intent(in) :: A
    real(kind=mireal), dimension(n,n), intent(out) :: Ainv
    real(kind=mireal), intent(out) :: detlog
    !
    integer :: info, i

    ! External procedures defined in LAPACK
    external DPOTRF
    external DPOTRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! https://www.math.utah.edu/software/lapack/lapack-d/dpotrf.html
    call DPOTRF('U',n,Ainv,n,info)

    !Now Ainv = L, such that A = L L*
    !The determinant of det A is then det A = det L x det L* = (det L)^2
    !and det L = mult_1^N L(i,i)

    !Since we are interested on log det A (for the likelihood calculation)
    !log(det(A)) = 2 (log (det (L)))

    detlog = 0.d0
    do i = 1, n
      detlog = detlog + log((Ainv(i,i)))
    end do
    detlog = 2.*detlog

    !DPOTRI calculates the inverse matrix using as input the output of DPOTRF
    !http://www.math.utah.edu/software/lapack/lapack-d/dpotri.html
    call DPOTRI('U',n,Ainv,n,info)

    do i = 2, n
  !    do j = 1, i
  !      Ainv(i,j) = Ainv(j,i)
        Ainv(i,1:i) = Ainv(1:i,i)
  !    end do
    end do

  end subroutine cholesky_inv_det

  subroutine mtops(M,V,n)
  use constants
  implicit none
    integer, intent(in) :: n
    real(kind=mireal), dimension(n,n), intent(in) :: M
    real(kind=mireal), dimension(n*(n+1)/2), intent(out) :: V
    !
    integer :: i, j

    do j = 1, n
      do i = 1, j
        V(i-1 + j*(j-1)/2) = M(i,j)
      end do
    end do

  end subroutine mtops


  subroutine pstom(V,M,n)
  use constants
  implicit none
    integer, intent(in) :: n
    real(kind=mireal), dimension(n*(n+1)/2), intent(in) :: V
    real(kind=mireal), dimension(n,n), intent(out) :: M
    !
    integer :: i, j

    M = 0.d0
    do j = 1, n
      do i = 1, j
        M(i,j) = V(i-1 + j*(j-1)/2)
        if (i == j) M(i,j) = M(i,j)/2.
      end do
    end do

    M = M + transpose(M)

  end subroutine pstom

!Taken from
!https://msvlab.hre.ntou.edu.tw/grades/now/inte/Inverse%20%26%20Border/border-LuTT.pdf
!This will be used to invert the first two terms of the mega matrix with
subroutine invert_square_block_matrices(a,b,d,xi,detlog,n)
use constants
  integer, intent(in) :: n
  real(kind=mireal), dimension(n,n), intent(in) :: a
  real(kind=mireal), dimension(n,n), intent(in) :: b
  real(kind=mireal), dimension(n,n), intent(in) :: d
  real(kind=mireal), dimension(2*n,2*n), intent(out) :: xi
  real(kind=mireal), intent(out) :: detlog
  !local variables
  real(kind=mireal), dimension(n,n) :: ai
  real(kind=mireal) :: detlog_a

  !Compute the inverse and determinant of the submatrix a
  call cholesky_inv_det(a,ai,detlog_a,n)
  !Compute the inverse and determinant of the whole matrix
  call invert_2x2block_matrix(b,d,ai,detlog_a,xi,detlog,n,n)

end subroutine invert_square_block_matrices

!This subroutine will compute the inverse of a matrix assumint it is a 2x2 block matrix
subroutine invert_2x2block_matrix(b,d,ai,detlog_a,xi,detlog,na,nd)
use constants
  integer, intent(in) :: na,nd
  real(kind=mireal), dimension(na,na), intent(in) :: ai
  real(kind=mireal), dimension(na,nd), intent(in) :: b
  real(kind=mireal), dimension(nd,nd), intent(in) :: d
  real(kind=mireal), intent(in) :: detlog_a !logarithmic value of the determinant of A
  real(kind=mireal), dimension(na+nd,na+nd), intent(out) :: xi
  real(kind=mireal), intent(out) :: detlog
  !local variables
  real(kind=mireal), dimension(nd,nd) :: dbtab, idbtab
  real(kind=mireal), dimension(nd,na) :: bt, btai
  real(kind=mireal), dimension(na,nd) :: aib
  real(kind=mireal) :: det

  bt = transpose(b)

  aib  = matmul(ai,b)
  btai = matmul(bt,ai)
  !Compute the D - B^T inv(A) B term
  dbtab = d - matmul(bt,matmul(ai,b))
  !invert the D - B^T inv(A) B term
  call cholesky_inv_det(dbtab,idbtab,det,nd)
  !Compute the determinant of the matrix
  !if X is a matrix ( A & B \\ B^T & D) then
  !Det(X) = Det(A) x Det( D - B^T inv(A) B )
  detlog = det + detlog_a

  !fill the matrices
  xi(1:na,1:na)             = ai + matmul(aib,matmul(idbtab,btai))
  xi(na+1:na+nd,1:na)       = - matmul(idbtab,btai)
  xi(1:na,na+1:na+nd)       = - matmul(aib,idbtab)
  xi(na+1:na+nd,na+1:na+nd) = idbtab

end subroutine invert_2x2block_matrix

!This subroutine will compute the inverse of a matrix assumint it is a 2x2 block matrix
subroutine invert_2x2block_matrix_new(b,d,ai,detlog_a,xi,detlog,na,nd)
use constants
  integer, intent(in) :: na,nd
  real(kind=mireal), dimension(na,na), intent(in) :: ai
  real(kind=mireal), dimension(na,nd), intent(in) :: b
  real(kind=mireal), dimension(nd,nd), intent(in) :: d
  real(kind=mireal), intent(in) :: detlog_a !logarithmic value of the determinant of A
  real(kind=mireal), dimension(na+nd,na+nd), intent(out) :: xi
  real(kind=mireal), intent(out) :: detlog
  !local variables
  real(kind=mireal), dimension(nd,nd) :: dbtab, idbtab
  real(kind=mireal), dimension(nd,na) :: bt, btai
  real(kind=mireal), dimension(na,nd) :: aib
  real(kind=mireal) :: det

  !Comparison between matmul and dgemm in
  !https://modelingguru.nasa.gov/docs/DOC-1762/diff

  aib  = matmul(ai,b)
  call DGEMM('n','n',na,nd,na,1.d0,ai,na,b,na,0.d0,aib,na)
  !print*,'1'
  btai = matmul(bt,ai)
  call DGEMM('n','n',nd,na,na,1.d0,bt,nd,ai,na,0.d0,btai,nd)
  !print*,'2'
  !Compute the D - B^T inv(A) B term
  dbtab = d - matmul(bt,aib)
  call DGEMM('n','n',nd,na,na,1.d0,bt,nd,aib,na,0.d0,dbtab,nd)
  !print*,'3'
  dbtab = d - dbtab
  !invert the D - B^T inv(A) B term
  call cholesky_inv_det(dbtab,idbtab,det,nd)
  !Compute the determinant of the matrix
  !if X is a matrix ( A & B \\ B^T & D) then
  !Det(X) = Det(A) x Det( D - B^T inv(A) B )
  detlog = det + detlog_a

  !fill the matrices
  xi(1:n,1:n)         = ai + matmul(aib,matmul(idbtab,btai))
  call DGEMM('n','n',nd,na,nd,1.d0,idbtab,nd,btai,nd,0.d0,xi(1:na,1:na),nd)
  !print*,'4'
  call DGEMM('n','n',na,na,nd,1.d0,aib,na,xi(1:na,1:na),nd,0.d0,xi(1:na,1:na),na)
  !print*,'5'
  xi(1:na,1:na) = ai + xi(1:na,1:na)
  xi(1:n,n+1:2*n)     = - matmul(aib,idbtab)
  call DGEMM('n','n',na,nd,nd,-1.d0,aib,na,idbtab,nd,0.d0,xi(na+1:na+nd,1:na),na)
  !print*,'6'
  xi(n+1:2*n,1:n)     = - matmul(idbtab,btai)
  call DGEMM('n','n',nd,na,nd,-1.d0,idbtab,nd,btai,nd,0.d0,xi(1:na,na+1:na+nd),nd)
  !print*,'7'
  xi(na+1:na+nd,na+1:na+nd) = idbtab
  !print*,'8'

end subroutine invert_2x2block_matrix_new

!Subroutine to invert the matrix for the multi-dimensional GP
!Input: Covariance matrix K, number of time-series nt
!Output: Inverse of the covariance matrix K, ln of matrix determinant K
!It follows !https://msvlab.hre.ntou.edu.tw/grades/now/inte/Inverse%20%26%20Border/border-LuTT.pdf
!computing iteratively the determinant of the matrix A, and the matrix D will always have dimension of
!the number of observations
subroutine invert_multi_gp_matrix(K,nt,Kinv,detlog,n)
   use constants
    integer, intent(in) :: n, nt
    real(kind=mireal), dimension(n,n), intent(in) :: K
    real(kind=mireal), dimension(n,n), intent(out) :: Kinv
    real(kind=mireal), intent(out) :: detlog
    !
    integer :: m, i, j
    real(kind=mireal)  :: detlog_a
    real(kind=mireal), dimension(n,n) :: Ainv

    !m is the number of data points per each time-series
    m = n/nt
    Kinv(:,:) = 0.d0
    Ainv(:,:) = 0.d0

    !if we only have one time-series then we compute it with a
    !cholesky method
    if (nt == 1) then
      call cholesky_inv_det(K,Kinv,detlog,n)
    !if we only have two time-series then we compute them as a 2x2 block matrix
    !inversion
    else if (nt == 2) then
      call invert_square_block_matrices(K(1:n/2,1:n/2),K(1:n/2,n/2+1:n),K(n/2+1:n,n/2+1:n),Kinv,detlog,n/2)
    !If we have more than two time-series, we need to peform an interive 2x2 block
!      call invert_square_block_matrices(K(1:n/2,1:n/2),K(n/2+1:n,1:n/2),K(n/2+1:n,n/2+1:n),Kinv,detlog,n/2)
    !matrix computation
    else if (nt > 2) then
      do i = 1, nt - 1
        j = i + 1
        if (i == 1) then
	  !This will compute the determinant of the first 2x2 Block, this will be the new matrix a and its determinant
          call invert_square_block_matrices &
            (K(1:m,1:m),K(1:m,m+1:2*m),K(m+1:2*m,m+1:2*m),Kinv(1:2*m,1:2*m),detlog,m)
            detlog_a = detlog
            Ainv = Kinv
          !Here Kinv(1:2*m,1:2*m) is the Ainv in the referred paper, we can use it  to compute the inverse interatively
        else
          !For i => 2, then we can use the subset of K and the inverse Kinv of the previous case as the new a and ai
          call invert_2x2block_matrix &
            (K(1:i*m,i*m+1:j*m),K(i*m+1:j*m,i*m+1:j*m),Ainv(1:i*m,1:i*m),detlog_a,Kinv(1:j*m,1:j*m),detlog,i*m,m)
            detlog_a = detlog
            Ainv = Kinv
          !This will create a new set of a and ai matrices that will be used in the next iteration
          !note that invert_block_matrices makes only one matrix inversion of a matrix of dimension mxm
        end if
      end do
      !When this do cycle ends, we have the inverse of the multi-GP matrix!
    end if

end subroutine invert_multi_gp_matrix