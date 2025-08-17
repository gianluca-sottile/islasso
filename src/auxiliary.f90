subroutine cov_from_hess_spd(p, hess, xtx, cov, hi, info)
  implicit none
  integer, intent(in)    :: p
  integer, intent(out)   :: info
  double precision, intent(in)  :: hess(p,p), xtx(p,p)
  double precision, intent(out) :: cov(p,p), hi(p)

  double precision :: R(p,p)
  integer :: k

  ! Copio H e fattorizzo Cholesky (UPLO='U': H = R^T R, con R triang. sup.)
  R = hess
  call DPOTRF('U', p, R, p, info)
  if (info /= 0) return

  ! Parto da XtX e faccio due solve "a sinistra" per ottenere H^{-1}*XtX
  cov = xtx
  call DTRSM('L','U','T','N', p, p, 1.0d0, R, p, cov, p)   ! cov := R^{-T} * cov
  call DTRSM('L','U','N','N', p, p, 1.0d0, R, p, cov, p)   ! cov := R^{-1} * cov = H^{-1}*XtX

  ! Salvo la diagonale dell’hat (come nel tuo codice: diag(H^{-1}*XtX))
  do k = 1, p
     hi(k) = cov(k,k)
  end do

  ! Ora post-moltiplico per H^{-1} = R^{-1} R^{-T} con due solve "a destra"
  call DTRSM('R','U','N','N', p, p, 1.0d0, R, p, cov, p)   ! cov := cov * R^{-1}
  call DTRSM('R','U','T','N', p, p, 1.0d0, R, p, cov, p)   ! cov := cov * R^{-T}

  ! Simmetrizzazione numerica (per eliminare errori d’arrotondamento)
  call symmetrize_upper(p, cov)
end subroutine cov_from_hess_spd

subroutine symmetrize_upper(n, A)
  implicit none
  integer, intent(in) :: n
  double precision, intent(inout) :: A(n,n)
  integer :: i, j
  do j = 1, n
     do i = j+1, n
        A(i,j) = A(j,i)
     end do
  end do
end subroutine symmetrize_upper

subroutine inv_posdef(n, a, ainv, info)
  implicit none
  integer, intent(in) :: n
  integer, intent(out) :: info
  double precision, intent(in) :: a(n,n)
  double precision, intent(out) :: ainv(n,n)

  ! Copia A in AINV
  ainv = a

  ! Fattorizzazione di Cholesky: A = L * L^T (UPLO = 'U' o 'L')
  call DPOTRF('U', n, ainv, n, info)
  if (info /= 0) return

  ! Inversione a partire dalla fattorizzazione di Cholesky
  call DPOTRI('U', n, ainv, n, info)
  if (info /= 0) return

  ! Ricostruzione della matrice completa simmetrica
  call copy_upper_to_lower(n, ainv)

end subroutine inv_posdef

!--------------------------------------------------------
subroutine copy_upper_to_lower(n, mat)
  implicit none
  integer, intent(in) :: n
  double precision, intent(inout) :: mat(n,n)
  integer :: i, j

  do j = 1, n
     do i = j+1, n
        mat(i,j) = mat(j,i)
     end do
  end do
end subroutine copy_upper_to_lower


! More efficient matrix inversion using LAPACK
subroutine inv_lapack(n, a, ainv, info, ipiv, work)
implicit none
integer :: n, info, ipiv(n)
double precision :: a(n,n), ainv(n,n), work(n*n)

! Make a copy of A in AINV
ainv = a

! LU decomposition
call DGETRF(n, n, ainv, n, ipiv, info)
if (info /= 0) return

! Matrix inversion
call DGETRI(n, ainv, n, ipiv, work, n*n, info)
end subroutine inv_lapack

subroutine standardize(X, xm, xse, n, p, intercept)
implicit none
integer :: n, p, intercept
double precision :: X(n,p), xm(p), xse(p)
integer :: j, start
double precision :: sum_x, sum_x2

xm = 0.d0
xse = 1.d0
start = 1
if((intercept == 1) .and. (p > 1)) start = 2

do j = start, p
    sum_x = sum(X(:, j))
    sum_x2 = dot_product(X(:, j), X(:, j))

    xm(j) = sum_x / n
    xse(j) = sqrt(sum_x2 / n - xm(j)**2)

    ! Avoid division by zero
    if(xse(j) < 1.0d-12) xse(j) = 1.0d0

    X(:, j) = (X(:, j) - xm(j)) / xse(j)
end do
end subroutine standardize

subroutine check_out(theta, cov, xm, xse, p, intercept)
    implicit none
    integer, intent(in) :: p, intercept
    double precision, intent(inout) :: theta(p), cov(p,p)
    double precision, intent(in) :: xm(p), xse(p)
    integer :: i, j
    double precision :: temp_vector(p)

    ! Vectorized division
    theta = theta / xse

    ! Handle intercept adjustment for theta
    if ((intercept == 1) .and. (p > 1)) then
        theta(1) = theta(1) - dot_product(theta(2:p), xm(2:p))
    end if

    ! Optimize the covariance matrix calculation
    ! First handle diagonal elements
    do i = 1, p
        cov(i,i) = cov(i,i) / (xse(i) * xse(i))
    end do

    ! Handle off-diagonal elements with a single loop
    ! This avoids redundant assignments and improves cache locality
    do i = 1, p-1
        do j = i+1, p
            cov(i,j) = cov(i,j) / (xse(i) * xse(j))
            cov(j,i) = cov(i,j)  ! Mirror value to avoid recalculation
        end do
    end do

    ! Handle intercept adjustment for covariance matrix
    if ((intercept == 1) .and. (p > 1)) then
        ! Store the result of matrix multiplication to avoid redundant calculation
        temp_vector = matmul(xm(2:p), cov(2:p,:))

        ! Update first row
        cov(1,:) = cov(1,:) - temp_vector

        ! Update first column (mirror of first row)
        cov(:,1) = cov(1,:)

        ! Update the (1,1) element
        cov(1,1) = cov(1,1) - dot_product(cov(1,2:p), xm(2:p))
    end if
end subroutine check_out
