subroutine islasso(X, y, n, p, theta, se, cov, lambda, alpha, pi, estpi, &
& itmax, itmaxse, tol, sigma2, trace, adaptive, offset, conv, stand, intercept, &
& eta, mu, varmu, mu_eta_val, w, res, dev, weights, hi, edf, xtw, xtx, grad, hess, &
& invH, pen)

implicit none

!external variables
integer :: n, p, estpi, itmax, itmaxse, trace, conv, adaptive, stand, intercept
double precision :: X(n,p), y(n), theta(p), se(p), cov(p,p), lambda(p), alpha
double precision :: pi(p), tol, sigma2, offset(n), eta(n), mu(n), varmu(n)
double precision :: mu_eta_val(n), w(n), res(n), dev, weights(n), hi(p), edf
double precision :: xtw(p,n), xtx(p,p), grad(p), hess(p,p), invH(p,p), pen(p)

! internal variables
integer :: i, j, k, info
double precision :: X_orig(n,p), xm(p), xse(p), xtwy(p), lmbd0, lambdaadapt(p)
double precision :: theta0(p), cov0(p,p), se0(p), s2, hat_matrix(p,p), redf, ind, ind2

! Add workspace arrays for LAPACK
double precision, allocatable :: work(:)
integer, allocatable :: ipiv(:)

! Allocate workspace for LAPACK routines
allocate(ipiv(p))
allocate(work(p*p))

info = 0

! Initialize variables
s2 = sigma2
conv = 0
varmu = 1.d0
mu_eta_val = 1.d0
w = 1.d0
hi = 1.d0

xm = 0.d0
xse = 1.d0

lmbd0 = MAXVAL(lambda)

X_orig = X
! standardizing
if(stand.eq.1) then
    call standardize(X, xm, xse, n, p, intercept)
    lambda = lambda / xse
end if
lambdaadapt = lambda

! Pre-compute X*weights once before the loop
do i = 1, p
    xtw(i,:) = X(:,i) * weights
end do

! Compute X'WX and X'Wy more efficiently
call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
call DGEMV('N', p, n, 1.d0, xtw, p, y - offset, 1, 0.d0, xtwy, 1)

cov0 = cov
se0 = se

! main body
do i = 1, itmaxse
    call rchkusr()
    if(trace.eq.9) call islasso_trace2_3(i)

    if((adaptive.eq.1).and.(i.gt.5)) then
      do k = 1, p
        hi(k) = max(min(hi(k), 1.0d0), 0.00001d0)
        lambdaadapt(k) = lambda(k) * (1 - hi(k)) / hi(k)
      end do
    end if

    theta0 = theta
    ! computing mixture parameters c
    if((estpi.eq.1).and.(i.gt.5)) then
        call logitlinkinv(abs(theta0 / se0), p, pi)
        pi = 0.75d0 * (2.d0 * pi - 1.d0) + 0.25d0
    end if

    do j = 1, itmax
        ! computing IWLS - use optimized versions of fn1
        call hessian_theta(theta0, se0, lambdaadapt, xtx, pi, p, hess, alpha)

        ! Replace solve with direct LAPACK call for better performance
        call DCOPY(p, xtwy, 1, work, 1)  ! Copy xtwy to work
        !call DGESV(p, 1, hess, p, ipiv, work, p, info)

        call DPOTRF('U', p, hess, p, info)
        call DPOTRS('U', p, 1, hess, p, work, p, info)

        if(info.ne.0) then
            conv = 2
            exit
        end if

        call DCOPY(p, work, 1, theta, 1)  ! Copy result back to theta
        theta = theta0 + 0.5d0 * (theta - theta0)

        ind2 = MAXVAL(abs(theta - theta0))
        if(ind2.le.tol) then
            exit
        end if

        theta0 = theta
    end do
    if((conv.eq.2).or.(itmaxse.eq.0)) exit

    ! Use more efficient matrix-vector operation
    call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
    eta = eta + offset
    mu = eta
    res = y - mu

    dev = sum(weights * (res**2))

    ! updating components for variance covariance matrix
    call hessian(theta, se0, lambdaadapt, xtx, pi, p, hess, alpha)

    call inv_posdef(p, hess, invH, info, ipiv, work)

    !! More efficient matrix inversion using LAPACK
    !call inv_lapack(p, hess, invH, info, ipiv, work)
    if(info.ne.0) then
        conv = 2
        exit
    end if

    call DGEMM('N', 'N', p, p, p, 1.d0, invH, p, xtx, p, 0.d0, hat_matrix, p)
    call DGEMM('N', 'N', p, p, p, 1.d0, hat_matrix, p, invH, p, 0.d0, cov, p)
    cov = cov0 + 0.1d0 * (cov - cov0)

    do k = 1, p
      hi(k) = hat_matrix(k,k)
    end do

    edf = sum(hi)
    redf = n - edf
    if (redf .le. 1.0d-12) redf = 1.0d0
    if(sigma2.le.0.0d0) s2 = dev / redf

    do k = 1, p
        se(k) = sqrt(s2 * cov(k,k))
    end do

    ! checking possible convergence criterion
    ind = MAXVAL(abs(se - se0))
    if(trace.eq.2) call islasso_trace1_2_2(tol, i, lmbd0, dev, redf, s2, ind, ind2)
    if(trace.eq.1) call islasso_trace1_7_2(tol, i, lmbd0, dev, redf, s2, ind, ind2)

    if(ind.le.(tol*10.d0)) then
        if((trace.eq.1).or.(trace.eq.2)) call islasso_trace1_8(1)
        exit
    end if

    ! conv = 1 if i >= itmax
    if(i.ge.itmaxse) then
        conv = 1
        exit
    end if

    cov0 = cov
    se0 = se
end do

! if standardized beta and se, then return to the original scale
if(stand.eq.1) then
    call check_out(theta, cov, xm, xse, p, intercept)
    lambdaadapt = lambdaadapt * xse

    do k = 1, p
        se(k) = sqrt(s2 * cov(k,k))
    end do

    X = X_orig
    do i = 1, p
        xtw(i,:) = X(:,i) * weights
    end do
    call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)

    call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
    eta = eta + offset
    mu = eta
    res = y - mu
    dev = sum(weights * (res**2))
end if

! updating components for variance covariance matrix
call gradient(theta, se, lambdaadapt, xtw, res, pi, n, p, grad, alpha)
call hessian(theta, se, lambdaadapt, xtx, pi, p, hess, alpha)
!call inv_lapack(p, hess, invH, info, ipiv, work)
call inv_posdef(p, hess, invH, info, ipiv, work)
if(info.ne.0) then
    conv = 2
end if
call penalty(theta, se, pi, p, pen, alpha)

! updating output components
itmax = i
lambda = lambdaadapt
tol = ind
sigma2 = s2

! Clean up
deallocate(ipiv)
deallocate(work)

end subroutine islasso

subroutine islasso_glm(X, y, n, p, theta, se, cov, lambda, alpha, pi, estpi, &
& itmax, itmaxse, tol, sigma2, trace, adaptive, offset, conv, stand, intercept, &
& eta, mu, varmu, mu_eta_val, w, res, dev, weights, hi, edf, xtw, xtx, grad, hess, &
& invH, pen, fam, link)

implicit none

!external variables
integer :: n, p, estpi, itmax, itmaxse, trace, conv, adaptive, stand, intercept, fam, link
double precision :: X(n,p), y(n), theta(p), se(p), cov(p,p), lambda(p), alpha
double precision :: pi(p), tol, sigma2, offset(n), eta(n), mu(n), varmu(n)
double precision :: mu_eta_val(n), w(n), res(n), dev, weights(n), hi(p), edf
double precision :: xtw(p,n), xtx(p,p), grad(p), hess(p,p), invH(p,p), pen(p)

! internal variables
integer :: i, j, k, info
double precision :: X_orig(n,p), xm(p), xse(p), xtwz(p), z(n), lmbd0, lambdaadapt(p)
double precision :: theta0(p), cov0(p,p), se0(p), s2, hat_matrix(p,p), redf, ind, ind2

! Add workspace arrays for LAPACK
integer, allocatable :: ipiv(:)
double precision, allocatable :: work(:)

! Allocate workspace for LAPACK routines
allocate(ipiv(p))
allocate(work(p*p))

info = 0

! Initialize variables
s2 = sigma2
conv = 0
hi = 1.d0

xm = 0.d0
xse = 1.d0

lmbd0 = MAXVAL(lambda)

X_orig = X
! standardizing
if(stand.eq.1) then
    call standardize(X, xm, xse, n, p, intercept)
    lambda = lambda / xse
end if
lambdaadapt = lambda

! Initial calculation of eta, mu, varmu, mu_eta_val
call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
eta = eta + offset
call family(fam, link, 2, eta, n, mu)           ! link inverse
call family(fam, link, 4, mu, n, varmu)         ! var(mu)
call family(fam, link, 3, eta, n, mu_eta_val)   ! dmu/deta
res = (y - mu) / mu_eta_val
w = weights * (mu_eta_val**2) / varmu
! Compute X'WX and X'Wz efficiently
do k = 1, p
    xtw(k,:) = X(:,k) * w
end do
call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
z = (eta - offset) + res
call DGEMV('N', p, n, 1.d0, xtw, p, z, 1, 0.d0, xtwz, 1)

cov0 = cov
se0 = se

do i = 1, itmaxse
    call rchkusr()

    ! Update adaptive lambda if enabled
    if((adaptive.eq.1) .and. (i.gt.5)) then
        do k = 1, p
            hi(k) = max(min(hi(k), 1.0d0), 0.00001d0)
            lambdaadapt(k) = lambda(k) * (1 - hi(k)) / hi(k)
        end do
    end if

    theta0 = theta
    ! Update mixture parameters if needed
    if((estpi.eq.1) .and. (i.gt.5)) then
        call logitlinkinv(abs(theta0 / se0), p, pi)
        pi = 0.75d0 * (2.d0 * pi - 1.d0) + 0.25d0
    end if

    ! Inner IWLS iteration
    do j = 1, itmax

        ! Computing IWLS
        call hessian_theta(theta0, se0, lambdaadapt, xtx, pi, p, hess, alpha)

        ! Solve linear system using LAPACK
        call DCOPY(p, xtwz, 1, work, 1)
        !call DGESV(p, 1, hess, p, ipiv, work, p, info)

        call DPOTRF('U', p, hess, p, info)
        call DPOTRS('U', p, 1, hess, p, work, p, info)

        if(info /= 0) then
            conv = 2
            exit
        end if
        call DCOPY(p, work, 1, theta, 1)
        theta = theta0 + 0.5d0 * (theta - theta0)

        ! Check convergence
        ind2 = MAXVAL(abs(theta - theta0))
        if(ind2.le.tol) then
            exit
        end if

        theta0 = theta
    end do
    if((conv == 2) .or. (itmaxse == 0)) exit

    ! Update linear predictor and other variables
    call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
    eta = eta + offset
    call family(fam, link, 2, eta, n, mu)
    call family(fam, link, 4, mu, n, varmu)
    call family(fam, link, 3, eta, n, mu_eta_val)
    res = (y - mu) / mu_eta_val
    w = weights * (mu_eta_val**2) / varmu
    do k = 1, p
        xtw(k,:) = X(:,k) * w
    end do
    call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)
    z = (eta - offset) + res
    call DGEMV('N', p, n, 1.d0, xtw, p, z, 1, 0.d0, xtwz, 1)

    dev = sum(w * (res**2))

    ! Update covariance matrix
    call hessian(theta, se0, lambdaadapt, xtx, pi, p, hess, alpha)
    !call inv_lapack(p, hess, invH, info, ipiv, work)
    call inv_posdef(p, hess, invH, info, ipiv, work)
    if(info.ne.0) then
        conv = 2
        exit
    end if

    call DGEMM('N', 'N', p, p, p, 1.d0, invH, p, xtx, p, 0.d0, hat_matrix, p)
    call DGEMM('N', 'N', p, p, p, 1.d0, hat_matrix, p, invH, p, 0.d0, cov, p)
    cov = cov0 + 0.1d0 * (cov - cov0)

    ! Extract hat matrix diagonal
    do k = 1, p
        hi(k) = hat_matrix(k,k)
    end do

    edf = sum(hi)
    redf = n - edf
    if(sigma2.le.0) s2 = dev / redf

    do k = 1, p
        se(k) = sqrt(s2 * cov(k,k))
    end do

    ! Check convergence of standard errors
    ind = MAXVAL(abs(se - se0))
    if(trace == 2) call islasso_trace2_2_2(tol, i, lmbd0, dev, redf, s2, ind, ind2)
    if(trace == 1) call islasso_trace2_7_2(tol, i, lmbd0, dev, redf, s2, ind, ind2)

    if(ind.le.(tol*10)) then
        if((trace == 1) .or. (trace == 2)) call islasso_trace1_8(1)
        if(trace == 9) call islasso_trace2_6(i)
        exit
    end if

    ! Check maximum iterations
    if(i >= itmaxse) then
        conv = 1
        exit
    end if

    cov0 = cov
    se0 = se
end do

! if standardized beta and se, then return to the original scale
if(stand.eq.1) then
    call check_out(theta, cov, xm, xse, p, intercept)
    lambdaadapt = lambdaadapt * xse

    do k = 1, p
        se(k) = sqrt(s2 * cov(k,k))
    end do

    X = X_orig
    call DGEMV('N', n, p, 1.d0, X, n, theta, 1, 0.d0, eta, 1)
    eta = eta + offset
    call family(fam, link, 2, eta, n, mu)
    call family(fam, link, 4, mu, n, varmu)
    call family(fam, link, 3, eta, n, mu_eta_val)
    res = (y - mu) / mu_eta_val
    w = weights * (mu_eta_val**2) / varmu
    do k = 1, p
        xtw(k,:) = X(:,k) * w
    end do
    call DGEMM('N', 'N', p, p, n, 1.d0, xtw, p, X, n, 0.d0, xtx, p)

    dev = sum(w * (res**2))
end if

! updating components for variance covariance matrix
call gradient(theta, se, lambdaadapt, xtw, res, pi, n, p, grad, alpha)
call hessian(theta, se, lambdaadapt, xtx, pi, p, hess, alpha)
!call inv_lapack(p, hess, invH, info, ipiv, work)
call inv_posdef(p, hess, invH, info, ipiv, work)
if(info.ne.0) then
    conv = 2
end if
call penalty(theta, se, pi, p, pen, alpha)

! updating output components
itmax = i
lambda = lambdaadapt
tol = ind
sigma2 = s2

! Clean up
deallocate(ipiv)
deallocate(work)

end subroutine islasso_glm
