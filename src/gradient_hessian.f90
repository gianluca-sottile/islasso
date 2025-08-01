subroutine hessian_theta(theta, se, lambda, xtx, pi, p, hess, alpha)
implicit none
integer :: p, i
double precision :: theta(p), se(p), lambda(p), xtx(p,p), pi(p), hess(p,p)
double precision :: theta2(p), pnm, temp1, alpha, threshold

! Copy xtx to hess directly
hess = xtx
threshold = 1.0d-6

do i = 1, p
    theta2(i) = theta(i)
    if(abs(theta2(i)) < threshold) theta2(i) = threshold

    temp1 = theta2(i) / se(i)
    hess(i,i) = hess(i,i) + lambda(i) * alpha * ( pi(i) * (2.d0 * pnm(temp1, 0.d0, 1.d0) - 1.d0) + &
        & (1.d0 - pi(i)) * (2.d0 * pnm(temp1, 0.d0, 0.00001d0) - 1.d0) ) / theta2(i) + &
        & lambda(i) * (1.d0 - alpha)
end do
end subroutine hessian_theta

subroutine hessian(theta, se, lambda, xtx, pi, p, hess, alpha)
implicit none
integer :: p, i
double precision :: theta(p), se(p), lambda(p), xtx(p,p), pi(p), hess(p,p)
double precision :: dnm, temp1, alpha

! Copy xtx to hess directly
hess = xtx

do i = 1, p
    temp1 = theta(i) / se(i)
    hess(i,i) = hess(i,i) + 2.d0 * lambda(i) * alpha * ( pi(i) * dnm(temp1, 0.d0, 1.d0) + &
        & (1.d0 - pi(i)) * dnm(temp1, 0.d0, 0.00001d0) ) / se(i) + (1.d0 - alpha) * lambda(i)
end do
end subroutine hessian

subroutine gradient(theta, se, lambda, xtw, res, pi, n, p, grad, alpha)
implicit none
integer :: n, p, i
double precision :: theta(p), se(p), lambda(p), xtw(p,n), res(n), pi(p), grad(p)
double precision :: pnm, temp1, alpha

! Compute xtw * res efficiently with BLAS
call DGEMV('N', p, n, -1.d0, xtw, p, res, 1, 0.d0, grad, 1)
! Note: using -1.d0 scalar to avoid the need for the negation after DGEMV

do i = 1, p
    temp1 = theta(i) / se(i)
    grad(i) = grad(i) + lambda(i) * alpha * ( &
              pi(i) * (2.d0 * pnm(temp1, 0.d0, 1.d0) - 1.d0) + &
              (1.d0 - pi(i)) * (2.d0 * pnm(temp1, 0.d0, 0.00001d0) - 1.d0) &
              ) + lambda(i) * (1.d0 - alpha) * theta(i)
end do
end subroutine gradient

subroutine penalty(theta, se, pi, p, pen, alpha)
implicit none
integer :: p, i
double precision :: theta(p), se(p), pi(p), pen(p)
double precision :: pnm, alpha, temp1

do i = 1, p
    temp1 = theta(i) / se(i)
    pen(i) = alpha * ( pi(i) * (2.d0 * pnm(temp1, 0.d0, 1.d0) - 1.d0) + &
              (1.d0 - pi(i)) * (2.d0 * pnm(temp1, 0.d0, 0.00001d0) - 1.d0) &
              ) + (1.d0 - alpha) * theta(i)
end do
end subroutine penalty
