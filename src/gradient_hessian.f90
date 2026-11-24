subroutine hessian_theta(theta, se, lambda, xtx, pi, p, hess, alpha)
implicit none
integer :: p, i
double precision, intent(in) :: theta(p), se(p), lambda(p), xtx(p,p), pi(p), alpha
double precision, intent(out) :: hess(p,p)
double precision :: theta2(p), temp1, val_pi, term_pnm, pnm
double precision, parameter :: threshold = 1.0d-6

! Copy xtx to hess directly
hess = xtx

do i = 1, p
    theta2(i) = theta(i)
    if(abs(theta2(i)) < threshold) theta2(i) = threshold

    temp1 = theta2(i) / se(i)
    val_pi = pi(i)

    term_pnm = val_pi * (2.d0 * pnm(temp1, 0.d0, 1.d0) - 1.d0) + &
                (1.d0 - val_pi) * (2.d0 * pnm(temp1, 0.d0, 1.d-5) - 1.d0)

    hess(i,i) = hess(i,i) + lambda(i) * (alpha * term_pnm / theta2(i) + (1.d0 - alpha))
end do
end subroutine hessian_theta

subroutine hessian(theta, se, lambda, xtx, pi, p, hess, alpha)
implicit none
integer :: p, i
double precision, intent(in) :: theta(p), se(p), lambda(p), xtx(p,p), pi(p), alpha
double precision, intent(out) :: hess(p,p)
double precision :: temp1, val_pi, term_dnm, dnm

! Copy xtx to hess directly
hess = xtx

do i = 1, p
    temp1 = theta(i) / se(i)
    val_pi = pi(i)

    term_dnm = val_pi * dnm(temp1, 0.d0, 1.d0) + &
                (1.d0 - val_pi) * dnm(temp1, 0.d0, 1.d-5)

    hess(i,i) = hess(i,i) + lambda(i) * (2.d0 * alpha * term_dnm / se(i) + &
                 (1.d0 - alpha))
end do
end subroutine hessian

subroutine gradient(theta, se, lambda, xtw, res, pi, n, p, grad, alpha)
implicit none
integer :: n, p, i
double precision, intent(in) :: theta(p), se(p), lambda(p), xtw(p,n), res(n), pi(p), alpha
double precision, intent(out) :: grad(p)
double precision :: temp1, val_pi, term_pnm, pnm

! Compute - xtw * res efficiently with BLAS
call DGEMV('N', p, n, -1.d0, xtw, p, res, 1, 0.d0, grad, 1)
! Note: using -1.d0 scalar to avoid the need for the negation after DGEMV

do i = 1, p
    temp1 = theta(i) / se(i)
    val_pi = pi(i)

    term_pnm = val_pi * (2.d0 * pnm(temp1, 0.d0, 1.d0) - 1.d0) + &
                (1.d0 - val_pi) * (2.d0 * pnm(temp1, 0.d0, 1.d-5) - 1.d0)

    grad(i) = grad(i) + lambda(i) * (alpha * term_pnm + (1.d0 - alpha) * theta(i))
end do
end subroutine gradient

subroutine penalty(theta, se, pi, p, pen, alpha)
implicit none
integer :: p, i
double precision, intent(in) :: theta(p), se(p), pi(p), alpha
double precision, intent(out) :: pen(p)
double precision :: temp1, val_pi, term_pnm, pnm

do i = 1, p
    temp1 = theta(i) / se(i)
    val_pi = pi(i)

    term_pnm = val_pi * (2.d0 * pnm(temp1, 0.d0, 1.d0) - 1.d0) + &
               (1.d0 - val_pi) * (2.d0 * pnm(temp1, 0.d0, 1.d-5) - 1.d0)

    pen(i) = alpha * term_pnm + (1.d0 - alpha) * theta(i)
end do
end subroutine penalty
