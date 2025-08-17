!===============================
! FAMILY DISPATCHER
!===============================
subroutine family(fam, link, func, x, n, y)
    implicit none
    integer, intent(in) :: fam, link, func, n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: y(n)

    ! fam: 1=Binomial, 2=Poisson, 3=Gamma
    ! link: dipende dalla famiglia (vedi sotto)
    ! func: 1=link, 2=linkinv, 3=mueta, 4=variance

    y = 0.0d0

    select case (fam)

    case (1)   ! Binomial
        select case (link)
        case (1)   ! logit
            select case (func)
            case (1); call logitlink(x, n, y)
            case (2); call logitlinkinv(x, n, y)
            case (3); call logitmueta(x, n, y)
            case (4); call binomial_variance(x, n, y)
            end select
        case (2)   ! probit
            select case (func)
            case (1); call probitlink(x, n, y)
            case (2); call probitlinkinv(x, n, y)
            case (3); call probitmueta(x, n, y)
            case (4); call binomial_variance(x, n, y)
            end select
        end select

    case (2)   ! Poisson
        select case (link)
        case (1)   ! log
            select case (func)
            case (1); call loglink(x, n, y)
            case (2); call loglinkinv(x, n, y)
            case (3); call logmueta(x, n, y)
            case (4); call poisson_variance(x, n, y)
            end select
        end select

    case (3)   ! Gamma
        select case (link)
        case (1)   ! inverse
            select case (func)
            case (1); call inverselink(x, n, y)
            case (2); call inverselinkinv(x, n, y)
            case (3); call inversemueta(x, n, y)
            case (4); call gamma_variance(x, n, y)
            end select
        case (2)   ! log
            select case (func)
            case (1); call loglink(x, n, y)
            case (2); call loglinkinv(x, n, y)
            case (3); call logmueta(x, n, y)
            case (4); call gamma_variance(x, n, y)
            end select
        case (3)   ! identity
            select case (func)
            case (1); call identitylink(x, n, y)
            case (2); call identitylinkinv(x, n, y)
            case (3); call identitymueta(x, n, y)
            case (4); call gamma_variance(x, n, y)
            end select
        end select
    end select
end subroutine family

!===============================
! BINOMIAL
!===============================
subroutine binomial_variance(x, n, varmu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: varmu(n)
    ! var(mu) = mu*(1-mu)
    varmu = x * (1.0d0 - x)
end subroutine binomial_variance

subroutine logitlink(x, n, mu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)   ! domain: x in (0,1)
    double precision, intent(out) :: mu(n)
    ! mu = log(x/(1-x))
    mu = log( x / (1.0d0 - x) )
end subroutine logitlink

subroutine logitlinkinv(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: eta(n)
    double precision, parameter :: thresh = 2.220446049250313D-16
    double precision, parameter :: invthresh = 4503599627370496.d0
    double precision :: tmp(n)

    ! tmp = exp(x) con clamp ai tuoi cut-off (-30,+30) e costanti fornite
    tmp = exp(x)
    where (x < -30.0d0) tmp = thresh
    where (x >  30.0d0) tmp = invthresh

    eta = tmp / (1.0d0 + tmp)
end subroutine logitlinkinv

subroutine logitmueta(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: eta(n)
    double precision, parameter :: thresh = 2.220446049250313D-16
    double precision :: ex(n)

    ex = exp(x)
    eta = ex / ( (1.0d0 + ex)**2 )
    where ((x < -30.0d0) .or. (x > 30.0d0)) eta = thresh
end subroutine logitmueta

subroutine probitlink(x, n, mu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: mu(n)
    double precision, external :: qnm
    integer :: i
    do i = 1, n
        mu(i) = qnm(x(i), 0.0d0, 1.0d0)
    end do
end subroutine probitlink

subroutine probitlinkinv(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: eta(n)
    double precision, external :: pnm
    double precision, parameter :: thresh = 8.12589066470190d6
    integer :: i
    double precision :: val
    do i = 1, n
        val = x(i)
        if (val <= -thresh) val = -thresh
        if (val >=  thresh) val =  thresh
        eta(i) = pnm(val, 0.0d0, 1.0d0)
    end do
end subroutine probitlinkinv

subroutine probitmueta(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: eta(n)
    double precision, external :: dnm
    double precision, parameter :: thresh = 2.220446049250313D-16
    integer :: i
    do i = 1, n
        eta(i) = dnm(x(i), 0.0d0, 1.0d0)
        if (eta(i) <= thresh) eta(i) = thresh
    end do
end subroutine probitmueta

!===============================
! POISSON
!===============================
subroutine poisson_variance(x, n, varmu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: varmu(n)
    varmu = x
end subroutine poisson_variance

subroutine loglink(x, n, mu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)   ! domain: x>0
    double precision, intent(out) :: mu(n)
    mu = log(x)
end subroutine loglink

subroutine loglinkinv(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: eta(n)
    double precision, parameter :: thresh = 2.220446049250313D-16
    eta = exp(x)
    where (eta <= thresh) eta = thresh
end subroutine loglinkinv

subroutine logmueta(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: eta(n)
    double precision, parameter :: thresh = 2.220446049250313D-16
    eta = exp(x)
    where (eta <= thresh) eta = thresh
end subroutine logmueta

!===============================
! GAMMA
!===============================
subroutine gamma_variance(x, n, varmu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: varmu(n)
    varmu = x**2
end subroutine gamma_variance

subroutine inverselink(x, n, mu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: mu(n)
    mu = 1.0d0 / x
end subroutine inverselink

subroutine inverselinkinv(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: eta(n)
    eta = 1.0d0 / x
end subroutine inverselinkinv

subroutine inversemueta(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: eta(n)
    eta = -1.0d0 / (x**2)
end subroutine inversemueta

subroutine identitylink(x, n, mu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: mu(n)
    mu = x
end subroutine identitylink

subroutine identitylinkinv(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: eta(n)
    eta = x
end subroutine identitylinkinv

subroutine identitymueta(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in)  :: x(n)
    double precision, intent(out) :: eta(n)
    eta = 1.0d0
end subroutine identitymueta
