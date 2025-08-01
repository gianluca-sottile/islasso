subroutine family(fam, link, func, x, n, y)
    implicit none
    integer, intent(in) :: fam, link, func, n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: y(n)

    ! Use direct calls based on input parameters instead of nested select cases
    if (fam == 1) then       ! Binomial
        if (link == 1) then  ! logit
            if (func == 1) then
                call logitlink(x, n, y)
            else if (func == 2) then
                call logitlinkinv(x, n, y)
            else if (func == 3) then
                call logitmueta(x, n, y)
            else if (func == 4) then
                call binomial_variance(x, n, y)
            end if
        else if (link == 2) then  ! probit
            if (func == 1) then
                call probitlink(x, n, y)
            else if (func == 2) then
                call probitlinkinv(x, n, y)
            else if (func == 3) then
                call probitmueta(x, n, y)
            else if (func == 4) then
                call binomial_variance(x, n, y)
            end if
        end if
    else if (fam == 2) then  ! Poisson
        if (link == 1) then  ! log
            if (func == 1) then
                call loglink(x, n, y)
            else if (func == 2) then
                call loglinkinv(x, n, y)
            else if (func == 3) then
                call logmueta(x, n, y)
            else if (func == 4) then
                call poisson_variance(x, n, y)
            end if
        end if
    else if (fam == 3) then  ! Gamma
        if (link == 1) then  ! inverse
            if (func == 1) then
                call inverselink(x, n, y)
            else if (func == 2) then
                call inverselinkinv(x, n, y)
            else if (func == 3) then
                call inversemueta(x, n, y)
            else if (func == 4) then
                call gamma_variance(x, n, y)
            end if
        else if (link == 2) then  ! log
            if (func == 1) then
                call loglink(x, n, y)
            else if (func == 2) then
                call loglinkinv(x, n, y)
            else if (func == 3) then
                call logmueta(x, n, y)
            else if (func == 4) then
                call gamma_variance(x, n, y)
            end if
        else if (link == 3) then  ! identity
            if (func == 1) then
                call identitylink(x, n, y)
            else if (func == 2) then
                call identitylinkinv(x, n, y)
            else if (func == 3) then
                call identitymueta(x, n, y)
            else if (func == 4) then
                call gamma_variance(x, n, y)
            end if
        end if
    end if
end subroutine family

! BINOMIAL
subroutine binomial_variance(x, n, varmu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: varmu(n)
    
    ! Vectorized operation
    varmu = x * (1.d0 - x)
end subroutine binomial_variance

subroutine logitlink(x, n, mu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: mu(n)
    
    ! Vectorized operation
    mu = log(x / (1.d0 - x))
end subroutine logitlink

subroutine logitlinkinv(x, n, eta)
    implicit none
    integer, intent(in) :: n
    integer :: i
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: eta(n)
    double precision :: tmp
    double precision, parameter :: thresh = 2.220446049250313D-16
    double precision, parameter :: invthresh = 4503599627370496.d0
    
    ! Loop optimization with branch prediction hints
    do i = 1, n
        if (x(i) < -30.d0) then
            tmp = thresh
        else if (x(i) > 30.d0) then
            tmp = invthresh
        else
            tmp = exp(x(i))
        end if
        eta(i) = tmp / (1.d0 + tmp)
    end do
end subroutine logitlinkinv

subroutine logitmueta(x, n, eta)
    implicit none
    integer, intent(in) :: n
    integer :: i
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: eta(n)
    double precision, parameter :: thresh = 2.220446049250313D-16
    double precision :: expx
    
    ! Loop optimization with temporary variable to avoid redundant calculations
    do i = 1, n
        if ((x(i) < -30.d0) .or. (x(i) > 30.d0)) then
            eta(i) = thresh
        else
            expx = exp(x(i))
            eta(i) = expx / ((1.d0 + expx)**2)
        end if
    end do
end subroutine logitmueta

subroutine probitlink(x, n, mu)
    implicit none
    integer, intent(in) :: n
    integer :: i
    double precision, intent(in) :: x(n)
    double precision, external :: qnm
    double precision, intent(out) :: mu(n)
    
    do i = 1, n
        mu(i) = qnm(x(i), 0.d0, 1.d0)
    end do
end subroutine probitlink

subroutine probitlinkinv(x, n, eta)
    implicit none
    integer, intent(in) :: n
    integer :: i
    double precision, intent(in) :: x(n)
    double precision, external :: pnm
    double precision, intent(out) :: eta(n)
    double precision, parameter :: thresh = 8.12589066470190d6
    double precision :: val
    
    ! Loop optimization with temporary variable
    do i = 1, n
        val = x(i)
        if (val <= -thresh) val = -thresh
        if (val >= thresh) val = thresh
        eta(i) = pnm(val, 0.d0, 1.d0)
    end do
end subroutine probitlinkinv

subroutine probitmueta(x, n, eta)
    implicit none
    integer, intent(in) :: n
    integer :: i
    double precision, intent(in) :: x(n)
    double precision, external :: dnm
    double precision, intent(out) :: eta(n)
    double precision, parameter :: thresh = 2.220446049250313D-16
    
    do i = 1, n
        eta(i) = dnm(x(i), 0.d0, 1.d0)
        if (eta(i) <= thresh) eta(i) = thresh
    end do
end subroutine probitmueta

! POISSON
subroutine poisson_variance(x, n, varmu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: varmu(n)
    
    ! Vectorized operation
    varmu = x
end subroutine poisson_variance

subroutine loglink(x, n, mu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: mu(n)
    
    ! Vectorized operation
    mu = log(x)
end subroutine loglink

subroutine loglinkinv(x, n, eta)
    implicit none
    integer, intent(in) :: n
    integer :: i
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: eta(n)
    double precision, parameter :: thresh = 2.220446049250313D-16
    
    ! Two-phase approach: vectorized operation first, then fix threshold values
    eta = exp(x)
    
    ! Only loop through elements that need checking
    do i = 1, n
        if (eta(i) <= thresh) eta(i) = thresh
    end do
end subroutine loglinkinv

subroutine logmueta(x, n, eta)
    implicit none
    integer, intent(in) :: n
    integer :: i
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: eta(n)
    double precision, parameter :: thresh = 2.220446049250313D-16
    
    ! Two-phase approach: vectorized operation first, then fix threshold values
    eta = exp(x)
    
    ! Only loop through elements that need checking
    do i = 1, n
        if (eta(i) <= thresh) eta(i) = thresh
    end do
end subroutine logmueta

! GAMMA
subroutine gamma_variance(x, n, varmu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: varmu(n)
    
    ! Vectorized operation
    varmu = x**2
end subroutine gamma_variance

subroutine inverselink(x, n, mu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: mu(n)
    
    ! Vectorized operation
    mu = 1.d0 / x
end subroutine inverselink

subroutine inverselinkinv(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: eta(n)
    
    ! Vectorized operation
    eta = 1.d0 / x
end subroutine inverselinkinv

subroutine inversemueta(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: eta(n)
    
    ! Vectorized operation
    eta = -1.d0 / (x**2)
end subroutine inversemueta

subroutine identitylink(x, n, mu)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: mu(n)
    
    ! Vectorized operation
    mu = x
end subroutine identitylink

subroutine identitylinkinv(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: eta(n)
    
    ! Vectorized operation
    eta = x
end subroutine identitylinkinv

subroutine identitymueta(x, n, eta)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(out) :: eta(n)
    
    ! Removed redundant operation x = x
    ! Vectorized operation
    eta = 1.d0
end subroutine identitymueta