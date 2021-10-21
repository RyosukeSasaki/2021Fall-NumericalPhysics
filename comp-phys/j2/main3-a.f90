program main
    use differential
    implicit none
    DOUBLE PRECISION, PARAMETER :: t_begin = 0.0, t_end = 8.0
    DOUBLE PRECISION, PARAMETER :: tau = 0.01, interval = 0.1
    DOUBLE PRECISION :: t = 0.0, x(2)
    x(1) = 0.0; x(2) = 1.0

    do while (t < t_end)
        x = runge_kutta(f, 2, x, t, t + interval, tau)
        t = t + interval
        WRITE(*, *) t, x(1), x_true(t)
    end do

    contains
    function f(tp, xp, n)
        implicit none
        INTEGER, INTENT(IN) :: n
        DOUBLE PRECISION :: f(n)
        DOUBLE PRECISION, INTENT(IN) :: xp(:), tp
        f(1) = xp(2)
        f(2) = -xp(1)
    end function f
    function x_true(tp)
        implicit none
        DOUBLE PRECISION :: x_true
        DOUBLE PRECISION, INTENT(IN) :: tp
        x_true = sin(tp)
    end function x_true
end program main