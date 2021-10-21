program main
    use differential
    implicit none
    DOUBLE PRECISION, PARAMETER :: t_begin = 0.0, t_end = 10.0, x0(1) = 1.0
    DOUBLE PRECISION, PARAMETER :: tau = 0.001, interval = 0.1
    DOUBLE PRECISION :: t = 0, x(1)
    x(:) = x0(:)
    
    do while (t < t_end)
        x = runge_kutta(f, 1, x, t, t + interval, tau)
        t = t + interval
        write(*, *) t, x(1), x_true(t)
    end do

    contains
    function f(tp, xp, n)
        implicit none
        INTEGER, INTENT(IN) :: n
        DOUBLE PRECISION :: f(n)
        DOUBLE PRECISION, INTENT(in) :: tp, xp(:)

        f(1) = xp(1)*cos(tp)
    end function f
    function x_true(t)
        implicit none
        DOUBLE PRECISION :: x_true
        DOUBLE PRECISION, INTENT(in) :: t
        
        x_true = exp(sin(t))
    end function x_true
end program main