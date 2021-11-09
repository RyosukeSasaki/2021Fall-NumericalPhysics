program main
    use differential
    implicit none
    DOUBLE PRECISION, PARAMETER :: t_begin = 0.0, t_end = 10.0, x0(1) = 1.0
    DOUBLE PRECISION, PARAMETER :: tau = 0.001, interval = 0.1
    DOUBLE PRECISION :: t = 0.0, x(1)
    x(:) = x0(:)
    
    do while (t < t_end)
        x = runge_kutta(f, 1, x, t, tau)
        t = t + interval
        WRITE(*, *) t, x(1), x_true(t)
    end do

    contains
    function f(tp, xp, n, const)
        implicit none
        INTEGER, INTENT(IN) :: n
        DOUBLE PRECISION :: f(n)
        DOUBLE PRECISION, INTENT(in) :: tp, xp(:)
        DOUBLE PRECISION, OPTIONAL :: const(:)

        f(1) = xp(1)*cos(tp)
    end function f
    function x_true(tp)
        implicit none
        DOUBLE PRECISION :: x_true
        DOUBLE PRECISION, INTENT(in) :: tp
        
        x_true = exp(sin(tp))
    end function x_true
end program main