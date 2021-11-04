program main
    use differential
    implicit none
    DOUBLE PRECISION :: t = 0d0, x(2), tau = 0.01d0, t_end = 10d0
    INTEGER :: n, i
    x(1) = 0d0; x(2) = 1d0

    !do n = 2, 12
    !    t = 0d0
    !    tau = 1d0 / 2d0**dble(n)
    !    x(1) = 0d0; x(2) = 1d0
    !    do while (t < t_end)
    !        x = runge_kutta(f, 2, x, t, tau)
    !    end do
    !    WRITE(*, *) tau, t, abs(x(1) - x_true(t)), abs(x(2) - v_true(t))
    !end do
    do while (t <= t_end)
        x = runge_kutta(f, 2, x, t, tau)
        i = i + 1
        if(mod(i, 10) == 0) then
            WRITE(*, *) t, x(1), x_true(t), x(2), v_true(t)
        end if
    end do
    

    contains
    function f(tp, xp, n, const)
        implicit none
        INTEGER, INTENT(IN) :: n
        DOUBLE PRECISION :: f(n)
        DOUBLE PRECISION, INTENT(IN) :: xp(:), tp
        DOUBLE PRECISION, OPTIONAL :: const(:)
        f(1) = xp(2)
        f(2) = -xp(1)
    end function f
    function x_true(tp)
        implicit none
        DOUBLE PRECISION :: x_true
        DOUBLE PRECISION, INTENT(IN) :: tp
        x_true = sin(tp)
    end function x_true
    function v_true(tp)
        implicit none
        DOUBLE PRECISION :: v_true
        DOUBLE PRECISION, INTENT(IN) :: tp
        v_true = cos(tp)
    end function v_true
end program main