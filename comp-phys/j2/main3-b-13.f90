program main
    use differential
    implicit none
    DOUBLE PRECISION, PARAMETER :: pi = 4d0*atan(1d0)
    DOUBLE PRECISION :: t_begin = 0d0, t_end
    DOUBLE PRECISION :: tau, t = 0d0, x(3) = 0d0
    DOUBLE PRECISION :: const(3)
    INTEGER :: interval = 1000, i = 0
    !const(:) : (Q, G, Omega)
    !xp(:) : (omega, theta, phi)
    const(1) = 2d0; const(3) = 2d0/3d0;
    t_end = 1050d0*2d0*pi/const(3)
    tau = 2d0*pi/const(3)*10d0**(-3d0)

    const(2) = 1.5d0
    do while (t <= 50d0*2d0*pi/const(3))
        x = runge_kutta(f, 3, x, t, tau, const, boundary)
    end do
    do while (t <= t_end)
        x = runge_kutta(f, 3, x, t, tau, const, boundary)
        if (mod(i, interval) == 0) then
            WRITE(*, *) x(3), x(2), x(1)
        end if
        i = i + 1
    end do
    
    contains
    function f(tp, xp, n, const)
        implicit none
        INTEGER, INTENT(IN) :: n
        DOUBLE PRECISION :: f(n)
        DOUBLE PRECISION, INTENT(IN) :: xp(:), tp
        DOUBLE PRECISION, OPTIONAL :: const(:)
        !xp(:) : (omega, theta, phi)
        !const(:) : (Q, G, Omega, pi)
        f(1) = -1.0/const(1)*xp(1) - sin(xp(2)) + const(2)*cos(xp(3))
        f(2) = xp(1)
        f(3) = const(3)
    end function f
    function boundary(xp, n)
        implicit none
        INTEGER, INTENT(IN) :: n
        DOUBLE PRECISION, INTENT(IN) :: xp(:)
        DOUBLE PRECISION :: boundary(n)
        !xp(:) : (omega, theta, phi)
        boundary(1) = xp(1)

        if (xp(2) > pi) then
            boundary(2) = xp(2) - 2*pi
        elseif (xp(2) <= -pi) then
            boundary(2) = xp(2) + 2*pi
        else
            boundary(2) = xp(2)
        end if
        if (xp(3) >= 2*pi) then
            boundary(3) = xp(3) - 2*pi
        elseif (xp(3) < 0) then
            boundary(3) = xp(3) + 2*pi
        else
            boundary(3) = xp(3)
        end if
    end function boundary
end program main