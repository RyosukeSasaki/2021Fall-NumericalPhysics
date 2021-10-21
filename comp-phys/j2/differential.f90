module differential
    implicit none
     
    contains
    function runge_kutta(arg, n, init, t_begin, t_end, tau)
        implicit none
        interface
            function arg(t, x, n)
                INTEGER, INTENT(IN) :: n
                DOUBLE PRECISION :: arg(n)
                DOUBLE PRECISION, INTENT(in) :: x(:), t
            end function arg
        end interface
        INTEGER :: i
        INTEGER, PARAMETER :: order = 4
        INTEGER, INTENT(IN) :: n
        DOUBLE PRECISION, INTENT(in) :: init(:), t_begin, t_end, tau
        DOUBLE PRECISION :: runge_kutta(n), x(n), t, s(n), delta(n)
        DOUBLE PRECISION :: a(order), b(order)
        a(1)=0.0; a(2)=0.5; a(3)=0.5; a(4)=1.0
        b(1)=1.0/6.0; b(2)=1.0/3.0; b(3)=1.0/3.0; b(4)=1.0/6.0
        x(:) = init(:)
        t = t_begin

        do while (t < t_end)
            s = 0; delta = 0
            do i = 1, order
                delta = arg(t+a(i)*tau, x(:) + a(i) * delta, n)*tau
                s(:) = s(:) + b(i) * delta
            end do
            x(:) = x(:) + s(:); t = t + tau
        end do
        runge_kutta(:) = x(:)
    end function runge_kutta
end module differential