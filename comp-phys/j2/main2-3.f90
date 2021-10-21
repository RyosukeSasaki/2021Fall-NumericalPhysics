program main
    use differential
    implicit none
    double precision, parameter :: t_begin = 0.0, t_end = 10.0, x0 = 1.0
    double precision, parameter :: tau = 0.01, interval = 0.1
    double precision :: t = 0, x = x0
    
    do while (t < t_end)
        x = runge_kutta(f, x, t, t + interval, tau)
        t = t + interval
        write(*, *) t, x
    end do

    contains
    function f(t, x)
        implicit none
        double precision :: f
        double precision, intent(in) :: t, x

        f = x*cos(t)
    end function f
end program main