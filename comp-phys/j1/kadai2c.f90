module integral
    implicit none

    contains
    function simpson(arg, a, b, n)
        implicit none
        interface
            function arg(x)
                double precision :: arg
                double precision, intent(in) :: x
            end function arg
        end interface 
        double precision :: simpson, h, sum, x
        double precision, intent(in) :: a, b
        integer :: i
        integer, intent(in) :: n
        h = (b - a) / n

        sum = arg(a) + arg(b)
        do i = 1, n - 1, 2
            x = a + i * h
            sum = sum + 4d0*arg(x)
        end do
        do i = 2, n - 1, 2
            x = a + i*h
            sum = sum + 2d0*arg(x)
        end do
        simpson = sum*h / 3d0
    end function simpson
end module integral

module funcs
    use integral
    implicit none
    double precision, parameter :: pi = 3.14159265358979323846264

    contains
    function debye(x, n)
        implicit none
        double precision :: debye
        double precision, intent(in) :: x
        integer, intent(in) :: n
        
        debye = 3d0*x**3d0*simpson(integrand, 0d0, 1d0 / x, n)
        
        contains
        function integrand(x)
            implicit none
            double precision :: integrand, f
            double precision, intent(in) :: x
            f = x**4d0*exp(x) / (exp(x) - 1)**2d0
            if (isnan(f)) then
                integrand = 0
            else
                integrand = f
            end if
        end function integrand
    end function debye

    function dulong_petit()
        implicit none
        double precision :: dulong_petit
        dulong_petit = 1d0
    end function dulong_petit

    function low_temperature_behavior(x)
        implicit none
        double precision :: low_temperature_behavior
        double precision, intent(in) :: x
        low_temperature_behavior = 4d0*pi**4d0*x**3d0 / 5d0
    end function low_temperature_behavior
end module funcs

program  main
    use funcs
    implicit none
    integer :: n, i
    integer, parameter :: num_samples = 1000
    double precision :: x, f_exact, f
    double precision, parameter :: x_start = 0d0, x_end = 2d0
    double precision, parameter :: interval = (x_end - x_start) / num_samples
    open(18, file='output2c-1', status='replace')
    do n = 1, num_samples
        x = x_start + interval * n
        write(18, *) x, debye(x, 2**10), dulong_petit(), low_temperature_behavior(x)
    end do
    close(18)
    open(18, file='output2c-2', status='replace')
    f_exact = debye(x_end, 2**10)
    do i = 3, 9
        n = 2**i
        write(18, *) (1d0 / x_end - 0d0)/n, (f_exact - debye(x_end, n))/f_exact
    end do 
    close(18)

end program  main