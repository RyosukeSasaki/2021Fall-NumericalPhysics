program main
    implicit none
    integer :: n
    integer, parameter :: num_samples = 10000
    double precision :: x
    double precision, parameter :: x_start = 0, x_end = 10
    double precision, parameter :: interval = x_end / num_samples

    do n = 1, num_samples
        x = x_start + interval * n
        write(*, *) x, planck(x), rayleigh_jeans(x), wien(x)
    end do
    
    contains
    function planck(x)
        implicit none
        double precision :: planck, x
        planck = x**3d0 / (exp(x) - 1)
    end
    
    function rayleigh_jeans(x)
        implicit none
        double precision :: rayleigh_jeans, x
        rayleigh_jeans = x**2d0
    end
    
    function wien(x)
        implicit none
        double precision :: wien, x
        wien = x**3d0 * exp(-x)
    end
end