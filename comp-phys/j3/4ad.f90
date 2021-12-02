program main
    implicit none
    DOUBLE PRECISION, PARAMETER :: dt = 1d0*10d0**(-5d0)
    DOUBLE PRECISION, PARAMETER :: x0 = -0.5d0, k0 = 20d0, sigma=0.1d0
    DOUBLE PRECISION :: xav, xs, t
    INTEGER, PARAMETER :: nmax = 100
    INTEGER :: n1max = int(0.05d0/dt/dble(nmax)), n
    write(*,'( ''#''/ ''# time'',12x,''  <x> '',12x, ''  <(x-<x>)**2>'' )')
    do n=0,nmax
        t = dt*n1max*n
        xav = -0.5d0 + k0 * t
        xs = sigma**2d0 +t**2d0 / 4d0 / sigma**2d0
        write(*,'(3e18.8e3)') t,xav,xs
    end do
end program main