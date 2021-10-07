program main
      implicit none
      integer, parameter :: nmax = 10
      integer :: n
      double precision :: a, b, c, d

      do n = 1, nmax
            a = dble(n)
            b = a**(1.0d0/3.0d0)
            c = a**(1.0d0/4.0d0)
            d = a**(1.0d0/5.0d0)
            call sleep(1)
      enddo

end
