program main
  implicit none
  INTEGER :: ndim,n,i
  DOUBLE PRECISION :: eigval
  DOUBLE PRECISION, PARAMETER :: pi = 4d0*atan(1d0)
  write(0,*) 'number of lattice points ?'
  read (*,*)  ndim
  
  n=int(dble(ndim)/2d0)
  do i=-n+1, n
    eigval = 2*abs(sin(pi*dble(i)/(2*dble(n))))
    write(*,*) i, eigval
  end do
end program main