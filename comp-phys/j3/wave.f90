!======================================================================-
!  Schroedinger equation
!    A free particle in one dimension
!    <x> and <(x-<x>)**2> for a Gaussian wave packet
!    originally written by H. Takano
!    modified for the f90 form by J.Yamauchi 2011-10-06
!======================================================================-
program wave
  !======================================================================-
  !  sengen
  !======================================================================-
  implicit none
  integer,parameter:: idim=4097
  real(kind=8)::  a(4), b(4), r(idim), x(idim)
  complex(kind=8)::  cp(idim), cs(idim), cps(idim), cdp(idim)
  !
  integer:: iok,jmax,nmax,n1max,j,n,n1,m
  real(kind=8):: rl,sig,rk,x0,dt,dx,w,xj,sum,rnf,xav,xs
  complex(kind=8):: ci
  !======================================================================-
  !  nyuuryoku
  !======================================================================-
  !iok=0
  !do while (iok.eq.0)
  !  write(0,*) ' input Jmax,' &
  !       ,' where Jmax must be even and 2<= Jmax <=',idim-1,'.'
  !  read(*,*) jmax
  !  if (jmax.le.0) then
  !    write(0,*)'error: Jmax must be positive.'
  !  else if (jmax.ge.idim) then
  !    write(0,*)'error: Jmax must be less than',idim,'.'
  !  else if (mod(jmax,2).ne.0) then
  !    write(0,*) 'error: Jmax must be even.'
  !  else
  !    iok=1
  !  endif
  !enddo
  !!
  !iok=0
  !do while (iok.eq.0)
  !  write(0,*) ' input Nmax(> 0).'
  !  read(*,*) nmax
  !  if(nmax.le.0) then
  !    write(0,*)'error: Nmax must be positive.'
  !  else
  !    iok=1
  !  endif
  !enddo
  !!
  !iok=0
  !do while (iok.eq.0)
  !  write(0,*) ' input N1max (> 0).'
  !  read(*,*) n1max
  !  if(n1max.le.0) then
  !    write(0,*)'error: N1max must be positive.'
  !  else
  !    iok=1
  !  endif
  !enddo
  !!
  !!
  !iok=0
  !do while (iok.eq.0)
  !  write(0,*)'  input L(>0).'
  !  read(*,*) rl
  !  if(rl.le.0.0d0) then
  !    write(0,*)'error: L must be positive.'
  !  else
  !    iok=1
  !  endif
  !enddo
  !!
  !iok=0
  !do while (iok.eq.0)
  !  write(0,*)'  input sig(>0).'
  !  read(*,*) sig
  !  if(sig.le.0.0d0) then
  !    write(0,*)'error: sig must be positive.'
  !  else
  !    iok=1
  !  endif
  !enddo
  !!
  !write(0,*)'  input k0.'
  !read(*,*) rk
  !write(0,*)'  input x0.'
  !read(*,*) x0
  !!
  !iok=0
  !do while (iok.eq.0)
  !  write(0,*)'  input dt.'
  !  read(*,*) dt
  !  if(dt.le.0.0d0) then
  !    write(0,*)'error: dt must be positive.'
  !  else
  !    iok=1
  !  endif
  !enddo
  !
  Nmax = 100
  rl = 4d0
  sig = 0.1d0
  rk = 20d0
  x0 = -0.5d0
  jmax = 2*int(rl)*64
  dt = 1d0*10d0**(-5d0)
  n1max = int(0.05d0/dt/dble(nmax))
  write(*,'(''# Jmax  = '',i12)') jmax
  write(*,'(''# Nmax  = '',i12)') nmax
  write(*,'(''# N1max = '',i12)') n1max
  write(*,'(''# L     = '',e18.8e3)') rl
  write(*,'(''# sig   = '',e18.8e3)') sig
  write(*,'(''# k0    = '',e18.8e3)') rk
  write(*,'(''# x0    = '',e18.8e3)') x0
  write(*,'(''# dt    = '',e18.8e3)') dt
  !======================================================================-
  ! junbi
  !======================================================================-
  a(1)=0.0d0
  a(2)=0.5d0
  a(3)=0.5d0
  a(4)=1.0d0
  b(1)=1.0d0/6.0d0
  b(2)=1.0d0/3.0d0
  b(3)=1.0d0/3.0d0
  b(4)=1.0d0/6.0d0
  dx=2.0d0*rl/dble(jmax)
  w =0.5d0/dx**2
  ci=(0.0d0,1.0d0)
  do j=1,jmax+1
    x(j)=dx*dble(j-1)-rl
  enddo
  !======================================================================-
  ! shoki-jooken
  !======================================================================-
  sum=0.0d0
  do j=2,jmax
    xj=x(j)
    cp(j)=cdexp(ci*rk*xj-((xj-x0)/(2.0d0*sig))**2)
    sum=sum+cp(j)*dconjg(cp(j))
  enddo
  cp(1)     =(0.d0,0.d0)
  cp(jmax+1)=(0.d0,0.d0)
  rnf=1.0d0/dsqrt(sum*dx)
  do j=2,jmax
    cp(j)=rnf*cp(j)
  enddo
  !
  n=0
  write(*,'( ''#''/ ''# time'',12x,''  <x> '',12x, ''  <(x-<x>)**2>'' )')
  !======================================================================-
  ! shuukei
  !======================================================================-
  do j=2,jmax
    r(j)=cp(j)*dconjg(cp(j))
  enddo
  xav=0.0d0
  do j=2,jmax
    xav=xav+x(j)*r(j)*dx
  enddo
  xs=0.0d0
  do j=2,jmax
    xs=xs+((x(j)-xav)**2)*r(j)*dx
  enddo
  write(*,'(3e18.8e3)') dt*n1max*n,xav,xs
  !
  cdp(:)=0.0d0
  !======================================================================-
  do n=1,Nmax
    !==========================================================-
    do n1=1,N1max
      !==============================================-
      ! 1 step sekibun
      !==============================================-
      do j=2,jmax
        cs(j)=(0.0d0,0.0d0)
      enddo
      do m=1,4
        do j=2,jmax
          cps(j)=cp(j)+a(m)*cdp(j)
        enddo
        do j=2,jmax
          cdp(j)=(cps(j+1)-2.0d0*cps(j)+cps(j-1))*w*ci*dt
        enddo
        do j=2,jmax
          cs(j)=cs(j)+b(m)*cdp(j)
        enddo
      enddo
      do j=2,jmax
        cp(j)=cp(j)+cs(j)
      enddo
    enddo
    !==========================================================-
    ! shuukei
    !==========================================================-
    do j=2,jmax
      r(j)=cp(j)*dconjg(cp(j))
    enddo
    xav=0.0d0
    do j=2,jmax
      xav=xav+x(j)*r(j)*dx
    enddo
    xs=0.0d0
    do j=2,jmax
      xs=xs+((x(j)-xav)**2)*r(j)*dx
    enddo
    write(*,'(3e18.8e3)') dt*n1max*n,xav,xs
  enddo
  !======================================================================-
end program wave
