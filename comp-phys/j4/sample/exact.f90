!======================================================================-
!
!       The internal energy and specific heat of the Ising model
!   on the L x L squrare lattice with the periodic boundary condition
!               exact calculation
!    originally written by Y. Saito.
!                  1995.06.07.
!
!    modified for the f90 form by J.Yamauchi 2011-10-20
!======================================================================-
!
program exact
  implicit none
  integer,parameter:: lmax=50
  integer:: lsize
  real(kind=8):: temp,rk,f,e,c
  !
1000 continue
  !
  write(0,*) 'linear size?'
  read (*,*,end=9000) lsize
  !
  write(0,*) 'temperature(T/J)?'
  read (*,*,end=9000) temp
  !
  if(lsize.le.1)  stop
  if(lsize.gt.lmax) stop
  !
  rk=1.0d0/temp
  !
  call calc(lsize,rk,f,e,c)
  !
  write(*,'(a10,i12)')     '   lsize= ',lsize
  write(*,'(a10,e24.16)')  '    temp= ',temp
  write(*,'(a10,e24.16)')  '       e= ',e
  write(*,'(a10,e24.16/)') '       c= ',c
  !
  goto 1000
  !
9000 continue
  !
  stop
end program exact
!
!======================================================================-
subroutine calc (n,rk,fn,en,cn)
  !======================================================================-
  !
  implicit none
  !
  integer,parameter:: lmax=50, lmax2=lmax*2
  integer:: n, n2, nn, i, i1, k, i20, i21 
  real(kind=8):: fn, en, cn, rk, rk2,ch, sh, th, pai, ccs, pn, sc &
       ,c0, c1, c2, dnh &
       ,g00, g01, g10, g11, g20, g21, gsum0, gsum1 &
       ,dc0, dc1, de0, de1, ds0, ds1, dt0, dt1 &
       ,edg, zt0, zt1, zt2 & 
       ,c(lmax2), g(lmax2), g1(lmax2), g2(lmax2), z(4), z1(4), z2(4)
  !
  !======================================================================-
  ! parameter check
  !======================================================================-
  !
  fn=0d0
  en=0d0
  cn=0d0
  !
  if(n.le.1)  return
  if(n.gt.lmax) return
  !
  !======================================================================-
  ! set constants
  !======================================================================-
  !
  rk2=rk*2d0
  n2 =n*2
  nn =n*n
  ch =dcosh(rk2)
  sh =dsinh(rk2)
  th =dtanh(rk2)
  pai=dacos(0d0)*2d0
  !
  !======================================================================-
  ! calculate c's
  !======================================================================-
  !
  ccs=ch*ch/sh
  pn=pai/n
  !
  do i=1,n2
    i1=i-1
    c(i)=ccs-dcos(i1*pn)
  enddo
  !
  c1=2d0*ch*(1d0-1d0/(sh*sh))
  c2=8d0*ch*ch/sh**3+4d0*(sh-1d0/sh)
  !
  !======================================================================-
  ! calculate g's
  !======================================================================-
  !
  dnh=dfloat(n)/2d0
  !
  g (1)=(rk2+dlog(tanh(rk)))*dnh
  g1(1)=(2d0*(1d0+1d0/sh)  )*dnh
  g2(1)=(-4d0*ch/(sh*sh)   )*dnh
  !
  do i=2,n2
    c0=c(i)
    sc=sqrt(c0*c0-1d0)
    g (i)=(dlog(c0+sc)         )*dnh
    g1(i)=(c1/sc               )*dnh
    g2(i)=(c2/sc-c1*c1*c0/sc**3)*dnh
  enddo
  !
  !======================================================================-
  ! calculate Z's
  !======================================================================-
  !
  do k=1,4
    z (k)=1d0
    z1(k)=0d0
    z2(k)=0d0
  enddo
  !
  gsum1=0d0
  gsum0=0d0
  !
  do i=1,n
    !
    i21=2*(i-1)+1+1
    i20=2*(i-1)+1
    !
    g01=g (i21)
    g00=g (i20)
    g11=g1(i21)
    g10=g1(i20)
    g21=g2(i21)
    g20=g2(i20)
    !
    dc1=dcosh(g01)
    dc0=dcosh(g00)
    !
    ds1=dsinh(g01)
    ds0=dsinh(g00)
    !
    dt1=dtanh(g01)
    dt0=dtanh(g00)
    !
    de1=dexp(-2d0*g01)
    de0=dexp(-2d0*g00)
    !
    gsum1=gsum1+g01
    gsum0=gsum0+g00
    !
    z (1)=z (1)*(1d0+de1)
    z (2)=z (2)*(1d0-de1)
    z (3)=z (3)*(1d0+de0)
    z (4)=z (4)*(1d0-de0)
    !
    z1(1)=z1(1)+g11*dt1
    z1(2)=z1(2)+g11/dt1
    z1(3)=z1(3)+g10*dt0
    z1(4)=z1(4)+g10/dt0
    !
    z2(1)=z2(1)+g21*dt1+g11*g11/(dc1*dc1)
    z2(2)=z2(2)+g21/dt1-g11*g11/(ds1*ds1)
    z2(3)=z2(3)+g20*dt0+g10*g10/(dc0*dc0)
    z2(4)=z2(4)+g20/dt0-g10*g10/(ds0*ds0)
    !
  enddo
  !
  do k=1,4
    z2(k)=z2(k)+z1(k)**2
  enddo
  !
  edg=dexp(gsum0-gsum1)
  z(3) =z (3)*edg
  z(4) =z (4)*edg
  !
  do k=1,4
    z1(k)=z1(k)*z (k)
    z2(k)=z2(k)*z (k)
  enddo
  !
  !======================================================================-
  !                               calculate f,e and c
  !======================================================================-
  !
  zt0=0d0
  zt1=0d0
  zt2=0d0
  !
  do k=1,4
    zt0=zt0+z (k)
    zt1=zt1+z1(k)
    zt2=zt2+z2(k)
  enddo
  !
  fn=0.5d0*dlog(2d0*sh)+(-dlog(2d0)+dlog(zt0)+gsum1)/nn
  !     fn=-fn/rk
  !
  en=-ch/sh-zt1/(zt0*nn)
  !
  cn=-2d0/(sh*sh)+(zt2/zt0-(zt1/zt0)**2)/nn
  cn=cn*rk*rk
  !
  return
end subroutine calc
