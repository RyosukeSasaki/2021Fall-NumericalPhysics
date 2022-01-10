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
  !!
  !write(0,*) 'linear size?'
  !read (*,*,end=9000) lsize
  !!
  !write(0,*) 'temperature(T/J)?'
  !read (*,*,end=9000) temp
  !!
  lsize = 8
  temp = 3.0
  
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
  !goto 1000
  !
9000 continue
  !
  stop
end program exact
