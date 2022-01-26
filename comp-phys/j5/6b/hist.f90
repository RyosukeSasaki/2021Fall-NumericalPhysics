!======================================================================-
! histogram
!    originally written by M. Eto,
!    modified for the f90 form by J.Yamauchi 2011-10-20
!======================================================================-
program hist
  !======================================================================-
  ! sengen
  !======================================================================-
  implicit none
  !
  integer,parameter:: nmax  =200
  !
  integer:: ndiv,n,no,icount(0:nmax)
  !
  real(kind=8):: xmin,xmax,xdiv,x
  !
  character*60 ifname
  !
  !======================================================================-
  ! parameter nyuryoku
  !======================================================================-
  !==========================================================-
  !   nyuryoku
  !==========================================================-
  write(0,*) 'saisho-chi?'
  read (*,*)  xmin
  write(0,*) 'saidai-chi ?'
  read (*,*)  xmax
  write(0,*) 'bunkatsu-su ?'
  read (*,*)  ndiv
  write(0,*) 'input-file-name ?'
  read (*,*)  ifname
  !==========================================================-
  !   kensa
  !==========================================================-
  if((ndiv.gt.nmax).or.(ndiv.lt.1)) then
    write(0,*) 'error in ndiv:',ndiv
    stop
  endif
  !======================================================================-
  ! junbi
  !======================================================================-
  xdiv=(xmax-xmin)/dble(ndiv)
  do n=0,ndiv
    icount(n)=0
  enddo
  no=0
  !======================================================================-
  ! shuukei
  !======================================================================-
  open(10,file=ifname)
  do while (.TRUE.)
    read(10,*,end=9000) x
    n=int((x-xmin)/xdiv)
    if((n.ge.0).and.(n.le.ndiv)) then
      icount(n)=icount(n)+1
    endif
    no=no+1
  enddo
9000 continue
  close(10)
  !======================================================================-
  ! shutsuryoku
  !======================================================================-
  do n=0,ndiv
    write(*,*) xmin+(dble(n)+0.5d0)*xdiv &
         ,dble(icount(n))/(dble(no)*xdiv)
  enddo
  !
  stop
end program hist
