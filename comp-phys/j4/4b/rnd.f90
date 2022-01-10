!======================================================================-
!     random number generator
!       based on Tausworthe's algorithm
!       (maximum-length linearly recurring sequence)
!       initialization by Fushimi's algorithm
!       (M. Fushimi, ``Ransu'',(Tokyo Daigaku Shuppan, 1989, Tokyo))
!======================================================================-
!     call rndini(iseed)
!       initializes random number generator
!       iseed: input, seed for initialization, integer*4.
!     call rndu(rnd)
!       generates a random number uniformly distributed in [0,1).
!       rnd: output, random number, real*8.
!     note:  rndini must be called before rndu is called.
!    modified for the f90 form by J.Yamauchi 2011-10-20
!======================================================================-
!
subroutine rndini(iseed)
  !
  implicit none
  !
  integer,parameter:: ip=521, iq=32, nq=ip-iq, nbit=32, nbit1=31 &
       ,ia=69069
  !
  integer:: iseed,iptr
  integer:: ir(ip),iw(ip)
  !
  integer:: i,j,ii,ij,mj,ih
  !
  common /rndpar/ir,iptr
  save   /rndpar/
  !
  !      if(iseed.le.0) then
  !        write(*,*) 'rndini: iseed must be positive.'
  !        stop
  !      end if
  !
  if(mod(iseed,2).eq.0) then
    write(*,*) 'rndini: iseed must be odd.'
    stop
  end if
  !
  do i=1,ip
    iseed=iseed*ia
    iw(i)=isign(1,iseed)
  enddo
  !
  do j=1,ip
    ih=mod((j-1)*nbit,ip)
    mj=0
    do i=1,nbit1
      ii=mod(ih+i-1, ip)+1
      mj=2*mj+(iw(ii)-1)/(-2)
      ij=mod(ii+nq-1,ip)+1
      iw(ii)=iw(ii)*iw(ij)
    enddo
    ir(j)=mj
    ii=mod(ih+nbit1,ip)+1
    ij=mod(ii+nq-1, ip)+1
    iw(ii)=iw(ii)*iw(ij)
  enddo
  !
  iptr=0
  !
  return
end subroutine rndini
!
!======================================================================-
!
subroutine rndu(rnd)
  !
  implicit none
  !
  integer,parameter:: ip=521, iq=32
  real(kind=8),parameter:: ra=2.0d0**(-31), rb=2.0d0**(-32)
  !
  integer::  iptr, ir(ip)
  integer:: jptr
  real(kind=8):: rnd
  !
  common /rndpar/ir,iptr
  save   /rndpar/
  !
  iptr=iptr+1
  if(iptr.gt.ip) iptr=1
  jptr=iptr-iq
  if(jptr.le.0) jptr=jptr+ip
  !
  ir(iptr)=ieor(ir(iptr),ir(jptr))
  rnd=dble(ir(iptr))*ra+rb
  !
  return
end subroutine rndu
