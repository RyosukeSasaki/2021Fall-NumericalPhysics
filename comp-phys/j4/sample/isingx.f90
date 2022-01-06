!======================================================================-
! Monte Carlo simulation of 2-D Ising model (no field)
!    originally written by Y. Saito,
!    modified for the f90 form by J.Yamauchi 2011-10-20
!======================================================================-
program ising
  !======================================================================-
  ! sengen
  !======================================================================-
  implicit none
  integer,parameter:: lmax=64, nnno=4
  integer:: lsize,inist,imcs,nmcs,nblock,iseed &
       ,ispin(0:lmax-1,0:lmax-1),ip(0:lmax-1),im(0:lmax-1)
  real(kind=8):: temp &
       ,trprob(-nnno:nnno)
  !
  integer:: no,lsize1,i,ix,iy,nm,n,nb,ie,iene,imag
  real(kind=8):: de,amag,amag2,aene,aene2,ave,avc,dene,avm,avx,dmag,rnd
  !============================================================
  !  for graphics
  !============================================================
  integer,parameter:: igxlen=640, igylen=480, igxmin=0, igymin=0 &
       ,igxmid=igxlen/2, igymid=igylen/2, igxmax=igxlen-1, igymax=igylen-1 &
       ,igcolp=4, igcolm=1
  integer:: igw,iglen,igx0,igy0 &
       ,igx(0:lmax-1), igy(0:lmax-1)
  character*10 wname
  integer::  igcol(-1:1)
  data igcol/igcolm,0,igcolp/
  !======================================================================-
  ! parameter nyuryoku
  !======================================================================-
  !   nyuryoku
  !==========================================================-
  write(0,*) 'linear size?'
  read (*,*) lsize
  !
  write(0,*) 'temperature?'
  read (*,*) temp
  !
  write(0,*) 'initial configuration?'
  write(0,*) ' -1: all down'
  write(0,*) '  0: random'
  write(0,*) '  1: all up'
  read (*,*) inist
  !
  write(0,*) 'number of initial Monte Carlo steps?'
  read (*,*) imcs
  write(0,*) 'number of Monte Carlo steps per block?'
  read (*,*) nmcs
  write(0,*) 'number of blocks?'
  read (*,*) nblock
  !
  write(0,*) 'seed of random numbers (odd positive integer)?'
  read (*,*) iseed
  !==========================================================-
  !   kensa
  !==========================================================-
  if((lsize.gt.lmax).or.(lsize.lt.1)) then
    write(*,*) 'error in lsize:',lsize
    stop
  endif
  if((inist.lt.-1).or.(inist.gt.1)) then
    write(0,*) 'error in inist:',inist
    stop
  endif
  if((iseed.le.0).or.(mod(iseed,2).eq.0)) then
    write(0,*) 'error in iseed:',iseed
    stop
  endif
  if(imcs.lt.0) then
    write(0,*) 'error in imcs:',imcs
    stop
  endif
  if(nmcs.lt.0) then
    write(0,*) 'error in nmcs:',nmcs
    stop
  endif
  if(nblock.lt.0) then
    write(0,*) 'error in nblock:',nblock
    stop
  endif
  !==========================================================-
  !   shutsuryoku
  !==========================================================-
  write(*,'(a10,i12)')    '   lsize= ',lsize
  write(*,'(a10,e24.16)') '    temp= ',temp
  write(*,'(a10,i12)')    '   inist= ',inist
  write(*,'(a10,i12)')    '    imcs= ',imcs
  write(*,'(a10,i12)')    '    nmcs= ',nmcs
  write(*,'(a10,i12)')    '  nblock= ',nblock
  write(*,'(a10,i12)')    '   iseed= ',iseed
  !======================================================================-
  ! shokisettei
  !======================================================================-
  !   ransu shokika
  !==========================================================-
  call rndini(iseed)
  !==========================================================-
  !   tonari no hyo shokika
  !==========================================================-
  no=lsize*lsize
  lsize1=lsize-1
  do i=0,lsize1
    ip(i)=i+1
    im(i)=i-1
  enddo
  ip(lsize1)=0
  im(0)=lsize1
  !==========================================================-
  !   spin haichi shokika
  !==========================================================-
  do ix=0,lsize1
    do iy=0,lsize1
      if(inist.ge.1) then
        ispin(ix,iy)=1
      else if(inist.le.-1) then
        ispin(ix,iy)=-1
      else
        call rndu(rnd)
        if(rnd.ge.0.5d0) then
          ispin(ix,iy)=1
        else
          ispin(ix,iy)=-1
        endif
      endif
    enddo
  enddo
  !==========================================================-
  !   sen'i-kakuritsu shokika
  !==========================================================-
  do n=-nnno,nnno
    de=2.0d0*dble(n)/temp
    if(de.gt.0.0d0) then
      trprob(n)=dexp(-de)
    else
      trprob(n)=1.0d0
    endif
  enddo
  !==========================================================-
  !  for graphics
  !==========================================================-
  igw=igylen/lsize
  if(igw.le.0) then
    write(*,*) 'error: lsize =',lsize,' > ',igylen
    stop
  endif
  iglen=igw*lsize
  igx0=igxmid-iglen/2
  igy0=igymid+iglen/2
  do ix=0,lsize1
    igx(ix)=igx0+igw*ix
  enddo
  do iy=0,lsize1
    igy(iy)=igy0-igw*(iy+1)
  enddo
  wname(1:5)='ising'
  call ginit(wname)
  call gcls
  do ix=0,lsize1
    do iy=0,lsize1
      call gbox(igx(ix),igy(iy),igw,igw,igcol(ispin(ix,iy)))
    enddo
  enddo
  call gdisp
  !======================================================================-
  ! shoki loop
  !======================================================================-
  do nm=1,imcs
    !==========================================================-
    !   1 Monte Carlo step
    !==========================================================-
    do n=1,no
      call rndu(rnd)
      ix=lsize*rnd
      call rndu(rnd)
      iy=lsize*rnd
      ie=ispin(ix,iy)*(ispin(ip(ix),   iy )+ispin(im(ix),   iy ) &
           +ispin(   ix ,ip(iy))+ispin(   ix ,im(iy)))
      call rndu(rnd)
      if(rnd.lt.trprob(ie)) then
        ispin(ix,iy)=-ispin(ix,iy)
        call gbox(igx(ix),igy(iy),igw,igw,igcol(ispin(ix,iy)))
        !            call gdisp
      endif
    enddo
    call gdisp
  enddo
  !======================================================================-
  ! block loop
  !======================================================================-
  do nb=1,nblock
    !==========================================================-
    !   block-wa shokika
    !==========================================================-
    amag =0.0d0
    amag2=0.0d0
    aene =0.0d0
    aene2=0.0d0
    !==========================================================-
    !   Monte Carlo step loop
    !==========================================================-
    do nm=1,nmcs
      !==============================================-
      !     1 Monte Carlo step
      !==============================================-
      do n=1,no
        call rndu(rnd)
        ix=lsize*rnd
        call rndu(rnd)
        iy=lsize*rnd
        ie=ispin(ix,iy)*(ispin(ip(ix),   iy )+ispin(im(ix),   iy ) &
             +ispin(   ix ,ip(iy))+ispin(   ix ,im(iy)))
        call rndu(rnd)
        if(rnd.lt.trprob(ie)) then
          ispin(ix,iy)=-ispin(ix,iy)
          call gbox(igx(ix),igy(iy),igw,igw,igcol(ispin(ix,iy)))
          !              call gdisp
        endif
      enddo
      call gdisp
      !==============================================-
      !     sokutei
      !==============================================-
      iene=0
      imag=0
      do ix=0,lsize1
        do iy=0,lsize1
          imag=imag+ispin(ix,iy)
          iene=iene+ispin(ix,iy)*(ispin(ip(ix),iy)+ispin(ix,ip(iy)))
        enddo
      enddo
      !==============================================-
      !     block-wa kasan
      !==============================================-
      aene =aene +dble(iene)
      aene2=aene2+dble(iene)*dble(iene)
      amag =amag +dble(imag)
      amag2=amag2+dble(imag)*dble(imag)
    enddo
    !==========================================================-
    !   block heikin
    !==========================================================-
    aene =aene /dble(nmcs)
    aene2=aene2/dble(nmcs)
    amag =amag /dble(nmcs)
    amag2=amag2/dble(nmcs)
    write(*,'(a10,i12)')    '   block= ',nb
    write(*,'(a10,e24.16)') '    aene= ',aene
    write(*,'(a10,e24.16)') '   aene2= ',aene2
    write(*,'(a10,e24.16)') '    amag= ',amag
    write(*,'(a10,e24.16)') '   amag2= ',amag2
    dene =aene2-aene*aene
    ave  =-aene/dble(no)
    avc  =(dene/dble(no))/(temp*temp)
    dmag =amag2-amag*amag
    avm  =amag/dble(no)
    avx  =dmag/dble(no)
    write(*,'(a10,e24.16)') '       e= ',ave
    write(*,'(a10,e24.16)') '       c= ',avc
    write(*,'(a10,e24.16)') '       m= ',avm
    write(*,'(a10,e24.16)') '       x= ',avx
  enddo
  !
  stop
end program ising
