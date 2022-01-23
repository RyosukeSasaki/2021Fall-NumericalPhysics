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
  integer:: no,lsize1,i,ix,iy,nm,n,nb,ie,iene,imag,j
  real(kind=8):: de,amag,amag2,aene,aene2,ave,avc,dene,avm,avx,dmag,rnd
  real(kind=8):: sume,sumc,summ,sumx
  real(kind=8):: rk,f,c,e
  !======================================================================-
  ! parameter nyuryoku
  !======================================================================-
  !   nyuryoku
  !==========================================================-
  !write(0,*) 'linear size?'
  !read (*,*) lsize
  !!
  !write(0,*) 'temperature?'
  !read (*,*) temp
  !!
  !write(0,*) 'initial configuration?'
  !write(0,*) ' -1: all down'
  !write(0,*) '  0: random'
  !write(0,*) '  1: all up'
  !read (*,*) inist
  !!
  !write(0,*) 'number of initial Monte Carlo steps?'
  !read (*,*) imcs
  !write(0,*) 'number of Monte Carlo steps per block?'
  !read (*,*) nmcs
  !write(0,*) 'number of blocks?'
  !read (*,*) nblock
  !!
  !write(0,*) 'seed of random numbers (odd positive integer)?'
  !read (*,*) iseed
  lsize = 8
  temp = 2.3
  inist = 0
  imcs = 10000
  nmcs = 10000
  nblock = 10
  
  !do j=1, 300
  !iseed = 1+2*j
  iseed = 601
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
  !write(*,'(a10,i12)')    '   lsize= ',lsize
  !write(*,'(a10,e24.16)') '    temp= ',temp
  !write(*,'(a10,i12)')    '   inist= ',inist
  !write(*,'(a10,i12)')    '    imcs= ',imcs
  !write(*,'(a10,i12)')    '    nmcs= ',nmcs
  !write(*,'(a10,i12)')    '  nblock= ',nblock
  !write(*,'(a10,i12)')    '   iseed= ',iseed
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
      endif
    enddo
  enddo
  !======================================================================-
  ! block loop
  !======================================================================-
  sume=0;sumc=0;summ=0;sumx=0;
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
        endif
      enddo
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
    !write(*,'(a10,i12)')    '   block= ',nb
    !write(*,'(a10,e24.16)') '    aene= ',aene
    !write(*,'(a10,e24.16)') '   aene2= ',aene2
    !write(*,'(a10,e24.16)') '    amag= ',amag
    !write(*,'(a10,e24.16)') '   amag2= ',amag2
    dene =aene2-aene*aene
    ave  =-aene/dble(no)
    avc  =(dene/dble(no))/(temp*temp)
    dmag =amag2-amag*amag
    avm  =amag/dble(no)
    avx  =dmag/dble(no)
    !write(*,'(a10,e24.16)') '       e= ',ave
    write(*,'(a10,e24.16)') '       c= ',avc
    !write(*,'(a10,e24.16)') '       m= ',avm
    !write(*,'(a10,e24.16)') '       x= ',avx
    sume = sume + ave
    sumc = sumc + avc
    summ = summ + avm
    sumx = sumx + avx
  enddo
  sume = sume / nblock
  sumc = sumc / nblock
  summ = summ / nblock
  sumx = sumx / nblock
  !write(*,'(a10,e24.16)') '       e= ',sume
  !write(*,'(a10,e24.16)') '       c= ',sumc
  !write(*,'(a10,e24.16)') '       m= ',summ
  !write(*,'(a10,e24.16)') '       x= ',sumx

  !write(*,'(i0, e24.16)') 1+2*j, sume

  !rk=1.0d0/temp
  !call calc(lsize,rk,f,e,c)
  !write(*,'(e24.16,e24.16,e24.16,e24.16,e24.16)') temp, sume, sumc, e, c
  !end do
  !
  stop
end program ising
