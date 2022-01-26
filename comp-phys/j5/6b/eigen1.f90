!======================================================================-
! 1-jigen koshi-shindo (doshu-ryushi)
!     originally written by M. Eto,
!     modified for the f90 form by J.Yamauchi 2011-10-20
!======================================================================-
program phonon
  !======================================================================-
  ! sengen
  !======================================================================-
  implicit none
  !
  integer,parameter:: nmax  =400, mmax  =nmax*(nmax+1)/2
  real(kind=8),parameter:: eps   =1.0d-8
  !
  integer:: ndim,ip,im,i,j,k,icon
  real(kind=8):: rmatrx(mmax),eigval(nmax),eigvec(nmax,nmax),gamma,bodd,beven,bi,bj
  ! For LAPACK
  real(kind=8):: work(nmax*3)
  character*1,parameter:: JOBZ='V',UPLO='U'
  !======================================================================-
  ! parameter nyuryoku
  !======================================================================-
  !==========================================================-
  !   nyuryoku
  !==========================================================-
  write(0,*) 'number of lattice points ?'
  read (*,*)  ndim
  write(0,*) 'gamma ?'
  read (*,*)  gamma
  !==========================================================-
  !   kensa
  !==========================================================-
  if((ndim.gt.nmax).or.(ndim.lt.1)) then
    write(*,*) 'error in ndim:',ndim
    stop
  endif
  if(gamma.eq.0d0) then
    write(*,*) 'error in gamma:',gamma
    stop
  endif
  bodd = 1d0
  beven = 1d0/sqrt(gamma)
  !======================================================================-
  ! gyoretsu-yoso no keisan
  !======================================================================-
  ! Make MATRIX A(j,i)=A(i*(i-1)/2+j)
  do i=1,ndim
    if(mod(i,2).eq.1) then
        bi = bodd
    else
        bi = beven
    endif

    if(i.eq.ndim) then
      ip=1
    else
      ip=i+1
    endif
    if(i.eq.1) then
      im=ndim
    else
      im=i-1
    endif
    do j=1,i
      if(mod(j,2).eq.1) then
        bj = bodd
      else
        bj = beven
      endif

      k=(i*(i-1)/2)+j
      if (j.eq.i) then
        rmatrx(k)= 2.0d0
      else if ((j.eq.ip).or.(j.eq.im)) then
        rmatrx(k)=-1.0d0
      else
        rmatrx(k)= 0.0d0
      end if
      rmatrx(k) = rmatrx(k)*bi*bj
    end do
  end do
  !======================================================================-
  ! gyoretsu no taikaku-ka (LAPACK)
  !======================================================================-
  !     call DSPEV(JOBZ,UPLO,N,AP,W,Z,LDZ,WORK,INFO)
  !               INFO = 0  : Normal end
  !                    > 0  : Not convergent
  !                    < 0  : (=-i) The value of the i-th argument is wrong
  call DSPEV(JOBZ,UPLO,ndim,rmatrx,eigval,eigvec,nmax,work,icon)
  !==========================================================-
  !   error check
  !==========================================================-
  if(icon.ne.0) then
    write(*,*) 'lapack(dspev) error: icon=', icon
  endif
  !======================================================================-
  ! shutsuryoku
  !======================================================================-
  do j=1,ndim
    if((eigval(j).lt.0.0d0).and.(eigval(j).gt.-eps)) then
      eigval(j)=0.0d0
    endif
    eigval(j)=dsqrt(eigval(j))
    write(*,*) eigval(j)
  enddo
  !
  stop
end program phonon
