subroutine self_ph(ePt,nLead,norbs,dimph,pheigv,coup,fock0,Lambda, &
                    miu,telec, cselfr, cselfl,greenl,cselfg)
!
!  first-order BA, the self-energy is constructed by bare electron GF
! 
!  output of greater self-energy is added
!
use parameters
use transmod,    only: dBlkSeq, bocc
implicit none

complex*16,   intent(in)    :: ePt
integer,      intent(in)    :: nLead,norbs,dimph
real*8,       intent(in)    :: pheigv(dimph)
complex*16,   intent(in)    :: coup(norbs,norbs,dimph)
real*8,       intent(in)    :: fock0(norbs,norbs)
type(dBlkSeq),intent(inout) :: Lambda(nLead)
real*8,       intent(in)    :: miu(nLead), telec
complex*16,   intent(out)   :: cselfr(norbs,norbs), cselfl(norbs,norbs)
complex*16,   intent(out)   :: greenl(norbs,norbs), cselfg(norbs,norbs)
!
integer,    allocatable :: ipiv(:)
real*8,     allocatable :: tlamda(:,:)
complex*16, allocatable :: hmat(:,:)
complex*16, allocatable :: ctmp1(:,:), ctmp2(:,:), ctmp3(:,:)
complex*16, allocatable :: greenr(:,:),greeng(:,:)
!
integer    :: i,ni, nj, nk, kd, istat,info
real*8     :: fermi(nLead),tmp1
complex*16 :: ztmp1,ePtq
!
real*8, external :: FermiDirac
real*8, external :: BoseEinstein
!
allocate(ipiv(norbs))
allocate(hmat(norbs,norbs),tlamda(norbs,norbs))
allocate(ctmp1(norbs,norbs),ctmp2(norbs,norbs),ctmp3(norbs,norbs), STAT=istat)
allocate(greenr(norbs,norbs),greeng(norbs,norbs)  )
!
cselfr = czero
cselfl = czero
cselfg = czero

!
do nk = 1, dimph
do kd = 1, 2
  if(kd==1) then
    ePtq = ePt + dcmplx(pheigv(nk),0.d0)
  else 
    ePtq = ePt - dcmplx(pheigv(nk),0.d0)
  endif

  do i=1,nLead
     fermi(i)=FermiDirac(telec,dble(ePtq),miu(i))
  enddo

  !------------------------
  !energy-dependent Gamma !
  !------------------------

  call self_lead(ePtq,nLead,norbs,miu,fock0,Lambda,tlamda)
  hmat=dcmplx(fock0,-tlamda)

  greenr=-hmat
  do ni = 1, norbs
     greenr(ni,ni) = ePtq + greenr(ni,ni)
  enddo
   
  call zgetrf(norbs,norbs, greenr, norbs, ipiv, info)
  call zgetri(norbs, greenr, norbs, ipiv, ctmp2,norbs*norbs, info)! G^r(epsilon+omega)
!
  ctmp2 = czero
  ctmp3 = czero
  do i=1,nLead
    ! sigma^< = i f(epsilon-miu) * Lamda  ! lesser self-energy
    ! greater self-energy
    ctmp2 = ctmp2 + 2.d0*eye*fermi(i)*Lambda(i)%ge
    ctmp3 = ctmp3 + 2.d0*eye*(fermi(i)-1.d0)*Lambda(i)%ge
  enddo
!
  call zgemm('n','n',norbs,norbs,norbs,cunity,greenr,norbs,ctmp2,norbs,czero,ctmp1,norbs)
  call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,greenr,norbs,czero,greenl,norbs)  ! lesser green's function
  
  call zgemm('n','n',norbs,norbs,norbs,cunity,greenr,norbs,ctmp3,norbs,czero,ctmp1,norbs)
  call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,greenr,norbs,czero,greeng,norbs)  ! greater green's function
!
  if(kd==1) then
    ctmp1 = bocc(nk)*greenr - 0.5d0*greenl
  else
    ctmp1 = (1.d0+bocc(nk))*greenr + 0.5d0*greenl
  endif
!
  !call dcmulti(norbs,norbs,coup(1:norbs,1:norbs,nk),ctmp1, cunity, ctmp2)
  !call cdmulti(norbs,norbs,ctmp2,coup(1:norbs,1:norbs,nk), cunity, ctmp3)
  call zgemm('n','n',norbs,norbs,norbs,cunity,coup(1:norbs,1:norbs,nk),norbs,ctmp1,norbs,czero,ctmp2, norbs)
  call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp2,norbs,coup(1:norbs,1:norbs,nk),norbs,czero,ctmp3, norbs)

!  retarded and lesser self-energy
  if(kd==1) then
    ztmp1 = dcmplx(1.d0 + bocc(nk), 0.d0)
  else 
    ztmp1 = dcmplx(bocc(nk), 0.d0)
  endif
  !call dcmulti(norbs,norbs,coup(1:norbs,1:norbs,nk),greenl,ztmp1,  ctmp2)
  !call cdmulti(norbs,norbs,ctmp2,coup(1:norbs,1:norbs,nk), cunity, ctmp1)

  call zgemm('n','n',norbs,norbs,norbs,ztmp1,coup(1:norbs,1:norbs,nk),norbs,greenl,norbs,czero,ctmp2, norbs)
  call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp2,norbs,coup(1:norbs,1:norbs,nk),norbs,czero,ctmp1, norbs)

  cselfr = cselfr + ctmp3
  cselfl = cselfl + ctmp1

  ! greater self-energy
  if(kd==1) then
    ztmp1 = dcmplx(bocc(nk), 0.d0)
  else 
    ztmp1 = dcmplx(bocc(nk)+1.d0, 0.d0)
  endif
  !call dcmulti(norbs,norbs,coup(1:norbs,1:norbs,nk),greeng,ztmp1,ctmp2)
  !call cdmulti(norbs,norbs,ctmp2,coup(1:norbs,1:norbs,nk),cunity,ctmp1)

  call zgemm('n','n',norbs,norbs,norbs,ztmp1,coup(1:norbs,1:norbs,nk),norbs,greeng,norbs,czero,ctmp2, norbs)
  call zgemm('n','n',norbs,norbs,norbs,cunity,ctmp2,norbs,coup(1:norbs,1:norbs,nk),norbs,czero,ctmp1, norbs)
  
  !cselfg(kd,kd) = cselfg(kd,kd) + ctmp1(kd,kd)  ! RWA

  cselfg = cselfg + ctmp1
  
enddo !kd
enddo !nk

cselfr = 0.5d0 * (cselfg-cselfl)    

deallocate(hmat,tlamda,ipiv)
deallocate(greenr, greeng) 
deallocate(ctmp1,ctmp2,ctmp3, STAT=istat)
end
