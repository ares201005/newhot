subroutine absorp_ph(ePt,nLead,norbs,dimph,pheigv,coup,fock0,Lambda,miu,telec,weight,absorb)
!
!  first-order BA, the self-energy is constructed by bare electron GF
! 
!  output of greater self-energy is added
!
use parameters
use transmod,    only: dBlkSeq,bocc
implicit none

complex*16,   intent(in)    :: ePt
integer,      intent(in)    :: nLead,norbs,dimph
real*8,       intent(in)    :: pheigv(dimph)
complex*16,   intent(in)    :: coup(norbs,norbs,dimph)
real*8,       intent(in)    :: fock0(norbs,norbs)
type(dBlkSeq),intent(inout) :: Lambda(nLead)
real*8,       intent(in)    :: miu(nLead), telec
real*8,       intent(in)    :: weight
real*8,     intent(inout)   :: absorb(dimph)
!
real*8,     allocatable :: tlamda(:,:)
complex*16, allocatable :: hmat(:,:)
complex*16, allocatable :: ctmp1(:,:), ctmp2(:,:), ctmp3(:,:)
complex*16, allocatable :: greenr(:,:),greeng(:,:),greenl(:,:)
complex*16, allocatable :: greeng0(:,:),greenl0(:,:)
!
integer    :: ipiv(norbs), info
integer    :: i,ni, nj, nk, kd, istat
real*8     :: fermi(nLead),tmp1
complex*16 :: ztmp1,ePtq
!
real*8, external :: FermiDirac
real*8, external :: BoseEinstein
!
allocate(hmat(norbs,norbs),tlamda(norbs,norbs))
allocate(ctmp1(norbs,norbs),ctmp2(norbs,norbs),ctmp3(norbs,norbs), STAT=istat)
allocate(greenr(norbs,norbs),greeng(norbs,norbs),greenl(norbs,norbs)  )
allocate(greeng0(norbs,norbs),greenl0(norbs,norbs)  )
!

!-------------------------
! get G^<(E) and G^>(E)  !
!-------------------------

call self_lead(ePt,nLead,norbs,miu,fock0,Lambda,tlamda)
hmat=dcmplx(fock0,-tlamda)

greenr=-hmat
do ni = 1, norbs
   greenr(ni,ni) = ePt + greenr(ni,ni)
enddo

call zgetrf(norbs,norbs, greenr, norbs, ipiv, info)
call zgetri(norbs, greenr, norbs, ipiv, ctmp2,norbs*norbs, info)! G^r(epsilon+omega)
!
ctmp2 = czero
ctmp3 = czero
do i=1,nLead
  fermi(i)=FermiDirac(telec,dble(ePt),miu(i))
  ctmp2 = ctmp2 + 2.d0*eye*fermi(i)*Lambda(i)%ge
  ctmp3 = ctmp3 + 2.d0*eye*(fermi(i)-1.d0)*Lambda(i)%ge
enddo

! G^<(E)
call zgemm('n','n',norbs,norbs,norbs,cunity,greenr,norbs,ctmp2,norbs,czero,ctmp1,norbs)
call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,greenr,norbs,czero,greenl0,norbs)

! G^>(E) 
call zgemm('n','n',norbs,norbs,norbs,cunity,greenr,norbs,ctmp3,norbs,czero,ctmp1,norbs)
call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,greenr,norbs,czero,greeng0,norbs)

!--------------------------------------
do nk = 1, dimph
  ePtq = ePt - dcmplx(pheigv(nk),0.d0)

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

  !---------------
  !  G^<(E-w)    !
  !---------------
  call zgemm('n','n',norbs,norbs,norbs,cunity,greenr,norbs,ctmp2,norbs,czero,ctmp1,norbs)
  call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,greenr,norbs,czero,greenl,norbs)
  
  !---------------
  !  G^>(E-w)    !
  !---------------
  call zgemm('n','n',norbs,norbs,norbs,cunity,greenr,norbs,ctmp3,norbs,czero,ctmp1,norbs)
  call zgemm('n','c',norbs,norbs,norbs,cunity,ctmp1,norbs,greenr,norbs,czero,greeng,norbs)

  !--------------------------------------------
  ! M[G^>(E) M G^<(E-w) - G^<(E) M G^>(E-w)]  !
  !--------------------------------------------
  ztmp1 = dcmplx(bocc(nk), 0.d0)
  !call dcmulti(norbs,norbs,coup(1:norbs,1:norbs,nk),greenl,cunity,ctmp1)
  !call dcmulti(norbs,norbs,coup(1:norbs,1:norbs,nk),greeng,cunity,ctmp2)
  call zgemm('n','n',norbs,norbs,norbs,cunity,coup(1:norbs,1:norbs,nk),norbs,greenl,norbs,czero,ctmp1,norbs)
  call zgemm('n','n',norbs,norbs,norbs,cunity,coup(1:norbs,1:norbs,nk),norbs,greeng,norbs,czero,ctmp2,norbs)

  call zgemm('n','n',norbs,norbs,norbs,cunity,greeng0,norbs,ctmp1,norbs,czero,ctmp3,norbs)
  call zgemm('n','n',norbs,norbs,norbs,-1.d0*cunity,greenl0,norbs,ctmp2,norbs,cunity,ctmp3,norbs)

  !call dcmulti(norbs,norbs,coup(1:norbs,1:norbs,nk),ctmp3,ztmp1,ctmp2)
  call zgemm('n','n',norbs,norbs,norbs,ztmp1,coup(1:norbs,1:norbs,nk),norbs,ctmp3,norbs,czero,ctmp2,norbs)

  do i=1,norbs
    absorb(nk)=absorb(nk) + dble(ctmp2(i,i))*weight
  enddo
  
enddo !nk


deallocate(hmat,tlamda)
deallocate(greenr, greeng,greenl,greenl0,greeng0) 
deallocate(ctmp1,ctmp2,ctmp3, STAT=istat)
end subroutine absorp_ph
