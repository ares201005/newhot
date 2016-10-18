!----------------------
! self-energy of lead !
!----------------------
subroutine self_lead(ePt,nLead,norbs,miu,fock0,Lambda,tlamda)
use variables, only : ns
use transmod, only: dBlkSeq,mps
implicit none
complex*16,    intent(in)    :: ePt
integer,       intent(in)    :: nLead
integer,       intent(in)    :: norbs
real*8,        intent(in)    :: miu(nLead)
real*8,        intent(in)    :: fock0(norbs,norbs)
type(dBlkSeq), intent(inout) :: Lambda(nLead)
real*8,        intent(out)   :: tlamda(norbs,norbs)
!
integer :: i,j
real*8  :: tmp1,tmp2

tlamda=0.d0

do i=1,nLead
  tmp1=Lambda(i)%width
  tmp1=tmp1*tmp1/((dble(ePt)-miu(i))*(dble(ePt)-miu(i))+tmp1*tmp1)

  tmp1=0.d0
  tmp2=0.d0
  if(dble(ePt)>=fock0(ns/2+1,ns/2+1)) then
    tmp2= sqrt(dble(ePt)-fock0(ns/2+1,ns/2+1))
  else if(dble(ePt)<=fock0(ns/2,ns/2)) then
    tmp1= sqrt(-dble(ePt)+fock0(ns/2,ns/2))
  endif

  Lambda(i)%ge = 0.d0
  Lambda(i)%ge(1:ns/2,1:ns/2)=Lambda(i)%g0(1:ns/2,1:ns/2)*tmp1
  Lambda(i)%ge(ns/2+1:ns,ns/2+1:ns)=Lambda(i)%g0(ns/2+1:ns,ns/2+1:ns)*tmp2

  if(norbs>ns) Lambda(i)%ge(ns+1:norbs,ns+1:norbs)=Lambda(i)%g0(ns+1:norbs,ns+1:norbs)

  tlamda=tlamda+Lambda(i)%ge

enddo



end subroutine self_lead
