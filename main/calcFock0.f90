subroutine calcFock0(ns,norbs,metal_ps,fock0)
use transmod, only: mps
use variables, only : telec
implicit none
!
integer,  intent(in)    :: ns
integer,  intent(in)    :: norbs
type(mps),intent(in)    :: metal_ps
real*8,   intent(inout) :: fock0(norbs,norbs)
!
integer :: i,j,k
real*8  :: dtmp

do i=1,metal_ps%nst
 j = i+ns
 fock0(j,j) = metal_ps%eig(i)
enddo

do i=1,ns
  do j=1,metal_ps%nst
    k = j+ns

    if(metal_ps%eig(j)<=0.d0.and.i>(ns/2)) cycle
    if(metal_ps%eig(j)>0.d0.and.i<(ns/2)) cycle

    dtmp = 0.d0
    if(metal_ps%eig(j)>=metal_ps%barrier_e) then
      dtmp=metal_ps%mscoup
    else if(metal_ps%eig(j)>=0.d0.and.metal_ps%eig(j)<metal_ps%barrier_e) then
      dtmp=metal_ps%mscoup*dexp((metal_ps%eig(j)-metal_ps%barrier_e)/telec)
    else if(metal_ps%eig(j)<0.d0.and.metal_ps%eig(j)>(-metal_ps%barrier_h)) then
      dtmp=metal_ps%mscoup*dexp((-metal_ps%eig(j)-metal_ps%barrier_h)/telec)
    else
      dtmp=metal_ps%mscoup
    endif

    fock0(k,i) = dtmp
    fock0(i,k) = dtmp
  enddo
enddo

end subroutine calcFock0
