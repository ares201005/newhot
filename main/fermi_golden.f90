subroutine fermi_golden(nst,freq)
implicit none
integer, intent(in) :: nst
real*8,  intent(in) :: freq

real*8,  allocatable :: eig(:), dipole(:,:)
!
character*10 :: chartmp
real*8  :: eshift, dtmp(2), delta,dE
integer :: i,j

allocate(eig(nst),dipole(nst,nst))
open(333,file="eig.dat",status='old')
read(333,*) eshift
do i=1,nst
  read(333,*) j, eig(i)
  eig(i) = eig(i) - eshift
  write(6,*) j, eig(i)
enddo

open(444,file="dipole-x.dat",status='old')
read(444,*) chartmp
read(444,*) chartmp
read(444,*) chartmp
read(444,*) chartmp
do i=1,nst
  read(444,*) dipole(1:nst,i)
enddo
close(444)

delta = 0.1
dtmp=0.d0
do i=1,nst
  do j=i,nst
   if(eig(j)>0.d0.and.eig(i)<0.d0) then
    dE = eig(j) - eig(i) - freq
    dtmp(2) = delta / (dE*dE + delta*delta)
    dtmp(1) = dtmp(1) + dtmp(2) * dipole(j,i)*dipole(j,i)
   endif
  enddo
enddo


write(6,'(A,e15.6)') "fermi-golden abs = ", dtmp(1)

end subroutine fermi_golden

