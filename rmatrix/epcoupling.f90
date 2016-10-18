subroutine ep_am15(ns,nst,mode,beta,freq,epcoup,pheigv,bocc,coup,tdielectric)
use parameters
implicit none
integer,   intent(in)  :: ns,nst, mode
real*8,    intent(in)  :: beta, freq, epcoup
real*8,    intent(out) :: pheigv(mode),bocc(mode)
complex*16,intent(out) :: coup(ns+nst,ns+nst,mode)
logical,   intent(in)  :: tdielectric
!
integer :: i,j,k
real*8  :: dE,tmp,flux,factor,power
complex*16 :: epsilonw,ztmp

!
pheigv= 0.d0
bocc  = 0.d0
coup  = cmplx(0.d0,0.d0)

dE=freq/dble(mode)

if(mode==1) dE=1.0d0

factor=0.d0
do k=1,mode
  pheigv(k)= dble(k-0.5d0)*dE
  if(mode==1) pheigv(k)=freq
  flux=pheigv(k)*pheigv(k)/(dexp(pheigv(k)*beta)-1.d0)
  factor=factor+flux*pheigv(k)*dE
enddo


write(6,*)
if(tdielectric) then
  write(6,*) ' node   energy        occ           flux      coupling-factor , (eps-1/(eps+2)-1'
else
  write(6,*) ' node   energy        occ           flux      coupling-factor '
endif
power=0.d0
do k=1,mode
  flux=pheigv(k)*pheigv(k)/(dexp(pheigv(k)*beta)-1.d0)/factor
  power=power+flux*pheigv(k)*dE

  bocc(k) = 1.d0/(dexp( beta * pheigv(k)) - 1.d0)
  !bocc(k) = 1.d0
  !dE is the weight of summation over photon mode
  !tmp=dsqrt(flux/bocc(k)/pheigv(k)*dE)  
  tmp=dsqrt(flux/bocc(k)*pheigv(k)*dE)  


  do i=1,nst/2
    do j=nst/2+1,nst
      coup(ns+j,ns+i,k)=dcmplx(epcoup*tmp,0.d0)
      coup(ns+i,ns+j,k)=dcmplx(epcoup*tmp,0.d0)
    enddo
  enddo

  if(tdielectric) then
    call dielectric(3,'Au',pheigv(k),epsilonw)
    ztmp = (epsilonw-1.d0)/(epsilonw+2.d0) - 1.d0
    coup(ns+1:ns+nst,ns+1:ns+nst,k) = coup(ns+1:ns+nst,ns+1:ns+nst,k) * ztmp
  endif

  if(tdielectric) then
    write(6,'(I4,2X, f9.4, 5e15.5)') k, pheigv(k), bocc(k), flux, tmp, ztmp
  else
    write(6,'(I4,2X, f9.4, 3e15.5)') k, pheigv(k), bocc(k), flux, tmp
  endif
  !
enddo

write(6,'(A,f9.4,A)') "  power density is:", power,'kW/m^2/ev'

end subroutine ep_am15
