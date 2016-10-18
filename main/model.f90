program model
use transmod
use variables
use parameters
implicit none

character*10 :: date,time

real*8, allocatable :: fock0(:,:)
real*8, allocatable :: miu(:)

integer :: ierror(10)=0
integer :: i, j, k, l
!
type(dBlkSeq), allocatable :: Lambda(:)

call date_and_time(date,time)
write(6,*)
write(6,1000) date(7:8),date(5:6),date(1:4),time(1:2),time(3:4),time(5:6)


call io

allocate( miu(nLead), Lambda(nLead) )

miu(1)=miu0+ampl
miu(2)=miu0+ampr
if(nLead>2) miu(3)=miu0+ampG

norbs=ns
if(tPlasmon) then
  norbs = ns + metal_ps%nst
endif

do i=1,nLead
  allocate(Lambda(i)%g0(norbs,norbs))
  allocate(Lambda(i)%ge(norbs,norbs))
  Lambda(i)%g0=0.d0
  Lambda(i)%ge=0.d0
enddo
!
allocate( fock0(norbs,norbs))
allocate(lamdaL(norbs,norbs), lamdaR(norbs,norbs) )
allocate(deltapt(norbs,norbs), deltah(norbs,norbs))
allocate(newhmat(norbs,norbs) )

fock0=0.d0
fock0(1,1) = -onsite
fock0(2,2) =  onsite

if(tPlasmon) call calcFock0(ns,norbs,metal_ps,Fock0)

!-------------------------
!  line-width function   !
!-------------------------

lamdaL = 0.d0
lamdaR = 0.d0
do i=1,nst
 lamdaL(i,i)=eImag
 lamdaR(i,i)=eImag
enddo

lamdaL(1,1) = elcoup;
lamdaR(2,2) = elcoup
if(ampr<onsite) then
  dtmp=exp(-(onsite-ampr)/telec)
else
  dtmp=1.d0
endif
lamdaL(2,2) = elcoup*dtmp
lamdaR(1,1) = elcoup*dtmp 

Lambda(1)%g0 = LamdaL
Lambda(2)%g0 = LamdaR

if(nLead>2) then
  Lambda(3)%g0(1,1)=elcoup
  Lambda(3)%g0(2,2)=elcoup
endif

do i=1,nLead
  Lambda(i)%width=half_w(i)
enddo


write(6,*)
write(6,'(A26,2e15.7)') "  fock0(1,1), fock0(2,2) :", fock0(1,1),fock0(2,2)
write(6,'(A26,2e15.7)') "  lamdaL(1,1),lamdaL(2,2):", lamdaL(1,1),lamdaL(2,2)
write(6,'(A26,2e15.7)') "  lamdaR(1,1),lamdaR(2,2):", lamdaR(1,1),lamdaR(2,2)

call landauer(nLead,norbs,miu,beta,fock0,Lambda)

if(lphoton.or.tplasmon) then
  allocate( coup(norbs,norbs,mode))
  allocate( pheigv(mode),bocc(mode) )
endif

if(lphoton) then
  call ep_am15(0,norbs,mode,pt_beta,phfreq,epcoup,pheigv,bocc,coup,.false.)
endif
!
if(tplasmon) then
  if(t_read_file) then
    !------------------------------
    !read dipole matrix from file !
    !------------------------------
   call read_coup(ns,nst,mode,pt_beta,phfreq,pheigv,bocc,coup,'dipole-x.dat',.true.)
   !call fermi_golden(nst,phfreq)
  else
    call ep_am15(ns,metal_ps%nst,mode,pt_beta,phfreq,eps_coup,pheigv,bocc,coup,.true.)
  endif
endif

!------------------------------
! get steady state current    !
!------------------------------

if(lphoton.or.tplasmon) then
  call landauerphlop(nLead,norbs,mode,fock0,Lambda,miu,telec,pheigv,coup)
endif

call date_and_time(date,time)
write(6,1014) date(7:8),date(5:6),date(1:4),time(1:2),time(3:4),time(5:6)

1000 format(/,'  Calculation started on ',/,&
            ' ',A2,"/",A2,"/",A4,"  ",A2,":",A2,":",A2,/)
1001 format(A,e15.7)
1002 format(A,f15.7)
1003 format(A,I5)

1014 format(/,'  Calculation Finished on ',/,&
            ' ',A2,"/",A2,"/",A4,"  ",A2,":",A2,":",A2)

!---------------
! deallocation !
!---------------
deallocate(deltapt,deltah,newhmat,lamdaL,lamdaR)
deallocate(miu)
do i=1,nLead
  deallocate(Lambda(i)%g0,Lambda(i)%ge)
enddo
deallocate(Lambda)

if(lphoton) then
  deallocate( pheigv, coup,bocc )
endif
!
end program model
