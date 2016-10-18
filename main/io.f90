subroutine io
use parameters
use transmod
use variables
implicit none

namelist /general/ tplasmon,lphoton,lscba
namelist /system/ nLead, ns, onsite, hop, elcoup,miu0,half_w
namelist /temperature/ telec
namelist /voltage/ ampl, ampr,ampG
namelist /photon/ mode, epcoup, phfreq, pt_telec
namelist /plasmon/ barrier_e,barrier_h,eps_coup,relax,nst,eig,mscoup,t_read_file
!
integer :: i

eImag = 0.d0

!---------------------------
lscba=.false.
lphoton = .false.
tplasmon= .false.
rewind(5)
read(5,nml=general,end=101)
101 continue
!
nLead = 2
ns    = 2
onsite = 0.d0
hop = 1.d0
elcoup = 1.d0
miu0 = 0.d0
half_w =0.1d0
rewind(5)
read(5,nml=system,end=102)
102 continue

write(6,1003) ' ns   ', ns
write(6,1003) ' nLead   ', nLead
write(6,1001) ' elcoup: ', elcoup
write(6,1001) ' hop     ', hop
write(6,1001) ' miu0    ', miu0
write(6,1001) ' onsite  ', onsite
write(6,'(A20,5f12.4)') ' half_w  ', half_w(1:5)
write(6,*)
!
telec = 300  
rewind(5)
read(5,nml=temperature,end=103)
103 continue

write(6,1001) ' temperature: ', telec

beta = ev2j / (boltz * telec)         ! 1/kT in 1/ev  
telec = (telec * boltz)/ev2j          ! Kelvin to ev.

write(6,1001) ' beta(ev^-1): ', beta
write(6,1001) ' telec(ev):   ', telec

ampl = 0.d0
ampr = 0.d0
ampG = 0.d0
rewind(5)
read(5,nml=voltage,end=104)
104 continue

write(6,1001) " left voltage:  ", ampl
write(6,1001) " right voltage: ", ampr
write(6,1001) " gate voltage:  ", ampG
write(6,*)

!----------------------
! get photon namelist !
!----------------------
mode    = 1
epcoup  = 1.d0
phfreq  = 0.05d0
pt_telec=6.d3
rewind(5)
read(5,photon, end=105)
105 continue

pt_beta = ev2j / (boltz * pt_telec)   ! 1/kT in 1/ev  

write(6,'(A,L2)') ' Is photon considered?  ', lphoton
write(6,'(A,L2)') ' Is plasmon considered? ', tplasmon
write(6,1003) ' number of mode:', mode
write(6,1001) ' photon freq:   ', phfreq
write(6,1001) ' e-p coupling:  ', epcoup
write(6,1001) ' photon telec:  ', pt_telec
write(6,*)

if(lphoton.and.tplasmon) then
  write(6,'(A)') "warning: only one of lphoton and tplasmon can be true!"
  stop
endif


if(tPlasmon) then
  barrier_e = 0.0
  barrier_h = 0.0
  eps_coup  = 1.d0
  mscoup    = 1.d0
  relax     = 1.d0
  nst       = 2
  eig(1)     = -1.d0
  eig(2)     = 1.d0
  t_read_file= .false.
  rewind(5)
  read(5,plasmon,end=106)
  106 continue
  metal_ps%barrier_e= barrier_e
  metal_ps%barrier_h= barrier_h
  metal_ps%mscoup   = mscoup
  metal_ps%relax    = relax
  metal_ps%nst      = nst
  write(6,1001) " barrier for electron:      ", metal_ps%barrier_e
  write(6,1001) " barrier for hole:          ", metal_ps%barrier_h
  write(6,1001) " metal-system coupling:     ", metal_ps%mscoup
  write(6,1001) " metal-plasmon coupling:    ", eps_coup
  write(6,1001) " relaxation time:           ", metal_ps%relax
  write(6,1003) " number of states in metal: ", metal_ps%nst
  !
  allocate(metal_ps%eig(nst))
  metal_ps%eig=0.d0
  if(t_read_file) then
    call read_eig(nst,metal_ps%eig,"eig.dat")
  else
    write(6,*) "states of metal"
    do i=1,nst
      if(nst>1) metal_ps%eig(i) = eig(1) + dble(i-1)*(eig(2)-eig(1))/dble(nst-1)
      write(6,'(I5,e15.7)') i,metal_ps%eig(i)
    enddo
  endif
  write(6,*)
  eImag = metal_ps%relax/2.d0
endif


1001 format(A,e15.7)
1003 format(A,I5)

end subroutine io
