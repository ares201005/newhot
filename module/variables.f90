module variables
use transmod
implicit none

logical :: lphoton
logical :: tplasmon
logical :: lscba
logical :: t_read_file

integer :: mode
integer :: ns, nLead,norbs
integer :: nst

real*8  :: onsite, hop, elcoup,half_w(5)
real*8  :: barrier_e,barrier_h,eps_coup,mscoup,relax,eig(2)
real*8  :: width,dtmp

real*8  :: pt_telec,pt_beta
real*8  :: telec,beta
real*8  :: ampL, ampR, ampG, miu0, ebot
real*8  :: epscoup,epcoup, phfreq         
real*8  :: eImag

type(mps):: metal_ps

end module variables
