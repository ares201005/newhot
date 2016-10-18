module parameters
implicit none

real*8, parameter :: pi=3.141592653589793238d0
real*8, parameter :: auev=27.2113845d0
real*8, parameter :: hbarinv=1.51924888335207073622d0
real*8, parameter :: hbar=0.65822d0
real*8, parameter :: qe = 1.6022d-19
integer, parameter :: mgauleg1=8

double precision, parameter :: boltz  = 1.3806506d-23
double precision, parameter :: ev2j   = 1.60210d-19
double precision, parameter :: dsmall = 1.D-8
!
real*8, parameter :: freqau=1.378999779d6
real*8, parameter :: aufs=2.418884326505d-2
complex*16, parameter :: cunity=(1.d0, 0.d0), czero=(0.d0, 0.d0), eye=(0.d0, 1.d0)

double precision :: dpt(mgauleg1), dw(mgauleg1)
data dpt    /                &
             -0.96028986d0,  & 
             -0.79666648d0,  &
             -0.52553241d0,  &
             -0.18343464d0,  &
              0.18343464d0,  &
              0.52553241d0,  &
              0.79666648d0,  &
              0.96028986d0 /
data dw     /                &
              0.10122854d0,  &  
              0.22238103d0,  & 
              0.31370665d0,  & 
              0.36268378d0,  & 
              0.36268378d0,  & 
              0.31370665d0,  & 
              0.22238103d0,  & 
              0.10122854d0 /


end module parameters
