module transmod
!
complex*16, allocatable, save :: newhmat(:,:)
real*8, allocatable, save     :: deltah(:,:)
real*8, allocatable, save     :: LamdaL(:,:), LamdaR(:,:)
real*8, allocatable, save     :: deltapt(:, :)

real*8, allocatable, save :: pheigv(:) !phonon frequency
real*8, allocatable, save :: bocc(:)     ! phonon occupation number

complex*16, allocatable, save :: coup(:,:,:)
complex*16, allocatable, save :: sigma0(:,:)
!
!* external field 
real*8, save :: freq(5)
real*8, save :: ex, ey, ez

!*
real*8, save :: charge

type :: dBlkSeq
  real*8, allocatable :: g0(:,:)
  real*8, allocatable :: ge(:,:)
  real*8 :: width
end type dBlkSeq
!

type :: mps
  real*8 :: barrier_e  ! barrier for electron
  real*8 :: barrier_h  ! barrier for hole
  real*8 :: mscoup     ! metal-system coupling
  real*8 :: relax      ! relaxation time of metal
  integer:: nst        ! number of states in the metallic contatcts
  !
  real*8, allocatable :: eig(:) ! energy of states
end type mps

end module transmod
