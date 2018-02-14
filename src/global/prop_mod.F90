!################################################################################
module prop_mod

  implicit none

  character(len=16) :: prop_type
  integer :: cyctot
  integer :: nstep
  integer :: totstep
  real(kind(0d0)) :: dstep
  real(kind(0d0)) :: dstep_ll, dstep_ul

  integer :: inistep
  complex(kind(0d0)) :: initime

  real(kind(0d0)) :: prop_tol, prop_safety, prop_ulscal
  real(kind(0d0)) :: rk5_a(2:6)
  real(kind(0d0)) :: rk5_b(1:5, 2:6)
  real(kind(0d0)) :: rk5_crk5(1:6)
  real(kind(0d0)) :: rk5_crk4(1:6)

end module prop_mod
!################################################################################
!old!################################################################################
!oldmodule prop_mod
!old
!old  implicit none
!old
!old  character(len=16) :: prop_type
!old  integer :: cyctot
!old  integer :: nstep
!old  integer :: totstep
!old  real(kind(0d0)) :: dstep
!old
!old  integer :: inistep
!old  real(kind(0d0)) :: initime
!old
!old  real(kind(0d0)) :: prop_tol, prop_safety, prop_ulscal
!old  real(kind(0d0)), parameter :: rk5_a(2:6) = &
!old       & (/ 1./5.,  &
!old       &    3./10., &
!old       &    3./5.,  &
!old       &    1.,     &
!old       &    7./8. /)
!old  real(kind(0d0)), parameter :: rk5_b(1:5, 2:6) = &
!old       & (/ 1./5.,       0.,        0.,            0.,           0., &
!old       &    3./40.,      9./40.,    0.,            0.,           0., &
!old       &    3./10.,     -9./10.,    6./5.,         0.,           0., &
!old       &  -11./54.,      5./2.,   -70./27.,       35./27,        0., &
!old       & 1631./55296., 175./512., 575./13824., 44275./110592., 253./4096. /)
!old
!old  real(kind(0d0)), parameter :: rk5_crk5(1:6) = &
!old       & (/ 37./378., &
!old       &     0.,      &
!old       &   250./621., &
!old       &   125./594., &
!old       &     0.,      &
!old       &   512./1771. /)
!old  real(kind(0d0)), parameter :: rk5_crk4(1:6) = &
!old       & (/ 2825./27648., &
!old       &       0.,        &
!old       &   18575./48384., &
!old       &   13525./55296., &
!old       &     277./14336., &
!old       &       1./4. /)
!old
!oldend module prop_mod
!old!################################################################################
