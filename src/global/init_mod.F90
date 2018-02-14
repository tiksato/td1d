!################################################################################
module init_mod

  implicit none

  integer :: maxcyc
  integer :: maxdav

  character(len=16) :: guess_type
  character(len=16) :: init_type
  character(len=16) :: prop_type_init
  character(len=16) :: diag_type

  logical :: broken_symmetry
  integer :: cic_11(1:2)
  real(kind(0d0)) :: distep, distepx, distepci, distepmo, distep1, distep2, distep3

  integer :: init_ah_couple(1:3)
  integer :: init_ah_xnew
  integer :: init_ah_maxvec
  integer :: init_ah_maxvec2(1:3)
  real(kind(0d0)) :: init_ah_alpha
  real(kind(0d0)) :: init_ah_trust
  logical :: init_ah_norm

  real(kind(0d0)) :: init_diffe
  logical :: guess_interpolate
  logical :: guess_alter

end module init_mod
!################################################################################
