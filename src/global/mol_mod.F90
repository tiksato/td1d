!################################################################################
module mol_mod

  implicit none

  integer, parameter :: maxatom = 100
  integer :: natom
  integer :: nae
  integer :: nbe
  integer :: ne(3)
  integer :: mult
  real(kind(0d0)) :: z(1:maxatom)
  real(kind(0d0)) :: c(1:maxatom)
  real(kind(0d0)) :: enen

end module mol_mod
!################################################################################
