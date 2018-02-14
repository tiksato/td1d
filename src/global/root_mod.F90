!################################################################################
module root_mod

  implicit none

  character(len = 256) :: name
  integer :: ical
  integer :: iprint
  integer :: isoftnuc
  integer :: isoftr12
  real(kind(0d0)) :: softnuc
  real(kind(0d0)) :: softr12
  integer :: icomp
  logical :: nocoulomb

  integer :: nprint_ene
  integer :: nprint_op1e
  integer :: nprint_op1x
  integer :: nprint_ionp
  integer :: nprint_rho1
  integer :: nprint_wfn

  logical :: root_debug
  logical :: use_blas

end module root_mod
!################################################################################
