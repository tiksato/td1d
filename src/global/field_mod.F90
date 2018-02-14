!################################################################################
module field_mod

  implicit none

  character(len=16) :: env_type
  character(len=1) :: gauge

  integer :: ncyc_on
  integer :: ncyc_off
  integer :: ncyc_flat
  integer :: sin2_ncyc
  real(kind(0d0)) :: fwhm
  real(kind(0d0)) :: tau
  real(kind(0d0)) :: cep
  real(kind(0d0)) :: fint
  real(kind(0d0)) :: famp
  real(kind(0d0)) :: wlen
  real(kind(0d0)) :: freq
  real(kind(0d0)) :: cyc0
  real(kind(0d0)) :: period
  real(kind(0d0)) :: rquiv
  real(kind(0d0)) :: pulse_t
  real(kind(0d0)) :: pulse_dt

end module field_mod
!################################################################################
