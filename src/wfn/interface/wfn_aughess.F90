!################################################################################
real(kind(0d0)) function wfn_aughess(lfield, accept, trad, wfn, imethod)

  implicit none
  !--------------------------------------------------------------------
  real(kind(0d0)), intent(in) :: lfield
  logical, intent(out) :: accept
  real(kind(0d0)), intent(inout) :: trad
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------

end function wfn_aughess
!################################################################################
