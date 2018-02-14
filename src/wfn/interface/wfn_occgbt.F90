!################################################################################
real(kind(0d0)) function wfn_occgbt(calene, lfield, wfn0, wfn, imethod)

  use const_mod, only : zero
  use root_mod, only : icomp
  use wfn_mod, only : ormas

  implicit none
  !--------------------------------------------------------------------
  logical, intent(in) :: calene
  real(kind(0d0)), intent(in) :: lfield
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in)  :: wfn0(*), wfn(*)
  !--------------------------------------------------------------------

end function wfn_occgbt
!################################################################################
