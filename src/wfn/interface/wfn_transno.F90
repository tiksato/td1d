!################################################################################
subroutine wfn_transno(wfn, imethod)

  use wfn_mod, only : ormas

  implicit none
  complex(kind(0d0)), intent(out) :: wfn(*)
  integer, intent(in) :: imethod

  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
  else if (imethod == 0) then
  end if

end subroutine wfn_transno
!################################################################################
