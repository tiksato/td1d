!################################################################################
subroutine wfn_print(iunit, wfn, imethod)

  use wfn_mod, only : ormas

  implicit none
  integer, intent(in) :: iunit, imethod
  complex(kind(0d0)), intent(in) :: wfn(*)

  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     call x2e_print(iunit, wfn)
  else if (imethod == 0) then
     call hf_print(iunit, wfn)
  end if

end subroutine wfn_print
!################################################################################
