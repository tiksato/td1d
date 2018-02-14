!################################################################################
subroutine wfn_ort(wfn, imethod)

  use wfn_mod, only : ormas

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout)  :: wfn(*)
  !--------------------------------------------------------------------
  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     call x2e_ort(wfn)
  else if (imethod == 0) then
     call hf_ort(wfn)
  end if

end subroutine wfn_ort
!################################################################################
