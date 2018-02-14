!################################################################################
subroutine wfn_ort_mo(wfn, imethod)

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
  else if (imethod == 0) then
     call hf_ort(wfn)
  end if

end subroutine wfn_ort_mo
!################################################################################
