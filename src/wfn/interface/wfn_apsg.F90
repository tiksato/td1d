!################################################################################
subroutine wfn_apsg(wfn, imethod)

  use wfn_mod, only : nblock, type_block

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  integer :: poscic
  integer, external :: wfn_poscic

  if (maxval(type_block) > 2) stop 'apsg: only geminals !'

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     stop 'x2e_apsg nyi.'
  else if (imethod == 0) then
     stop 'hf_apsg nyi.'
  end if

end subroutine wfn_apsg
!################################################################################
