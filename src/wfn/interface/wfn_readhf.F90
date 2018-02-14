!################################################################################
subroutine wfn_readhf(wfnin, wfn, imethod)

  use wfn_mod, only : ormas

  implicit none
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: wfnin(*)
  complex(kind(0d0)), intent(out) :: wfn(*)

  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     call x2e_readhf(wfnin, wfn)
  else if (imethod == 0) then
     call hf_readhf(wfnin, wfn)
  end if

end subroutine wfn_readhf
!################################################################################
