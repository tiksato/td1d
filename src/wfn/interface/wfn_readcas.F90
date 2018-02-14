!################################################################################
subroutine wfn_readcas(wfnin, wfn, imethod)

  implicit none
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: wfnin(*)
  complex(kind(0d0)), intent(out) :: wfn(*)

  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     stop 'sorry, x2e_readcas nyi.'
  else if (imethod == 0) then
     stop 'hf_readcas nyi.'
  end if

end subroutine wfn_readcas
!################################################################################
