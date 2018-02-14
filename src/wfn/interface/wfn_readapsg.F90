!################################################################################
subroutine wfn_readapsg(wfnin, wfn, imethod)

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: wfnin(*)
  complex(kind(0d0)), intent(out) :: wfn(*)
  !--------------------------------------------------------------------
  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     stop 'x2e_readapsg nyi.'
  else if (imethod == 0) then
     stop 'hf_readapsg nyi.'
  end if

end subroutine wfn_readapsg
!################################################################################
