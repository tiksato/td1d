!################################################################################
subroutine wfn_exps2(wfn, imethod)

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout)  :: wfn(*)
  !--------------------------------------------------------------------
  integer :: poscic
  integer, external :: wfn_poscic

end subroutine wfn_exps2
!################################################################################
