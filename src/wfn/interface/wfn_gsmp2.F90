!######################################################################
subroutine wfn_gsmp2(wfn, imethod)

  use const_mod, only : zero
  use grid_mod, only : ngrid
  use wfn_mod, only : nspin

  implicit none
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)

end subroutine wfn_gsmp2
!######################################################################
