!################################################################################
subroutine wfn_dealloc

  use wfn_mod, only : wfn

  implicit none

  deallocate(wfn)

end subroutine wfn_dealloc
!################################################################################
