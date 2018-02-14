!################################################################################
subroutine wfn_final

  use wfn_mod, only : imethod
  use wfn_mod, only : ormas

  implicit none

  call wfn_dealloc

end subroutine wfn_final
!################################################################################
