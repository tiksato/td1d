!################################################################################
subroutine wfn_init

  use wfn_mod, only : wfn, imethod
  use wfn_mod, only : ormas

  implicit none

  call wfn_alloc(imethod)
  call wfn_clear(wfn, imethod)

end subroutine wfn_init
!################################################################################
