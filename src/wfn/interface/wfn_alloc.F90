!################################################################################
subroutine wfn_alloc(imethod)

  use wfn_mod, only : wfn, nfun

  implicit none
  integer, intent(in) :: imethod
  integer, external :: wfn_size
  integer :: len

  len = wfn_size(imethod)
  allocate(wfn(len))

end subroutine wfn_alloc
!################################################################################
