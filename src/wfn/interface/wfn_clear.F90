!################################################################################
subroutine wfn_clear(wfn, imethod)

  use const_mod, only : czero

  implicit none
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(out) :: wfn(*)

  integer :: len
  integer, external :: wfn_size

  len = wfn_size(imethod)
  call util_zcopy(len, czero, 0, wfn, 1)

end subroutine wfn_clear
!################################################################################
