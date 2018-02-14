!######################################################################
subroutine x2e_clear(wfn)

  use const_mod, only : czero

  implicit none
  complex(kind(0d0)), intent(out) :: wfn(*)

  integer, external :: x2e_size
  integer :: len

  len = x2e_size()
  wfn(1:len) = czero

end subroutine x2e_clear
!######################################################################
