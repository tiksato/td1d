!######################################################################
subroutine hf_clear(wfn)

  use const_mod, only : czero

  implicit none
  complex(kind(0d0)), intent(out) :: wfn(*)

  integer, external :: hf_size
  integer :: len

  len = hf_size()
  wfn(1:len) = czero

end subroutine hf_clear
!######################################################################
