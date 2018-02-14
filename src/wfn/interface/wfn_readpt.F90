!######################################################################
subroutine wfn_readpt(fname, wfn, imethod)

  use io_mod, only : ior

  implicit none
  character(len = *), intent(in) :: fname
  complex(kind(0d0)), intent(out) :: wfn(*)
  integer, intent(in) :: imethod

  integer :: pospt1
  integer, external :: wfn_pospt1

  open(unit = ior, file = trim(fname), status = 'old', form = 'formatted')

  pospt1 = wfn_pospt1(imethod)
  call x2e_read(ior, wfn(pospt1))

  close(unit = ior)

end subroutine wfn_readpt
!######################################################################
