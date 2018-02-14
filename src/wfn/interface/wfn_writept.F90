!######################################################################
subroutine wfn_writept(fname, wfn, imethod)

  use io_mod, only : iow

  implicit none
  character(len = *), intent(in) :: fname
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: wfn(*)

  integer :: posx
  integer, external :: wfn_pospt1

  open(unit = iow, file = trim(fname), status = 'unknown', form = 'formatted')

  posx = wfn_pospt1(imethod)
  call x2e_write(iow, wfn(posx))

  close(unit = iow)

end subroutine wfn_writept
!######################################################################
