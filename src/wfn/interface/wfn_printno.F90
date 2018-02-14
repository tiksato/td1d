!######################################################################
subroutine wfn_printno(fname, wfn, imethod)

  use io_mod, only : iow

  implicit none
  character(len = *), intent(in) :: fname
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: wfn(*)

  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)
  open(unit = iow, file = trim(fname), status = 'unknown', form = 'formatted')

  if (imethod == -1) then
     call x2e_printno(iow, wfn)
  else if (imethod == 0) then
  end if

  close(unit = iow)

end subroutine wfn_printno
!######################################################################
