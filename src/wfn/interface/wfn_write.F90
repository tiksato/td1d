!################################################################################
subroutine wfn_write(fname, wfn, imethod)

  use io_mod, only : iow
  use wfn_mod, only : ormas

  implicit none
  character(len = *), intent(in) :: fname
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: wfn(*)

  integer :: poscic
  integer, external :: wfn_poscic

  open(unit = iow, file = trim(fname), status = 'unknown', form = 'formatted')

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     call x2e_write(iow, wfn)
  else if (imethod == 0) then
     call hf_write(iow, wfn)
  end if

  close(unit = iow)

end subroutine wfn_write
!################################################################################
