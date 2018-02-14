!######################################################################
subroutine wfn_write_orb(fname, wfn, imethod)

  use io_mod, only : iow
  use wfn_mod, only : ormas

  implicit none
  character(len = *), intent(in) :: fname
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: wfn(*)

  open(unit = iow, file = trim(fname), status = 'unknown', form = 'formatted')

  if (imethod == -1) then
  else if (imethod == 0) then
     call hf_write_orb(iow, wfn)
  end if

  close(unit = iow)

end subroutine wfn_write_orb
!######################################################################
