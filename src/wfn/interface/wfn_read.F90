!################################################################################
subroutine wfn_read(fname, wfn, imethod)

  use io_mod, only : ior
  use wfn_mod, only : ormas

  implicit none
  character(len = *), intent(in) :: fname
  complex(kind(0d0)), intent(out) :: wfn(*)
  integer, intent(in) :: imethod

  integer :: poscic
  integer, external :: wfn_poscic

  open(unit = ior, file = trim(fname), status = 'old', form = 'formatted')

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     call x2e_read(ior, wfn)
  else if (imethod == 0) then
     call hf_read(ior, wfn)
  end if

  close(unit = ior)

end subroutine wfn_read
!################################################################################
subroutine wfn_read_fcore(fname, wfn, imethod)

  use io_mod, only : ior
  use wfn_mod, only : ormas, nfcore, ncore

  implicit none
  character(len = *), intent(in) :: fname
  complex(kind(0d0)), intent(out) :: wfn(*)
  integer, intent(in) :: imethod

end subroutine wfn_read_fcore
!################################################################################
