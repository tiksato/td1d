!################################################################################
subroutine guess_x2e(wfn, imethod)

  ! read CAS and ...
  ! X2E: nyi,
  ! HF: nyi,
  ! GVB: nyi,
  ! CAS: nyi,
  ! APSG: nyi,
  ! CIS: nyi,
  ! MRMP: nyi,

  use root_mod, only : name

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(out) :: wfn(*)
  !--------------------------------------------------------------------
  character(len = 256) :: fname
  integer, parameter :: ix2e = -1
  integer :: len
  integer, external :: wfn_size
  complex(kind(0d0)), allocatable :: wfnin(:)

  len = wfn_size(ix2e)
  allocate(wfnin(1:len))
  call wfn_clear(wfnin, ix2e)

  fname = trim(name)//".wfnin"
  call wfn_read(fname, wfnin, ix2e)
  call wfn_readx2e(wfnin, wfn, imethod)

  deallocate(wfnin)

end subroutine guess_x2e
!################################################################################
