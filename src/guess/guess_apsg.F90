!################################################################################
subroutine guess_apsg(wfn, imethod)

  ! read apsg and ...
  ! X2E: nyi,
  ! HF: nyi,
  ! GVB: nyi,
  ! CAS: nyi,
  ! APSG: nyi,
  ! CIS: nyi.
  ! MRMP: nyi.

  use root_mod, only : name

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(out) :: wfn(*)
  !--------------------------------------------------------------------
  character(len = 256) :: fname
  integer, parameter :: iapsg = 3
  integer :: len
  integer, external :: wfn_size
  complex(kind(0d0)), allocatable :: wfnin(:)

  len = wfn_size(iapsg)
  allocate(wfnin(1:len))
  call wfn_clear(wfnin, iapsg)

  fname = trim(name)//".wfnin"
  call wfn_read(fname, wfnin, iapsg)
  call wfn_readapsg(wfnin, wfn, imethod)

  deallocate(wfnin)

end subroutine guess_apsg
!################################################################################
