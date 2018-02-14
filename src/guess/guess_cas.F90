!################################################################################
subroutine guess_cas(wfn, imethod)

  ! read CAS and ...
  ! X2E: expanded to 2e wfn,
  ! HF: nyi,
  ! GVB: nyi,
  ! CAS: do nothing,
  ! APSG: diagonalize CAS 1RDM to make natural expansion,
  ! CIS: nyi,
  ! MRMP: nyi,

  use root_mod, only : name

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(out) :: wfn(*)
  !--------------------------------------------------------------------
  character(len = 256) :: fname
  integer, parameter :: icas = 2
  integer :: len
  integer, external :: wfn_size
  complex(kind(0d0)), allocatable :: wfnin(:)

  len = wfn_size(icas)
  allocate(wfnin(1:len))
  call wfn_clear(wfnin, icas)

  fname = trim(name)//".wfnin"
  call wfn_read(fname, wfnin, icas)
  call wfn_readcas(wfnin, wfn, imethod)

  deallocate(wfnin)

end subroutine guess_cas
!################################################################################
