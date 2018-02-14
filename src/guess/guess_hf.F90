!################################################################################
subroutine guess_hf(wfn, imethod)

  use root_mod, only : name

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(out) :: wfn(*)
  !--------------------------------------------------------------------
  character(len = 256) :: fname
  integer, parameter :: ihf = 0
  integer :: len
  integer, external :: wfn_size
  complex(kind(0d0)), allocatable :: wfnin(:)

  len = wfn_size(ihf)
  allocate(wfnin(1:len))
  call wfn_clear(wfnin, ihf)

  fname = trim(name)//".wfnin"
  call wfn_read(fname, wfnin, ihf)
  call wfn_readhf(wfnin, wfn, imethod)

  deallocate(wfnin)

end subroutine guess_hf
!################################################################################
