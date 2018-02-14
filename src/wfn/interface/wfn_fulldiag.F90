!################################################################################
subroutine wfn_fulldiag(wfn, imethod)

  use const_mod, only : zero
  use wfn_mod, only : nspin

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  integer :: poscic, lenx
  integer, external :: wfn_poscic, wfn_size
  complex(kind(0d0)), allocatable :: work(:)

  lenx = wfn_size(-1)
  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     stop 'x2e_fulldiag nyi'
  else if (imethod == 0) then
     allocate(work(lenx*nspin))
     call hf_fulldiag(zero, wfn, work)
     deallocate(work)
  end if

end subroutine wfn_fulldiag
!################################################################################
