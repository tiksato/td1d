!######################################################################
subroutine hf_write(iunit, wfn)

  use grid_mod, only : ngrid, x
  use wfn_mod, only : nfun, nspin

  implicit none
  integer, intent(in) :: iunit
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nfun, 1:nspin)

  write(iunit, "( I25)") ngrid
  write(iunit, "( E25.15)") x(0:ngrid)
  write(iunit, "(2E25.15)") wfn(0:ngrid, 1:nfun, 1:nspin)

end subroutine hf_write
!######################################################################
