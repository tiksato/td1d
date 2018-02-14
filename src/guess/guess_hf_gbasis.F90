!################################################################################
subroutine guess_hf_gbasis(wfn, imethod)

  use grid_mod, only : ngbas,g09out
  use wfn_mod, only : nfun,nspin

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(out) :: wfn(1:*)
  !--------------------------------------------------------------------
  integer, external :: wfn_size
  integer, external :: wfn_poscic

  call guess_read_gbasis(wfn, imethod)

end subroutine guess_hf_gbasis
!################################################################################
