!################################################################################
subroutine hf_readhf(wfnin, wfn)

  use init_mod, only : guess_type
  use grid_mod, only : ngrid, x
  use wfn_mod, only : nfun, nspin

  implicit none
  complex(kind(0d0)), intent(in) :: wfnin(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: wfn(0:ngrid, 1:nfun, 1:nspin)

  call util_zcopy((ngrid+1)*nfun*nspin, wfnin, 1, wfn, 1)

end subroutine hf_readhf
!################################################################################
