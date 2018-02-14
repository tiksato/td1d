!######################################################################
subroutine hf_v2prod(dtime, wfn, hwfn)

  use omp_mod
  use mol_mod, only : ne
  use const_mod, only :czero
  use grid_mod, only : ngrid, gll, gul
  use wfn_mod, only : nfcore, nfun, nspin, hf_doproj

  implicit none
  real(kind(0d0)), intent(in) :: dtime
  complex(kind(0d0)), intent(in)  :: wfn (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: hwfn(0:ngrid, 1:nfun, 1:nspin)

  hwfn = czero
  if(ne(3) == 1) return

  !$omp parallel default(shared)
  call omp_mod_thread(gll, gul)
  call hf_hprod_meanfield(wfn, wfn, hwfn, ng0, ng1)
  !$omp end parallel
  if (hf_doproj) call hf_hprod_proj(dtime, wfn, hwfn)

end subroutine hf_v2prod
!######################################################################
