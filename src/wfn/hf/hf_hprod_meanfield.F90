!######################################################################
subroutine hf_hprod_meanfield(wfnin, wfn, hwfn, ng0, ng1)

  use mol_mod, only : ne
  use root_mod, only : icomp, nocoulomb
  use grid_mod, only : ngrid
  use wfn_mod, only : nfun, nspin, saex, crapx, hcore, ldacx

  implicit none
  complex(kind(0d0)), intent(in)    :: wfnin(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in)    :: wfn  (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(inout) :: hwfn (0:ngrid, 1:nfun, 1:nspin)
  integer, intent(in) :: ng0, ng1

  if (icomp < 0 .or. ne(3) <= 1 .or. nocoulomb .or. hcore) return

  if (saex .or. crapx) then
     call hf_hprod_meanfield_nosi(wfnin, wfn, hwfn, ng0, ng1)
  else if (ldacx) then
     stop 'lda not yet implemented.'
!nyi     call hf_hprod_meanfield_def(wfnin, wfn, hwfn, ng0, ng1)
!nyi     call hf_hprod_meanfield_lda(wfnin, wfn, hwfn, ng0, ng1)
  else
     call hf_hprod_meanfield_def(wfnin, wfn, hwfn, ng0, ng1)
  end if

end subroutine hf_hprod_meanfield
!######################################################################
