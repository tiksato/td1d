!######################################################################
subroutine hf_setfroz(wfn0, wfn, wfnin)

  use root_mod, only : icomp
  use grid_mod, only : ngrid
  use const_mod, only : czero
  use wfn_mod, only : nfun, nspin, sae, crapola, saex, crapx

  implicit none
  complex(kind(0d0)), intent(in)  :: wfn0 (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in)  :: wfn  (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: wfnin(0:ngrid, 1:nfun, 1:nspin)

  saex = sae .and. icomp == 1
  crapx = crapola .and. icomp == 1

  if (saex) then
     wfnin(0:ngrid, 1:nfun, 1) = wfn0(0:ngrid, 1:nfun, 1)
     wfnin(0:ngrid, 1:nfun, 2) = wfn (0:ngrid, 1:nfun, 2)
  else if (crapx) then
     wfnin(0:ngrid, 1:nfun, 1) = wfn0(0:ngrid, 1:nfun, 1)
     wfnin(0:ngrid, 1:nfun, 2) = wfn (0:ngrid, 1:nfun, 2)
  else
     wfnin(0:ngrid, 1:nfun, 1:nspin) = wfn(0:ngrid, 1:nfun, 1:nspin)
  end if

  return

end subroutine hf_setfroz
!######################################################################
