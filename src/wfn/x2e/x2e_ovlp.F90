!######################################################################
subroutine x2e_ovlp(rmax, wfnl, wfnr, ovlp)

  use grid_mod, only : ngrid
  use wfn_mod, only : nfun

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfnl(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(in) :: wfnr(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(out) :: ovlp(1:nfun, 1:nfun)

  complex(kind(0d0)), external :: x2e_iprod

  ovlp(1, 1) = x2e_iprod(rmax, wfnl, wfnr)
end subroutine x2e_ovlp
!######################################################################
