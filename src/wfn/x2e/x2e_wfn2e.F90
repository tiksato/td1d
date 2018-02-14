!######################################################################
subroutine x2e_wfn2e(wfn, wfn2e)

  use grid_mod, only : ngrid

  implicit none
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(out) :: wfn2e(0:ngrid, 0:ngrid)

  wfn2e(0:ngrid, 0:ngrid) = wfn(0:ngrid, 0:ngrid)

end subroutine x2e_wfn2e
!######################################################################
