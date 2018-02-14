!######################################################################
subroutine x2e_ort(wfn)

  use const_mod, only : one
  use grid_mod, only : ngrid, dgrid

  implicit none
  complex(kind(0d0)), intent(inout) :: wfn(0:ngrid, 0:ngrid)

  real(kind(0d0)) :: norm, fac
  complex(kind(0d0)), external :: x2e_iprod

  norm = sqrt(real(x2e_iprod(-one, wfn, wfn)))
  fac = one / norm
  wfn(0:ngrid, 0:ngrid) = wfn(0:ngrid, 0:ngrid) * fac

end subroutine x2e_ort
!######################################################################
