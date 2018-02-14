!######################################################################
subroutine x2e_capped(rmax, wfn, wfnp, q0, q1, q2)

  use grid_mod, only : ngrid

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  real(kind(0d0)), intent(inout) :: q0, q1, q2
  complex(kind(0d0)), intent(in) :: wfn (0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(in) :: wfnp(0:ngrid, 0:ngrid)

  real(kind(0d0)) :: p0, p1, p2, ptot  
  real(kind(0d0)), external :: x2e_norm

  ptot = x2e_norm(rmax, wfn, p0, p1, p2)
  q0 = q0 + p0
  q1 = q1 + p1
  q2 = q2 + p2

  ptot = x2e_norm(rmax, wfnp, p0, p1, p2)
  q0 = q0 - p0
  q1 = q1 - p1
  q2 = q2 - p2

end subroutine x2e_capped
!######################################################################
