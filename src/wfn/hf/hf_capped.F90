!######################################################################
subroutine hf_capped(rmax, wfn, wfnp, q0, q1, q2)

  use grid_mod, only : ngrid
  use wfn_mod, only : nfun, nspin

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  real(kind(0d0)), intent(inout) :: q0, q1, q2
  complex(kind(0d0)), intent(in) :: wfn (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in) :: wfnp(0:ngrid, 1:nfun, 1:nspin)

  real(kind(0d0)) :: p0, p1, p2, ptot  
  real(kind(0d0)), external :: hf_norm

  ptot = hf_norm(rmax, wfn, p0, p1, p2)
  q0 = q0 + p0
  q1 = q1 + p1
  q2 = q2 + p2

  ptot = hf_norm(rmax, wfnp, p0, p1, p2)
  q0 = q0 - p0
  q1 = q1 - p1
  q2 = q2 - p2

end subroutine hf_capped
!######################################################################
