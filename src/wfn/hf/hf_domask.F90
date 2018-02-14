!######################################################################
subroutine hf_domask(rmax, dt, wfn, q0, q1, q2)

  use omp_mod
  use grid_mod, only : ngrid, x, xmask, mask, gll, gul
  use wfn_mod, only : nfun, nspin
  use const_mod, only : half

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: dt
  real(kind(0d0)), intent(inout) :: q0, q1, q2
  complex(kind(0d0)), intent(inout) :: wfn(0:ngrid, 1:nfun, 1:nspin)

  integer :: igrid, ifun, ispin
  real(kind(0d0)) :: p0, p1, p2, ptot  
  complex(kind(0d0)) :: dt2
  real(kind(0d0)), external :: hf_norm

  dt2 = dt * half

  ptot = hf_norm(rmax, wfn, p0, p1, p2)
  q0 = q0 + p0
  q1 = q1 + p1
  q2 = q2 + p2

!$omp parallel default(shared) private(igrid, ispin, ifun)

  call omp_mod_thread(gll, gul)

  do igrid = ng0, ng1
     if(abs(x(igrid)) > xmask) then
        do ispin = 1, nspin
           do ifun = 1, nfun
!              wfn(igrid, ifun, ispin) = wfn(igrid, ifun, ispin) * mask(igrid) * dt2
              wfn(igrid, ifun, ispin) = wfn(igrid, ifun, ispin) * mask(igrid)
           end do
        end do
     end if
  end do

!$omp end parallel

  ptot = hf_norm(rmax, wfn, p0, p1, p2)
  q0 = q0 - p0
  q1 = q1 - p1
  q2 = q2 - p2

end subroutine hf_domask
!######################################################################
