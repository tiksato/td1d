!######################################################################
subroutine hf_exp_potential(dstep, lfield, wfn, expwfn)

  use omp_mod
  use const_mod, only : iunit
  use root_mod, only : icomp
  use grid_mod, only : ngrid, x, v1, docap, mask, gll, gul
  use wfn_mod, only : nfun, nspin

  implicit none
  complex(kind(0d0)), intent(in)  :: dstep
  real(kind(0d0)), intent(in)     :: lfield
  complex(kind(0d0)), intent(in)  :: wfn   (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: expwfn(0:ngrid, 1:nfun, 1:nspin)

  integer :: igrid, ifun, ispin
  complex(kind(0d0)), allocatable :: exp_veff(:)

  allocate(exp_veff(0:ngrid))

!$omp parallel default(shared) private(igrid, ifun, ispin)

  call omp_mod_thread(gll, gul)

  exp_veff(ng0:ng1) = v1(ng0:ng1) - lfield * x(ng0:ng1)
  if (docap .and. icomp == 1) then
     exp_veff(ng0:ng1) = exp_veff(ng0:ng1) - iunit * mask(ng0:ng1)
  end if

  do igrid = ng0, ng1
     exp_veff(igrid) = exp(-iunit * dstep * exp_veff(igrid))
  end do

  do ispin = 1, nspin
     do ifun = 1, nfun
        do igrid = ng0, ng1
           expwfn(igrid, ifun, ispin) = exp_veff(igrid) * wfn(igrid, ifun, ispin)
        end do
     end do
  end do

!$omp end parallel

  deallocate(exp_veff)

end subroutine hf_exp_potential
!######################################################################
