!######################################################################
subroutine hf_full_potential(lfield, fock, ng0, ng1)

  use const_mod, only : iunit
  use root_mod, only : icomp
  use grid_mod, only : ngrid, x, v1, docap, mask
  use wfn_mod, only : nspin

  integer, intent(in) :: ng0, ng1
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(inout) :: fock(0:ngrid, 0:ngrid, 1:nspin)

  integer :: igrid, ispin
  complex(kind(0d0)), allocatable :: veff(:)

  allocate(veff(ng0:ng1))

  veff(ng0:ng1) = v1(ng0:ng1) - lfield * x(ng0:ng1)
  if (docap .and. icomp == 1) then
     veff(ng0:ng1) = veff(ng0:ng1) - iunit * mask(ng0:ng1)
  end if

  do ispin = 1, nspin
     do igrid = ng0, ng1
        fock(igrid, igrid, ispin) &
    & = fock(igrid, igrid, ispin) &
    & + veff(igrid)
        end do
  end do

  deallocate(veff)

end subroutine hf_full_potential
!######################################################################
