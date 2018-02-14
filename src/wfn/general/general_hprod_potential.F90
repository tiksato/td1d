!######################################################################
subroutine general_hprod_potential(lfield, wfn, hwfn, ng0, ng1)

  use const_mod, only : iunit
  use field_mod, only : gauge
  use root_mod, only : icomp
  use grid_mod, only : ngrid, x, v1, docap, mask
  use wfn_mod, only : nfcore, nfun

  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in)    :: wfn (0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:nfun)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, ifun
  complex(kind(0d0)), allocatable :: veff(:)

  allocate(veff(ng0:ng1))

  ! nuclear attraction energy
  veff(ng0:ng1) = v1(ng0:ng1)

  ! length gauge field
  if (trim(gauge) == 'L') then
     veff(ng0:ng1) = veff(ng0:ng1) - lfield * x(ng0:ng1)
  end if

  ! complex absorbing potential
  if (docap .and. icomp == 1) then
     veff(ng0:ng1) = veff(ng0:ng1) - iunit * mask(ng0:ng1)
  end if

  do ifun = nfcore + 1, nfun
     do igrid = ng0, ng1
        hwfn(igrid, ifun) &
    & = hwfn(igrid, ifun) &
    & +  wfn(igrid, ifun) * veff(igrid)
     end do
  end do

  deallocate(veff)

end subroutine general_hprod_potential
!######################################################################
subroutine general_hprod_potential2(h1, wfn, hwfn, ng0, ng1)

  use const_mod, only : iunit
  use field_mod, only : gauge
  use root_mod, only : icomp
  use grid_mod, only : ngrid
  use wfn_mod, only : nfcore, nfun

  real(kind(0d0)), intent(in)    :: h1(0:ngrid)
  complex(kind(0d0)), intent(in)    :: wfn(0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:nfun)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, ifun

  do ifun = 1, nfun
     do igrid = ng0, ng1
        hwfn(igrid, ifun) = hwfn(igrid, ifun) + wfn(igrid, ifun) * h1(igrid)
     end do
  end do

end subroutine general_hprod_potential2
!######################################################################
