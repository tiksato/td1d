!######################################################################
subroutine general_hwfn_potential(lfield, nfun, wfn, hwfn, ng0, ng1)

  use const_mod, only : iunit
  use root_mod, only : icomp
  use grid_mod, only : ngrid, x, v1, docap, mask

  real(kind(0d0)), intent(in) :: lfield
  integer, intent(in) :: nfun
  complex(kind(0d0)), intent(in) :: wfn (0:ngrid, 1:*)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:*)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, ifun
  complex(kind(0d0)), allocatable :: veff(:)
!debug
!  write(6, "('WARNING: skip 1e-potential operator!')")
!  return
!debug

  allocate(veff(ng0:ng1))

  do igrid = ng0, ng1
     veff(igrid) = v1(igrid) - lfield * x(igrid)
  end do
  if (docap .and. icomp == 1) then
     do igrid = ng0, ng1
        veff(igrid) = veff(igrid) - iunit * mask(igrid)
     end do
  end if

!f90  veff(ng0:ng1) = v1(ng0:ng1) - lfield * x(ng0:ng1)
!f90  if (docap .and. icomp == 1) then
!f90     veff(ng0:ng1) = veff(ng0:ng1) - iunit * mask(ng0:ng1)
!f90  end if

  do ifun = 1, nfun
     do igrid = ng0, ng1
        hwfn(igrid, ifun) &
    & = hwfn(igrid, ifun) &
    & +  wfn(igrid, ifun) * veff(igrid)
     end do
  end do

  deallocate(veff)

end subroutine general_hwfn_potential
!######################################################################
