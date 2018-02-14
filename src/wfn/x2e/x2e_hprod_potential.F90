!######################################################################
subroutine x2e_hprod_potential(lfield, wfn, hwfn, ng0, ng1)

  use const_mod, only : iunit
  use field_mod, only : gauge
  use root_mod, only : icomp, nocoulomb
  use grid_mod, only : ngrid, x, v1, v2, docap, mask

  implicit none
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in) :: wfn (0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 0:ngrid)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, jgrid
  complex(kind(0d0)) :: v12
  complex(kind(0d0)), allocatable :: veff(:)

  allocate(veff(0:ngrid))

  if (trim(gauge) == 'L') then
     veff(0:ngrid) = v1(0:ngrid) - lfield * x(0:ngrid)
  else
     veff(0:ngrid) = v1(0:ngrid)
  end if

  if (docap .and. icomp == 1) then
     veff(0:ngrid) = veff(0:ngrid) - iunit * mask(0:ngrid)
  end if

  if (nocoulomb) then
     do igrid = ng0, ng1
        do jgrid = 1, igrid
           v12 = veff(igrid) + veff(jgrid)
           hwfn(jgrid, igrid) = hwfn(jgrid, igrid) + v12 * wfn(jgrid, igrid)
        end do
     end do
  else
     do igrid = ng0, ng1
        do jgrid = 1, igrid
           v12 = veff(igrid) + veff(jgrid) + v2(jgrid, igrid)
           hwfn(jgrid, igrid) = hwfn(jgrid, igrid) + v12 * wfn(jgrid, igrid)
        end do
     end do
  end if

  deallocate(veff)

end subroutine x2e_hprod_potential
!######################################################################
