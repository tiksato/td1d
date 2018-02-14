!################################################################################
subroutine general_meanfield(norb1, norb2, orb1, orb2, veff, ng0, ng1)

  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, v2

  implicit none
  integer, intent(in) :: norb1, norb2, ng0, ng1
  complex(kind(0d0)), intent(in) :: orb1(0:ngrid, 1:*)
  complex(kind(0d0)), intent(in) :: orb2(0:ngrid, 1:*)
  complex(kind(0d0)), intent(inout) :: veff(0:ngrid, 1:norb1, 1:norb2)

  integer :: iorb, jorb, igrid, jgrid
  complex(kind(0d0)) :: rho

  do jorb = 1, norb2
     do iorb = 1, norb1
        do jgrid = 0, ngrid
           rho = conjg(orb1(jgrid, iorb)) * orb2(jgrid, jorb) * dgrid
           if (abs(rho) > thrwfn) then
              do igrid = ng0, ng1
                 veff(igrid, iorb, jorb) = veff(igrid, iorb, jorb) + rho * v2(igrid, jgrid)
              end do
           end if
        end do
     end do
  end do

end subroutine general_meanfield
!################################################################################
