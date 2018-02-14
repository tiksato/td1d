!######################################################################
subroutine general_full_kinetic(fock, ng0, ng1)

  use const_mod, only : czero
  use grid_mod, only : ngrid, dgrid, gll, gul
  use grid_mod, only : fd_coeff_in, fd_ohalf

  implicit none
  integer, intent(in) :: ng0, ng1
  complex(kind(0d0)), intent(inout) :: fock(0:ngrid, 0:ngrid)

  integer :: igrid, jgrid, sigma

  do igrid = ng0, ng1
     do jgrid = gll, gul
        sigma = igrid - jgrid
        if (abs(sigma) <= fd_ohalf) then
           fock(jgrid, igrid) = fock(jgrid, igrid) + fd_coeff_in(sigma)
        end if
     end do
  end do

end subroutine general_full_kinetic
!######################################################################
