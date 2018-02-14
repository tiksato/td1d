!######################################################################
subroutine x2e_hprod_kinetic(lfield, wfn, hwfn, ng0, ng1)

  use const_mod, only : two, half
  use grid_mod, only : ngrid, dgrid

  implicit none
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in)    :: wfn (0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 0:ngrid)
  integer, intent(in) :: ng0, ng1

!old  integer :: igrid, jgrid
!old  complex(kind(0d0)) :: tmp
!old  real(kind(0d0)) :: fac

  call general_hprod_kinetic_2d(lfield, wfn, hwfn, ng0, ng1)

!old  fac = -half / (dgrid * dgrid)
!old
!old  do igrid = ng0, ng1
!old     do jgrid = 1, igrid
!old        tmp = (wfn(jgrid - 1, igrid) &
!old          & -  wfn(jgrid,     igrid) * two &
!old          & +  wfn(jgrid + 1, igrid) &
!old          & +  wfn(jgrid, igrid - 1) &
!old          & -  wfn(jgrid, igrid    ) * two &
!old          & +  wfn(jgrid, igrid + 1)) * fac
!old        hwfn(jgrid, igrid) = hwfn(jgrid, igrid) + tmp
!old     end do
!old  end do

end subroutine x2e_hprod_kinetic
!######################################################################
