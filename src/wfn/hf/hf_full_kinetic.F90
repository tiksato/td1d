!######################################################################
subroutine hf_full_kinetic(fock, ng0, ng1)

  use grid_mod, only : ngrid
  use wfn_mod, only : nspin

  implicit none
  integer, intent(in) :: ng0, ng1
  complex(kind(0d0)), intent(inout) :: fock(0:ngrid, 0:ngrid, 1:nspin)
  integer :: ispin

  do ispin = 1, nspin
     call general_full_kinetic(fock(0,0,ispin), ng0, ng1)
  end do

!old  fac = -half / ( dgrid * dgrid )
!old
!old  do ispin = 1, nspin
!old     do ifun = 1, nfun
!old        do igrid = ng0, ng1
!old           hwfn(igrid, ifun, ispin) = hwfn(igrid, ifun, ispin) &
!old    &    +      (wfn(igrid-1, ifun, ispin) &
!old    &    - two * wfn(igrid,   ifun, ispin) &
!old    &    +       wfn(igrid+1, ifun, ispin)) * fac
!old        end do
!old     end do
!old  end do

end subroutine hf_full_kinetic
!######################################################################
