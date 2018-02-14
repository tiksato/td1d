!######################################################################
subroutine general_hprod_op1e(iop, wfn, hwfn, ng0, ng1)

  use const_mod, only : iunit
  use root_mod, only : icomp
  use grid_mod, only : ngrid, x, v1, gv1
  use wfn_mod, only : nfcore, nfun

  integer, intent(in) :: iop
  complex(kind(0d0)), intent(in)    :: wfn (0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:nfun)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, ifun
  complex(kind(0d0)), allocatable :: veff(:)

  allocate(veff(ng0:ng1))

  if (iop == 0) then
     ! dipole
     veff(ng0:ng1) = x(ng0:ng1)
  else if (iop == 1) then
     ! acceleration
     veff(ng0:ng1) = gv1(ng0:ng1)
  end if

  do ifun = 1, nfun
     do igrid = ng0, ng1
        hwfn(igrid, ifun) = hwfn(igrid, ifun) + wfn(igrid, ifun) * veff(igrid)
     end do
  end do

  deallocate(veff)

end subroutine general_hprod_op1e
!######################################################################
