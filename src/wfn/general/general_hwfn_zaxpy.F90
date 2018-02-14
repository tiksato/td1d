!######################################################################
subroutine general_hwfn_zaxpy(nfun, fac, h1wfn, h2wfn, ng0, ng1)

  use const_mod, only : iunit
  use root_mod, only : icomp
  use grid_mod, only : ngrid, x, v1, docap, mask

  implicit none
  integer, intent(in) :: nfun
  complex(kind(0d0)), intent(in) :: fac
  complex(kind(0d0)), intent(in) :: h1wfn(0:ngrid, 1:*)
  complex(kind(0d0)), intent(inout) :: h2wfn(0:ngrid, 1:*)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, ifun

  do ifun = 1, nfun
     do igrid = ng0, ng1
        h2wfn(igrid, ifun) &
    & = h2wfn(igrid, ifun) &
    & + h1wfn(igrid, ifun) * fac
     end do
  end do

end subroutine general_hwfn_zaxpy
!######################################################################
