!######################################################################
subroutine general_hprod_proj(wfn, hwfn)

  use const_mod, only : one
  use grid_mod, only : ngrid
  use wfn_mod, only : nfcore, nfun

  implicit none
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:nfun)

  integer :: ifun, jfun, nbas
  complex(kind(0d0)) :: fac
  complex(kind(0d0)), allocatable :: wfn_dot_hwfn(:,:)
  

  nbas = ngrid + 1
  allocate(wfn_dot_hwfn(1:nfun, 1:nfun))
  call general_ovlp(.true., -one, wfn, hwfn, wfn_dot_hwfn)

  do ifun = nfcore + 1, nfun
     do jfun = 1, nfun
        fac = - wfn_dot_hwfn(jfun, ifun)
        call zaxpy(nbas, fac, wfn(0, jfun), 1, hwfn(0, ifun), 1)
!old        hwfn(0:ngrid, ifun) &
!old    & = hwfn(0:ngrid, ifun) &
!old    & -  wfn(0:ngrid, jfun) * wfn_dot_hwfn(jfun, ifun)
     end do
  end do

  deallocate(wfn_dot_hwfn)

end subroutine general_hprod_proj
!######################################################################
