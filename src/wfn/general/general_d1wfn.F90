!######################################################################
subroutine general_d1wfn(wfn, d1wfn, ng0, ng1)

  use grid_mod, only : ngrid, fd_ohalf, fd_coeff1_in, fd_coeff1_l, fd_coeff1_r
  use wfn_mod, only : nfun

  implicit none
  integer, intent(in) :: ng0, ng1
  complex(kind(0d0)), intent(in)  :: wfn  (0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(out) :: d1wfn(0:ngrid, 1:nfun)

  call general_diff(wfn, d1wfn, fd_ohalf, fd_coeff1_in, fd_coeff1_l, fd_coeff1_r, ng0, ng1)

!old  fac = -half / ( dgrid * dgrid )
!old
!old  do ispin = 1, nspin
!old     do ifun = nfroz(ispin) + 1, nfun
!old        do igrid = ng0, ng1
!old           hwfn(igrid, ifun, ispin) = hwfn(igrid, ifun, ispin) &
!old    &    +      (wfn(igrid-1, ifun, ispin) &
!old    &    - two * wfn(igrid,   ifun, ispin) &
!old    &    +       wfn(igrid+1, ifun, ispin)) * fac
!old        end do
!old     end do
!old  end do

end subroutine general_d1wfn
!######################################################################
