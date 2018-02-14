!######################################################################
subroutine general_hwfn_weff(nfun1, nfun2, wfn1, wfn2, weff, ng0, ng1)

  use const_mod, only : czero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid, dgrid, v2, gll, gul

  implicit none
  integer, intent(in) :: nfun1, nfun2
  complex(kind(0d0)), intent(in) :: wfn1(0:ngrid, 1:*)
  complex(kind(0d0)), intent(in) :: wfn2(0:ngrid, 1:*)
  complex(kind(0d0)), intent(inout) :: weff(0:ngrid, 1:nfun1, 1:nfun2)
  integer, intent(in) :: ng0, ng1

  integer :: igrid, jgrid, ifun, jfun
  complex(kind(0d0)) :: rho

  do ifun = 1, nfun2
     do jfun = 1, nfun1
        do jgrid = gll, gul
           rho = conjg(wfn1(jgrid, jfun)) &
                   & * wfn2(jgrid, ifun) * dgrid
           if (abs(rho) > thrwfn) then
              do igrid = ng0, ng1
                 weff(igrid, jfun, ifun) = &
               & weff(igrid, jfun, ifun) + rho * v2(igrid, jgrid)
              end do
           end if
        end do
     end do
  end do

end subroutine general_hwfn_weff
!######################################################################
