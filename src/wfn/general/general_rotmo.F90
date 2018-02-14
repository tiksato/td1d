!######################################################################
subroutine general_rotmo(norb, rmat, wfn, twfn, ng0, ng1)
!
! MO rotation
!
  use const_mod, only : czero
  use thresh_mod, only : thrwfn
  use grid_mod, only : ngrid

  implicit none
  integer, intent(in) :: norb, ng0, ng1
  complex(kind(0d0)), intent(in) :: rmat(1:norb, 1:norb)
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:*)
  complex(kind(0d0)), intent(inout) :: twfn(0:ngrid, 1:*)

  integer :: ifun, jfun, igrid
  complex(kind(0d0)) :: fac

  if (norb == 0) return

  do ifun = 1, norb
     do jfun = 1, norb
        fac = rmat(jfun, ifun)
        if (abs(fac) > thrwfn) then
           do igrid = ng0, ng1
              twfn(igrid, ifun) = twfn(igrid, ifun) + wfn(igrid, jfun) * fac
           end do
!blas           call util_zaxpy(nbno, fac, wfn(ng0,jfun), 1, twfn(ng0,ifun), 1)
        end if
     end do
  end do

end subroutine general_rotmo
!######################################################################
