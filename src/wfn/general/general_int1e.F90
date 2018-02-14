!######################################################################
subroutine general_int1e(nfun1, nfun2, wfn, h1wfn, int1e, ng0, ng1)

  use const_mod, only : czero
  use grid_mod, only : ngrid, dgrid

  implicit none
  integer, intent(in) :: nfun1, nfun2
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:*)
  complex(kind(0d0)), intent(in) :: h1wfn(0:ngrid, 1:*)
  complex(kind(0d0)), intent(inout) :: int1e(1:nfun1, 1:nfun2)
  integer, intent(in) :: ng0, ng1

  integer :: ifun, jfun, igrid
  complex(kind(0d0)) :: tmp

  do ifun = 1, nfun2
     do jfun = 1, nfun1
        tmp = czero
        do igrid = ng0, ng1
           tmp = tmp + conjg(wfn(igrid, jfun)) * h1wfn(igrid, ifun)
        end do
        int1e(jfun, ifun) = int1e(jfun, ifun) + tmp * dgrid
     end do
  end do

end subroutine general_int1e
!######################################################################
