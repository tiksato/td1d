!################################################################################
subroutine general_op1e(norb1, norb2, opx, orb1, orb2, opmat)

  use omp_mod
  use const_mod, only : czero
  use grid_mod, only : ngrid, dgrid

  implicit none

  integer, intent(in) :: norb1, norb2
  complex(kind(0d0)), intent(in) :: orb1(0:ngrid, 1:*)
  complex(kind(0d0)), intent(in) :: orb2(0:ngrid, 1:*)
  real(kind(0d0)), intent(in) :: opx(0:ngrid)
  complex(kind(0d0)), intent(out) :: opmat(1:norb1, 1:norb2)

  integer :: ifun, jfun, igrid
  complex(kind(0d0)) :: tmp

  !$omp parallel default(shared) private(tmp)
  !$omp do
  do jfun = 1, norb2
     do ifun = 1, norb1
        tmp = czero
        do igrid = 0, ngrid
           tmp = tmp + conjg(orb1(igrid, ifun)) * opx(igrid) * orb2(igrid, jfun)
        end do
        opmat(ifun, jfun) = tmp * dgrid
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine general_op1e
!################################################################################
