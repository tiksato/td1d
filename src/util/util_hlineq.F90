!################################################################################
subroutine util_hlineq(n, maxcyc, a, b, c)
!
! solve a*c = b by an iterative method
!
  use const_mod, only : czero, runit

  implicit none
  integer, intent(in) :: n, maxcyc
  complex(kind(0d0)), intent(in) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(in) :: b(1:n)
  complex(kind(0d0)), intent(out) :: c(1:n)

  ! initial guess
  do i = 1, n
     x(i, 1) = b(i) / a(i, i)
  end do

  do icyc = 1, maxcyc

     ! y = a * x
     call zgemm('n', 'n', n, 1, n, runit, a, n, x(1,icyc), n, czero, y(1,icyc), n)

     ! r.h.s in the current basis
     b0(icyc) = util_zdotc(n, x(1, icyc), 1, b, 1)

     ! coefficients in the current basis
     a0(icyc, icyc) = util_zdotc(n, x(1, icyc), 1, y(1, icyc), 1)
     do jcyc = 1, icyc - 1
        a0(jcyc, icyc) = util_zdotc(n, x(1, jcyc), 1, y(1, icyc), 1)
        a0(icyc, jcyc) = conjg(a0(jcyc, icyc))
     end do

     call util_hlineq_solve(icyc, n, a0, b0, c0)
     call util_hlineq_resid(icyc, y, b0, resid)
     call util_hlineq_new(icyc+1, resid, x)

  end do

end subroutine util_lineq
!################################################################################
