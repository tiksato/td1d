!################################################################################
subroutine util_lineq_solve(n, dima, fac, a, b, c)
!
! solve a*c = b by direct inversion
!
  use const_mod, only : czero, runit

  implicit none
  integer, intent(in) :: n, dima
  complex(kind(0d0)), intent(in) :: fac
  complex(kind(0d0)), intent(in) :: a(1:dima, 1:*)
  complex(kind(0d0)), intent(in) :: b(1:*)
  complex(kind(0d0)), intent(out) :: c(1:*)

  integer :: i, j
  complex(kind(0d0)), allocatable :: aa(:,:), ra(:,:), bb(:)

  allocate(aa(1:n,1:n))
  allocate(ra(1:n,1:n))
  allocate(bb(1:n))

  do i = 1, n
     do j = 1, n
        aa(j, i) = fac * a(j, i)
     end do
     bb(i) = fac * b(i)
  end do

  call util_matinv(n, aa, ra)
  call zgemm('n', 'n', n, 1, n, runit, ra, n, bb, n, czero, c, n)

  deallocate(bb)
  deallocate(ra)
  deallocate(aa)

end subroutine util_lineq_solve
!################################################################################
!################################################################################
subroutine util_lineq_resid(len, nvec, coeff, yvec, bvec, resid)
!
! have residue
!
  use const_mod, only : czero, runit

  implicit none
  integer, intent(in) :: nvec, len
  complex(kind(0d0)), intent(in) :: coeff(1:*)
  complex(kind(0d0)), intent(in) :: yvec(1:len, 1:*)
  complex(kind(0d0)), intent(in) :: bvec(1:len)
  complex(kind(0d0)), intent(out) :: resid(1:len)
  integer :: ivec

  call util_zcopy(len, czero, 0, resid, 1)
  do ivec = 1, nvec
     call util_zaxpy(len, coeff(ivec), yvec(1,ivec), 1, resid, 1)
  end do
  call util_zaxpy(len, -runit, bvec, 1, resid, 1)

end subroutine util_lineq_resid
!################################################################################
!################################################################################
subroutine util_lineq_value(len, nvec, coeff, xvec, sol)
!
! have residue
!
  use const_mod, only : czero

  implicit none
  integer, intent(in) :: nvec, len
  complex(kind(0d0)), intent(in) :: coeff(1:*)
  complex(kind(0d0)), intent(in) :: xvec(1:len, 1:*)
  complex(kind(0d0)), intent(out) :: sol(1:len)
  integer :: ivec

  call util_zcopy(len, czero, 0, sol, 1)
  do ivec = 1, nvec
     call util_zaxpy(len, coeff(ivec), xvec(1,ivec), 1, sol, 1)
  end do

end subroutine util_lineq_value
!################################################################################
!################################################################################
subroutine util_lineq_new(len, nvec, norm, resid, xvec)
!
! resid is orthogonalized against xvec(1:nvec-1) and copied to xvec(nvec)
!

  use const_mod, only : runit
  use thresh_mod, only : throcc

  implicit none
  integer, intent(in) :: len, nvec
  integer, intent(inout) :: norm
  complex(kind(0d0)), intent(in) :: resid(1:len)
  complex(kind(0d0)), intent(inout) :: xvec(1:len, 1:*)

  integer :: ivec
  complex(kind(0d0)) :: fac
  complex(kind(0d0)), external :: util_zdotc

  call util_zcopy(len, resid, 1, xvec(1,nvec), 1)
  do ivec = 1, nvec - 1
     fac = util_zdotc(len, xvec(1,ivec), 1, resid, 1)
     call util_zaxpy(len, -fac, xvec(1,ivec), 1, xvec(1,nvec), 1)
  end do

  if (norm >= 0) then
     fac = util_zdotc(len, xvec(1,nvec), 1, xvec(1,nvec), 1)
     if (norm > 0 .or. abs(fac) > throcc) then
        !debug
        write(6, "('lineq_new: xvec is normalized. <x|x>  = ', e15.5)") dble(fac)
        !debug
        fac = sqrt(runit / fac)
        call util_zscal(len, fac, xvec(1,nvec), 1)
        norm = 1
     end if
  end if

end subroutine util_lineq_new
!################################################################################
!!################################################################################
!subroutine util_lineq_check(len, nvec, resid, xvec)
!!
!! resid is orthogonalized against xvec(1:nvec-1) and copied to xvec(nvec)
!!
!  implicit none
!  integer, intent(in) :: len, nvec
!  complex(kind(0d0)), intent(in) :: resid(1:len)
!  complex(kind(0d0)), intent(inout) :: xvec(1:len, 1:*)
!
!  integer :: ivec
!  complex(kind(0d0)) :: fac
!
!  call util_zcopy(len, resid, 1, xvec(1,nvec), 1)
!  do ivec = 1, nvec - 1
!     fac = util_zdotc(len, xvec(1,ivec), 1, resid, 1)
!     call util_zaxpy(len, -fac, xvec(1,ivec), 1, xvec(1,nvec), 1)
!  end do
!
!end subroutine util_lineq_check
!!################################################################################
