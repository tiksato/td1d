!################################################################################
subroutine util_gmatinv(n, thresh, a, inva)

  use const_mod, only : czero, one

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: thresh
  complex(kind(0d0)), intent(in) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: inva(1:n, 1:n)

  integer :: i, j, k
!  complex(kind(0d0)), allocatable :: aa(:,:), raa(:,:)
!  complex(kind(0d0)), parameter :: czero = (0.d+0, 0.d+0)
  real(kind(0d0)) :: reig
  real(kind(0d0)), allocatable :: aeig(:)
  complex(kind(0d0)), allocatable :: uvec(:,:)
  complex(kind(0d0)), allocatable :: vvec(:,:)

  allocate(aeig(1:n))
  allocate(uvec(1:n, 1:n))
  allocate(vvec(1:n, 1:n))

  call util_svd_comp(n, n, a, aeig, uvec, vvec)
  inva(1:n, 1:n) = czero
  do k = 1, n
     reig = aeig(k) / (aeig(k) * aeig(k) + thresh * thresh)
!    reig = one / (aeig(k) + thresh * exp(-aeig(k)/thresh))
     do i = 1, n
        do j = 1, n
           inva(j, i) = inva(j, i) + vvec(j, k) * reig * conjg(uvec(i, k))
        end do
     end do
  end do

  deallocate(vvec)
  deallocate(uvec)
  deallocate(aeig)

!  allocate(aa(1:n,1:n))
!  allocate(raa(1:n,1:n))
!
!  aa(1:n, 1:n) = czero
!  do i = 1, n
!     do j = 1, n
!        do k = 1, n
!           aa(j, i) = aa(j, i) + a(j, k) * conjg(a(i, k))
!        end do
!     end do
!  end do
!  call util_matinv_reg(n, thresh, aa, raa)
!
!  inva(1:n, 1:n) = czero
!  do i = 1, n
!     do j = 1, n
!        do k = 1, n
!           inva(j, i) = inva(j, i) + conjg(a(k, j)) * raa(k, i)
!        end do
!     end do
!  end do
!
!  deallocate(raa)
!  deallocate(aa)

end subroutine util_gmatinv
!################################################################################
