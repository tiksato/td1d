!################################################################################
subroutine util_matinv(n, a, inva)

  implicit none
  integer, intent(in) :: n
  complex(kind(0d0)), intent(in) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: inva(1:n, 1:n)

  integer :: info, lwork
  complex(kind(0d0)) :: clwork
  integer, allocatable :: ipiv(:)
  complex(kind(0d0)), allocatable :: ca(:,:)
  complex(kind(0d0)), allocatable :: work(:)

  allocate(ipiv(1:n))
  allocate(ca(1:n,1:n))

  ca(1:n, 1:n) = a(1:n, 1:n)

  call zgetrf(n, n, ca, n, ipiv, info)

  call zgetri(n, ca, n, ipiv, clwork, -1, info)
  lwork = int(clwork)
  allocate(work(lwork))
  call zgetri(n, ca, n, ipiv, work, lwork, info)

  inva(1:n, 1:n) = ca(1:n, 1:n)

  deallocate(work)
  deallocate(ca)
  deallocate(ipiv)

end subroutine util_matinv
!################################################################################
!################################################################################
subroutine util_matinv_real(n, a, inva)

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: a(1:n, 1:n)
  real(kind(0d0)), intent(out) :: inva(1:n, 1:n)

  integer :: info, lwork
  real(kind(0d0)) :: rlwork
  integer, allocatable :: ipiv(:)
  real(kind(0d0)), allocatable :: ca(:,:)
  real(kind(0d0)), allocatable :: work(:)

  allocate(ipiv(1:n))
  allocate(ca(1:n,1:n))

  ca(1:n, 1:n) = a(1:n, 1:n)

  call dgetrf(n, n, ca, n, ipiv, info)

  call dgetri(n, ca, n, ipiv, rlwork, -1, info)
  lwork = int(rlwork)
  allocate(work(lwork))
  call dgetri(n, ca, n, ipiv, work, lwork, info)

  inva(1:n, 1:n) = ca(1:n, 1:n)

  deallocate(work)
  deallocate(ca)
  deallocate(ipiv)

end subroutine util_matinv_real
!################################################################################
!################################################################################
subroutine util_matinv2_real(m, n, dim, a)

  implicit none
  integer, intent(in) :: m, n, dim
  real(kind(0d0)), intent(inout) :: a(1:dim, 1:dim)

  integer :: i, j, k
  real(kind(0d0)), allocatable :: aa(:,:)
  real(kind(0d0)), allocatable :: atmp(:,:)

  allocate(aa(1:n,1:n))
  allocate(atmp(1:n,1:n))

  aa(1:n, 1:n) = 0.d+0
  do i = 1, n
     do j = 1, n
        do k = 1, m
           aa(i, j) = aa(i, j) + a(k, i) * a(k, j)
        end do
     end do
  end do

  atmp(1:n, 1:n) = 0.d+0
  call util_matinv_real(n, aa, atmp)

  a(1:dim, 1:dim) = 0.d+0
  do i = 1, n
     do j = 1, n
        a(j, i) = atmp(j, i)
     end do
  end do

  deallocate(atmp)
  deallocate(aa)

end subroutine util_matinv2_real
!################################################################################
