!################################################################################
subroutine util_hermc(n, a)

  use const_mod, only : chalf

  implicit none
  integer, intent(in) :: n
  complex(kind(0d0)), intent(inout) :: a(1:n, 1:n)

  integer :: i, j
  complex(kind(0d0)), allocatable :: ca(:,:)

  allocate(ca(1:n,1:n))

  ca(1:n, 1:n) = a(1:n, 1:n)
  do i = 1, n
     do j = 1, n
        a(i, j) = (a(i, j) + conjg(ca(j, i))) * chalf
     end do
  end do

  deallocate(ca)

end subroutine util_hermc
!################################################################################
!################################################################################
subroutine util_aherm(n, a)

  implicit none
  integer, intent(in) :: n
  complex(kind(0d0)), intent(inout) :: a(1:n, 1:n)

  integer :: i, j
  complex(kind(0d0)), allocatable :: ca(:,:)

  allocate(ca(1:n,1:n))

  ca(1:n, 1:n) = a(1:n, 1:n)
  do i = 1, n
     do j = 1, n
        a(i, j) = ca(i, j) - conjg(ca(j, i))
     end do
  end do

  deallocate(ca)

end subroutine util_aherm
!################################################################################
