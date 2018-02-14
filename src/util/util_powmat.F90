!################################################################################
subroutine util_powmat(n, pow, a)
!
! b = a ** pow; a and b are hermetican
!
  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: pow
  complex(kind(0d0)), intent(inout) :: a(1:n, 1:n)

  complex(kind(0d0)), allocatable :: b(:), u(:,:)
  integer :: i, j, k

  allocate(b(1:n))
  allocate(u(1:n,1:n))

  ! diagonalize a
  call util_diag_comp(.true., n, a, u)

  ! b**pow
  do i = 1, n
     b(i) = a(i, i) ** pow
  end do

  ! a**pow = u * b**pow * u+
  do i = 1, n
     do j = 1, i
        a(i, j) = (0.d+0, 0.d+0)
        do k = 1, n
           a(i, j) = a(i, j) + u(i, k) * b(k) * conjg(u(j, k))
        end do
     end do
  end do

  ! hermitization
  do i = 1, n
     do j = i + 1, n
        a(i, j) = conjg(a(j, i))
     end do
  end do

  deallocate(u)
  deallocate(b)

end subroutine util_powmat
!################################################################################
