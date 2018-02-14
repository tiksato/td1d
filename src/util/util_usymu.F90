!################################################################################
subroutine util_usymu(n, u, f, g)
  !
  !  g = u^T * f * u, f and g are symmetric, u is unitary.
  !  
  use const_mod, only : zero, czero
  
  implicit none
  integer, intent(in) :: n
  complex(kind(0d0)), intent(in) :: u(1:n, 1:n)
  complex(kind(0d0)), intent(in) :: f(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: g(1:n, 1:n)
  integer :: i, j, k
  complex(kind(0d0)) :: tmp
  complex(kind(0d0)), allocatable :: work(:,:)

  allocate(work(1:n, 1:n))
  
  do i = 1, n
     do j = 1, n
        tmp = czero
        do k = 1, n
           tmp = tmp + f(i, k) * u(k, j)
        end do
        work(i, j) = tmp
     end do
  end do

  do i = 1, n
     do j = 1, i
        tmp = czero
        do k = 1, n
           tmp = tmp + u(k, i) * work(k, j)
        end do
        g(i, j) = tmp
     end do
  end do

  do i = 1, n
     do j = i + 1, n
        g(j, i) = g(i, j)
     end do
  end do
  
  deallocate(work)

end subroutine util_usymu
!################################################################################
