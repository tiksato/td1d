!################################################################################
subroutine util_asyma(n, a, f, g)
  !
  !  g = a * f * a; f, a, and g are hermitian
  !  
  use const_mod, only : zero, czero
  
  implicit none
  integer, intent(in) :: n
  complex(kind(0d0)), intent(in) :: a(1:n, 1:n)
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
           tmp = tmp + f(i, k) * a(k, j)
        end do
        work(i, j) = tmp
     end do
  end do

  do i = 1, n
     do j = 1, i
        tmp = czero
        do k = 1, n
           tmp = tmp + a(i, k) * work(k, j)
        end do
        g(i, j) = tmp
     end do
  end do

  do i = 1, n
     do j = i + 1, n
        g(j, i) = conjg(g(i, j))
     end do
  end do
  
  deallocate(work)

end subroutine util_asyma
!################################################################################
