!################################################################################
subroutine util_exp_asym(n, a, expa)

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(inout) :: a(1:n, 1:n)
  real(kind(0d0)), intent(out) :: expa(1:n, 1:n)

  integer :: i, j, k, l, ii, no2, num1
  real(kind(0d0)) :: tmp, norm
  real(kind(0d0)), parameter :: tiny = 1.D-10
  real(kind(0d0)), allocatable :: aa(:,:)
  real(kind(0d0)), allocatable :: uv(:,:)

!  if (mod(n, 2) /= 0) stop 'util_exp_asym: n is not even.'

  allocate(aa(1:n, 1:n))
  allocate(uv(1:n, 1:n))

  aa(1:n, 1:n) = 0.d+0
  do i = 1, n
     do j = 1, n
        do k = 1, n
           aa(j, i) = aa(j, i) + a(j, k) * a(i, k)
        end do
     end do
  end do

!debug  write(6, "('a:')")
!debug  do i = 1, n
!debug     do j = 1, n
!debug        write(6, "(f12.5)", advance = 'no') a(j, i)
!debug     end do
!debug     write(6, *)
!debug  end do
!debug  write(6, "('aa+:')")
!debug  do i = 1, n
!debug     do j = 1, n
!debug        write(6, "(f12.5)", advance = 'no') aa(j, i)
!debug     end do
!debug     write(6, *)
!debug  end do

! call util_diag_real(.false., n, aa, uv)
  call util_diag_real(.true., n, aa, uv)

!debug  write(6, "('eigenvalues of aa+:')")
!debug  do j = 1, n
!debug     write(6, "(f12.5)", advance = 'no') aa(j, j)
!debug  end do
!debug  write(6, *)

!  write(6, "('eigenvectors of aa+:')")
  num1 = 0
  do i = 1, n
     norm = 0.d+0
     do j = 1, n
        norm = norm + uv(j, i) * uv(j, i)
     end do
     norm = sqrt(norm)
!     write(6, "(f12.5)", advance = 'no') norm
     if (abs(aa(i,i)) > tiny) num1 = num1 + 1
  end do
!  write(6, *)

!  do j = 1, n
!     do i = 1, n
!        write(6, "(f12.5)", advance = 'no') uv(j, i)
!     end do
!     write(6, *)
!  end do

  if (mod(num1, 2) == 0) then
     aa(1:n, 1:n) = 0.d+0
     do i = 1, num1
        do j = 1, num1
           do k = 1, n
              do l = 1, n
                 aa(j, i) = aa(j, i) + uv(k, j) * a(k, l) * uv(l, i)
              end do
           end do
        end do
     end do
     
!debug     write(6, "('eigen-factorization of a:')")
!debug     do i = 1, n
!debug        do j = 1, n
!debug           write(6, "(f12.5)", advance = 'no') aa(j, i)
!debug        end do
!debug        write(6, *)
!debug     end do
     
     no2 = n / 2
     do ii = 1, no2
        i = 2 * ii - 1
        tmp = aa(i, i+1)
        aa(i, i)     =  cos(tmp)
        aa(i, i+1)   =  sin(tmp)
        aa(i+1, i)   = -sin(tmp)
        aa(i+1, i+1) =  cos(tmp)
     end do
     do i = 2 * no2 + 1, n
        aa(i, i) = 1.d+0
     end do
     
     expa(1:n, 1:n) = 0.d+0
     do i = 1, n
        do j = 1, n
           do ii = 1, no2
              k = 2 * ii - 1
              l = 2 * ii - 1
              expa(j, i) = expa(j, i) + uv(j, k) * aa(k, l) * uv(i, l)
     
              k = 2 * ii - 1
              l = 2 * ii
              expa(j, i) = expa(j, i) + uv(j, k) * aa(k, l) * uv(i, l)
     
              k = 2 * ii
              l = 2 * ii - 1
              expa(j, i) = expa(j, i) + uv(j, k) * aa(k, l) * uv(i, l)
     
              k = 2 * ii
              l = 2 * ii
              expa(j, i) = expa(j, i) + uv(j, k) * aa(k, l) * uv(i, l)
           end do
           do k = 2 * no2 + 1, n
              expa(j, i) = expa(j, i) + uv(j, k) * aa(k, k) * uv(i, k)
           end do
        end do
     end do
     
!debug     write(6, "('exp-a in the spectral form:')")
!debug     do i = 1, n
!debug        do j = 1, n
!debug           write(6, "(f12.5)", advance = 'no') aa(j, i)
!debug        end do
!debug        write(6, *)
!debug     end do
!debug     
!debug     write(6, "('exp-a:')")
!debug     do i = 1, n
!debug        do j = 1, n
!debug           write(6, "(f12.5)", advance = 'no') expa(j, i)
!debug        end do
!debug        write(6, *)
!debug     end do

!  else if (num1 == n) then
!
!     expa(1:n, 1:n) = 0.d+0
!     do i = 1, n
!        expa(i, i) = 1.d+0
!     end do
!
!     write(6, "('exp-a:')")
!     do i = 1, n
!        do j = 1, n
!           write(6, "(f12.5)", advance = 'no') expa(j, i)
!        end do
!        write(6, *)
!     end do
  else
     stop 'bad num1 in util_exp_asym'
  end if

  deallocate(uv)
  deallocate(aa)

end subroutine util_exp_asym
!################################################################################
