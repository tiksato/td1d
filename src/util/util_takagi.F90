!################################################################################
subroutine util_takagi_ceig(dim, amat, aeig, uvec)
!
  use root_mod, only : iprint
!
! complex skew-symmetric matrix is factorized to the 
! direct sum of 2x2 real skew-symmetric factors.
! input  n: dimension
!        a: skew-symmetric matrix
! output a: diaga
!           diaga(2i-1,2i) =  sqrt(diag(aa+))
!           diaga(2i,2i-1) = -sqrt(diag(aa+))
!        u: unitary transformation matrix
!           diaga = u^T * a * u
!           a = v * diaga * v^T, v = conjg(u)
!
  implicit none
  integer, intent(in) :: dim
  complex(kind(0d0)), intent(in) :: amat(1:dim, 1:dim)
  complex(kind(0d0)), intent(out) :: aeig(1:*)
  complex(kind(0d0)), intent(out) :: uvec(1:dim, 1:dim)

  integer :: i, j, k, l, ii, dim2
  complex(kind(0d0)) :: tmp
  complex(kind(0d0)), allocatable :: atmp(:,:)
  complex(kind(0d0)), allocatable :: vvec(:,:)

  allocate(atmp(1:dim, 1:dim))
  allocate(vvec(1:dim, 1:dim))

  vvec(1:dim, 1:dim) = (0.d+0, 0.d+0)
  atmp(1:dim, 1:dim) = (0.d+0, 0.d+0)
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           atmp(i, j) = atmp(i, j) + amat(i, k) * conjg(amat(j, k))
        end do
     end do
  end do
  call util_diag_comp(.true., dim, atmp, vvec)

  do i = 1, dim
     do j = 1, dim
        uvec(i, j) = conjg(vvec(i, j))
     end do
  end do

  ! adiag = U^T a U
  atmp(1:dim, 1:dim) = (0.d+0, 0.d+0)
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           do l = 1, dim
              atmp(i, j) = atmp(i, j) + uvec(k, i) * amat(k, l) * uvec(l, j)
           end do
        end do
     end do
  end do

  !debug
  if (iprint > 4) then
     write(6, "('Takagi decomposition of a. real-part:')")
     do i = 1, dim
        do j = 1, dim
           write(6, "(f12.5)", advance = 'no') dble(atmp(i, j))
        end do
        write(6, *)
     end do
     write(6, "('Takagi decomposition of a. imag-part:')")
     do i = 1, dim
        do j = 1, dim
           write(6, "(f12.5)", advance = 'no') aimag(atmp(i, j))
        end do
        write(6, *)
     end do
  end if
  !debug

  dim2 = dim / 2
  do ii = 1, dim2
     i = 2 * ii - 1
     aeig(ii) = atmp(i, i+1)
  end do

  write(6, "('util_takagi_ceig. eigenvalues:')")
  do ii = 1, dim2
     write(6, "(i10, 2f16.8, e16.4)") ii, aeig(ii), abs(aeig(ii))
  end do
!debug  write(6, "('util_takagi. eigenvectors:')")
!debug  do ii = 1, dim2
!debug     i = 2 * ii - 1
!debug     write(6, "(i16)", advance = 'no') i
!debug  end do
!debug  write(6, *)
!debug  do j = 1, dim
!debug     do ii = 1, dim2
!debug        i = 2 * ii - 1
!debug        write(6, "(e16.4)", advance = 'no') dble(uvec(j, i))
!debug     end do
!debug     write(6, *)
!debug  end do
!debug  do ii = 1, dim2
!debug     i = 2 * ii
!debug     write(6, "(i16)", advance = 'no') i
!debug  end do
!debug  write(6, *)
!debug  do j = 1, dim
!debug     do ii = 1, dim2
!debug        i = 2 * ii
!debug        write(6, "(e16.4)", advance = 'no') dble(uvec(j, i))
!debug     end do
!debug     write(6, *)
!debug  end do

  deallocate(vvec)
  deallocate(atmp)

!debug
!  stop "for debug in util_takagi"
!debug

end subroutine util_takagi_ceig
!################################################################################
!################################################################################
subroutine util_takagi_reig(dim, amat, aeig, uvec)
!
! complex skew-symmetric matrix is factorized to the 
! direct sum of 2x2 real skew-symmetric factors.
! input  n: dimension
!        a: skew-symmetric matrix
! output a: diaga
!           diaga(2i-1,2i) =  sqrt(diag(aa+))
!           diaga(2i,2i-1) = -sqrt(diag(aa+))
!        u: unitary transformation matrix
!           diaga = u^T * a * u
!           a = v * diaga * v^T, v = conjg(u)
!
  implicit none
  integer, intent(in) :: dim
  complex(kind(0d0)), intent(in) :: amat(1:dim, 1:dim)
  complex(kind(0d0)), intent(out) :: aeig(1:*)
  complex(kind(0d0)), intent(out) :: uvec(1:dim, 1:dim)

  integer :: i, j, k, l, ii, dim2
  complex(kind(0d0)) :: tmp
  complex(kind(0d0)), allocatable :: atmp(:,:)
  complex(kind(0d0)), allocatable :: vvec(:,:)

  allocate(atmp(1:dim, 1:dim))
  allocate(vvec(1:dim, 1:dim))

  vvec(1:dim, 1:dim) = (0.d+0, 0.d+0)
  atmp(1:dim, 1:dim) = (0.d+0, 0.d+0)
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           atmp(i, j) = atmp(i, j) + amat(i, k) * conjg(amat(j, k))
        end do
     end do
  end do
  call util_diag_comp(.true., dim, atmp, vvec)

  do i = 1, dim
     do j = 1, dim
        uvec(i, j) = conjg(vvec(i, j))
     end do
  end do

  ! adiag = U^T a U
  atmp(1:dim, 1:dim) = (0.d+0, 0.d+0)
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           do l = 1, dim
              atmp(i, j) = atmp(i, j) + uvec(k, i) * amat(k, l) * uvec(l, j)
           end do
        end do
     end do
  end do

  dim2 = dim / 2
  do ii = 1, dim2
     i = 2 * ii - 1
     tmp = abs(atmp(i, i+1))
     tmp = atmp(i, i+1) / tmp
     tmp = (1.d+0, 0.d+0) / sqrt(tmp)
     do j = 1, dim
        uvec(j, i)   = uvec(j, i)   * tmp
        uvec(j, i+1) = uvec(j, i+1) * tmp
     end do
  end do

!debug
!  vvec(1:dim, 1:dim) = (0.d+0, 0.d+0)
!  do i = 1, dim
!     do j = 1, dim
!        do k = 1, dim
!           do l = 1, dim
!              vvec(i, j) = vvec(i, j) + uvec(k, i) * amat(k, l) * uvec(l, j)
!           end do
!        end do
!     end do
!  end do
!  write(6, "('Takagi decomposition of a. real-part:')")
!  do i = 1, dim
!     do j = 1, dim
!        write(6, "(f12.5)", advance = 'no') dble(vvec(i, j))
!     end do
!     write(6, *)
!  end do
!  write(6, "('Takagi decomposition of a. imag-part:')")
!  do i = 1, dim
!     do j = 1, dim
!        write(6, "(f12.5)", advance = 'no') aimag(vvec(i, j))
!     end do
!     write(6, *)
!  end do
!debug

  do ii = 1, dim2
     i = 2 * ii - 1
     aeig(ii) = abs(atmp(i, i+1))
  end do

  write(6, "('util_takagi_reig. eigenvalues:')")
  do ii = 1, dim2
     write(6, "(i10, 2f16.8, e16.4)") ii, aeig(ii), abs(aeig(ii))
  end do

  deallocate(vvec)
  deallocate(atmp)

end subroutine util_takagi_reig
!################################################################################
!################################################################################
subroutine util_takagi_real(dim, amat, aeig, uvec, vvec)
!
! real symmetric matrix is factorized into the diagonal form
! input   dim: dimension
!        amat: symmetric matrix
! output aeig: diaga(i) = sqrt(diag(aa+))
!        uvec: unitary transformation matrix
!        vvec: unitary transformation matrix
!              diaga = u^T * a * v
!              a = u * diaga * v^T
!
  implicit none
  integer, intent(in) :: dim
  real(kind(0d0)), intent(in) :: amat(1:dim, 1:dim)
  real(kind(0d0)), intent(out) :: aeig(1:*)
  real(kind(0d0)), intent(out) :: uvec(1:dim, 1:dim)
  real(kind(0d0)), intent(out) :: vvec(1:dim, 1:dim)

  integer :: i, j, k, l
  real(kind(0d0)), allocatable :: atmp(:,:)

  allocate(atmp(1:dim, 1:dim))

  vvec(1:dim, 1:dim) = 0.d+0
  atmp(1:dim, 1:dim) = 0.d+0
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           atmp(i, j) = atmp(i, j) + amat(i, k) * amat(j, k)
        end do
     end do
  end do
  call util_diag_real(.true., dim, atmp, vvec)

  do i = 1, dim
     do j = 1, dim
        uvec(i, j) = vvec(i, j)
     end do
  end do

  write(6, "('util_takagi_real: A')")
  do i = 1, dim
     do j = 1, dim
        write(6, "(f20.10)", advance = 'no') amat(i, j)
     end do
     write(6, *)
  end do

  write(6, "('util_takagi_real. eigenvalues of aa^t:')")
  do i = 1, dim
     write(6, "(i10, f20.10)") i, atmp(i, i)
  end do

  ! adiag = U^T a U
  atmp(1:dim, 1:dim) = 0.d+0
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           do l = 1, dim
              atmp(i, j) = atmp(i, j) + uvec(k, i) * amat(k, l) * uvec(l, j)
           end do
        end do
     end do
  end do

  write(6, "('util_takagi_real: a = Vt * A * V')")
  do i = 1, dim
     do j = 1, dim
        write(6, "(f20.10)", advance = 'no') atmp(i, j)
     end do
     write(6, *)
  end do

  do i = 1, dim
     if (atmp(i, i) < 0.d+0) then
        do j = 1, dim
           uvec(j, i)   = -uvec(j, i)
        end do
     end if
  end do

  ! adiag = U^T a V
  atmp(1:dim, 1:dim) = 0.d+0
  do i = 1, dim
     do j = 1, dim
        do k = 1, dim
           do l = 1, dim
              atmp(i, j) = atmp(i, j) + uvec(k, i) * amat(k, l) * vvec(l, j)
           end do
        end do
     end do
  end do

  write(6, "('takagi_real: a = Ut * A * V <=> A = U * a * Vt')")
  do i = 1, dim
     do j = 1, dim
        write(6, "(f20.10)", advance = 'no') atmp(i, j)
     end do
     write(6, *)
  end do

  do i = 1, dim
     aeig(i) = atmp(i, i)
  end do

  write(6, "('util_takagi_real. eigenvalues:')")
  do i = 1, dim
     write(6, "(i10, f20.10)") i, aeig(i)
  end do

  deallocate(atmp)

end subroutine util_takagi_real
!################################################################################
