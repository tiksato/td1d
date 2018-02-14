!################################################################################
subroutine util_gdiag_comp_new(dsc, n, a, uleft, uright)

  implicit none
  logical, intent(in) :: dsc
  integer, intent(in) :: n
  complex(kind(0d0)), intent(inout) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: uleft(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: uright(1:n, 1:n)

  integer :: info, lwork, i, j, k, len
  complex(kind(0d0)) :: clwork, test

  integer, allocatable :: index(:)
  complex(kind(0d0)), allocatable :: tmpa(:, :)
  complex(kind(0d0)), allocatable :: eig(:), work(:), vecl(:,:), vecr(:,:), plr(:)
  real(kind(0d0)), allocatable :: rwork(:), deig(:)
!debug  complex(kind(0d0)) :: norm
!debug  complex(kind(0d0)) :: sij
!debug  complex(kind(0d0)), external :: util_zdotc, util_zdotu

  len = 2 * n
  allocate(index(n))
  allocate(rwork(len))
  allocate(tmpa(n, n))
  allocate(vecl(n, n))
  allocate(vecr(n, n))
  allocate(eig(n))
  allocate(deig(n))
  allocate(plr(n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

!debug
!debugwrite(6, "('A: A')")
!debugwrite(6, "(2f20.10)") tmpa(1:n, 1:n)
!debug
  call zgeev("v", "v", n, tmpa, n, eig, vecl, n, vecr, n, clwork, -1, rwork, info)
  lwork = int(clwork)

  allocate(work(lwork))
  call zgeev("v", "v", n, tmpa, n, eig, vecl, n, vecr, n, work, lwork, rwork, info)

  if (info /= 0) then
     write(6, "('util_gdiag_comp: info = ', i20)") info
     stop 'error in util_gdiag_comp.'
  else
     uleft (1:n, 1:n) = (0.d+0, 0.d+0)
     uright(1:n, 1:n) = (0.d+0, 0.d+0)
     a     (1:n, 1:n) = (0.d+0, 0.d+0)
!     if (dsc) then
!        do i = 1, n
!           a(i, i) = eig(n - i + 1)
!           uleft (1:n, i) = vecl(1:n, n - i + 1)
!           uright(1:n, i) = vecr(1:n, n - i + 1)
!        end do
!     else
!        do i = 1, n
!           a(i, i) = eig(i)
!           uleft (1:n, i) = vecl(1:n, i)
!           uright(1:n, i) = vecr(1:n, i)
!        end do
!     end if

     deig(1:n) = dble(eig(1:n))
     call util_index_double(dsc, n, deig, index)
     do i = 1, n
        j = index(i)
        a(i, i) = eig(j)
!debug        norm = util_zdotc(n, vecl(1,j), 1, vecr(1,j), 1)
!debug        uleft (1:n, i) = vecl(1:n, j)
!debug        uright(1:n, i) = vecr(1:n, j) / norm
        uleft (1:n, i) = vecl(1:n, j)
        uright(1:n, i) = vecr(1:n, j)
     end do
  end if

  ! renormalize eigenvector to hold biorthonormality
  ! write(6, "('gdiag_comp_new: biorthogonality check')")
  do i = 1, n
     do j = 1, n
        test = (0.d+0, 0.d+0)
        do k = 1, n
           test = test + conjg(uleft(k, i)) * uright(k, j)
        end do
        if (i == j) plr(i) = test
        ! write(6, "(2i5, 2f20.10)") i, j, test
     end do
  end do
  do i = 1, n
     do j = 1, n
        test = sqrt(plr(j))
        uright(i, j) = uright(i, j) / test
        uleft(i, j) = uleft(i, j) / conjg(test)
     end do
  end do

  deallocate(work)
  deallocate(plr)
  deallocate(deig)
  deallocate(eig)
  deallocate(vecr)
  deallocate(vecl)
  deallocate(tmpa)
  deallocate(rwork)
  deallocate(index)

end subroutine util_gdiag_comp_new
!################################################################################
subroutine util_gdiag_comp(dsc, n, a, uleft, uright)

  implicit none
  logical, intent(in) :: dsc
  integer, intent(in) :: n
  complex(kind(0d0)), intent(inout) :: a(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: uleft(1:n, 1:n)
  complex(kind(0d0)), intent(out) :: uright(1:n, 1:n)

  integer :: info, lwork, i, j, len
  complex(kind(0d0)) :: clwork

  integer, allocatable :: index(:)
  complex(kind(0d0)), allocatable :: tmpa(:, :)
  complex(kind(0d0)), allocatable :: eig(:), work(:), vecl(:,:), vecr(:,:)
  real(kind(0d0)), allocatable :: rwork(:), deig(:)
!debug  complex(kind(0d0)) :: norm
!debug  complex(kind(0d0)) :: sij
!debug  complex(kind(0d0)), external :: util_zdotc, util_zdotu

  len = 2 * n
  allocate(index(n))
  allocate(rwork(len))
  allocate(tmpa(n, n))
  allocate(vecl(n, n))
  allocate(vecr(n, n))
  allocate(eig(n))
  allocate(deig(n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

!debug
!debugwrite(6, "('A: A')")
!debugwrite(6, "(2f20.10)") tmpa(1:n, 1:n)
!debug
  call zgeev("v", "v", n, tmpa, n, eig, vecl, n, vecr, n, clwork, -1, rwork, info)
  lwork = int(clwork)

  allocate(work(lwork))
  call zgeev("v", "v", n, tmpa, n, eig, vecl, n, vecr, n, work, lwork, rwork, info)

!debug
!debugwrite(6, "('A: A = U * a * U+')")
!debugdo i = 1, n
!debug   do j = 1, n
!debug      a(i, j) = (0.d+0, 0.d+0)
!debug      do k = 1, n
!debug         a(i, j) = a(i, j) + tmpa(i, k) * eig(k) * conjg(tmpa(j, k))
!debug      end do
!debug   end do
!debugend do
!debugwrite(6, "(2f20.10)") a(1:n, 1:n)
!debug

  if (info /= 0) then
     write(6, "('util_gdiag_comp: info = ', i20)") info
     stop 'error in util_gdiag_comp.'
  else
     uleft (1:n, 1:n) = (0.d+0, 0.d+0)
     uright(1:n, 1:n) = (0.d+0, 0.d+0)
     a     (1:n, 1:n) = (0.d+0, 0.d+0)
!     if (dsc) then
!        do i = 1, n
!           a(i, i) = eig(n - i + 1)
!           uleft (1:n, i) = vecl(1:n, n - i + 1)
!           uright(1:n, i) = vecr(1:n, n - i + 1)
!        end do
!     else
!        do i = 1, n
!           a(i, i) = eig(i)
!           uleft (1:n, i) = vecl(1:n, i)
!           uright(1:n, i) = vecr(1:n, i)
!        end do
!     end if

     deig(1:n) = dble(eig(1:n))
     call util_index_double(dsc, n, deig, index)
     do i = 1, n
        j = index(i)
        a(i, i) = eig(j)
!debug        norm = util_zdotc(n, vecl(1,j), 1, vecr(1,j), 1)
!debug        uleft (1:n, i) = vecl(1:n, j)
!debug        uright(1:n, i) = vecr(1:n, j) / norm
        uleft (1:n, i) = vecl(1:n, j)
        uright(1:n, i) = vecr(1:n, j)
     end do
  end if

!debug
!write(6,"('gdiag: eigenvalues')")
!do i = 1, 3
!   write(6,"(4F20.10)") eig(i), a(i,i)
!end do
!write(6,"('gdiag: bi-orthonormality check:')")
!i = 1; j = 1
!write(6,"(i10,i10,4f20.10)") i, j, util_zdotc(n, uleft(1,i), 1, uright(1,j), 1), util_zdotu(n, uleft(1,i), 1, uright(1,j), 1)
!i = 1; j = 2
!write(6,"(i10,i10,4f20.10)") i, j, util_zdotc(n, uleft(1,i), 1, uright(1,j), 1), util_zdotu(n, uleft(1,i), 1, uright(1,j), 1)
!i = 2; j = 1
!write(6,"(i10,i10,4f20.10)") i, j, util_zdotc(n, uleft(1,i), 1, uright(1,j), 1), util_zdotu(n, uleft(1,i), 1, uright(1,j), 1)
!i = 2; j = 2
!write(6,"(i10,i10,4f20.10)") i, j, util_zdotc(n, uleft(1,i), 1, uright(1,j), 1), util_zdotu(n, uleft(1,i), 1, uright(1,j), 1)
!
!write(6,"('gdiag: orthonormality check:')")
!i = 1; j = 1
!write(6,"(i10,i10,4f20.10)") i, j, util_zdotc(n, uleft(1,i), 1, uleft(1,j), 1), util_zdotu(n, uright(1,i), 1, uright(1,j), 1)
!i = 1; j = 2
!write(6,"(i10,i10,4f20.10)") i, j, util_zdotc(n, uleft(1,i), 1, uleft(1,j), 1), util_zdotu(n, uright(1,i), 1, uright(1,j), 1)
!i = 2; j = 1
!write(6,"(i10,i10,4f20.10)") i, j, util_zdotc(n, uleft(1,i), 1, uleft(1,j), 1), util_zdotu(n, uright(1,i), 1, uright(1,j), 1)
!i = 2; j = 2
!write(6,"(i10,i10,4f20.10)") i, j, util_zdotc(n, uleft(1,i), 1, uleft(1,j), 1), util_zdotu(n, uright(1,i), 1, uright(1,j), 1)
!debug

  deallocate(work)
  deallocate(deig)
  deallocate(eig)
  deallocate(vecr)
  deallocate(vecl)
  deallocate(tmpa)
  deallocate(rwork)
  deallocate(index)

end subroutine util_gdiag_comp
!################################################################################
