!################################################################################
complex(kind(0d0)) function util_det(n, thresh, a)

  use const_mod, only : one, runit, czero

  implicit none
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: thresh
  complex(kind(0d0)), intent(in) :: a(1:n, 1:n)

  logical :: zero_triv
  integer :: npiv, i, j
  complex(kind(0d0)) :: test
  complex(kind(0d0)), external :: util_zdotu
  complex(kind(0d0)), allocatable :: tmpa(:,:), norm(:)

  ! check the trivial zero determinant
  ! **********************************
  zero_triv = .false.

  allocate(norm(1:n))
  do i = 1, n
     test = util_zdotu(n, a(1, i), 1, a(1, i), 1)
     norm(i) = sqrt(test)
     if (abs(test) < thresh) then
        zero_triv = .true.
        exit
     end if
  end do

  do i = 1, n
     if (zero_triv) exit
     do j = 1, i - 1
        test = util_zdotu(n, a(1, j), 1, a(1, i), 1) / (norm(j) * norm(i))
        if (abs(runit - test) < thresh .or. abs(-runit - test) < thresh) then
           zero_triv = .true.
           exit
        end if
     end do
  end do
  deallocate(norm)

  if (zero_triv) then
     util_det = czero
     return
  end if
  ! **********************************

  allocate(tmpa(1:n, 1:n))
  tmpa(1:n, 1:n) = a(1:n, 1:n)

!debug
!  write(6, "('before:')")
!  do i = 1, n
!     write(6, "(3f20.10)") real(tmpa(i, 1:n))
!  end do
!debug

  call util_trimat(n, thresh, npiv, tmpa)
  util_det = (-one) ** npiv
  do i = 1, n
     util_det = util_det * tmpa(i, i)
  end do

!debug
!  write(6, "('after:')")
!  do i = 1, n
!     write(6, "(3f20.10)") real(tmpa(i, 1:n))
!  end do
!debug

  deallocate(tmpa)

end function util_det
!################################################################################
