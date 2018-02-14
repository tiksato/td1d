!################################################################################
subroutine util_matoutc(iout, n, mat)

  implicit none
  integer, intent(in) :: iout, n
  complex(kind(0d0)), intent(in) :: mat(1:n, 1:n)
  !--------------------------------------------------------------------
  integer :: i, j


!old  do i = 1, n
!old     do j = 1, n
!old        write(iout, "(e12.4)", advance='no') dble(mat(j, i))
!old     end do
!old     write(iout, "(' | ')", advance='no')
!old     do j = 1, n
!old        write(iout, "(e12.4)", advance='no') aimag(mat(j, i))
!old     end do
!old     write(iout, *)
!old  end do

  write(iout, "(' mat-real')")
  do i = 1, n
     do j = 1, n
        write(iout, "(e16.8)", advance='no') dble(mat(i, j))
     end do
     write(iout, *)
  end do


  write(iout, "(' mat-imag')")
  do i = 1, n
     do j = 1, n
        write(iout, "(e16.8)", advance='no') aimag(mat(i, j))
     end do
     write(iout, *)
  end do

end subroutine util_matoutc
!################################################################################
!################################################################################
subroutine util_matoutc2(iout, n, m, mat)

  implicit none
  integer, intent(in) :: iout, n, m
  complex(kind(0d0)), intent(in) :: mat(1:n, 1:m)
  !--------------------------------------------------------------------
  integer :: i, j


  write(iout, "(' mat-real')")
  do i = 1, m
     do j = 1, n
        write(iout, "(e16.8)", advance='no') dble(mat(i, j))
     end do
     write(iout, *)
  end do


  write(iout, "(' mat-imag')")
  do i = 1, m
     do j = 1, n
        write(iout, "(e16.8)", advance='no') aimag(mat(i, j))
     end do
     write(iout, *)
  end do

end subroutine util_matoutc2
!################################################################################
!################################################################################
subroutine util_linoutc(iout, n, mat)

  implicit none
  integer, intent(in) :: iout, n
  complex(kind(0d0)), intent(in) :: mat(1:*)
  !--------------------------------------------------------------------
  integer :: i, j, ji


  do i = 1, n
     do j = 1, i
        ji = (i * (i - 1)) / 2 + j
        write(iout, "(e16.8)", advance='no') dble(mat(ji))
     end do
     do j = i + 1, n
        write(iout, "(12x)", advance='no')
     end do
     write(iout, "(' | ')", advance='no')
     do j = 1, i
        ji = (i * (i - 1)) / 2 + j
        write(iout, "(e16.8)", advance='no') aimag(mat(ji))
     end do
     write(iout, *)
  end do

end subroutine util_linoutc
!################################################################################
