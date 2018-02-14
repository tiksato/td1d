!######################################################################
subroutine util_itochar(i, n, i_char)

  implicit none
  integer, intent(in) :: i, n
  character(len=16), intent(out) :: i_char
  integer :: nzero, izero

  write(i_char, "(i16)") i

  nzero = aint(log10(dble(max(1,n)))) - aint(log10(dble(max(1,i))))
  do izero = 1, nzero
     i_char = '0'//trim(adjustl(i_char))
  end do

  i_char = adjustl(i_char)

end subroutine util_itochar
!######################################################################
