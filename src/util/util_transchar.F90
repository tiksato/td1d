!################################################################################
subroutine util_transchar(a)
!
! transform lower-case letters to upper-case ones
!
  implicit none
  character(len=*), intent(inout) :: a
!
  integer :: i
!
  do i = 1, len_trim(a)
     if ((a(i:i) >= 'a') .and. (a(i:i) <= 'z')) then
        a(i:i) = char(ichar(a(i:i)) - (ichar('a')-ichar('A')))
     end if
  end do
end subroutine util_transchar
!################################################################################
