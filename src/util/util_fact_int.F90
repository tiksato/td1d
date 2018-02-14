!################################################################################
integer function util_ifact(n)

  implicit none
  integer, intent(in) :: n
  integer i, m

  m = 1
  do i = 2, n
     m = m * i
  end do
  
  util_ifact = m

end function util_ifact
!################################################################################
