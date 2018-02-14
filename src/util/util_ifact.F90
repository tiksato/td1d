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
!################################################################################
integer function util_ifact2(n, k)

  implicit none
  integer, intent(in) :: n, k
  integer i, m

  m = 1
  do i = n, n - k + 1, -1
     m = m * i
  end do
  
  util_ifact2 = m

end function util_ifact2
!################################################################################
