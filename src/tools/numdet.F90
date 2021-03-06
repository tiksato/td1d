!################################################################################
program numdet

  integer :: ne, no

end program numdet
!################################################################################
!################################################################################
integer function util_bicoeff(n, k)

  implicit none
  integer, intent(in) :: n, k
  integer, external :: util_ifact2
  integer :: kx, num, denom

!old  nk = n - k
!old  num = util_ifact(n)
!old  denom1 = util_ifact(nk)
!old  denom2 = util_ifact(k)
!old  util_bicoeff = num / (denom1 * denom2)
  kx = min(k, n - k)
  num = util_ifact2(n, kx)
  denom = util_ifact2(kx, kx)
  util_bicoeff = num / denom

!debugwrite(6, "('util_bicoeff: n = ', i5)"), n
!debugwrite(6, "('util_bicoeff: k = ', i5)"), kx
!debugwrite(6, "('util_bicoeff: num = ', i20)"), num
!debugwrite(6, "('util_bicoeff: de1 = ', i20)"), denom

end function util_bicoeff
!################################################################################
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
