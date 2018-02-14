!################################################################################
integer(8) function util_maxloc(ifabs, n, a)

  implicit none
  logical, intent(in) :: ifabs
  integer, intent(in) :: n
  real(kind(0d0)), intent(in) :: a(1:n)

  integer :: locmax, i
  real(kind(0d0)) :: valmax

  locmax = 1
  valmax = a(1)

  if (ifabs) then
     do i = 2, n
        if (abs(a(i)) > valmax) then
           locmax = i
           valmax = abs(a(i))
        end if
     end do
  else
     do i = 2, n
        if (a(i) > valmax) then
           locmax = i
           valmax = a(i)
        end if
     end do
  end if

  util_maxloc = locmax

end function util_maxloc
!################################################################################
