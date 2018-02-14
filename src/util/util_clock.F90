!################################################################################
real(kind(0d0)) function util_clock(init, t1)

  use const_mod, only : zero

  implicit none
  logical, intent(in) :: init
  integer, intent(inout) :: t1
  integer :: t2, t_rate, t_max, diff

  if (init) then
     call system_clock(t1)
     util_clock = zero
  else
     call system_clock(t2, t_rate, t_max)
     if (t2 < t1) then
        diff = t_max - t1 + t2
     else
        diff = t2 - t1
     end if
     util_clock = diff / dble(t_rate)
  end if

end function util_clock
!################################################################################
