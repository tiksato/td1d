!################################################################################
real(kind(0d0)) function afield_sin2(time)

! sin2 envelope

  use const_mod, only : zero, pi
  use field_mod, only : tau

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  !--------------------------------------------------------------------
  real(kind(0d0)) :: t0, t
  real(kind(0d0)), external :: int_field_sin2

  t0 = zero
  t = min(dble(time), tau * pi)

  afield_sin2 = int_field_sin2(t, 1) &
            & - int_field_sin2(t0, 1)

  return

end function afield_sin2
!################################################################################
