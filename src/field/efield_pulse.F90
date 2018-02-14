!################################################################################
real(kind(0d0)) function efield_pulse(time)

  use const_mod, only : one, zero
  use field_mod, only : famp, pulse_t, pulse_dt

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  !--------------------------------------------------------------------
  real(kind(0d0)) :: t, env

  t = dble(time)
  if (abs(t - pulse_t) < pulse_dt) then
     env = one
  else
     env = zero
  end if

  efield_pulse = famp * env
  return

end function efield_pulse
!################################################################################
