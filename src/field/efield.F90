!################################################################################
real(kind(0d0)) function efield(time)

  use const_mod, only : zero
  use field_mod, only : env_type

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: efield_trape, efield_sin2, efield_flat, efield_pulse

  if (trim(env_type) == 'ZERO') then
     efield = zero
  else if (trim(env_type) == 'TRAPEZOIDAL') then
     efield = efield_trape(time)
  else if (trim(env_type) == 'SIN2') then
     efield = efield_sin2(time)
  else if (trim(env_type) == 'FLAT') then
     efield = efield_flat(time)
  else if (trim(env_type) == 'PULSE') then
     efield = efield_pulse(time)
  else
     stop 'bad env_type.'
  end if

  return

end function efield
!################################################################################
