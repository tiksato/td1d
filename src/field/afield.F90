!################################################################################
real(kind(0d0)) function afield(time)

  use const_mod, only : zero
  use field_mod, only : env_type

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: afield_trape, afield_sin2, afield_flat

  if (trim(env_type) == 'ZERO') then
     afield = zero
  else if (trim(env_type) == 'TRAPEZOIDAL') then
     afield = afield_trape(time)
  else if (trim(env_type) == 'SIN2') then
     afield = afield_sin2(time)
  else if (trim(env_type) == 'FLAT') then
     afield = afield_flat(time)
  else
     stop 'bad env_type.'
  end if

  return

end function afield
!################################################################################
