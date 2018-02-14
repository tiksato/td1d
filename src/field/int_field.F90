!################################################################################
real(kind(0d0)) function int_field(time, nint)

  use const_mod, only : zero
  use field_mod, only : env_type

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  integer, intent(in) :: nint
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: int_field_trape, int_field_sin2, int_field_flat

  int_field = zero
  if (env_type == 'ZERO') then
     int_field = zero
  else if (env_type == 'TRAPEZOIDAL') then
!nyi     int_field = int_field_trape(time, nint)
  else if (env_type == 'SIN2') then
     int_field = int_field_sin2(time, nint)
  else if (env_type == 'FLAT') then
     int_field = int_field_flat(time, nint)
  end if

  return

end function int_field
!################################################################################
