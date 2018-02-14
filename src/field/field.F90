!################################################################################
real(kind(0d0)) function field(time)

  use field_mod, only : gauge

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: efield, afield

  if (aimag(time) > 1D-10) then
     field = 0d0
     write(6,"('field: imaginary time -> field is zero...')")
  else if (trim(gauge) == 'L') then
     field = efield(time)
  else if (trim(gauge) == 'V') then
     field = afield(time)
  else
     stop 'bad gauge.'
  end if

  return

end function field
!################################################################################
