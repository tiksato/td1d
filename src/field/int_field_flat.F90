!################################################################################
real(kind(0d0)) function int_field_flat(time, nint)

! analytical integral expression for flat envelope field

  use field_mod, only : cep, famp, freq

  implicit none
  !--------------------------------------------------------------------
  real(kind(0d0)), intent(in) :: time
  integer, intent(in) :: nint
  !--------------------------------------------------------------------
  real(kind(0d0)) :: wt_cep

  wt_cep = freq * time + cep

  if (nint == 0) then
     int_field_flat = sin(wt_cep)
  else if (nint == 1) then
     int_field_flat = - cos(wt_cep) / freq
  else if (nint == 2) then
     int_field_flat = - sin(wt_cep) / (freq * freq)
  else
     stop 'bad nint'
  end if

  int_field_flat = int_field_flat * famp

  return

end function int_field_flat
!################################################################################
