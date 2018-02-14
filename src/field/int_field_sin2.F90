!################################################################################
real(kind(0d0)) function int_field_sin2(time, nint)
!
! analytical integral expression for sin2 envelope field
!
  use const_mod, only : zero, two, half, quart, pi
  use field_mod, only : tau, cep, famp, freq

  implicit none
  !--------------------------------------------------------------------
  real(kind(0d0)), intent(in) :: time
  integer, intent(in) :: nint
  !--------------------------------------------------------------------
  real(kind(0d0)) :: wplus, wmins, wt_cep, wpt_cep, wmt_cep

  wplus = freq + two / tau
  wmins = freq - two / tau
  wt_cep = freq * time + cep
  wpt_cep = wplus * time + cep
  wmt_cep = wmins * time + cep

  if (nint == 0) then
     int_field_sin2 = half * sin(wt_cep) &
    &  - quart * (sin(wpt_cep) &
    &           + sin(wmt_cep))
  else if (nint == 1) then
     int_field_sin2 = - half * cos(wt_cep) / freq &
    &  + quart * (cos(wpt_cep) / wplus &
    &           + cos(wmt_cep) / wmins)
  else if (nint == 2) then
     int_field_sin2 = - half * sin(wt_cep) / (freq * freq) &
    &  + quart * (sin(wpt_cep) / (wplus * wplus) &
    &           + sin(wmt_cep) / (wmins * wmins))
  else
     stop 'bad nint'
  end if

  int_field_sin2 = int_field_sin2 * famp

  return

end function int_field_sin2
!################################################################################
