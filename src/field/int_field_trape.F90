!################################################################################
real(kind(0d0)) function int_field_trape1(time, nint)

! int_field_trape during linear turning on

  use const_mod, only : zero, one, two
  use field_mod, only : cep, famp, freq, period, ncyc_on, ncyc_off, ncyc_flat

  implicit none
  !--------------------------------------------------------------------
  real(kind(0d0)), intent(in) :: time
  integer, intent(in) :: nint
  !--------------------------------------------------------------------
  real(kind(0d0)) :: wt_cep, sint, cost, t0, t1, d, wd, w2d, w3d

  wt_cep = freq * time + cep
  sint = sin(wt_cep)
  cost = cos(wt_cep)

  t0 = zero
  t1 = t0 + dble(ncyc_on) * period

  d = t1 - t0
  wd = freq * d
  w2d = freq * wd
  w3d = freq * w2d

  if (nint == 0) then
     int_field_trape1 = (time - t0) * sint / d
  else if (nint == 1) then
     int_field_trape1 = sint / w2d - (time - t0) * cost / wd
  else if (nint == 2) then
     int_field_trape1 = - two * cost / w3d - (time - t0) * sint / w2d
  else
     stop 'bad nint'
  end if

  int_field_trape1 = int_field_trape1 * famp

  return

end function int_field_trape1
!################################################################################
!################################################################################
real(kind(0d0)) function int_field_trape2(time, nint)
!
! int_field_trape during flat maximum
!
  use const_mod, only : zero, one
  use field_mod, only : cep, famp, freq

  implicit none
  !--------------------------------------------------------------------
  real(kind(0d0)), intent(in) :: time
  integer, intent(in) :: nint
  !--------------------------------------------------------------------
  real(kind(0d0)) :: wt_cep, sint, cost, w, w2

  wt_cep = freq * time + cep
  sint = sin(wt_cep)
  cost = cos(wt_cep)
  w = freq
  w2 = w * w

  if (nint == 0) then
     int_field_trape2 = sint
  else if (nint == 1) then
     int_field_trape2 = -cost / w
  else if (nint == 2) then
     int_field_trape2 = -sint / w2
  else
     stop 'bad nint'
  end if

  int_field_trape2 = int_field_trape2 * famp

  return

end function int_field_trape2
!################################################################################
!################################################################################
real(kind(0d0)) function int_field_trape3(time, nint)
!
! int_field_trape during linear turning off
!
  use const_mod, only : zero, one, two
  use field_mod, only : cep, famp, freq, period, ncyc_on, ncyc_off, ncyc_flat

  implicit none
  !--------------------------------------------------------------------
  real(kind(0d0)), intent(in) :: time
  integer, intent(in) :: nint
  !--------------------------------------------------------------------
  real(kind(0d0)) :: wt_cep, sint, cost, t0, t1, t2, t3, d, wd, w2d, w3d

  wt_cep = freq * time + cep
  sint = sin(wt_cep)
  cost = cos(wt_cep)

  t0 = zero
  t1 = t0 + dble(ncyc_on) * period
  t2 = t1 + dble(ncyc_flat) * period
  t3 = t2 + dble(ncyc_off) * period

  d = t3 - t2
  wd = freq * d
  w2d = freq * wd
  w3d = freq * w2d

  if (nint == 0) then
     int_field_trape3 = (t3 - time) * sint / d
  else if (nint == 1) then
     int_field_trape3 = - sint / w2d - (t3 - time) * cost / wd
  else if (nint == 2) then
     int_field_trape3 = two * cost / w3d - (t3 - time) * sint / w2d
  else
     stop 'bad nint'
  end if

  int_field_trape3 = int_field_trape3 * famp

  return

end function int_field_trape3
!################################################################################
