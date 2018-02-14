!################################################################################
real(kind(0d0)) function afield_trape(time)
!
! vector potential of trapezoidal envelope
!
  use const_mod, only : zero
  use field_mod, only : period, ncyc_on, ncyc_off, ncyc_flat

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  !--------------------------------------------------------------------
  real(kind(0d0)) :: t, t0, t1, t2, t3
  real(kind(0d0)), external :: int_field_trape1, int_field_trape2, int_field_trape3

  t = dble(time)
  t0 = zero
  t1 = t0 + dble(ncyc_on) * period
  t2 = t1 + dble(ncyc_flat) * period
  t3 = t2 + dble(ncyc_off) * period

  afield_trape = zero
  if (t < t0) then
     afield_trape = zero
  else if (t < t1) then
     afield_trape = int_field_trape1(t, 1) &
                & - int_field_trape1(t0, 1)
  else if (t < t2) then
     afield_trape = int_field_trape1(t1, 1) &
                & - int_field_trape1(t0, 1) &
                & + int_field_trape2(t, 1) &
                & - int_field_trape2(t1, 1)
  else if (t < t3) then
     afield_trape = int_field_trape1(t1, 1) &
                & - int_field_trape1(t0, 1) &
                & + int_field_trape2(t2, 1) &
                & - int_field_trape2(t1, 1) &
                & + int_field_trape3(t, 1) &
                & - int_field_trape3(t2, 1)
  else if (t > t3) then
     afield_trape = int_field_trape1(t1, 1) &
                & - int_field_trape1(t0, 1) &
                & + int_field_trape2(t2, 1) &
                & - int_field_trape2(t1, 1) &
                & + int_field_trape3(t3, 1) &
                & - int_field_trape3(t2, 1)
  end if

  return

end function afield_trape
!################################################################################
