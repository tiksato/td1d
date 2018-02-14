!################################################################################
real(kind(0d0)) function efield_trape(time)
!
! trapezoidal envelope
!
  use const_mod, only : one, zero
  use field_mod, only : cep, famp, freq, period, ncyc_on, ncyc_off, ncyc_flat

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  !--------------------------------------------------------------------
  real(kind(0d0)) :: env, wt, t, t0, t1, t2, t3

  wt = freq * real(time)
  t0 = zero
  t1 = t0 + dble(ncyc_on) * period
  t2 = t1 + dble(ncyc_flat) * period
  t3 = t2 + dble(ncyc_off) * period
  t = dble(time)

  env = zero
  if (t < t0) then
     env = zero
  else if (t < t1) then
     env = t / t1
  else if (t < t2) then
     env = one
  else if (t < t3) then
     env = (t3 - t) / (t3 - t2)
  else if (t > t3) then
     env = zero
  end if

  efield_trape = famp * env * sin(wt + cep)
  return

end function efield_trape
!################################################################################
