!################################################################################
real(kind(0d0)) function efield_sin2(time)

! sin2 envelope

  use const_mod, only : one, zero, pi
  use field_mod, only : tau, cep, famp, freq

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  !--------------------------------------------------------------------
  real(kind(0d0)) :: env, t, wt

  t = dble(time)
  wt = freq * t

  if (t < tau * pi) then
     env = sin(t / tau)
     env = env * env
  else
     env = zero
  end if

  efield_sin2 = famp * env * sin(wt + cep)
  return

end function efield_sin2
!################################################################################
