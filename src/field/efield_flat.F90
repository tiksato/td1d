!################################################################################
real(kind(0d0)) function efield_flat(time)

  ! flat envelope

  use field_mod, only : cep, famp, freq

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  !--------------------------------------------------------------------
  real(kind(0d0)) :: wt

  wt = freq * real(time)
  efield_flat = famp * sin(wt + cep)
  return

end function efield_flat
!################################################################################
