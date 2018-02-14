!################################################################################
real(kind(0d0)) function afield_flat(time)

! vector potential of flat envelope

  use const_mod, only : zero

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(in) :: time
  !--------------------------------------------------------------------
  real(kind(0d0)) :: time0, time1
  real(kind(0d0)), external :: int_field_flat

  time0 = zero
  time1 = dble(time)
  afield_flat = int_field_flat(time1, 1) - int_field_flat(time0, 1)

  return

end function afield_flat
!################################################################################
