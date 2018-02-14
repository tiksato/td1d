!################################################################################
subroutine prop_new(time, dt, len, wfn0, wfn, hwfn, imethod, prop_type)

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: len, imethod
  character(len=*), intent(in) :: prop_type
  complex(kind(0d0)), intent(in) :: time
  complex(kind(0d0)), intent(inout) :: dt       ! overwritten in 'VRK5'
  complex(kind(0d0)), intent(in) :: wfn0(1:*)
  complex(kind(0d0)), intent(inout) :: wfn(1:*)
  complex(kind(0d0)), intent(inout) :: hwfn(1:*)
  !--------------------------------------------------------------------

  if (trim(prop_type) == 'RK4') then
     ! fixed step-size fourth order runge-kutta (rk4)
     call prop_rk4_new(time, dt, len, wfn0, wfn, hwfn, imethod)
  else if (trim(prop_type) == 'SPLIT_RK4') then
     ! fixed step-size fourth order runge-kutta (rk4)
     call prop_split_rk4(time, dt, len, wfn0, wfn, hwfn, imethod)
!nyi  else if (trim(prop_type) == 'ARK4') then
!nyi     ! fixed step-size fourth order runge-kutta (rk4)
!nyi     call prop_rk4(time, dt, len, wfn0, wfn, hwfn, imethod)
  else if (trim(prop_type) == 'VRK5') then
     ! variable step-size fifth order runge-kutta (vrk5)
     call prop_vrk5_new(time, dt, len, wfn0, wfn, hwfn, imethod)
!vyi  else if (trim(prop_type) == 'SPLIT_EXP') then
!vyi     ! split t and v, spectral for t, and exponential (2nd order magnus) for v
!vyi     if (imethod > 0) stop 'prop_split_exp, only for x2e.'
!vyi     call prop_split_exp(time, dt, len, wfn, hwfn, imethod)
  else if (trim(prop_type) == 'EULER') then
     ! First-order Euler method
     call prop_euler_new(len, wfn, hwfn)
  else
     write(6, "('bad prop_type in prop_new.')")
     stop
  end if

end subroutine prop_new
!################################################################################
subroutine prop(time, dt, len, wfn0, wfn, hwfn, imethod, prop_type)

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: len, imethod
  character(len=*), intent(in) :: prop_type
  complex(kind(0d0)), intent(in) :: time
  complex(kind(0d0)), intent(inout) :: dt       ! overwritten in 'VRK5'
  complex(kind(0d0)), intent(in) :: wfn0(1:*)
  complex(kind(0d0)), intent(inout) :: wfn(1:*)
  complex(kind(0d0)), intent(inout) :: hwfn(1:*)
  !--------------------------------------------------------------------

  if (trim(prop_type) == 'RK4') then
     ! fixed step-size fourth order runge-kutta (rk4)
     call prop_rk4(time, dt, len, wfn0, wfn, hwfn, imethod)
  else if (trim(prop_type) == 'ARK4') then
     ! fixed step-size fourth order runge-kutta (rk4)
     call prop_rk4(time, dt, len, wfn0, wfn, hwfn, imethod)
  else if (trim(prop_type) == 'VRK5') then
     ! variable step-size fifth order runge-kutta (vrk5)
     call prop_vrk5(time, dt, len, wfn0, wfn, hwfn, imethod)
  else if (trim(prop_type) == 'SPLIT_RK4') then
     ! split t and v, spectral for t, and rk4 for v
     call prop_split_rk4(time, dt, len, wfn0, wfn, hwfn, imethod)
  else if (trim(prop_type) == 'SPLIT_EXP') then
     ! split t and v, spectral for t, and exponential (2nd order magnus) for v
     if (imethod > 0) stop 'prop_split_exp, only for x2e.'
     call prop_split_exp(time, dt, len, wfn, hwfn, imethod)
  else if (trim(prop_type) == 'EULER') then
     ! First-order Euler method
     call prop_euler(dt, len, wfn, hwfn)
  else
     ! fixed step-size fourth order runge-kutta (rk4)
     call prop_rk4(time, dt, len, wfn0, wfn, hwfn, imethod)
  end if

end subroutine prop
!################################################################################
