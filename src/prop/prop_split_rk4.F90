!################################################################################
subroutine prop_split_rk4(time, dt, len, wfn0, wfn, hwfn, imethod)

  use mol_mod, only : ne
  use const_mod, only : half

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: len, imethod
  complex(kind(0d0)), intent(in) :: time, dt
  complex(kind(0d0)), intent(in)    :: wfn0(1:len)
  complex(kind(0d0)), intent(inout) :: wfn (1:len)
  complex(kind(0d0)), intent(out)   :: hwfn(1:len)
  !--------------------------------------------------------------------

  complex(kind(0d0)) :: dt2, ttime
  real(kind(0d0)), external :: wfn_v2prod

  dt2 = dt * half
  ! first step
  ttime = time
  !call wfn_exp_kinetic_fft(dt2, wfn, hwfn, imethod)
  call general_lanczos(ttime, dt2, wfn)

  ! second step
  if (ne(3) > 1) call prop_v2rk4(dt, len, wfn, hwfn, imethod)

  ! third step
  ttime = time + dt2
  !call wfn_exp_kinetic_fft(dt2, hwfn, wfn, imethod)
  call general_lanczos(ttime, dt2, wfn)

end subroutine prop_split_rk4
!################################################################################
