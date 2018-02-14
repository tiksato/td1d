!################################################################################
subroutine prop_split_exp(time, dt, len, wfn, hwfn, imethod)

  use const_mod, only : half

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: len, imethod
  complex(kind(0d0)), intent(in) :: time, dt
  complex(kind(0d0)), intent(inout) :: wfn (1:len)
  complex(kind(0d0)), intent(out)   :: hwfn(1:len)
  !--------------------------------------------------------------------
  real(kind(0d0)) :: lfield
  complex(kind(0d0)) :: ttime, dt1, dt2
  real(kind(0d0)), external :: field

  dt1 = dt
  dt2 = dt * half

  ! first step
  call wfn_exp_kinetic_fft(dt2, wfn, hwfn, imethod)

  ! second step
  ttime = time + dt2
  lfield = field(ttime)
  call wfn_exp_potential(dt1, lfield, hwfn, hwfn, imethod)

  ! third step
  call wfn_exp_kinetic_fft(dt2, hwfn, wfn, imethod)

end subroutine prop_split_exp
!################################################################################
