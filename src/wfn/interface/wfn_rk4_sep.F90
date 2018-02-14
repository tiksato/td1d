!################################################################################
subroutine wfn_rk4_sep(time, dt, wfn, hwfn, moci)
!
! RK4 propagation for mo (moci = 1) or ci (moci = 2).
!
  use const_mod, only : zero, one, two, three, six, half, iunit

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: moci
  complex(kind(0d0)), intent(in) :: time, dt
  complex(kind(0d0)), intent(inout) :: wfn(1:*)
  complex(kind(0d0)), intent(inout) :: hwfn(1:*)
  !--------------------------------------------------------------------
  real(kind(0d0)) :: lfield, ene, dtime, dinit
  complex(kind(0d0)) :: ttime, oo1, oo2, oo3, oo6
  integer :: len, pos
  real(kind(0d0)), external :: field
  complex(kind(0d0)), allocatable :: wfnini(:)
  complex(kind(0d0)), allocatable :: wfntmp(:)
  integer, external :: wfn_size, wfn_poscic
end subroutine wfn_rk4_sep
!################################################################################
