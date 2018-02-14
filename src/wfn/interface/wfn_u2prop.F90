!################################################################################
subroutine wfn_u2prop(time, dt, wfn, hwfn, imethod)

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: time, dt
  complex(kind(0d0)), intent(inout) :: wfn(1:*)
  complex(kind(0d0)), intent(inout) :: hwfn(1:*)
  !--------------------------------------------------------------------
end subroutine wfn_u2prop
!################################################################################
