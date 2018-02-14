!################################################################################
subroutine wfn_u1prop(time, dt, wfn, hwfn, imethod)

  use wfn_mod, only : domo, doci, nelact
  use const_mod, only : czero, zero, one, two, three, six, half, iunit

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: time, dt
  complex(kind(0d0)), intent(inout) :: wfn(1:*)
  complex(kind(0d0)), intent(inout) :: hwfn(1:*)
  !--------------------------------------------------------------------
end subroutine wfn_u1prop
!################################################################################
