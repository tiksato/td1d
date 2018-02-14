!################################################################################
subroutine wfn_adjdt(dt, dt0, wfn, hwfn, imethod)

  use wfn_mod, only : ormas

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(inout) :: dt
  complex(kind(0d0)), intent(in) :: dt0
  complex(kind(0d0)), intent(in) :: wfn(*)
  complex(kind(0d0)), intent(in) :: hwfn(*)
  integer, intent(in) :: imethod
  !--------------------------------------------------------------------

  if (imethod == -1) then
  else if (imethod == 0) then
  end if

end subroutine wfn_adjdt
!################################################################################
