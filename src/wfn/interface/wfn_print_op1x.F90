!################################################################################
subroutine wfn_print_op1x(iout, istep, time, wfn, imethod)

  use wfn_mod, only : ormas

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: iout, istep, imethod
  complex(kind(0d0)), intent(in) :: time
  complex(kind(0d0)), intent(in)  :: wfn(*)
  !--------------------------------------------------------------------
end subroutine wfn_print_op1x
!################################################################################
