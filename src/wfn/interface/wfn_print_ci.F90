!################################################################################
subroutine wfn_print_ci(iout, istep, time, wfn, imethod)

  use const_mod, only : zero, two
  use grid_mod, only : ngrid, x

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: iout, istep, imethod
  complex(kind(0d0)), intent(in) :: time
  complex(kind(0d0)), intent(in) :: wfn(*)
  !--------------------------------------------------------------------
  integer :: poscic, pospt1
  integer, external :: wfn_poscic, wfn_pospt1

  poscic = wfn_poscic(imethod)
  pospt1 = wfn_pospt1(imethod)

  if (imethod == -1) then
  else if (imethod == 0) then
  end if

end subroutine wfn_print_ci
!################################################################################
