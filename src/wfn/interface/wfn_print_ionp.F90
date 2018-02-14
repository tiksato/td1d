!################################################################################
subroutine wfn_print_ionp(iout, istep, time, rmax, wfn0, wfn, imethod)

  use wfn_mod, only : ormas

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: iout, istep, imethod
  complex(kind(0d0)), intent(in) :: time
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in)  :: wfn0(*), wfn(*)
  !--------------------------------------------------------------------
  integer :: poscic, pospt1
  integer, external :: wfn_poscic, wfn_pospt1

  poscic = wfn_poscic(imethod)
  pospt1 = wfn_pospt1(imethod)

  if (imethod == -1) then
  else if (imethod == 0) then
  end if

end subroutine wfn_print_ionp
!################################################################################
