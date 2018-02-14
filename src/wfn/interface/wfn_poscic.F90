!################################################################################
integer function wfn_poscic(imethod)

  use wfn_mod, only : ormas

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  !--------------------------------------------------------------------

  if (imethod == -1) then
     wfn_poscic = 1
  else if (imethod == 0) then
     wfn_poscic = 1
  end if

end function wfn_poscic
!################################################################################
