!################################################################################
integer function wfn_pospt1(imethod)

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  !--------------------------------------------------------------------

  if (imethod == -1) then
     wfn_pospt1 = 1
  else if (imethod == 0) then
     wfn_pospt1 = 1
  end if

end function wfn_pospt1
!######################################################################
