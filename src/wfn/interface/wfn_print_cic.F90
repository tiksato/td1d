!######################################################################
subroutine wfn_print_cic(iunit, wfn, imethod)

  implicit none
  integer, intent(in) :: iunit, imethod
  complex(kind(0d0)), intent(in) :: wfn(*)

  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
  else if (imethod == 0) then
  end if

end subroutine wfn_print_cic
!######################################################################
