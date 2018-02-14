!######################################################################
subroutine wfn_print_den1(io, wfn, imethod)

  implicit none
  integer, intent(in) :: io, imethod
  complex(kind(0d0)), intent(in) :: wfn(*)

  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
  else if (imethod == 0) then
  end if

end subroutine wfn_print_den1
!######################################################################
