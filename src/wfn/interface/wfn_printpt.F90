!######################################################################
subroutine wfn_printpt(iunit, wfn, imethod)

  implicit none
  integer, intent(in) :: iunit, imethod
  complex(kind(0d0)), intent(in) :: wfn(*)

  integer :: pospt1
  integer, external :: wfn_pospt1

  pospt1 = wfn_pospt1(imethod)
  call x2e_print(iunit, wfn(pospt1))

end subroutine wfn_printpt
!######################################################################
