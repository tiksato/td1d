!######################################################################
subroutine wfn_print_ovlp(io, rmax, wfn, norm, imethod)

  implicit none
  integer, intent(in) :: io, imethod
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: wfn(*)
  real(kind(0d0)), intent(out) :: norm

  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
  else if (imethod == 0) then
     call hf_print_ovlp(io, rmax, wfn, norm)
  end if

end subroutine wfn_print_ovlp
!######################################################################
