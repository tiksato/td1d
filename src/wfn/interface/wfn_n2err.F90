!######################################################################
real(kind(0d0)) function wfn_n2err(wfn1, wfn2, imethod)

  use wfn_mod, only : ormas

  implicit none
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in)  :: wfn1(*)
  complex(kind(0d0)), intent(in)  :: wfn2(*)

  integer :: poscic
  integer, external :: wfn_poscic
  real(kind(0d0)), external :: x2e_n2err, hf_n2err

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     wfn_n2err = x2e_n2err(wfn1, wfn2)
  else if (imethod == 0) then
     wfn_n2err = hf_n2err(wfn1, wfn2)
  end if

end function wfn_n2err
!######################################################################
