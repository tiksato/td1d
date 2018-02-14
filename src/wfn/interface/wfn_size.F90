!######################################################################
integer function wfn_size(imethod)

  use root_mod, only : icomp
  use wfn_mod, only : ormas

  implicit none
  integer, intent(in) :: imethod

  integer, external :: x2e_size, hf_size

  if (imethod == -1) then
     wfn_size = x2e_size()
  else if (imethod == 0) then
     wfn_size = hf_size()
  end if

end function wfn_size
!######################################################################
