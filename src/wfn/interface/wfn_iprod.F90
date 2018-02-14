!######################################################################
complex(kind(0d0)) function wfn_iprod(rmax, wfnl, wfnr, imethod)

  use const_mod, only : czero
  use root_mod, only : icomp
  use wfn_mod, only : ormas

  implicit none
  real(kind(0d0)), intent(in) :: rmax
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in)  :: wfnl(*)
  complex(kind(0d0)), intent(in)  :: wfnr(*)

  integer :: poscic, pospt1
  integer, external :: wfn_poscic, wfn_pospt1
  complex(kind(0d0)), external :: x2e_iprod, hf_iprod

  poscic = wfn_poscic(imethod)
  pospt1 = wfn_pospt1(imethod)

  if (imethod == -1) then
     wfn_iprod = x2e_iprod(rmax, wfnl, wfnr)
  else if (imethod == 0) then
     wfn_iprod = hf_iprod(rmax, wfnl, wfnr)
  end if

end function wfn_iprod
!######################################################################
