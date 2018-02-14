!######################################################################
real(kind(0d0)) function wfn_vprod(calene, lfield, wfn0, wfn, hwfn, imethod)

  use const_mod, only : zero
  use root_mod, only : icomp

  implicit none
  logical, intent(in) :: calene
  real(kind(0d0)), intent(in) :: lfield
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in)  :: wfn0(*), wfn(*)
  complex(kind(0d0)), intent(out) :: hwfn(*)

  integer :: poscic, pospt1
  integer, external :: wfn_poscic, wfn_pospt1
  real(kind(0d0)), external :: x2e_hprod, hf_hprod

  poscic = wfn_poscic(imethod)
  pospt1 = wfn_pospt1(imethod)

  if (imethod == -1) then
     wfn_vprod = x2e_hprod(.false., calene, lfield, wfn, hwfn)
  else if (imethod == 0) then
     wfn_vprod = hf_hprod(.false., calene, lfield, wfn0, wfn, hwfn)
  end if

end function wfn_vprod
!######################################################################
