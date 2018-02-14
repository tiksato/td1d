!######################################################################
real(kind(0d0)) function wfn_op1e(iop, rmax, wfn, imethod)

  use const_mod, only : zero
  use root_mod, only : icomp

  implicit none
  integer, intent(in) :: iop, imethod
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in)  :: wfn(*)

  integer :: poscic, pospt1
  integer, external :: wfn_poscic, wfn_pospt1
  real(kind(0d0)), external :: x2e_op1e, hf_op1e

  poscic = wfn_poscic(imethod)
  pospt1 = wfn_pospt1(imethod)  

  if (imethod == -1) then
     wfn_op1e = x2e_op1e(iop, rmax, wfn)
  else if (imethod == 0) then
     wfn_op1e = hf_op1e(iop, rmax, wfn)
  end if

end function wfn_op1e
!######################################################################
