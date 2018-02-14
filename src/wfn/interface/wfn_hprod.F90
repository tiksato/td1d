!################################################################################
real(kind(0d0)) function wfn_hprod(calene, lfield, dtime, wfn0, wfn, hwfn, imethod)

  use const_mod, only : zero
  use root_mod, only : icomp
  use wfn_mod, only : ormas, tdcc
  use wfn_mod, only : nfcore, ncore, sep_fc

  implicit none
  !--------------------------------------------------------------------
  logical, intent(in) :: calene
  real(kind(0d0)), intent(in) :: lfield
  integer, intent(in) :: imethod
  real(kind(0d0)), intent(in) :: dtime
  complex(kind(0d0)), intent(in)  :: wfn0(*), wfn(*)
  complex(kind(0d0)), intent(out) :: hwfn(*)
  !--------------------------------------------------------------------
  integer :: poscic, pospt1
  integer, external :: wfn_poscic, wfn_pospt1
  real(kind(0d0)), external :: hf_hprod

  poscic = wfn_poscic(imethod)
  pospt1 = wfn_pospt1(imethod)

  if (imethod == -1) then
!nyi     wfn_hprod = x2e_hprod(.true., calene, lfield, wfn, hwfn)
  else if (imethod == 0) then
     wfn_hprod = hf_hprod(.true., calene, lfield, dtime, wfn0, wfn, hwfn)
  end if

end function wfn_hprod
!################################################################################
