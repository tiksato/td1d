!################################################################################
real(kind(0d0)) function wfn_prop_p2(lfield, dstep, wfn, hwfn, imethod)

  use const_mod, only : zero
  use root_mod, only : icomp
  use wfn_mod, only : ormas

  implicit none
  !--------------------------------------------------------------------
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in) :: dstep
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout)  :: wfn(*)
  complex(kind(0d0)), intent(out) :: hwfn(*)
  !--------------------------------------------------------------------

end function wfn_prop_p2
!################################################################################
