!################################################################################
subroutine wfn_v2prod(dtime, wfn, hwfn, imethod)

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  real(kind(0d0)), intent(in) :: dtime
  complex(kind(0d0)), intent(in)  :: wfn(*)
  complex(kind(0d0)), intent(out) :: hwfn(*)
  !--------------------------------------------------------------------
  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == 0) then
     call hf_v2prod(dtime, wfn, hwfn)
  end if

end subroutine wfn_v2prod
!################################################################################
