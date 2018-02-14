!################################################################################
real(kind(0d0)) function wfn_diagp(calene, lfield, wfn0, wfn, imethod)

  implicit none
  !--------------------------------------------------------------------
  logical, intent(in) :: calene
  real(kind(0d0)), intent(in) :: lfield
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: wfn0(*)
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     stop 'x2e_diagp nyi.'
  else if (imethod == 0) then
     stop 'hf_diagp nyi.'
  end if

end function wfn_diagp
!################################################################################
