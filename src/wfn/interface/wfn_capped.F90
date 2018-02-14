!################################################################################
subroutine wfn_capped(rmax, wfn, wfnp, q0, q1, q2, imethod)

  implicit none
  !--------------------------------------------------------------------
  real(kind(0d0)), intent(in) :: rmax
  integer, intent(in) :: imethod
  real(kind(0d0)), intent(inout) :: q0, q1, q2
  complex(kind(0d0)), intent(in)  :: wfn(*)
  complex(kind(0d0)), intent(in)  :: wfnp(*)
  !--------------------------------------------------------------------
  integer :: poscic
  integer, external :: wfn_poscic

  poscic = wfn_poscic(imethod)

  if (imethod == -1) then
     call x2e_capped(rmax, wfn, wfnp, q0, q1, q2)
  else if (imethod == 0) then
     call hf_capped(rmax, wfn, wfnp, q0, q1, q2)
  end if

end subroutine wfn_capped
!################################################################################
