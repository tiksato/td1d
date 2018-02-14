!################################################################################
subroutine wfn_domask(rmax, dt, wfn, q0, q1, q2, imethod)

  use root_mod, only : icomp

  implicit none
  !--------------------------------------------------------------------
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in) :: dt
  integer, intent(in) :: imethod
  real(kind(0d0)), intent(inout) :: q0, q1, q2
  complex(kind(0d0)), intent(inout)  :: wfn(*)
  !--------------------------------------------------------------------
  integer :: poscic, pospt1
  integer, external :: wfn_poscic, wfn_pospt1

  poscic = wfn_poscic(imethod)
  pospt1 = wfn_pospt1(imethod)

  if (imethod == -1) then
     call x2e_domask(rmax, dt, wfn, q0, q1, q2)
  else if (imethod == 0) then
     call hf_domask(rmax, dt, wfn, q0, q1, q2)
  end if

end subroutine wfn_domask
!################################################################################
