!################################################################################
subroutine tdse_summary(wfn, imethod)

  use io_mod, only : iostdo
  use root_mod, only : iprint

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(in) :: wfn(*)
  !--------------------------------------------------------------------

  if (iprint > 2) then
     write(iostdo, "('# final wavefunction')")
     call wfn_print(iostdo, wfn, imethod)
  end if

end subroutine tdse_summary
!################################################################################
