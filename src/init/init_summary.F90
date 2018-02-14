!################################################################################
subroutine init_summary(wfn, imethodx)

  use io_mod, only : iostdo
  use root_mod, only : iprint

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethodx
  complex(kind(0d0)), intent(in) :: wfn(*)
  !--------------------------------------------------------------------

  call wfn_exps2(wfn, imethodx)

  if (iprint > 2) then
     write(iostdo, "('# initial wavefunction')")
     call wfn_print(iostdo, wfn, imethodx)
  end if

end subroutine init_summary
!################################################################################
