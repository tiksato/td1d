!################################################################################
subroutine wfn_fft_destroy_plan(imethod)

  use fft_mod, only : fft_planf, fft_planb, fft_wfnr, fft_wfnk

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  !--------------------------------------------------------------------

#ifdef FFT_DISABLED
  stop 'wfn_fft_destroy_plan: fft disabled.'
#else
  call dfftw_destroy_plan(fft_planb)
  call dfftw_destroy_plan(fft_planf)
  deallocate(fft_wfnk)
  deallocate(fft_wfnr)
#endif

end subroutine wfn_fft_destroy_plan
!################################################################################
