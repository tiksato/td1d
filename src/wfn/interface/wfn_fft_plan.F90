!################################################################################
subroutine wfn_fft_plan(imethod)

  use grid_mod, only : ngrid
  use fft_mod, only : fft_planf, fft_planb, fft_wfnr, fft_wfnk
  use omp_mod, only : omp_get_num_threads, omp_get_thread_num

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  !--------------------------------------------------------------------
#ifdef FFT_DISABLED
  stop 'wfn_fft_plan: fft disabled.'
#else
  include 'fftw3.f'

  if (imethod == -1) then
     allocate(fft_wfnr(ngrid * ngrid))
     allocate(fft_wfnk(ngrid * ngrid))
     call dfftw_plan_dft_2d(fft_planf, ngrid, ngrid, fft_wfnr, fft_wfnk, fftw_forward, fftw_measure)
     call dfftw_plan_dft_2d(fft_planb, ngrid, ngrid, fft_wfnk, fft_wfnr, fftw_backward, fftw_measure)
 !    call dfftw_plan_dft_2d(fft_planf, ngrid-1, ngrid-1, fft_wfnr, fft_wfnk, fftw_forward, fftw_measure)
 !    call dfftw_plan_dft_2d(fft_planb, ngrid-1, ngrid-1, fft_wfnk, fft_wfnr, fftw_backward, fftw_measure)
  else
     allocate(fft_wfnr(ngrid))
     allocate(fft_wfnk(ngrid))
     call dfftw_plan_dft_1d(fft_planf, ngrid, fft_wfnr, fft_wfnk, fftw_forward, fftw_measure)
     call dfftw_plan_dft_1d(fft_planb, ngrid, fft_wfnk, fft_wfnr, fftw_backward, fftw_measure)
 !    call dfftw_plan_dft_1d(fft_planf, ngrid-1, fft_wfnr, fft_wfnk, fftw_forward, fftw_measure)
 !    call dfftw_plan_dft_1d(fft_planb, ngrid-1, fft_wfnk, fft_wfnr, fftw_backward, fftw_measure)
  end if
#endif

end subroutine wfn_fft_plan
!################################################################################
