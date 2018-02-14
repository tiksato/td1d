!######################################################################
subroutine general_hwfn_kinetic_fft(nfun, wfn, hwfn)

  use const_mod, only : one, two
  use grid_mod, only : ngrid, efree
  use fft_mod, only : fft_planf, fft_planb, fft_wfnr, fft_wfnk

  implicit none
  integer, intent(in) :: nfun
  complex(kind(0d0)), intent(in)    :: wfn (0:ngrid, 1:*)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:*)

  integer :: kgrid, ifun
  real(kind(0d0)) :: norm
  complex(kind(0d0)), allocatable :: eigt(:)

#ifdef FFT_DISABLED
  stop 'general_hwfn_kinetic_fft: fft disabled.'
#else
!debug
!  write(6, "('WARNING: skip kinetic operator!')")
!  return
!debug
  allocate(eigt(1:ngrid))

  norm = one / dble(ngrid)
  do kgrid = 1, ngrid
     eigt(kgrid) = efree(kgrid - 1) * norm
  end do

  do ifun = 1, nfun

     call util_zcopy(ngrid, wfn(1, ifun), 1, fft_wfnr, 1)
     call dfftw_execute(fft_planf, fft_wfnr, fft_wfnk)

     do kgrid = 1, ngrid
        fft_wfnk(kgrid) = fft_wfnk(kgrid) * eigt(kgrid)
     end do

     call dfftw_execute(fft_planb, fft_wfnk, fft_wfnr)
     call util_zcopy(ngrid, fft_wfnr, 1, hwfn(1, ifun), 1)
        
  end do

  deallocate(eigt)
#endif

end subroutine general_hwfn_kinetic_fft
!######################################################################
