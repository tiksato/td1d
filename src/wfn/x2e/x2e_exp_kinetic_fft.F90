!######################################################################
subroutine x2e_exp_kinetic_fft(dstep, wfn, expwfn, wfnr, wfnk)

  use omp_mod
  use const_mod, only : one, iunit
  use grid_mod, only : ngrid, efree2
  use fft_mod, only : fft_planf, fft_planb

  implicit none
  complex(kind(0d0)), intent(in)  :: dstep
  complex(kind(0d0)), intent(in)    :: wfn   (0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(inout) :: expwfn(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(out) :: wfnr(1:ngrid, 1:ngrid)
  complex(kind(0d0)), intent(out) :: wfnk(1:ngrid, 1:ngrid)

  integer :: kgrid, lgrid
  real(kind(0d0)) :: fac
  complex(kind(0d0)) :: expt

#ifdef FFT_DISABLED
  stop 'x2e_exp_kinetic_fft: fft disabled.'
#else
  fac = one / dble(ngrid * ngrid)

  wfnr(1:ngrid, 1:ngrid) = wfn(1:ngrid, 1:ngrid)
  call dfftw_execute(fft_planf, wfnr, wfnk)

  !$omp parallel default(shared) private(kgrid, lgrid)
  call omp_mod_thread(0, ngrid)

  do kgrid = ng0, ng1
     do lgrid = 1, kgrid - 1
        expt = exp(-iunit * dstep * efree2(lgrid - 1, kgrid - 1)) * fac
        wfnk(lgrid, kgrid) = wfnk(lgrid, kgrid) * expt
        wfnk(kgrid, lgrid) = wfnk(lgrid, kgrid)
     end do
     expt = exp(-iunit * dstep * efree2(kgrid - 1, kgrid - 1)) * fac
     wfnk(kgrid, kgrid) = wfnk(kgrid, kgrid) * expt
  end do
  !$omp end parallel

  call dfftw_execute(fft_planb, wfnk, wfnr)
  expwfn(1:ngrid, 1:ngrid) = wfnr(1:ngrid, 1:ngrid)
#endif

end subroutine x2e_exp_kinetic_fft
!######################################################################
