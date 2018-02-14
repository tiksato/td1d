!######################################################################
subroutine general_exp_kinetic_fft(dstep, nfun, nfroz, ndim, wfn, expwfn)

  use const_mod, only : one, iunit
  use grid_mod, only : ngrid, efree
  use fft_mod, only : fft_planf, fft_planb, fft_wfnr, fft_wfnk

  implicit none
  complex(kind(0d0)), intent(in)  :: dstep
  integer, intent(in) :: nfun, nfroz, ndim
  complex(kind(0d0)), intent(in)  :: wfn   (0:ngrid, 1:nfun, 1:ndim)
  complex(kind(0d0)), intent(out) :: expwfn(0:ngrid, 1:nfun, 1:ndim)

  integer :: kgrid, ifun, idim
  real(kind(0d0)) :: norm
  complex(kind(0d0)), allocatable :: expt(:)

#ifdef FFT_DISABLED
  stop 'general_exp_kinetic_fft: fft disabled.'
#else
  allocate(expt(1:ngrid))

  norm = one / dble(ngrid)
  do kgrid = 1, ngrid
     expt(kgrid) = exp(-iunit * dstep * efree(kgrid - 1)) * norm
  end do

  do idim = 1, ndim
     do ifun = nfroz + 1, nfun

        call util_zcopy(ngrid, wfn(1, ifun, idim), 1, fft_wfnr, 1)
        call dfftw_execute(fft_planf, fft_wfnr, fft_wfnk)

        do kgrid = 1, ngrid
           fft_wfnk(kgrid) = fft_wfnk(kgrid) * expt(kgrid)
        end do

        call dfftw_execute(fft_planb, fft_wfnk, fft_wfnr)
        call util_zcopy(ngrid, fft_wfnr, 1, expwfn(1, ifun, idim), 1)
           
     end do
  end do

  deallocate(expt)
#endif

end subroutine general_exp_kinetic_fft
!######################################################################
