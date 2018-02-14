!######################################################################
subroutine general_hprod_kinetic_fft(nfun, nfroz, ndim, lfield, wfn, hwfn)

  use root_mod, only : icomp
  use const_mod, only : one, two
  use field_mod, only : gauge
  use grid_mod, only : ngrid, efree, p
  use fft_mod, only : fft_planf, fft_planb, fft_wfnr, fft_wfnk

  implicit none
  integer, intent(in) :: nfun, nfroz, ndim
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in)    :: wfn (0:ngrid, 1:nfun, 1:ndim)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:nfun, 1:ndim)

  integer :: kgrid, ifun, idim
  real(kind(0d0)) :: norm
  complex(kind(0d0)), allocatable :: eigt(:)

#ifdef FFT_DISABLED
  stop 'general_hprod_kinetic_fft: fft disabled.'
#else
  allocate(eigt(1:ngrid))

  norm = one / dble(ngrid)
  if (icomp /= 1 .or. trim(gauge) == 'L') then
     ! kinetic energy
     do kgrid = 1, ngrid
        eigt(kgrid) = efree(kgrid - 1) * norm
     end do
  else
     ! kinetic energy + velocity gauge field
     do kgrid = 1, ngrid
        eigt(kgrid) = (efree(kgrid - 1) + lfield * p(kgrid - 1)) * norm
     end do
  end if

  do idim = 1, ndim
     do ifun = nfroz + 1, nfun

        call util_zcopy(ngrid, wfn(1, ifun, idim), 1, fft_wfnr, 1)
        call dfftw_execute(fft_planf, fft_wfnr, fft_wfnk)

        do kgrid = 1, ngrid
           fft_wfnk(kgrid) = fft_wfnk(kgrid) * eigt(kgrid)
        end do

        call dfftw_execute(fft_planb, fft_wfnk, fft_wfnr)
        call util_zcopy(ngrid, fft_wfnr, 1, hwfn(1, ifun, idim), 1)
           
     end do
  end do

  deallocate(eigt)
#endif

end subroutine general_hprod_kinetic_fft
!######################################################################
