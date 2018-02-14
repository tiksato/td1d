!######################################################################
subroutine x2e_hprod_kinetic_fft(lfield, wfn, hwfn, wfnr, wfnk)

  use root_mod, only : icomp
  use const_mod, only : one, two, runit, czero
  use field_mod, only : gauge
  use grid_mod, only : ngrid, efree2, p2
  use fft_mod, only : fft_planf, fft_planb
  use mol_mod, only : mult

  implicit none
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in)    :: wfn (0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 0:ngrid)
  complex(kind(0d0)), intent(out) :: wfnr(1:ngrid, 1:ngrid)
  complex(kind(0d0)), intent(out) :: wfnk(1:ngrid, 1:ngrid)

  integer :: kgrid, lgrid
  real(kind(0d0)) :: fac
  complex(kind(0d0)) :: t12, mfac

!nyi
  if (mult == 1) stop 'x2e_hprod_kinetic_fft has bug for triplet systems.'
!nyi
 
#ifdef FFT_DISABLED
  stop 'x2e_hprod_kinetic_fft: fft disabled.'
#else
  fac = one / dble(ngrid * ngrid)
  mfac = (-runit) ** mult

  wfnr(1:ngrid, 1:ngrid) = wfn(1:ngrid, 1:ngrid)
  call dfftw_execute(fft_planf, wfnr, wfnk)

  if (icomp /= 1 .or. trim(gauge) == 'L') then
     ! kinetic energy
     do kgrid = 1, ngrid
        do lgrid = 1, kgrid - 1
           t12 = efree2(lgrid - 1, kgrid - 1) * fac
           wfnk(lgrid, kgrid) = wfnk(lgrid, kgrid) * t12
           wfnk(kgrid, lgrid) = wfnk(lgrid, kgrid) * mfac
        end do
        if (mult == 0) then
           t12 = efree2(kgrid - 1, kgrid - 1) * fac
           wfnk(kgrid, kgrid) = wfnk(kgrid, kgrid) * t12
        else
           wfnk(kgrid, kgrid) = czero
        end if
     end do
  else
     ! kinetic energy + velocity gauge field
     do kgrid = 1, ngrid
        do lgrid = 1, kgrid - 1
           t12 = (efree2(lgrid - 1, kgrid - 1) + p2(lgrid - 1, kgrid - 1) * lfield) * fac
           wfnk(lgrid, kgrid) = wfnk(lgrid, kgrid) * t12
           wfnk(kgrid, lgrid) = wfnk(lgrid, kgrid) * mfac
        end do
        if (mult == 0) then
           t12 = (efree2(kgrid - 1, kgrid - 1) + p2(kgrid - 1, kgrid - 1) * lfield) * fac
           wfnk(kgrid, kgrid) = wfnk(kgrid, kgrid) * t12
        else
           wfnk(kgrid, kgrid) = czero
        end if
     end do
  end if

  call dfftw_execute(fft_planb, wfnk, wfnr)
  hwfn(1:ngrid, 1:ngrid) = wfnr(1:ngrid, 1:ngrid)
#endif

end subroutine x2e_hprod_kinetic_fft
!######################################################################
