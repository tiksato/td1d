!######################################################################
subroutine general_hprod_kinetic(lfield, wfn, hwfn, ng0, ng1)

  use root_mod, only : icomp
  use const_mod, only : czero, iunit
  use field_mod, only : gauge
  use grid_mod, only : ngrid, dgrid
  use grid_mod, only : fd_order, fd_ohalf, fd_coeff_in, fd_coeff_l, fd_coeff_r, &
       & fd_coeff1_in, fd_coeff1_l, fd_coeff1_r
  use wfn_mod, only : nfun, nfcore

  implicit none
  real(kind(0d0)), intent(in) :: lfield
  complex(kind(0d0)), intent(in)    :: wfn (0:ngrid, 1:nfun)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:nfun)
  integer, intent(in) :: ng0, ng1

  complex(kind(0d0)) :: tmp
  complex(kind(0d0)) :: coeff_in(-fd_ohalf:fd_ohalf)
  complex(kind(0d0)) :: coeff_l(0:(fd_order + 1), 0:fd_ohalf)
  complex(kind(0d0)) :: coeff_r(0:(fd_order + 1), 0:fd_ohalf)
  integer :: k, kgrid, igrid, ifun, icoeff, kmax
!
  kmax = 2 * fd_ohalf + 1
  if (icomp /= 1 .or. trim(gauge) == 'L') then
     ! kinetic energy
     coeff_in(-fd_ohalf:fd_ohalf) = fd_coeff_in(-fd_ohalf:fd_ohalf)
  else
     ! kinetic energy + velocity gauge field
     coeff_in(-fd_ohalf:fd_ohalf) = fd_coeff_in (-fd_ohalf:fd_ohalf) &
               & - iunit * lfield * fd_coeff1_in(-fd_ohalf:fd_ohalf)
  end if

  do ifun = nfcore + 1, nfun

!debug     do igrid = ng0, ng1
!debug        tmp = wfn(igrid, ifun) * coeff_in(0)
!debug        do k = 1, fd_ohalf
!debug           if(igrid - k >= 0)     tmp = tmp + wfn(igrid - k, ifun) * coeff_in(-k)
!debug           if(igrid + k <= ngrid) tmp = tmp + wfn(igrid + k, ifun) * coeff_in(+k)
!debug        end do
!debug        hwfn(igrid, ifun) = hwfn(igrid, ifun) + tmp
!debug     end do

     ! inside ===============================================
     !         q <= igrid <= ngrid - q (ngrid - 2q + 1 terms)
     ! igrid - q <= kgrid <= igrid + q (        2q + 1 terms)
     ! ======================================================
!old     do igrid = max(ng0, fd_ohalf), min(ng1, ngrid - fd_ohalf)
!old        tmp = wfn(igrid, ifun) * coeff_in(0)
!old        do k = 1, fd_ohalf
!old           tmp = tmp + wfn(igrid - k, ifun) * coeff_in(-k) &
!old                   & + wfn(igrid + k, ifun) * coeff_in(+k)
!old        end do
!old        hwfn(igrid, ifun) = hwfn(igrid, ifun) + tmp
!old     end do
     do igrid = ng0, ng1
        tmp = wfn(igrid, ifun) * coeff_in(0)
        do k = 1, fd_ohalf
           if (igrid - k >= 0) then
              tmp = tmp + wfn(igrid - k, ifun) * coeff_in(-k)
           end if
           if (igrid + k <= ngrid) then
              tmp = tmp + wfn(igrid + k, ifun) * coeff_in(+k)
           end if
        end do
        hwfn(igrid, ifun) = hwfn(igrid, ifun) + tmp
     end do

!BOUNDARY
!     ! left boundary =====================
!     ! 0 <= igrid <=  q - 1 (     q terms)
!     ! 0 <= kgrid <= 2q + 1 (2q + 2 terms)
!     ! ===================================
!     do igrid = ng0, min(ng1, fd_ohalf - 1)
!        icoeff = igrid
!        tmp = czero
!        do k = 0, kmax
!           kgrid = k
!           tmp = tmp + wfn(kgrid, ifun) * coeff_l(k, icoeff)
!!           tmp = tmp + wfn(kgrid, ifun) * fd_coeff_l(k, icoeff)
!        end do
!        hwfn(igrid, ifun) = hwfn(igrid, ifun) + tmp
!     end do
!
!     ! right boundary ================================
!     ! ngrid -  q + 1 <= igrid <= ngrid (     q terms)
!     ! ngrid - 2q - 1 <= kgrid <= ngrid (2q + 2 terms)
!     ! ===============================================
!     do igrid = max(ng0, ngrid - fd_ohalf + 1), ng1
!        icoeff = ngrid - igrid
!        tmp = czero
!        do k = 0, kmax
!           kgrid = ngrid - kmax + k
!           tmp = tmp + wfn(kgrid, ifun) * coeff_r(k, icoeff)
!!           tmp = tmp + wfn(kgrid, ifun) * fd_coeff_r(k, icoeff)
!        end do
!        hwfn(igrid, ifun) = hwfn(igrid, ifun) + tmp
!     end do
!BOUNDARY

  end do

end subroutine general_hprod_kinetic
!######################################################################
