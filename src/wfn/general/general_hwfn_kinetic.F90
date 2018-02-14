!######################################################################
subroutine general_hwfn_kinetic(nfun, wfn, hwfn, ng0, ng1)

  use const_mod, only : czero
  use grid_mod, only : ngrid, dgrid
  use grid_mod, only : fd_ohalf, fd_coeff_in, fd_coeff_l, fd_coeff_r

  implicit none
  integer, intent(in) :: nfun
  complex(kind(0d0)), intent(in)    :: wfn (0:ngrid, 1:*)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:*)
  integer, intent(in) :: ng0, ng1

  complex(kind(0d0)) :: tmp
  integer :: k, kgrid, igrid, ifun, icoeff, kmax
!debug
!  write(6, "('WARNING: skip kinetic operator!')")
!  return
!debug
!
  kmax = 2 * fd_ohalf + 1
!
  do ifun = 1, nfun

     ! inside ===============================================
     !         q <= igrid <= ngrid - q (ngrid - 2q + 1 terms)
     ! igrid - q <= kgrid <= igrid + q (        2q + 1 terms)
     ! ======================================================
     do igrid = max(ng0, fd_ohalf), min(ng1, ngrid - fd_ohalf)
        tmp = wfn(igrid, ifun) * fd_coeff_in(0)
        do k = 1, fd_ohalf
           tmp = tmp + wfn(igrid - k, ifun) * fd_coeff_in(-k) &
                   & + wfn(igrid + k, ifun) * fd_coeff_in(+k)
        end do
        hwfn(igrid, ifun) = hwfn(igrid, ifun) + tmp
     end do

!     write(6, "('general_hwfn_kinetic: boundaries!')")
!     ! left boundary =====================
!     ! 0 <= igrid <=  q - 1 (     q terms)
!     ! 0 <= kgrid <= 2q + 1 (2q + 2 terms)
!     ! ===================================
!     do igrid = ng0, min(ng1, fd_ohalf - 1)
!        icoeff = igrid
!        tmp = czero
!        do k = 0, kmax
!           kgrid = k
!           tmp = tmp + wfn(kgrid, ifun) * fd_coeff_l(k, icoeff)
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
!           tmp = tmp + wfn(kgrid, ifun) * fd_coeff_r(k, icoeff)
!        end do
!        hwfn(igrid, ifun) = hwfn(igrid, ifun) + tmp
!     end do

  end do

end subroutine general_hwfn_kinetic
!######################################################################
