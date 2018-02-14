!################################################################################
module fft_mod

  implicit none

  logical :: dofft
  integer :: fft_planf(8)
  integer :: fft_planb(8)
  integer :: fft_nthreads

  complex(kind(0d0)), allocatable :: fft_wfnr(:)
  complex(kind(0d0)), allocatable :: fft_wfnk(:)

end module fft_mod
!################################################################################
