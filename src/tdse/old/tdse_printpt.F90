!######################################################################
subroutine tdse_printpt(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, nprint, imethodx)

  use root_mod, only : name, nprintwfn, printwfn
  use const_mod, only : one, four, runit
  use grid_mod, only : rmax, x0, xmask
  use field_mod, only : period
  use prop_mod, only : totstep

  implicit none
  integer, intent(in) :: istep
  complex(kind(0d0)), intent(in) :: time
  real(kind(0d0)), intent(inout) :: lfield, ene, q0, q1, q2
  complex(kind(0d0)), intent(in) :: wfn0(*)
  complex(kind(0d0)), intent(in) :: wfn(*)
  integer, intent(inout) :: nprint
  integer, intent(in) :: imethodx

  integer :: lenx
  integer, external :: wfn_size
!debug  integer :: poscic
!debug  integer, external :: cas_poscic
  complex(kind(0d0)), allocatable :: wfn2e0(:)
  complex(kind(0d0)), allocatable :: wfn2e(:)

  lenx = wfn_size(-1)
  allocate(wfn2e0(lenx))
  allocate(wfn2e(lenx))

  call wfn_wfn2e(wfn0, wfn2e0, imethodx)
  call wfn_wfn2e(wfn, wfn2e, imethodx)
!debug  poscic = cas_poscic()
!debug  call cas_wfn2e(wfn0, wfn0(poscic), wfn2e0)
!debug  call cas_wfn2e(wfn, wfn(poscic), wfn2e)

!old  wfac = runit / sqrt(runit + norm0)
!old  call wfn_wfn2e(wfn0, wfn2e0, imethodx)
!old  call zaxpy(lenx, runit, wfnpt0, 1, wfn2e0, 1)
!old  call zscal(lenx, wfac, wfn2e0, 1)
!old  wfac = runit / sqrt(runit + norm)
!old  call wfn_wfn2e(wfn, wfn2e, imethodx)
!old  call zaxpy(lenx, runit, wfnpt, 1, wfn2e, 1)
!old  call zscal(lenx, wfac, wfn2e, 1)  

  call tdse_print(istep, time, lfield, ene, q0, q1, q2, wfn2e0, wfn2e, nprint, -1)

  deallocate(wfn2e)
  deallocate(wfn2e0)

end subroutine tdse_printpt
!######################################################################
