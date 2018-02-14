!######################################################################
!nyisubroutine tdse_printpt(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, nprint, imethod)
!nyi
!nyi  use root_mod, only : name
!nyi  use const_mod, only : one, four, runit
!nyi  use grid_mod, only : rmax, x0, xmask
!nyi  use field_mod, only : period
!nyi  use prop_mod, only : totstep
!nyi
!nyi  implicit none
!nyi  integer, intent(in) :: istep
!nyi  complex(kind(0d0)), intent(in) :: time
!nyi  real(kind(0d0)), intent(inout) :: lfield, ene, q0, q1, q2
!nyi  complex(kind(0d0)), intent(in) :: wfn0(*)
!nyi  complex(kind(0d0)), intent(in) :: wfn(*)
!nyi  integer, intent(inout) :: nprint
!nyi  integer, intent(in) :: imethod
!nyi
!nyi  integer :: lenx
!nyi  integer, external :: wfn_size
!nyi!debug  integer :: poscic
!nyi!debug  integer, external :: cas_poscic
!nyi  complex(kind(0d0)), allocatable :: wfn2e0(:)
!nyi  complex(kind(0d0)), allocatable :: wfn2e(:)
!nyi
!nyi  lenx = wfn_size(-1)
!nyi  allocate(wfn2e0(lenx))
!nyi  allocate(wfn2e(lenx))
!nyi
!nyi  call wfn_wfn2e(wfn0, wfn2e0, imethod)
!nyi  call wfn_wfn2e(wfn, wfn2e, imethod)
!nyi!debug  poscic = cas_poscic()
!nyi!debug  call cas_wfn2e(wfn0, wfn0(poscic), wfn2e0)
!nyi!debug  call cas_wfn2e(wfn, wfn(poscic), wfn2e)
!nyi
!nyi!old  wfac = runit / sqrt(runit + norm0)
!nyi!old  call wfn_wfn2e(wfn0, wfn2e0, imethod)
!nyi!old  call zaxpy(lenx, runit, wfnpt0, 1, wfn2e0, 1)
!nyi!old  call zscal(lenx, wfac, wfn2e0, 1)
!nyi!old  wfac = runit / sqrt(runit + norm)
!nyi!old  call wfn_wfn2e(wfn, wfn2e, imethod)
!nyi!old  call zaxpy(lenx, runit, wfnpt, 1, wfn2e, 1)
!nyi!old  call zscal(lenx, wfac, wfn2e, 1)  
!nyi
!nyi  call tdse_print(istep, time, lfield, ene, q0, q1, q2, wfn2e0, wfn2e, nprint, -1)
!nyi
!nyi  deallocate(wfn2e)
!nyi  deallocate(wfn2e0)
!nyi
!nyiend subroutine tdse_printpt
!######################################################################
