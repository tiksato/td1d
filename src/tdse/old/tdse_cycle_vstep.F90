!######################################################################
subroutine tdse_cycle_vstep(wfn0, wfn, imethodx, prop_type)

  use const_mod, only : zero, czero, runit
  use grid_mod, only : rmax, domask, docap
  use field_mod, only : period
  use prop_mod, only : cyctot, period, dstep
  use wfn_mod, only : domo, doci
  use fft_mod, only : dofft, fft_wfnr, fft_wfnk

  implicit none
  character(len=*), intent(in) :: prop_type
  integer, intent(in) :: imethodx
  complex(kind(0d0)), intent(in) :: wfn0(*)
  complex(kind(0d0)), intent(inout) :: wfn(*)

  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod
  complex(kind(0d0)), allocatable :: hwfn(:)
  complex(kind(0d0)), allocatable :: wfnp(:)
  integer :: len, istep, nprint
  real(kind(0d0)) :: lfield, ene, q0, q1, q2
  complex(kind(0d0)) :: dt, time

  len = wfn_size(imethodx)
  allocate(hwfn(len))
  allocate(wfnp(len))

  domo = .true.
  doci = .true.
  dt = dstep * runit
  nprint = 0

  q0 = zero
  q1 = zero
  q2 = zero

  if (dofft) call wfn_fft_plan(imethodx)

  istep = 0
  time = czero

  do

     lfield = field(time)

     ! this is not necessary for split operator propagators...
     ene = wfn_hprod(.true., lfield, wfn0, wfn, hwfn, imethodx)

     ! output
     call tdse_print(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, nprint, imethodx)

     ! propagation from t to t + dt
     call prop(time, dt, len, wfn0, wfn, hwfn, imethodx, prop_type)

     if (domask) then
        ! mask
        call wfn_domask(rmax, dt, wfn, q0, q1, q2, imethodx)
     else if (docap) then
!nyi        ! cap
!nyi        call wfn_capped(rmax, dt, wfn, wfnp, q0, q1, q2, imethodx)
     end if

     if (dble(time) > cyctot*period) exit
     istep = istep + 1
     time = time + dt

  end do

  if (dofft) call wfn_fft_destroy_plan(imethodx)

  deallocate(wfnp)
  deallocate(hwfn)

end subroutine tdse_cycle_vstep
!######################################################################
