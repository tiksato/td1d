!################################################################################
subroutine init_cycle_diag(wfn, imethod)

  use mol_mod, only : enen
  use io_mod, only : iostdo
  use const_mod, only : iunit, zero, czero
  use thresh_mod, only : threne, thrgbt
  use init_mod, only : maxcyc, distepmo, prop_type_init
  use wfn_mod, only : cionly, doci, domoq, domop1, domop2, domop, domo, &
       & postapsg, energy, max_grad, noci, nomo, nolag1, nolag2, nolag3

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_diagp, wfn_hprod
  complex(kind(0d0)), allocatable :: hwfn(:)

  integer :: len, istep
  real(kind(0d0)) :: ene, ene1, lfield, test
  complex(kind(0d0)) :: dt, time

  len = wfn_size(imethod)
  allocate(hwfn(len))

  time = czero
  lfield = field(time)

  dt = - distepmo * iunit
  ene1 = zero

  max_grad(1) = zero
  max_grad(2) = zero
  max_grad(3) = zero

  do istep = 1, maxcyc

     ! augmented Hessian
     doci = .not. noci
     domo = .not. nomo
     domop1 = .not. nolag1
     domop2 = .not. nolag2
     domop = domop1 .or. domop2
     domoq = .not. nolag3
     ene = wfn_diagp(.true., lfield, wfn, wfn, imethod)
     call wfn_ort(wfn, imethod)
!     ene = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)

     if (cionly) exit

!     ! mo update
!     doci = .false.
!     domo = .not. nomo
!     domop1 = .false.
!     domop2 = .false.
!     domop  = .false.
!     domoq  = .not. nolag3
!     ene_ci = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
!     call wfn_ort(wfn, imethod)

     test = ene - ene1

     if (test < zero) then
        write(iostdo, "('# cycle ', i10, f20.10, 4e20.10)") &
             & istep, ene, test, max_grad(1:3)
     else
        write(iostdo, "('# cycle ', i10, f20.10, 4e20.10, ' <-- WARNING')") &
             & istep, ene, test, max_grad(1:3)
     end if

     if(abs(test) < threne .and. &
          & (nolag1 .or. max_grad(1) < thrgbt) .and. &
          & (nolag2 .or. max_grad(2) < thrgbt) .and. &
          & (nolag3 .or. max_grad(3) < thrgbt)) exit

     ene1 = ene

  end do

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy

  deallocate(hwfn)

end subroutine init_cycle_diag
!################################################################################
