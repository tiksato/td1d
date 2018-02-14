!################################################################################
subroutine init_cycle_diagprop(wfn, imethod)

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
  real(kind(0d0)), external :: wfn_diag, wfn_hprod
  complex(kind(0d0)), allocatable :: hwfn(:)

  integer :: len, istep
  real(kind(0d0)) :: ene, ene1, ene_ci, ene_mo, lfield, test, test1, test2
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

     ! ci update
     doci = .not. noci
     domo = .false.
     domop1 = .false.
     domop2 = .false.
     domop = .false.
     domoq = .false.
     ene_mo = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     ene_ci = wfn_diag(.true., lfield, wfn, wfn, imethod)
     call wfn_ort(wfn, imethod)
     if (cionly) then
        ene = ene_ci
        exit
     end if

     ! mo update
     doci = .false.
     domo = .true.
     domop1 = .not. nolag1
     domop2 = .not. nolag2
     domoq  = .not. nolag3
     domop  = domop1 .or. domop2
     ene_ci = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
     call wfn_ort(wfn, imethod)

     ene = ene_ci
     test1 = ene - ene1
     test2 = ene_ci - ene_mo
     test = test1

     if (test1 < zero .and. test2 < zero) then
        write(iostdo, "('# cycle ', i10, 2f20.10, 4e20.10)") &
             & istep, ene_mo, ene_ci, test, max_grad(1:3)
     else
        write(iostdo, "('# cycle ', i10, 2f20.10, 4e20.10, ' <-- WARNING')") &
             & istep, ene_mo, ene_ci, test, max_grad(1:3)
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

end subroutine init_cycle_diagprop
!################################################################################
