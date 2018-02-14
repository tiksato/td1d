!################################################################################
subroutine init_cycle_propcip2p1q(wfn, imethod)

  use mol_mod, only : enen
  use io_mod, only : iostdo
  use const_mod, only : iunit, zero, czero
  use thresh_mod, only : threne, thrgbt
  use init_mod, only : maxcyc, distepci, distep1, distep2, distep3, prop_type_init
  use wfn_mod, only : doci, domoq, domop1, domop2, domop, domo, &
       & postapsg, energy, max_grad, noci, nomo, nolag1, nolag2, nolag3

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod
  complex(kind(0d0)), allocatable :: hwfn(:)

  integer :: len, istep
  real(kind(0d0)) :: ene, ene1, ene_ci, ene_mo, lfield, test, test1, test2
  complex(kind(0d0)) :: dt, dt1, dt2, time

  len = wfn_size(imethod)
  allocate(hwfn(len))

  time = czero
  lfield = field(time)

  dt1 = - min(distepci, distep2) * iunit
  dt2 = - min(distep1, distep3) * iunit
  ene1 = zero

  max_grad(1) = zero
  max_grad(2) = zero
  max_grad(3) = zero

  do istep = 1, maxcyc

     ! ci and a-a rotations update
     dt = dt1
     doci = .not. noci
     domo = .true.
     domop1 = .false.
     domop2 = .not. nolag2
     domop  = domop2
     domoq  = .false.
     ene_mo = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
     call wfn_ort(wfn, imethod)
!debug
!write(6, "(' after ci update. ')")
!debug

     ! c-a and v-o rotations update
     dt = dt2
     doci = .false.
     domo = .true.
     domop1 = .not. nolag1
     domop2 = .false.
     domop  = domop1
     domoq  = .not. nolag3
     ene_ci = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
     call wfn_ort(wfn, imethod)
!debug
!write(6, "(' after mo update. ')")
!debug

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

  if (postapsg) then
     call wfn_apsg(wfn, imethod)
     ene = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
  end if

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy

  deallocate(hwfn)

end subroutine init_cycle_propcip2p1q
!################################################################################
