!################################################################################
subroutine init_cycle_diagcip2p1q(wfn, imethod)

  use mol_mod, only : enen
  use io_mod, only : iostdo
  use const_mod, only : czero, iunit, zero, one
  use thresh_mod, only : threne, thrgbt
  use init_mod, only : maxcyc, distep1, distep2, distep3, prop_type_init
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

  logical :: needp, needq
  integer :: len, istep
  real(kind(0d0)) :: ene, enep, ene1, ene2, ene3, lfield, test, test1, test2, test3
  complex(kind(0d0)) :: dt, dt1, dt2, time

  len = wfn_size(imethod)
  allocate(hwfn(len))

  time = dt * (istep - 1)
  lfield = field(time)

  dt1 = - distep2 * iunit
  dt2 = - min(distep1, distep3) * iunit
  enep = zero

  needp = .true.
  needq = .true.
  test1 = one
  test2 = one
  test3 = one
  max_grad(1) = zero
  max_grad(2) = zero
  max_grad(3) = zero

  do istep = 1, maxcyc

!nyi     ! ci and a-a rotations update
!nyi     doci = .not. noci
!nyi     domo = .true.
!nyi     domop1 = .not. nolag1
!nyi     domop2 = .not. nolag2
!nyi     domop = domop1 .or. domop2
!nyi     domoq = .false.
!nyi     ene_mo = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!nyi     ene_ci = wfn_diag(.true., lfield, wfn, wfn, imethod)
!nyi     call wfn_ort(wfn, imethod)
!nyi     if (cionly) then
!nyi        ene = ene_ci
!nyi        exit
!nyi     end if
     ! ci update
     if (needp) then
        doci = .not. noci
        domo = .false.
        domop1 = .false.
        domop2 = .false.
        domop = .false.
        domoq = .false.
        ene1 = wfn_diag(.true., lfield, wfn, wfn, imethod)
        call wfn_ort(wfn, imethod)
        if (cionly) then
           ene = ene1
           exit
        end if
   
        ! a-a rotations update
        dt = dt1
        doci = .false.
        domo = .true.
        domop1 = .false.
        domop2 = .not. nolag2
        domop  = domop2
        domoq  = .false.
        ene1 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
        call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
        call wfn_ort(wfn, imethod)
     else
        ene1 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
        ene2 = ene1
     end if
!nyi

     ! c-a and v-o rotations update
     dt = dt2
     doci = .false.
     domo = .true.
     domop1 = .not. nolag1
     domop2 = .false.
     domop  = domop1
     domoq  = .not. nolag3
     ene2 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
     call wfn_ort(wfn, imethod)
     ene3 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)

     ene = ene3
     test = ene - enep
     test1 = ene1 - enep
     test2 = ene2 - ene1
     test3 = ene3 - ene2
!nyi     needp = abs(test1) > abs(test3) .or. &
!nyi           & abs(test1) < abs(test3) .and. abs(test1) > abs(test3) * 100
!nyi     needq = abs(test3) > abs(test1) .or. &
!nyi           & abs(test3) < abs(test1) .and. abs(test3) > abs(test1) * 100

!     if (test < zero .and. test1 < zero .and. test2 < zero .and. test3 < zero) then
        write(iostdo, "('# cycle ', i10, 3f20.10, 6e15.5)") &
             & istep, ene1, ene2, ene3, test1, test2, test3, max_grad(1:3)
!     else
!        write(iostdo, "('# cycle ', i10, 2f20.10, 4e20.10, ' <-- WARNING')") &
!             & istep, ene_mo, ene_ci, test, test1, test2, test3, max_grad(1:3)
!     end if

     if(abs(test) < threne .and. &
          & (nolag1 .or. max_grad(1) < thrgbt) .and. &
          & (nolag2 .or. max_grad(2) < thrgbt) .and. &
          & (nolag3 .or. max_grad(3) < thrgbt)) exit

     enep = ene

  end do

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy

  deallocate(hwfn)

end subroutine init_cycle_diagcip2p1q
!################################################################################
