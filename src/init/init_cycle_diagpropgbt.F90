!################################################################################
subroutine init_cycle_diagpropgbt(wfn, imethod)

  use mol_mod, only : enen
  use io_mod, only : iostdo
  use thresh_mod, only : threne
  use const_mod, only : iunit, zero
  use init_mod, only : maxcyc, distep, prop_type_init
  use wfn_mod, only : cionly, doci, domoq, domop1, domop2, domop, domo, &
       & postapsg, energy, max_grad, noci, nomo, nolag1, nolag2, nolag3

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_diag, wfn_hprod, wfn_occgbt
  complex(kind(0d0)), allocatable :: hwfn(:)

  integer :: len, istep
  real(kind(0d0)) :: ene, ene1, ene_ci, ene_mo1, ene_mo2, lfield, test
  complex(kind(0d0)) :: dt, time

  len = wfn_size(imethod)
  allocate(hwfn(len))

  dt = - distep * iunit
  ene1 = zero

  max_grad(1) = zero
  max_grad(2) = zero
  max_grad(3) = zero

  do istep = 1, maxcyc

     time = dt * (istep - 1)
     lfield = field(time)

     ! ci update
     doci = .true.
     domo = .false.
     domop1 = .false.
     domop2 = .false.
     domop = .false.
     domoq = .false.
     ene_ci = wfn_diag(.true., lfield, wfn, wfn, imethod)
     call wfn_ort(wfn, imethod)
     ene_ci = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     if (cionly) then
        ene = ene_ci
        exit
     end if

     ! mo update: occ-vir  #####old     nolagx = .true. #####
     doci = .false.
     domo = .true.
     domop1 = .false.
     domop2 = .false.
     domop  = .false.
     domoq  = .not. nolag3
     ene_mo1 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
     call wfn_ort(wfn, imethod)
     ene_mo1 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)

     ! mo update: occ-occ  #####old     nolagx = .false. #####
     doci = .false.
     domo = .true.
     domop1 = .not. nolag1
     domop2 = .not. nolag2
     domop  = domop1 .or. domop2
     domoq  = .false.
     ene_mo2 = wfn_occgbt(.true., lfield, wfn, wfn, imethod)
     call wfn_ort(wfn, imethod)
     ene_mo2 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)

     ene = ene_mo2
     test = ene - ene1
     write(iostdo, "('# cycle ', i10, 3f20.10, e20.10)") istep, ene_ci, ene_mo1, ene_mo2, test
     if(abs(test) < threne) then
        exit
     else if (test > zero) then
        write(iostdo, "('# WARNING: energy raised.')")
     end if
     ene1 = ene

  end do

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy

  deallocate(hwfn)

end subroutine init_cycle_diagpropgbt
!################################################################################
!nyi!################################################################################
!nyisubroutine init_cycle_diagpropgbt2(wfn, imethod)
!nyi
!nyi  use mol_mod, only : enen
!nyi  use io_mod, only : iostdo
!nyi  use thresh_mod, only : threne
!nyi  use const_mod, only : iunit, zero
!nyi  use init_mod, only : maxcyc, distep, prop_type_init
!nyi  use wfn_mod, only : domo, doci, postapsg, cionly, energy, nolagx
!nyi
!nyi  implicit none
!nyi  !--------------------------------------------------------------------
!nyi  integer, intent(in) :: imethod
!nyi  complex(kind(0d0)), intent(inout) :: wfn(*)
!nyi  !--------------------------------------------------------------------
!nyi  real(kind(0d0)), external :: field
!nyi  integer, external :: wfn_size
!nyi  real(kind(0d0)), external :: wfn_hprod, wfn_diag, wfn_occgbt
!nyi  complex(kind(0d0)), allocatable :: hwfn(:)
!nyi
!nyi  integer :: len, istep, jstep
!nyi  real(kind(0d0)) :: ene, ene1, ene_ci, ene_mov, ene_mo, lfield, test
!nyi  complex(kind(0d0)) :: dt, time
!nyi
!nyi  dt = - distep * iunit
!nyi  len = wfn_size(imethod)
!nyi  allocate(hwfn(len))
!nyi
!nyi  ene = zero
!nyi  do istep = 1, maxcyc
!nyi
!nyi     ! ci update
!nyi     domo = .false.
!nyi     doci = .true.
!nyi     ene_ci = wfn_diag(.true., zero, wfn, wfn, imethod)
!nyi     call wfn_ort(wfn, imethod)
!nyi     ene_ci = wfn_hprod(.true., zero, dt, wfn, wfn, hwfn, imethod)
!nyi
!nyi     ! mo update: occ-vir
!nyi     ene1 = ene_ci
!nyi     domo = .true.
!nyi     doci = .false.
!nyi     nolagx = .true.
!nyi     do jstep = 1, maxcyc
!nyi        time = dt * (jstep - 1)
!nyi        lfield = field(time)
!nyi
!nyi        ene_mov = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!nyi        call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
!nyi        call wfn_ort(wfn, imethod)
!nyi        ene_mov = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!nyi
!nyi        test = ene_mov - ene1
!nyi        ene1 = ene_mov
!nyi        write(iostdo, "('     # vir-micro-cycle:', i10, f20.10, e20.10)") jstep, ene_mov, test
!nyi        if(abs(test) < threne) then
!nyi           write(iostdo, "('     # vir-micro-cycle converged:', i10, f20.10)") jstep, ene_mov
!nyi           exit
!nyi        else if (test > zero) then
!nyi           write(iostdo, "('     # WARNING: energy raised.')")
!nyi        end if
!nyi     end do
!nyi
!nyi     ! mo update: occ-occ
!nyi     ene1 = ene_mov
!nyi     domo = .true.
!nyi     doci = .false.
!nyi     nolagx = .false.
!nyi     do jstep = 1, maxcyc
!nyi        time = dt * (jstep - 1)
!nyi        lfield = field(time)
!nyi
!nyi        ene_mo = wfn_occgbt(.true., lfield, wfn, wfn, imethod)
!nyi        call wfn_ort(wfn, imethod)
!nyi        ene_mo = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!nyi
!nyi        test = ene_mo - ene1
!nyi        ene1 = ene_mo
!nyi        write(iostdo, "('     # occ-micro-cycle:', i10, f20.10, e20.10)") jstep, ene_mo, test
!nyi        if(abs(test) < threne) then
!nyi           write(iostdo, "('     # occ-micro-cycle converged:', i10, f20.10)") jstep, ene_mo
!nyi           exit
!nyi        else if (test > zero) then
!nyi           write(iostdo, "('     # WARNING: energy raised.')")
!nyi        end if
!nyi     end do
!nyi
!nyi!     test = ene_mo - ene
!nyi     test = ene_mo - ene_mov
!nyi     ene = ene_mo
!nyi     write(iostdo, "('# cycle:', i10, 3f20.10, e20.10)") istep, ene_ci, ene_mov, ene_mo, test
!nyi     if(abs(test) < threne) then
!nyi        exit
!nyi     else if (test > zero) then
!nyi        write(iostdo, "('# WARNING: energy raised.')")
!nyi     end if
!nyi
!nyi  end do
!nyi
!nyi!  energy = ene + enen
!nyi!  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy
!nyi
!nyi  deallocate(hwfn)
!nyi
!nyiend subroutine init_cycle_diagpropgbt2
!nyi!################################################################################
!nyi!################################################################################
!nyisubroutine init_cycle_diagpropgbt3(wfn, imethod)
!nyi
!nyi  use mol_mod, only : enen
!nyi  use io_mod, only : iostdo
!nyi  use thresh_mod, only : threne
!nyi  use const_mod, only : iunit, zero, one
!nyi  use init_mod, only : maxcyc, distep, prop_type_init
!nyi  use wfn_mod, only : domo, doci, postapsg, cionly, energy, nolagx
!nyi
!nyi  implicit none
!nyi  !--------------------------------------------------------------------
!nyi  integer, intent(in) :: imethod
!nyi  complex(kind(0d0)), intent(inout) :: wfn(*)
!nyi  !--------------------------------------------------------------------
!nyi  real(kind(0d0)), external :: field
!nyi  integer, external :: wfn_size
!nyi  real(kind(0d0)), external :: wfn_hprod, wfn_diag, wfn_occgbt
!nyi  complex(kind(0d0)), allocatable :: hwfn(:)
!nyi
!nyi  integer :: len, istep, jstep
!nyi  real(kind(0d0)) :: ene, ene1, ene_ci, ene_vir, ene_occ, lfield, test
!nyi  complex(kind(0d0)) :: dt, time
!nyi
!nyi  dt = - distep * iunit
!nyi  len = wfn_size(imethod)
!nyi  allocate(hwfn(len))
!nyi
!nyi  ene = zero
!nyi  ene1 = one
!nyi  do istep = 1, maxcyc
!nyi
!nyi     do jstep = 1, maxcyc
!nyi        ! ci update
!nyi        domo = .false.
!nyi        doci = .true.
!nyi        ene_ci = wfn_diag(.true., zero, wfn, wfn, imethod)
!nyi        call wfn_ort(wfn, imethod)
!nyi        ene_ci = wfn_hprod(.true., zero, dt, wfn, wfn, hwfn, imethod)
!nyi
!nyi        ! mo update: occ-vir
!nyi        domo = .true.
!nyi        doci = .false.
!nyi        nolagx = .true.
!nyi        time = dt * (jstep - 1)
!nyi        lfield = field(time)
!nyi        ene_vir = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!nyi        call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
!nyi        call wfn_ort(wfn, imethod)
!nyi        ene_vir = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!nyi
!nyi        test = ene_vir - ene_ci
!nyi        write(iostdo, "('     # vir-micro-cycle:', i10, 2f20.10, e20.10)") jstep, ene_ci, ene_vir, test
!nyi        if(abs(test) < threne) then
!nyi           write(iostdo, "('     # vir-micro-cycle: converged', i20, f20.10)") jstep, ene_vir
!nyi           exit
!nyi        else if (test > zero) then
!nyi           write(iostdo, "('     # WARNING: energy raised.')")
!nyi        end if
!nyi     end do
!nyi
!nyi     do jstep = 1, maxcyc
!nyi        ! ci update
!nyi        domo = .false.
!nyi        doci = .true.
!nyi        ene_ci = wfn_diag(.true., zero, wfn, wfn, imethod)
!nyi        call wfn_ort(wfn, imethod)
!nyi        ene_ci = wfn_hprod(.true., zero, dt, wfn, wfn, hwfn, imethod)
!nyi
!nyi        ! mo update: occ-occ
!nyi        domo = .true.
!nyi        doci = .false.
!nyi        nolagx = .false.
!nyi        time = dt * (jstep - 1)
!nyi        lfield = field(time)
!nyi        ene_occ = wfn_occgbt(.true., lfield, wfn, wfn, imethod)
!nyi        call wfn_ort(wfn, imethod)
!nyi        ene_occ = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!nyi
!nyi        test = ene_occ - ene_ci
!nyi        write(iostdo, "('     # occ-micro-cycle:', i10, 2f20.10, e20.10)") jstep, ene_ci, ene_occ, test
!nyi        if(abs(test) < threne) then
!nyi           write(iostdo, "('     # occ-micro-cycle: converged', i20, f20.10)") jstep, ene_occ
!nyi           exit
!nyi        else if (test > zero) then
!nyi           write(iostdo, "('     # WARNING: energy raised.')")
!nyi        end if
!nyi     end do
!nyi
!nyi     ene = ene_occ
!nyi     test = ene - ene1
!nyi     ene1 = ene
!nyi     write(iostdo, "('# macro-cycle:', i10, 2f20.10, e20.10)") istep, ene_vir, ene_occ, test
!nyi     if(abs(test) < threne) then
!nyi        exit
!nyi     else if (test > zero) then
!nyi        write(iostdo, "('# WARNING: energy raised.')")
!nyi     end if
!nyi
!nyi  end do
!nyi
!nyi!  energy = ene + enen
!nyi!  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy
!nyi
!nyi  deallocate(hwfn)
!nyi
!nyiend subroutine init_cycle_diagpropgbt3
!nyi!################################################################################
