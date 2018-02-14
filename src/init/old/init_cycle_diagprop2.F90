!################################################################################
subroutine init_cycle_diagprop2(wfn, imethod)

  use mol_mod, only : enen
  use io_mod, only : iostdo
  use thresh_mod, only : threne
  use const_mod, only : iunit, zero
  use init_mod, only : maxcyc, distep, prop_type_init
  use wfn_mod, only : domo, doci, postapsg, cionly, energy

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod, wfn_diag
  complex(kind(0d0)), allocatable :: hwfn(:)

  integer :: len, istep
  real(kind(0d0)) :: ene, ene1, ene_ci, ene_mo, lfield, test
  complex(kind(0d0)) :: dt, time

  len = wfn_size(imethod)
  allocate(hwfn(len))

  dt = - distep * iunit
  ene1 = zero

  do istep = 1, maxcyc

     time = dt * (istep - 1)
     lfield = field(time)

     ! ci update
     domo = .false.
     doci = .true.
     ene_ci = wfn_diag(.true., lfield, wfn, wfn, imethod)
     call wfn_ort(wfn, imethod)
     ene_ci = wfn_hprod(.true., lfield, wfn, wfn, hwfn, imethod)
     if (cionly) then
        ene = ene_ci
        exit
     end if

!!     ! mo update
!!     domo = .true.
!!     doci = .false.
!!     ene_mo = wfn_hprod(.true., lfield, wfn, wfn, hwfn, imethod)
!!     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
!!     call wfn_ort(wfn, imethod)
!!     ene_mo = wfn_hprod(.true., lfield, wfn, wfn, hwfn, imethod)

     call init_cycle1(wfn, imethod, .true., .false., .true.)
     ene_mo = energy - enen

     ene = ene_mo
     test = ene - ene1
     write(iostdo, "('# cycle ', i10, 2f20.10, e20.10)") istep, ene_ci, ene_mo, test
     if(abs(test) < threne) then
        exit
     else if (test > zero) then
        write(iostdo, "('# WARNING: energy raised.')")
     end if
     ene1 = ene

  end do

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' macro-cycles:', 2f20.10)") istep, ene, energy

  deallocate(hwfn)

end subroutine init_cycle_diagprop2
!################################################################################
