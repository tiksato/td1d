!################################################################################
subroutine init_cycle_prop1gbt(wfn, imethod)

  use mol_mod, only : enen
  use io_mod, only : iostdo
  use const_mod, only : iunit, zero
  use thresh_mod, only : threne, thrgbt
  use init_mod, only : maxcyc, distep, prop_type_init
  use wfn_mod, only : domo, doci, energy, nolagx, max_grad

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod, wfn_diag, wfn_occgbt
  complex(kind(0d0)), allocatable :: hwfn(:)

  integer :: len, istep
  real(kind(0d0)) :: ene, ene1, ene_Q, ene_P, lfield, test
  complex(kind(0d0)) :: dt, time

  len = wfn_size(imethod)
  allocate(hwfn(len))

  dt = - distep * iunit
  ene1 = zero
  max_grad(1) = zero
  max_grad(2) = zero

  do istep = 1, maxcyc

     time = dt * (istep - 1)
     lfield = field(time)

     ! update cic's and occ-vir rotations
     domo = .true.
     doci = .true.
     nolagx = .true.
     ene_Q = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
     call wfn_ort(wfn, imethod)
     ene_Q = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)

     ! mo update: occ-occ
     domo = .true.
     doci = .false.
     nolagx = .false.
     ene_P = wfn_occgbt(.true., lfield, wfn, wfn, imethod)
     call wfn_ort(wfn, imethod)
     ene_P = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)

     ene = ene_P
     test = ene - ene1
     write(iostdo, "('# cycle ', i10, 2f17.10, 3e15.5)") &
          & istep, ene_Q, ene_P, test, max_grad(1:2)
     if (test > zero) write(iostdo, "('# WARNING: energy raised.')")

     if(abs(test) < threne .and. &
          & max_grad(1) < thrgbt .and. &
          & max_grad(2) < thrgbt) exit
     ene1 = ene

  end do

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") &
       & istep, ene, energy

  deallocate(hwfn)

end subroutine init_cycle_prop1gbt
!################################################################################
