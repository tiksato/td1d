!######################################################################
subroutine init_cycle_diagdiag(wfn, imethod)

  use io_mod, only : iostdo
  use thresh_mod, only : threne
  use const_mod, only : zero, iunit
  use init_mod, only : maxcyc, distep, prop_type_init
  use wfn_mod, only : domo, doci, cionly, energy

  implicit none
  complex(kind(0d0)), intent(inout) :: wfn(*)
  integer, intent(in) :: imethod

  integer :: istep, len
  real(kind(0d0)) :: ene, ene1, ene_ci, ene_mo, test
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_diag, wfn_diagmo, wfn_hprod
  complex(kind(0d0)) :: dt, time
  complex(kind(0d0)), allocatable :: hwfn(:)

  len = wfn_size(imethod)
  allocate(hwfn(len))

  dt = - distep * iunit
  ene1 = 1.d+10

  domo = .false.
  doci = .true.
  ene = wfn_diag(.true., zero, wfn, wfn, imethod)
!no_orth  call wfn_ort(wfn, imethod)
  if (cionly) then
     energy = ene
     write(iostdo, "('# total energy after ', i10, ' cycles:', f20.10)") 1, energy
     return
  else
     write(iostdo, "('# init  ', 10x, f20.10)") ene
  end if

  ene1 = ene

  do istep = 1, maxcyc

     time = dt * (istep - 1)

     ! mo update
     domo = .true.
     doci = .false.
     ene_ci = wfn_diagmo(.true., zero, wfn, wfn, imethod)
!no_orth     call wfn_ort(wfn, imethod)

     ! ci update
     domo = .false.
     doci = .true.
!debug
     ene_mo = wfn_diag(.true., zero, wfn, wfn, imethod)
!     ene_mo = wfn_hprod(.true., zero, wfn, wfn, hwfn, imethod)
!     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
!debug
     call wfn_ort(wfn, imethod)

     ene = ene_mo

     test = ene - ene1
     write(iostdo, "('# cycle ', i10, 2f20.10, e20.10)") istep, ene_ci, ene_mo, test
     if(abs(test) < threne) exit

     if (test > zero) then
        write(iostdo, "('# WARNING: energy raised.')")
     end if

     ene1 = ene

  end do

  energy = ene
  write(iostdo, "('# total energy after ', i10, ' cycles:', f20.10)") istep, energy

  deallocate(hwfn)

end subroutine init_cycle_diagdiag
!######################################################################
