!################################################################################
subroutine init_cycle2(wfn, imethod)

  use mol_mod, only : enen
  use io_mod, only : iostdo
  use const_mod, only : iunit, zero
  use thresh_mod, only : threne, thrgbt
  use init_mod, only : maxcyc, distep, prop_type_init
  use wfn_mod, only : domo, doci, postapsg, energy, max_grad, nolag1, nolag2

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod
  real(kind(0d0)), external :: wfn_occgbt
  complex(kind(0d0)), allocatable :: pwfn(:)
  complex(kind(0d0)), allocatable :: hwfn(:)

  integer :: len, istep
  real(kind(0d0)) :: ene, ene1, ene_ci, ene_mo, lfield, test, test1, test2
  complex(kind(0d0)) :: dt, time

  len = wfn_size(imethod)
  allocate(pwfn(len))
  allocate(hwfn(len))

  dt = - distep * iunit
  ene1 = 1.d+10
  max_grad(1) = zero
  max_grad(2) = zero
  pwfn(1:len) = wfn(1:len)

  do istep = 1, maxcyc

     time = dt * (istep - 1)
     lfield = field(time)

     ! ci update
     domo = .false.
     doci = .true.
     ene_mo = wfn_hprod(.true., lfield, wfn, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
     call wfn_ort(wfn, imethod)

     ! mo update
     domo = .true.
     doci = .false.
     ene_ci = wfn_hprod(.true., lfield, wfn, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
     call wfn_ort(wfn, imethod)

     ene = ene_ci

     test1 = ene - ene1
     if (test1 > zero) then
        write(iostdo, "('# WARNING: energy raised (ene > ene_1).')")
     end if
     test2 = ene_ci - ene_mo
     if (test2 > zero) then
        write(iostdo, "('# WARNING: energy raised (ene_mo > ene_ci).')")
     end if

!strict     if (abs(test1) > abs(test2)) then
!strict        test = test1
!strict     else
!strict        test = test2
!strict     end if
     test = test1
     write(iostdo, "('# cycle ', i10, 2f20.10, 3e20.10)") &
          & istep, ene_mo, ene_ci, test, max_grad(1:2)

     if(abs(test) < threne .and. &
          & (nolag1 .or. max_grad(1) < thrgbt) .and. &
          & (nolag2 .or. max_grad(2) < thrgbt)) exit

     ene1 = ene
     pwfn(1:len) = wfn(1:len)

  end do

  if (postapsg) then
     call wfn_apsg(wfn, imethod)
     ene = wfn_hprod(.true., lfield, wfn, wfn, hwfn, imethod)
  end if

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy

  deallocate(hwfn)
  deallocate(pwfn)

end subroutine init_cycle2
!################################################################################
