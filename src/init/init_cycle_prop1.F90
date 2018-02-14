!################################################################################
subroutine init_cycle_prop1(wfn, imethod)

  use mol_mod, only : enen
  use root_mod, only : name
  use io_mod, only : iostdo
  use const_mod, only : iunit, zero
  use thresh_mod, only : threne, thrgbt
  use init_mod, only : maxcyc, distep, prop_type_init
  use grid_mod, only : gbasis
  use wfn_mod, only : doci, domoq, domop1, domop2, domop, domo, &
       postapsg, energy, max_grad, noci, nomo, nolag1, nolag2, nolag3, &
       dpsi_test
  use cc_mod, only : fock,cc_solve_itr

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer,external :: wfn_poscic
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod
  complex(kind(0d0)), allocatable :: hwfn(:)

  integer :: len, istep
  real(kind(0d0)) :: ene, ene1, dtime, lfield, test
  complex(kind(0d0)) :: dt, time
!debug
  character(len = 16) :: cistep
  character(len = 16) :: fname
!debug

  if (prop_type_init == 'ABM') then
     call init_cycle_prop1_abm(wfn, imethod)
     return
  end if

  len = wfn_size(imethod)
  allocate(hwfn(len))

  dt = - distep * iunit
  ene1 = 1.d+10
  max_grad(1) = zero
  max_grad(2) = zero
  max_grad(3) = zero
  dpsi_test = threne * 1000

  ! gaussian case: make MO integrals 
  if (gbasis .and. nomo) then
     doci   = .true.
     domo   = .true.
     domop1 = .false.
     domop2 = .false.
     domoq  = .true.
     domop  = .false.
     ene = wfn_hprod(.true., lfield, dtime, wfn, wfn, hwfn, imethod)
  end if

  doci   = .not. noci
  domo   = .not. nomo
  domop1 = .not. nolag1
  domop2 = .not. nolag2
  domoq  = .not. nolag3
  domop  = domop1 .or. domop2

  do istep = 0, maxcyc

     dtime = abs(dt)
     time = dt * (istep - 1)
     lfield = field(time)
     ene = wfn_hprod(.true., lfield, dtime, wfn, wfn, hwfn, imethod)
     energy = ene + enen
!debug
!     call util_itochar(istep-1, 100, cistep)
!     fname = trim(trim(adjustl(name))//".orb.cyc"//trim(adjustl(cistep)))
!     call wfn_write_orb(fname, wfn, imethod)
!debug

     test = energy - ene1
     dpsi_test = test
     if (test < zero) then
        write(iostdo, "('# cycle ', i10, f20.10, 4e20.10)") &
             & istep, energy, test, max_grad(1:3)
     else
        write(iostdo, "('# cycle ', i10, f20.10, 4e20.10, ' <-- WARNING')") &
             & istep, energy, test, max_grad(1:3)
     end if

     if(abs(test) < threne .and. &
          & (nolag1 .or. max_grad(1) < thrgbt) .and. &
          & (nolag2 .or. max_grad(2) < thrgbt) .and. &
          & (nolag3 .or. max_grad(3) < thrgbt)) exit
     ene1 = energy

     !STOP "for debug @ init_cycle_prop1."

     if (trim(prop_type_init) == 'DIAG') then
        call wfn_fulldiag(wfn, imethod)
        call wfn_ort(wfn, imethod)
     else if (trim(prop_type_init) == 'U2PROP') then
        call wfn_u2prop(time, dt, wfn, hwfn, imethod)
        call wfn_ort(wfn, imethod)
     else
        call prop_new(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
        call wfn_ort(wfn, imethod)
     end if
     !if (istep == 100) stop 'for debug @ init_cycle_prop1.'

!DEBUG
!     stop 'for debug @ init_cycle_prop1'
!DEBUG
  end do

  if (postapsg) then
     call wfn_apsg(wfn, imethod)
     ene = wfn_hprod(.true., lfield, dtime, wfn, wfn, hwfn, imethod)
     energy = ene + enen
  end if

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy
  if (imethod == 0) call hf_print_eig(wfn)

  deallocate(hwfn)

end subroutine init_cycle_prop1
!################################################################################
!################################################################################
subroutine init_cycle_prop1_abm(wfn, imethod)

  use mol_mod, only : enen
  use root_mod, only : name
  use io_mod, only : iostdo
  use const_mod, only : iunit, zero
  use thresh_mod, only : threne, thrgbt
  use init_mod, only : maxcyc, distep, prop_type_init
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
  complex(kind(0d0)), allocatable :: hwfn(:,:)
  complex(kind(0d0)), allocatable :: twfn(:)

  integer :: len, istep
  real(kind(0d0)) :: ene, ene1, lfield, test
  complex(kind(0d0)) :: dt, time

  len = wfn_size(imethod)
  allocate(hwfn(1:len, 0:3))
  allocate(twfn(1:len))

  doci   = .not. noci
  domo   = .not. nomo
  domop1 = .not. nolag1
  domop2 = .not. nolag2
  domoq  = .not. nolag3
  domop  = domop1 .or. domop2

  dt = - distep * iunit
  ene1 = 1.d+10
  max_grad(1) = zero
  max_grad(2) = zero
  max_grad(3) = zero

  ! generate first three points by RK4 method
  do istep = 1, 3

     time = dt * (istep - 1)
     lfield = field(time)

     ene = wfn_hprod(.true., lfield, dt, wfn, wfn, twfn, imethod)
     hwfn(1:len, 3-istep+1) = twfn(1:len)

     test = ene - ene1
     if (test < zero) then
        write(iostdo, "('# cycle ', i10, f20.10, 4e20.10)") &
             & istep, ene, test, max_grad(1:3)
     else
        write(iostdo, "('# cycle ', i10, f20.10, 4e20.10, ' <-- WARNING')") &
             & istep, ene, test, max_grad(1:3)
     end if
     ene1 = ene

     call prop(time, dt, len, wfn, wfn, twfn, imethod, 'RK4')
     call wfn_ort(wfn, imethod)

  end do

  do istep = 4, maxcyc

     time = dt * (istep - 1)
     lfield = field(time)
     ene = wfn_hprod(.true., lfield, dt, wfn, wfn, twfn, imethod)
     hwfn(1:len, 0) = twfn(1:len)

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

     call prop_abm4(prop_type_init, time, dt, len, wfn, wfn, hwfn, imethod)
     call wfn_ort(wfn, imethod)

  end do

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy
  if (imethod == 0) call hf_print_eig(wfn)

  deallocate(twfn)
  deallocate(hwfn)

end subroutine init_cycle_prop1_abm
!################################################################################
