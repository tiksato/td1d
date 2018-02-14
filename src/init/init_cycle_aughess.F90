!################################################################################
subroutine init_cycle_aughess(wfn, imethod)

  use mol_mod, only : enen
  use root_mod, only : name
  use io_mod, only : iostdo
  use const_mod, only : iunit, zero, one, czero
  use thresh_mod, only : thrwfn, threne, thrgrd
  use init_mod, only : maxcyc, distepmo, prop_type_init, init_diffe, init_ah_trust
  use wfn_mod, only : cionly, doci, domoq, domop1, domop2, domop, domo, &
       & postapsg, energy, rms_grad, noci, nomo, nolag1, nolag2, nolag3

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_aughess, wfn_hprod
  complex(kind(0d0)), allocatable :: hwfn(:)

  integer :: len, istep
  complex(kind(0d0)) :: dt, time
  logical :: accept
  real(kind(0d0)) :: trad, ene, ene1, lfield, test
!debug
  character(len = 16) :: cistep
  character(len = 16) :: fname
!debug

  len = wfn_size(imethod)
  allocate(hwfn(len))

  time = czero
  lfield = field(time)

  dt = - distepmo * iunit
  init_diffe = -one

  ! initial trust radius
  trad = init_ah_trust

  doci = .not. noci
  domo = .not. nomo
  domop1 = .not. nolag1
  domop2 = .not. nolag2
  domop = domop1 .or. domop2
  domoq = .not. nolag3
!old  doci = .not. noci
!old  domo = .false.
!old  domop1 = .false.
!old  domop2 = .false.
!old  domop = .false.
!old  domoq = .false.
  ene1 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!debug  ene1 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!debug  ene1 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
!debug  ene1 = wfn_aughess(lfield, accept, trad, wfn, imethod)
  write(iostdo, "('# cycle ', i6, f20.12)") 0, ene1       

  do istep = 1, maxcyc

!debug
     call util_itochar(istep-1, 100, cistep)
     fname = trim(trim(adjustl(name))//".orb.cyc"//trim(adjustl(cistep)))
     call wfn_write_orb(fname, wfn, imethod)
!debug
     ! augmented Hessian
     doci = .not. noci
     domo = .not. nomo
     domop1 = .not. nolag1
     domop2 = .not. nolag2
     domop = domop1 .or. domop2
     domoq = .not. nolag3
!debug     if (istep == 1) then
!debug        domo = .false.
!debug        domop1 = .false.
!debug        domop2 = .false.
!debug        domop = .false.
!debug        domoq = .false.
!debug     end if
     ene = wfn_aughess(lfield, accept, trad, wfn, imethod)
     call wfn_ort(wfn, imethod)

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

     if (.not. accept) then
        write(iostdo, "('# cycle ', i6, 2f20.12, 3f16.8, ' <-- Rejected')") &
             & istep, ene, test, rms_grad(1:3)
     else if (test > zero .and. test > threne * 10.d+0) then
        write(iostdo, "('# cycle ', i6, 2f20.12, 3f16.8, ' <-- WARNING')") &
             & istep, ene, test, rms_grad(1:3)
     else
        write(iostdo, "('# cycle ', i6, 2f20.12, 3f16.8)") &
             & istep, ene, test, rms_grad(1:3)
     end if

     if(abs(test) < thrwfn) exit
     if(abs(test) < threne .and. &
          rms_grad(1) < thrgrd .and. &
          & (nolag1 .or. rms_grad(2) < thrgrd) .and. &
          & (nolag2 .or. rms_grad(2) < thrgrd) .and. &
          & (nolag3 .or. rms_grad(3) < thrgrd)) exit

     ene1 = ene
     init_diffe = test

  end do

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy

  deallocate(hwfn)

end subroutine init_cycle_aughess
!################################################################################
