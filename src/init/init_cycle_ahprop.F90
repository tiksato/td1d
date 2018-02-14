!################################################################################
subroutine init_cycle_ahprop(wfn, imethod)

  use mol_mod, only : enen
  use root_mod, only : name
  use io_mod, only : iostdo
  use const_mod, only : iunit, zero, one, czero
  use thresh_mod, only : thrwfn, threne, thrgrd, thrgbt
  use init_mod, only : maxcyc, maxdav, distepmo, prop_type_init, init_diffe, init_ah_trust
  use wfn_mod, only : cionly, doci, domoq, domop1, domop2, domop, domo, &
       & postapsg, energy, rms_grad, max_grad, noci, nomo, nolag1, nolag2, nolag3

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_aughess, wfn_hprod
  complex(kind(0d0)), allocatable :: pwfn(:)
  complex(kind(0d0)), allocatable :: hwfn(:)

  integer :: len, istep, icyc, max_prop
  complex(kind(0d0)) :: dt0, dt, time
  logical :: accept
  real(kind(0d0)) :: trad, ene_ci, ene, ene1, lfield, test
!debug
  character(len = 16) :: cistep
  character(len = 16) :: fname
!debug

  len = wfn_size(imethod)
  allocate(pwfn(len))
  allocate(hwfn(len))

  rms_grad(1) = zero
  rms_grad(2) = zero
  rms_grad(3) = zero
  max_grad(1) = zero
  max_grad(2) = zero
  max_grad(3) = zero

! max_prop = maxdav
  max_prop = 10

  time = czero
  lfield = field(time)

  dt0 = - distepmo * iunit
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

     ! augmented Hessian
     doci = .not. noci
     domo = .not. nomo
     domop1 = .not. nolag1
     domop2 = .not. nolag2
     domop = domop1 .or. domop2
     domoq = .false.
!debug     if (istep == 1) then
!debug        domo = .false.
!debug        domop1 = .false.
!debug        domop2 = .false.
!debug        domop = .false.
!debug        domoq = .false.
!debug     end if
     ene = wfn_aughess(lfield, accept, trad, wfn, imethod)
     ene_ci = ene
     call wfn_ort(wfn, imethod)

     ! mo update
     doci = .not. noci
!    doci = .false.

     domo = .not. nomo

     domop1 = .not. nolag1  ! for rough convergence
!    domop1 = .false.       ! for final convergence

     domop2 = .false.
!    domop2 = .not. nolag2

     domop = domop1 .or. domop2
     domoq  = .not. nolag3

     if (doci .or. domo) then
        dt = dt0
        pwfn(1:len) = wfn(1:len)
        do icyc = 1, max_prop
           ene = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
           call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
           call wfn_ort(wfn, imethod)
           ene = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)

           test = ene - ene_ci
           write(iostdo, "('# q-loop ', i6, f10.5, 2f20.12)") icyc, abs(dt), ene, test
           if (icyc == max_prop .or. ene < ene_ci + threne) then
              exit
           else
              dt = dt * 0.5d+0
              wfn(1:len) = pwfn(1:len)
           end if
        end do
     end if

     test = ene - ene1
     if (test > zero .and. test > threne * 10.d+0) then
        write(iostdo, "('# cycle ', i6, 3f20.12, 3f16.8, ' <-- WARNING')") &
             & istep, ene_ci, ene, test, rms_grad(1:2), max_grad(3)
     else
        write(iostdo, "('# cycle ', i6, 3f20.12, 3f16.8)") &
             & istep, ene_ci, ene, test, rms_grad(1:2), max_grad(3)
     end if

     if(abs(test) < thrwfn) exit
     if(abs(test) < threne .and. &
          rms_grad(1) < thrgrd .and. &
          & (nolag1 .or. rms_grad(2) < thrgrd) .and. &
          & (nolag2 .or. rms_grad(2) < thrgrd) .and. &
          & (nolag3 .or. max_grad(3) < thrgbt)) exit

     ene1 = ene
     init_diffe = test

  end do

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy

  deallocate(hwfn)
  deallocate(pwfn)

end subroutine init_cycle_ahprop
!################################################################################
