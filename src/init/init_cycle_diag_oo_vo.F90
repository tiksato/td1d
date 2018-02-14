!################################################################################
subroutine init_cycle_diag_oo_vo(wfn, imethod)

  use mol_mod, only : enen
  use io_mod, only : iostdo
  use const_mod, only : iunit, zero, czero
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
  complex(kind(0d0)), allocatable :: pwfn(:)

  integer :: len, istep
  real(kind(0d0)) :: ene, ene_prev, ene1, ene2, ene3, lfield, test, test1, test2, test3
  complex(kind(0d0)) :: dt, dtp, dtq, time

  len = wfn_size(imethod)
  allocate(hwfn(len))
  allocate(pwfn(len))
  call util_zcopy(len, wfn, 1, pwfn, 1)

  time = czero
  lfield = field(time)

  dtp = - min(distep1, distep2) * iunit
  dtq = - distep3 * iunit
  ene_prev = zero

  max_grad(1) = zero
  max_grad(2) = zero
  max_grad(3) = zero

  do istep = 1, maxcyc

     ! ci update
     doci = .not. noci
     domo = .false.
     domop1 = .false.
     domop2 = .false.
     domop = .false.
     domoq = .false.
     ene1 = wfn_diag(.true., lfield, wfn, wfn, imethod)
     test1 = ene1 - ene_prev
     call wfn_ort(wfn, imethod)
!     if (test1 > threne) then
!        call util_zcopy(len, pwfn, 1, wfn, 1)
!        ene1 = ene_prev
!        test1 = zero
!     else
        call util_zcopy(len, wfn, 1, pwfn, 1)
!     end if
     if (cionly) then
        ene = ene1
        exit
     end if

     ! o-o rotations update
     dt = dtp
     doci = .false.
     domo = .true.
     domop1 = .not. nolag1
     domop2 = .not. nolag2
     domop  = domop1 .or. domop2
     domoq  = .false.
     ene1 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
     ene2 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     test2 = ene2 - ene1
!     if (test2 > threne) then
!        call util_zcopy(len, pwfn, 1, wfn, 1)
!        ene2 = ene1
!        test2 = zero
!     else
        call util_zcopy(len, wfn, 1, pwfn, 1)
!     end if


     ! v-o rotations update
     dt = dtq
     doci = .false.
     domo = .true.
     domop1 = .false.
     domop2 = .false.
     domop  = .false.
     domoq  = .not. nolag3
     ene2 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
     call wfn_ort(wfn, imethod)
     ene3 = wfn_hprod(.true., lfield, dt, wfn, wfn, hwfn, imethod)
     test3 = ene3 - ene2
!     if (test3 > threne) then
!        call util_zcopy(len, pwfn, 1, wfn, 1)
!        ene3 = ene2
!        test3 = zero
!     else
        call util_zcopy(len, wfn, 1, pwfn, 1)
!     end if

     ene = ene3
     test = ene - ene_prev

     write(iostdo, "('# cycle ', i10, 3f20.10, 4e20.10)", advance = 'no') &
          & istep, ene1, ene2, ene3, test, max_grad(1:3)
     if (test1 < threne .and. test2 < threne .and. test3 < threne) then
        write(iostdo, *)
     else 
        write(iostdo, "(' <-- ')", advance = 'no')
        if (test1 > threne) write(iostdo, "('CI ')", advance = 'no')
        if (test2 > threne) write(iostdo, "('OO ')", advance = 'no')
        if (test3 > threne) write(iostdo, "('VO ')", advance = 'no')
        write(iostdo, "('raise energy.')")
     end if

     if(abs(test) < threne .and. &
          & (nolag1 .or. max_grad(1) < thrgbt) .and. &
          & (nolag2 .or. max_grad(2) < thrgbt) .and. &
          & (nolag3 .or. max_grad(3) < thrgbt)) exit

     ene_prev = ene

  end do

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy

  deallocate(pwfn)
  deallocate(hwfn)

end subroutine init_cycle_diag_oo_vo
!################################################################################
