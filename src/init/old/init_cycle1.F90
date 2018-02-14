!################################################################################
subroutine init_cycle1(wfn, imethod, docix, domoqx, domop1x, domop2x)

  use mol_mod, only : enen
  use io_mod, only : iostdo
  use const_mod, only : iunit, zero
  use thresh_mod, only : threne, thrgbt
  use init_mod, only : maxcyc, distep, prop_type_init
  use wfn_mod, only : doci, domoq, domop1, domop2, domop, domo, &
       & postapsg, energy, max_grad, nolag1, nolag2
!debug
!  use wfn_mod, only : nfun
!debug

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  logical, intent(in) :: docix, domoqx, domop1x, domop2x
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  real(kind(0d0)), external :: field
  integer, external :: wfn_size
  real(kind(0d0)), external :: wfn_hprod
  complex(kind(0d0)), allocatable :: hwfn(:)
  integer :: len, istep
  real(kind(0d0)) :: ene, ene1, lfield, test
  complex(kind(0d0)) :: dt, time

  len = wfn_size(imethod)
  allocate(hwfn(len))

  doci = docix
  domoq = domoqx
  domop1 = domop1x .and. .not. nolag1
  domop2 = domop2x .and. .not. nolag2
  domop = domop1 .or. domop2
  domo = domoq .or. domop
  dt = - distep * iunit
  ene1 = 1.d+10
  max_grad(1) = zero
  max_grad(2) = zero

  do istep = 1, maxcyc

     time = dt * (istep - 1)
     lfield = field(time)
     ene = wfn_hprod(.true., lfield, wfn, wfn, hwfn, imethod)

     test = ene - ene1
     if (test < zero) then
        write(iostdo, "('# cycle ', i10, f20.10, 3e20.10)") &
             & istep, ene, test, max_grad(1:2)
     else
        write(iostdo, "('# cycle ', i10, f20.10, 3e20.10, ' <-- WARNING')") &
             & istep, ene, test, max_grad(1:2)
     end if

     if(abs(test) < threne .and. &
          & (nolag1 .or. max_grad(1) < thrgbt) .and. &
          & (nolag2 .or. max_grad(2) < thrgbt)) exit
     ene1 = ene

     if (trim(prop_type_init) /= 'DIAG') then
        call prop(time, dt, len, wfn, wfn, hwfn, imethod, prop_type_init)
!debug
!  write(6, "('init_cycle1: force g- or u-symmetry!')")
!  call general_symm(nfun, wfn)
!debug
        call wfn_ort(wfn, imethod)
     else
        call wfn_fulldiag(wfn, imethod)
        call wfn_ort(wfn, imethod)
     end if
!debug
!     if (istep == 1) then
!        write(6, "('stop for debug in init_cycle1', f20.10)") ene
!        stop
!     end if
!debug

  end do

  if (postapsg) then
     call wfn_apsg(wfn, imethod)
     ene = wfn_hprod(.true., lfield, wfn, wfn, hwfn, imethod)
  end if

  energy = ene + enen
  write(iostdo, "('# total energy after ', i10, ' cycles:', 2f20.10)") istep, ene, energy
  if (imethod == 0) call hf_print_eig(wfn)

  deallocate(hwfn)

end subroutine init_cycle1
!################################################################################
