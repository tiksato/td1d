!################################################################################
subroutine guess_lcao(wfn, imethod)

  use mol_mod, only : natom, c
  use const_mod, only : two, four, pi
  use grid_mod, only : ngrid, x
  use wfn_mod, only : nfun, nspin
  use init_mod, only : broken_symmetry

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(out) :: wfn(0:ngrid, 1:nfun, 1:nspin, *)
  !--------------------------------------------------------------------
  integer :: igrid, ifun, ispin, iatom
  real(kind(0d0)) :: xn, x2, alpha, gauss

  call hf_clear(wfn)
  alpha = four * log(two)

  ! canonical rhf
  if (nspin == 1) then
     do iatom = 1, natom
        do igrid = 0, ngrid
           xn = x(igrid) - c(iatom)
           x2 = xn ** two
           do ifun = 1, nfun
              gauss = exp(-alpha * x2) * xn ** (ifun - 1)
              wfn(igrid, ifun, 1, 1) &
          & = wfn(igrid, ifun, 1, 1) + exp(-alpha * x2) * x(igrid) ** (ifun - 1)
           end do
        end do
     end do
  else
     if (broken_symmetry) then
! only for hydrogen-chain like systems, i.e., 
! number of electrons (sum of up- and down-spin orbitals) should be equal to the number of nuclei.
        do iatom = 1, natom
           do igrid = 0, ngrid
              xn = x(igrid) - c(iatom)
              x2 = xn ** two
              do ifun = 1, nfun
                 gauss = exp(-alpha * x2) * xn ** (ifun - 1)
                 wfn(igrid, ifun, iatom, 1) = exp(-alpha * x2) * x(igrid) ** (ifun - 1)
              end do
           end do
        end do
     else
        do iatom = 1, natom
           do igrid = 0, ngrid
              xn = x(igrid) - c(iatom)
              x2 = xn ** two
              do ifun = 1, nfun
                 gauss = exp(-alpha * x2) * xn ** (ifun - 1)
                 do ispin = 1, 2
                    wfn(igrid, ifun, ispin, 1) &
                & = wfn(igrid, ifun, ispin, 1) + exp(-alpha * x2) * x(igrid) ** (ifun - 1)
                 end do
              end do
           end do
        end do
     end if
  end if

  if (imethod >= 10) then
!!bug, why?     call util_zcopy((ngrid+1)*nfun*nspin, wfn(0, 1, 1, 1), 1, wfn(0, 1, 1, 2), 1)
!!     wfn(0:ngrid, 1:nfun, 1:nspin, 2) = wfn(0:ngrid, 1:nfun, 1:nspin, 1)
!     do ispin = 1, nspin
!        do ifun = 1, nfun
!           do igrid = 0, ngrid
!              wfn(igrid, ifun, ispin, 2) = wfn(igrid, ifun, ispin, 1)
!           end do
!        end do
!     end do
  end if
!debug
!write(6, "('init: guess wavefunction')")
!call wfn_print(6, wfn, imethod)
!debug

  call hf_ort(wfn)

 end subroutine guess_lcao
!######################################################################
