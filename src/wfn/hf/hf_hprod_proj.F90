!######################################################################
subroutine hf_hprod_proj(dtime, wfn, hwfn)

  use root_mod, only : icomp
  use mol_mod, only : ne
  use grid_mod, only : ngrid
  use wfn_mod, only : nfcore, nvir, nfun, nspin, hf_doproj
  use const_mod, only : one, czero, runit, iunit

  implicit none
  real(kind(0d0)), intent(in) :: dtime
  complex(kind(0d0)), intent(in) :: wfn(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(inout) :: hwfn(0:ngrid, 1:nfun, 1:nspin)

  integer :: ispin, ifun, jfun, nbas, noccx
  complex(kind(0d0)) :: fac
  complex(kind(0d0)), allocatable :: wfn_dot_hwfn(:,:)
  

  nbas = ngrid + 1

  if (hf_doproj) then
     allocate(wfn_dot_hwfn(1:nfun, 1:nfun))
     do ispin = 1, nspin
   !bug noccx = nocc
        noccx = ne(ispin)
   
        call hf_ovlp(.true., -one, wfn(0, 1, ispin), hwfn(0, 1, ispin), wfn_dot_hwfn)
   
        ! occupied orbitals
        do ifun = nfcore + 1, noccx
           do jfun = 1, noccx
              fac = - wfn_dot_hwfn(jfun, ifun)
              call zaxpy(nbas, fac, wfn(0, jfun, ispin), 1, hwfn(0, ifun, ispin), 1)
           end do
        end do
   
        ! virtual orbitals
        do ifun = noccx + 1, nfun
           do jfun = 1, noccx
              fac = conjg(wfn_dot_hwfn(ifun, jfun))
              call util_zcopy(nbas, czero, 0, hwfn(0, ifun, ispin), 1)
              call zaxpy(nbas, fac, wfn(0, jfun, ispin), 1, hwfn(0, ifun, ispin), 1)
           end do
        end do
     end do
     deallocate(wfn_dot_hwfn)
  end if

  ! multipy time step
  if (icomp == 1) then
     fac = -iunit * dtime
  else
     fac = -runit * dtime
  end if
  call util_zscal(nbas*nfun*nspin, fac, hwfn, 1)

end subroutine hf_hprod_proj
!######################################################################
