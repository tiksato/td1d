!######################################################################
real(kind(0d0)) function hf_hprod_sum(calene, wfn, h1wfn, h2wfn, hwfn, ng0, ng1)

  use mol_mod, only : ne
  use root_mod, only : icomp
  use const_mod, only : zero
  use grid_mod, only : ngrid
  use wfn_mod, only : nfun, nspin

  implicit none
  logical, intent(in) :: calene
  complex(kind(0d0)), intent(in) :: wfn  (0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in) :: h1wfn(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(in) :: h2wfn(0:ngrid, 1:nfun, 1:nspin)
  complex(kind(0d0)), intent(out) :: hwfn(0:ngrid, 1:nfun, 1:nspin)
  integer, intent(in) :: ng0, ng1

  real(kind(0d0)) :: ene1, ene2
  real(kind(0d0)), external :: hf_enec, hf_ene2

  hf_hprod_sum = zero

  if(icomp < 0 .or. ne(3) == 1) then

     hwfn(ng0:ng1, 1:nfun, 1:nspin) = h1wfn(ng0:ng1, 1:nfun, 1:nspin)
     if (calene) hf_hprod_sum = hf_enec(wfn, h1wfn, ng0, ng1)
  else

     hwfn(ng0:ng1, 1:nfun, 1:nspin) = h1wfn(ng0:ng1, 1:nfun, 1:nspin) &
                                  & + h2wfn(ng0:ng1, 1:nfun, 1:nspin)
     if (calene) then
        ene1 = hf_enec(wfn, h1wfn, ng0, ng1)
        ene2 = hf_ene2(wfn, h2wfn, ng0, ng1)
        hf_hprod_sum = ene1 + ene2
!debug
!write(6, "(' ene1 = ', f20.10)") ene1
!write(6, "(' ene2 = ', f20.10)") ene2
!debug
     end if
  end if

  return

end function hf_hprod_sum
!######################################################################
