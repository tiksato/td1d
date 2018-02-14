!################################################################################
subroutine init_cycle(wfn, imethodx)

  use fft_mod, only : dofft
  use cc_mod, only : cc_solve,fock
  use wfn_mod, only : tdcc,semicanonical
  use init_mod, only : init_type, maxcyc

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethodx
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  integer,external :: wfn_poscic

  if (maxcyc <= 0) return
  if (dofft) call wfn_fft_plan(imethodx)

  if (trim(init_type) == 'FULLDIAG') then
     call init_cycle_fulld(wfn, imethodx)
  else if (trim(init_type) == 'AHPROP') then
     call init_cycle_ahprop(wfn, imethodx)
  else if (trim(init_type) == 'AUGHESS') then
     call init_cycle_aughess(wfn, imethodx)
  else if (trim(init_type) == 'AUGHESS_PROP') then
     call init_cycle_aughess_prop(wfn, imethodx)
  else if (trim(init_type) == 'DIAG') then
     call init_cycle_diag(wfn, imethodx)
  else if (trim(init_type) == 'PROP1') then
     call init_cycle_prop1(wfn, imethodx)
  else if (trim(init_type) == 'PROP2') then
     call init_cycle_prop2(wfn, imethodx)
  else if (trim(init_type) == 'PROP3') then
     call init_cycle_prop3(wfn, imethodx)
  else if (trim(init_type) == 'PROP_CIP_Q') then
     call init_cycle_propcipq(wfn, imethodx)
  else if (trim(init_type) == 'DIAGP_PROPQ') then
     call init_cycle_diagp_propq(wfn, imethodx)
  else if (trim(init_type) == 'DIAG_PROP') then
     call init_cycle_diagprop(wfn, imethodx)
  else if (trim(init_type) == 'DIAG_OO_VO') then
     call init_cycle_diag_oo_vo(wfn, imethodx)
  else if (trim(init_type) == 'DIAG_CIP_Q') then
     call init_cycle_diagcipq(wfn, imethodx)
  else if (trim(init_type) == 'PROP_CIP2_P1Q') then
     call init_cycle_propcip2p1q(wfn, imethodx)
  else if (trim(init_type) == 'DIAG_CIP2_P1Q') then
     call init_cycle_diagcip2p1q(wfn, imethodx)
  else if (trim(init_type) == 'PROP1_GBT') then
     call init_cycle_prop1gbt(wfn, imethodx)
  else if (trim(init_type) == 'DIAG_PROP_GBT') then
!nyi     call init_cycle_diagpropgbt(wfn, imethodx)
  else if (trim(init_type) == 'PROP_GBT') then
!nyi     call init_cycle_proppropgbt(wfn, imethodx)
  else if (trim(init_type) == 'PROP_PROP') then
!nyi     call init_cycle_proppropprop(wfn, imethodx)
  else
     stop 'init_cycle: bad init_type.'
  end if

  if (dofft) call wfn_fft_destroy_plan(imethodx)

end subroutine init_cycle
!################################################################################
