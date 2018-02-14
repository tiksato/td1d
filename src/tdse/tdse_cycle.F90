!################################################################################
subroutine tdse_cycle(wfn, imethodx)

  use fft_mod, only : dofft
  use prop_mod, only : prop_type

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethodx
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------

  if (dofft) call wfn_fft_plan(imethodx)

  if (trim(prop_type(1:3)) == 'VRK') then
     call tdse_cycle_vstep(wfn, imethodx)
  else if (trim(prop_type(1:3)) == 'ARK') then
     call tdse_cycle_astep(wfn, imethodx)
  else if (trim(prop_type(1:2)) == 'U1') then
     call tdse_cycle_u1prop(wfn, imethodx)
  else if (trim(prop_type(1:2)) == 'U2') then
     call tdse_cycle_u2prop(wfn, imethodx)
  else if (trim(prop_type(1:3)) == 'ARN') then
     call tdse_cycle_arnoldi(wfn, imethodx)
  else if (trim(prop_type(1:3)) == 'RK4') then
     call tdse_cycle_cstep(wfn, imethodx)
  else if (trim(prop_type(1:9)) == 'SPLIT_RK4') then
     call tdse_cycle_cstep(wfn, imethodx)
  else if (trim(prop_type(1:5)) == 'EULER') then
     call tdse_cycle_cstep(wfn, imethodx)
  else
     stop 'tdse_cycle: Bad prop_type.'
  end if

  if (dofft) call wfn_fft_destroy_plan(imethodx)

end subroutine tdse_cycle
!################################################################################
