!################################################################################
subroutine init_cycle_old(wfn, imethodx)

  use init_mod, only : init_type
  use wfn_mod, only : cionly
  use fft_mod, only : dofft

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethodx
  complex(kind(0d0)), intent(inout) :: wfn(*)
  !--------------------------------------------------------------------
  integer, external :: wfn_size

  if (dofft) call wfn_fft_plan(imethodx)

  if (imethodx <= 1 .or. imethodx == 4 .or. imethodx == 5) then
     ! tdse, hf, gvb, td-cis
     call init_cycle1(wfn, imethodx, .true., .true., .true., .true.)

  else
     ! cas, apsg, mrmp
     if (cionly) then
        if (trim(init_type) == 'DIAG_PROP') then
           call init_cycle_diagprop(wfn, imethodx)
        else
           call init_cycle1(wfn, imethodx, .true., .false., .false., .false.)
        end if
     else
        if (trim(init_type) == 'DIAG_PROP_GBT') then
           call init_cycle_diagpropgbt(wfn, imethodx)
!           call init_cycle_diagpropgbt2(wfn, imethodx)
!           call init_cycle_diagpropgbt3(wfn, imethodx)
        else if (trim(init_type) == 'PROP_GBT') then
           call init_cycle1(wfn, imethodx, .true., .false., .false., .false.)
           call init_cycle_propgbt(wfn, imethodx)
        else if (trim(init_type) == 'PROP1_GBT') then
           call init_cycle_prop1gbt(wfn, imethodx)
        else if (trim(init_type) == 'PROP2_GBT') then
           call init_cycle_prop2gbt(wfn, imethodx)
        else if (trim(init_type) == 'DIAG_PROP') then
           call init_cycle_diagprop(wfn, imethodx)
        else if (trim(init_type) == 'DIAG_PROP2') then
           call init_cycle_diagprop2(wfn, imethodx)
        else if (trim(init_type) == 'PROPO_PROPQ') then
           call init_cycle1(wfn, imethodx, .true., .true., .false.)
           call init_cycle1(wfn, imethodx, .true., .false., .true.)
        else if (trim(init_type) == 'PROP1') then
           call init_cycle1(wfn, imethodx, .false., .true., .false.)
           call init_cycle1(wfn, imethodx, .true., .true., .true.)
        else if (trim(init_type) == 'PROP1_NOCYC1') then
           call init_cycle1(wfn, imethodx, .true., .true., .true.)
        else if (trim(init_type) == 'PROP2') then
           call init_cycle1(wfn, imethodx, .false., .true., .false.)
           call init_cycle2(wfn, imethodx)
        else if (trim(init_type) == 'PROP2_NOCYC1') then
           call init_cycle2(wfn, imethodx)
        else if (trim(init_type) == 'PROP3') then
           call init_cycle1(wfn, imethodx, .false., .true., .false.)
           call init_cycle3(wfn, imethodx)
        else
           stop 'init_cycle: bad init_type.'
        end if
     end if

  end if

  if (dofft) call wfn_fft_destroy_plan(imethodx)

end subroutine init_cycle_old
!################################################################################
