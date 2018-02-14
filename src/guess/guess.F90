!################################################################################
subroutine guess

  use io_mod, only : iostdo
  use grid_mod, only : gbasis
  use wfn_mod, only : wfn, imethod
  use root_mod, only : name, iprint
  use init_mod, only : guess_type, guess_alter

  implicit none
  real(kind(0d0)), external :: util_clock
  character(len = 256) :: fname
  integer :: itime
  real(kind(0d0)) :: rtime

  rtime = util_clock(.true., itime)

  if (trim(guess_type) == 'NONE') then
     call wfn_clear(wfn, imethod)
  else if (trim(guess_type) == 'READ') then
     call guess_read(wfn, imethod)
  else if (trim(guess_type) == 'READG') then
     if (gbasis) then
        call guess_read_gbasis(wfn, imethod)
     else
        stop 'READG not supported for grid simulation.'
     end if
  else if (trim(guess_type) == 'READOPT') then
     call guess_read(wfn, imethod)
  else if (trim(guess_type) == 'HF') then
     call guess_hf(wfn, imethod)
  else if (trim(guess_type) == 'HFG') then
     if (gbasis) then
        call guess_read_gbasis(wfn, imethod)
        call init_cycle_fulld(wfn,0)
     else
        stop 'HFG not supported for grid simulation.'
     end if
  else if (trim(guess_type) == 'NO') then
     call guess_no(wfn, imethod)
  else if (trim(guess_type) == 'SUBNO') then
     call guess_no(wfn, imethod)
  else if (trim(guess_type) == 'CAS') then
     call guess_cas(wfn, imethod)
  else if (trim(guess_type) == 'APSG') then
     call guess_apsg(wfn, imethod)
  else if (trim(guess_type) == 'X2E') then
     call guess_x2e(wfn, imethod)
  else if (trim(guess_type) == 'LCAO') then
     if (gbasis) then
        call guess_core(wfn, imethod)
     else
        call guess_lcao(wfn, imethod)
     end if
!  else if (trim(guess_type) == 'READ_LCAO') then
!     call guess_lcao(wfn, imethod)
!     call guess_read(wfn, imethod)
  else
     stop 'guess: bad guess_type.'
  end if

  if (guess_alter) call guess_order(wfn)

  fname = trim(name)//".wfng"
  call wfn_write(fname, wfn, imethod)

  if (iprint > 2) then
     write(iostdo, "('# guess wavefunction')")
     call wfn_print(iostdo, wfn, imethod)
  end if

  rtime = util_clock(.false., itime)
  write(6, "('# wall clock time for guess:   ', E12.5)") rtime

!debug
!stop 'for debug in guess.'
!debug

!nyi  call initpt(wfn, imethod)

end subroutine guess
!################################################################################
