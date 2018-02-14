!######################################################################
subroutine tdse(wfn, imethodx)

  use root_mod, only : name, ical, icomp, iprint
  use const_mod, only : iostdo, czero
  use prop_mod, only : prop_type
  use wfn_mod, only : dopt

  implicit none
  integer, intent(in) :: imethodx
  complex(kind(0d0)), intent(inout) :: wfn(*)

  integer, external :: wfn_size
  real(kind(0d0)), external :: util_clock
  complex(kind(0d0)), allocatable :: wfn0(:)
  character(len = 256) :: fname
  integer :: len
  integer :: itime
  real(kind(0d0)) :: rtime
!debug
  integer :: poscic, pospt1
  real(kind(0d0)) :: ene
  integer, external :: cas_poscic, cas_pospt1
  real(kind(0d0)), external :: mrmp_ene_mp2, mrmp_ene_pt1
!debug

  if (ical < 0) return
  icomp = 1

  rtime = util_clock(.true., itime)

!debug
  if (dopt) then
     poscic = cas_poscic()
     pospt1 = cas_pospt1()
!clear wfnpt
     len = wfn_size(-1)
     call util_zcopy(len, czero, 0, wfn(pospt1), 1)
!clear wfnpt
     ene = mrmp_ene_mp2(wfn, wfn(poscic), wfn(pospt1))
     write(iostdo, "('# second order energy:                 ', f20.10)") ene
     ene = mrmp_ene_pt1(wfn, wfn(poscic), wfn(pospt1))
     write(iostdo, "('# energy of 1st order wavefunction:    ', f20.10)") ene
  end if
!debug

  len = wfn_size(imethodx)
  allocate(wfn0(1:len))
  wfn0(1:len) = wfn(1:len)

  if (trim(prop_type) == 'VRK5') then
     call tdse_cycle_vstep(wfn0, wfn, imethodx, prop_type)
  else
     call tdse_cycle(wfn0, wfn, imethodx, prop_type)
  end if

  fname = trim(name)//".wfn1"
  call wfn_write(fname, wfn, imethodx)
  if (dopt) then
     fname = trim(name)//".wfnpt1"
     call wfn_writept(fname, wfn, imethodx)
  end if

  if (iprint > 2) then
     write(iostdo, "('# final wavefunction')")
     call wfn_print(iostdo, wfn, imethodx)
  end if

  deallocate(wfn0)

  rtime = util_clock(.false., itime)
  write(6, "('# wall clock time for tdse:   ', E12.5)") rtime

end subroutine tdse
!######################################################################
