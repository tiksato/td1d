!######################################################################
subroutine tdse_print(istep, time, lfield, ene, q0, q1, q2, wfn0, wfn, nprint, imethodx)

  use root_mod, only : name, nprintwfn, printwfn
  use const_mod, only : one, four
  use grid_mod, only : rmax, x0, xmask
  use field_mod, only : period
  use prop_mod, only : totstep

  implicit none
  integer, intent(in) :: istep
  complex(kind(0d0)), intent(in) :: time
  real(kind(0d0)), intent(inout) :: lfield, ene, q0, q1, q2
  complex(kind(0d0)), intent(in) :: wfn0(*)
  complex(kind(0d0)), intent(in) :: wfn(*)
  integer, intent(inout) :: nprint
  integer, intent(in) :: imethodx

  character(len = 16) :: istep_char = ''
  character(len = 256) :: fname
  real(kind(0d0)), external :: wfn_op1e
  real(kind(0d0)), external :: wfn_norm
  complex(kind(0d0)), external :: wfn_iprod
  real(kind(0d0)) :: optcyc, p0, p1, p2, ptot, dip, acc
  complex(kind(0d0)) :: c0

  optcyc = real(time) / period
  c0 = wfn_iprod(x0, wfn0, wfn, imethodx)
  dip = wfn_op1e(0, -one, wfn, imethodx)
  acc = wfn_op1e(1, -one, wfn, imethodx)
  ptot = wfn_norm(rmax, wfn, p0, p1, p2, imethodx)

!debug dip = wfn_op1e(0, rmax, wfn, imethodx)
!debug acc = wfn_op1e(1, rmax, wfn, imethodx)
!debug call wfn_print_cic(6, wfn, imethodx)

  write(6,"(' step ', i10, 14E15.6)") &
&   istep, real(time), optcyc, lfield, c0, p0, p1, p2, q0, q1, q2, dip, acc, ene

!no write  if (istep == 0 .or. &
!no write    & istep == totstep .or. &
!no write    & printwfn(nprint) == istep .or. &
!no write    & nprintwfn < 0 .and. mod(istep, abs(nprintwfn)) == 0) then
!no write
!no write     nprint = nprint + 1
!no write
!no write     call util_itochar(istep, totstep, istep_char)
!no write     fname = trim(name)//".Step"//trim(istep_char)//".wfn"
!no write     call wfn_write(fname, wfn, imethodx)
!no write  end if

end subroutine tdse_print
!######################################################################
