!######################################################################
subroutine hf_print_op1e(iout, istep, time, rmax, wfn)

  use omp_mod
  use mol_mod, only : ne
  use root_mod, only : iprint
  use field_mod, only : period
  use const_mod, only : zero, one, two, czero
  use grid_mod, only : ngrid, dgrid, x, gv1, gll, gul
  use wfn_mod, only : nfun, nspin

  implicit none
  integer, intent(in) :: iout, istep
  complex(kind(0d0)), intent(in) :: time
  real(kind(0d0)), intent(in) :: rmax
  complex(kind(0d0)), intent(in)  :: wfn(0:ngrid, 1:nfun, 1:nspin)

  real(kind(0d0)), allocatable  :: oppno(:,:,:,:)
  integer igrid, ifun, ispin, llgrid, ulgrid
  real(kind(0d0)) :: hffac, rho, optcyc, dip, acc, occd

  allocate(oppno(1:nfun, 1:2, 1:2, 0:(nproc-1)))

  if (nspin == 2 .or. ne(3) == 1) then
     hffac = one
  else
     hffac = two
  end if

  dip = zero
  acc = zero
  oppno(1:nfun, 1:2, 1:2, 0:(nproc-1)) = zero
  call get_irmax(rmax, llgrid, ulgrid, gll, gul)

!$omp parallel default(shared) private(ifun,igrid,ispin,rho)

  call omp_mod_thread(llgrid, ulgrid)

  do ispin = 1, nspin
     do ifun = 1, ne(ispin)
        do igrid = ng0, ng1
           rho = dble(conjg(wfn(igrid, ifun, ispin)) * wfn(igrid, ifun, ispin))
           oppno(ifun, ispin, 1, iproc) = oppno(ifun, ispin, 1, iproc) + rho * x(igrid)
           oppno(ifun, ispin, 2, iproc) = oppno(ifun, ispin, 2, iproc) + rho * gv1(igrid)
        end do
     end do
  end do

!$omp end parallel

  do iproc = 1, nproc - 1
     do ispin = 1, nspin
        do ifun = 1, ne(ispin)
           oppno(ifun, ispin, 1, 0) = oppno(ifun, ispin, 1, 0) + oppno(ifun, ispin, 1, iproc)
           oppno(ifun, ispin, 2, 0) = oppno(ifun, ispin, 2, 0) + oppno(ifun, ispin, 2, iproc)
        end do
     end do
  end do

  do ispin = 1, nspin
     do ifun = 1, ne(ispin)
        occd = hffac * dgrid
        oppno(ifun, ispin, 1, 0) = oppno(ifun, ispin, 1, 0) * occd
        oppno(ifun, ispin, 2, 0) = oppno(ifun, ispin, 2, 0) * occd
        dip = dip + oppno(ifun, ispin, 1, 0)
        acc = acc + oppno(ifun, ispin, 2, 0)
     end do
  end do

  optcyc = dble(time) / period
  write(iout, "(' op1e: ', i10, 3E15.6)") istep, optcyc, dip, acc

  if (iprint > 1) then
     write(iout, "(' dip-no: ', i10, E15.6)", advance = 'no') istep, optcyc
     do ispin = 1, nspin
        do ifun = 1, ne(ispin)
           write(iout, "(E15.6)", advance = 'no') oppno(ifun, ispin, 1, 0)
        end do
     end do
     write(iout, *)

     write(iout, "(' acc-no: ', i10, E15.6)", advance = 'no') istep, optcyc
     do ispin = 1, nspin
        do ifun = 1, ne(ispin)
           write(iout, "(E15.6)", advance = 'no') oppno(ifun, ispin, 2, 0)
        end do
     end do
     write(iout, *)
  end if

  deallocate(oppno)

end subroutine hf_print_op1e
!######################################################################
