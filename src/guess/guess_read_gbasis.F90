!################################################################################
subroutine guess_read_gbasis(wfn, imethod)

  use grid_mod, only : ngbas,g09out
  use wfn_mod, only : nfun,nspin,tdcc

  implicit none
  !--------------------------------------------------------------------
  integer, intent(in) :: imethod
  complex(kind(0d0)), intent(out) :: wfn(1:*)
  !--------------------------------------------------------------------
  integer, external :: wfn_size
  integer, external :: wfn_poscic

  wfn(1:wfn_size(imethod)) = 0d0
  if (imethod == 2 .and. .not.tdcc) wfn(wfn_poscic(imethod))=1d0
  call guess_read_gbasis_orb(wfn)
  call wfn_ort_mo(wfn, imethod)

end subroutine guess_read_gbasis
!################################################################################
subroutine guess_read_gbasis_orb(wfn)

  use grid_mod, only : ngbas,g09out,goinv
  use wfn_mod, only : nfun,nspin

  implicit none
  !--------------------------------------------------------------------
  complex(kind(0d0)), intent(out) :: wfn(1:ngbas, 1:nfun, 1:nspin)
  !--------------------------------------------------------------------
  character(len=256) :: line
  integer :: mu,nu,ifun,ispin,nitr,iitr,ill,ioerr
  real(kind(0d0)), allocatable :: go12(:,:)
  complex(kind(0d0)), allocatable :: wtmp(:)
  real(kind(0d0)) :: wread(1:nfun)

  if (nspin == 2) stop 'guess_core only supports nspin = 1.'
  nitr = nfun/5 + min(mod(nfun,5),1)

  open(unit=1,file=trim(g09out),status='old',form='formatted')
  rewind(unit=1)
  ioerr=0

  write(6,"('Searching for MO coefficients....  ')",advance='no')
  do while (ioerr==0)
     read(1,"(A)",iostat=ioerr) line
     if (line(1:36)=='     Molecular Orbital Coefficients:') then
        write(6,"('FOUND!')")
        exit
     end if
  end do
  write(6,*)
  if (ioerr.ne.0) stop 'Failed to find MO coefficients.'

!debug
!write(6,"('ngbas =',i5)") ngbas
!write(6,"('nfun  =',i5)") nfun
!write(6,"('nitr  =',i5)") nitr
!debug

  do iitr=1, nitr
     read(1,"(A)") line; read(1,"(A)") line; read(1,"(A)") line
     ill=5*(iitr-1)+1
     do mu=1, ngbas
        read(1,"(A)") line
   !bug read(line(22:71),*) wfn(mu,ill:min(nfun,ill+4),1)
        read(line(22:71),*) wread(ill:min(nfun,ill+4))
        wfn(mu,ill:min(nfun,ill+4),1) = wread(ill:min(nfun,ill+4))
     end do
  end do
  close(unit=1)

!debug
!stop 'I am here!'
!debug

  !DEBUG
  !write(6,"('guess_read_gbasis: MO coefficients before transformation.')")
  !do mu=1, ngbas
  !   write(6,"(i5)",advance='no') mu
  !   do ifun = 1, nfun
  !      write(6,"(2f14.10)",advance='no') wfn(mu,ifun,1)
  !   end do
  !   write(6,*)
  !end do
  !stop
  !DEBUG

  !transform to the representation in orthonormalized basis
  allocate(wtmp(1:ngbas))
  allocate(go12(1:ngbas,1:ngbas))
  call util_matinv_real(ngbas,goinv,go12)
  do ifun=1, nfun
     wtmp = 0d0
     do mu=1, ngbas
        do nu=1, ngbas
           wtmp(mu)=wtmp(mu)+go12(mu,nu)*wfn(nu,ifun,1)
        end do
     end do
     wfn(:,ifun,1)=wtmp
  end do
  deallocate(go12)
  deallocate(wtmp)

  !DEBUG
  !write(6,"('guess_read_gbasis: MO coefficients after transformation.')")
  !do mu=1, ngbas
  !   write(6,"(i5)",advance='no') mu
  !   do ifun = 1, nfun
  !      write(6,"(2f14.10)",advance='no') wfn(mu,ifun,1)
  !   end do
  !   write(6,*)
  !end do
  !stop
  !DEBUG

end subroutine guess_read_gbasis_orb
!######################################################################
