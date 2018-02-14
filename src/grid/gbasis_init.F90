!#########################################################
!program read_gaussian
!  use grid_mod
!  implicit none
!  character(len=256) :: fname
!  fname=trim('He.out')
!  call gbasis_init(fname)
!  call gbasis_ovlp(fname)
!  call gbasis_core(fname)
!  call gbasis_eris(fname)
!  call gbasis_final
!end program read_gaussian
!#########################################################
subroutine gbasis_init(fname)
  use root_mod, only : iprint
  use grid_mod, only : ngbas,ngeris,govlp,goinv,gcore,geris,gmap2,gmap4
  implicit none

  character(len=*) :: fname
  integer :: mu,nu,mu2,nu2,mu3,nu3,I4,J4,munu,munu2
  real(kind(0d0)) :: s12r
  real(kind(0d0)),allocatable :: xmat(:,:),umat(:,:),ERI1(:,:,:,:),ERI2(:,:,:,:)
  
  ! read G09 information
  call gbasis_alloc(fname)
  call gbasis_read_ovlp(fname)
  call gbasis_read_core(fname)
  call gbasis_read_eris(fname)

  ! orthonormalize AOs
  allocate(xmat(1:ngbas,1:ngbas))
  allocate(umat(1:ngbas,1:ngbas))
  allocate(ERI1(1:ngbas,1:ngbas,1:ngbas,1:ngbas))
  allocate(ERI2(1:ngbas,1:ngbas,1:ngbas,1:ngbas))
  xmat = govlp
  call util_diag_real(.true.,ngbas,xmat,umat)
  goinv = 0d0
  do mu2=1, ngbas
     s12r = 1d0/sqrt(xmat(mu2,mu2))
     do mu=1, ngbas
        do nu=1, ngbas
           goinv(mu,nu)=goinv(mu,nu)+umat(mu,mu2)*s12r*umat(nu,mu2)
        end do
     end do
  end do
!DEBUG
!  umat = matmul(goinv,govlp)
!  xmat = matmul(umat,goinv)
!  write(6,"('1/sqrt(S)*S/sqrt(S):')")
!  do mu=1, ngbas
!     do nu=1, ngbas
!        write(6,"(2i5,E20.12)") mu,nu,xmat(mu,nu)
!     end do
!  end do
!DEBUG

  ! transformed core hamiltonian
  xmat = 0d0
  do mu=1, ngbas
     do nu=1, ngbas
        do mu2=1, ngbas
           xmat(mu,nu)=xmat(mu,nu)+gcore(mu,mu2)*goinv(mu2,nu)
        end do
     end do
  end do
  gcore = 0d0
  do mu=1, ngbas
     do nu=1, ngbas
        do mu2=1, ngbas
           gcore(mu,nu)=gcore(mu,nu)+goinv(mu2,mu)*xmat(mu2,nu)
        end do
     end do
  end do

  ! transformed ERIs
  ERI1 = 0d0
  do mu=1, ngbas
     do nu=1, ngbas
        do mu2=1, ngbas
           do nu2=1, ngbas
              !I4=gmap4(gmap2(mu,nu),gmap2(mu2,nu2))
              do mu3=1, ngbas
                 J4=gmap4(gmap2(mu,nu),gmap2(mu3,nu2))
                 ERI1(mu,nu,mu2,nu2)=ERI1(mu,nu,mu2,nu2)+geris(J4)*goinv(mu3,mu2)
              end do
           end do
        end do
     end do
  end do
  ERI2 = 0d0
  do mu=1, ngbas
     do nu=1, ngbas
        do mu2=1, ngbas
           do nu2=1, ngbas
              !I4=gmap4(gmap2(mu,nu),gmap2(mu2,nu2))
              do nu3=1, ngbas
                 !J4=gmap4(gmap2(mu,nu),gmap2(mu2,nu3))
                 ERI2(mu,nu,mu2,nu2)=ERI2(mu,nu,mu2,nu2)+ERI1(mu,nu,mu2,nu3)*goinv(nu3,nu2)
              end do
           end do
        end do
     end do
  end do
  ERI1 = 0d0
  do mu=1, ngbas
     do nu=1, ngbas
        do mu2=1, ngbas
           do nu2=1, ngbas
              !I4=gmap4(gmap2(mu,nu),gmap2(mu2,nu2))
              do mu3=1, ngbas
                 !J4=gmap4(gmap2(mu3,nu),gmap2(mu2,nu2))
                 ERI1(mu,nu,mu2,nu2)=ERI1(mu,nu,mu2,nu2)+goinv(mu3,mu)*ERI2(mu3,nu,mu2,nu2)
              end do
           end do
        end do
     end do
  end do
  ERI2 = 0d0
  do mu=1, ngbas
     do nu=1, ngbas
        do mu2=1, ngbas
           do nu2=1, ngbas
              !I4=gmap4(gmap2(mu,nu),gmap2(mu2,nu2))
              do nu3=1, ngbas
                 !J4=gmap4(gmap2(mu,nu3),gmap2(mu2,nu2))
                 ERI2(mu,nu,mu2,nu2)=ERI2(mu,nu,mu2,nu2)+goinv(nu3,nu)*ERI1(mu,nu3,mu2,nu2)
              end do
           end do
        end do
     end do
  end do
  geris = 0d0
  do mu=1, ngbas
     do nu=1, mu
        munu=(mu*(mu-1))/2+nu
        do mu2=1, ngbas
           do nu2=1, mu2
              munu2=(mu2*(mu2-1))/2+nu2
              if (munu >= munu2) then
                 geris(gmap4(munu,munu2))=ERI2(mu,nu,mu2,nu2)
              end if
           end do
        end do
     end do
  end do
  deallocate(ERI1)
  deallocate(ERI2)
  deallocate(umat)
  deallocate(xmat)

  if (iprint > 4) then
     write(6,"('transformed CORE:')")
     do mu=1, ngbas
        do nu=1, ngbas
           write(6,"(2i5,E20.12)") mu,nu,gcore(mu,nu)
        end do
     end do
     write(6,"('transformed ERIs')")
     do mu=1, ngbas
        do nu=1, ngbas
           do mu2=1, ngbas
              do nu2=1, ngbas
                 write(6,"(4i5,E20.12)") mu,nu,mu2,nu2,geris(gmap4(gmap2(mu,nu),gmap2(mu2,nu2)))
              end do
           end do
        end do
     end do
     !stop 'for debug @ gbasis_init.'
  end if
  
end subroutine gbasis_init
!#########################################################
subroutine gbasis_alloc(fname)
  use mol_mod, only : enen
  use grid_mod, only : ngbas,ngb2,ngeris,govlp,goinv,gcore,geris,gmap2,gmap4
  implicit none
  character(len=*) :: fname
  character(len=256) :: line
  integer :: i,j,k,l,ioerr,ij,kl,ngbasx

  open(unit=1,file=trim(fname),status='old',form='formatted')
  ! read nuclear repulsion energy
  rewind(unit=1)
  ioerr=0
  do while (ioerr==0)
     read(1,"(A)",iostat=ioerr) line
     if (line(1:31)=='       nuclear repulsion energy') then
        read(line(32:52),*) enen
        exit
     end if
  end do
  write(6,"('enen  = ',f20.10)") enen

  ! read number of basis functions
  rewind(unit=1)
  ioerr=0
  do while (ioerr==0)
     read(1,"(A)",iostat=ioerr) line
     if (line(1:12)=='    NBasis =') then
        read(line(13:16),*) ngbasx
        exit
     end if
  end do
  close(unit=1)
  if (ngbas .ne. ngbasx) then
     write(6,"('ngbas  = ',i5)") ngbas
     write(6,"('ngbasx = ',i5)") ngbasx
     stop
  end if

  ngb2=(ngbas*(ngbas+1))/2
  allocate(gmap2(1:ngbas,1:ngbas))
  allocate(gmap4(1:ngb2,1:ngb2))
  do i=1, ngbas
     do j=1, ngbas
        if (i >= j) then
           ij=(i*(i-1))/2+j
        else
           ij=(j*(j-1))/2+i
        end if
        gmap2(i,j) = ij
     end do
  end do

  ngeris=0
  do i=1, ngbas
     do j=1, i
        ij=(i*(i-1))/2+j
        do k=1, ngbas
           do l=1, k
              kl=(k*(k-1))/2+l
              if (ij >= kl) then
                 ngeris=ngeris+1
                 gmap4(ij,kl) = ngeris
              end if
           end do
        end do
     end do
  end do
  do ij=2, ngb2
     do kl=1, ij-1
        gmap4(kl,ij) = gmap4(ij,kl)
     end do
  end do
  allocate(govlp(1:ngbas,1:ngbas))
  allocate(goinv(1:ngbas,1:ngbas))
  allocate(gcore(1:ngbas,1:ngbas))
  allocate(geris(1:ngeris))
!  write(6,"('gmap2:')")
!  do i=1, ngbas
!     do j=1, ngbas
!        write(6,"(i5)",advance='no') gmap2(i,j)
!     end do
!     write(6,*)
!  end do
!
!  write(6,"('gmap4:')")
!  do i=1, ngb2
!     do j=1, ngb2
!        write(6,"(i5)",advance='no') gmap4(i,j)
!     end do
!     write(6,*)
!  end do
!
!  write(6,"('gmap from quartets to serial ERI orders:')")
!  do i=1, ngbas
!     do j=1, ngbas
!        do k=1, ngbas
!           do l=1, ngbas
!              write(6,"(8i5)") i,j,k,l,gmap2(i,j),gmap2(k,l),gmap4(gmap2(i,j),gmap2(k,l))
!           end do
!        end do
!     end do
!  end do
end subroutine gbasis_alloc
!#########################################################
subroutine gbasis_read_eris(fname)
  use root_mod, only : iprint
  use grid_mod, only : ngbas,ngeris,geris,gmap2,gmap4
  implicit none
  character(len=*) :: fname
  character(len=256) :: line
  character(len=10) :: cxx
  integer :: i,j,k,l,ij,kl,ioerr,ngerisx,ixx
  real(kind(0d0)) :: val

  geris=0d+0
  open(unit=1,file=trim(fname),status='old',form='formatted')
  rewind(unit=1)
  ioerr=0
  do while (ioerr==0)
     read(1,"(A)",iostat=ioerr) line
     if (line(1:39)==' *** Dumping Two-Electron integrals ***') exit
  end do
  read(1,*); read(1,*); read(1,*); read(1,*); read(1,*);
  read(1,*) cxx,ixx,cxx,ngerisx
  do ixx = 1, ngerisx
     read(1,*) cxx,i,cxx,j,cxx,k,cxx,l,cxx,val
     geris(gmap4(gmap2(i,j),gmap2(k,l))) = val
  end do
  close(unit=1)
  if (iprint > 4) then
     write(6,"('ERIs',4i10)") ngbas,ngerisx,ngeris,ngbas**4
     do i=1, ngbas
        do j=1, ngbas
           do k=1, ngbas
              do l=1, ngbas
                 write(6,"(4i5,E20.12)") i,j,k,l,geris(gmap4(gmap2(i,j),gmap2(k,l)))
              end do
           end do
        end do
     end do
  end if
end subroutine gbasis_read_eris
!#########################################################
subroutine gbasis_read_ovlp(fname)
  use root_mod, only : iprint
  use grid_mod, only : ngbas,govlp
  implicit none
  character(len=*) :: fname
  character(len=256) :: line
  character(len=10) :: cxx
  integer :: i,j,ioerr,ixx,nitr,iitr,ill
  real(kind(0d0)) :: val

  govlp=0d+0
!old  goinv=0d+0
  nitr = ngbas/5 + min(mod(ngbas,5),1)

  open(unit=1,file=trim(fname),status='old',form='formatted')
  rewind(unit=1)
  ioerr=0
  do while (ioerr==0)
     read(1,"(A)",iostat=ioerr) line
     if (line(1:16)==' *** Overlap ***') exit
  end do
  do iitr=1, nitr
     read(1,"(A)") line
     ill=5*(iitr-1)+1
     do i=ill, ngbas
        read(1,*) ixx, govlp(i,ill:min(i,ill+4))
     end do
  end do
  close(unit=1)
  do i=1, ngbas-1
     do j=i+1, ngbas
        govlp(i,j) = govlp(j,i)
     end do
  end do

!old  ! inverse overlap
!old  call util_matinv_real(ngbas,govlp,goinv)

  if (iprint > 4) then
     write(6,"('OVLP:')")
     do i=1, ngbas
        do j=1, ngbas
           write(6,"(2i5,E20.12)") i,j,govlp(i,j)
        end do
     end do
!old     write(6,"('OINV:')")
!old     do i=1, ngbas
!old        do j=1, ngbas
!old           write(6,"(2i5,E20.12)") i,j,goinv(i,j)
!old        end do
!old     end do
     !STOP "for debug @ gbasis_read_ovlp."
  end if

end subroutine gbasis_read_ovlp
!#########################################################
subroutine gbasis_read_core(fname)
  use root_mod, only : iprint
  use grid_mod, only : ngbas,gcore
  implicit none
  character(len=*) :: fname
  character(len=256) :: line
  character(len=10) :: cxx
  integer :: i,j,ioerr,ixx,nitr,iitr,ill
  real(kind(0d0)) :: val

  gcore=0d+0
  nitr = ngbas/5 + min(mod(ngbas,5),1)

  open(unit=1,file=trim(fname),status='old',form='formatted')
  rewind(unit=1)
  ioerr=0
  do while (ioerr==0)
     read(1,"(A)",iostat=ioerr) line
     if (line(1:31)==' ****** Core Hamiltonian ******') exit
  end do
  do iitr=1, nitr
     read(1,*)
     ill=5*(iitr-1)+1
     do i=ill, ngbas
        read(1,*) ixx, gcore(i,ill:min(i,ill+4))
     end do
  end do
  close(unit=1)
  do i=1, ngbas-1
     do j=i+1, ngbas
        gcore(i,j) = gcore(j,i)
     end do
  end do
  if (iprint > 4) then
     write(6,"('CORE:')")
     do i=1, ngbas
        do j=1, ngbas
           write(6,"(2i5,E20.12)") i,j,gcore(i,j)
        end do
     end do
  end if
end subroutine gbasis_read_core
!#########################################################
