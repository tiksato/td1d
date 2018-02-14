!######################################################################
subroutine hf_read(iunit, wfn)

  use const_mod, only : czero
  use grid_mod, only : ngrid, x
  use wfn_mod, only : nfun, nspin
  use init_mod, only : guess_interpolate

  implicit none
  integer, intent(in) :: iunit
  complex(kind(0d0)), intent(out) :: wfn(0:ngrid, 1:nfun, 1:nspin)

  integer :: ngrid0, nval, nshift, igrid0, igrid, ifun, jfun, ispin
  real(kind(0d0)), allocatable :: x1(:)
  complex(kind(0d0)) :: unorm
  complex(kind(0d0)), external :: util_zdotc
  complex(kind(0d0)), allocatable :: wfn1(:,:,:), wfn2(:,:,:)
  complex(kind(0d0)), allocatable :: uvec(:,:,:)

  read(iunit, "(I25)") ngrid0
  allocate(x1(0:ngrid0))
  allocate(wfn1(0:ngrid0, 1:nfun, 1:nspin))

  read(iunit, "( E25.15)") x1(0:ngrid0)
  read(iunit, "(2E25.15)") wfn1(0:ngrid0, 1:nfun, 1:nspin)

  if (guess_interpolate) then
     nval = nfun * nspin
     call util_proj(ngrid0, ngrid, nval, wfn1, wfn)
     call hf_ort(wfn)
  else
     if (ngrid >= ngrid0) then
        nshift = (ngrid - ngrid0) / 2
        do igrid0 = 0, ngrid0
           igrid = igrid0 + nshift
           wfn(igrid, 1:nfun, 1:nspin) = wfn1(igrid0, 1:nfun, 1:nspin)
        end do
     else
        nshift = (ngrid0 - ngrid) / 2
        do igrid = 0, ngrid
           igrid0 = igrid + nshift
           wfn(igrid, 1:nfun, 1:nspin) = wfn1(igrid0, 1:nfun, 1:nspin)
        end do
     end if
  end if

  deallocate(wfn1)
  deallocate(x1)

end subroutine hf_read
!######################################################################
