!################################################################################
module cc_mod

  implicit none

  public
  logical :: tonly, cc_solve
  integer :: cc_solve_itr, cc_maxcyc
  integer :: cc_rank, cc_read_rank, cc_xij_xab
  logical :: occd, brueckner, optbrueckner, bccd, obccd, obccd2, occsd, ccdl1, ccsdl1
  logical :: cc_sreal
  logical :: cc_l1_disc, cc_l2_disc
  logical :: cc_read_ort, cc_read_donly, cc_nonredundant
  integer :: len_tcc1, len_gcc1, len_tcc2, len_gcc2, len_tcc3, len_gcc3
  integer :: ind_tcc1, ind_gcc1, ind_tcc2, ind_gcc2, ind_tcc3, ind_gcc3

  real(kind(0d0)) :: thrtamp
  real(kind(0d0)) :: thrgamp

  complex(kind(0d0)), allocatable :: cc_i0(:)
  complex(kind(0d0)), allocatable :: cc_i1(:)
  complex(kind(0d0)), allocatable :: cc_i2(:)
  complex(kind(0d0)), allocatable :: cc_i3(:)

  complex(kind(0d0)), allocatable :: fock(:,:,:)
  complex(kind(0d0)), allocatable :: int2x(:,:,:,:,:)
  complex(kind(0d0)), allocatable :: tamp1(:,:,:)
  complex(kind(0d0)), allocatable :: gamp1(:,:,:)
  complex(kind(0d0)), allocatable :: tamp2(:,:,:,:,:)
  complex(kind(0d0)), allocatable :: gamp2(:,:,:,:,:)
  complex(kind(0d0)), allocatable :: tamp3(:,:,:,:,:,:,:)
  complex(kind(0d0)), allocatable :: gamp3(:,:,:,:,:,:,:)
  complex(kind(0d0)), allocatable :: dtamp1(:,:,:)
  complex(kind(0d0)), allocatable :: dgamp1(:,:,:)
  complex(kind(0d0)), allocatable :: dtamp2(:,:,:,:,:)
  complex(kind(0d0)), allocatable :: dgamp2(:,:,:,:,:)
  complex(kind(0d0)), allocatable :: dtamp3(:,:,:,:,:,:,:)
  complex(kind(0d0)), allocatable :: dgamp3(:,:,:,:,:,:,:)

  complex(kind(0d0)), allocatable :: den1s(:,:,:)
  complex(kind(0d0)), allocatable :: den2s(:,:,:,:,:)
  complex(kind(0d0)), allocatable :: den1_noref(:,:,:)
  complex(kind(0d0)), allocatable :: den2_noref(:,:,:,:,:)

  complex(kind(0d0)), allocatable :: itm_hh(:,:)
  complex(kind(0d0)), allocatable :: itm_pp(:,:)
  complex(kind(0d0)), allocatable :: itm_hp(:,:)
  complex(kind(0d0)), allocatable :: itm_ph(:,:)
  complex(kind(0d0)), allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)), allocatable :: itm_pppp(:,:,:,:)
  complex(kind(0d0)), allocatable :: itm_pphp(:,:,:,:)
  complex(kind(0d0)), allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)), allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)), allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)), allocatable :: itm_hpph(:,:,:,:)

  complex(kind(0d0)), allocatable :: itm_hhh(:,:,:)
  complex(kind(0d0)), allocatable :: itm_hph(:,:,:)
  complex(kind(0d0)), allocatable :: itm_pph(:,:,:)
  complex(kind(0d0)), allocatable :: itm_ppp(:,:,:)
  complex(kind(0d0)), allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)), allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)), allocatable :: itm_phpp(:,:,:,:)
  complex(kind(0d0)), allocatable :: itm_ppph(:,:,:,:)
  complex(kind(0d0)), allocatable :: itm_hhphhh(:,:,:,:,:,:)

end module cc_mod
!################################################################################
