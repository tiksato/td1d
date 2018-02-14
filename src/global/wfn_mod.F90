!################################################################################
module wfn_mod

  implicit none

  public
  complex(kind(0d0)), allocatable :: wfn(:)

  logical :: hcore, sae, crapola, ehf
  logical :: cionly, moonly
  logical :: apsg, postapsg
  logical :: ldacx
  logical :: ooinit_minus
  logical :: normci
  logical :: semicanonical

  character(len=16) :: method
  integer :: imethod
  integer :: nspin
  integer :: nfun
  integer :: nocc
  integer :: n2ov
  integer :: nvir, nvira, nvirx
  integer :: nelcore(3), nelact(3)
  integer :: nfcore, ndcore, ncore, nact, nfact, nopen
  integer :: sep_fc
  logical :: fc_exc

  integer :: nblock
  integer, parameter :: max_nblock = 10
  integer :: type_block(max_nblock)
  integer :: nfun_block(max_nblock)

  logical :: doci, domo, domoq, domop, domop1, domop2, docip2
  logical :: saex
  logical :: crapx
  integer :: nfroz(2)

  logical :: dopt
  logical :: tcham

  real(kind(0d0)) :: exps2
  real(kind(0d0)) :: energy
  real(kind(0d0)) :: eneref ! ene0 + ene1
  real(kind(0d0)) :: ene0   ! 0th order energy
  real(kind(0d0)) :: ene1   ! 1st order energy
  real(kind(0d0)) :: ene2   ! 2nd order energy
  real(kind(0d0)) :: ene3   ! energy expectation value of 0+1th order wavefunction
  real(kind(0d0)) :: norm   ! squared norm of 1st order wavefunction correction
  logical :: mrmp_h0_nofield
  character(len=16) :: mrmp_pt0
  character(len=16) :: mrmp_proj_type

  integer :: max_krylov
  integer :: cmf_maxorb
  integer :: cmf_maxcic
  integer :: cmf_nddt
  integer :: cmf_ncorr
  character(len=16) :: cmf_type

  ! ormas parameters
  logical :: ormas, ormas0, rasscf, qcas, donly, ormas_sd1, ormas_allowed
  logical :: noci, nomo, nolag1, nolag2, nolag3, qcouple
  logical :: s2zero
  integer, parameter :: nsub_max = 10
  integer :: nstate, target_state
  integer :: nsub
  integer :: norb_sub(1:nsub_max), lorb_sub(1:2, 1:nsub_max)
  integer :: nel_sub_alph(1:nsub_max), nel_sub_beta(1:nsub_max)
  integer :: min_sub(1:nsub_max), min_sub_alph(1:nsub_max), min_sub_beta(1:nsub_max)
  integer :: max_sub(1:nsub_max), max_sub_alph(1:nsub_max), max_sub_beta(1:nsub_max)
  integer :: ndist, ndist_alph, ndist_beta
  integer :: ndet, nstr_alph, nstr_beta
  integer, allocatable :: sub_orb(:)
  integer, allocatable :: dist(:,:), dist_alph(:,:), dist_beta(:,:)
  integer, allocatable :: det_allowed(:,:), ndet_dist(:), mapf_det(:,:), mapr_det(:,:)
  integer, allocatable :: nstr_alph_dist(:), lstr_alph_dist(:,:), nstr_alph_dist_sub(:,:)
  integer, allocatable :: nstr_beta_dist(:), lstr_beta_dist(:,:), nstr_beta_dist_sub(:,:)
  ! ##### obsolete #####
  integer, allocatable :: arcwgt(:,:)
  ! ##### obsolete #####
  integer, allocatable :: arc_alph(:,:)
  integer, allocatable :: arc_beta(:,:)
  integer, allocatable :: dist_str_alph(:,:), substr_alph(:,:)
  integer, allocatable :: dist_str_beta(:,:), substr_beta(:,:)
  integer, allocatable :: onv_alph(:,:), orb_alph(:,:)
  integer, allocatable :: onv_beta(:,:), orb_beta(:,:)
  integer, allocatable :: n1x_alph(:,:), p1x_alph(:,:), h1x_alph(:,:), eq1x_alph(:,:), &
       & sgn1x_alph(:,:), n1xr_alph(:,:), r1xr_alph(:,:,:), l1xr_alph(:,:,:), &
       & sgn1xr_alph(:,:,:)
  integer, allocatable :: n1x_beta(:,:), p1x_beta(:,:), h1x_beta(:,:), eq1x_beta(:,:), &
       & sgn1x_beta(:,:), n1xr_beta(:,:), r1xr_beta(:,:,:), l1xr_beta(:,:,:), &
       & sgn1xr_beta(:,:,:)
  integer :: nrotoo, nrotca, nrotaa, nrotov
  integer, allocatable :: rotoo_mapf(:,:), rotoo_mapb(:,:)
  integer, allocatable :: rotca_mapf(:,:), rotca_mapb(:,:)
  integer, allocatable :: rotaa_mapf(:,:), rotaa_mapb(:,:)
  integer, allocatable :: rotov_mapf(:,:), rotov_mapb(:,:)
  real(kind(0d0)) :: max_grad(1:3)
  real(kind(0d0)) :: rms_grad(0:3)

  ! ras parameters
  logical :: doras
  integer :: nact1, nact2, nact3
  integer :: max_hole, max_hole1, max_hole2
  integer :: max_elec, max_elec1, max_elec2
  integer :: max_type, ntypea, ntypeb
  integer :: maxdpsi
  logical :: dpsi_itr, dpsi_nod2x, dpsi_nocp
  integer :: dpsi_pene, dpsi_ncut, dpsi_reg
  real(kind(0d0)) :: dpsi_smin
  real(kind(0d0)) :: dpsi_damp, dpsi_dsv, dpsi_dyreg0, dpsi_dyreg1, dpsi_dyreg2, dpsi_maxang, dpsi_test
  integer :: nrota
  integer, allocatable :: map2f(:,:), map2b(:,:)
  integer, allocatable :: nstra_type(:), llstra(:), ulstra(:)
  integer, allocatable :: nstrb_type(:), llstrb(:), ulstrb(:)
  integer :: npair
  integer, allocatable :: nfunpp(:), pp2fun(:), fun2pp(:)
  complex(kind(0d0)), allocatable :: occn(:), coeffa(:,:), coeffb(:,:), coeffc(:,:)

  logical :: reg_occd
  integer :: reg_type
  integer :: max_ipd, max_ipx

! method specific options
  logical :: hf_doproj
  logical :: apsg_sepphase

  integer :: hcic_type
  integer :: den1_type
  integer :: den2_type

! Time-dependent Coupled-Cluster
  logical :: tdcc
  integer :: norb1, norb2
  integer :: nbiort

! For gaussian-basis calculation
  complex(kind(0d0)), allocatable :: int1g(:,:), int2g(:,:,:,:)

!OBSOLETE
  logical :: nolag, nolagx, nolag2x
  integer :: nstra                       ! number of alpha string
  integer, allocatable :: wgta(:,:)      ! arc weight for alpha string graph
  integer, allocatable :: onva(:,:)      ! alpha occupation number vector
  integer, allocatable :: orba(:,:)      ! alpha orbital index vector
  integer, allocatable :: n1xa(:)        ! number of nonzero internal single excitations (1x)
  integer, allocatable :: p1xa(:,:)      ! particle indecies for each 1x of from each string
  integer, allocatable :: h1xa(:,:)      ! hole indecies for each 1x of from each string
  integer, allocatable :: eq1xa(:,:)     ! indexing the 1x-resultant string from each string
  integer, allocatable :: sgn1xa(:,:)    ! phase connecting eq1x and original string
  integer, allocatable :: ord1xa(:,:)    ! two-dimensional ordering of each excitation
  integer, allocatable :: n1xra(:,:)     ! number of strings surviving excitation j -> i.
  integer, allocatable :: r1xra(:,:,:)   ! index of string before applying Eij
  integer, allocatable :: l1xra(:,:,:)   ! index of string after applying Eij
  integer, allocatable :: sgn1xra(:,:,:) ! phase connecting r1xr and l1xr
  integer, allocatable :: ord1xra(:,:,:) ! two-dimensional ordering of each excitation

  integer :: nstrb
  integer, allocatable :: wgtb(:,:)
  integer, allocatable :: onvb(:,:)
  integer, allocatable :: orbb(:,:)
  integer, allocatable :: n1xb(:)
  integer, allocatable :: p1xb(:,:)
  integer, allocatable :: h1xb(:,:)
  integer, allocatable :: eq1xb(:,:)
  integer, allocatable :: sgn1xb(:,:)
  integer, allocatable :: ord1xb(:,:)
  integer, allocatable :: n1xrb(:,:)
  integer, allocatable :: r1xrb(:,:,:)
  integer, allocatable :: l1xrb(:,:,:)
  integer, allocatable :: sgn1xrb(:,:,:)
  integer, allocatable :: ord1xrb(:,:,:)
!OBSOLETE

end module wfn_mod
!################################################################################
