!################################################################################
subroutine input_wfn(ioin)

  use mol_mod, only : ne
  use const_mod, only : zero, one, two
  use wfn_mod, only : max_nblock, imethod
  use wfn_mod, only : method, nfun, fc_exc, sep_fc, nfcore, ndcore, ncore, nact, nfact, nopen, nocc, nvira, nvirx, &
       & nvir, n2ov, nelcore, nelact, nspin, nblock, type_block, nfun_block, hcore, sae, crapola, &
       & ehf, cionly, moonly, apsg, postapsg, dopt, tcham, mrmp_pt0, mrmp_proj_type, mrmp_h0_nofield, &
       & npair, noci, nomo, nolag1, nolag2, nolag3, qcouple, ldacx, reg_occd, reg_type, max_ipd, max_ipx, &
       & apsg_sepphase, doras, max_hole, max_elec, nact1, nact2, nact3, dpsi_itr, dpsi_nod2x, dpsi_reg, &
       & dpsi_nocp, dpsi_pene, dpsi_ncut, dpsi_damp, dpsi_dsv, dpsi_dyreg0, dpsi_dyreg1, dpsi_dyreg2, &
       & dpsi_maxang, maxdpsi, max_hole1, max_elec1, ormas, ormas0, rasscf, qcas, nsub_max, nsub, &
       & norb_sub, lorb_sub, min_sub, max_sub, nel_sub_alph, nel_sub_beta, nstate, target_state, s2zero, &
       & hf_doproj, donly, ormas_sd1, ormas_allowed, ooinit_minus, max_krylov, cmf_maxorb, cmf_maxcic, normci, hcic_type, &
       & den1_type, den2_type, cmf_nddt, cmf_ncorr, cmf_type, tdcc, norb1, norb2, nbiort, semicanonical
  use cc_mod, only : cc_rank, cc_read_rank, len_tcc1, len_gcc1, len_tcc2, len_gcc2, len_tcc3, len_gcc3, &
       ind_tcc1, ind_gcc1, ind_tcc2, ind_gcc2, ind_tcc3, ind_gcc3, cc_l1_disc, cc_l2_disc, &
       cc_sreal, cc_read_ort, cc_read_donly, ccdl1, ccsdl1, occd, brueckner, optbrueckner, bccd, obccd, obccd2, occsd, tonly, cc_solve, &
       cc_solve_itr, thrtamp, thrgamp, cc_maxcyc, cc_nonredundant, cc_xij_xab
!OBSOLETE
  use wfn_mod, only : nolag
!OBSOLETE

  implicit none
  integer, intent(in) :: ioin
  integer :: ioerr, nfun_tot, iblock, nfcorex, ndcorex, ncorex, nactx, nopenx, &
       & nvira1, nvirx1, npair1, isub, ntota, ntotb
  namelist /wfn/ method, nfun, fc_exc, sep_fc, nfcore, ndcore, ncore, nact, nfact, nopen, nocc, nvira, nvirx, nvir, n2ov, &
       & nelcore, nelact, nspin, nblock, type_block, nfun_block, hcore, sae, crapola, ehf, &
       & cionly, moonly, apsg, postapsg, dopt, tcham, mrmp_pt0, mrmp_proj_type, mrmp_h0_nofield, &
       & npair, noci, nomo, nolag1, nolag2, nolag3, qcouple, ldacx, reg_occd, reg_type, max_ipd, max_ipx, &
       & apsg_sepphase, doras, max_hole, max_elec, nact1, nact2, nact3, dpsi_itr, dpsi_nod2x, dpsi_reg, &
       & dpsi_nocp, dpsi_pene, dpsi_ncut, dpsi_damp, dpsi_dsv, dpsi_dyreg0, dpsi_dyreg1, dpsi_dyreg2, &
       & dpsi_maxang, maxdpsi, max_hole1, max_elec1, ormas, ormas0, rasscf, qcas, nsub, &
       & norb_sub, min_sub, max_sub, nel_sub_alph, nel_sub_beta, nstate, target_state, s2zero, &
       & hf_doproj, donly, ormas_sd1, ormas_allowed, ooinit_minus, max_krylov, cmf_maxorb, cmf_maxcic, normci, hcic_type, &
       & den1_type, den2_type, cmf_nddt, cmf_ncorr, cmf_type, tdcc, cc_l1_disc, cc_l2_disc, cc_sreal, nbiort, cc_read_ort, &
       & cc_read_donly, ccdl1, ccsdl1, occd, brueckner, optbrueckner, bccd, obccd, obccd2, occsd, tonly, cc_read_rank, semicanonical, cc_solve, &
       & cc_solve_itr, thrtamp, thrgamp, cc_maxcyc, cc_nonredundant, cc_xij_xab

  ! essential
  method = ''                   ! HF, APSG, APSG-CI, CASSCF, CASCI, (2E: X2E, GVB, CIS)
  nblock = -2                   ! number of distinct spaces of orbitals
  type_block(1:max_nblock) = -2 ! type of the orbital spaces, -1(fc)/0(dc)/1(act)/2(va)/3(vx)
  nfun_block(1:max_nblock) = -2 ! number of orbitals in each block

  ! optional
  nfun = -1                 ! number of orbitals
  fc_exc = .true.           ! do frozen-core exchange
  sep_fc = 0                ! independent frozen-core
  nfcore = -1               ! number of frozen core orbitals
  ndcore = -1               ! number of dynamical core orbitals
  ncore = -1                ! nfcore + ndcore
  nact = -1                 ! number of active orbitals
  nfact = 0                 ! number of frozen active orbitals
  nopen = -1                ! number of singly occupied orbitals
  nocc = -1                 ! ncore + nact + nopen
  nvira = -1                ! number of primary virtual orbitals
  nvirx = -1                ! number of secondary virtual orbitals
  nvir = -1                 ! nvira + nvirx
  npair = -1                ! number of gvb geminals
  nelcore(1:2) = -1         ! number of alpha/beta core electrons
  nelact(1:2) = -1          ! number of alpha/beta active electrons
  nspin =  1                ! restricted or unrestricted for hartree-fock
  noci = .false.            ! ci's are fixed (same as moonly = .true.)
  nomo = .false.            ! mo's are fixed (same as cionly = .true.)
  nolag1 = .false.          ! mcscf: core-active rotations are fixed
  nolag2 = .false.          ! mcscf: active-active rotations are fixed
  nolag3 = .false.          ! mcscf: occupied-virtual rotations are fixed
  qcouple = .true.          ! mcscf: couple q and other variables in aughess method
  hcore = .false.           ! only core hamiltonian
  sae = .false.             ! single active electron approximation
  crapola = .false.         ! "crapola" model
  ehf = .false.             ! extended Hartree-Fock instead of GVB
  cionly = .false.          ! fixed-orbital CI variation
  moonly = .false.          ! fixed-CI orbital variation
  apsg = .false.            ! transform to NO after each iteration of CAS
  postapsg = .false.        ! transform to NO after CAS optimization completed
  dopt = .false.            ! many-body perturbation theory
  tcham = .false.           ! transcorrelated Hamiltonian
  mrmp_pt0 = 'FULL'         ! FULL, VV
  mrmp_proj_type = 'PSI0'   ! PSI0, OV, VV
  mrmp_h0_nofield = .false. ! zero-th order Hamiltonian without laser field
  ldacx = .false.           ! local-density approximation of frozen-core exchange
  reg_occd = .true.         ! regularize the occupation number denominator in td-casscf
  reg_type = 0              ! type of the hermitian matrix regularization: 0, 1, or 2.
  max_ipd = 2               ! to which-tuple ionization probabilites estimated from RDM
  max_ipx = 0               ! to which-tuple ionization probabilites computed exactly
  max_krylov = 100          ! maximum krylov subspace dimension
  cmf_maxorb = 1            ! maximum krylov subspace dimension for orbitals
  cmf_maxcic = 1            ! maximum krylov subspace dimension for CI
  cmf_nddt = 1              ! number of subgrids of dt
  cmf_ncorr = 0             ! number of corrector steps
  cmf_type = 'RK4'          ! orbital propagation method in the cmf method. RK4 or EXP
  tdcc = .false.            ! time-dependent coupled-cluster
  cc_l1_disc = .true.       ! include disconnected (but linked) diagrams in the lambda-S equation
  cc_l2_disc = .true.       ! include disconnected (but linked) diagrams in the lambda-D equation
  cc_sreal = .true.         ! real or complex action
  nbiort = 1                ! 1: orthonormal, 2: biorthonormal
  cc_read_ort = .false.     ! read same orbitals for ket and bra
  cc_read_donly = .false.   ! read T2 and G2 from CCD output
  ccdl1 = .false.           ! CCD with L1 amplitudes
  ccsdl1 = .false.          ! CCSD with Delta^a_i = -delta T^a_i
  occd = .false.            ! OCCD, equivalent to DOnly.and.NBiOrt==1
  brueckner = .false.       ! Brueckner coupled-cluster theory
  optbrueckner = .false.    ! Optimized Brueckner coupled-cluster theory
  bccd = .false.            ! Brueckner CCD
  obccd = .false.           ! Optimized Brueckner CCD
  obccd2 = .false.          ! Optimized Brueckner CCD
  occsd = .false.           ! Optimized CCSD, SHOULD BE equivalent to CCSD
  tonly = .false.           ! Optimize t-amplitude only fixed lambda's and orbitals
  cc_read_rank = 0          ! Rank of reading CC amplitudes
  cc_solve = .false.        ! Iteratively solve T (and L) amplitude equation
  cc_solve_itr = 0          ! Solve T (and L) equation during imag. prop.
  semicanonical = .false.   ! Fock of occ and vir spaces are separately diagonalized
  thrtamp = 1.D-12          ! Threshold for T amplitude equation
  thrgamp = 1.D-12          ! Threshold for L amplitude equation
  cc_maxcyc = 20
  cc_nonredundant = .false.
  cc_xij_xab = 0

  hf_doproj = .true.        ! do or donot project against occupied orbitals
  apsg_sepphase = .true.    ! apply (1 - \gem(p)\gem^+(p)) to i\partial_t gem(p).
  ooinit_minus = .true.
  normci = .false.          ! normalize ci coefficients in real time propagation
  hcic_type = 2
  den1_type = 1
  den2_type = 2

  ! ormas settings
  s2zero = .false.              ! spin-unpolarized singlet
  ormas = .false.                ! using ormas program
  ormas0 = .false.              ! using ormas0 constraint within ormas program
  rasscf = .false.              ! using ras constraint within ormas program
  qcas = .false.                ! using qcas constraint within ormas program
  donly = .false.               ! only double transitions (miyagi and madsen)
  ormas_sd1 = .false.           ! only double transitions (miyagi and madsen)
  ormas_allowed = .false.       ! direct input of det_allowed from Name.allowed
  nsub = 1                      ! number of subactive spaces
  norb_sub(1:nsub_max) = -1     ! number of orbitals in each subactive spaces
  min_sub(1:nsub_max) = -1      ! minimum occupations
  max_sub(1:nsub_max) = -1      ! maximum occupations
  nel_sub_alph(1:nsub_max) = -1 ! numbers of alpha electrons (QCAS)
  nel_sub_beta(1:nsub_max) = -1 ! numbers of beta electrons (QCAS)
  nstate = 1                    ! number of electronic states
  target_state = 0              ! target electronic state

  ! ras settings
  doras = .false.           ! using rasci program
  max_hole = -1             ! maximum allowed number of holes in ras1
  max_elec = -1             ! maximum allowed number of elecs in ras3
  max_hole1 = -1            ! maximum allowed number of alpha/beta holes in ras1
  max_elec1 = -1            ! maximum allowed number of alpha/beta elecs in ras3
  nact1 = -1                ! number of ras1 orbitals
  nact2 = -1                ! number of ras2 orbitals
  nact3 = -1                ! number of ras3 orbitals
  maxdpsi = 100             ! maximum iteration numbers for ras_hprod_dpsidt
  dpsi_itr = .false.        ! iterative solution of active-active rotations
  dpsi_nod2x = .false.      ! neglect den2x contributions in the a-a rot equation
  dpsi_nocp = .false.       ! neglect the ci-orbital coupling
  dpsi_pene = 1             ! Pseudo energy W: 0=0, 1=<H>, 2=<H-R>
  dpsi_reg = 1              ! 2: x+d*exp(-x/d), 1: x/(x**2+d**2), 2: 1/(x+d), -1: 1/x
  dpsi_damp = 1.D-10
  dpsi_dsv = 1.D-15
  dpsi_dyreg0 = 1.D-10      ! threshold for D inverse
  dpsi_dyreg1 = 1.D-10      ! threshold for (2-D) inverse
  dpsi_dyreg2 = 1.D-10      ! threshold for A inverse
  dpsi_ncut = 0             ! 1 <= i <= nvar - dpsi_ncut: thresh1, otherwise: thresh2
  dpsi_maxang = -one

  rewind(ioin)
  read(unit=ioin, nml=wfn, iostat=ioerr)
  if(ioerr /= 0) stop "error in namelist wfn."

  call util_transchar(method)
  call util_transchar(mrmp_pt0)
  call util_transchar(mrmp_proj_type)

  if (trim(method) == 'TDSE') then
     imethod = -1
  else if (trim(method) == 'HF') then
     imethod = 0
  else if (trim(method) == 'CORE') then
     imethod = 0
     hcore = .true.
  else if (trim(method) == 'GVB') then
     imethod = 1
  else if (trim(method) == 'CAS') then
     imethod = 2
  else if (trim(method) == 'MCSCF') then
     imethod = 2
  else if (trim(method) == 'TDCC') then
     imethod = -2
  else if (trim(method) == 'CASSCF') then
     imethod = 2
  else if (trim(method) == 'CASCI') then
     imethod = 2
     cionly = .true.
  else if (trim(method) == 'APSG') then
     imethod = 3
  else if (trim(method) == 'APSG-CI') then
     imethod = 3
     cionly = .true.
  else if (trim(method) == 'TDCIS') then
     imethod = 4
  else if (trim(method) == 'CIS') then
     imethod = -4
  else if (trim(method) == 'TCHF') then
     stop 'TCHF nyi.'
  else if (trim(method) == 'TCCAS') then
     stop 'TCCAS nyi.'
  else if (trim(method) == 'MP') then
     stop 'MP nyi.'
  else if (trim(method) == 'MRMP') then
     stop 'MRMP nyi.'
  else
     stop 'wrong method.'
  end if

  if (nblock < 0) stop "error in wfn.nblock."
  do iblock = 1, nblock
     if (type_block(iblock) < -1) stop "error in wfn.type_block."
     if (nfun_block(iblock) < 0) stop "error in wfn.nfun_block."
  end do

  npair1 = 0
  nfcorex = 0
  ndcorex = 0
  ncorex = 0
  nactx = 0
  nopenx = 0
  nvira1 = 0
  nvirx1 = 0
  nfun_tot = 0
  do iblock = 1, nblock
     nfun_tot = nfun_tot + nfun_block(iblock)
     if (type_block(iblock) == -1) then
        nfcorex = nfcorex + nfun_block(iblock)
        ncorex = ncorex + nfun_block(iblock)
     else if (type_block(iblock) == 0) then
        ndcorex = ndcorex + nfun_block(iblock)
        ncorex = ncorex + nfun_block(iblock)
     else if (type_block(iblock) == 1) then
        nactx = nactx + nfun_block(iblock)
        if (imethod == 3) then
           npair1 = npair1 + 1
        end if
!nyi     else if (type_block(iblock) == 2) then
!nyi        nopenx = nopenx + nfun_block(iblock)
     else if (type_block(iblock) == 2) then
        nvira1 = nvira1 + nfun_block(iblock)
     else if (type_block(iblock) == 3) then
        nvirx1 = nvirx1 + nfun_block(iblock)
     else
        stop 'bad type_block for cas wavefunction.'
     end if
  end do
  if (nfun > 0 .and. nfun_tot /= nfun) stop 'wfn.nfun_tot .ne. wfn.nfun.'

  if (nfcore < 0) nfcore = nfcorex
  if (ndcore < 0) ndcore = ndcorex
  if (ncore < 0) ncore = ncorex
  if (nact < 0) nact = nactx
  if (nopen < 0) nopen = nopenx
  if (nocc < 0) nocc = ncore + nactx + nopenx
  if (nvira < 0) nvira = nvira1
  if (nvirx < 0) nvirx = nvirx1
  if (nvir < 0) nvir = nvira1 + nvirx1
  if (nfun < 0) nfun = nocc + nvir
  if (npair < 0) npair = npair1
  n2ov = nocc * nvir

  if (nelcore(1) < 0) nelcore(1) = ncore
  if (nelcore(2) < 0) nelcore(2) = ncore
  nelcore(3) = nelcore(1) + nelcore(2)

  if (nelact(1) < 0) nelact(1) = ne(1) - nelcore(1)
  if (nelact(2) < 0) nelact(2) = ne(2) - nelcore(2)
  nelact(3) = nelact(1) + nelact(2)


  ! check ormas parameters
  if (qcas) then
     ormas = .true.
     dpsi_nocp = .true.
     ntota = 0
     ntotb = 0
     do isub = 1, nsub
        if (nel_sub_alph(isub) < 0) stop 'error in wfn.nel_sub_alph.'
        if (nel_sub_beta(isub) < 0) stop 'error in wfn.nel_sub_alph.'
        ntota = ntota + nel_sub_alph(isub)
        ntotb = ntotb + nel_sub_beta(isub)
        min_sub(isub) = nel_sub_alph(isub) + nel_sub_beta(isub)
        max_sub(isub) = min_sub(isub)
     end do
     if (ntota /= nelact(1)) stop 'error in wfn.nel_sub_alph.'
     if (ntotb /= nelact(2)) stop 'error in wfn.nel_sub_beta.'
  else if (ormas0) then
     ormas = .true.
     dpsi_nocp = .true.
     do isub = 1, nsub
        if (min_sub(isub) /= max_sub(isub)) stop 'error in wfn.mix/max_sub for ormas0.'
     end do
  else if (rasscf) then
     ormas = .true.
     if (nsub > 3) stop 'error in wfn.nsub for rasscf.'
  end if
  !##### DONLY #####
!  if (donly) then
!     dpsi_nocp = .true.
!  end if
  !##### DONLY #####
  if (ormas) then
     s2zero = s2zero .and. nelact(1) == nelact(2)
     if (nsub == 0) stop 'error in wfn.nsub.'
     if (sum(norb_sub(1:nsub)) /= nact) stop 'error in wfn.norb_sub.'
     if (sum(norb_sub(1:nsub))*2 < nelact(3)) stop 'error in wfn.norb_sub.'
     if (sum(min_sub(1:nsub)) > nelact(3)) stop 'error in wfn.min_sub.'
     if (sum(max_sub(1:nsub)) < nelact(3)) stop 'error in wfn.max_sub.'
     do isub = 1, nsub
        if (min_sub(isub) > max_sub(isub)) stop 'error in wfn.min/max_sub.'
        if (max_sub(isub) > min(nelact(3), 2*norb_sub(isub))) stop 'error in wfn.min/max_sub.'
     end do
     lorb_sub(1, 1) = 1
     lorb_sub(2, 1) = norb_sub(1)
     do isub = 2, nsub
        lorb_sub(1, isub) = lorb_sub(2, isub-1) + 1
        lorb_sub(2, isub) = lorb_sub(2, isub-1) + norb_sub(isub)
     end do

     ! tdcc settings
     norb1 = norb_sub(1)
     norb2 = norb_sub(2)
     cc_rank = max_sub(2)
     if (optbrueckner) brueckner = .true.
     if (donly .and. brueckner) stop 'ccd and brueckner simultaneously.'
     if (donly) then
        len_tcc1 = 0
        len_gcc1 = 0
     else if (brueckner) then
        len_tcc1 = 0
        len_gcc1 = norb1 * norb2
     else
        len_tcc1 = norb1 * norb2
        len_gcc1 = norb1 * norb2
     end if
     len_tcc2 = (norb1 * norb2)**2 * 2
     len_gcc2 = (norb1 * norb2)**2 * 2
     if (cc_rank == 3) then
        len_tcc3 = (norb1*norb2)**3 * 2
        len_gcc3 = (norb1*norb2)**3 * 2
     end if
     if (cc_read_rank <= 0) cc_read_rank = cc_rank

     ind_tcc1 = 1
     ind_gcc1 = ind_tcc1 + len_tcc1
     ind_tcc2 = ind_gcc1 + len_gcc1
     ind_gcc2 = ind_tcc2 + len_tcc2
     ind_tcc3 = ind_gcc2 + len_gcc2
     ind_gcc3 = ind_tcc3 + len_tcc3
     if (occd .and. bccd) stop "choose either occd or bccd"
     if (occd .and. obccd) stop "choose either occd or obccd"
     if (occd .and. obccd2) stop "choose either occd or obccd2"
     if (bccd .and. obccd) stop "choose either bccd or obccd"
     if (bccd .and. obccd2) stop "choose either bccd or obccd2"
     if (obccd .and. obccd2) stop "choose either obccd or obccd2"

     ! orbital optimization
     noci = noci .or. moonly
     nomo = nomo .or. cionly
     nolag1 = nomo .or. nolag1 .or. ncore == 0 .or. nact == 0
     nolag2 = nomo .or. nolag2 .or. nsub == 1
     nolag3 = nomo .or. nolag3
     nomo = nomo .or. (nolag1 .and. nolag2 .and. nolag3)
  end if

  ! check ras parameters
  if (doras) then
     if (max_hole < 0) max_hole = nelact(1) + nelact(2)
     if (max_elec < 0) max_elec = nelact(1) + nelact(2)
     if (max_hole1 < 0) max_hole1 = max_hole
     if (max_elec1 < 0) max_elec1 = max_elec
     if (nact1 < 0) nact1 = 0
     if (nact3 < 0) nact3 = 0
     if (nact2 < 0) nact2 = nact - nact1 - nact3
!     doras = (max_hole > 0 .and. nact1 > 0) &
!      & .or. (max_elec > 0 .and. nact3 > 0)

     ! occ-occ orbital rotations
     nolag = nolag1 .and. nolag2
  end if

  if (nfcore .ne. nfcorex) stop "error in wfn.nfcore."
  if (ndcore .ne. ndcorex) stop "error in wfn.ndcore."
  if (ncore .ne. nfcore + ndcore) stop "error in wfn.ncore."
  if (nact .ne. nactx) stop "error in wfn.nact."
  if (nopen .ne. nopenx) stop "error in wfn.nact."
  if (nocc .ne. ncore + nact + nopen) stop "error in wfn.nocc."
  if (nvira < 0) stop "error in wfn.nvira."
  if (nvirx < 0) stop "error in wfn.nvirx."
  if (nvir < 0) stop "error in wfn.nvir."
  if (nfun .ne. nocc + nvir) stop "error in wfn.nfun."
  if (nelcore(1) .ne. ncore) stop "error in wfn.nelcore(1)."
  if (nelcore(2) .ne. ncore) stop "error in wfn.nelcore(2)."
  if (nelact(1) .ne. ne(1) - nelcore(1)) stop "error in wfn.nelact(1)."
  if (nelact(2) .ne. ne(2) - nelcore(2)) stop "error in wfn.nelact(2)."
  if (nelact(3) .ne. ne(3) - nelcore(3)) stop "error in wfn.nelact(2)."
  if (nspin /= 1 .and. nspin /= 2) stop "error in wfn.nspin."
  if (nblock < 0 .or. nblock > max_nblock) stop "error in wfn.nblock."
  if (nstate <= 0) stop "error in wfn.nstate"
  if (target_state < 0) stop "error in wfn.target_state"
  if (nstate < target_state) stop "error in wfn.nstate/target_state"
  if (.not. reg_occd) stop 'bad reg_occd in mcscf_dpsidt_dir.'
  if (dpsi_damp < zero) stop 'bad dpsi_damp in mcscf_dpsidt_dir.'
  if (dpsi_dyreg1 < zero) stop 'bad dpsi_dyreg1 in mcscf_dpsidt_dir.'
  if (dpsi_dyreg2 < zero) stop 'bad dpsi_dyreg2 in mcscf_dpsidt_dir.'

  write(6, nml=wfn)

!debug
!  stop 'stop for debug in input_wfn.'
!debug

end subroutine input_wfn
!################################################################################
