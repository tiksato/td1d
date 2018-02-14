program laser

  implicit none
  character(len=256) :: argc
  real(kind(0d0)) :: int, wlen, ip
  real(kind(0d0)) :: amp, ene, period, rquiv, epon, keldish
  real(kind(0d0)) :: ene_ev, period_fs, rquiv_ang, epon_ev

  real(kind(0d0)), parameter :: two = 2.0d+0
  real(kind(0d0)), parameter :: four = 4.0d+0
  real(kind(0d0)), parameter :: pi = 3.1415926535d+0
  real(kind(0d0)), parameter :: int2amp = 3.50944506d+16
  real(kind(0d0)), parameter :: lam2ene = 45.5633524d+0
  real(kind(0d0)), parameter :: toev = 27.2113845d+0
  real(kind(0d0)), parameter :: tofs = 0.0241889d+0
  real(kind(0d0)), parameter :: toang = 0.5291772086d+0

  call getarg(1, argc); read(argc, *) int
  call getarg(2, argc); read(argc, *) wlen
  call getarg(3, argc); read(argc, *) ip

  amp = sqrt(int / int2amp)
  ene = lam2ene / wlen
  period = two * pi /ene
  rquiv = amp / (ene * ene)
  epon = (amp * amp) / (four * ene * ene)
  keldish = sqrt(ip / toev / two / epon)

  ene_ev = ene * toev
  period_fs = period * tofs
  rquiv_ang = rquiv * toang
  epon_ev = epon * toev

  write(6, "(' field amplitide         :', f20.10, ' au')") amp
  write(6, "(' photone energy          :', f20.10, ' au ', f20.10, ' eV')") ene, ene_ev
  write(6, "(' period                  :', f20.10, ' au ', f20.10, ' fs')") period, period_fs
  write(6, "(' quiver-motion amplitide :', f20.10, ' au ', f20.10, ' ang')") rquiv, rquiv_ang
  write(6, "(' ponderomotive energy    :', f20.10, ' au ', f20.10, ' eV')") epon, epon_ev
  write(6, "(' keldish parameter       :', f20.10)") keldish

end program laser
