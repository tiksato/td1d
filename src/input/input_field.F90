!################################################################################
subroutine input_field(ioin)

! input field parameters

  use const_mod, only : zero, one, two, pi, au2fs
  use field_mod, only : gauge, env_type, ncyc_on, ncyc_off, ncyc_flat, sin2_ncyc, &
       & fwhm, tau, cep, fint, famp, wlen, freq, cyc0, period, rquiv, pulse_t, pulse_dt

  implicit none
  integer, intent(in) :: ioin
  integer :: ioerr
  namelist /field/ gauge, env_type, ncyc_on, ncyc_off, ncyc_flat, sin2_ncyc, &
       & fwhm, tau, cep, fint, famp, wlen, freq, cyc0, period, rquiv, pulse_t, pulse_dt

  gauge = 'L'       ! Length or Velocity gauge
  env_type = 'FLAT' ! TRAPEZOIDAL, GAUSS, SIN2, ZERO, PULSE
  ncyc_on = -1      ! (trapozoidal) number of cycles for linear turning on
  ncyc_off = -1     ! (trapozoidal) number of cycles for linear turning off
  ncyc_flat = -1    ! (trapozoidal) number of cycles for maximum amplitude
  sin2_ncyc = -1    ! (sin2) number of fundamental cycles in a sin2 window
  fwhm = -one       ! (sin2) full width half maximum of envelope in fs, then transformed to au
  tau = -one        ! (sin2) tau = 2 * fwhm / pi = half-cycle of sin2 envelop, in fs, then transformed to au
  cep = zero        ! carrier-envelope phase in degree, then transformed to radian
  fint = -one       ! field intensity (W/cm^2)
  famp = -one       ! field amplicude (au)
  wlen = -one       ! field wavelength (nm)
  freq = -one       ! field frequency (photon energy) (au)
  cyc0 = zero       ! electron's ejection time in unit of period (for 3-step model simulation)
  pulse_t = -one
  pulse_dt = -one

  rewind(ioin)
  read(unit=ioin, nml=field, iostat=ioerr)
  if(ioerr /= 0) stop "error in namelist field."
!
  call util_transchar(gauge)
  call util_transchar(env_type)

  if(wlen < zero .and. freq < zero) stop "error in field.wlen or field.freq."
  if (freq < zero) then
     freq = 1239.84190d+0 / (wlen * 27.2113845d+0)
  else
     wlen = 1239.84190d+0 / (freq * 27.2113845d+0)
  end if

  if(fint < zero .and. famp < zero) stop "error in field.fint or field.famp."
  if (famp < zero) then
     famp = sqrt(fint / 3.50944506d+16)
  else
     fint = 3.50944506d+16 * fint * fint
  end if

  period = two * pi / freq;
  rquiv = famp / (freq * freq)

  if (env_type == 'ZERO') then
  else if (env_type == 'FLAT') then
  else if (env_type == 'PULSE') then
     if (pulse_t < zero) stop "error in field.pulse_t."
     if (pulse_dt < zero) stop "error in field.pulse_dt."
  else if (env_type == 'TRAPEZOIDAL') then
     if(ncyc_on < 0) stop "error in field.ncyc_on."
     if(ncyc_off < 0) stop "error in field.ncyc_off."
     if(ncyc_flat < 0) stop "error in field.ncyc_flat."
  else if (env_type == 'SIN2') then
     if(sin2_ncyc <= 0 .and. fwhm < zero .and. tau < zero) &
          & stop "error in field.sin2_ncyc or field.fwhm or field.tau."
     if (sin2_ncyc > 0) then
        tau = two * sin2_ncyc / freq
        fwhm = tau * pi / two
     else if(fwhm < zero) then
        fwhm = tau * pi / two
        fwhm = fwhm / au2fs
        tau = tau / au2fs
     else
        tau = fwhm / pi * two
        fwhm = fwhm / au2fs
        tau = tau / au2fs
     end if
  else if (env_type == 'GAUSS') then
     stop "env_type == GAUSS, nyi."
  else
     stop "bad env_type."
  end if

  cep = cep / 180.D+0 * pi

  write(6, nml=field)

end subroutine input_field
!################################################################################
