module fb
  use iso_fortran_env, only: error_unit, real64
  use ieee_arithmetic, only: ieee_is_nan
  use param, only: CSPEED, EP, EE, PI

  implicit none

  real(kind=real64), parameter :: OMGP1 = 574745     !in 1/min
  real(kind=real64), parameter :: OMGE1 = 1055307413 !in 1/min
contains
  real(kind=real64) function fs0(tacc, xpk, bm, u1, densw, ob, amach, rsh, vth, pinj, pc)
    real(kind=real64), intent(in) :: xpk(6)
    real(kind=real64), intent(in) :: tacc, bm, u1
    !   densw in cm^-3, fs in cm^-3 (GeV/c)^-3,
    real(kind=real64), intent(in) :: densw
    real(kind=real64), intent(in) :: ob, amach, rsh
    !   vth thermal speed of protons with downstream Maxwellian temperature
    real(kind=real64), intent(in) :: vth
    real(kind=real64), intent(in out) :: pinj, pc

    !   fampb is amplification of upstream magnetic field
    !   ratkp is parallel to perpendicular diff ceof ratio
    real(kind=real64) :: ecrmax, ecrmin, ob0, obw, etinj, fampb, ratkp
    common /acc/ecrmax,ecrmin,ob0,obw,etinj,fampb,ratkp

    real(kind=real64) :: rnz, rnm
    common /specie/rnz,rnm

    real(kind=real64) :: ob2, dob, obinj, sgm, facob, pth, e1c, f0, omegae
    real(kind=real64) :: p1, omegap, p1inj

    ob2 = atan(rsh*sin(ob)/cos(ob))
    dob = amach / sqrt(40 * (1 + rsh*rsh))
    obinj = max(cos(ob), 0.25d0) !min(ob2,ob2-dob)
    etinj = 2.5 / obinj !injection at 3*u1
    !   shock slope
    !rpw = rsh * (amach - 1) / (amach + 1.0*sqrt(rsh))
    !if (rpw < 1.00001) rpw = 1.00001
    sgm = 3*rsh/(rsh-1)
    !if (sgm > 5.0) then
    !  skp = sgm/2 - 1
    !else
    !  skp = 1.5
    !end if

    !skp = 30.0

    !   obliuity factor for normal diffusion coeffiecient
    facob = sqrt(ratkp)
    facob = facob*cos(ob)**2 + sin(ob)**2/facob

    if (rnm > 0.5) then
      pth = EP * vth / CSPEED ! * sqrt((skp-1.5)/skp) ! to kappa
      !   assume all species have the same thermal speed
      !    or temperature proportional to mass (Jacco Vink et al 2015,A&A)
      pth = rnm / rnz * pth !in rigidity
      p1 = EP * u1 / CSPEED
      omegap = OMGP1 * bm * fampb
      p1inj = etinj * p1
      pinj = rnm / rnz * p1inj
      e1c = hypot(EP, p1inj) + p1*p1/EP*(1-1/rsh)*tacc*omegap/facob*rnz/rnm
      pc = sqrt(e1c*e1c - EP*EP) * rnm / rnz
    else
      ! same thermal speed as protons
      ! no thermalization between e and p
      pth = EE * vth / CSPEED

      p1 = EE * u1 / CSPEED
      omegae = OMGE1*bm*fampb
      p1inj = etinj*p1
      pinj = p1inj
      e1c = hypot(EE, p1inj) + p1*p1/EE*(1-1/rsh)*tacc*omegae/facob/100
      pc = sqrt(e1c*e1c-EE*EE)
    end if

    !ecr = (ecrmax + ecrmin - (ecrmax-ecrmin) * tanh((ob-ob0)/obw)) / 2

    !   only inject particles reached shock acceleration equilibrium
    if (pc > xpk(4)) then
      !   matching distribution at pinj
      f0 = densw*rsh/(PI*pth*pth)**1.5/exp((pinj/pth)**2) !Maxwellian
      !f0 = densw * rsh / (PI*pth*pth)**1.5 * &
      !   exp(log_gamma(skp+1) - log_gamma(skp-0.5)) / skp**1.5 / &
      !   (1 + (pinj/pth)**2 / skp)**(skp+1)
      fs0 = f0 * sgm * u1 * (1 - 1/rsh) / 3 / (xpk(4)/pinj)**sgm
      if (ieee_is_nan(fs0)) then
        !write(*,*) fs0, densw, rsh, pth, pinj, pc, sgm
        fs0 = 0.0
      end if
    else
      fs0 = 0.0
    end if

  end function


  subroutine preparefb
    integer             :: nfbconst
    real(kind=real64)   :: rb0, rmax, rk, deltat, tc, tl, tmodel0
    common/fbcnst/rb0,rmax,rk,deltat,tc,tl,tmodel0,nfbconst
    real(kind=real64)   :: sclat, sclong, scanw
    namelist/fb/rb0,rmax,rk,deltat,tc,tl,tmodel0,nfbconst,sclat,sclong,scanw
    character(len=256)  :: dir
    common/dir/dir
    real(kind=real64)   :: sp0, gp, ap, sp, h0
    real(kind=real64)   :: trgtfs(4)
    common /srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0
    integer             :: nfb

    open(newunit=nfb, file=trim(dir)//'inputfb.nml', status='old')
    read(nfb, nml=fb)
    close(nfb)
    trgtfs(1) = sin(sclat*PI/180)
    if (trgtfs(1) == 0) trgtfs(1) = 1e-27
    trgtfs(2) = cos(sclat*PI/180)
    trgtfs(3) = sin(sclong*PI/180)
    trgtfs(4) = cos(sclong*PI/180)
    scanw = scanw*PI/180
  end subroutine

  real(kind=real64) function fb0(torg, rpb)
    real(kind=real64), intent(in)  :: torg, rpb(5)
    !integer             :: nfbconst
    !real(kind=real64)              :: rb0, rmax, rk, deltat, tc, tl, tmodel0
    !common/fbcnst/rb0,rmax,rk,deltat,tc,tl,tmodel0,nfbconst
    !real(kind=real64)              :: sp0, gp, ap, sp, scanw, h0
    !real(kind=real64)              :: trgtfs(4)
    !common /srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0

    !rp = rpb(4)
    !if (torg < 0 .or. torg > deltat) then
    !  fb0 = 0.0
    !else
    !  if (nfbconst == 0) then
    !    tmodel = torg + tmodel0
    !    fbtimes = exp(-tc/tmodel - tmodel/tl) / tmodel
    !    ad = cos(rpb(2)) * trgtfs(2) + sin(rpb(2)) * trgtfs(1) * &
    !        (cos(rpb(3)) * trgtfs(4) + sin(rpb(3)) * trgtfs(3))
    !    if (ad < cos(scanw)) fbtimes = 0.0
    !  else
    !    fbtimes = 1.0
    !  end if
    !    fb0 = rp2e(rp)**rk / (rp*rp) * fbtimes
    !end if

    fb0 = 0.0

  end function


  !real(kind=real64) function getdeltat()
  !  real(kind=real64)  :: rb0, rmax, rk, deltat, tc, tl, tmodel0
  !  integer :: nfbconst
  !  common/fbcnst/rb0,rmax,rk,deltat,tc,tl,tmodel0,nfbconst
  !  getdeltat = deltat
  !end function


  !subroutine getfbcnst(rk0, deltat0, tc0, tl0, tmodel00)
  !  real(kind=real64), intent(out) :: rk0, deltat0, tc0, tl0, tmodel00
  !  integer                        :: nfbconst
  !  real(kind=real64)              :: rb0, rmax, rk, deltat, tc, tl, tmodel0
  !  common/fbcnst/rb0,rmax,rk,deltat,tc,tl,tmodel0,nfbconst
  !  rk0 = rk
  !  deltat0 = deltat
  !  tc0 = tc
  !  tl0 = tl
  !  tmodel00 = tmodel0
  !end subroutine

end module fb
