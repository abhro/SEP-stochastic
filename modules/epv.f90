module epv
  use param, only: EP, EE, cspeed
  use iso_fortran_env, only: real64
  implicit none

  ! Quantities:
  ! beta - speed in units of c (speed of light)
  ! gamma - lorentz factor
  ! v - velocity
  ! rp - ??? some sort of mass quantity (in GeV/c^2?)
  ! e -
  !
  ! rnm - ???
  ! rnz - ???


contains
  real(kind=real64) function e2p(e)
    real(kind=real64), intent(in)  :: e
    real(kind=real64)              :: rnz, rnm
    common/specie/rnz,rnm
    if (rnm > 0.5) then
      e2p = dsqrt(2.0*EP*e + e*e) * rnm / rnz
    else
      e2p = dsqrt(2.0*EE*e + e*e)
    end if
  end function

  real(kind=real64) function rp2e(rp)
    real(kind=real64), intent(in)  :: rp
    real(kind=real64)              :: rnz, rnm, rp1
    common/specie/rnz,rnm
    if (rnm > .5) then
      rp1 = rp*rnz/rnm
      rp2e = hypot(EP, rp1) - EP
    else
      rp1 = rp
      rp2e = hypot(EE, rp1) - EE
    end if
  end function

  real(kind=real64) function rp2beta(rp)
    real(kind=real64), intent(in)  :: rp
    real(kind=real64)              :: rnz, rnm, rp1
    common/specie/rnz,rnm
    if (rnm > 0.5) then
      rp1 = rp * rnz / rnm
      rp2beta = rp1 / hypot(EP, rp1)
    else
      rp1 = rp
      rp2beta = rp1 / hypot(EE, rp1)
    end if
  end function

  real(kind=real64) function rp2v(rp)
    real(kind=real64), intent(in)  :: rp
    rp2v = rp2beta(rp) * cspeed
  end function

  real(kind=real64) function e2beta(e)
    real(kind=real64), intent(in)  :: e
    real(kind=real64)              :: rp
    rp = e2p(e)
    e2beta = rp2beta(rp)
  end function

  real(kind=real64) pure function beta2gamma(beta)
    !! convert beta (velocity in units of c) to Lorentz factor gamma
    real(kind=real64), intent(in) :: beta
    beta2gamma = 1.0 / dsqrt(1.0 - beta*beta)
  end function

  real(kind=real64) function e2gamma(e)
    real(kind=real64), intent(in)  :: e
    real(kind=real64)              :: beta
    beta = e2beta(e)
    e2gamma = 1.0 / dsqrt(1.0 - beta*beta)
  end function

  real(kind=real64) function e2v(e)
    real(kind=real64), intent(in)  :: e
    e2v = e2beta(e) * cspeed
  end function

  real(kind=real64) function v2p(v)
    real(kind=real64), intent(in)  :: v
    v2p = beta2p(v / cspeed)
  end function

  real(kind=real64) function beta2p(beta)
    real(kind=real64), intent(in)  :: beta
    real(kind=real64)              :: rnz, rnm
    common/specie/rnz,rnm
    if (rnm > 0.5) then
      beta2p = beta*EP/sqrt(1.0-beta*beta)*rnm/rnz
    else
      beta2p = beta*EE/sqrt(1.0-beta*beta)
    end if
  end function

end module epv
