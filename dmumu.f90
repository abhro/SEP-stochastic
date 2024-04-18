module dmumu
  use iso_fortran_env, only: real64
  use param, only: PI, CSPEED

  implicit none

contains
  subroutine preparedmumu(ndmumu)
    integer, intent(in) :: ndmumu
    if (ndmumu == 0) then
      call read_dmumuAH
    else if (ndmumu == 1) then
      call read_dmumuBK
    else if (ndmumu == 2) then
      call read_dmumuDG
    else
      stop
    end if
  end subroutine

  subroutine cofdu(p, pa, beta, bm, du, ddu, ndmumu)
    !   subroutine to calculate pitch angle diffusion coefficient
    !   and its derivative with respect to pa = cos().

    real(kind=real64), intent(in)     :: p
    real(kind=real64), intent(in out) :: pa
    real(kind=real64), intent(in)     :: beta, bm
    real(kind=real64), intent(out)    :: du, ddu
    integer, intent(in)               :: ndmumu

    !real(kind=real64)                :: g0(3), b(3), a(3)

    if (ndmumu == 0) then
      du = dmumuAH(p, pa, bm, ddu)
    else if (ndmumu == 1) then
      du = dmumuBK(p, pa, beta*CSPEED, bm, ddu)
    else if (ndmumu == 2) then
      du = dmumuDG(beta*CSPEED, pa, ddu)
    else
      stop
    end if
  end subroutine

  real(kind=real64) function dmumuAH(p, pa, bm, ddu)
    real(kind=real64), intent(in)  :: p, pa, bm
    real(kind=real64), intent(out) :: ddu

    real(kind=real64)              :: du0, b1, a1
    common/dmumuAHcnst/du0,b1,a1

    dmumuAH = du0 * p**b1 / bm**a1 * (1-pa*pa)

    if (dmumuAH < 0) dmumuAH = 0.0

    ddu = -2 * du0 * p**b1 / bm**a1 * pa
  end function

  real(kind=real64) function dmumuBK(p, pa, v, bm, ddu)
    real(kind=real64), intent(in)  :: p, pa, v, bm
    real(kind=real64), intent(out) :: ddu

    real(kind=real64)              :: absmu, sgn, du00, absmuq1sigma

    real(kind=real64)              :: du0, b1, a1, qindx, sigma, h
    common /dmumuBKcnst/du0,b1,a1,qindx,sigma,h

    absmu = abs(pa)
    sgn = sign_my(pa)
    absmuq1sigma = absmu**(qindx-1.0) * (1.0 - sgn*sigma)
    !du00 = du0 * p**(qindx-2) * v / bm**(-qindx)
    du00 = du0 * p**b1 * v / bm**a1

    dmumuBK = du00 * (1-pa*pa) * (absmuq1sigma+h)
    if (dmumuBK < 0) dmumuBK = 0.0

    if (pa == 0) then
      ddu = 0.0
    else
      ddu = (-2.0*(absmuq1sigma+h)*pa+sgn*(qindx-1.0)*absmuq1sigma/&
                 absmu*(1.0-pa*pa))*du00
    end if
  end function

  real(kind=real64) pure elemental function sign_my(x)
    real(kind=real64), intent(in) :: x
    if (x > 0.0) then
      sign_my = 1.0
    else if (x < 0.0) then
      sign_my = -1.0
    else
      sign_my = 0.0
    end if
  end function

  real(kind=real64) function dmumuDG1(v, rmu, ddmumuDG1)
    real(kind=real64), intent(in)     :: v
    real(kind=real64), intent(in out) :: rmu
    real(kind=real64), intent(out)    :: ddmumuDG1

    real(kind=real64)                 :: du0, rnu, delta, va
    common/dmumuDGcnst/du0,rnu,delta,va

    real(kind=real64)                 :: vova, q, rr, theta

    vova = v/va
    q = rnu

    if (rmu >  1.0) rmu =  1.0
    if (rmu < -1.0) rmu = -1.0
    rr = hypot(delta, rmu*vova)
    theta = pi - dacos(rmu*vova/rr)

    if (rmu >= 1.0 .or. rmu <= -1.0) then
      dmumuDG1 = 0.0
    else
      dmumuDG1 = -du0 * (1.0-rmu*rmu) * rr**(q - 1.0) * dsin((1-q)*theta)
    end if
    if (dmumuDG1 < 0) then
      print *, 'mu =', rmu, ' dmumuDG1 =', dmumuDG1
      dmumuDG1 = 0.0
    end if

    if (1.0-rmu*rmu == 0.0) then
      ddmumuDG1 = 2.0 * rmu * du0 * rr ** (q - 1.0) * dsin((1-q)*theta)
    else
      ddmumuDG1 = ( &
          - (2.*rmu) / (1.-rmu*rmu) &
          - (1.-q) / (rr*rr) * rmu * vova * vova &
        ) * dmumuDG1 &
        + du0*(1.0-rmu*rmu)*rr**(-(1.0-q)) * dcos((1.0-q)*theta)*vova*delta/(rr*rr)
    end if
  end function

  real(kind=real64) function dmumuDG(v, rmu0, ddmumuDG)
    real(kind=real64), intent(in)     :: v
    real(kind=real64), intent(in out) :: rmu0
    real(kind=real64), intent(out)    :: ddmumuDG
    real(kind=real64)                 :: rmu1
    real(kind=real64)                 :: ddmumuDG1, ddmumuDG2
    real(kind=real64)                 :: d1, d2

    d1 = dmumuDG1(v, rmu0, ddmumuDG1)
    rmu1 = -rmu0
    d2 = dmumuDG1(v, rmu1, ddmumuDG2)
    rmu0 = rmu1
    dmumuDG = d1 + d2
    ddmumuDG = ddmumuDG1 - ddmumuDG2
  end function

  subroutine set_du0AH(du0rt, du0)
    real(kind=real64), intent(in) :: du0rt
    real(kind=real64), intent(in out) :: du0
    du0 = du0 * du0rt
  end subroutine

  subroutine set_du0BK(du0rt, du0)
    real(kind=real64), intent(in) :: du0rt
    real(kind=real64), intent(in out) :: du0
    du0 = du0 * du0rt
  end subroutine

  subroutine set_du0DG(du0rt, du0)
    real(kind=real64), intent(in) :: du0rt
    real(kind=real64), intent(in out) :: du0
    du0 = du0 * du0rt
  end subroutine

  subroutine read_dmumuDG
    real(kind=real64)   :: du0, rnu, delta, va
    common/dmumuDGcnst/du0,rnu,delta,va

    namelist /inputdmumuDG/du0,rnu,delta,va

    character(len=256)  :: dir
    common/dir/dir

    integer             :: fileunit

    open(newunit=fileunit, file=trim(dir)//'inputdmumu.nml', status='old')
    read(fileunit, nml=inputdmumuDG)
    close(fileunit)

    !     convert Alfven speed (km/s) to (AU/DAY)
    va = va * 5.76e-4
  end subroutine

  subroutine read_dmumuAH
    real(kind=real64)   :: du0, b1, a1
    common/dmumuAHcnst/du0,b1,a1

    namelist /inputdmumuAH/du0,b1,a1

    character(len=256)  :: dir
    common/dir/dir

    integer             :: fileunit

    open(newunit=fileunit, file=trim(dir)//'inputdmumu.nml', status='old')
    read(fileunit, nml=inputdmumuAH)
    close(fileunit)

  end subroutine

  subroutine read_dmumuBK
    integer             :: nlambda

    real(kind=real64)   :: du0, b1, a1, qindx, sigma, h, rlambda
    common/dmumuBKcnst/du0,b1,a1,qindx,sigma,h

    namelist /inputdmumuBK/du0,rlambda,nlambda,b1,a1,qindx,sigma,h

    character(len=256)  :: dir
    common/dir/dir

    integer             :: fileunit

    open(newunit=fileunit, file=trim(dir)//'inputdmumu.nml', status='old')
    read(fileunit, nml=inputdmumuBK)
    close(fileunit)
  end subroutine
end module dmumu
