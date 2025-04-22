module dxx
  use iso_fortran_env, only: real64
  use param, only: CSPEED
  use epv, only: e2p, rp2beta
  use dmumu, only: cofdu, set_du0DG, set_du0BK, set_du0AH
  implicit none
contains
  subroutine preparedxx(ndxx)
    integer, intent(in) :: ndxx

    if (ndxx == 1) call read_dxx(ndxx)
  end subroutine

  subroutine cofm(r, p, pa, beta, bv, bm, cvtu, gbm, dbbds, b1s, gb1s, g, dg)
    !!   calculate diffusion coeficients in magnetic field coordinate
    !!   and derivatives of g with respect to r, theta, phi
    !!   2 perpendicular coeficients must equal to the poles
    real(kind=real64), intent(in)  :: p, pa, beta, bm, dbbds, b1s
    real(kind=real64), intent(in)  :: r(3), bv(3), cvtu(3), gbm(3), gb1s(3)
    real(kind=real64), intent(out) :: g(3), dg(3)
    real(kind=real64)              :: ggper(3), dbb(3)
    real(kind=real64)              :: ab1s, vpl, dvpl, bdotggper

    real(kind=real64)              :: g0(3), b(3), a(3)
    common/dxxcnst/g0,b,a

    integer                        :: ndxx
    common/ndxx/ndxx

    real(kind=real64)              :: densw0, vsw, k4ok2, k6ok2, omega, b1au, vom, facip
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip


    ab1s = abs(b1s)

    if (ndxx == 1) then

      !  g in magnetic field - curvature coordinates
      g(1) = 0
      !  adhoc GCR formula
      !g(2) = g0(2) * beta * p**b(2) / bm**a(2)
      !  particle following FLRW \beta*c*B(1rs)<dx^2>/(2*Vpl*B*<dt>)
      !   Vpl-plasma speed in km/s, use Leblanc et al (1998) model
      !   supergranulation dx=30000 km in dt=24hr.
      vpl = vsw / (1 + k4ok2/r(1)**2 + k6ok2/r(1)**4)
      dvpl = vpl * (2*k4ok2/r(1)**3 + 4*k6ok2/r(1)**5) / &
                (1 + k4ok2/r(1)**2 + k6ok2/r(1)**4)
      g(2) = g0(2) * beta * ab1s * CSPEED / (2 * vpl) !g0(2)=<dx^2>/<dt>
      !g(2) = g(2) * sqrt(pa*pa)
      g(3) = g(2)
      !  g tensor = gper*(1-bb)
      !  divergence of g tensor in spherical coordinates b/c of all other divergence
      !   adhoc GCR
      !ggper = -a(2) * g(2) * gbm / bm
      !   FLRW
      ggper(1) = g(2) * gb1s(1) / ab1s - g(2) * dvpl / vpl
      ggper(2) = g(2) * gb1s(2) / ab1s
      ggper(3) = g(2) * gb1s(3) / ab1s
      bdotggper = (bv(1)*ggper(1) + bv(2)*ggper(2) + bv(3)*ggper(3)) / bm
      dbb = cvtu - dbbds*bv/bm
      dg = ggper - bdotggper*bv/bm - g(2)*dbb
    else
      g(:) = 0.0
      dg(:) = 0.0
    end if
    return
  end subroutine

  subroutine read_dxx(ndxx)
    integer, intent(in) :: ndxx

    real(kind=real64)   :: g0(3), b(3), a(3)
    common/dxxcnst/g0,b,a

    real(kind=real64)   :: gperCMS, bper, aper
    namelist /inputdxx/gperCMS,bper,aper

    character(len=256)  :: dir
    common/dir/dir


    integer             :: fileunit

    if (ndxx == 1) then
      open(newunit=fileunit, file=trim(dir)//'inputdxx.nml', status='old')
      read(fileunit, nml=inputdxx)
      close(fileunit)
      g0(1) = 0.0
      !  for p=1GV B=1G,beta=1,c*rg/3=1E+17 cm^2/s
      g0(2) = gperCMS*1.2386e-20 !to rs^2/min
      g0(3) = g0(2)

      b = [1d0, bper, bper]

      a = [1d0, aper, aper]
    end if
  end subroutine

  subroutine set_rlambdax(e0)
    real(kind=real64), intent(in) :: e0
    integer                       :: nlambdax, ndxx
    common/nlambdax/nlambdax
    real(kind=real64)             :: rlambdax, rlambday
    common/rlambdax/rlambdax,rlambday
    common/ndxx/ndxx
    real(kind=real64)             :: rxl0, ryl0, g0rt
    if (ndxx == 1) then
      rxl0 = get_rlambdax0(e0)
      ryl0 = get_rlambday0(e0)
      if (nlambdax == 0) then
        rlambdax = rxl0
        rlambday = ryl0
      else
        g0rt = rxl0 / rlambdax
        !call set_dxx(g0rt)
        g0rt = ryl0 / rlambday
        !call set_dyy(g0rt)
      end if
    end if
    return
  end subroutine

  real(kind=real64) function get_rlambdax()
    real(kind=real64)  :: rlambdax, rlambday
    common/rlambdax/rlambdax,rlambday

    integer :: ndxx
    common/ndxx/ndxx

    if (ndxx == 1) then
      get_rlambdax = rlambdax
    else
      get_rlambdax = 0.0
    end if
  end

  real(kind=real64) function get_rlambda0(e0)
    real(kind=real64), intent(in)  :: e0

    integer, parameter             :: num = 20000

    real(kind=real64)              :: rp, beta, v, bm
    real(kind=real64)              :: ddd, dmu, rmu, du, ddu, dmumu0

    integer                        :: ndmumu
    common/ndmumu/ndmumu

    integer                        :: i

    rp = e2p(e0)
    beta = rp2beta(rp)
    v = beta*CSPEED
    bm = 1.0

    ddd = 0.0
    dmu = 1.0 / num
    rmu = 0.0
    do i = 1, num
      call cofdu(rp, rmu, beta, bm, du, ddu, ndmumu)
      dmumu0 = du

      if (dmumu0 /= 0) ddd = ddd + (1.0 - rmu*rmu)**2 / (2 * dmumu0)
      rmu = rmu+dmu
    end do
    get_rlambda0 = 1.5 * v * ddd * dmu
  end function

  real(kind=real64) function get_rlambdax0(e0)
    real(kind=real64), intent(in) :: e0
    real(kind=real64) :: re(3), bv(3), cvtu(3), gbm(3), gb1s(3), g(3), dg(3)
    real(kind=real64) :: rp, beta, v, bm, pa, p, dbbds, b1s
    rp = e2p(e0)
    beta = rp2beta(rp)
    v = beta * CSPEED
    bm = 1.0
    pa = 0.5

    call cofm(re, p, pa, beta, bv, bm, cvtu, gbm, dbbds, b1s, gb1s, g, dg)
    get_rlambdax0 = 3.0 * g(2) / v
  end function

  real(kind=real64) function get_rlambday0(e0)
    real(kind=real64), intent(in)  :: e0
    real(kind=real64) :: re(3), bv(3), cvtu(3), gbm(3), gb1s(3), g(3), dg(3)
    real(kind=real64) :: rp, beta, v, bm, pa
    real(kind=real64) :: p, dbbds, b1s
    rp = e2p(e0)
    beta = rp2beta(rp)
    v = beta*CSPEED
    bm = 1.0
    pa = 0.5

    call cofm(re, p, pa, beta, bv, bm, cvtu, gbm, dbbds, b1s, gb1s, g, dg)
    get_rlambday0 = 3.0 * g(3) / v
  end function

  subroutine set_rlambda(e0)
    real(kind=real64), intent(in) :: e0
    real(kind=real64)             :: du0rt, rl0

    integer                       :: ndmumu
    common /ndmumu/ndmumu

    integer                       :: nlambda
    common /nlambda/nlambda

    real(kind=real64)             :: rlambda
    common /rlambda/rlambda

    ! only du0 needed
    real(kind=real64)             :: du0, rnu, delta, va
    common /dmumuDGcnst/du0,rnu,delta,va

    rl0 = get_rlambda0(e0)
    if (nlambda == 0) then
      rlambda = rl0
    else
      du0rt = rl0/rlambda

      if (ndmumu == 0) then
        call set_du0AH(du0rt, du0)
      else if (ndmumu == 1) then
        call set_du0BK(du0rt, du0)
      else if (ndmumu == 2) then
        call set_du0DG(du0rt, du0)
      else
        stop
      end if

    end if
  end subroutine

  real(kind=real64) function get_rlambda()
    real(kind=real64) :: rlambda
    common/rlambda/rlambda
    get_rlambda = rlambda
  end function

end module dxx
