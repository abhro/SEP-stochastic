module sim3d_utils
  use iso_fortran_env, only: real64

  use param, only: PI, TWOPI, GAMMA_CS, NSPMAX, nseedmax, CSPEED, &
                   N_R, N_THETA, N_PHI, bgrid, gbgrid, b1rsgrid
  use fb, only: fb0
  use epv, only: rp2beta
  use mtrx, only: mrtx, mbtr, trilinear, trilineardif

  implicit none

contains

  subroutine f0mod(r, pa, f0, df0, ddf0, df0dmu, ddf0dmu2)

    !  function modification factor for distribution
    !    f=f0mod*g -- g is the new function to be solve the stochastic

    real(kind=real64), intent(in)  :: r(3), pa
    real(kind=real64), intent(out) :: f0, df0(3), ddf0(3,3), df0dmu, ddf0dmu2

    real(kind=real64)  :: rs(3), sintheta, costheta, sinphi, cosphi

    real(kind=real64)  :: sp, sp0, gp, ap, trgtfs(4), scanw, h0
    common /srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0

    real(kind=real64)  :: densw0, vsw, k4ok2, k6ok2, omega, b1au, vom, facip
    common /bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip

    real(kind=real64)  :: rb0, rmax, rk, deltat, tc, tl, tmodel0
    integer :: nfbconst
    common /fbcnst/rb0,rmax,rk,deltat,tc,tl,tmodel0,nfbconst

    real(kind=real64)  :: dcddt, dcddp, ddcddtt, ddcddtp, ddcddpp
    real(kind=real64)  :: daddt, daddp, ddaddtt, ddaddtp, ddaddpp
    real(kind=real64)  :: ad, cd, sd, sdp

    rs(1) = rb0
    rs(2) = r(2)
    rs(3) = r(3) + (r(1) - rs(1)) / vom
    sintheta = sin(rs(2))
    if (sintheta == 0) sintheta = 1e-27
    costheta = cos(rs(2))
    sinphi = sin(rs(3))
    cosphi = cos(rs(3))

    cd = costheta * trgtfs(2) + sintheta * trgtfs(1) * (cosphi * trgtfs(4) + sinphi * trgtfs(3))
    dcddt = -sintheta * trgtfs(2) + costheta * trgtfs(1) * (cosphi * trgtfs(4) + sinphi * trgtfs(3))
    dcddp = -sintheta * trgtfs(1) * (sinphi * trgtfs(4) - cosphi * trgtfs(3))
    ddcddtt = -cd
    ddcddtp = -costheta * trgtfs(1) * (sinphi * trgtfs(4) - cosphi * trgtfs(3))
    ddcddpp = -sintheta * trgtfs(1) * (cosphi * trgtfs(4) + sinphi * trgtfs(3))
    ad = acos(cd)
    sd = sin(ad)
    daddt = -1 / sd * dcddt
    daddp = -1 / sd * dcddp
    ddaddtt = -cd / sd * daddt * daddt - 1 / sd * ddcddtt
    ddaddtp = -cd / sd * daddp * daddt - 1 / sd * ddcddtp
    ddaddpp = -cd / sd * daddp * daddp - 1 / sd * ddcddpp

    !f0 = (1+pa/sp) * exp(-ad*ad/ap/ap/2)
    f0 = log(1+pa/sp) - ad*ad/ap/ap/2   !used as exp(f0)

    ! \frac{df_0}{d\mu}
    df0dmu = 1/(sp+pa)   !normalized by f0
    ! \frac{d^2 f_0}{d\mu^2}
    ddf0dmu2 = 0 !gp * (gp-1) / (sp + pa)**2

    !  df0(1-3) = derivative in sphere coordinates /f0
    !    first order
    if (sd == 0) then
      df0(2) = 0.0
      df0(3) = 0.0
    else
      df0(2) = -ad / ap / ap * daddt
      df0(3) = -ad / ap / ap * daddp
    end if
    df0(1) = df0(3)/vom
    !    second order
    if (sd == 0) then
      if (ad == 0) then
        ddf0(2,2) = df0(2) * df0(2) - 1 / ap / ap
        ddf0(2,3) = df0(2) * df0(3) + 0
        ddf0(3,3) = df0(3) * df0(3)
      else
        sdp = 1e-10
        ddf0(2,2) = df0(2) * df0(2) + PI / ap / ap / sdp
        ddf0(2,3) = df0(2) * df0(3) + 0
        ddf0(3,3) = df0(3) * df0(3) + PI / ap / ap / sdp * trgtfs(1)**2
      end if
    else
      !ddsd = (sd-ad*cd)/sd/sd/ap/ap
      !ddf0(2,2) = df0(2)*df0(2) - ad*cd/ap/ap/sd + ddsd*dcddt*dcddt
      !ddf0(2,3) = df0(2)*df0(3) + ad/ap/ap/sd*ddcddtp + ddsd*dcddt*dcddp
      !ddf0(3,3) = df0(3)*df0(3) + ad/ap/ap/sd*ddcddpp + ddsd*dcddp*dcddp
      ddf0(2,2) = -ad/ap/ap*ddaddtt + df0(2)*df0(2)*(1-ap*ap/ad/ad)
      ddf0(2,3) = -ad/ap/ap*ddaddtp + df0(2)*df0(3)*(1-ap*ap/ad/ad)
      ddf0(3,3) = -ad/ap/ap*ddaddpp + df0(3)*df0(3)*(1-ap*ap/ad/ad)
    end if
    !ddf0(1,1) = ddf0(3,3)/vom**2
    ddf0(1,1) = -(ad*ddaddpp+daddp*daddp)/ap/ap/vom/vom + df0(1)**2
    !ddf0(1,2) = ddf0(2,3)/vom
    ddf0(1,2) = -(ad*ddaddtp+daddt*daddp)/ap/ap/vom + df0(1)*df0(2)
    !ddf0(1,3) = ddf0(3,3)/vom
    ddf0(1,3) = -(ad*ddaddpp+daddp*daddp)/ap/ap/vom + df0(1)*df0(3)
    !   gradient-f0 /f0 in spherical
    df0(2) = df0(2) / r(1)
    df0(3) = df0(3) / (r(1)*sintheta)
    !   grdaient-gradient-f0 /f0
    ddf0(2,1) = ddf0(1,2)/r(1) - df0(2)/r(1)/r(1)
    ddf0(3,1) = (ddf0(1,3) - df0(3)/r(1))/r(1)/sintheta
    ddf0(2,2) = (ddf0(2,2)/r(1) + df0(1))/r(1)
    ddf0(3,2) = (ddf0(2,3)-costheta/sintheta*df0(3))/r(1)/r(1)/sintheta
    ddf0(3,3) = (ddf0(3,3)/r(1)/sintheta + df0(1)*sintheta + &
      df0(2)/r(1)*costheta)/r(1)/sintheta
    ddf0(1,2) = ddf0(2,1)
    ddf0(1,3) = ddf0(3,1)
    ddf0(2,3) = ddf0(3,2)

  end subroutine


  real(kind=real64) function solarwindtemp(r) result(temp)
    real(kind=real64), intent(in)  :: r(3)
    !  empirical model in the corona from Withbroe (ApJ 325,442,1988) 10.1086/166015

    if (r(1) > 10.0) temp = 1.4d6 * (10/r(1))**1.3333

    if (r(1) <= 10.0 .and. r(1) >= 1.125) temp = 1.4d6

    if (r(1) < 1.125) then
      temp = (1.0d5**3.5 + 25.33965d6 * (1/1.0287 - 1/r(1))) ** (1/3.5)
    end if

  end function


  !   use bisection to search for root of shock adiabatic equation
  !   calculate shock compression ratio of oblique MHD shock
  function compress(amach, smach, ob) result(rsh)
    real(kind=real64), intent(in) :: amach, smach, ob
    real(kind=real64)             :: rsh(3)
    !real(kind=real64)             :: xr, y1, y2, xmin, xmax, froot
    real(kind=real64)             ::  rsh0
    real(kind=real64)             :: sintheta2 !, costheta2, tbn
    real(kind=real64)             :: pbeta
    real(kind=real64)             :: a0, a1, a2, a3, b0, b1, del
    complex(kind=real64)          :: wth, pr3, root(3)
    real(kind=real64)             :: amach2, amach4, amach6
    integer                       :: k

    pbeta = (amach/smach)**2
    sintheta2 = sin(ob)**2
    amach2 = amach * amach
    amach4 = amach2 * amach2
    amach6 = amach4 * amach2

    a3 = -2*pbeta + pbeta*sintheta2*4 - pbeta*sintheta2*sintheta2*2 - amach2*sintheta2&
          - amach2*GAMMA_CS + amach2 + amach2*GAMMA_CS*sintheta2
    a2 = -amach2 * (&
      -GAMMA_CS - pbeta*4 + GAMMA_CS*sintheta2 + pbeta*sintheta2*4 + sintheta2&
      -amach2*GAMMA_CS*2 + amach2*2 + amach2*GAMMA_CS*sintheta2 - 1)
    a1 = -amach4 * (&
      GAMMA_CS*2 + pbeta*2 - GAMMA_CS*sintheta2 - sintheta2*2 + amach2*GAMMA_CS - amach2 + 2)
    a0 = amach6 * GAMMA_CS + amach6

    rsh0 = -a2 / (3*a3)
    b1 = (3*a3*a1-a2*a2) / (3*a3*a3)
    b0 = (2*a2*a2*a2-9*a3*a2*a1+27*a3*a3*a0) / (27*a3*a3*a3)
    del = 4*b1*b1*b1 + 27*b0*b0


    if (del < 0) then
      pr3 = (-b1/3)**0.5
      wth = 3*b0/(2*b1)*(-3/b1)**0.5
      wth = acos(wth) / 3
      do k = 0, 2
        root(k+1) = 2*pr3*cos(wth-2*PI*k/3) + rsh0
        if (abs(aimag(root(k+1))) < 1.0d-10) then
          rsh(k+1) = real(root(k+1))
        else
          rsh(k+1) = 1
        end if
      end do
    else
      if (b1 < 0) then
        wth = -3*abs(b0)/(2*b1)*(-3/b1)**0.5
        wth = acosh(wth)/3
        !wth = log(wth + sqrt(wth + 1) * sqrt(wth - 1)) / 3 ! this is acosh
        root(1) = -2*abs(b0)/b0*(-b1/3)**0.5*cosh(wth) + rsh0
      else
        wth = 3 * b0/(2*b1)*(3/b1)**0.5
        wth = asinh(wth)/3
        !wth = log(wth + sqrt(wth * wth + 1)) / 3 ! this is asinh
        root(1) = -2*(b1/3)**0.5*sinh(wth) + rsh0
      end if
      root(2) = 1
      root(3) = 1
      rsh(:) = real(root(:))
    end if
  end function

  recursive subroutine split(&
      iseed, rpb, ck, fs, t, nsplvl, dnsk, bv0, flx, dflx, walk3d, nsts)

    integer, intent(in)               :: iseed
    real(kind=real64), intent(in)     :: rpb(5)
    real(kind=real64), intent(in out) :: ck, fs
    real(kind=real64), intent(in)     :: t
    integer, intent(in out)           :: nsplvl
    real(kind=real64), intent(in out) :: dnsk, bv0(3), flx, dflx
    external                          :: walk3d
    integer, intent(in)               :: nsts

    integer           :: nodr(NSPMAX)
    real(kind=real64) :: t0sv(2**(NSPMAX+1)), cksv(2**(NSPMAX+1))
    real(kind=real64) :: rpbsv(5,2**(NSPMAX+1))
    real(kind=real64) :: fssv(2**(NSPMAX+1))
    real(kind=real64) :: dnsksv(2**(NSPMAX+1))
    real(kind=real64) :: bv0sv(3,2**(NSPMAX+1))
    common /svsp/nodr,t0sv,cksv,rpbsv,fssv,dnsksv,bv0sv

    real(kind=real64) :: rps0(5)
    !real(kind=real64) :: rp0org(5)
    real(kind=real64) :: rb(3)
    real(kind=real64) :: fb_, f1

    real(kind=real64)  :: t0org, te, tdl, dmapjul, tcme0
    common /tmprm/t0org,te,tdl,dmapjul,tcme0

    real(kind=real64)  :: rb0, rmax
    common /rb0max/rb0,rmax

    integer :: nlambdaconst
    common /nlambdaconst/nlambdaconst

    real(kind=real64) :: densw0, vsw, k4ok2, k6ok2, omega, b1au, vom, facip
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip
    !common /param5/p1,beta1,pa1,pol1

    integer :: ndpdt
    common/ndpdt/ndpdt

    !real(kind=real64) :: trgtf(4)
    real(kind=real64) :: sp, sp0, gp, ap
    real(kind=real64) :: trgtfs(4)
    real(kind=real64) :: scanw, h0
    common /srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0
    real(kind=real64) :: df0(3),ddf0(3,3)
    real(kind=real64) :: t0, tsp, hb, df0dmu, ddf0dmu2, pab, rate, tb
    integer           :: nsv, i, ns
    character(len=*), parameter :: writefmt = "(i1,i3,8(1pe13.5),i3)"

    nsplvl = nsplvl + 1

    !  save data for splitting

    nsv = 2**nsplvl - 2
    ! nsv = 0
    ! do i = 1, nsplvl-1
    !   nsv = nsv + 2**i
    ! end do
    do i = 1, nsplvl-1
      nsv = nsv + 2*2**i*nodr(i)
    end do
    t0sv(nsv+1) = t
    cksv(nsv+1) = ck
    fssv(nsv+1) = fs
    rpbsv(:, nsv+1) = rpb(:)
    dnsksv(nsv+1) = dnsk
    bv0sv(:, nsv+1) = bv0(:)
    t0sv(nsv+2) = t
    cksv(nsv+2) = ck
    fssv(nsv+2) = fs
    rpbsv(:, nsv+2) = rpb(:)
    dnsksv(nsv+2) = dnsk
    bv0sv(:, nsv+2) = bv0(:)

    do i = 1,2
      nodr(nsplvl) = i-1
      t0 = t0sv(nsv+i)
      ck = cksv(nsv+i)
      fs = fssv(nsv+i)
      rps0(1:5) = rpbsv(1:5,nsv+i)
      dnsk = dnsksv(nsv+i)
      bv0(1:3) = bv0sv(1:3,nsv+i)

      tsp = (nsplvl+1)*tdl
      if (tsp > te) tsp = te
      call walk3d(iseed, rps0, rpb, ck, fs, t0, t, tsp, ns, dnsk, bv0, nsplvl)
      if (ns == -1 .and. tsp < te) then
        call split(iseed, rpb, ck, fs, t, nsplvl, dnsk, bv0, flx, dflx, walk3d, nsts)
      else
        rb(1:3) = rpb(1:3)
        pab = rpb(5)
        call f0mod(rb, pab, hb, df0, ddf0, df0dmu, ddf0dmu2)
        rate = exp(ck - hb + h0)
        if (ns >= 0) ns = 1
        fb_ = fb0(tb, rpb) * rate
        f1 = (fs+fb_) / 2**nsplvl
        !if (fs > 0.) write(nsts,writefmt) 0,nsplvl,fs
        if (fs > 1.d-30) write(nsts,writefmt) 0,nsplvl,fs
        !if (fb_ > 0.) write(nsts,writefmt) 2,nsplvl,fb_,t0org-t,rate,rpb,ns
        !call flush(nsts)
        flush(nsts)
        !OMP CRITICAL sum
        flx = flx + f1
        dflx = dflx + f1*f1
        !OMP END CRITICAL sum
      end if
    end do
    nsplvl = nsplvl - 1
  end subroutine

  subroutine vfunc(t, xpk, dxpkdt, du, gxw2, gxw3, bv0, densw, vpl, gper, b1s)
    real(kind=real64), intent(in) :: t, xpk(6)
    real(kind=real64), intent(out) :: dxpkdt(6)
    real(kind=real64)             :: du
    real(kind=real64), intent(out) :: gxw2(3), gxw3(3)
    real(kind=real64)             :: densw
    real(kind=real64), intent(out) :: bv0(3)
    real(kind=real64)             :: vpl(3)
    real(kind=real64), intent(out) :: gper
    real(kind=real64)             :: b1s

    real(kind=real64)             :: rnz, rnm
    common/specie/rnz,rnm

    real(kind=real64)             :: densw0, vsw, k4ok2, k6ok2, omega, b1au, vom, facip
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip

    real(kind=real64)             :: r(3), b2r(3,3), r2x(3,3), b2x(3,3)
    real(kind=real64)             :: bv(3), vd(3)
    real(kind=real64)             :: avr(3), avx(3)
    real(kind=real64)             :: cvtu(3), gbmag(3), bxgb2(3), pol, gb1s(3)
    real(kind=real64)             :: culper(3), culpar(3)
    real(kind=real64)             :: uax1(3), uax2(3), uax3(3)
    real(kind=real64)             :: gb(3), gr(3,3), dgr(3)!, dgx(3)
    real(kind=real64)             :: gvpl(3,3)
    real(kind=real64)             :: df0(3), ddf0(3,3)
    real(kind=real64), parameter  :: rg1 = 4.79259e-5
    real(kind=real64)             :: trgtfs(4)
    real(kind=real64)             :: sp, sp0, gp, ap, scanw, h0
    common /srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0
    real(kind=real64)             :: vx(3)
    integer                       :: lpb!, lpb2
    real(kind=real64)             :: p, pa
    real(kind=real64)             :: beta, sinmu, cosp
    real(kind=real64)             :: sintheta, costheta, sinphi, cosphi
    real(kind=real64)             :: adpadt, bm, um, rgmu2, bbgv, dcs, f0, dbbds
    real(kind=real64)             :: ddu, divv, dpadt, fk, vpcl, df0dmu, ddf0dmu2

    integer                       :: ndmumu
    common/ndmumu/ndmumu

    integer                       :: i, j

    r(1) = norm2(xpk(1:3))
    r(2) = acos(xpk(3)/r(1))
    r(3) = atan2(xpk(2),xpk(1))
    p = xpk(4)
    beta = rp2beta(p)
    pa = xpk(5)

    sintheta = sin(r(2))
    if (sintheta == 0) sintheta = 1e-27
    costheta = cos(r(2))
    sinphi = sin(r(3))
    cosphi = cos(r(3))

    call drvbmag(r, bv, bm, cvtu, gbmag, bxgb2, dbbds, pol, b1s, gb1s)
    bv0 = bv
    lpb = pol
    !  local coordinates
    !   unit outward field line direction: uax1
    uax1 = lpb*bv/bm
    !   unit perdicular directions
    uax2(1:3) = [0d0, -lpb*bv(3), lpb*bv(2)]
    um = hypot(uax2(2), uax2(3))
    if (um == 0) then
      uax2(2) = 1
      uax2(3) = 0
    else
      uax2(2) = uax2(2)/um
      uax2(3) = uax2(3)/um
    end if
    ! take cross product and normalize
    uax3(1) = uax1(2)*uax2(3) - uax1(3)*uax2(2)
    uax3(2) = uax1(3)*uax2(1) - uax1(1)*uax2(3)
    uax3(3) = uax1(1)*uax2(2) - uax1(2)*uax2(1)
    uax3 = uax3 / norm2(uax3(1:3))


    !   focusing along outward field line
    dbbds = lpb*dbbds

    !  calculate matrix from polar spheric to xyz coordinates
    r2x = mrtx(sintheta, costheta, sinphi, cosphi)
    !  calculate matrix from magnetic to polar spheric coordinates
    b2r = mbtr(uax1, uax2, uax3)
    !  calculate matrix from magnetic to xyz coordinates
    b2x = matmul(r2x, b2r)

    !  diffusion tensor in magnetic coordinates
    call cofm(r, p, pa, beta, bv, bm, cvtu, gbmag, dbbds, b1s, gb1s, gb, dgr)
    !  calculate diffusion coeficients
    !     in polar spheric coordinates
    gper = gb(2) !gb(2)=gb(3)
    do i = 1, 3
      do j = 1, 3
        gr(i,j) = sum(b2r(i,:) * gb(:) * b2r(j,:))
      end do
    end do

    !   1. parallel particle speed + solar wind speed
    vpcl = beta * CSPEED
    !   2. drift velocity in spherical coordinates
    dcs = 1.0d0 / norm2(bxgb2)
    sinmu = 1 - pa*pa
    if (sinmu < 0) then
      sinmu = 0
    else
      sinmu = sqrt(sinmu)
    end if
    rgmu2 = rg1 * p / bm * sinmu * 1.414 !sqrt(2),linear; 2,step
    if (dcs > rgmu2) then !regular drift
      culpar = dot_product(bv(1:3), cvtu(1:3))/bm/bm*bv
      culper = cvtu - culpar
      vd = rg1 * p * vpcl / bm * (pa*pa*culper/bm &
                                  + (1-pa*pa)/2/bm*culpar + (1+pa*pa)/2*bxgb2)
    else ! current sheet drift (square delta function)
      vd = -vpcl * sinmu / 2.828 * bxgb2 * dcs
    end if
    if (rnz < 0) vd = -vd
    !   3. solar wind speed
    call solarwind(r, vpl, gvpl, densw)

    !   4 artificial drift
    call f0mod(r, pa, f0, df0, ddf0, df0dmu, ddf0dmu2)
    avr(1:3) = 2 * matmul(gr, df0)
    avx(1:3) = matmul(r2x, avr)

    !   total
    do i = 1,3
      vx(i) = dot_product(r2x(i,1:3), vpl)
      dxpkdt(i) = -b2x(i,1)*vpcl*pa - vx(i) - dot_product(r2x(i,1:3), vd)
      !  add dg and artifical drift
      dxpkdt(i) = dxpkdt(i) + avx(i) + dot_product(r2x(i,1:3), dgr)
    end do

    !  get pitch angle diffusion coefficient
    call cofdu(p, pa, beta, bm, du, ddu, ndmumu)
    if (r(1) > 20.0) then
      cosp = bv(1) / bm
      du = du * cosp * cosp / facip
      ddu = ddu * cosp * cosp / facip
    end if
    !   pitch angle drift term
    !   focusing
    dpadt = -vpcl * dbbds * (1 - pa*pa) / 2.0
    !   cooling
    bbgv = 0.0
    do i = 1, 3
      do j = 1, 3
        bbgv = bv(i)/bm*gvpl(i,j)*bv(j)/bm + bbgv
      end do
    end do
    divv = gvpl(1,1)+gvpl(2,2)+gvpl(3,3)
    !divv1 = 2*vpl(1)/r(1) + gvpl(1,1)  ! for test
    dpadt = dpadt + pa * (1-pa*pa) / 2 * (divv-3*bbgv)
    !   derivative of pitch-angle diffusion
    !   artificial due to modification of function
    adpadt = 2.0 * du * df0dmu
    dxpkdt(5) = -dpadt + ddu + adpadt

    !   momentum loss
    dxpkdt(4) = ((1-pa*pa)/2*(divv-bbgv) + pa*pa*bbgv) * p

    !   killing factor
    fk = (-dpadt+ddu)*df0dmu + du*ddf0dmu2
    ! increment fk by trace(matmul(gr,  ddf0))
    fk = fk + ( &
      gr(1,1)*ddf0(1,1) + gr(1,2)*ddf0(2,1) + gr(1,3)*ddf0(3,1) + &
      gr(2,1)*ddf0(1,2) + gr(2,2)*ddf0(2,2) + gr(2,3)*ddf0(3,2) + &
      gr(3,1)*ddf0(1,3) + gr(3,2)*ddf0(2,3) + gr(3,3)*ddf0(3,3) )
    fk = fk + dot_product(dgr, df0) &
            - dot_product(vd,  df0) &
            - dot_product(vpl, df0)
    dxpkdt(6) = fk

    !  diffusion vector applied to dw2 and dw3
    gxw2(:) = b2x(1:3,2) * sqrt(2*gb(2))
    gxw3(:) = b2x(1:3,3) * sqrt(2*gb(3))

    !write(54,"(f14.7,13(1pe14.5))") t, xpk, dxpkdt,bm
    !call flush(54)
  end subroutine

  subroutine solarwind(r, vpl, gvpl, densw)
    !!   calculate solar wind velocity in the corotating frame
    !!     and its gradient
    !!    in spherical coordinate system
    real(kind=real64), intent(in)  :: r(3)
    real(kind=real64), intent(out) :: vpl(3), gvpl(3,3), densw
    !     use leBalnc(1998) model
    real(kind=real64)              :: densw0, vsw, k4ok2, k6ok2
    real(kind=real64)              :: omega, b1au, vom, facip
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip

    real(kind=real64)              :: sintheta, sinthetap, costheta, rpl

    sintheta = sin(r(2))
    if (sintheta == 0) then
      sinthetap = 1e-20
    else
      sinthetap = sintheta
    end if
    costheta = cos(r(2))
    rpl = 1 + k4ok2/r(1)**2 + k6ok2/r(1)**4

    if (r(1) < 2.5) then
      vpl = 0.0
      gvpl = 0.0
    else
      vpl(1) = vsw/rpl
      vpl(2) = 0
      vpl(3) = -omega*(r(1)-2.5)*sintheta !rss=2.5Rs for pfss

      gvpl(1,1) = vpl(1)*(2*k4ok2/r(1)**3+4*k6ok2/r(1)**5)/&
        (1+k4ok2/r(1)**2+k6ok2/r(1)**4)
      gvpl(1,2) = 0.d0
      gvpl(1,3) = -omega*sintheta
      gvpl(2,1) = 0.0
      gvpl(2,2) = vpl(1)/r(1)
      gvpl(2,3) = -omega*(r(1)-2.5)*costheta/r(1)
      gvpl(3,1) = -vpl(3)/r(1)
      gvpl(3,2) = -costheta/sinthetap*vpl(3)/r(1)
      gvpl(3,3) = vpl(1)/r(1)
    end if
    densw = densw0 * rpl / r(1)**2
  end subroutine

  function solarwind1(r) result(vpl)
    !   calculate solar wind velocity in the corotating frame
    !     and its gradient in sphereical coordinate system
    real(kind=real64), intent(in) :: r(3)
    real(kind=real64)             :: vpl(3)
    !     use leBalnc(1998) model
    real(kind=real64)             :: densw0, vsw, k4ok2, k6ok2
    real(kind=real64)             :: omega, b1au, vom, facip
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip

    if (r(1) < 2.5) then
      vpl = 0.0
    else
      vpl(1) = vsw / (1 + k4ok2/r(1)**2 + k6ok2/r(1)**4)
      vpl(2) = 0.d0
      vpl(3) = -omega * (r(1)-2.5) * sin(r(2)) !rss=2.5Rs for pfss
    end if
  end function

  subroutine drvbmag(r1, b, bmag, cvtu, gbmag, bxgb2, dbbds, pol, b1rs, gb1rs)
    real(kind=real64), intent(in)  :: r1(3)
    real(kind=real64), intent(out) :: b(3), bmag, cvtu(3), gbmag(3), bxgb2(3), dbbds, pol
    !  cvtu is now curl B
    real(kind=real64), intent(out) :: b1rs, gb1rs(3)
    real(kind=real64), parameter   :: rss = 2.5d0
    real(kind=real64)              :: r(3)
    real(kind=real64)              :: phic(2,2,2), px(3), phi
    real(kind=real64)              :: df(3)
    !real(kind=real64)  :: b1rsgrid(0:N_R, 0:N_THETA, 0:N_PHI)
    integer                        :: irr, itheta, iphi, m
    integer                        :: i1, i2, i3
    real(kind=real64)              :: bx, by, bz, gbx, gby, gbz, bm
    real(kind=real64)              :: vr, vt, tpsw
    real(kind=real64)              :: densw0, vsw, k4ok2, k6ok2, omega, b1au, vom, facip
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip
    real(kind=real64)              :: rsovr, plbs, tanp, oneptan2, sroneptan2
    real(kind=real64)              :: vpl(3), dvr1, dvf1, dvf2
    real(kind=real64)              :: gbmags(3), bxgb2s(3)
    real(kind=real64)              :: bplus(3), bminus(3)
    integer                        :: nplus, nminus
    integer                        :: m1, m2, m3
    real(kind=real64)              :: sintheta, costheta

    r = r1
    if (r(1) > rss) then
      tpsw = ((r(1) - k4ok2/r(1) - k6ok2/r(1)**3/3) &
        - (rss - k4ok2/rss - k6ok2/rss**3/3) ) / vsw
    else
      tpsw = 0.0
    end if
    r(2) = acos(cos(r(2)))
    r(3) = r(3) + omega*tpsw
    r(3) = atan2(sin(r(3)), cos(r(3)))
    if (r(3) < 0.0) r(3) = r(3) + 2*pi
    if (r(1) > rss) r(1) = rss
    if (r(1) < 1.0) r(1) = 1.0
    sintheta = sin(r(2))
    if (sintheta == 0.0) sintheta = 1e-20

    !  find the grid cell
    irr = floor((r(1)-1) / (rss-1) * N_R)
    if (irr >= N_R) irr = N_R - 1

    itheta = floor(r(2) / pi * N_THETA)
    if (itheta >= N_THETA) itheta = N_THETA - 1

    iphi = floor(r(3) / pi / 2 * N_PHI)
    if (iphi >= N_PHI) iphi = N_PHI-1

    !  relative displacement from lower grids
    px(1) = (r(1)-1.0) / (rss-1.0) * N_R - irr
    px(2) = r(2)/pi* N_THETA - itheta
    px(3) = r(3)/TWOPI* N_PHI - iphi

    do m = 1, 3
      phic(1,1,1) = bgrid(irr,   itheta,   iphi,   m)
      phic(2,1,1) = bgrid(irr+1, itheta,   iphi,   m)
      phic(1,2,1) = bgrid(irr,   itheta+1, iphi,   m)
      phic(2,2,1) = bgrid(irr+1, itheta+1, iphi,   m)
      phic(1,1,2) = bgrid(irr,   itheta,   iphi+1, m)
      phic(2,1,2) = bgrid(irr+1, itheta,   iphi+1, m)
      phic(1,2,2) = bgrid(irr,   itheta+1, iphi+1, m)
      phic(2,2,2) = bgrid(irr+1, itheta+1, iphi+1, m)
      phi = trilinear(phic, px)
      b(m) = phi
      ! phic(1,1,1) = cvgrid(irr,     itheta,     iphi,     m)
      ! phic(2,1,1) = cvgrid(irr + 1, itheta,     iphi,     m)
      ! phic(1,2,1) = cvgrid(irr,     itheta + 1, iphi,     m)
      ! phic(2,2,1) = cvgrid(irr + 1, itheta + 1, iphi,     m)
      ! phic(1,1,2) = cvgrid(irr,     itheta,     iphi + 1, m)
      ! phic(2,1,2) = cvgrid(irr + 1, itheta,     iphi + 1, m)
      ! phic(1,2,2) = cvgrid(irr,     itheta + 1, iphi + 1, m)
      ! phic(2,2,2) = cvgrid(irr + 1, itheta + 1, iphi + 1, m)
      ! phi = trilinear(phic, px)
      ! cvtu(m) = phi
      phic(1,1,1) = gbgrid(irr,   itheta,   iphi,   m)
      phic(2,1,1) = gbgrid(irr+1, itheta,   iphi,   m)
      phic(1,2,1) = gbgrid(irr,   itheta+1, iphi,   m)
      phic(2,2,1) = gbgrid(irr+1, itheta+1, iphi,   m)
      phic(1,1,2) = gbgrid(irr,   itheta,   iphi+1, m)
      phic(2,1,2) = gbgrid(irr+1, itheta,   iphi+1, m)
      phic(1,2,2) = gbgrid(irr,   itheta+1, iphi+1, m)
      phic(2,2,2) = gbgrid(irr+1, itheta+1, iphi+1, m)
      phi = trilinear(phic, px)
      gbmag(m) = phi
    end do


    do i1 = 0, 1
      do i2 = 0, 1
        do i3 = 0, 1
          bx  =  bgrid(irr+i1, itheta+i2, iphi+i3, 1)
          by  =  bgrid(irr+i1, itheta+i2, iphi+i3, 2)
          bz  =  bgrid(irr+i1, itheta+i2, iphi+i3, 3)
          bm  = sqrt(bx*bx+by*by+bz*bz)
          gbx = gbgrid(irr+i1, itheta+i2, iphi+i3, 1)
          gby = gbgrid(irr+i1, itheta+i2, iphi+i3, 2)
          gbz = gbgrid(irr+i1, itheta+i2, iphi+i3, 3)
          phic(i1+1,i2+1,i3+1) = (by*gbz-bz*gby)/bm/bm
        end do
      end do
    end do
    phi = trilinear(phic,px)
    bxgb2(1) = phi

    do i1 = 0, 1
      do i2 = 0, 1
        do i3 = 0, 1
          bx = bgrid(irr+i1, itheta+i2, iphi+i3, 1)
          by = bgrid(irr+i1, itheta+i2, iphi+i3, 2)
          bz = bgrid(irr+i1, itheta+i2, iphi+i3, 3)
          bm = sqrt(bx*bx+by*by+bz*bz)
          gbx = gbgrid(irr+i1, itheta+i2, iphi+i3, 1)
          gby = gbgrid(irr+i1, itheta+i2, iphi+i3, 2)
          gbz = gbgrid(irr+i1, itheta+i2, iphi+i3, 3)
          phic(i1+1,i2+1,i3+1) = (bz*gbx-bx*gbz)/bm/bm
        end do
      end do
    end do
    phi = trilinear(phic,px)
    bxgb2(2) = phi

    do i1 = 0, 1
      do i2 = 0, 1
        do i3 = 0, 1
          bx = bgrid(irr+i1, itheta+i2, iphi+i3, 1)
          by = bgrid(irr+i1, itheta+i2, iphi+i3, 2)
          bz = bgrid(irr+i1, itheta+i2, iphi+i3, 3)
          bm = sqrt(bx*bx + by*by + bz*bz)
          gbx = gbgrid(irr+i1, itheta+i2, iphi+i3, 1)
          gby = gbgrid(irr+i1, itheta+i2, iphi+i3, 2)
          gbz = gbgrid(irr+i1, itheta+i2, iphi+i3, 3)
          phic(i1+1,i2+1,i3+1) = (bx*gby-by*gbx)/bm/bm
        end do
      end do
    end do
    phi = trilinear(phic,px)
    bxgb2(3) = phi


    phic(1,1,1) = abs(b1rsgrid(irr,   itheta,   iphi))
    phic(2,1,1) = abs(b1rsgrid(irr+1, itheta,   iphi))
    phic(1,2,1) = abs(b1rsgrid(irr,   itheta+1, iphi))
    phic(2,2,1) = abs(b1rsgrid(irr+1, itheta+1, iphi))
    phic(1,1,2) = abs(b1rsgrid(irr,   itheta,   iphi+1))
    phic(2,1,2) = abs(b1rsgrid(irr+1, itheta,   iphi+1))
    phic(1,2,2) = abs(b1rsgrid(irr,   itheta+1, iphi+1))
    phic(2,2,2) = abs(b1rsgrid(irr+1, itheta+1, iphi+1))
    phi = trilinear(phic,  px)
    df = trilineardif(phic, px)
    b1rs = phi
    gb1rs(1) = df(1) / ((rss-1)/N_R)
    gb1rs(2) = df(2) / (PI/N_THETA)/r(1)
    gb1rs(3) = df(3) / (TWOPI/N_PHI)/r(1)/sintheta

    nplus = 0
    nminus = 0
    bplus = 0.0
    bminus = 0.0
    do m3 = 0, 1
      do m2 = 0, 1
        do m1 = 0, 1
          if (b1rsgrid(irr+m1, itheta+m2, iphi+m3) > 0) then
            nplus = nplus + 1
            bplus(1) = bplus(1) + bgrid(irr+m1,itheta+m2,iphi+m3,1)
            bplus(2) = bplus(2) + bgrid(irr+m1,itheta+m2,iphi+m3,2)
            bplus(3) = bplus(3) + bgrid(irr+m1,itheta+m2,iphi+m3,3)
          else
            nminus = nminus + 1
            bminus(1) = bminus(1) + bgrid(irr+m1,itheta+m2,iphi+m3,1)
            bminus(2) = bminus(2) + bgrid(irr+m1,itheta+m2,iphi+m3,2)
            bminus(3) = bminus(3) + bgrid(irr+m1,itheta+m2,iphi+m3,3)
          end if
        end do
      end do
    end do
    if ((nplus > 0) .and. (nminus > 0)) then
      bplus = bplus / nplus
      bminus = bminus / nminus
      if (dot_product(bplus, b) > dot_product(bminus, b)) then
        pol = 1.0
      else
        pol = -1.0
      end if
    else
      if (nplus == 0) pol = -1.0
      if (nminus == 0) pol = 1.0
    end if

    bmag = norm2(b(1:3))
    dbbds = dot_product(b(1:3), gbmag(1:3)) / (bmag * bmag)
    cvtu = 0.0 ! potential field

    if (r1(1) >= rss) then
      rsovr = rss / r1(1)
      if (b(1) >= 0) then
        plbs = 1
      else
        plbs = -1
      end if
      b(1) = b(1) * rsovr**2
      vpl = solarwind1(r1)
      tanp = vpl(3) / vpl(1)
      oneptan2 = 1 + tanp**2
      sroneptan2 = sqrt(oneptan2)
      b(3) = b(1) * tanp
      bmag = norm2(b(1:3))
      dvr1 = vsw * (2*k4ok2/r1(1)**3 + 4*k6ok2/r1(1)**5) / &
        (1 + k4ok2/r1(1)**2 + k6ok2/r1(1)**4)**2
      dvr1 = dvr1/vpl(1)
      dvf1 = -omega*sintheta / vpl(3)
      dvf2 = -omega*(r1(1)-rss)*cos(r1(2))/r1(1) / vpl(3)
      dvf2 = dvf2/vpl(3)
      gbmags = gbmag
      gbmag(1) = -2*bmag/r1(1) + bmag*tanp**2/oneptan2*(dvf1-dvr1)
      gbmag(2) = gbmags(2) * rsovr**3 * sroneptan2&
        + bmag * tanp**2 / oneptan2 * dvf2
      gbmag(3) = gbmags(3)*rsovr**3*sroneptan2
      bxgb2s = bxgb2
      bxgb2(1) = - bxgb2s(3)*rsovr*tanp/sroneptan2&
                 - b(1)*tanp**3/oneptan2*dvf2/bmag
      bxgb2(2) = bxgb2(2)*rsovr/sroneptan2&
                 - 2*b(1)/r1(1)*tanp/bmag&
                 + b(1)*tanp**3/oneptan2*(dvf1-dvr1)/bmag
      bxgb2(3) = bxgb2(3)*rsovr/sroneptan2&
                 + b(1)*tanp**2/oneptan2*dvf2/bmag
      cvtu(1) = b(3) * (dvf2 + cos(r1(2))/sintheta/r1(1)) &
        + tanp * rsovr**3 * plbs * gbmags(2)
      cvtu(2) = plbs * gbmags(3) * rsovr**3 &
        + b(3) * (1/r1(1) - dvf1 + dvr1)
      cvtu(3) = -plbs * gbmags(2) * rsovr**3
      dbbds = dot_product(b, gbmag) / bmag / bmag
      b1rs = b1rs / rsovr**2 / sroneptan2
      gb1rs(1) = 2 * b1rs / r1(1) - b1rs * tanp**2 / oneptan2 * (dvf1-dvr1)
      gb1rs(2) = gb1rs(2) / rsovr / sroneptan2 - b1rs * tanp**2 / oneptan2 * dvf2
      gb1rs(3) = gb1rs(3) / rsovr / sroneptan2
    end if
  end subroutine ! drvbmag
end module
