program sim3d
  !  Source injection at shocks
  !  time backward simulation to calculate fluxes at locations in IP
  !  pfss magnetic field model
  !  calculation is done in corotation reference frame
  !  pitch angle with outward magnetic field line
  !  pitch angle diffusion (symmetric D_{\mu\mu})
  !  perpendicular diffusion added

  use iso_fortran_env, only: real64, output_unit
  use datetime_utils, only: caldate
  use param, only: PI, NSPMAX, NFMAX, N_R, N_THETA, N_PHI, CSPEED, GAMMA_CS, &
                   nseedmax, bgrid, gbgrid, b1rsgrid, epsilon
  use cme_cross, only: inorout, preparecme
  use sim3d_utils, only: f0mod, compress, solarwindtemp, split, drvbmag, vfunc
  use epv, only: rp2e, e2p
  use fb, only: fb0, preparefb, fs0
  use mtrx, only: norm2, mrtx
  use loadptcl, only: prepareptcl
  use dmumu, only: preparedmumu
  use dxx, only: preparedxx, set_rlambda, set_rlambdax
  use file_op, only: read_maggrid, read_b1rs, read_param, record_nodes, fl_open, write_head
  use random, only: gasdev

  implicit none

  include 'omp_lib.h'

  integer, parameter  :: NM1 = 16, NMXID = 40

  real(kind=real64)   :: rpb(5), rp0(5)
  !real(kind=real64)   :: rp0org(5)
  real(kind=real64)   :: r0(3), rb(3)
  real(kind=real64)   :: x0(6)

  real(kind=real64)   :: rnz, rnm
  common/specie/rnz,rnm

  integer             :: npp, num(3)
  real(kind=real64)   :: t0, t, tnp

  integer             :: nsucmin
  common/nsucmin/nsucmin

  integer             :: npmax
  common/npmax/npmax

  real(kind=real64)   :: rb0, rmax, rk, deltat, tc, tl, tmodel0
  integer             :: nfbconst
  common/fbcnst/rb0,rmax,rk,deltat,tc,tl,tmodel0,nfbconst

  integer             :: ndpdt
  common /ndpdt/ndpdt

  integer             :: nodr(NSPMAX)
  real(kind=real64)   :: t0sv(2**(NSPMAX+1)), cksv(2**(NSPMAX+1))
  real(kind=real64)   :: rpbsv(5,2**(NSPMAX+1))
  common /svsp/nodr,t0sv,cksv,rpbsv

  real(kind=real64)   :: t0org, te, tdl, dmapjul, tcme0
  common /tmprm/t0org,te,tdl,dmapjul,tcme0

  real(kind=real64), allocatable :: tf(:), rf(:,:), ef(:), rmuf(:)
  integer             :: np, nf
  !real(kind=real64)   :: tf(NFMAX), rf(3,NFMAX), ef(NFMAX), rmuf(NFMAX)
  !common /ldptcl/tf,rf,ef,rmuf,np,nf

  character(len=256)  :: dir
  common /dir/dir

  integer             :: nodes, chunk
  integer             :: id

  character(len=2)    :: rankstr
  common/rankstr/rankstr

  integer             :: iseed
  integer, allocatable:: nseeds(:)

  real(kind=real64)   :: densw0, vsw, k4ok2, k6ok2, omega, b1au, vom, facip
  common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip

  real(kind=real64)   :: trgtfs(4), sp, sp0, gp, ap, scanw, h0
  common /srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0

  integer             :: nsplvl
  real(kind=real64)   :: df0(3), ddf0(3,3)
  real(kind=real64)   :: b1s, gb1s(3)
  real(kind=real64)   :: bv0(3), bm, cvtu(3), gbmag(3), bxgb2(3), dbbds, pol
  real(kind=real64)   :: dnsk0, dnsk00
  real(kind=real64)   :: vsk
  !real(kind=real64)   :: vnr(3)
  real(kind=real64)   :: vnx(3)

  real(kind=real64)   :: vsksw, tauf, tauc1_0, tauc2, tauc2_0, vcme0kmPs, vcme0 !LC
  common /vsksw/vsksw,tauf,tauc1_0,tauc2,tauc2_0,vcme0kmPs,vcme0 !LC

  integer             :: i, j, lsp, ns

  real(kind=real64)   :: df0dmu, ddf0dmu2
  real(kind=real64)   :: flux, dflux, flx, dflx
  real(kind=real64)   :: doy
  real(kind=real64)   :: e0, f1, pa0, ck
  real(kind=real64)   :: fs, pab, hb, rate, tsp, tb
  real(kind=real64)   :: fb_
  real(kind=real64)   :: rdf, treal, tod
  integer             :: itjul, iyear, iyday

  character(len=*), parameter :: writefmt = "(i1,i3,8(1pe13.5),i3)"

  integer :: nfl, nsts

  integer :: ndmumu
  common /ndmumu/ndmumu

  integer :: ndxx
  common /ndxx/ndxx

  !real(kind=real64) :: b1rsgrid(0:N_R, 0:N_THETA, 0:N_PHI) ! dummy storage

  densw0 = 166410.0 !332820.d0
  k4ok2 = 12.4242d0
  k6ok2 = 242.4242d0

  nodes = 1 ! default, will be read in readparm from 'dir.dat'
  if (nodes > NMXID) stop 'nodes > NMXID'
  chunk = 1

  ! USES environment variables PARAM_OUTDIR_PATH PARAM_NODES
  call read_param(nodes, nseeds)

  ! USES environment variables MAGGRID_FILE
  call read_maggrid(bgrid, gbgrid)
  ! USES environment variables B1RS_FILE
  call read_b1rs(b1rsgrid)

  ! USES environment variables PTCL_FILE
  ! this also allocates each of the arrays
  call prepareptcl(tf, rf, ef, rmuf, np, nf)
  print *, "nf = ", nf
  print *, "tf = ", tf
  print *, "np = ", np

  print "(a,i0,a)", "Preparing dmumu (ndmumu = ", ndmumu,")"
  call preparedmumu(ndmumu)
  print *, "Finished preparing dmumu"

  call preparedxx(ndxx)

  call preparefb

  call record_nodes(nodes)
  print *, 'nodes =',nodes

  !   normalize mean free path at 1 GV
  e0 = rp2e(1.0d0)
  call set_rlambda(e0)
  call set_rlambdax(e0)

  ! USES environment variable CME_DATA_FILE
  call preparecme
  print *, "tcme0 =", tcme0

  num = 1
  iseed = -1
  te = 0.0
  id = 0
  write(rankstr,"(i0.2)") id
  call fl_open(nfl, nsts)
  call write_head(nfl, np)

  ! __LINE__ is compiler dependent, and preprocessing is required
  !print *, __FILE__, ":", __LINE__, ": Finished preloop preparation"

  print *, "nf = ", nf

  do i = 1, nf
    print *, "i = ", i
    flx = 0.0
    dflx = 0.0
    tnp = 0.0
    rp0(:) = [rf(1,i), rf(2,i)*PI/180, rf(3,i)*PI/180, e2p(ef(i)), rmuf(i)]
    t0org = tf(i) ! comes from ptcl
    te = t0org - tcme0 !time from cme onset
    print *, "te = ", te
    print *, "tdl = ", tdl
    lsp = floor(te/tdl)
    print *, "lsp = ", lsp
    npp = ceiling(1.d0 * np * (lsp+1) / (2.0**(lsp+1) - 1))
    print *, "npp = ", npp
    flush(output_unit)
    if (npp <= 0) npp = 1

    sp = sp0 + (gp-sp0) * log(1.0 + te/180.0)
    !call setasgp(t0org,rp0)
    r0(1:3) = rp0(1:3)
    pa0 = rp0(5)
    call f0mod(r0, pa0, h0, df0, ddf0, df0dmu, ddf0dmu2)
    call drvbmag(r0, bv0, bm, cvtu, gbmag, bxgb2, dbbds, pol, b1s, gb1s)

    ! Convert spherical coordinates to Cartesian
    x0(1) = r0(1) * sin(r0(2)) * cos(r0(3))
    x0(2) = r0(1) * sin(r0(2)) * sin(r0(3))
    x0(3) = r0(1) * cos(r0(2))

    call inorout(t0org, x0, dnsk00, vsk, vnx)
    !print *, __FILE__, ":", __LINE__, ":", "Finished calling inorout"

    write(nsts,*) 'For flux at point', i
    write(nsts,*) 'Time,postion,energy/n,\mu'
    write(nsts,"(f12.5,3f11.4,e13.5,f9.5,2(1pe13.5))") tf(i),rf(:,i),ef(i),rmuf(i)
    write(nsts,"(' sp = ',1pe12.4,'; ap = ',1pe12.4)") sp, ap
    write(nsts,*) 'number of original particles npp = ', npp
    write(nsts,*) 'Nonzero source or boundary value list (f1)'
    write(nsts,*) 'Indicator: 0=summed,1=source,2=boundary value'
    write(nsts,*) 'List may be longer than npp due to split'
    !call flush(nsts)
    flush(nsts)

    !$OMP  PARALLEL NUM_THREADS(nodes) DEFAULT(firstprivate)&
    !$OMP& SHARED(rp0,h0,chunk,np,nseeds,bgrid,gbgrid,flx,dflx)
    !!$OMP& SHARED(cvgrid)
    id = OMP_GET_THREAD_NUM()
    iseed = nseeds(id+1)
    !$OMP DO SCHEDULE(DYNAMIC,chunk) REDUCTION(+:flx,dflx)
    do j = 1, npp
      print *, "Running step i = ", i, ", j = ", j
      ck = 0.0
      fs = 0.0
      nsplvl = 0
      t0 = 0.0
      tsp = t0 + tdl
      ns = 0
      dnsk0 = dnsk00
      if (tsp > te) tsp = te

      call walk3d(iseed, rp0, rpb, ck, fs, t0, t, tsp, ns, dnsk0, bv0, nsplvl)

      if (ns == -1 .and. tsp < te) then
        call split(iseed, rpb, ck, fs, t, nsplvl, dnsk0, bv0, flx, dflx, walk3d, nsts)
      else
        rb(1:3) = rpb(1:3)
        pab = rpb(5)
        call f0mod(rb, pab, hb, df0, ddf0, df0dmu, ddf0dmu2)
        rate = exp(ck - hb + h0)
        if (ns >= 0) ns = 1
        fb_ = fb0(tb,rpb) * rate
        f1 = (fs+fb_) / 2 ** nsplvl
        if (fs > 1.d-30) write(nsts,writefmt) 0,nsplvl,fs
        !if (fb_ > 1.d-30) write(nsts,writefmt) 2,nsplvl,fb_,t0org-t,rate,rpb,ns
        !call flush(nsts)
        flush(nsts)
        !OMP CRITICAL sum
        flx = flx + f1
        dflx = dflx + f1*f1
        !OMP END CRITICAL sum
      end if
    end do
    !$OMP END DO
    !$OMP BARRIER
    !$OMP END PARALLEL
    !write(nsts,*) 'end of list'
    !call flush(nsts)
    flush(nsts)
    flux = flx/npp
    dflux = sqrt(dflx)/npp
    rdf = dflux/flux
    flux = flux*rp0(4)**2*3e7 !flux in 1/(cm^2 s sr MeV/n)
    dflux = dflux*rp0(4)**2*3e7
    if (rnm > 0.5) then
      flux = flux * rnm / rnz
      dflux = dflux * rnm / rnz
    end if
    treal = dmapjul + tf(i) / 1440.0
    itjul = floor(treal)
    call caldate(treal, iyear, iyday)
    tod = treal - itjul
    doy = iyday + tod
    write(nfl,"(i4,f12.7,7(1pe12.4))") iyear, doy, rf(:,i), ef(i), rmuf(i), flux, dflux
    !call flush(nfl)
    flush(nfl)
    print *, ""
  end do
  !print *, __FILE__, ":", __LINE__, ": Exited loop with", i, "runs"
  close(nfl)
  close(nsts)

contains


  !  random walk of energetic particles in magnetic
  !  variables: t, xp(5)
  !  x - spatial coordinators
  !  p - momentum
  !  pa - pitch angle

  subroutine walk3d(iseed, rp0, rpb, ck, fs, t0, t, tsp, ns, dnsk0, bv0, nsplvl)

    ! rp0,rpb = (r,theta, phi, p, pa)  theta=const, phi follows Parker spiral
    ! initial or boundary value
    integer, intent(in out)           :: iseed
    real(kind=real64), intent(in)     :: rp0(5)
    real(kind=real64), intent(out)    :: rpb(5)
    real(kind=real64), intent(in out) :: ck, fs
    real(kind=real64), intent(in)     :: t0, tsp
    real(kind=real64), intent(out)    :: t
    integer, intent(out)              :: ns
    real(kind=real64), intent(in out) :: dnsk0, bv0(3)
    integer, intent(in)               :: nsplvl

    real(kind=real64), parameter :: rdpmax = 100

    !integer           :: id
    !integer           :: is(3)
    real(kind=real64) :: xpk(6)
    !real(kind=real64) :: x(3), vx(3), bv(3)
    real(kind=real64) :: r(3)
    real(kind=real64) :: dxpkdt(6), gxw2(3), gxw3(3)
    !real(kind=real64) :: dxpkdt1(6), gxw2s(3), gxw3s(3)
    !real(kind=real64) :: dxpk1(6), dxpk2(6), xpk1(6)
    real(kind=real64) :: dxpk(6)
    real(kind=real64) :: cvtu(3), cvtu0(3), b1s
    !real(kind=real64) :: gb1s(3)
    !real(kind=real64) :: gbmag(3), bxgb2(3)
    !real(kind=real64) :: culper(3), culpar(3)
    !real(kind=real64) :: uax1(3), uax2(3), uax3(3), uax20(3)
    real(kind=real64) :: e(2), sqrte(2)
    real(kind=real64) :: dw(3)
    !real(kind=real64) :: bw(3)
    !real(kind=real64) :: gb(3), gr(3,3), dgr(3), dgx(3)
    !real(kind=real64) :: b2x(3,3), b2r(3,3)
    real(kind=real64) :: r2x(3,3)
    real(kind=real64) :: vpl(3)!, gvpl(3,3)

    !real(kind=real64) :: vd(3)

    integer           :: nlambdaconst
    common/nlambdaconst/nlambdaconst

    real(kind=real64) :: densw0, vsw, k4ok2, k6ok2, vom, facip, b1au, omega
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip

    integer           :: ndpdt
    common/ndpdt/ndpdt

    real(kind=real64) :: p, pa, p0, pa0, pas, hs

    real(kind=real64) :: t0org, te, tdl, dmapjul, tcme0
    common /tmprm/t0org,te,tdl,dmapjul,tcme0

    real(kind=real64) :: rnz, rnm
    common/specie/rnz,rnm

    real(kind=real64) :: rb0, rmax, rk, deltat, tc, tl, tmodel0

    integer           :: nfbconst
    common/fbcnst/rb0,rmax,rk,deltat,tc,tl,tmodel0,nfbconst

    real(kind=real64) :: trgtfs(4), sp, sp0, gp, ap, scanw, h0
    common/srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0

    integer           :: iw, n, isp, kf
    real(kind=real64) :: df0(3), ddf0(3,3)
    real(kind=real64) :: vsk, vnr(3), vnx(3)
    real(kind=real64) :: rsh(3), fs1, dt, srdt, du, densw, gper, tsh, ob, vptc
    real(kind=real64) :: dtmax, dtmin1, dtmin2, dtmin3, dtmin4, dtmin5
    real(kind=real64) :: bm, amach, va, u1, vswn, costheta2, sintheta2, df0dmu, ddf0dmu2, dlt
    real(kind=real64) :: vs, tempsw, dtsk, dxdtn, g2n, pinj, pc, rate, rbf, rprs, rshf
    real(kind=real64) :: tempds, vthds, smach, tacc, dnsk

    ! e = [5e-6, 0.005]
    e = [5e-8, 0.005]

    r(1:3) = rp0(1:3)
    xpk(1) = r(1) * dsin(r(2)) * dcos(r(3))
    xpk(2) = r(1) * dsin(r(2)) * dsin(r(3))
    xpk(3) = r(1) * dcos(r(2))
    p0 = rp0(4)
    pa0 = rp0(5)
    p = p0
    pa = pa0
    fs1 = 0.0
    xpk(4) = p
    xpk(5) = pa
    xpk(6) = ck
    dnsk = dnsk0

    t = t0 !0.0
    ns = 0

    !print *, "t=", t
    !ck = ck0 !0.0
    n = 0
    iw = 1
    isp = 1
    sqrte(1:2) = sqrt(e(1:2))
    dt = e(isp)
    srdt = sqrte(isp)

    do while (iw == 1)

      !write(59,"(f14.6,12(1pe14.5))") t,xpk,fs1,fs
      call vfunc(t, xpk, dxpkdt, du, gxw2, gxw3, bv0, densw, vpl, gper, b1s)

      dtmax = epsilon(2) * r(1) / CSPEED
      dt = dtmax
      dtmin1 = epsilon(1)**2/(2*du)
      if (dt > dtmin1) dt = dtmin1
      dtmin2 = abs(epsilon(1)/dxpkdt(5))
      if (dt > dtmin2) dt = dtmin2
      dtmin3 = (epsilon(2)*r(1))**2/(2.*gper)
      if (dt > dtmin3) dt = dtmin3
      !vptc = norm2(dxpkdt(1:3))
      dtmin4 = epsilon(2) * r(1) / norm2(dxpkdt(1:3))
      if (dt > dtmin4) dt = dtmin4
      dtmin5 = epsilon(4) / abs(dxpkdt(6))
      if (dt > dtmin5) dt = dtmin5
      !   slow down at shock
      !ddif2 = 16*gper*gper/vptc/vptc
      !if (dnsk0*dnsk0 < ddif2) then
      !  dtmin6 = 0.125*gper/vptc**2
      !  if (dt > dtmin6) dt = dtmin6
      !end if
      dt = dt/2
      srdt = sqrt(dt)

      !  use EM scheme
      dw(1) = gasdev(iseed)
      dxpk = dxpkdt*dt
      dxpk(5) = dxpk(5) + sqrt(2*du)*dw(1)*srdt
      xpk = xpk + dxpk
      if (xpk(5) > 1.0) xpk(5) = 0.99999
      if (xpk(5) < -1.0) xpk(5) = -0.99999

      !    g1,2,3 ---- kappa z,x,y
      !   step forward stochastic differential equation (Euler scheme)
      !   calculate Wiener noise for spatial diffusion
      dw(2) = gasdev(iseed)
      dw(3) = gasdev(iseed)

      !   calculate increments due to spatial diffusion term
      xpk(1:3) = xpk(1:3) + (gxw2(1:3)*dw(2) + gxw3(1:3)*dw(3)) * srdt

      t = t + dt
      n = n + 1
      !   change the position to spheric coordinates
      r(1) = norm2(xpk(1:3))
      r(2) = dacos(xpk(3)/r(1))
      r(3) = datan2(xpk(2),xpk(1))
      !  sum source
      tsh = t0org - t
      call inorout(tsh, xpk, dnsk, vsk, vnx)
      if (dnsk0 * dnsk < 0) then
        bm = norm2(bv0)
        r2x = mrtx(sin(r(2)), cos(r(2)), sin(r(3)), cos(r(3)))
        vnr(1:3) = vnx(1)*r2x(1,1:3) + vnx(2)*r2x(2,1:3) + vnx(3)*r2x(3,1:3)
        ob = acos(abs(dot_product(bv0, vnr)) / bm)
        vswn = sum(vpl * vnr)
        u1 = vsk - vswn
        !if (u1 < 0) write(*,*) 'NaN1'
        va = 187.8*bm / sqrt(densw)
        amach = u1/va
        tempsw = solarwindtemp(r)
        vs = 7.83e-6 * sqrt(GAMMA_CS * tempsw)
        smach = u1/vs
        if (amach > 1.0) then
          rsh = compress(amach, smach, ob)
          costheta2 = cos(ob)**2
          sintheta2 = 1-costheta2
          rshf = 1
          do kf = 1, 3
            rbf = rsh(kf) * (amach*amach-costheta2) / (amach*amach-rsh(kf)*costheta2)
            if (rbf > 1.0000001d0 .and. rsh(kf) > 1.0000001d0) rshf = rsh(kf)
          end do
          rprs = 1 + GAMMA_CS*smach*smach*(rshf-1)/rshf*(1-rshf*sintheta2*((rshf+1)&
            * amach * amach - 2 * rshf * costheta2) / 2 / (amach*amach-rshf*costheta2)**2)
          tempds = tempsw / rshf * rprs
          vthds = 7.83e-6 * sqrt(2*tempds)
          !    vthds ~ Vsh due to shock heating (mainly protons)
          if (dxpkdt(4) > 0) then
            tacc = xpk(4) / dxpkdt(4)
          else
            tacc = 1e10
          end if
          dtsk = tsh - tcme0
          if (tacc > dtsk) tacc = dtsk
          fs1 = fs0(tacc, xpk, bm, u1, densw, ob, amach, rshf, vthds, pinj, pc)
          pas = xpk(5)
          ck = xpk(6)
          call f0mod(r, pas, hs, df0, ddf0, df0dmu, ddf0dmu2)
          rate = exp(ck - hs + h0)
          g2n = gper * sin(ob)**2
          dxdtn = dxpkdt(1)*vnx(1) + dxpkdt(2)*vnx(2) + dxpkdt(3)*vnx(3)
          dlt = abs(dnsk) / (g2n + dxdtn**2*dt/2)
          fs = fs + fs1 * dlt * rate
          !   output to analyze source injection profiles
          !if (fs1 > 0.0) then
          !  write(nsts,"(i1,i3,8(1pe13.5),i3)") 1,nsplvl,dtsk,xpk,fs1,rate,dlt,u1,amach,smach,ob
          !end
          !   if local shock acceleration injection cutoff, add to p instead of source
          if (pc < xpk(4)) xpk(4) = xpk(4) * (1 - u1*(1-1/rshf)/3*dlt)
        end if
      end if

      dnsk0 = dnsk

      if ((xpk(4) < rp0(4)/rdpmax) .or. (xpk(4) > rdpmax*rp0(4))) then
        iw = 0
        ns = -2
      end if
      if (r(1) > rmax) then
        iw = 0
        ns = -3
      end if
      if (r(1) < rb0) then
        iw = 0
        ns = n
      end if
      if (t > tsp) then
        iw = 0
        ns = -1
      end if
      if (t > te) then
        iw = 0
        ns = -4
      end if
    end do

    rpb(1:3) = r(1:3)
    rpb(4:5) = xpk(4:5)
    ck = xpk(6)
  end subroutine

end program
