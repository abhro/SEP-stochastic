program shockfront
  !  Source injection at shocks
  !  time backward simulation to calculate fluxes at locations in IP
  !  pfss magnetic field model
  !  calculation is done in corotation reference frame
  !  pitch angle with outward magnetic field line
  !  pitch angle diffusion (symmetric D_{\mu\mu})
  !  perpendicular diffusion added
  use iso_fortran_env, only: real64
  use datetime_utils, only: caldate
  use param, only: PI, GAMMA_CS, NSPMAX, NFMAX, NSEEDMAX, &
                   N_R, N_THETA, N_PHI, epsilon, bgrid, gbgrid, b1rsgrid!, cvgrid
  use epv, only: rp2e, e2p
  use cme_cross, only: inorout, preparecme, locate
  use fb, only: fb0, preparefb
  use sim3d_utils, only: f0mod, compress, solarwindtemp, drvbmag, vfunc, split
  use mtrx, only: mrtx, mxptr, norm2
  use file_op, only: read_maggrid, read_b1rs, read_param, record_nodes, fl_open, write_head
  use loadptcl, only: prepareptcl
  use dmumu, only: preparedmumu
  use dxx, only: preparedxx, set_rlambda, set_rlambdax
  use random, only: gasdev

  implicit none

  include 'omp_lib.h'

  integer, parameter  :: NM1 = 16, NMXID = 40

  real(kind=real64)   :: rpb(5), rp0(5), rp0org(5)
  real(kind=real64)   :: r0(3), rb(3)
  real(kind=real64)   :: x0(6)

  real(kind=real64)   :: rnz, rnm
  common/specie/rnz,rnm

  integer             :: nsucmin
  common/nsucmin/nsucmin

  integer             :: npmax
  common/npmax/npmax

  integer             :: npp, n1, num(3)

  real(kind=real64)   :: rb0, rmax, rk, deltat, tc, tl, tmodel0
  integer             :: nfbconst
  common/fbcnst/rb0,rmax,rk,deltat,tc,tl,tmodel0,nfbconst

  real(kind=real64)   :: t0, t, tnp

  integer             :: ndpdt
  common/ndpdt/ndpdt

  real(kind=real64)   :: t0sv(2**(NSPMAX+1)), cksv(2**(NSPMAX+1))
  real(kind=real64)   :: rpbsv(5,2**(NSPMAX+1))
  integer             :: nodr(NSPMAX)
  common /svsp/nodr,t0sv,cksv,rpbsv

  real(kind=real64)   :: t0org, te, tdl, dmapjul, tcme0
  common /tmprm/t0org,te,tdl,dmapjul,tcme0

  real(kind=real64), allocatable :: tf(:), rf(:,:), ef(:), rmuf(:)
  !real(kind=real64)   :: tf(NFMAX), rf(3, NFMAX), ef(NFMAX), rmuf(NFMAX)
  integer             :: np,nf
  !common/ldptcl/tf,rf,ef,rmuf,np,nf

  character(len=256)  :: dir
  common /dir/dir

  integer             :: nodes, chunk, id

  character(len=2)    :: rankstr
  common/rankstr/rankstr

  integer             :: iseed

  integer, allocatable:: nseeds(:)

  real(kind=real64)   :: densw0, vsw, k4ok2, k6ok2, omega, b1au, vom, facip
  common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip

  real(kind=real64)   :: trgtfs(4), ap, h0, gp, sp, sp0, scanw
  common /srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0

  integer :: ndmumu
  common /ndmumu/ndmumu

  integer :: ndxx
  common /ndxx/ndxx

  real(kind=real64)   :: df0(3), ddf0(3,3)
  integer             :: nsplvl
  real(kind=real64)   :: bv0(3), bm, cvtu(3), gbmag(3), bxgb2(3), dbbds
  real(kind=real64)   :: pol, b1s, gb1s(3)
  real(kind=real64)   :: dnsk0, vsk, vnr(3), vnx(3)

  integer             :: i, j, lsp, ns
  real(kind=real64)   :: df0dmu, ddf0dmu2, dflux, dflx, doy, flx, flux
  real(kind=real64)   :: e0, f1, pa0, ck
  real(kind=real64)   :: fs, pab, hb, rate, tsp, tb
  real(kind=real64)   :: fb_
  real(kind=real64)   :: rdf, treal, tod
  integer             :: itjul, iyear, iyday
  character(len=*), parameter :: writefmt = "(i1,i3,8(1pe13.5),i3)"
  integer :: nfl, nsts

  !real(kind=real64) :: bgrid(0:N_R, 0:N_THETA, 0:N_PHI, 3)
  !real(kind=real64) :: gbgrid(0:N_R, 0:N_THETA, 0:N_PHI, 3)
  !real(kind=real64) :: b1rsgrid(0:N_R, 0:N_THETA, 0:N_PHI) ! dummy storage

  densw0 = 66410.0 !332820.d0
  k4ok2 = 12.4242d0
  k6ok2 = 242.4242d0

  nodes = 1 ! default, will be read in readparm from 'dir.dat'
  if (nodes > NMXID) stop 'nodes > NMXID'
  chunk = 1
  call read_param(nodes, nseeds)
  nodes = 1!!! force nodes = 1
  call read_maggrid(bgrid, gbgrid)
  call read_b1rs(b1rsgrid)
  call prepareptcl(tf, rf, ef, rmuf, np, nf)
  call preparedmumu(ndmumu)
  call preparedxx(ndxx)
  call preparefb
  call record_nodes(nodes)
  write(*,*) 'nodes =', nodes
  !   normalize mean free path at 1 GV
  e0 = rp2e(1.0d0)
  call set_rlambda(e0)
  call set_rlambdax(e0)
  call preparecme

  num = 1
  iseed = -1
  te = 0.0
  id = 0
  write(rankstr,"(i0.2)") id
  call fl_open(nfl, nsts)
  call write_head(nfl, np)

  ! ========== force 1 particle and no split =============================
  np = 1
  tdl = 1.e10
  ! =======================================

  do i = 1, nf
    flx = 0.0
    dflx = 0.0
    tnp = 0.0
    rp0(1) = rf(1,i)
    rp0(2) = rf(2,i)*PI/180
    rp0(3) = rf(3,i)*PI/180
    rp0(4) = e2p(ef(i))
    rp0(5) = rmuf(i)
    t0org = tf(i)
    te = t0org-tcme0 !time from cme onset
    lsp = floor(te/tdl)
    npp = ceiling(1.d0*np*(lsp+1)/(2.0**(lsp+1)-1))
    if (npp <= 0) npp = 1

    sp = sp0+(gp-sp0)*log(1.+te/180.)
    !call setasgp(t0org,rp0)
    r0(1:3) = rp0(1:3)
    pa0 = rp0(5)
    call f0mod(r0,pa0,h0,df0,ddf0,df0dmu,ddf0dmu2)
    call drvbmag(r0,bv0,bm,cvtu,gbmag,bxgb2,dbbds,pol,b1s,gb1s)
    x0(1) = r0(1) * sin(r0(2)) * cos(r0(3))
    x0(2) = r0(1) * sin(r0(2)) * sin(r0(3))
    x0(3) = r0(1) * cos(r0(2))
    call inorout(t0org,x0,dnsk0,vsk,vnx)

    write(nsts,*) 'For flux at point', i
    write(nsts,*) 'Time,postion,energy/n,\mu'
    write(nsts,"(f12.5,3f11.4,e13.5,f9.5,2(1pe13.5))") tf(i),rf(1:3,i),ef(i),rmuf(i)
    write(nsts,"(' sp = ',1pe12.4,'; ap = ',1pe12.4)") sp,ap
    write(nsts,*) 'number of original particles npp = ', npp
    write(nsts,*) 'Nonzero source or boundary value list (f1)'
    write(nsts,*) 'Indicator: 0=summed,1=source,2=boundary value'
    write(nsts,*) 'List may be longer than npp due to split'
    !call flush(nsts)
    flush(nsts)

    !$OMP  PARALLEL NUM_THREADS(nodes) DEFAULT(firstprivate)&
    !$OMP& SHARED(rp0,h0,chunk,np,nseeds,bgrid,gbgrid,flx,dflx)
    !! SHARED(cvgrid)
    id = OMP_GET_THREAD_NUM()
    iseed = nseeds(id+1)
    !$OMP DO SCHEDULE(DYNAMIC,chunk) REDUCTION(+:flx,dflx)
    do n1 = 1,npp
      ck = 0.0
      fs = 0.0
      nsplvl = 0
      t0 = 0.0
      tsp = t0+tdl
      ns = 0

      if (tsp > te) tsp = te
      call walk3d(iseed, rp0, rpb, ck, fs, t0, t, tsp, ns, dnsk0, bv0, nsplvl)
      if (ns == -1 .and. tsp < te) then
        call split(iseed, rpb, ck, fs, t, nsplvl, dnsk0, bv0, flx, dflx, walk3d, nsts)
      else
        rb(1:3) = rpb(1:3)
        pab = rpb(5)
        call f0mod(rb, pab, hb, df0, ddf0, df0dmu, ddf0dmu2)
        rate = exp(ck-hb+h0)
        if (ns >= 0) ns = 1
        fb_ = fb0(tb,rpb) * rate
        f1 = (fs+fb_) / 2 ** nsplvl
        if (fs > 1.d-40) write(nsts,writefmt) 0,nsplvl,fs
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
    flux = flx / npp
    dflux = sqrt(dflx) / npp
    rdf = dflux / flux
    flux = flux * rp0(4)**2 * 3e7 !flux in 1/(cm^2 s sr MeV/n)
    dflux = dflux * rp0(4)**2 * 3e7
    if (rnm > 0.5) then
      flux = flux * rnm / rnz
      dflux = dflux * rnm / rnz
    end if
    treal = dmapjul + tf(i)/1440.0
    itjul = floor(treal)
    call caldate(treal,iyear,iyday)
    tod = treal-itjul
    doy = iyday+tod
    write(nfl,"(i4,f12.7,7(1pe12.4))") iyear, doy, rf(:,i), ef(i), rmuf(i), flux, dflux
    !call flush(nfl)
    flush(nfl)
  end do
  close(nfl)
  close(nsts)
  close(54)
contains


  !>  random walk of energetic particles in magnetic
  !>  variables: t, xp(5)
  !>  - x - spatial coordinators
  !>  - p - momentum
  !>  - pa - pitch angle
  subroutine walk3d(iseed, rp0, rpb, ck, fs, t0, t, tsp, ns, dnsk0, bv0, nsplvl)

    ! rp0,rpb = (r,theta, phi, p, pa)  theta=const, phi follows Parker spiral
    ! initial or boundary value
    use param, only: CSPEED
    use fb, only: fs0
    implicit none
    real(kind=real64) :: ck, fs, t0, t, tsp
    integer           :: ns, nsplvl
    integer           :: id
    real(kind=real64), parameter :: rdpmax = 100
    integer           :: is(3)
    real(kind=real64) :: rp0(5),rpb(5)
    real(kind=real64) :: xpk(6),x(3),r(3),vx(3),bv(3),bv0(3)
    real(kind=real64) :: dxpkdt(6),gxw2(3),gxw3(3)
    real(kind=real64) :: dxpkdt1(6),gxw2s(3),gxw3s(3)
    real(kind=real64) :: dxpk1(6),dxpk2(6),dxpk(6),xpk1(6)
    real(kind=real64) :: cvtu(3),cvtu0(3),gbmag(3),bxgb2(3),b1s,gb1s(3)
    real(kind=real64) :: culper(3),culpar(3)
    real(kind=real64) :: uax1(3),uax2(3),uax3(3),uax20(3)
    real(kind=real64) :: e(2),sqrte(2)
    real(kind=real64) :: dw(3),bw(3)
    real(kind=real64) :: gb(3), gr(3,3), dgr(3),dgx(3)
    real(kind=real64) :: b2x(3,3),r2x(3,3),b2r(3,3)
    real(kind=real64) :: vpl(3),gvpl(3,3)
    integer           :: iseed
    real(kind=real64) :: vd(3)
    integer           :: nlambdaconst
    common/nlambdaconst/nlambdaconst
    real(kind=real64) :: densw0,vsw,k4ok2, k6ok2, vom, facip, b1au, omega
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip
    integer           :: ndpdt
    common/ndpdt/ndpdt
    !data e/5e-6, 0.005/
    data e/5e-8, 0.005/
    real(kind=real64) :: p,pa,p0,pa0,pas,hs
    real(kind=real64) :: te, tdl, tcme0, dmapjul, t0org
    common /tmprm/t0org,te,tdl,dmapjul,tcme0
    real(kind=real64) :: rnz, rnm
    common/specie/rnz,rnm
    real(kind=real64) :: rb0, rmax, rk, deltat, tc, tl, tmodel0
    integer           :: nfbconst
    common/fbcnst/rb0,rmax,rk,deltat,tc,tl,tmodel0,nfbconst
    real(kind=real64) :: trgtfs(4), sp, sp0, gp, ap, scanw, h0
    common/srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0
    integer           :: iw
    real(kind=real64) :: df0(3),ddf0(3,3)
    real(kind=real64) :: dnsk0,dnsk,vsk,vnr(3),vnx(3)
    real(kind=real64) :: rsh(3)
    real(kind=real64) :: fs1, dt, srdt, tsh, vptc
    integer           :: n, isp, i, kf
    real(kind=real64) :: du, densw, gper, dtmax, dtmin1, dtmin2, dtmin3, dtmin4, dtmin5
    real(kind=real64) :: bm, ob, amach, va, u1, vswn, costheta2, sintheta2, df0dmu, ddf0dmu2, dlt
    real(kind=real64) :: vs, tempsw, dtsk, dxdtn, g2n, pinj, pc, rate, rbf, rprs, rshf
    real(kind=real64) :: smach, tacc, tempds, vthds

    r(1:3) = rp0(1:3)
    xpk(1) = r(1)*dsin(r(2))*dcos(r(3))
    xpk(2) = r(1)*dsin(r(2))*dsin(r(3))
    xpk(3) = r(1)*dcos(r(2))
    p0 = rp0(4)
    pa0 = rp0(5)
    p = p0
    pa = pa0
    fs1 = 0.
    xpk(4) = p
    xpk(5) = pa
    xpk(6) = ck
    dnsk = dnsk0

    t = t0 !0.0
    ns = 0

    !write(*,*)'t =',t
    !ck = ck0 !0.0
    n = 0
    iw = 1
    isp = 1
    sqrte(1:2) = sqrt(e(1:2))
    dt = e(isp)
    srdt = sqrte(isp)

    do while (iw == 1)

      if (n == 64440) then
        close(69)
        stop
      end if
      t = t0org + 90.0 - 0.05*n
      tsh = t0org - t
      call sksurface(tsh, xpk, mod(n, 2))

      !write(59,"(f14.6,12(1pe14.5))") t,xpk,fs1,fs
      call vfunc(t,xpk,dxpkdt,du,gxw2,gxw3,bv0,densw,vpl,gper,b1s)

      dtmax = epsilon(2)*r(1)/CSPEED
      dt = dtmax
      dtmin1 = epsilon(1)**2/(2*du)
      if(dt > dtmin1) dt = dtmin1
      dtmin2 = abs(epsilon(1)/dxpkdt(5))
      if(dt > dtmin2) dt = dtmin2
      dtmin3 = (epsilon(2)*r(1))**2 / (2.*gper)
      if(dt > dtmin3) dt = dtmin3
      vptc = sqrt(dxpkdt(1)**2+dxpkdt(2)**2+dxpkdt(3)**2)
      dtmin4 = epsilon(2)*r(1)/vptc
      if(dt > dtmin4) dt = dtmin4
      dtmin5 = epsilon(4)/abs(dxpkdt(6))
      if(dt > dtmin5) dt = dtmin5
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
      dxpk = dxpkdt * dt
      dxpk(5) = dxpk(5) + sqrt(2*du)*dw(1)*srdt
      xpk = xpk+dxpk
      if (xpk(5) >  1.0) xpk(5) =  0.99999
      if (xpk(5) < -1.0) xpk(5) = -0.99999

      !    g1,2,3 ---- kappa z,x,y
      !   step forward stochastic differential equation (Euler scheme)
      !   calculate Wiener noise for spatial diffusion
      dw(2) = gasdev(iseed)
      dw(3) = gasdev(iseed)

      !   calculate increments due to spatial diffusion term
      xpk(1:3) = xpk(1:3)+(gxw2(1:3)*dw(2)+gxw3(1:3)*dw(3))*srdt

      t = t+dt
      n = n+1
      !   change the position to spheric coordinates
      r(1) = norm2(xpk(1:3))
      r(2) = dacos(xpk(3)/r(1))
      r(3) = datan2(xpk(2),xpk(1))
      !  sum source
      tsh = t0org-t
      call inorout(tsh, xpk, dnsk, vsk, vnx)
      if(dnsk0*dnsk<0) then
        bm = norm2(bv0(1:3))
        r2x = mrtx(sin(r(2)), cos(r(2)), sin(r(3)), cos(r(3)))
        vnr(1:3) = vnx(1)*r2x(1,1:3)+vnx(2)*r2x(2,1:3)+vnx(3)*r2x(3,1:3)
        ob = acos(abs(dot_product(bv0(1:3), vnr(1:3)))/bm)
        vswn = dot_product(vpl(1:3), vnr(1:3))
        u1 = vsk-vswn
        ! if (u1 < 0) write(*,*) 'NaN1'
        va = 187.8*bm/sqrt(densw)
        amach = u1/va
        tempsw = solarwindtemp(r)
        vs = 7.83e-6*sqrt(GAMMA_CS*tempsw)
        smach = u1/vs
        if (amach > -1.0d20) then
          rsh = compress(amach, smach, ob)
          costheta2 = cos(ob)**2
          sintheta2 = 1-costheta2
          rshf = 1
          do kf = 1,3
            rbf = rsh(kf)*(amach*amach-costheta2)/(amach*amach-rsh(kf)*costheta2)
            if ((rbf > 1.0000001d0) .and. (rsh(kf) > 1.0000001d0)) rshf = rsh(kf)
          end do
          rprs = 1+GAMMA_CS*smach*smach*(rshf-1)/rshf*(1-rshf*sintheta2*((rshf+1)&
            *amach*amach-2*rshf*costheta2)/2/(amach*amach-rshf*costheta2)**2)
          tempds = tempsw/rshf*rprs
          vthds = 7.83e-6*sqrt(2*tempds)
          !    vthds ~ Vsh due to shock heating (mainly protons)
          if(dxpkdt(4)>0) then
            tacc =(xpk(4)/dxpkdt(4))
          else
            tacc = 1e10
          end if
          dtsk = tsh-tcme0
          if (tacc > dtsk) tacc = dtsk
          fs1 = fs0(tacc,xpk,bm,u1,densw,ob,amach,rshf,vthds,pinj,pc)
          write(69,"(20e16.6e3)") xpk(1:3),ob,rshf,vthds,pinj,pc,u1,amach,smach,b1s,fs1,tsh
          pas = xpk(5)
          ck = xpk(6)
          call f0mod(r,pas,hs,df0,ddf0,df0dmu,ddf0dmu2)
          rate = exp(ck-hs+h0)
          g2n = gper*sin(ob)**2
          dxdtn = dxpkdt(1)*vnx(1)+dxpkdt(2)*vnx(2)+dxpkdt(3)*vnx(3)
          dlt = abs(dnsk)/(g2n+dxdtn**2*dt/2)
          fs = fs+fs1*dlt*rate
          !   output to analyze source injection profiles
          ! if (fs1 > 0.0) write(nsts,"(i1,i3,16(1pe13.5))") 1,nsplvl,dtsk,xpk,fs1,rate,dlt,u1,amach,smach,ob
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

  subroutine sksurface(tsh, x, n)
    real(kind=real64), intent(in) :: tsh
    real(kind=real64), intent(out):: x(6)
    integer, intent(in)           :: n


    real(kind=real64) :: vskf0, acsk, tska(20), pska(20,8)
    integer           :: nsk
    common /cmesk/vskf0,acsk,tska,pska,nsk

    real(kind=real64) :: t

    real(kind=real64) :: rskc, drskc, rskmax, rskmin, ask(3), dask(3)
    real(kind=real64) :: rskf, vskf
    real(kind=real64) :: thetaskc, phiskc, gmsk, dthetaskc, dphiskc, dgmsk
    real(kind=real64) :: sinthetask, costhetask, sinphisk, cosphisk
    real(kind=real64) :: r2xsk(3,3)
    real(kind=real64) :: xp2r(3,3), xp2x(3,3)
    real(kind=real64) :: xskc(3), xc(3), xpc(3)

    real(kind=real64) :: densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip

    real(kind=real64) :: t0org,te,tdl,dmapjul,tcme0
    common /tmprm/t0org,te,tdl,dmapjul,tcme0

    real(kind=real64) :: vsksw,tauf,tauc1_0,tauc2,tauc2_0,vcme0kmPs,vcme0 !LC
    common /vsksw/vsksw,tauf,tauc1_0,tauc2,tauc2_0,vcme0kmPs,vcme0

    integer           :: jt
    real(kind=real64) :: theta, phi, ep, rac, rbc, rhclf

    i = floor(n / 360.0)
    j = mod(n, 360)

    theta = (i+0.5) * PI/180.0
    phi   = (j+0.5) * PI/180.0

    if (tsh < tska(nsk)) then
      jt = locate(tska,nsk,tsh)
      if (jt == 0) jt = 1
      !  maximum shock radial distance for quick search
      drskc = (pska(jt+1,1)-pska(jt,1)) / (tska(jt+1)-tska(jt))
      rskc = pska(jt,1)+(tsh-tska(jt))*drskc
      dask(3) = (pska(jt+1,6)-pska(jt,6))/(tska(jt+1)-tska(jt))
      ask(3) = pska(jt,6)+(tsh-tska(jt))*dask(3)
    else
      jt = nsk-1
      vcme0 = vskf0
      if (tsh < tauc2_0) then
        rskf = pska(nsk,1)+pska(nsk,6)+vcme0*(tsh-tska(nsk)) !LC
        vskf = vcme0
      else
        if (tauc2_0>tska(nsk)) then
          rskf = 3.0/2.0*(vcme0-vsw)*tauc2*(((tsh-tska(1))/tauc2)&
            **(2.0/3.0)-1.0d0)+vsw*tauc2*(((tsh-tska(1))/tauc2)-1.0d0)&
            + pska(nsk,1)+pska(nsk,6)+vcme0*(tauc2_0-tska(nsk)) !LC
        else
          rskf = 3.0/2.0*(vcme0-vsw)*tauc2*(((tsh-tska(1))/tauc2)**(2.0/3.0&
            )-((tska(nsk)-tska(1))/tauc2)**(2.0/3.0))+vsw*(tsh-tska(nsk))&
            + pska(nsk,1)+pska(nsk,6)
        end if
        vskf = (vcme0-vsw)*(((tsh-tska(1))/tauc2)**(-1.0d0/3.0d0)) + vsw
      end if

      rhclf = pska(nsk,6)/(pska(nsk,1)+pska(nsk,6))
      rskc = (1-rhclf)*rskf
      drskc = (1-rhclf)*vskf
      ask(3) = rhclf*rskf
      dask(3) = rhclf*vskf
    end if
    dthetaskc = (pska(jt+1,2)-pska(jt,2)) / (tska(jt+1)-tska(jt))
    thetaskc = pska(jt,2) + (tsh-tska(jt))*dthetaskc
    dphiskc = (pska(jt+1,3)-pska(jt,3)) / (tska(jt+1)-tska(jt))
    phiskc = pska(jt,3) + (tsh-tska(jt))*dphiskc
    dgmsk = (pska(jt+1,7)-pska(jt,7)) / (tska(jt+1)-tska(jt))
    gmsk = pska(jt,7) + (tsh-tska(jt))*dgmsk
    if (tsh < tska(nsk)) then
      dask(1) = (pska(jt+1,4) - pska(jt,4)) / (tska(jt+1) - tska(jt))
      dask(2) = (pska(jt+1,5) - pska(jt,5)) / (tska(jt+1) - tska(jt))
      ask(1) = pska(jt,4) + (tsh-tska(jt)) * dask(1)
      ask(2) = pska(jt,5) + (tsh-tska(jt)) * dask(2)
    else
      rac = pska(nsk,4)/pska(nsk,6)
      ask(1) = rhclf*rskf*rac
      dask(1) = rhclf*vskf*rac
      rbc = pska(nsk,5)/pska(nsk,6)
      ask(2) = rhclf*rskf*rbc
      dask(2) = rhclf*vskf*rbc
    end if

    !  thetaskc is colatitude
    sinthetask = sin(thetaskc)
    costhetask = cos(thetaskc)
    sinphisk = sin(phiskc)
    cosphisk = cos(phiskc)
    r2xsk = mrtx(sinthetask, costhetask, sinphisk, cosphisk)
    !  r2xsk matrix to convert r-the-phi of shock ellipsoid center to xyz of
    !    of HEEQ. r=z theta=x, phi=y
    xp2r = mxptr(gmsk)

    xp2x = matmul(r2xsk, xp2r)

    if (mod(n,2) == 0) then
      ep = -1e-6
    else
      ep = 1e-6
    end if

    xpc(1) = (ep+ask(1)) * sin(theta) * cos(phi)
    xpc(2) = (ep+ask(2)) * sin(theta) * sin(phi)
    xpc(3) = (ep+ask(3)) * cos(theta)

    xc(1:3) = matmul(xp2x, xpc)

    xskc(1) = rskc * sinthetask * cosphisk
    xskc(2) = rskc * sinthetask * sinphisk
    xskc(3) = rskc * costhetask

    x(1:3) = xc(1:3) + xskc(1:3)

    print "(4a, 9e16.6e3)", __FILE__, ":", __LINE__, ":", xskc(:), thetaskc, phiskc, x(1:3), tsh
    return
  end subroutine
end program
