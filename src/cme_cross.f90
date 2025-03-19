module cme_cross
  use iso_fortran_env, only: real64
  use param, only: PI, TWOPI
  use sim3d_utils, only: solarwindtemp
  use mtrx, only: drvbmag, solarwind, mrtx, mxptr, dmxptr, norm2, dmrtx
  use datetime_utils, only: modjulianday
  use file_op, only: open_file_from_environment
  implicit none
contains
  !  read in parameters for CME shock ellipsoid in Kwon's model
  subroutine preparecme

    real(kind=real64)  :: pexsk(0:144,8), vskf0, acsk, tska(20), pska(20,8)
    integer :: nsk
    common /cmesk/pexsk,vskf0,acsk,tska,pska,nsk

    real(kind=real64)  :: vskfska(20)
    integer :: jt

    character(len=256) :: dir
    common/dir/dir

    real(kind=real64)  :: prm(9)

    real(kind=real64)  :: densw0, vsw, k4ok2, k6ok2, omega, b1au, vom, facip
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip

    real(kind=real64)  :: t0org, te, tdl, dmapjul, tcme0
    common /tmprm/t0org,te,tdl,dmapjul,tcme0

    real(kind=real64)  :: vsw_real, rfnsk

    real(kind=real64)  :: vsksw, tauf, tauc1_0, tauc2, tauc2_0, vcme0kmPs, vcme0 !LC
    common /vsksw/vsksw,tauf,tauc1_0,tauc2,tauc2_0,vcme0kmPs,vcme0 !LC

    real(kind=real64)  :: a_tau, c_tau, rc1, rc1sh_3d(3), rc1_3d(3), r1au(3)
    real(kind=real64)  :: bv(3), bm, cvtu(3), gbmag(3), bxgb2(3)
    real(kind=real64)  :: dbbds, pol, b1s, gb1s(3)
    real(kind=real64)  :: vpl(3), gvpl(3,3), densw
    integer :: ntauc2
    real(kind=real64)  :: tsh1, ca1, cs1, ca2, cs2, tpsw1
    real(kind=real64), parameter :: RS_PER_MIN_TO_KM_PER_SEC = 6.96340d5 / 60.0
    real(kind=real64), parameter :: gamma_cs = 1.6666666666666666666d0
    character(len=*), parameter :: &
      writefmt = "(10e18.10)", &
      !           YYYY - MM -  DD  T   HH :  MM :  SS.SSS
      readfmt  = "(i4,a1,i2,a1,i2, a1, i2,a1,i2,a1,f6.3,  3e15.5,4e18.7,f15.6,f15.3)", &
      readfmt2 = "(i4,a1,i2,a1,i2, a1, i2,a1,i2,a1,f6.3)"

    ! the following are dummy variables. they are to read in the hyphens and dashes
    ! that a datetime string might have, since fortran throws a hissy fit about
    ! having constant strings in read input formats
    character :: hyphen1, hyphen2, bigT, colon1, colon2


    integer :: iyr, imon, iday, ihr, imin, i, ii, jj
    real(kind=real64)  :: sec, fracday, hwinf, tmap0, tcme, alinf
    real(kind=real64)  :: cm1sq, costhetahbv2, dal
    real(kind=real64)  :: dphiskc, drskf, dso, dthetaskc, phic, phiskc, rskf
    real(kind=real64)  :: tauc1, tcme1
    real(kind=real64)  :: tempsw, thc, theta_skc, tsh
    real(kind=real64)  :: vthc, vphic, vphic0, vthc0, vskf, vsw0
    integer :: cme_fileunit, cpm_fileunit, ecme_fileunit

    call open_file_from_environment("CME_DATA_FILE", cme_fileunit, "FORMATTED")
    !open(newunit=cme_fileunit, file=trim(dir)//'cme.dat')

    read(cme_fileunit,readfmt2) iyr, hyphen1, imon,hyphen2, iday, bigT, ihr, colon1, imin,colon2, sec
    fracday = (ihr + imin/60.0 + sec/3600.0) / 24.0
    dmapjul = modjulianday(iyr, imon, iday, fracday)
    tmap0 = dmapjul * 1440
    read(cme_fileunit, *) nsk, tauf, hwinf !tauf rising phase time;minutes case by case
    read(cme_fileunit, *)
    hwinf = hwinf/180.0*PI
    do i = 1, nsk
      if (nsk > 20) stop 'nsk too large'
      read(cme_fileunit, readfmt) &
        iyr, hyphen1, imon, hyphen2, iday, bigT, &
        ihr, colon1, imin, colon2, sec, prm
      fracday = (ihr+imin/60.+sec/3600.) / 24.0
      tcme = modjulianday(iyr, imon, iday, fracday)
      tcme = tcme * 1440
      tska(i) = tcme - tmap0

      pska(i,1) = prm(4)/6.96340d8
      pska(i,2) = acos(prm(3) / norm2(prm))
      pska(i,3) = atan2(prm(2), prm(1))

      ! transform to HEEQ+60 longitude at map time
      pska(i,3) = pska(i,3) - (omega-2.0d0*PI/5.256e5)*tska(i) + PI/3
      pska(i,4) = (prm(5)+prm(6))/2/6.96340d8
      pska(i,5) = pska(i,4)
      pska(i,6) = prm(7)/6.96340d8
      pska(i,7) = 0.0
      pska(i,8) = prm(9)/180*PI
    end do
    close(cme_fileunit)

    write(*,*)'tska1,tska2',tska(1:2)
    write(*,*)'pska',pska(1,6),pska(2,6),pska(1,1),pska(2,1)
    do i = 2, nsk
      vskfska(i) = (pska(i,6) + pska(i,1) - pska(i-1,6) - pska(i-1,1)) / &
                   (tska(i) - tska(i-1))
      vskfska(i) = vskfska(i) * RS_PER_MIN_TO_KM_PER_SEC !km/s
    end do
    !  find the zero time of CME shock size when pska(6) is extrapolated to zero
    tcme0 = tska(1) - pska(1,6) / (pska(2,6)-pska(1,6)) * (tska(2)-tska(1))
    !  find the zero time of CME front when pska(6)+pska(1) is one
    tcme1 = tska(1) - (pska(1,6) + pska(1,1) - 1.0) / &
                      (pska(2,6) + pska(2,1) - pska(1,6)  -pska(1,1)) * &
                      (tska(2) - tska(1))
    if (tcme0 < tcme1) tcme0 = tcme1
    !  cme shock slow down last the last frame
    !  cme shock goes to vsw asymptotically
    !  cme shock pass 1 AU t1au since tska(nsk)
    rfnsk = pska(nsk,6)+pska(nsk,1) !LC last shock front
    !vsw_real = vsw/(1+k4ok2/rfnsk**2+k6ok2/rfnsk**4)!LC
    ! average shock front speed
    vskf0 = (pska(nsk,6)+pska(nsk,1)-pska(1,6)-pska(1,1)) / (tska(nsk)-tska(1))

    !LC
    !    1.sum(vskfska(2:nsk))/(nsk-1) average speed from cme.dat 2.vskfska(nsk) last data point 3. observarion
    vcme0kmPs = sum(vskfska(2:nsk))/(nsk-1)!vskf0*6.96340d5/60.!991.!km/s case by case
    vcme0 = vcme0kmPs/RS_PER_MIN_TO_KM_PER_SEC!Rs/min
    c_tau = 3380.*(tauf/60.)**(-1.14)/(vcme0kmPs-vsw*RS_PER_MIN_TO_KM_PER_SEC)
    a_tau = vcme0/vsw+1.0/sqrt(c_tau)*(vcme0/vsw-1.)
    tauc1 = a_tau*(1.0+sqrt(c_tau))/(a_tau-1.)*tauf !min !critical time for CME leading edge
    rc1 = pska(1,6)+pska(1,1)+vcme0*tauc1 !critical r for CME leading edge
    !rc1 = pska(1,1)+vcme0*tauc1 !critical r for CME leading edge
    if (tauc1 + tska(1) < tska(nsk)) then
      jt = locate(tska, nsk, tauc1 + tska(1))
    else
      jt = nsk-1
    end if
    drskf = (pska(jt+1,1)+pska(jt+1,6)-pska(1,1)-pska(1,6))/&
      (tska(jt+1)-tska(1))
    rskf = pska(jt,1)+pska(jt,6)+(tauc1+tska(1)-tska(jt))*drskf


    !dso1 = pska(nsk,1)-pska(nsk,6) !standoff distance ! temporary
    dthetaskc = (pska(jt+1,2)-pska(1,2))/(tska(jt+1)-tska(1))
    theta_skc = pska(jt,2)+(tauc1+tska(1)-tska(jt))*dthetaskc
    dphiskc = pska(jt+1,3) - pska(jt,3)
    if (dphiskc >  PI) dphiskc = dphiskc - TWOPI
    if (dphiskc < -PI) dphiskc = dphiskc + TWOPI
    dphiskc = dphiskc / (tska(jt+1)-tska(jt))
    phiskc = pska(jt,3) + (tauc1+tska(1)-tska(jt))*dphiskc


    !!=======if using data or observation at 1au to calculate Va and Vs=======
    !!bm densw,tempsw could be replaced with those from observation
    !r1au(1) = 215.032
    !r1au(2) = theta_skc

    !tpsw = ((r1au(1)-k4ok2/r1au(1)-k6ok2/r1au(1)**3/3)
    !&   -(rc1-k4ok2/rc1-k6ok2/rc1**3/3))/vsw
    ! r1au(3) = phiskc-omega*tpsw
    ! call drvbmag(r1au, bv, bm, cvtu, gbmag, bxgb2, dbbds, pol, b1s, gb1s)
    ! call solarwind(r1au, vpl, gvpl, densw)
    ! va1au = 187.8*bm/sqrt(densw)
    ! tempsw = solarwindtemp(r1au)
    ! vs1au = 7.83e-6 * sqrt(gamma_cs*tempsw)
    ! ca1 = va1au*r1au(1)/rc1*
    !&  sqrt((vsw+(rc1-2.)*omega)/(vsw+(r1au(1)-2.)*omega))
    ! cs1 = vs1au*(r1au(1)/rc1)**(gamma_cs-1)
    ! costhetahbv2 = vsw**2/(vsw**2+((r1au(1)-2.)*omega)**2)
    !!=======if using data at rc1 to calculate Va and Vs directly=====
    vsw0 = vsw

    rc1_3d(1) = rc1
    rc1_3d(2) = theta_skc
    rc1_3d(3) = phiskc
    call drvbmag(rc1_3d, bv, bm, cvtu, gbmag, bxgb2, dbbds, pol, b1s, gb1s)
    call solarwind(rc1_3d, vpl, gvpl, densw)
    ca1 = 187.8*bm/sqrt(densw)
    tempsw = solarwindtemp(rc1_3d)
    cs1 = 7.83e-6*sqrt(gamma_cs*tempsw)
    ! cm1sq: square of mach number at rc1
    !!magnetosonic speed (Vms): Vms**2 = 0.5 * {Va**2 + Vs**2 + [(Va**2 + Vs**2)**2 -
    !! 4 * Va**2 * Vs**2 * (cos (theta))**2]**0.5}
    costhetahbv2 = vpl(1)**2 / (vpl(1)**2+vpl(2)**2)
    cm1sq = 2. * (vcme0-vsw)**2 / (ca1**2+cs1**2+&
      sqrt((ca1**2 + cs1**2)**2 - 4 * ca1**2 * cs1**2 * costhetahbv2))
    r1au(1) = 215.032
    dso = 0.264*r1au(1)*((gamma_cs-1)*cm1sq+2)/((gamma_cs+1)*(cm1sq-1))*&
      (rc1/r1au(1))**0.78 !standoff distance

    rc1sh_3d(1) = rc1+dso/2!sheath position at CME leading edge critical time
    rc1sh_3d(2) = theta_skc
    tpsw1 = ((rc1sh_3d(1)-k4ok2/rc1sh_3d(1)-k6ok2/rc1sh_3d(1)**3/3)&
      -(rc1-k4ok2/rc1-k6ok2/rc1**3/3))/vsw
    rc1sh_3d(3) = phiskc-omega*tpsw1
    call drvbmag(rc1sh_3d, bv, bm, cvtu, gbmag, bxgb2, dbbds, pol, b1s, gb1s)
    call solarwind(rc1sh_3d, vpl, gvpl, densw)
    ca2 = 187.8*bm/sqrt(densw)
    tempsw = solarwindtemp(rc1sh_3d)
    cs2 = 7.83e-6*sqrt(gamma_cs*tempsw)! bm densw tempsw from model (or observation)
    !cm1sq = (vcme0-vsw)**2/cs2**2

    tauc2 = dso/sqrt(ca2**2+cs2**2)+tauc1 !min !critical time for CME shock front
    tauc2_0 = tauc2 + tska(1)
    tauc1_0 = tauc1 + tska(1)
    !!====================================================
    call open_file_from_environment(&
      "CME_PROP_MODEL_FILE", cpm_fileunit, "FORMATTED")
    !open (newunit=cpm_fileunit, file = trim(dir)//'cmepropmodel.dat')
    write(cpm_fileunit,'(A)')'tauf(min): rising phase time '
    write(cpm_fileunit,'(A)')'tcme0(min): zero time of CME shock'
    write(cpm_fileunit,'(A)')'tauc1_0:critical time of CME leading edge'
    write(cpm_fileunit,'(A)')'tauc2_0:critical time of CME shock front'
    write(cpm_fileunit,'(A)')'rc1:position at critical time of CME leading edge'
    write(cpm_fileunit,'(A)')'dso(RS): standoff distance'
    write(cpm_fileunit,'(A)')'ca1,cs1(Rs/min): Alfven and sound speed at rc1 '
    write(cpm_fileunit,'(A)')'ca2,cs2(Rs/min): Alfven and sound speed at the &
      & shealth rc1+dso/2'
    write(cpm_fileunit,'(A)')'=============================================='
    write(cpm_fileunit,'(A)')'Calculated variables:'
    write(cpm_fileunit,'(A,F5.1)')'tauf(min)=',tauf
    write(cpm_fileunit,'(A,3F8.1)')'rc1(Rs,deg,deg)=',&
      rc1_3d(1),rc1_3d(2)/PI*180.,rc1_3d(3)/PI*180.
    write(cpm_fileunit,'(A,F5.1)')'dso(Rs)=',dso
    write(cpm_fileunit,'(3(A,F15.1,2x))') &
      'tcme0=',tcme0,', tauc1_0=',tauc1_0,', tauc2_0=',tauc2_0
    write(cpm_fileunit,'(4(A,F8.4,2x))')'CA1=',&
      ca1,',CS1=',cs1,',CA2=',ca2,',CS2=',cs2
    write(cpm_fileunit,'(A,F5.1)')'Mach number at rc1=', sqrt(cm1sq)
    write(cpm_fileunit,'(A,F7.2)')'Half Angle at 21.5 Rs= ', hwinf*57.29578
    write(cpm_fileunit,'(A)')'=============================================='
    !==============comment this part if just using tauc2 from the new model====
    !----------------find tauc2 from cme.dat----------------------
    ntauc2 = nsk
    do i = (nsk-2), 2, -1
      if (vskfska(i+1) < vskfska(i) .and. vskfska(i+2) < vskfska(i+1)) then
        ntauc2 = i
      end if
    end do
    !!----------choose smaller one-----------------------------------
    ! if (tauc2_0>tska(ntauc2)) then
    !     write(cpm_fileunit,'(2(A,F8.1,2x))')&
    !        'Calculated tauc2_0=', tauc2_0, &
    !        ' is larger than tauc2_0 from cme.dat', tska(ntauc2)
    !     write(cpm_fileunit,'(A)')'tauc2_0 from cme.dat is adopted'
    !     tauc2_0 = tska(ntauc2)+0.01
    ! else
    !     write(cpm_fileunit,'(A)')'tauc2_0 from shock model is adopted'
    !  end if
    ! write(cpm_fileunit,'(A,F8.3)')'adopted tauc2_0 is',tauc2_0
    !---ignore cme.data where tska>tauc2_0 if tauc2_0<tska(nsk)----------------
    !!if (tauc2_0<tska(nsk)) then
    !!write(cpm_fileunit,'(A)')'=============================================='
    !!write(cpm_fileunit,'(A)')'tauc2_0 is less than tska(nsk).'
    !!write(cpm_fileunit,'(A,i3)')'original nsk=',nsk
    !!jt = locate(tska,nsk,tauc2_0)
    !!if (jt>1) then
    !!  nsk = jt !update nsk
    !!else
    !!  nsk = 2
    !!  tauc2_0 = tska(2)+0.01   !update nsk and tauc2_0
    !!end if
    !!write(cpm_fileunit,'(A,i3)')'ignore cme.dat where nsk>',nsk
    !!write(cpm_fileunit,'(A,i3)')'updated nsk=',nsk
    !!write(cpm_fileunit,'(A)')'=============================================='
    !!end if
    !===========================================================
    !---------------- calculate t,vskf at 1au------------------
    rskf = 0.0
    vskf0 = (pska(nsk,6)+pska(nsk,1)-pska(1,6)-pska(1,1)) / (tska(nsk)-tska(1))
    vcme0 = vskf0
    tsh1 = tauc2_0

    do while (rskf < 215.032)
      if (tauc2_0 > tska(nsk)) then
        rskf = 3.0/2.0*(vcme0-vsw0)*tauc2*(((tsh1-tska(1))/tauc2)&
          **(2.0/3.0)-1.0d0) + vsw0*tauc2*(((tsh1-tska(1))/tauc2)-1.0d0)&
          + pska(nsk,1)+pska(nsk,6)+vcme0*(tauc2_0-tska(nsk)) !LC
      else
        rskf = 3.0/2.0*(vcme0-vsw0)*tauc2*(((tsh1-tska(1))/tauc2)**(2.0/3.0&
          )-((tska(nsk)-tska(1))/tauc2)**(2.0/3.0))+vsw0*(tsh1-tska(nsk))&
          + pska(nsk,1)+pska(nsk,6)
      end if
      vskf = (vcme0-vsw0) * (((tsh1-tska(1))/tauc2)**(-1.0d0/3.0d0)) + vsw0
      tsh1 = tsh1+0.5
    end do
    write(cpm_fileunit,'(A,4F8.1)')'t1au(min),vcme0,vskf1au,vsw0(km/s)=', &
      tsh1, vcme0kmPs, vskf*RS_PER_MIN_TO_KM_PER_SEC, &
      vsw0*RS_PER_MIN_TO_KM_PER_SEC

    !call flush(cpm_fileunit)
    flush(cpm_fileunit)
    close(cpm_fileunit)

    !------------------------------------------------------------------------
    !  Use the propagation model to calculate the radial, lat and long
    !  motion and save to an array up to 3days
    pexsk(0,:) = pska(nsk,:)
    thc = pexsk(0,2)
    phic = 0.0 !relative to pska(nsk,3)
    dthetaskc = (pska(nsk,2)-pska(1,2)) / (tska(nsk)-tska(1))
    dphiskc = pska(nsk,3) - pska(nsk-1,3)
    if (dphiskc >  PI) dphiskc = dphiskc - TWOPI
    if (dphiskc < -PI) dphiskc = dphiskc + TWOPI
    dphiskc = dphiskc / (tska(nsk)-tska(nsk-1))
    vthc0 = (pska(nsk,1)+pska(nsk,6)) * dthetaskc
    vphic0 = (pska(nsk,1)+pska(nsk,6)) * (dphiskc+omega) * sin(thc)
    ! limiting shock size to hwinf (after 21.5 Rs)
    !    hwinf = pska(nsk,6) * sin(alinf) / (pska(nsk,6)*cos(alinf)+pska(nsk,1))
    alinf = hwinf
    dal = 1.0
    do while (dal*dal > 0.0001)
      dal = (hwinf*(pska(nsk,6)*cos(alinf)+pska(nsk,1))-pska(nsk,6)*sin(alinf))&
            / (pska(nsk,6)*cos(alinf)+hwinf*pska(nsk,6)*sin(alinf))
      alinf = alinf + dal
    end do
    do i = 1, 43200
      tsh1 = tska(nsk) + 0.1*i
      if (tsh1 < tauc2_0) then
        rskf = pska(nsk,1)+pska(nsk,6)+vcme0*(tsh1-tska(nsk))
      else
        if (tauc2_0 > tska(nsk)) then
          rskf = 3.0/2.0*(vcme0-vsw0)*tauc2*(((tsh1-tska(1))/tauc2)&
            **(2.0/3.0)-1.0d0) + vsw0*tauc2*(((tsh1-tska(1))/tauc2)-1.0d0)&
            + pska(nsk,1)+pska(nsk,6)+vcme0*(tauc2_0-tska(nsk)) !LC
        else
          rskf = 3.0/2.0*(vcme0-vsw0)*tauc2*(((tsh1-tska(1))/tauc2)**(2.0/3.0&
            )-((tska(nsk)-tska(1))/tauc2)**(2.0/3.0))+vsw0*(tsh-tska(nsk))&
            + pska(nsk,1)+pska(nsk,6)
        end if
      end if
      if (tsh1 < tauc2_0) then
        vthc = vthc0
        vphic = vphic0
      else
        vthc = (vthc0-0)*((tsh1-tska(1))/tauc2)**(-1.0/3.0)
        vphic = (vphic0-0)*((tsh1-tska(1))/tauc2)**(-1.0/3.0)
      end if
      thc = thc+vthc/rskf*0.1  !integration time step is 0.1 min
      phic = phic+vphic/rskf/sin(thc)*0.1 !fixed frame
      if (mod(i,300) == 0) then
        ii = i/300
        pexsk(ii,1) = pska(nsk,1)/(pska(nsk,1)+pska(nsk,6))*rskf
        pexsk(ii,2) = thc
        pexsk(ii,3) = phic+pska(nsk,3)-omega*(tsh1-tska(nsk))
        pexsk(ii,4) = pska(nsk,4)/(pska(nsk,1)+pska(nsk,6))*rskf
        pexsk(ii,5) = pska(nsk,5)/(pska(nsk,1)+pska(nsk,6))*rskf
        pexsk(ii,6) = pska(nsk,6)/(pska(nsk,1)+pska(nsk,6))*rskf
        pexsk(ii,7) = 0.
        if (pexsk(ii,6)*cos(alinf) + pexsk(ii,1) > 21.5) then
          pexsk(ii,8) = alinf
        else
          jj = ii
        end if
      end if
    end do
    do ii = 1,jj
      pexsk(ii,8) = pexsk(0,8) + (alinf-pexsk(0,8))/jj*ii
    end do

    call open_file_from_environment(&
      "EXTENDED_CME_DATA_FILE", ecme_fileunit, "FORMATTED")
    !open(newunit=ecme_fileunit, file=trim(dir)//'extended_cme.dat')
    write(ecme_fileunit,'(A)') 'Time Maptime, 8 perameters in Rs as in cme.dat'
    do i = 1, nsk
      write(ecme_fileunit,writefmt) tska(i), pska(i,1:8)
    end do
    do i = 0, 144
      write(ecme_fileunit,writefmt) tska(nsk)+i*30, pexsk(i,:)
    end do
    !call flush(ecme_fileunit)
    flush(ecme_fileunit)
    close(ecme_fileunit)
  end subroutine

  ! find which side of CME shock a particle is
  !  if crossed, calculate distance to shock, shock speed and normal
  subroutine inorout(tsh, x, dnsk, vsk, vnx)
    real(kind=real64), intent(in)     :: tsh, x(6)
    real(kind=real64), intent(in out) :: dnsk
    real(kind=real64), intent(out)    :: vsk, vnx(3)

    real(kind=real64) :: ra
    real(kind=real64) :: dnsk0 !copy of input
    !  dnsk = distance to shock outside +; inside -
    !  vsk = shock speed,
    !  vnx = normal to shock in xyz
    real(kind=real64) :: pexsk(0:144,8), vskf0, acsk, tska(20), pska(20,8)
    integer           :: nsk
    common /cmesk/pexsk,vskf0,acsk,tska,pska,nsk

    real(kind=real64) :: rskc, drskc, rskmax, rskmin, ask(3), dask(3)
    real(kind=real64) :: theta_skc, phiskc, gmsk, dthetaskc, dphiskc, dgmsk
    real(kind=real64) :: sinthetask, costhetask, sinphisk, cosphisk
    real(kind=real64) :: r2xsk(3,3), dr2xsk(3,3)
    real(kind=real64) :: xp2r(3,3), xp2x(3,3), dxp2r(3,3), dxp2x(3,3)
    real(kind=real64) :: xskc(3), dxskc(3), dxpskc(3), xc(3), xpc(3), rho, vnm, vnxp(3)
    real(kind=real64) :: grf1(3), grf2(3), grf3(3), grf2m, drxc(3)
    real(kind=real64) :: dot1, dot2, dot3, dot4
    real(kind=real64) :: densw0, vsw, k4ok2, k6ok2, omega, b1au, vom, facip
    common/bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip
    real(kind=real64) :: t0org, te, tdl, dmapjul, tcme0
    common /tmprm/t0org,te,tdl,dmapjul,tcme0
    real(kind=real64) :: vsksw, tauf, tauc1_0, tauc2, tauc2_0, vcme0kmPs, vcme0 !LC
    common /vsksw/vsksw,tauf,tauc1_0,tauc2,tauc2_0,vcme0kmPs,vcme0 !LC
    integer           :: i, j, jt, jt1
    real(kind=real64) :: dtzskmin, tzskmin, rac, rbc, rhclf, rskf, vskf, tsh1

    dnsk0 = dnsk

    if (tsh < tska(nsk)) then
      jt = locate(tska,nsk,tsh)
      if (jt == 0) jt = 1
      !  maximum shock radial distance for quick search
      drskc = (pska(jt+1,1)-pska(jt,1))/(tska(jt+1)-tska(jt))
      rskc = pska(jt,1)+(tsh-tska(jt))*drskc
      dask(3) = (pska(jt+1,6)-pska(jt,6))/(tska(jt+1)-tska(jt))
      ask(3) = pska(jt,6)+(tsh-tska(jt))*dask(3)
      !vskf = vskf0 !only for fort.58; useless
      !rskf = pska(jt,1)+pska(jt,6)!only for fort.58; useless
      dthetaskc = (pska(jt+1,2)-pska(jt,2))/(tska(jt+1)-tska(jt))
      theta_skc = pska(jt,2)+(tsh-tska(jt))*dthetaskc
      dphiskc = (pska(jt+1,3)-pska(jt,3))
      if (dphiskc >  PI) dphiskc = dphiskc - TWOPI
      if (dphiskc < -PI) dphiskc = dphiskc + TWOPI
      dphiskc = dphiskc / (tska(jt+1) - tska(jt))
      phiskc = pska(jt,3) + (tsh-tska(jt)) * dphiskc
    else
      jt = nsk-1
      vcme0 = vskf0
      if (tsh < tauc2_0) then
        rskf = pska(nsk,1)+pska(nsk,6)+vcme0*(tsh-tska(nsk)) !LC
        vskf = vcme0
      else
        if (tauc2_0>tska(nsk)) then
          rskf = 3.0/2.0 * (vcme0-vsw)*tauc2*(((tsh-tska(1))/tauc2)**(2.0/3.0)&
            -1.0d0)+vsw*tauc2*(((tsh-tska(1))/tauc2)-1.0d0)&
            + pska(nsk,1)+pska(nsk,6)+vcme0*(tauc2_0-tska(nsk)) !LC
        else
          rskf = 3.0/2.0 * (vcme0-vsw)*tauc2*(((tsh-tska(1))/tauc2)**(2.0/3.0)&
            -((tska(nsk)-tska(1))/tauc2)**(2.0/3.0)) + vsw*(tsh-tska(nsk))&
            + pska(nsk,1)+pska(nsk,6)
        end if
        vskf = (vcme0-vsw) * (((tsh-tska(1))/tauc2)**(-1.0d0/3.0d0)) + vsw
      end if

      rhclf = pska(nsk,6) / (pska(nsk,1) + pska(nsk,6))
      rskc = (1-rhclf) * rskf
      drskc = (1-rhclf) * vskf
      if ((pska(nsk,1) - pska(nsk-1,1)) < 1.e-3) then
        drskc = 0.0
        rskc = pska(nsk,1)
      end if
      ask(3) = rhclf * rskf
      dask(3) = rhclf * vskf
    end if
    rskmax = rskc + ask(3)
    rskmin = rskc - ask(3)
    ra = norm2(x(1:3))
    !write(58,"(f14.8,3(1pe12.4))") tsh,ra,rskf,vskf
    if ((rskmax < ra*.99) .or. (rskmin > ra*1.01)) then
      dnsk = 0.01*ra  ! use step size as a minimum
      if (dnsk0 * dnsk >= 0.0) then
        vsk = 1d-6
        vnx = 0.57735026918962576450d0
        return
      end if
    end if

    if (tsh < tska(nsk)) then
      dthetaskc = (pska(jt+1,2)-pska(jt,2)) / (tska(jt+1)-tska(jt))
      theta_skc = pska(jt,2)+(tsh-tska(jt))*dthetaskc
      dphiskc = pska(jt+1,3) - pska(jt,3)
      if (dphiskc >  PI) dphiskc = dphiskc - TWOPI
      if (dphiskc < -PI) dphiskc = dphiskc + TWOPI
      dphiskc = dphiskc/(tska(jt+1)-tska(jt))
      phiskc = pska(jt,3)+(tsh-tska(jt))*dphiskc
      dgmsk = (pska(jt+1,7)-pska(jt,7))/(tska(jt+1)-tska(jt))
      gmsk = pska(jt,7)+(tsh-tska(jt))*dgmsk
      dask(1) = (pska(jt+1,4)-pska(jt,4))/(tska(jt+1)-tska(jt))
      ask(1) = pska(jt,4)+(tsh-tska(jt))*dask(1)
      dask(2) = (pska(jt+1,5)-pska(jt,5))/(tska(jt+1)-tska(jt))
      ask(2) = pska(jt,5)+(tsh-tska(jt))*dask(2)
      dtzskmin = (pska(jt+1,8)-pska(jt,8)) / (tska(jt+1)-tska(jt))
      tzskmin = pska(jt,8)+(tsh-tska(jt))*dtzskmin
    else
      tsh1 = tsh - tska(nsk)
      jt1 = floor(tsh1/30)
      dthetaskc = (pexsk(jt1+1,2)-pexsk(jt1,2))/30
      theta_skc = pexsk(jt1,2)+(tsh1-30*jt1)*dthetaskc
      dphiskc = pexsk(jt1+1,3)-pexsk(jt1,3)
      if (dphiskc >  PI) dphiskc = dphiskc - TWOPI
      if (dphiskc < -PI) dphiskc = dphiskc + TWOPI
      dphiskc = dphiskc / 30
      phiskc = pexsk(jt1,3) + (tsh1-30*jt1) * dphiskc
      dgmsk = (pexsk(jt1+1,7)-pexsk(jt1,7)) / 30
      gmsk = pexsk(jt1,7) + (tsh1-30*jt1)*dgmsk
      rac = pska(nsk,4) / pska(nsk,6)
      ask(1) = rhclf * rskf * rac
      dask(1) = rhclf * vskf * rac
      rbc = pska(nsk,5) / pska(nsk,6)
      ask(2) = rhclf * rskf * rbc
      dask(2) = rhclf * vskf * rbc
      dtzskmin = (pexsk(jt1+1,8)-pexsk(jt1,8)) / 30
      tzskmin = pexsk(jt1,8) + (tsh1-30*jt1) * dtzskmin
    end if
    tzskmin = cos(tzskmin)


    !tzskmin = -0.866
    ! write(*,*) rskc, drskc, ask, dask
    !  theta_skc is colatitude
    sinthetask = sin(theta_skc)
    costhetask = cos(theta_skc)
    sinphisk = sin(phiskc)
    cosphisk = cos(phiskc)
    r2xsk = mrtx(sinthetask, costhetask, sinphisk, cosphisk)
    !  r2xsk matrix to convert r-the-phi of shock ellipsoid center to xyz of
    !    of HEEQ. r = z theta = x, phi = y
    call mxptr(gmsk, xp2r)

    !xp2x = matmul(r2xsk, xp2r)
    do i = 1, 3
      xp2x(i,1:3) = r2xsk(i,1)*xp2r(1,1:3) + r2xsk(i,2)*xp2r(2,1:3) + &
          r2xsk(i,3)*xp2r(3,1:3)
    end do

    xskc(1) = rskc * sinthetask * cosphisk
    xskc(2) = rskc * sinthetask * sinphisk
    xskc(3) = rskc * costhetask

    xc(1:3) = x(1:3) - xskc(1:3)
    xpc(1:3) = xc(1)*xp2x(1,1:3)+xc(2)*xp2x(2,1:3)+xc(3)*xp2x(3,1:3)

    !  distance to the shock
    grf1(1:3) = xpc(1:3)/ask(1:3)
    rho = norm2(grf1(1:3))
    grf2(1:3) = grf1(1:3) / ask(1:3)
    grf3(1:3) = grf1(1:3) * grf1(1:3)
    grf2m = norm2(grf2(1:3))
    dnsk = (rho-1) / grf2m

    if (dnsk0*dnsk < 0.) then ! crossed shock
      if (grf1(3) < tzskmin) then !No shock on the backside
        vsk = 1d-6
        vnx = 0.57735026918962576450d0
        return
      end if
      !   shock normal
      vnxp = grf2/grf2m
      vnx = matmul(xp2x, vnxp)
      !  calculate shock velocity (normal)
      dr2xsk = dmrtx(sinthetask, costhetask, sinphisk, cosphisk, dthetaskc, dphiskc)
      dxp2r = dmxptr(gmsk, dgmsk)
      do i = 1, 3
        do j = 1, 3
          dxp2x(i,j) = dot_product(dr2xsk(i,1:3), xp2r(1:3,j)) + &
                       dot_product(r2xsk(i,1:3), dxp2r(1:3,j))
        end do
      end do
      dot1 = dot_product(grf2, vnxp)
      dot2 = dot_product(grf3, dask(:)/ask(:))
      dxskc(1) = drskc * sinthetask * cosphisk &
               + rskc * costhetask * dthetaskc * cosphisk &
               - rskc * sinthetask * sinphisk * dphiskc
      dxskc(2) = drskc * sinthetask * sinphisk &
               + rskc * costhetask * dthetaskc * sinphisk &
               + rskc * sinthetask * cosphisk * dphiskc
      dxskc(3) = drskc * costhetask - rskc * sinthetask * dthetaskc
      do i = 1, 3
        dxpskc(i) = dot_product(dxskc, xp2x(:,i))
      end do
      dot3 = dot_product(grf2, dxpskc)
      do i = 1, 3
        drxc(i) = dot_product(xc, dxp2x(1:3,i))
      end do
      dot4 = dot_product(grf2(1:3), drxc(1:3))
      vsk = (dot2 + dot3 - dot4) / dot1
    else
      vsk = 1.d-6
      vnx = 0.57735026918962576450d0
    end if
    !write(60,"(f14.6,11(1pe12.4))") tsh, ra, rskc, vsk, vnx, ask, rskf, vskf
    return
  end subroutine

  integer function locate(xx, n, x) result(j)
    integer, intent(in)           :: n
    real(kind=real64), intent(in) :: xx(n)
    real(kind=real64), intent(in) :: x

    integer              :: jl, jm, ju

    jl = 0
    ju = n + 1

10  if (ju - jl > 1) then
      jm = (ju+jl) / 2
      if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
        jl = jm
      else
        ju = jm
      end if
      goto 10
    end if
    j = jl
  end function

end module
