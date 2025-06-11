program mapb2s
  use iso_fortran_env, only: real64, int64
  use param, only: PI, N_R, N_THETA, N_PHI, TWOPI, DEG_TO_RAD
  use mtrx, only: norm2, trilinear
  use file_op, only: read_maggrid, write_b1rs
  use rksolvers, only: rk4

  implicit none

  include 'omp_lib.h'

  real(kind=real64), parameter :: RSS = 2.5d0

  real(kind=real64) :: magfieldgrid(0:N_R, 0:N_THETA, 0:N_PHI, 3)
  real(kind=real64) :: gbgrid(0:N_R, 0:N_THETA, 0:N_PHI, 3) ! dummy storage
  real(kind=real64) :: b1rs(0:N_R, 0:N_THETA, 0:N_PHI, 2)
  real(kind=real64) :: pol
  real(kind=real64) :: bmag, bmag0
  !real(kind=real64) :: v(3)
  real(kind=real64) :: t, dt
  real(kind=real64) :: r0(3), r1(3), r(3), b(3)
  !real(kind=real64) :: sintheta
  !real(kind=real64) :: nan(6), nan1
  !real(kind=real64) :: dis
  real(kind=real64) :: dr
  !real(kind=real64) :: bm
  integer           :: n, i, j, k
  integer           :: chunk
  !integer           :: id
  !real(kind=real64) :: ratio
  integer           :: lr, ltheta, lphi
  real(kind=real64) :: rmin(3), rmin1(3)
  integer           :: map(0:N_R, 0:N_THETA, 0:N_PHI)

  integer           :: n1
  !real(kind=real64) :: r_mhd, theta_mhd, phi_mhd

  t = 0d0
  dr = 0.01d0
  dt = 0.005d0
  chunk = 1

  !  load magnetic field grid
  call read_maggrid(magfieldgrid, gbgrid)

  b1rs = 0.d0
  map = 0

  print "('N_R = ',i0,', N_THETA = ',i0,', N_PHI = ',i0)", N_R, N_THETA, N_PHI

  !$OMP PARALLEL NUM_THREADS(20) DEFAULT(firstprivate) SHARED(chunk,magfieldgrid,b1rs,map)
  !$OMP DO SCHEDULE(static,chunk)
  !  trace field lines from surface surface back to 1 Rs

  do k = 0, N_R
    r0(1) = 1.00001d0 + k*0.01d0
    do i = 0, N_THETA
      r0(2) = i * DEG_TO_RAD
      do j = 0, N_PHI
        r0(3) = j * DEG_TO_RAD
        r = r0
        b = magfield(r)
        if (b(1) > 0) then
          pol = 1.0
        else
          pol = -1.0
        end if
        n = 0

        bmag0 = norm2(b)
        b1rs(k,i,j,2) = bmag0
        bmag = bmag0
        if (k == 1 .and. i == 1 .and. j == 1) write(*,*) b
        rmin(1) = 100000.d0
        do while (r(1) > 1.000009d0)
          r(2) = acos(cos(r(2)))
          r(3) = atan2(sin(r(3)), cos(r(3)))
          if (r(3) < 0) r(3) = r(3) + TWOPI
          r1 = r
          b = magfield(r)
          bmag = norm2(b)
          !v = -pol*b/bmag
          !v(2) = v(2)/r1(1)
          !sintheta = sin(r1(2))
          !if (sintheta == 0) sintheta = 1.0d-6
          !v(3) = v(3) / (r1(1)*sintheta)
          dt = dr * 0.1
          !r = rk4(r1, v, 3, t, dt, vfunc, pol)
          r = rk4(vfunc, t, r1, dt, pol)
          n = n + 1

          if (r(2) < 0) then
            r(2) = -r(2)
            r(3) = r(3) + pi
            if (r(3) > twopi) r(3) = r(3) - twopi
          end if

          if (r(1) < rmin(1)) rmin = r
          if ((r(1) - rmin(1)) > 2.5 .or. n > 10000) then !wrong direction?
            pol = -pol !try back with opposite polarity
            n1 = 0
            r = r0
            rmin1(1) = 100000.0
            do while (r(1) > 1.000009d0)
              r(2) = acos(cos(r(2)))
              r(3) = atan2(sin(r(3)),cos(r(3)))
              if (r(3) < 0) r(3) = r(3) + twopi
              r1 = r
              b = magfield(r)
              bmag = norm2(b)
              !v = -pol * b / bmag
              !v(2) = v(2) / r1(1)
              !sintheta = sin(r1(2))
              !if (sintheta == 0) sintheta = 1.0d-6
              !v(3) = v(3) / (r1(1)*sintheta)
              dt = dr * 0.1
              !r = rk4(r1, v, 3, t, dt, vfunc, pol)
              r = rk4(vfunc, t, r1, dt, pol)
              n1 = n1 + 1
              if (r(2) < 0) then
                r(2) = -r(2)
                r(3) = r(3)+pi
                if (r(3) > twopi) r(3) = r(3) - twopi
              end if
              if (r(1) < rmin1(1)) rmin1 = r

              ! double open field or stuck
              if ((r(1) - rmin1(1)) > 2.5 .or. n1 > 10000) then
                if (rmin(1) < rmin1(1)) then
                  pol = -pol
                  r = rmin ! use rmin to remap
                else
                  pol = pol
                  r = rmin1
                end if
                lr = floor((r(1)-1.0d0) / 0.01d0)
                ltheta = floor(r(2) / pi * 180.d0)
                lphi = floor(r(3) / pi * 180.d0)
                if (((lr*(N_THETA+1)+ltheta) * (N_PHI+1) + lphi) < &
                    ((k*(N_THETA+1) + i) * (N_PHI+1) + j)) then
                  map(k,i,j) = (lr*(N_THETA+1)+ltheta) * (N_PHI+1) + lphi
                else
                  write(*,*) "New x point at", lr, ltheta, lphi
                  map(k,i,j) = ((lr-1)*(N_THETA+1)+ltheta) * (N_PHI+1) + lphi
                end if
                bmag = 0.0
                exit
              end if
            end do
            exit
          end if
        end do

        !nan1 = -3.0
        !nan = sqrt(nan1)
        b1rs(k,i,j,1) = pol * bmag
        map(k,i,j) = pol * map(k,i,j)
      end do
    end do
  end do

  !$OMP BARRIER
  !$OMP END PARALLEL

  call write_b1rs(b1rs, map)

contains

  ! velocity function. find dr/dt given r and t (t not really needed)
  ! v = -sgn(B_r) / |B| (B_r r_hat + B_theta / r theta_hat + B_phi / (r sin theta) phi_hat)
  function vfunc(t, r, pol) result(v)
    real(kind=real64), intent(in)  :: t, r(3), pol
    real(kind=real64)              :: v(3)
    real(kind=real64)              :: b(3), bmag, sintheta

    b = magfield(r)
    bmag = norm2(b)
    if (bmag == 0) then
      v = 0
      return
    end if
    v = -pol * b / bmag
    v(2) = v(2) / r(1)
    sintheta = sin(r(2))
    if (sintheta == 0) sintheta = 1.0d-6
    v(3) = v(3) / (r(1) * sintheta)
  end function

  ! get \vec{B} = \vec{B}(\vec{r})
  function magfield(rvec) result(b)
    real(kind=real64), intent(in)  :: rvec(3)
    real(kind=real64)              :: b(3)

    real(kind=real64)              :: r, theta, phi ! input copy
    real(kind=real64)              :: fc(2,2,2)
    real(kind=real64)              :: px(3)

    integer(kind=int64)            :: irr, itheta, iphi ! grid positon
    integer(kind=int64)            :: m ! iteration variable

    b(:) = 0

    r = rvec(1)
    theta = rvec(2)
    phi = rvec(3)
    theta = acos(cos(theta))
    phi = atan2(sin(phi), cos(phi))
    if (phi < 0.0) phi = phi + twopi

    !  find the grid cell
    irr = floor((r-1) / (RSS-1) * N_R)
    if (irr >= N_R) irr = N_R - 1

    itheta = floor(theta / pi * N_THETA)
    if (itheta >= N_THETA) itheta = N_THETA - 1

    iphi = floor(phi / twopi * N_PHI)
    if (iphi >= N_PHI) iphi = N_PHI - 1

    if (irr <   0) then
      ! going into the sun, stop at the surface
      b = magfieldgrid(0, itheta, iphi, :)
      return
    end if

    !  relative displacement from lower grids
    px(1) = (r-1.0) / (RSS-1.0) * N_R - irr
    px(2) = theta / pi * N_THETA - itheta
    px(3) = phi / twopi * N_PHI - iphi

    !print "(a,':',i3,': ', ' irr = ', i3, ', itheta = ', i3, ', iphi = ', i3)", &
    !  __FILE__, __LINE__, irr, itheta, iphi

    do m = 1, 3
      fc(1,1,1) = magfieldgrid(irr,   itheta,   iphi,   m)
      fc(2,1,1) = magfieldgrid(irr+1, itheta,   iphi,   m)
      fc(1,2,1) = magfieldgrid(irr,   itheta+1, iphi,   m)
      fc(2,2,1) = magfieldgrid(irr+1, itheta+1, iphi,   m)
      fc(1,1,2) = magfieldgrid(irr,   itheta,   iphi+1, m)
      fc(2,1,2) = magfieldgrid(irr+1, itheta,   iphi+1, m)
      fc(1,2,2) = magfieldgrid(irr,   itheta+1, iphi+1, m)
      fc(2,2,2) = magfieldgrid(irr+1, itheta+1, iphi+1, m)

      ! interpolate the magnetic field
      b(m) = trilinear(fc, px)
    end do
  end function ! magfield

end program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Modification log
!!7/12/21:
!     corona data interpolated e.g. l246
!     adaptive dt based on dr e.g. l70
!     use b1rs from ir-1 if reconnection occurs at ir e.g. l81
!     Failed tracing overhaul 7/2/21
!     output polarity information in sign of b1rs
