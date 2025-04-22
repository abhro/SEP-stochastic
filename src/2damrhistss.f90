program histogram
  !!   build histogram

  use iso_fortran_env, only: real64

  implicit none

  real(kind=real64), parameter :: twopi = 2.d0 * 3.14159265358979323846d0

  ! km
  real(kind=real64), parameter :: rss = 2.5
  real(kind=real64), parameter :: rs = 6.96e5 !! solar radius in km

  ! raw km/s converted to Rs/min
  !> solar wind speed in Rs/min
  real(kind=real64), parameter :: vsw = 400.0 / rs * 60

  !> angular frequency of sun's rotation
  real(kind=real64), parameter :: omega = twopi / (27.27 * 24 * 60)
  real(kind=real64), parameter :: k4ok2 = 12.4242, k6ok2 = 242.4242

  real(kind=real64)            :: al(8)
  real(kind=real64)            :: ra, theta, phi
  real(kind=real64)            :: hist(91,46), hs(91,46)
  real(kind=real64)            :: tpsw
  real(kind=real64)            :: xmu
  !integer                     :: aa(2)
  !integer                     :: specfile
  integer                      :: n, nn, ns, nz
  integer                      :: i, j, ii, jj, jjj, iy, ix


  !open(newunit=specfile, file='spec.dat_e00')
  !do i = 1,10
  !  read(newunit=specfile,*)
  !end do
  !close(newunit=specfile)

  n = 0
  hist = 0.0
  do i = 1, 1569000
    read(*, *, end=99) al, nn
    if (al(1) == -1000.00) exit
    ra = al(3)
    tpsw = &
      (ra - rss) &
      - rss * log(ra/rss) &
      + k4ok2 * (1/rss/2 - 1/ra + rss/ra**2/2) &
      + k6ok2 * (1/rss**3/12 - 1/ra**3/3 + rss/ra**4/4)
    tpsw = tpsw / vsw
    theta = al(4)
    phi = al(5) + omega*tpsw
    phi = atan2(sin(phi), cos(phi))
    write(12,*) theta, phi
    if (phi < 0) phi = phi + twopi
    n = n + 1
    if (ra >= 2.50) then
      xmu = cos(theta)
      iy = floor(theta/twopi*90) + 1
      ix = floor(phi/twopi*90) + 1
      hist(ix,iy) = hist(ix,iy) + 1
    end if
  end do

  99 continue
  !  smoothing
  hs = 0
  do i = 1, 45
    do j = 1, 90
      if (hist(j,i) < 10) then ! use 5x5 box smoothing
        ns = 0
        nz = 0
        do ii = i-2, i+2
          do jj = j-2, j+2
            if (ii > 0 .and. ii < 46) then
              jjj = jj
              if (jj < 1) jjj = 90 + jj
              if (jj > 90) jjj = jj - 90
              ns = ns+1
              if (hist(jjj,ii) /= 0) nz = nz + 1
              hs(j,i) = hs(j,i) + hist(jjj,ii)
            end if
          end do
        end do
      end if
      if(3*nz>ns) then
        hs(j,i) = hs(j,i) / ns
      else
        hs(j,i) = hist(j,i)
      end if
      if (hist(j,i) > 10 .and. hist(j,i) < 100) then ! use 3x3 box
        ns = 0
        do ii = i-1, i+1
          do jj = j-1, j+1
            if (ii > 0 .and. ii < 46) then
              jjj = jj
              if (jj < 1) jjj = 90 + jj
              if (jj > 90) jjj = jj - 90
              ns = ns + 1
              if (hist(jjj,ii) /= 0) nz = nz + 1
              hs(j,i) = hs(j,i) + hist(jjj,ii)
            end if
          end do
        end do
      end if
      if (3*nz > ns) then
        hs(j,i) = hs(j,i) / ns
      else
        hs(j,i) = hist(j,i)
      end if
      if (hist(j,i) >= 100) hs(j,i) = hist(j,i)
    end do
  end do

  hist = hs

  do i = 1, 45
    do j = 1, 90
      theta = (i-0.5) / 90 * twopi
      phi = (j-0.5) / 90 * 360
      write(*,"(2f8.3,e12.4,f6.0)") theta*360/twopi, phi, hist(j,i)/n/sin(theta), hist(j,i)
    end do
  end do
end program
