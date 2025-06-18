program maggrid
  ! USES environment variables
  !     SHTC_FILE
  !     MAGGRID_FILE

  use iso_fortran_env, only: real64, output_unit

  use fgsl, only: fgsl_sf_legendre_Plm

  use param, only: pi, N_R, N_THETA, N_PHI, DEG_TO_RAD !, bgrid, gbgrid, cvgrid
  use mtrx, only: norm2
  use file_op, only: read_shtc, write_maggrid

  implicit none

  include 'omp_lib.h'

  integer, parameter :: NS = 40
  real(kind=real64) :: rvec(3), b(3), cvtu(3), gbmag(3)
  real(kind=real64) :: g(0:NS,0:NS), h(0:NS,0:NS)
  real(kind=real64) :: srfctrtl(NS,NS)
  real(kind=real64) :: bmag, dbbds
  integer           :: id
  integer           :: chunk
  real(kind=real64) :: starttime, stoptime
  integer           :: i, j, k

  ! bgrid and gbgrid must be manually specified, using the one in param segfaults
  real(kind=real64)  :: bgrid(0:N_R,0:N_THETA,0:N_PHI,3)
  real(kind=real64)  :: gbgrid(0:N_R,0:N_THETA,0:N_PHI,3)

  !>Vector of lower-triangular matrices, where each matrix is of the form
  !>```math
  !>\begin{bmatrix}
  !>P_0^0    & \cdot    & \cdot    & \cdots & \cdot       & \cdot       & \cdot  \\
  !>P_1^0    & P_1^1    & \cdot    & \cdots & \cdot       & \cdot       & \cdot  \\
  !>P_2^0    & P_2^1    & P_2^2    & \cdots & \cdot       & \cdot       & \cdot  \\
  !>\vdots   & \vdots   & \vdots   & \ddots & \vdots      & \vdots      & \vdots \\
  !>P_{38}^0 & P_{38}^1 & P_{38}^2 & \cdots & P_{38}^{38} & \cdot       & \cdot  \\
  !>P_{39}^0 & P_{39}^1 & P_{39}^2 & \cdots & P_{39}^{38} & P_{39}^{39} & \cdot  \\
  !>P_{40}^0 & P_{40}^1 & P_{40}^2 & \cdots & P_{40}^{38} & P_{40}^{39} & P_{40}^{40}
  !>\end{bmatrix}
  !>```
  !>where each ``P_ℓ^m`` is a matrix evaluated at an ``x`` given by the vector's
  !>index, and each ``⋅`` signifies a zero entry in the matrix.
  !> DALP_TABLE is the first derivative
  !> D2ALP_TABLE is the second derivative
  real(kind=real64), dimension(0:N_THETA, 0:NS, 0:NS) :: ALP_TABLE, DALP_TABLE, D2ALP_TABLE

  chunk = 1

  !  load spherical harmonic transform coefficients
  ! USES environment variable SHTC_FILE
  call read_shtc(g, h, NS)

  call factorialratio(srfctrtl)

  starttime = omp_get_wtime()

  !$OMP PARALLEL DEFAULT(private) SHARED(chunk,bgrid,gbgrid,g,h,srfctrtl)
  id = omp_get_thread_num()
  print *, __FILE__, ":", __LINE__, ": I am thread", id, ", Thread start time", starttime
  flush(output_unit)
  !$OMP DO SCHEDULE(STATIC,chunk)
  do i = 0, N_R
    rvec(1) = 1.0 + i/100.0
    do j = 0, N_THETA
      rvec(2) = j * pi/180.0
      do k = 0, N_PHI
        rvec(3) = k * pi/180.0
        call magfield(rvec, b, bmag, cvtu, gbmag, dbbds, srfctrtl, g, h)
        bgrid(i,j,k,:) = b(:)
        !  curvature not needed
        !cvgrid(i,j,k,:) = cvtu(:)
        gbgrid(i,j,k,:) = gbmag(:)
      end do
    end do
  end do
  !$OMP END PARALLEL

  call write_maggrid(bgrid, gbgrid)

  stoptime = omp_get_wtime()
  print "('elapsed time: ',f0.5,'s, my id: ',i0)", stoptime - starttime, id

contains

  subroutine factorialratio(fac_table)
    real(kind=real64), intent(out) :: fac_table(NS,NS)

    integer :: l, m
    ! we don't need the factor for l = 0 (which only has m = 0) or any m = 0
    do l = 1, NS
      do m = 1, l
        ! say, S = srfctrtl
        ! then, S_{l,m} = sqrt(2 * (l-m)! / (l+m)!)
        ! but using gamma directly would cause an overflow
        ! LaTeX: S_{\ell,m} = \sqrt{2 \frac{(\ell - m)!}{(\ell + m)!}}
        fac_table(l,m) = sqrt(2.d0 * exp( &
            ! the actual ratio part (turns into subtraction for ln)
            log_gamma(real(l-m+1, kind=real64)) - log_gamma(real(l+m+1, kind=real64)) &
          ))
      end do
    end do
  end subroutine

  subroutine magfield(rvec, b, bmag, cvtu, gbmag, dbbds, srfctrtl, g, h)
    real(kind=real64), intent(in)  :: rvec(3) ! position vector
    real(kind=real64), intent(out) :: b(3), bmag, cvtu(3), gbmag(3), dbbds
    real(kind=real64), intent(in)  :: srfctrtl(NS,NS)
    real(kind=real64), intent(in)  :: g(0:NS,0:NS), h(0:NS,0:NS)

    real(kind=real64), parameter   :: RSS = 2.5d0

    ! spherical coordinates of position vector (radial, colatitude/inclination, azimuthal)
    real(kind=real64)              :: r, theta, phi
    real(kind=real64)              :: ab(3), db(3,3) !derivative of components only

    ! legendre polynomials
    real(kind=real64)              :: alp(0:NS+2,0:NS+2)
    real(kind=real64)              :: schmt, dalp, d2alp

    real(kind=real64)              :: costheta, sintheta
    real(kind=real64)              :: cosmphi, sinmphi
    real(kind=real64)              :: r0rl2, rrs2l1, r0rs2l1, ghmf, dghmf
    integer                        :: l, m

    r = rvec(1)
    theta = rvec(2)
    phi = rvec(3)
    if (r > RSS) r = RSS
    costheta = cos(theta)
    sintheta = sin(theta)
    if (sintheta <= 1.745d-4) then
      sintheta = 1.745d-4
      costheta = sqrt(1.0 - sintheta**2)
      if (theta > 1.57) costheta = -costheta
    end if
    !  Associated Legendre Polynomial
    alp = 0.0
    do l = 0, NS + 2
      do m = 0, l
        alp(l,m) = plgndr(l, m, costheta)
      end do
    end do

    b = 0.0
    db = 0
    do l = 1, NS
      r0rl2 = r**(-2-l)
      rrs2l1 = (r/RSS)**(2*l+1)
      r0rs2l1 = RSS**(-2*l-1)*l + 1 + l
      do m = 0, l
        !  convert to Schmidt (semi-)normalized Associated Legendre Polynomial
        if (m == 0) then
          schmt = 1
          dalp = alp(l,1)
          d2alp = (-(l+1)*l*alp(l,0) + alp(l,2)) / 2.0
        else
          schmt = (-1)**m * srfctrtl(l,m)
          dalp = (-(l+m)*(l-m+1)*alp(l,m-1) + alp(l,m+1)) / 2.0
          if (m == 1) then
            d2alp = (alp(l,3) - (3*l*l-3*l-2) * alp(l,1)) / 4.0
          else
            d2alp = (&
              alp(l,m+2) &
              - ((l+m+1)*(l-m)+(l+m)*(l-m+1)) * alp(l,m) &
              + (l+m)*(l-m+1)*(l+m-1)*(l-m+2) * alp(l,m-2) &
              )/4.0
          end if
        end if
        cosmphi = cos(m*phi)
        sinmphi = sin(m*phi)
        ghmf = g(l,m)*cosmphi + h(l,m)*sinmphi
        dghmf = m * (-g(l,m)*sinmphi + h(l,m)*cosmphi)
        ab(1) =  r0rl2 * ghmf  * schmt * alp(l,m) * (l+1+l*rrs2l1) / r0rs2l1
        ab(2) = -r0rl2 * ghmf  * schmt * dalp * (1-rrs2l1) / r0rs2l1
        ab(3) = -r0rl2 * dghmf * schmt * alp(l,m) / sintheta * (1-rrs2l1) / r0rs2l1

        b(:) = ab(:) + b(:)

        if (ab(1) /= 0) then
          db(1,1) = db(1,1) - ab(1) * ((l+2)*(l+1)-l*(l-1)*rrs2l1) / (l+1+l*rrs2l1) / r
          db(2,1) = db(2,1) + ab(1) / alp(l,m) * dalp
          db(3,1) = db(3,1) + ab(1) / ghmf * dghmf
        end if
        if (ab(2) /= 0) then
          db(1,2) = db(1,2) - ab(2) * (l+2+(l-1)*rrs2l1) / (1-rrs2l1) / r
          db(2,2) = db(2,2) + ab(2) / dalp * d2alp
          db(3,2) = db(3,2) + ab(2) / ghmf * dghmf
        end if
        if (ab(3) /= 0) then
          db(1,3) = db(1,3) - ab(3) * (l+2+(l-1)*rrs2l1) / (1-rrs2l1) / r
          db(2,3) = db(1,3) + ab(3) / alp(l,m) * dalp - ab(3) * costheta / sintheta
          db(3,3) = db(1,3) + ab(3) / dghmf * m * m * ghmf
        end if
      end do
    end do
    bmag = norm2(b)

    gbmag(1) = sum(b(1:3) * db(1,1:3)) / bmag
    gbmag(2) = sum(b(1:3) * db(2,1:3)) / (bmag * r)
    gbmag(3) = sum(b(1:3) * db(3,1:3)) / (bmag * r * sintheta)
    dbbds = sum(b(1:3) * gbmag(1:3)) / (bmag * bmag)

    cvtu(1) = ( &
        b(1)*db(1,1) + b(2)*db(2,1)/r &
        + b(3)*db(3,1)/r/sintheta - (b(2)*b(2)+b(3)*b(3))/r &
      ) / bmag
    cvtu(2) = ( &
        b(1)*db(1,2) + b(2)*db(2,2)/r &
        + b(3)*db(3,2)/r/sintheta &
        + (b(2)*b(1)-b(3)*b(3)*costheta/sintheta)/r &
      ) / bmag
    cvtu(3) = ( &
        b(1)*db(1,3) + b(2)*db(2,3)/r &
        + b(3)*db(3,3)/r/sintheta &
        + (b(3)*b(1)+b(3)*b(2)*costheta/sintheta)/r &
      ) / bmag

    cvtu = (cvtu - dbbds*b) / bmag  ! curvature perpendicular to B

    if (rvec(1) > RSS) then
      b(1) = b(1) * (RSS/rvec(1))**2
      bmag = norm2(b)
      cvtu = 0.0
      gbmag(1) = -2 * bmag / rvec(1)
      gbmag(2) = gbmag(2) * (RSS/rvec(1))**3
      gbmag(3) = gbmag(3) * (RSS/rvec(1))**3
      dbbds = b(1) * gbmag(1) / bmag / bmag
    end if

  end subroutine

  real(kind=real64) function plgndr(l, m, x)
    integer, intent(in)           :: l, m
    real(kind=real64), intent(in) :: x

    integer           :: i     ! iteration variable
    real(kind=real64) :: pll
    real(kind=real64) :: somx2 ! square root of 1 - x^2
    real(kind=real64) :: pmm   ! legendre polynomial of order (m, m)
    real(kind=real64) :: pmmp1 ! legendre polynomial of order (m+1, m)

    if (m < 0 .or. m > l .or. abs(x) > 1.0) then
      print *, 'bad arguments in plgndr'
      stop 1
    end if

    pmm = 1.0 ! compute P_m^m = (-1)^m (1-x^2)^{m/2} \prod_{i=1}^m (2i+1)
    if (m > 0) then
      somx2 = sqrt(1.0 - x*x)
      ! compute \prod_{i=1}^m (2i+1) = (2m-1)!!
      do i = 1, m
        pmm = pmm * (2*i + 1)
      end do
      pmm = pmm * somx2**m
      if (mod(m, 2) == 1) pmm = -pmm ! factor of (-1)^m
    end if

    if (l == m) then
      plgndr = pmm
    else
      pmmp1 = x * (2*m + 1) * pmm ! compute P_{m+1}^m = x (2m + 1) P_m^m
      if (l == m+1) then
        plgndr = pmmp1
      else
        do i = m+2, l
          pll = (x * (2*i-1) * pmmp1 - (i+m-1) * pmm) / (i - m)
          pmm = pmmp1
          pmmp1 = pll
        end do
        plgndr = pll
      end if
    end if
  end function


  subroutine init_aplm(alp_arr, dalp_arr, d2alp_arr)
    !> populate alp, dalp, d2alp with the associated legendre polynomials,
    !> and its first and second derivatives respectively
    real(kind=real64), intent(out) ::   alp_arr(0:N_THETA,0:NS,0:NS)
    real(kind=real64), intent(out) ::  dalp_arr(0:N_THETA,0:NS,0:NS)
    real(kind=real64), intent(out) :: d2alp_arr(0:N_THETA,0:NS,0:NS)
    integer :: j
    integer :: l, m
    real(kind=real64) :: theta, costheta
    !alp(0,:,:) = 1
    do j = 0, N_THETA
      costheta = cos(j * DEG_TO_RAD)
      do l = 0, NS
        do m = 0, NS
          alp_arr(j, l, m) = fgsl_sf_legendre_Plm(l, m, costheta)
        end do
      end do
    end do
  end subroutine

end program
