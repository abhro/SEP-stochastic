module random
  use iso_fortran_env, only: real64
  implicit none
contains
  !real(kind=real64) function bindev(idum)
  !  integer, intent(in out) :: idum
  !  !U    USES ran2
  !  if (ran2(idum) > 0.5d0) then
  !    bindev = 1.0
  !  else
  !    bindev = -1.0
  !  end if
  !end function

  real(kind=real64) impure function gasdev(idum)
    !! This function uses Box-Muller transform to generate two normally
    !! distributed random numbers at a time
    integer, intent(in out) :: idum
    real(kind=real64)       :: fac, rsq, v1, v2
    !U    USES ran2
    logical, save           :: gaus_stored = .false.
    real(kind=real64), save :: gset = 0
    !$omp threadprivate(gset,gaus_stored)
    if (gaus_stored) then
      gasdev = gset
      gaus_stored = .false.
    else
      ! choose a point in the [-1, 1] x [-1, 1] square in the cartesian plane
1     v1 = 2.0 * ran2(idum) - 1.0
      v2 = 2.0 * ran2(idum) - 1.0

      ! if the point is also not within the unit circle, reject it and
      ! choose another point that will hopefully be better
      rsq = v1**2 + v2**2
      if (rsq >= 1.0 .or. rsq == 0.0) goto 1

      ! the actual Box-Muller part
      fac = sqrt(-2.0 * log(rsq)/rsq)

      gasdev = v1 * fac

      ! save this for the next time the function is called
      gset = v2 * fac
      gaus_stored = .true.
    end if
  end function


  real(kind=real64) function ran2(idum)
    ! this code is from WH Press's Numerical Recipes book

    integer, intent(in out) :: idum ! state of random number generator

    integer, parameter      :: &
      IM1 = 2147483563, IM2 = 2147483399, IMM1 = IM1-1, &
      IA1 = 40014, IA2 = 40692, &
      IQ1 = 53668, IQ2 = 52774, &
      IR1 = 12211, IR2 = 3791, &
      NTAB = 32, NDIV = 1 + IMM1/NTAB

    real(kind=real64), parameter :: AM = 1d0 / IM1, EPS = 1.2d-7
    real(kind=real64), parameter :: RNMX = 1.d0 - EPS

    integer, save           :: idum2 = 123456789, iv(NTAB) = 0, iy = 0

    integer                 :: j, k

    !$omp threadprivate(idum2,iv,iy)
    if (idum <= 0) then     ! initialize
      idum = max(-idum, 1)  ! prevent idum = 0
      idum2 = idum
      do j = NTAB+8,1,-1    ! load shuffle table after 8 warm-ups
        k = idum / IQ1
        idum = IA1*(idum-k*IQ1) - k*IR1
        if (idum < 0) idum = idum+IM1
        if (j <= NTAB) iv(j) = idum
      end do
      iy = iv(1)
    end if

    k = idum / IQ1
    idum = IA1 * (idum - k*IQ1) - k * IR1
    if (idum < 0) idum = idum + IM1

    k = idum2 / IQ2
    idum2 = IA2 * (idum2 - k*IQ2) - k * IR2
    if (idum2 < 0) idum2 = idum2 + IM2

    j = 1 + iy/NDIV
    iy = iv(j) - idum2
    iv(j) = idum
    if (iy < 1) iy = iy + IMM1
    ran2 = min(AM * iy, RNMX)
  end function
end module random
