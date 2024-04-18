program average
  use ieee_arithmetic, only: ieee_is_nan
  use iso_fortran_env, only: real64

  implicit none

  real(kind=real64) :: t, fl, fl2
  real(kind=real64) :: ar(9)
  integer           :: na, i
  na = 0

  do i = 1,6
    read(*,*)
  enddo

  do while (t > -10)
    read(*,*) t, ar
    !if (ar(7) * ar(7) >= 0) then ! nan testing
    if (.not. ieee_is_nan(ar(7))) then
      na = na + 1
      fl = ar(7) + fl
      fl2 = ar(7)**2 + fl2
    end if
  end do
  print *, fl/na, sqrt(fl2)/na
  stop
end program
