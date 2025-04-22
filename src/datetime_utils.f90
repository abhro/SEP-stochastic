module datetime_utils
  use iso_fortran_env, only: real64
  implicit none

  type gregorian_date
    integer :: iday  ! day of the month
    integer :: month ! month of the year
    integer :: iyear ! year
    integer :: iyday ! year day
  end type

contains
  pure logical function is_leap_year(year)
    !! returns if a function is a Gregorian leap year
    integer, intent(in) :: year
    is_leap_year = &
      mod(year, 400) == 0 .or. &
      (mod(year, 4) == 0 .and. mod(year, 100) /= 0)
  end function

  function seconds_of_year() result(total_seconds)
    !! return how many seconds have passed since the year started

    integer :: total_seconds

    ! how many days have passed since the start of the year at
    ! the start of this month. assumes leap year
    integer, parameter :: days_in_month(12) = [&
      0,    31,  59,  90, 120, 151, &
      181, 212, 243, 273, 304, 334]

    integer :: datetime(8)
    integer :: cur_year, cur_dom, cur_doy
    integer :: cur_hour, cur_min, cur_sec

    ! get current date and time and extract values from the array
    call date_and_time(values=datetime)
    cur_hour = datetime(5)
    cur_min  = datetime(6)
    cur_sec  = datetime(7)
    cur_year = datetime(1)
    cur_dom  = datetime(3)
    cur_doy  = cur_dom + days_in_month(datetime(2))
    ! if this is not a leap year take off the leap day
    if (.not. is_leap_year(cur_year)) cur_doy = cur_doy - 1

    total_seconds = cur_doy*86400 + cur_hour*3600 + cur_min*60 + cur_sec
  end function


  subroutine caldate(jday, iyear, iyday)
    !!  This routine takes the modified Julian date and
    !!  converts it to a date and time string.
    !!
    !!  Calls:  GREGORIAN

    !--------------------------------------------------
    !  Define local data.
    !--------------------------------------------------
    !
    real(kind=real64), intent(in) :: jday   !!  modified Julian day (integer)

    ! older output arguments
    !character(len=11)             :: dchar !!  date string (character)
    !character(len=11)             :: tchar !!  time string (character)
    integer, intent(out)          :: iyday  !!  year day (integer)
    integer, intent(out)          :: iyear  !!  year (integer)
    integer                       :: iday   !!  day of the month (integer)
    integer                       :: month  !!  month of the year (integer)

    integer                       :: ihour, imin, isec, julian
    real(kind=real64)             :: fjulian, hour, min
    character(len=3), parameter   :: mchar(12) = [ &
      'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
      'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    !            original from HOPS model

    integer, parameter :: offset = 2440000

    type(gregorian_date) :: gdtoday

    !==================================================
    !  Begin executable code.
    !==================================================
    !
    !  Add offset to get true Julian date.

    julian = offset + int(jday)

    !fjulian = jday - int(jday) + 0.5
    fjulian = abs(jday - int(jday))

    if (fjulian > 1.0) then
      julian = julian + 1
      fjulian = fjulian - 1.0
    end if

    !  Compute Gregorian date.
    gdtoday = gregorian(julian)
    iday = gdtoday%iday
    month = gdtoday%month
    iyear = gdtoday%iyear
    iyday = gdtoday%iyday

    !--------------------------------------------------
    !  Form date and time strings.
    !--------------------------------------------------

    hour = fjulian * 24.0
    ihour = int(hour)
    min = (hour - float(ihour)) * 60.0
    imin = int(min)
    isec = (min-imin) * 60.0

    !write(dchar,"(a3,i3,i5)") mchar(month), iday, iyear
    !write(tchar,"(i2.2,':',i2.2,':',i2.2)") ihour, imin, isec
  end subroutine

  pure type(gregorian_date) function gregorian(julian)
    !!  This routine converts Julian day number to calendar (Gregorian) date.


    !----------------------------------------------
    !  Define local data.
    !----------------------------------------------

    integer, intent(in) :: julian !!     Julian day (integer)

    integer, parameter  :: IYD(13) = [&
      1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
    integer, parameter  :: IYDL(13) = [&
      1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336, 367]
    integer, parameter  :: IGREG = 2299161

    ! the following parameters are from wikipedia, which is from
    ! Richards 2013: Explanatory Supplement to the Astronomical Almanac, 3rd ed
    ! and Richards 1998: Mapping Time: The Calendar and its History
    integer, parameter  :: YGREG = 4716, JGREG = 1401, MGREG = 2
    integer, parameter  :: NGREG = 12, RGREG = 4, PGREG = 1461
    integer, parameter  :: VGREG = 3, UGREG = 5, SGREG = 153
    integer, parameter  :: WGREG = 2, BGREG = 274277, CGREG = -38

    integer             :: iday, month, iyear, iyday

    integer :: f, e, g, h

    !==============================================
    !  Begin executable code.
    !==============================================

    ! all divisions here are meant to be integer divisions
    f = julian + JGREG + (((4 * julian + BGREG) / 146097) * 3) / 4 + CGREG
    e = RGREG * f + VGREG
    g = mod(e, PGREG) / RGREG
    h = UGREG * g + WGREG

    iday = mod(h, SGREG) / UGREG + 1
    month = mod(h / SGREG + MGREG, NGREG) + 1
    iyear = (e / PGREG) - YGREG + (NGREG + MGREG - month) / NGREG

    !if (julian > IGREG) then
    !  jalpha = int((julian - 1867216 - 0.25) / 36524.25)
    !  ja = julian + 1 + jalpha - int(0.25 * jalpha)
    !else
    !  ja = julian
    !end if
    !jb = ja + 1524
    !jc = int(6680.0 + (jb - 2439870 - 122.1) / 365.25)
    !jd = 365*jc + int(0.25*jc)
    !je = int((jb-jd) / 30.6001)
    !iday = jb - jd - int(30.6001*je)
    !month = je - 1
    !if (month > 12) month = month - 12
    !iyear = jc - 4715
    !if (month > 2) iyear = iyear - 1
    !if (iyear < 0) iyear = iyear - 1
    if (is_leap_year(iyear)) then
      iyday = IYDL(month) + iday - 1
    else
      iyday = IYD(month) + iday - 1
    end if

    gregorian%iday = iday
    gregorian%month = month
    gregorian%iyear = iyear
    gregorian%iyday = iyday
  end function


  real(kind=real64) function modjulianday(year, month, day, fracday)
    !! calculate the julian day from day, month, year and fraction of a day

    real(kind=real64), intent(in) :: fracday
    integer, intent(in)           :: day, month, year
    integer, parameter            :: offset = 2440000
    !                    original from HOPS model

    modjulianday = julday(month, day, year) - offset + fracday
  end function


  integer function julday(mm, id, iyyy)
    integer, intent(in)     :: mm, id, iyyy
    integer, parameter      :: IGREG = 15 + 31 * (10 + 12*1582)
    integer                 :: iyy
    integer                 :: jy, jm, ja
    if (iyyy == 0) stop 'there is no year zero.'
    iyy = iyyy
    if (iyy <  0) iyy = iyy + 1
    if (mm > 2) then
      jy = iyy
      jm = mm + 1
    else
      jy = iyy - 1
      jm = mm + 13
    end if
    julday = int(365.25*jy) + int(30.6001*jm) + id + 1720995
    if (id + 31*(mm+12*iyy) >= IGREG) then
      ja = int(0.01*jy)
      julday = julday + 2 - ja + int(0.25*ja)
    end if
  end function
end module datetime_utils
