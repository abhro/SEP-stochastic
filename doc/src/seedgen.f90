program seedgen

  ! USES environment variables
  !     SEEDS_FILENAME

  use datetime_utils, only: seconds_of_year
  use random, only: ran2
  use file_op, only: open_file_from_environment
  use param, only: NSEEDMAX

  implicit none

  integer             :: iseed ! store random state
  integer             :: i, j ! iteration variables
  integer             :: nran(NSEEDMAX) ! store all seeds
  integer             :: seeds_outfileunit

  integer             :: iran

  ! The random generator seed is -1 * Seconds of Year
  iseed = -1 * seconds_of_year()
  print *, "iseed = ", iseed


  call open_file_from_environment(&
    "SEEDS_FILE", seeds_outfileunit, 'FORMATTED')

  do i = 1, NSEEDMAX

    ! get a new random number and store it
1   iran = int(ran2(iseed) * 2147483647 - 1)
    !print *,  "iran = ", iran
    nran(i) = -iran

    ! make sure the random number hasn't already been generated
    do j = 1, i-1
      if (nran(i) == nran(j)) goto 1
    end do

    write(seeds_outfileunit,*)  nran(i)
  end do
  close(seeds_outfileunit)
end program
