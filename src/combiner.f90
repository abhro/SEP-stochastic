program combon

  ! USES environment variables
  !     COMBINER_B1RS_INFILENAME
  !     COMBINER_MAGGRID_INFILENAME
  !     COMBINER_MAGGRID_OUTFILENAME

  use iso_fortran_env, only: real64
  use file_op, only: open_file_from_environment

  implicit none

  real(kind=real64) :: b(3), cv(3), gb(3), b1rs
  integer :: i, j, k
  integer :: b1rs_infileunit
  integer :: maggrid_infileunit, maggrid_outfileunit

  call open_file_from_environment(&
    "COMBINER_B1RS_INFILENAME", b1rs_infileunit, 'FORMATTED')
  call open_file_from_environment(&
    "COMBINER_MAGGRID_INFILENAME",  maggrid_infileunit, 'FORMATTED')
  call open_file_from_environment(&
    "COMBINER_MAGGRID_OUTFILENAME", maggrid_outfileunit, 'FORMATTED')

  do k = 0, 150
    do i = 0, 180
      do j = 0, 360
        read(b1rs_infileunit,*) b1rs
        read(maggrid_infileunit,*) b(1:3), gb(1:3)
        write(maggrid_outfileunit,*) b(1:3), gb(1:3), b1rs
        !read(f2) b(1:3), cv(1:3), gb(1:3)
        !write(f3) b(1:3), cv(1:3), gb(1:3), b1rs(1:2)
      end do
    end do
  end do
  close(b1rs_infileunit)
  close(maggrid_infileunit)
  close(maggrid_outfileunit)
end program combon
