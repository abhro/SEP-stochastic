program arriv_time
  use iso_fortran_env, only: real64
  use epv, only: e2v
  implicit none
  real(kind=real64) :: rnz, rnm
  common/specie/rnz,rnm

  real(kind=real64) :: e, v

  rnz = 2
  rnm = 4
  print *, 'e = ?(GeV)'
  read (*,*) e
  v = e2v(e)
  print *, 'v = ', v, 'AU/day'
  print *, 'arriving time=', 1.0/v, 'days'
end program arriv_time
