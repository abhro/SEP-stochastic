module param
  use iso_fortran_env, only: real64, int64

  implicit none

  ! physical constants

  real(kind=real64), parameter :: PI = 3.14159265358979323846d0, TWOPI = 2.d0*PI

  real(kind=real64), parameter :: &
    CSPEED = 25.8441774d0, &    !! speed of light in R_sun/min
    QoMSI = 9.57883376d7, &     !! proton charge-to-mass ratio in coulomb/kg
    EP = 0.938d0, &             !! proton rest energy in GeV
    EE = 0.510998918d-3, &      !! electron rest energy in GeV
    GAMMA_CS = 1.6666666666666666666d0 !! heat capacity ratio

  ! program parameters

  integer, parameter :: NSEEDMAX = 2001, NFMAX = 200

  integer(kind=int64), parameter :: NBASE = 40, NSPMAX = 20

  ! for PFSS
  real(kind=real64), parameter :: RSS = 2.5d0
  integer, parameter :: N_R = 150, N_THETA = 180, N_PHI = 360 !! grid size
  real(kind=real64)  :: bgrid(0:N_R, 0:N_THETA, 0:N_PHI, 3)     !! magnetic field at grid points
  !real(kind=real64) :: cvgrid(0:N_R, 0:N_THETA, 0:N_PHI, 3)    !! curvature?
  real(kind=real64)  :: gbgrid(0:N_R, 0:N_THETA, 0:N_PHI, 3)    !! grad B
  real(kind=real64)  :: b1rsgrid(0:N_R, 0:N_THETA, 0:N_PHI)     !! mapping to solar surface

  real(kind=real64), parameter :: epsilon(4) = [0.04, 0.01, 0.04, 0.003]

  ! conversion factors
  real(kind=real64), parameter :: RAD_TO_DEG = 180d0 / pi
  real(kind=real64), parameter :: DEG_TO_RAD = pi / 180d0
  real(kind=real64), parameter :: RS_PER_MIN_TO_KM_PER_SEC = 6.96340d5 / 60.0
end module param
