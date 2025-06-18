MOdule rksolvers
  use iso_fortran_env, only: real64
  implicit none

    interface
      function odefun(x, y, p) result(v)
        use iso_fortran_env, only: real64
        implicit none
        real(kind=real64) :: v(3)
        real(kind=real64), intent(in) :: x
        real(kind=real64), intent(in) :: y(3)
        real(kind=real64), intent(in) :: p
      end function
    end interface

contains

  !function rk4(y, dydx, n, x, h, f, param) result(yout)
  !  integer, intent(in)            :: n
  !  real(kind=real64), intent(in)  :: h, x, dydx(n), y(n)
  !  procedure(odefun)              :: f
  !  real(kind=real64)              :: yout(n)
  !  !real(kind=real64), dimension(3), external :: odefun
  !  integer                        :: i
  !  real(kind=real64)              :: hh, xh, dym(n), dyt(n), yt(n)
  !  real(kind=real64)              :: param

  !  hh = h * 0.5
  !  xh = x + hh
  !  yt(:) = y(:) + hh * dydx(:)
  !  dyt = f(xh, yt, param)
  !  yt(:) = y(:) + hh * dyt(:)
  !  dym = f(xh, yt, param)
  !  yt(:) = y(:) + h*dym(:)
  !  dym(:) = dyt(1:n) + dym(:)
  !  dyt = f(x+h, yt, param)
  !  yout(1:n) = y(:) + (h/6.0) * (dydx(:) + dyt(:) + 2.0*dym(:))
  !end function

  !> single step of RK4
  function rk4(f, x0, y0, h, odefun_param) result(yout)
    procedure(odefun) :: f
    real(kind=real64) :: x0
    real(kind=real64) :: y0(:)
    real(kind=real64) :: h
    real(kind=real64) :: odefun_param
    real(kind=real64), allocatable :: yout(:)

    integer :: n

    real(kind=real64), allocatable :: k1(:), k2(:), k3(:), k4(:)

    n = size(y0)

    allocate(k1(n), k2(n), k3(n), k4(n), yout(n))

    ! basically just lifted from wikipedia
    k1 = f(x0,       y0,            odefun_param)
    k2 = f(x0 + h/2, y0 + h/2 * k1, odefun_param)
    k3 = f(x0 + h/2, y0 + h/2 * k2, odefun_param)
    k4 = f(x0 + h,   y0 + h   * k3, odefun_param)

    yout = y0 + h/6 * (k1 + 2*k2 + 2*k3 + k4)
  end function
end module
