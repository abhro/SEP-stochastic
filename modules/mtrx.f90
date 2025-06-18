module mtrx
  use iso_fortran_env, only: real64
  use param, only: CSPEED, N_R, N_THETA, N_PHI, pi, twopi
  use epv, only: rp2beta
  use dxx, only: cofm
  use dmumu, only: cofdu
  use file_op, only: open_file_from_environment

  implicit none

contains

  ! delete this function when the fortran 2008 standards are supported
  real(kind=real64) function norm2(arr)
    real(kind=real64), intent(in) :: arr(:)
    norm2 = sqrt(sum(arr * arr))
  end function norm2

  function cartesian2spherical(x) result(r)
    real(kind=real64), intent(in) :: x(3)
    real(kind=real64), intent(out) :: r(3)

    r(1) = norm2(x)
    r(2) = acos(x(3) / r(1))  ! acos(z/r)
    r(3) = atan2(x(2), x(1))  ! atan(y/x)
  end function

  function spherical2cartesian(r) result(x)
    real(kind=real64), intent(in) :: r(3)
    real(kind=real64), intent(out) :: x(3)

    cartesian(1) = r(1) * sin(r(2)) * cos(r(3)) ! r sin(theta) cos(phi)
    cartesian(2) = r(1) * sin(r(2)) * sin(r(3)) ! r sin(theta) sin(phi)
    cartesian(3) = r(1) * cos(r(2))             ! r cos(theta)
  end function


  function mxptr(gm) result(xptr)
    !!  calculate matrix for xyz' ellipsoid coordinate Kwon to spheric
    real(kind=real64), intent(in) :: gm
    real(kind=real64)             :: xptr(3,3)
    real(kind=real64)             :: sing, cosg

    sing = sin(gm)
    cosg = cos(gm)

    xptr(1,:) = [0.d0,  0.d0, 1.d0]
    xptr(2,:) = [cosg, -sing, 0.d0]
    xptr(3,:) = [sing,  cosg, 0.d0]
  end function


  function dmxptr(gm, dgm)
    !!  calculate matrix for xyz' ellipsoid coordinate Kwon to spheric
    real(kind=real64), intent(in) :: gm, dgm
    real(kind=real64)             :: dmxptr(3,3)
    real(kind=real64)             :: sing, cosg

    sing = sin(gm)
    cosg = cos(gm)

    dmxptr(1,:) = [     0.d0,      0.d0, 0.d0]
    dmxptr(2,:) = [-sing*dgm, -cosg*dgm, 0.d0]
    dmxptr(3,:) = [ cosg*dgm, -sing*dgm, 0.d0]
  end function


  function mbtr(uax1, uax2, uax3) result(b2r)
    !!  calculate martix from magnetic to polar spheric coordinates
    real(kind=real64), intent(in)  :: uax1(3), uax2(3), uax3(3)
    real(kind=real64)              :: b2r(3,3)

    b2r(:,1) = uax1(:)
    b2r(:,2) = uax2(:)
    b2r(:,3) = uax3(:)
  end function


  function mrtx(sintheta, costheta, sinphi, cosphi)
    !!  calculate martix from polar spheric to xyz coordinates
    real(kind=real64), intent(in) :: sintheta, costheta, sinphi, cosphi
    real(kind=real64)             :: mrtx(3,3)

    mrtx(1,:) = [sintheta*cosphi, costheta*cosphi, -sinphi]
    mrtx(2,:) = [sintheta*sinphi, costheta*sinphi,  cosphi]
    mrtx(3,:) = [       costheta,       -sintheta,    0.d0]

  end function

  function dmrtx(sintheta, costheta, sinphi, cosphi, dtheta, dphi)
    !!  calculate martix from polar spheric to xyz coordinates
    real(kind=real64), intent(in) :: sintheta, costheta, dtheta
    real(kind=real64), intent(in) :: sinphi, cosphi, dphi
    real(kind=real64)             :: dmrtx(3,3)

    dmrtx(1,1) =  costheta*dtheta*cosphi - sintheta*sinphi*dphi
    dmrtx(1,2) = -sintheta*dtheta*cosphi - costheta*sinphi*dphi
    dmrtx(1,3) = -cosphi*dphi
    dmrtx(2,1) =  costheta*dtheta*sinphi + sintheta*cosphi*dphi
    dmrtx(2,2) = -sintheta*dtheta*sinphi + costheta*cosphi*dphi
    dmrtx(2,3) = -sinphi*dphi
    dmrtx(3,1) = -sintheta*dtheta
    dmrtx(3,2) = -costheta*dtheta
    dmrtx(3,3) = 0.0
  end function

  !> trilinear interpolation
  function trilinear(phic, x) result(phi)

    !>  the value of phi at the corner of cubic box of side 1
    real(kind=real64), intent(in)  :: phic(2,2,2)
    !>  location inside the cube (0<=x<=1) or outside x<0 x>1
    real(kind=real64), intent(in)  :: x(3)
    !>  interpolated phi
    real(kind=real64)              :: phi

    phi = phic(1,1,1) * (1-x(1)) * (1-x(2)) * (1-x(3)) &
        + phic(2,1,1) *    x(1)  * (1-x(2)) * (1-x(3)) &
        + phic(1,2,1) * (1-x(1)) *    x(2)  * (1-x(3)) &
        + phic(1,1,2) * (1-x(1)) * (1-x(2)) *    x(3)  &
        + phic(2,1,2) *    x(1)  * (1-x(2)) *    x(3)  &
        + phic(1,2,2) * (1-x(1)) *    x(2)  *    x(3)  &
        + phic(2,2,1) *    x(1)  *    x(2)  * (1-x(3)) &
        + phic(2,2,2) *    x(1)  *    x(2)  *    x(3)

  end function

  function trilineardif(phic, x) result(dphi)
    !>  phic the value of phi at the corner of cubic box of side 1
    real(kind=real64), intent(in)  :: phic(2,2,2)
    !>  location inside the cube (0<=x<=1) or outside x<0 x>1
    real(kind=real64), intent(in)  :: x(3)
    !> dphi(1)=dphi/dx, dphi(2)=dphi/dy, dphi(3)=dphi/dz
    real(kind=real64)              :: dphi(3)

    dphi(1) = (phic(2,1,1)-phic(1,1,1)) * (1-x(2)) * (1-x(3)) &
            + (phic(2,2,1)-phic(1,2,1)) *    x(2)  * (1-x(3)) &
            + (phic(2,1,2)-phic(1,1,2)) * (1-x(2)) *    x(3)  &
            + (phic(2,2,2)-phic(1,2,2)) *    x(2)  *    x(3)
    dphi(2) = (phic(1,2,1)-phic(1,1,1)) * (1-x(1)) * (1-x(3)) &
            + (phic(2,2,1)-phic(2,1,1)) *    x(1)  * (1-x(3)) &
            + (phic(1,2,2)-phic(1,1,2)) * (1-x(1)) *    x(3)  &
            + (phic(2,2,2)-phic(2,1,2)) *    x(1)  *    x(3)
    dphi(3) = (phic(1,1,2)-phic(1,1,1)) * (1-x(1)) * (1-x(2)) &
            + (phic(2,1,2)-phic(2,1,1)) *    x(1)  * (1-x(2)) &
            + (phic(1,2,2)-phic(1,2,1)) * (1-x(1)) *    x(2)  &
            + (phic(2,2,2)-phic(2,2,1)) *    x(1)  *    x(2)
  end function
end module mtrx
