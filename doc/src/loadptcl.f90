module loadptcl
  use iso_fortran_env, only: real64
  !use param,           only: CSPEED!, NFMAX
  !use epv,             only: rp2e, rp2beta
  !use dxx,             only: get_rlambda0, cofm
  !use mtrx,            only: drvbmag
  use file_op,         only: open_file_from_environment

  implicit none
contains
  subroutine prepareptcl(tf, rf, ef, rmuf, np, nf)
    ! USES environment variables
    !     PTCL_FILE
    !real(kind=real64) :: tf(NFMAX), rf(3,NFMAX), ef(NFMAX), rmuf(NFMAX)
    real(kind=real64), intent(out), allocatable :: tf(:), rf(:,:), ef(:), rmuf(:)
    integer           :: np, nf
    !common/ldptcl/tf,rf,ef,rmuf,np,nf

    real(kind=real64) :: t0org, te, tdl, dmapjul, tcme0
    common /tmprm/t0org,te,tdl,dmapjul,tcme0

    character(len=256) :: dir
    common /dir/dir

    real(kind=real64) :: rnz, rnm
    common /specie/rnz,rnm

    real(kind=real64) :: sp0, gp, ap
    real(kind=real64) :: trgtfs(4)
    real(kind=real64) :: sp, scanw, h0
    common /srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0

    integer           :: i
    integer           :: ptcl_fileunit

    call open_file_from_environment("PTCL_FILE", ptcl_fileunit, "FORMATTED")
    !open(newunit=ptcl_fileunit, file=trim(dir)//'loadptcl.dat')

    read(ptcl_fileunit,*) ! sp0, gp, ap header line
    read(ptcl_fileunit,*) sp0, gp, ap
    read(ptcl_fileunit,*) ! header line for tdl
    read(ptcl_fileunit,*) tdl
    read(ptcl_fileunit,*) ! another header line (rnz, rnm)
    read(ptcl_fileunit,*) rnz, rnm
    read(ptcl_fileunit,*) ! you guessed it, header (nf, np)
    read(ptcl_fileunit,*) nf, np

    allocate(tf(nf), rf(3,nf), ef(nf), rmuf(nf))

    read(ptcl_fileunit,*) ! headers yay (time, r(3), ...)

    do i = 1, nf
      read(ptcl_fileunit,*) tf(i), rf(:,i), ef(i), rmuf(i)
    end do
    close(ptcl_fileunit)
  end subroutine

  !subroutine setasgp(t0,rp0)
  !  implicit none
  !  real(kind=real64) :: rp0(5), t0
  !  real(kind=real64) :: r(3), r2(3)
  !  real(kind=real64) :: bv(3), g2(3), dg2(3,3)
  !  real(kind=real64) :: cvtu(3), gbmag(3), bxgb2(3), pol
  !  real(kind=real64) :: sp0, gp, ap
  !  real(kind=real64) :: trgtfs(4)
  !  real(kind=real64) :: sp, scanw, h0
  !  common /srcmod/sp,sp0,gp,ap,trgtfs,scanw,h0

  !  real(kind=real64) :: p2, beta2
  !  real(kind=real64) :: tanp, sqr1p, dl
  !  real(kind=real64) :: bm, dbbds, b1s, gb1s(3)
  !  real(kind=real64) :: et1, rlambda1, rd, pa2, bm2
  !  ! set up parameters in f0mod to adjust artificial drift
  !  !  gp: order of pitch-angle anisotropy (typically arround 1.0)
  !  !  sp: inverse of pitch-angle anisotropy amplitude (>1.0)
  !  !  ap: lateral (lat-long) angular width of distribution

  !  gp = 1.0
  !  !  use radial diffusion to estimate pitch-angle anisotropy amplitude

  !  r(1:3) = rp0(1:3)
  !  p2 = rp0(4)
  !  beta2 = rp2beta(p2)
  !  call drvbmag(r, bv, bm, cvtu, gbmag, bxgb2, dbbds, pol, b1s, gb1s)
  !  tanp = abs(bv(3)/bv(1))
  !  sqr1p = sqrt(1+tanp*tanp)
  !  dl = r(1) / 2 * (sqr1p + log(tanp+sqr1p) / tanp)
  !  sp = 9*dl / (t0*beta2*CSPEED)
  !  if (sp > 0.9) sp = 0.9
  !  sp = 1/sp
  !  sp = 4.0

  !  !  estimate the radial diffusion
  !  et1 = rp2e(p2)
  !  rlambda1 = get_rlambda0(et1)
  !  rd = sqrt(1.5 * t0 * rlambda1 * beta2 * CSPEED/3.)

  !  !  perpendicular diffusion at half rd
  !  r2(1) = rd
  !  r2(2) = r(2)
  !  r2(3) = r(3)
  !  pa2 = rp0(5)
  !  call drvbmag(r, bv, bm2, cvtu, gbmag, bxgb2, dbbds, pol, b1s, gb1s)
  !  call cofm(r2, p2, pa2, beta2, bv, bm2, cvtu, gbmag, dbbds, b1s, gb1s, g2, dg2)
  !  !  estimate angular radius of perpendicular diffusion
  !  ap = sqrt(1.5*t0*g2(2)/(rd)**2)
  !  tanp = abs(bv(3)/bv(1))
  !  sqr1p = sqrt(1+tanp*tanp)
  !  ap = ap*sqr1p
  !  ap = 1e20
  !end subroutine
end module
