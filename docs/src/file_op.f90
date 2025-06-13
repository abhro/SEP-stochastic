module file_op
  use iso_fortran_env, only: real64, int64, error_unit
  use param, only: TWOPI, NFMAX, NSEEDMAX, N_R, N_THETA, N_PHI, pi
  implicit none


  ! USES environment variables
  !     SEEDS_FILENAME
  !     PARAM_OUTDIR_PATH
  !     PARAM_NODES


contains
  ! opens file pclfil and assigns to nsts and
  ! opens file anfil  and assigns to nfl
  ! pclfil and anfil is set in readparam
  subroutine fl_open(nfl, nsts)
    integer, intent(out)   :: nfl, nsts

    character(len=256)    :: pclfil, finfil, sucfil, nodesfil, anfil, anstfil
    common/filnm/pclfil,finfil,sucfil,nodesfil,anfil,anstfil

    character(len=2)      :: pclfil_n

    integer               :: nfinish
    common/nfinish/nfinish

    character(len=2)      :: rankstr
    common/rankstr/rankstr
    nfinish = 0
    open(newunit=nsts, file=trim(pclfil))
    open(newunit=nfl, file=trim(anfil))
  end subroutine


  subroutine record_nodes(nodes)
    integer, intent(in) :: nodes
    integer             :: nnds
    character(len=256)  :: pclfil, finfil, sucfil, nodesfil, anfil, anstfil
    common/filnm/pclfil,finfil,sucfil,nodesfil,anfil,anstfil
    open(newunit=nnds, file=trim(nodesfil))
    write(nnds,*) nodes
    close(nnds)
  end subroutine


  !subroutine fl_close(nfl, nsts)
  !  integer, intent(in) :: nfl, nsts
  !  close(nfl)
  !  close(nsts)
  !end subroutine


  subroutine write_head(nfl, nf)
    integer, intent(in)           :: nfl
    integer, intent(in)           :: nf

    write(nfl,'(a,i0)') 'nf = ', nf
    write(nfl,'(a)') 'Fluxes are calculated at the locations:'
    write(nfl,'(a)') 'Time,postion(r theta phi),energy/n,\mu,flux,dflux'

  end subroutine


  !subroutine update_finish
  !  integer, parameter  :: NSTEP = 10

  !  character(len=256)  :: pclfil, finfil, sucfil, nodesfil, anfil, anstfil
  !  common/filnm/pclfil,finfil,sucfil,nodesfil,anfil,anstfil
  !  integer             :: np, nmu, nloc, nfinish
  !  common/ldptclnum/np,nmu,nloc
  !  common/nfinish/nfinish
  !  character(len=2)    :: rankstr
  !  common/rankstr/rankstr
  !  integer             :: nfin

  !  nfinish = nfinish + 1

  !  open(newunit=nfin, file=trim(finfil)//rankstr(1:2))
  !  write(nfin,*) nfinish, 'of', nmu*nloc
  !  close(nfin)
  !end subroutine

  subroutine open_file_from_environment(&
      env_var_name, fileunit, file_form, file_status)
    use iso_fortran_env, only: error_unit
    ! get a file name from the environment, open it, and return the open unit,
    ! taking care of errors along the way (stopping the program)
    character(len=*), intent(in)  :: env_var_name, file_form
    character(len=*), intent(in), optional :: file_status
    integer, intent(out)          :: fileunit
    character(len=256)            :: filename
    integer                       :: stat

    call get_environment_variable(env_var_name, filename, status=stat)
    if (stat /= 0) then
      write (error_unit, '(a)') &
        "Error reading $" // env_var_name // " from the environment"
      stop 1
    end if

    if (present(file_status)) then
      open(newunit=fileunit, file=trim(filename), &
           iostat=stat, form=file_form, status=file_status)
    else
      open(newunit=fileunit, file=trim(filename), &
           iostat=stat, form=file_form)
    end if
    if (stat /= 0) then
      write (error_unit, '(a)') &
        "Error opening " // trim(filename) // &
        " (filename obtained from $" // env_var_name // ")"
      stop 1
    end if
  end subroutine

  ! THIS IS THE MOST EVIL FUNCTION YOU HAVE EVER SEEN. COMMON BLOCKS GALORE.
  ! YOU'LL NEED TO REFACTOR IT, BUT IT WILL NOT BE AN EASY TASK. HAVE FUN.
  ! AND GOOD LUCK. YOU'LL NEED IT.
  subroutine read_param(nodes, nseeds)

    ! USES environment variables
    !     PARAM_OUTDIR_PATH
    !     PARAM_NODES


    !integer, parameter :: timenum = 10
    ! read in parameters to characterize particle species
    ! and solar wind and magnetic field condition

    integer, intent(out) :: nodes

    !common /magfield/b1au
    !common /heliosphere/r0,rau,vsw,omega

    real(kind=real64)   :: densw0, vsw, k4ok2, k6ok2
    real(kind=real64)   :: omega, b1au, vom, facip
    common /bpark/densw0,vsw,k4ok2,k6ok2,omega,b1au,vom,facip

    integer, allocatable, intent(out):: nseeds(:)

    real(kind=real64)   :: vswSI, rsSI, bemNT, TsunDAY
    namelist /input/ vswSI,rsSI,bemNT,TsunDAY,facip

    character(len=256)  :: dir
    common /dir/dir

    character(len=256)  :: f1n, f2n
    namelist/io/f1n,f2n

    character(len=256)  :: pclfil, finfil, sucfil, nodesfil, anfil, anstfil
    common /filnm/pclfil,finfil,sucfil,nodesfil,anfil,anstfil

    real(kind=real64)   :: rlambda
    common/rlambda/rlambda
    real(kind=real64)   :: rlambdax, rlambday
    common /rlambdax/rlambdax,rlambday

    integer             :: nlambda
    common /nlambda/nlambda
    integer             :: nlambdaconst
    common /nlambdaconst/nlambdaconst
    integer             :: nlambdax
    common /nlambdax/nlambdax

    integer             :: ndmumu
    common /ndmumu/ndmumu
    integer             :: ndxx
    common /ndxx/ndxx

    integer             :: ndpdt
    common /ndpdt/ndpdt
    namelist /dpdt/ndpdt

    namelist /dmumu/rlambda,nlambda,ndmumu,nlambdaconst

    namelist /dxx/rlambdax,rlambday,ndxx,nlambdax

    real(kind=real64)  :: ecrmax, ecrmin, ob0, obw, etinj, fampb, ratkp
    namelist /shkacc/ecrmax,ecrmin,ob0,obw,etinj,fampb,ratkp
    common /acc/ecrmax,ecrmin,ob0,obw,etinj,fampb,ratkp

    integer             :: nfldfile

    ! Default constants
    !
    character(len=*), parameter :: NMLFL = 'inputfld.nml' ! TODO ini this
    character(len=*), parameter :: DIRPATH_ENV = "PARAM_OUTDIR_PATH"
    character(len=*), parameter :: NODE_ENV = "PARAM_NODES"
    real(kind=real64), parameter:: daySI = 86400.0

    ! locals
    !
    integer             :: staterr
    real(kind=real64)   :: rsperdaySI
    character(len=256)  :: tmp_env_buf

    ! there are two things in dir.dat (used in old implementation)
    ! the dir where to put other data files
    ! and the number of nodes.
    ! guess what they can be. that's right, environment variables.
    ! chocolate for you!

    call get_environment_variable(DIRPATH_ENV, value=dir, status=staterr)
    if (staterr /= 0) then
      write (error_unit, '(a)') "Error reading $" // DIRPATH_ENV // " from environment"
      stop 1
    end if
    !print *, dir

    call get_environment_variable(NODE_ENV, value=tmp_env_buf, status=staterr)
    if (staterr /= 0) then
      write (error_unit, '(a)') "Error reading $" // NODE_ENV // " from environment"
      stop 1
    end if
    ! parse into int
    read (tmp_env_buf, *) nodes

    call read_seeds(nodes, nseeds)

    open(newunit=nfldfile, file=trim(dir)//NMLFL, status='old', iostat=staterr)
    read(nfldfile, nml=input)
    read(nfldfile, nml=dmumu)
    read(nfldfile, nml=dxx)
    read(nfldfile, nml=shkacc)
    read(nfldfile, nml=dpdt)
    read(nfldfile, nml=io) ! NOTE this is where pclfil and nfil comes from
    close(nfldfile)
    !
    ! magnetic field at 1AU in the equator (nT)
    !

    !write(stdout, nml=input)
    !write(stdout, nml=seeds)
    !Tsun = 26.27
    !
    ! Va=30 km/s, Rs=6.96e8 m
    !
    rsperdaySI = rsSI / daySI
    vsw = vswSI / (rsperdaySI * 24 * 60)
    omega = twopi / (TsunDAY * 24 * 60)
    vom = vsw / omega ! distance in Rs now, time in min
    if (facip == 0) facip = 1.0

    b1au = bemNT / sqrt(1.0 + 1.0/vom**2)

    pclfil   = trim(dir) // trim(f1n)
    finfil   = trim(dir) // 'finished_'
    sucfil   = trim(dir) // 'success_'
    nodesfil = trim(dir) // 'nodes'
    anfil    = trim(dir) // trim(f2n)
    anstfil  = trim(dir) // 'anst.dat'
  end subroutine

  subroutine read_seeds(n, nseeds)
    ! USES environment variables
    !     SEEDS_FILENAME

    integer, intent(in) :: n
    integer, intent(out), allocatable :: nseeds(:)


    integer             :: i, seeds_fileunit

    allocate(nseeds(n))

    call open_file_from_environment("SEEDS_FILE", seeds_fileunit, 'FORMATTED')
    do i = 1, n
      read (seeds_fileunit, *) nseeds(i)
    end do
    close(seeds_fileunit)
    return
  end subroutine


  subroutine read_b1rs(b1rsgrid)
    ! USES environment variables
    !     MAGGRID_PROCESSED_INFILENAME

    real(kind=real64), intent(out) :: b1rsgrid(0:N_R, 0:N_THETA, 0:N_PHI)

    integer :: i, j, k
    integer :: b1rs_file_unit

    call open_file_from_environment("B1RS_FILE", b1rs_file_unit, 'FORMATTED')

    do i = 0, N_R
      do j = 0, N_THETA
        do k = 0, N_PHI
          read(b1rs_file_unit,*) b1rsgrid(i,j,k)
        end do
      end do
    end do

    close(b1rs_file_unit)

  end subroutine


  subroutine read_maggrid(magfieldgrid, gbgrid)
    ! USES environment variables
    !     MAPB2S_MAGGRID_INFILENAME

    real(kind=real64), intent(out) :: magfieldgrid(0:N_R,0:N_THETA,0:N_PHI,3)
    real(kind=real64), intent(out) :: gbgrid(0:N_R,0:N_THETA,0:N_PHI,3)

    integer :: i, j, k ! iteration variables
    integer :: maggrid_fileunit
    real(kind=real64) :: r, theta, phi

    call open_file_from_environment(&
      "MAGGRID_FILE", maggrid_fileunit, 'FORMATTED')

    do i = 0, N_R
      do j = 0, N_THETA
        do k = 0, N_PHI
          read(maggrid_fileunit,*) r, theta, phi, magfieldgrid(i,j,k,:), gbgrid(i,j,k,:)
        end do
      end do
    end do

    close(maggrid_fileunit)

  end subroutine

  subroutine read_shtc(g, h, n)
    ! populate g(:,:), h(:,:)
    ! n is the size of g and h (see declaration)

    ! USES environment variables
    !     SHTC_FILE

    integer, intent(in)            :: n
    real(kind=real64), intent(out) :: g(0:n,0:n), h(0:n,0:n)

    integer :: infileunit
    integer :: l, m ! iteration variables
    integer :: ll, mm ! dummy storage

    call open_file_from_environment("SHTC_FILE", infileunit, 'FORMATTED')

    do l = 0, n
      do m = 0, l
        read (infileunit,*) ll, mm, g(l,m), h(l,m)
      end do
    end do
    close(infileunit)
  end subroutine

  subroutine write_maggrid(bgrid, gbgrid)
    ! USES environment variables
    !     MAGGRID_OUTFILENAME

    ! maggrid_omp uses its own version of bgrid and gbgrid
    real(kind=real64), intent(in) :: bgrid(0:N_R,0:N_THETA,0:N_PHI,3)
    real(kind=real64), intent(in) :: gbgrid(0:N_R,0:N_THETA,0:N_PHI,3)

    integer           :: i, j, k
    integer           :: outfileunit
    real(kind=real64) :: r, theta, phi

    call open_file_from_environment("MAGGRID_FILE", outfileunit, 'FORMATTED')
    do i = 0, N_R
      r = 1.0 + i/100.0
      do j = 0, N_THETA
        theta = j * pi/180.0
        do k = 0, N_PHI
          phi = k * pi/180.0
          write(outfileunit,*) r, theta, phi, bgrid(i,j,k,:), gbgrid(i,j,k,:)
          !write(outfileunit,*) bgrid(i,j,k,:), cvgrid(i,j,k,:), gbgrid(i,j,k,:)
          !print "(3i4,6e14.6)", i, j, k, bgrid(i,j,k,:), gbgrid(i,j,k,:)
        end do
      end do
    end do
    close(outfileunit)
  end subroutine

  subroutine write_b1rs(b1rs, map)

    ! USES environment variables
    !     MAPB2S_B1RS_OUTFILENAME

    real(kind=real64), intent(in out) :: b1rs(0:N_R, 0:N_THETA, 0:N_PHI, 2)
    integer, intent(in)               :: map(0:N_R, 0:N_THETA, 0:N_PHI)

    integer :: i, j, k ! iteration variables
    integer :: lr, ltheta, lphi ! index storage
    integer :: b1rs_outfile
    integer :: nmap

    call open_file_from_environment("B1RS_FILE", b1rs_outfile, "FORMATTED")
    do i = 0, N_R
      do j = 0, N_THETA
        do k = 0, N_PHI

          if (b1rs(i,j,k,1) == 0.0) then

            nmap = abs(map(i,j,k))

            lphi = mod(nmap, N_PHI+1)

            nmap = floor(real(nmap) / (N_PHI+1))
            ltheta = mod(nmap, N_THETA+1)

            lr = floor(real(nmap) / (N_THETA+1))

            b1rs(i,j,k,1) = sign(1, map(i,j,k)) * abs(b1rs(lr,ltheta,lphi,1))
          end if

          write (b1rs_outfile,*) b1rs(i,j,k,1)/b1rs(i,j,k,2)
        end do
      end do
    end do
    close(b1rs_outfile)
  end subroutine

end module file_op
