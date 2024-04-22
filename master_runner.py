#!/usr/bin/env python

import sys
# rely on short circuit evaluation so that we don't trip on Python 4.3 or such
if sys.version_info.major < 3 or sys.version_info.minor < 6:
    sys.exit("Requires Python 3.6+ to run")
import os
import datetime
import argparse
import pprint

# for submitting slurm jobs
import subprocess

# for initializing dir
import uuid
import shutil

# for mail opts
import getpass
import grp

SBATCH_DEFAULT_EXECUTABLE = '/usr/bin/sbatch'

def main():

    # get time range or manual data sources
    args = get_parser().parse_args()

    # Get SHTC files
    if args.shtc_files:
        filenames = args.shtc_files
    elif args.start_time and args.end_time:
        # get shtc from web
        filenames = download_shtc_data(args.start_time, args.end_time)
    # one of start and end time was specified, but not both
    elif args.start_time or args.end_time:
        print("-s/--start-time and -e/--end-time must be used together",
              file=sys.stderr)
        sys.exit(1)
    else:
        # error shtc_files or (start_time, end_time) not specified
        print("Neither --shtc-files or "
              "(-s/--start-time, -e/--end-time) were specified",
              file=sys.stderr)
        sys.exit(1)


    # XXX not necessary anymore
    ## make sure you're given a url to get shtc from web
    ## shtc_url = os.getenv('SHTC_URL')
    ## if shtc_url is None:
    ##     sys.exit("SHTC_URL not provided")
    ## shtc_filename = os.getenv('SHTC_FILENAME', default=f'{dir_name}/shtc.dat')


    for shtc_fn in filenames:

        # create and initialize directory
        dir_name = init_dir(shtc_fn, skeleton_root=args.skeleton_dir)

        # Generate mag data
        # run maggrid_omp to generate maggrid.dat
        if args.maggrid_file:
            maggrid_file = args.maggrid_file
   	    else:
            maggrid_file = f'{dir_name}/maggrid.dat'
            maggrid_jobid = submit_sbatch_job(
                'maggrid.sbatch',
                dir_name=dir_name,
                sbatch_executable=args.sbatch_exec,
                env_vars={
                    'SHTC_FILE': dir_name + '/' + os.path.basename(shtc_fn),
                    'MAGGRID_FILE': maggrid_file,
                }
            )
            print(f"Running maggrid_omp as job {maggrid_jobid}")

        # run mapb2s to get b1rs.dat
        mapb2s_jobid = submit_sbatch_job(
                'mapb2s.sbatch',
                dir_name=dir_name,
                sbatch_executable=args.sbatch_exec,
                dependency=maggrid_jobid,
                env_vars={
                    'MAGGRID_FILE': maggrid_file,
                    'B1RS_FILE': f'{dir_name}/b1rs.dat'
                })
        print(f"Queued mapb2s as job {mapb2s_jobid}")

        # run combiner to combine the two
        # NOTE combiner is obsolete
        #combiner_jobid = submit_sbatch_job(
        #        'combiner.sbatch', dependency=mapb2s_jobid,
        #        dir_name=dir_name,
        #        sbatch_executable=args.sbatch_exec,
        #        env_vars={
        #            'COMBINER_B1RS_INFILENAME': f'{dir_name}/b1rs.dat',
        #            'COMBINER_MAGGRID_INFILENAME': maggrid_file,
        #            'COMBINER_MAGGRID_OUTFILENAME': f'{dir_name}/maggrid_combined.dat'
        #        }
        #)
        #print(f"Queued combiner as job {combiner_jobid}")

        # Generation of mag data over

        # Generate seeds
        # Run main code
        tt1s = [0, 60, 360, 1440]
        for time in tt1s:
            sim3d_env = {
                'B1RS_FILE': f'{dir_name}/b1rs.dat',
                'MAGGRID_FILE': maggrid_file,
                'ROOT_DIR': dir_name,
                'tt1': str(time)
            }

            if os.getenv('CME_DATA_FILE'):
                sim3d_env['CME_DATA_FILE'] = os.getenv('CME_DATA_FILE')

            sim3d_jobid = submit_sbatch_job(
                    'sim3d.sbatch',
                    #dependency=combiner_jobid,
                    dependency=mapb2s_jobid,
                    dir_name=dir_name,
                    sbatch_executable=args.sbatch_exec,
                    env_vars=sim3d_env
            )
            print(f"Submitted sim3d as job {sim3d_jobid}")


def get_parser():
    root_parser = argparse.ArgumentParser()

    dl_opts_group = root_parser.add_argument_group('download options')

    dl_opts_group.add_argument(
            '-s', '--start-time',
            help='Starting period of search duration. Requires --end-time. '
                 'Passed to sunpy.net.attrs.Time, so the string must '
                 'conform to its initializer.')

    dl_opts_group.add_argument(
            '-e', '--end-time',
            help='Ending period of search duration. Requires --start-time. '
                 'Passed to sunpy.net.attrs.Time, so the string must '
                 'conform to its initializer.')


    root_parser.add_argument(
            '--shtc-files', nargs='+',
            help='Specify the location of the shtc files. '
                 'To be used if download through sunpy is not possible. '
                 'If set, --start-time and --end-time are ignored.')

    root_parser.add_argument(
            '--sbatch-exec', default=SBATCH_DEFAULT_EXECUTABLE,
            help='Alternative sbatch binary executable')

    root_parser.add_argument(
            '--skeleton-dir', default='./skeleton-dir/',
            help='Alternative path to skeleton directory')

    root_parser.add_argument(
        	'--maggrid-file',
        	help='Pre-generated maggrid file (skips running maggrid)')

    return root_parser

def submit_sbatch_job(script_fn: str,
                      dir_name: str,
                      dependency: int = None,
                      env_vars: {str: str} = None,
                      sbatch_executable: str = SBATCH_DEFAULT_EXECUTABLE) -> int:

    mail_options = get_mail_opts()
    log_options = get_logging_opts(dir_name)

    args = [sbatch_executable, '--parsable', *mail_options, *log_options]
    if dependency is not None:
        args.append(f"--dependency=afterok:{dependency}")
    args.append(script_fn)

    if env_vars is None:
        env_vars = {"PATH": os.getenv("PATH")}
    else:
        env_vars["PATH"] = os.getenv("PATH")

    #print("submit_sbatch_job:{args = }")
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, env=env_vars)

    jobid, _ = proc.communicate()

    return int(jobid.decode().strip())

def download_shtc_data(start_time, end_time) -> [str]:

    # for downloading data
    try:
        from sunpy.net import Fido, attrs as a
        import gong_shtc
    except ImportError as ex:
        print("Could not import sunpy:", file=sys.stderr)
        print(ex, file=sys.stderr)
        print("Please provide the SHTC data files manually with the "
              "--shtc-files flag", file=sys.stderr)
        raise

    # search for the data
    results = Fido.search(a.Time(start_time, end_time),
                          a.Instrument('GONG_SHTC'))

    print(f'Found {results.file_num} in time range {start_time} to {end_time}')
    # download the data
    downloaded_files = Fido.fetch(results)

    return downloaded_files

def init_dir(shtc_path: str, skeleton_root: str) -> str:
    # each run gets a different id to run under
    # sth sth uuid
    #data_dir = 'soho'
    # everything is contained in a directory named (time now + random hex string)
    dir_name = (
            datetime.datetime.now().strftime('%Y%m%d-T%H%M%S')
            + '-'
            + uuid.uuid4().hex)

    print(f'Initializing directory {dir_name} for {shtc_path}')
    # make directory for this shtc file to live in
    os.mkdir(dir_name)
    print(f'Created directory {dir_name} for {shtc_path}')

    shutil.copy(shtc_path, dir_name)

    # copy over contents from skeleton directory:
    # skeleton directory
    # |- shtc file
    # |- ...
    # |- run-1/             # (for sim3d)
    # | |- flux.dat
    # | |- files.nml
    # | |- loadptcl.dat
    # | `- other stuff
    # |- run-2/             # (for sim3d)
    # | |- flux.dat
    # | |- files.nml
    # | |- loadptcl.dat
    # | `- other stuff
    # `- run-3/             # (for sim3d)
    # | |- flux.dat
    # | |- files.nml
    # | |- loadptcl.dat
    #   `- other stuff
    shutil.copytree(src=skeleton_root, dst=dir_name, dirs_exist_ok=True)

    print(f'Finished initializing directory "{dir_name}" for {shtc_path}')
    return dir_name


def get_mail_opts() -> (str, str):
    username = getpass.getuser()

    #student_gid = grp.getgrnam('student').gr_gid # bottleneck
    #faculty_gid = grp.getgrnam('faculty').gr_gid # bottleneck
    #staff_gid = grp.getgrnam('staff').gr_gid     # bottleneck
    student_gid = grp.getgrnam('student@fit.edu').gr_gid # bottleneck
    faculty_gid = grp.getgrnam('faculty@fit.edu').gr_gid # bottleneck
    staff_gid = grp.getgrnam('staff@fit.edu').gr_gid     # bottleneck

    pgid = os.getgid()

    if pgid == student_gid: # user is a student
        domain = 'my.fit.edu'
    elif pgid == faculty_gid or pgid == staff_gid: # user is faculty or staff
        domain = 'fit.edu'

    username = username.replace('@fit.edu','')

    email_address = f"{username}@{domain}"
    mail_frequency = "all,time_limit,array_tasks"

    mail_options = (f"--mail-type={mail_frequency}",
                    f"--mail-user={email_address}")
    return mail_options

def get_logging_opts(dir_name) -> (str, str):
    return ('--error=' + dir_name + '/%x.%J.err.txt',
            '--output=' + dir_name + '/%x.%J.out.txt')

if __name__ == '__main__':
    main()
