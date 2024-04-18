# SEP-stochastic
Fortran numerical integrator for SEP stochastic model

## Setting up virtual env

To set up virtual env

    $ module load python/3.9.0
    $ python -m pip install virtualenv
    $ # if venv-python doesn't already exist
    $ python -m virtualenv venv-python
    $ source venv-python/bin/activate

## Compiling
To compile, run `make all`, or just `make`. If there are some dependency errors,
you can use `make modules` to compile the modules first before redoing `make all`.
If errors still persist, use `make clean` to start from scratch. Current
versions of Fortran compilers don't support automatic dependency generation, so
some dependencies, like for mtrx.o and dxx.o have to be defined manually.

To compile using GFortran/GCC, set the environment variable `$toolset` to `gcc`.
By default the Makefile will assme a PGI toolset is being used.

    $ export toolset=gcc
    $ make clean
    $ make modules
    $ make all

## Running

INVALID. Does not work for the python script.

First, run `export SHTC_URL=<url from GONG here>` to set the environment
variable to the URL, the URL from which the spherical harmonic transform
coefficients will be downloaded. Then, run `./master-runner.sh` and sit back and
wait for it to email you. That is,

    $ export SHTC_URL=https://gong.nso.edu/data/magmap/QR/bqc/202004/mrbqc200421/mrbqc200421t0054c2230_348.dat
    $ ./master-runner.sh

This can also be combined into one line by setting the variable in the same line
right before the script invocation, like so:

    $ SHTC_URL=https://gong.nso.edu/data/magmap/QR/bqc/202004/mrbqc200421/mrbqc200421t0054c2230_348.dat ./master-runner.sh
where a SHTC file from 2020 April is used as an example.

If you are Florida Tech faculty/staff, you will want to change the line
containing the variable `$MAIL_ADDRESS` from

    MAIL_ADDRESS="$USER@my.fit.edu"
to

    MAIL_ADDRESS="$USER@fit.edu"
i.e., remove the 'my.' domain, to ensure you receive emails. This *could* be
automated by checking if `$HOME` contains `/udrive/faculty`, `/udrive/student`,
or `/udrive/staff`, but the author is not familiar enough with the email address
protocols of Florida Tech to reliably automate that in a future-proof fashion.

## Cleaning
There are three cleaning targets within the Makefile, `clean`, `clean-logs`, and
`clean-objs`. `make clean` will remove all object files, module files, and
compiled executables. `make clean-objs` will only remove object (.o) files,
leaving the executables and module (.mod) files intact. Finally,
`make clean-logs` will only remove the output and error files within the logs/
directory.

## TODO
- [ ] Get rid of control through environment variables and implement a conf file
  reader, a la [Finer](https://github.com/szaghi/FiNeR)
- [ ] Figure out what 2damrhistss is trying to do, and redirect I/O accordingly.
