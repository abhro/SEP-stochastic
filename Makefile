
SHELL = /bin/bash

ifeq ($(toolset),pgi)
  #MPIBIN = /usr/local/mpich/1.2.5/ip/up/pgi/ssh/bin/
  # Requires the pgi module
  FC = $(MPIBIN)pgf95

  #FFLAGS += -ansi -O -w4
  #FFLAGS += -mp
  FFLAGS += -g -mp -mcmodel=medium

  #LDFLAGS += -Mmpi=mpich -g -O0
  #LDFLAGS += -O0 -g -mp

# default assume gcc
else
  MPIBIN =
  FC = $(MPIBIN)gfortran
  #FC = f95

  #FFLAGS += -O
  FFLAGS += -fopenmp -ggdb3
  FFLAGS += -cpp
  FFLAGS += -fcheck=all
  #FFLAGS += -Wall -Wno-unused-variable
  FFLAGS += -fbacktrace
  FFLAGS += -std=f2008
  FFLAGS += -g -O0 -fvar-tracking -fopenmp
  #FFLAGS += -g -O0 -ffpe-trap=zero,overflow,underflow -fopenmp

  #LDFLAGS += -O
  #LDFLAGS += -g -O0 -fvar-tracking
  #LDFLAGS += -g -O0 -ffpe-trap=zero,overflow,underflow
  #LDFLAGS += -O0 -g

endif

FGSL_OPTS_COMP = -I/usr/local/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-9.4.0/gsl-2.7-x2znrtoned4nhaqpluk6oh4jma6tjixe/include -I/home1/arahman2021/.local/include/fgsl

FFLAGS += $(FGSL_OPTS_COMP)

FGSL_OPTS_LINK = -L/usr/local/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-9.4.0/gsl-2.7-x2znrtoned4nhaqpluk6oh4jma6tjixe/lib \
            -L/home1/arahman2021/.local/lib \
            -lfgsl -lgsl -lgslcblas -lm

LDFLAGS += $(FGSL_OPTS_LINK)



COMPILE.f90 = $(FC) $(FFLAGS) $(TARGET_ARCH) -c
LINK.f90 = $(FC) $(FFLAGS) $(LDFLAGS) $(TARGET_ARCH)

MODULE_SRCS = param.f90 file_op.f90 epv.f90 fb.f90 dmumu.f90 \
              dxx.f90 random.f90 datetime_utils.f90 mtrx.f90 \
              loadptcl.f90 sim3d_utils.f90 cme_cross.f90 rksolvers.f90
MODULE_OBJS = $(subst .f90,.o,$(MODULE_SRCS))
MODULE_MODS = $(subst .f90,.mod,$(MODULE_SRCS))

PROG_SRCS = 2damrhistss.f90 maggrid.f90 mapb2s.f90 \
            combiner.f90 seedgen.f90 sim3d.f90 shockfront.f90
#PROG_SRCS += sim3d_em.f90
PROG_OBJS = $(subst .f90,.o,$(PROG_SRCS))
PROG_BINS = $(subst .f90,,$(PROG_SRCS))

.SUFFIXES: .o .for .f .f90 .f95 .f03
.PHONY: modules executables all clean clean-logs clean-objs
.DEFAULT_GOAL := all

# how to compile any fortran source into an object file
#
# the $(OUTPUT_OPTION) is copies to mimic make's internal
# database, even though it's never set anywhere in this file
.f90.o:
	$(COMPILE.f90) $(OUTPUT_OPTION) $<
%.o: %.f90
	$(COMPILE.f90) $< $(OUTPUT_OPTION)
%.mod: %.f90
	$(COMPILE.f90) $<

# how to compile (link) any object file to make an executable, fortran specific
%: %.o
	@#$(LINK.f90) $^ $(MODULE_OBJS) $(LDLIBS) -o $@
	@#$(LINK.f90) $^ $(MODULE_MODS) $(LDLIBS) -o $@
	$(LINK.f90) $^ $(LDLIBS) -o $@

all: $(MODULE_OBJS) executables

# these targets create executables
executables: $(MODULE_OBJS) $(PROG_BINS)

clean:
	rm -vf $(PROG_BINS) *.o *.mod

clean-logs:
	rm -vf logs/*.out.txt logs/*.err.txt

clean-objs:
	rm -vf *.o

# these targets create modules: object files (.o) and
# module files (.mod)
modules: $(MODULE_OBJS) $(MODULE_MODS)

mtrx.o: param.o epv.o dxx.o dmumu.o sim3d_utils.o
dxx.o: dmumu.o

# this depedency is only relevant when making object files
# because they need the .mod files for function signatures
$(PROG_OBJS): $(MODULE_MODS)

# link together the object file with all the modules.
# this should ideally be done with automatic dependency
# generation, since not all programs use all modules.
# alas, we don't have that in fortran
$(PROG_BINS): $(MODULE_OBJS)
